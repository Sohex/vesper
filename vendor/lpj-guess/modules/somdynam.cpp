///////////////////////////////////////////////////////////////////////////////////////
/// \file somdynam.cpp
/// \brief Soil organic matter dynamics
///
/// \author Ben Smith (LPJ SOM dynamics, CENTURY), David Wårlind (CENTURY), Mateus Dantas (P-Cycle)
/// $Date: 2022-09-13 10:47:57 +0200 (Tue, 13 Sep 2022) $
///
///////////////////////////////////////////////////////////////////////////////////////


// WHAT SHOULD THIS FILE CONTAIN?
// Module source code files should contain, in this order:
//   (1) a "#include" directive naming the framework header file. The framework header
//       file should define all classes used as arguments to functions in the present
//       module. It may also include declarations of global functions, constants and
//       types, accessible throughout the model code;
//   (2) other #includes, including header files for other modules accessed by the
//       present one;
//   (3) type definitions, constants and file scope global variables for use within
//       the present module only;
//   (4) declarations of functions defined in this file, if needed;
//   (5) definitions of all functions. Functions that are to be accessible to other
//       modules or to the calling framework should be declared in the module header
//       file.
//
// PORTING MODULES BETWEEN FRAMEWORKS:
// Modules should be structured so as to be fully portable between models (frameworks).
// When porting between frameworks, the only change required should normally be in the
// "#include" directive referring to the framework header file.

#include "config.h"
#include "somdynam.h"
#include "ntransform.h"
#include "driver.h"
#include <assert.h>
#include <bitset>
#include <vector>

///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL CONSTANTS

// Turnover times for litter and SOM fractions at 10 deg C with ample moisture
// (Meentemeyer 1978; Foley 1995). In EARTH years: these are decomposition times
// measured on Earth, and decomposition does not know this world's orbit. They
// are used only by som_dynamics_lpj, the ifcentury 0 path; the live CENTURY path
// takes its K_MAX from Parton et al. (2010) already on a daily basis, which is
// absolute time and needs no conversion. See
// biosphere/notes/time-base-unit-contract.md.

static const double TAU_LITTER=2.85; // Thonicke, Sitch, pers comm, 26/11/01
static const double TAU_SOILFAST=33.0;
static const double TAU_SOILSLOW=1000.0;

static const double TILLAGE_FACTOR = 33.0 / 17.0;
	// (Inverse of "tillage factor" in Chatskikh et al. 2009, Value selected by T.Pugh)

static const double FASTFRAC=0.985;
	// fraction of litter decomposition entering fast SOM pool
static const double ATMFRAC=0.7;
	// fraction of litter decomposition entering atmosphere

// Corresponds to the amount of soil available nitrogen where SOM C:N ratio reach
// their minimum (nitrogen saturation) (Parton et al 1993, Fig. 4)
// Comment: NMASS_SAT is too high when considering BNF - Zaehle
static const double NMASS_SAT = 0.002 * 0.05;
// Corresponds to the nitrogen concentration in litter where SOM C:N ratio reach
// their minimum (nitrogen saturation) (Parton et al 1993, Fig. 4)
static const double NCONC_SAT = 0.02;

// The phosphorus saturation pair. Each is the fmax argument of setptoc(), the
// value of the driving quantity at which a SOM pool's C:P reaches its minimum.
// One of them is now derived and the other is not, and they fail differently.
//
// PMASS_SAT IS Parton, Stewart and Cole (1988), read directly off Fig. 3,
// p. 115: that figure plots the C:P ratio of the passive, slow and active soil
// organic matter against soil labile P over an axis running 0 to 2.0 gP/m2,
// and every line reaches its minimum at the right-hand end of that axis.
// 2.0 gP/m2 is 0.002 kgP/m2, the value here. The three (ctop_max, ctop_min)
// pairs passed to setptoc() below are that figure's three lines exactly: slow
// 200 to 90, passive 200 to 20, active 80 to 30. setptoc's linear
// interpolation is the figure's straight lines, and it sets the P:C of the
// pool RECEIVING carbon, which is the paper's own statement of how the organic
// P flows are computed (p. 115). So the constant is derived, the functional
// form is the source's, and the resemblance to the 0.002 inside NMASS_SAT is a
// coincidence of two unrelated readings.
//
// What is NOT derived is the pool this threshold is applied to. Parton's
// labile P is orthophosphate that is isotopically exchangeable or
// anion-exchange-resin extractable, over 0 to 20 cm (p. 115). This fork's
// soil.pmass_labile is the Hedley-labile pool of Wang et al. (2010) and Yang
// et al. (2013) -- labile inorganic PLUS labile organic P -- and Dantas de
// Paula et al. (2025) adopt that wider definition deliberately, to stand in
// for the biomineralisation pathway the model does not carry. Their own global
// numbers separate the two: 2.11 PgP simulated labile P and 3.6 PgP
// Hedley-labile, against 0.319 PgP by Olsen extraction over 0 to 20 cm. Over
// 1.3e14 m2 of ice-free land those are 16, 28 and 2.4 gP/m2, and Parton's
// 2.0 gP/m2 sits on the last of them. The fork's own sorption parameters say
// the same thing from the inside: kplab, the Langmuir half-saturation for
// labile P, is 10 to 78 gP/m2 by soil order (Wang et al. 2010, Table A1), so
// Parton's unconverted axis maximum is 5 to 39 times below the level at which
// the fork's own isotherm expects labile P to be halfway to saturating.
//
// Unconverted, fac exceeds fmax nearly everywhere and the slow, passive and
// soil microbial pools sit at their MINIMUM C:P always, which is their most
// phosphorus-rich end. That is not a wrong constant. It is a right constant
// reading a pool its source did not define. The one reader of what setptoc
// writes is transferdecomp()'s pinc, the phosphorus that rides carbon into a
// receiving pool, so this decides the phosphorus content of every organic
// transfer in the model.
//
// DECLARED DIVERGENCE FROM MAINLINE: pmass_sat_labile_currency, owner
// WORLD-Z01O. Vendored LPJ-GUESS-CNP declares
//     static const double PMASS_SAT = 0.002;
// which is Parton's axis maximum carried across unconverted. This project
// converts the THRESHOLD into the fork's own labile-P currency by the factor
// 6.6 applied below, giving 0.0132 kgP/m2. biosphere/config/somdynam.yaml is
// the register, on the ntransform.yaml pattern, and
// biosphere/scripts/somdynam_gate.py checks that this comment records the
// fork's line, that the fork's line is not what runs, and that the converted
// one is. The reference point there is the vendored subtree commit and not a
// release: LPJ-GUESS 4.1.1 has no phosphorus. The argument is
// biosphere/notes/phosphorus-cycle-parameterisation.md.
//
// The threshold and not the driver, because PMASS_SAT does a SECOND job:
// somfluxes() ends by pinning soil.pmass_labile to PMASS_SAT whenever !ifplim,
// and that pin is DEFINED as holding labile P at the value where the C:P ramp
// stops responding. The two uses are one number rather than two that coincide,
// so converting the threshold and letting the pin follow is one definition
// propagating and not a side effect. Splitting them into two constants would
// invent a second number with no independent derivation and let the two drift.
// The alternative -- leaving the threshold alone and driving setptoc() with a
// resin-equivalent FRACTION of soil.pmass_labile at the call sites only --
// would leave the pin behind and is not taken for that reason.
//
// 6.6 and not the middle of the bracket. The ratio between this fork's
// pmass_labile and the resin-extractable orthophosphate Fig. 3 was drawn
// against is bracketed 6.6 to 11.3 by Dantas de Paula et al. (2025)'s own
// global numbers: 6.6 is the model's SIMULATED labile P, 2.11 PgP, over
// Olsen-extractable 0.319 PgP; 11.3 is the OBSERVATIONAL Hedley-labile
// estimate, 3.6 PgP, over the same Olsen figure. setptoc() consumes the
// simulated pool and not the observation, so 6.6 is the ratio between the two
// quantities that are actually wired together. 11.3 would additionally carry
// the model's 41 per cent under-prediction of its own observational target,
// counting that error twice. The bracket cannot be tightened further from
// these numbers, because LPJ-GUESS-CNP's soil organic matter is a bulk pool
// with no depth and Parton's figure is per 0 to 20 cm, so the two sides of the
// ratio do not share a support; narrowing it needs a run of this fork.
//
// What moves. Under ifplim 0, the configuration this project runs, the pin
// moves with the threshold, fac still arrives at fmax, and every soil organic
// C:P ratio is unchanged: the divergence is INERT for the carbon and
// phosphorus flows. What does change is the reported stock, because
// commonoutput.cpp writes the pinned pmass_labile into PO4_mass and availp --
// the reported labile P moves from 2.0 to 13.2 gP/m2, which is the same order
// as the fork's own simulated 16 gP/m2 instead of an order below it. It is
// still a constant meaning "not limiting" and still not a simulated stock.
// Under ifplim 1 the divergence bites for real: the emergent labile P now has
// a threshold it can sit below, so the three pools ramp instead of saturating.
//
// PCONC_SAT has no phosphorus source at all. It carries NCONC_SAT's 0.02
// exactly, and Parton, Stewart and Cole (1988) contains no counterpart to it:
// that model has no surface microbial pool and no C:P ramp driven by a litter
// concentration. Its only C:P ramps are the three soil pools above, driven by
// labile P; litter P is set by a fixed structural C:P of 500 with the
// remainder going to the metabolic pool (p. 115). So both the ramp PCONC_SAT
// belongs to and its value are nitrogen's, carried across.
//
// It is compared against litter_pmass / (litter_cmass * 2), a phosphorus
// fraction of litter dry mass, whose whole attainable range is 3.1e-4 to
// 7.6e-4 for senesced-litter C:P of 1596 to 660 by mass (McGroddy et al. 2004,
// Table 1). 0.02 is 26 to 64 times above the richest litter the model can
// make, so the surface microbial pool sits at its MAXIMUM C:P of 80 always.
// That bounds any replacement from above; nothing in the cited source anchors
// it from below. WORLD-PIDX.
//
// PCONC_SAT is not changed here, because it is not settled by arithmetic: it
// has no phosphorus source at all, and parameters.cpp goes on refusing
// ifplim 1 while that stands. biosphere/notes/phosphorus-cycle-parameterisation.md
// carries the evidence, the arithmetic and what each is worth in the model's
// own reported stocks.
//
// 0.002 kgP/m2 is Fig. 3's axis maximum; 6.6 converts it into this fork's
// labile-P currency, as argued above.
static const double PMASS_SAT = 0.002 * 6.6;
static const double PCONC_SAT = 0.02;

// Phosphorus sorption rate constants, Wang et al. (2007) as restated by
// Wang et al. (2010) Appendix D.
//
// USORB and USSORB are EQUAL IN THE SOURCE, not by a copy here: Wang et al.
// (2010) state that the rate constants for the sorbed and strongly sorbed P
// pools "both are equal to 0.0067 year-1". Do not read the equality as a defect
// and do not split the two without a source that measures them apart. What
// follows from it is that Eq. D10, dPssb/dt = USORB*Psorb - USSORB*Pssb, drives
// the strongly sorbed pool to exactly the size of the sorbed pool and holds it
// there, so the strongly sorbed pool is a stock and not a sink. Nothing drains
// it: this model has no terminal occlusion, which is a declared gap and not a
// slow process. biosphere/notes/phosphorus-cycle-parameterisation.md argues the
// decision, and parameters.cpp refuses ifplim 1 while it stands.
//
// ABSOLUTE-RATE, per EARTH year: sorption is chemistry and does not know this
// world's orbit, so the divisor is the Earth year and not the simulation year.
// Dividing by date.year_length() delivered an Earth year of sorption every
// orbit, which is close to twice the published rate per unit absolute time.
// biosphere/notes/time-base-unit-contract.md.
//
// The reader is somfluxes(), which is called once per absolute day on the live
// daily path. equilsom() also calls it twelve times per model year with
// monthly-aggregated decay rates, so there these two constants act at 12/365 of
// their intended speed. That is a rate of approach and not an equilibrium:
// USORB == USSORB fixes the equilibrium at Pssb = Psorb whatever the constants
// are, and equilsom's 40000 model years reach it either way.
static const double USORB = 0.0067 / VESPER_EARTH_YEAR_DAYS;
static const double USSORB = 0.0067 / VESPER_EARTH_YEAR_DAYS;

///////////////////////////////////////////////////////////////////////////////////////
// FILE SCOPE GLOBAL VARIABLES

// Exponential decay constants for litter and SOM fractions
// Values set from turnover times (constants above) on first call to decayrates

static double k_litter10;
static double k_soilfast10;
static double k_soilslow10;

static bool firsttime=true;
	// indicates whether function decayrates has been called before

int dummy = 1;

///////////////////////////////////////////////////////////////////////////////////////
// SETCONSTANTS
// Internal function (do not call directly from framework)

void setconstants() {

	// DESCRIPTION
	// Calculate exponential decay constants (annual basis) for litter and
	// SOM fractions first time function decayrates is called

	k_litter10=1.0/TAU_LITTER;
	k_soilfast10=1.0/TAU_SOILFAST;
	k_soilslow10=1.0/TAU_SOILSLOW;
	firsttime=false;
}


///////////////////////////////////////////////////////////////////////////////////////
// BALANCE LABILE AND SORBED P POOLS
// Internal function (do not call directly from framework)

void pmass_add(Soil &soil, double delta) {

	double a = -1.0;
	double b = -soil.pmass_labile - soil.soiltype.kplab - soil.soiltype.spmax + soil.pmass_sorbed + delta;
	double c = soil.pmass_sorbed * soil.pmass_labile + soil.pmass_sorbed  * soil.soiltype.kplab + delta * soil.pmass_labile + delta * soil.soiltype.kplab - soil.pmass_labile * soil.soiltype.spmax;

	double bha = std::max(0.0, std::pow(b, 2.0) - 4.0 * a * c);

	double labile_inc = (-b - std::sqrt(bha)) / 2.0 * a;

	double sorbed_inc = delta - labile_inc;

		soil.pmass_labile += labile_inc;
	
		soil.pmass_sorbed += sorbed_inc;

}

/// Fraction of a stock surviving one daily removal operator.
/** The accelerator composes these fractions multiplicatively. Averaging daily
 *  removal fractions and applying the mean once is not the daily operator.
 */
double stock_survival(double before, double after) {
	if (negligible(before))
		return 1.0;
	return min(1.0, max(0.0, after / before));
}

/// Apply a composed survival fraction to the exchangeable P stock.
/** Every P mutation goes through pmass_add(), so the labile--sorbed isotherm
 *  remains true after accelerated uptake and leaching just as on the daily
 *  path.
 */
void pmass_apply_survival(Soil& soil, double survival) {
	const double exchangeable = soil.pmass_labile + soil.pmass_sorbed;
	pmass_add(soil, -exchangeable * (1.0 - survival));
}


//void pmass_add(Soil &soil, double delta) {
//
//	if (delta == 0.0)
//		return;
//
//	double k = soil.soiltype.kplab;
//	double s = soil.soiltype.spmax;
//
//	double total_p = soil.pmass_labile + soil.pmass_sorbed + delta;
//	
//	if (total_p < 1e-20) {
//		soil.pmass_labile = 0.0;
//		soil.pmass_sorbed = 0.0;
////		soil.patch.fluxes.report_flux(Fluxes::P_SOIL, soil.pmass_labile + soil.pmass_sorbed - delta);
//	}
//	else {
//
//		double b = -(1 - k / total_p - s / total_p);
//		double c = -k / total_p;
//		double bha = b * b - 4.0 * c;
//
//		double f11 = (-b + std::sqrt(bha)) / 2.0;
//
//		if ((1 - f11) * total_p < soil.soiltype.spmax) {
//			soil.pmass_labile = f11 * total_p;
//			soil.pmass_sorbed = (1 - f11) * total_p;
//		}
//		else {
//			soil.pmass_labile += delta;
//			soil.pmass_sorbed = soil.soiltype.spmax;
//		}
//
//
//		//if (soil.pmass_labile > PMASS_SAT) {
//		//	soil.pmass_labile = PMASS_SAT;
//		//	soil.patch.fluxes.report_flux(Fluxes::P_SOIL, soil.pmass_labile - PMASS_SAT);
//		//}
//
//		//if (soil.pmass_sorbed > (PMASS_SAT * soil.soiltype.spmax) / (soil.soiltype.kplab + PMASS_SAT)) {
//		//	soil.pmass_sorbed = (PMASS_SAT * soil.soiltype.spmax) / (soil.soiltype.kplab + PMASS_SAT);
//		//	soil.patch.fluxes.report_flux(Fluxes::P_SOIL, soil.pmass_sorbed - (PMASS_SAT * soil.soiltype.spmax) / (soil.soiltype.kplab + PMASS_SAT));
//		//}
//
//	}
//
//}



///////////////////////////////////////////////////////////////////////////////////////
// DECAYRATES
// Internal function (do not call directly from framework)
// used by som_dynamic_lpj()

void decayrates(double wcont,double gtemp_soil,double& k_soilfast,double& k_soilslow,
	double& fr_litter,double& fr_soilfast,double& fr_soilslow, bool tillage) {

	// DESCRIPTION
	// Calculation of fractional decay amounts for litter and fast and slow SOM
	// fractions given current soil moisture and temperature

	// INPUT PARAMETERS
	// wcont       = water content of upper soil layer (fraction of AWC)
	// gtemp_soil  = respiration temperature response incorporating damping of Q10
	//               response due to temperature acclimation (Eqn 11, Lloyd & Taylor
	//               1994)

	// OUTPUT PARAMETERS
	// k_soilfast  = adjusted daily decay constant for fast SOM fraction
	// k_soilslow  = adjusted daily decay constant for slow SOM fraction
	// fr_litter   = litter fraction remaining following today's decomposition
	// fr_soilfast = fast SOM fraction remaining following today's decomposition
	// fr_soilslow = slow SOM fraction remaining following today's decomposition

	double moist_response; // moisture modifier of decomposition rate

	// On first call only: set exponential decay constants

	if (firsttime) setconstants();

	// Calculate response of soil respiration rate to moisture content of upper soil layer
	// Foley 1995 Eqn 19

	moist_response=0.25+0.75*wcont;

	// Calculate litter and SOM fractions remaining following today's decomposition
	// (Sitch et al 2000 Eqn 71) adjusting exponential decay constants by moisture and
	// temperature responses and converting from annual to daily basis
	// NB: Temperature response (gtemp; Lloyd & Taylor 1994) set by framework

	// Annual to daily on ABSOLUTE time: the TAU constants above are Earth years,
	// so the divisor is the Earth year. Dividing by the simulation year instead
	// ran every litter and SOM pool through an Earth year of decomposition each
	// orbit, which on this calendar is close to twice as fast per unit absolute
	// time and left equilibrium soil carbon at about half what the same litter
	// input supports on Earth.
	k_soilfast=k_soilfast10*gtemp_soil*moist_response/VESPER_EARTH_YEAR_DAYS;
	if (tillage) {
		k_soilfast *= TILLAGE_FACTOR; // Increased HR for crops (tillage)
	}
	k_soilslow=k_soilslow10*gtemp_soil*moist_response/VESPER_EARTH_YEAR_DAYS;

	fr_litter=exp(-k_litter10*gtemp_soil*moist_response/VESPER_EARTH_YEAR_DAYS);
	fr_soilfast=exp(-k_soilfast);
	fr_soilslow=exp(-k_soilslow);
}


///////////////////////////////////////////////////////////////////////////////////////
// DECAYRATES
// Should be called by framework on last day of simulation year, following call to
// som_dynamics, once annual litter production and vegetation PFT composition are close
// to their long term equilibrium (typically 500-1000 simulation years).
// NB: should be called ONCE ONLY during simulation for a particular grid cell

void equilsom_lpj(Soil& soil) {

	// DESCRIPTION
	// Analytically solves differential flux equations for fast and slow SOM pools
	// assuming annual litter inputs close to long term equilibrium

	// INPUT PARAMETER (class defined in framework header file)
	// soil = current soil status

	double nyear;
		// number of years over which decay constants and litter inputs averaged

	nyear=soil.soiltype.solvesom_end-soil.soiltype.solvesom_begin+1;

	soil.decomp_litter_mean/=nyear;
	soil.k_soilfast_mean/=nyear;
	soil.k_soilslow_mean/=nyear;

	soil.cpool_fast=(1.0-ATMFRAC)*FASTFRAC*soil.decomp_litter_mean/
		soil.k_soilfast_mean;
	soil.cpool_slow=(1.0-ATMFRAC)*(1.0-FASTFRAC)*soil.decomp_litter_mean/
		soil.k_soilslow_mean;
}


///////////////////////////////////////////////////////////////////////////////////////
// SOM DYNAMICS
// To be called each simulation day for each modelled area or patch, following update
// of soil temperature and soil water.

void som_dynamics_lpj(Patch& patch, bool tillage) {

	// DESCRIPTION
	// Calculation of soil decomposition and transfer of C between litter and soil
	// organic matter pools.

	double k_soilfast; // adjusted daily decay constant for fast SOM fraction
	double k_soilslow; // adjusted daily decay constant for slow SOM fraction
	double fr_litter;
		// litter fraction remaining following one day's/one month's decomposition
	double fr_soilfast;
		// fast SOM fraction remaining following one day's/one month's decomposition
	double fr_soilslow;
		// slow SOM fraction remaining following one day's/one month's decomposition
	double decomp_litter; // litter decomposition today/this month (kgC/m2)
	double cflux; // accumulated C flux to atmosphere today/this month (kgC/m2)
	int p;

	// Obtain reference to Soil object
	Soil& soil=patch.soil;

	// Calculate decay constants and rates given today's soil moisture and
	// temperature

	decayrates(soil.get_soil_water_upper(),soil.gtemp,k_soilfast,k_soilslow,fr_litter,
		fr_soilfast,fr_soilslow, tillage);

	// From year soil.solvesom_begin, update running means for later solution
	// (at year soil.solvesom_end) of equilibrium SOM pool sizes

	if (date.year>=soil.soiltype.solvesom_begin) {
		soil.k_soilfast_mean+=k_soilfast;
		soil.k_soilslow_mean+=k_soilslow;
	}

	// Reduce litter and SOM pools, sum C flux to atmosphere from decomposition
	// and transfer correct proportions of litter decomposition to fast and slow
	// SOM pools

	// Reduce individual litter pools and calculate total litter decomposition
	// for today/this month

	decomp_litter=0.0;

	// Loop through PFTs

	for (p=0;p<npft;p++) {	//NB. also inactive pft's

		// For this PFT ...

		decomp_litter+=(patch.pft[p].cmass_litter_leaf+
			patch.pft[p].cmass_litter_root+
			patch.pft[p].cmass_litter_sap+
			patch.pft[p].cmass_litter_heart+
			patch.pft[p].cmass_litter_repr)*(1.0-fr_litter);

		patch.pft[p].cmass_litter_leaf*=fr_litter;
		patch.pft[p].cmass_litter_root*=fr_litter;
		patch.pft[p].cmass_litter_sap*=fr_litter;
		patch.pft[p].cmass_litter_heart*=fr_litter;
		patch.pft[p].cmass_litter_repr*=fr_litter;
	}

	if (date.year>=soil.soiltype.solvesom_begin)
		soil.decomp_litter_mean+=decomp_litter;

	// Partition litter decomposition among fast and slow SOM pools
	// and flux to atmosphere

	// flux to atmosphere
	cflux=decomp_litter*ATMFRAC;

	// remaining decomposition - goes to ...
	decomp_litter-=cflux;

	// ... fast SOM pool ...
	soil.cpool_fast+=decomp_litter*FASTFRAC;

	// ... and slow SOM pool
	soil.cpool_slow+=decomp_litter*(1.0-FASTFRAC);

	// Increment C flux to atmosphere by SOM decomposition
	cflux+=soil.cpool_fast*(1.0-fr_soilfast)+soil.cpool_slow*(1.0-fr_soilslow);

	// Reduce SOM pools
	soil.cpool_fast*=fr_soilfast;
	soil.cpool_slow*=fr_soilslow;

	// Updated soil fluxes. In wetlands, some portion of cflux can be emitted as CH4, so save cflux as dcflux_soil until Soil::methane() is called.  
	if (patch.stand.landcover != PEATLAND)
		patch.fluxes.report_flux(Fluxes::SOILC, cflux);
	else 
		soil.dcflux_soil=cflux;

	// Solve SOM pool sizes at end of year given by soil.solvesom_end

	if (date.year==soil.soiltype.solvesom_end && date.islastmonth && date.islastday)
		equilsom_lpj(soil);

}

/////////////////////////////////////////////////
// CENTURY SOM DYNAMICS

/// Data type representing a selection of SOM pools
/** A selection of SOM pools is represented by a bitset,
 *  the selected pools have their corresponding bits switched on.
 */
typedef std::bitset<NSOMPOOL> SomPoolSelection;


/// Reduce decay rates to keep the daily nitrogen balance in the soil
/** Only a selected subset of the SOM pools (as specified by the caller),
 *  are considered for reducion of decay rates.
 *  Usually the first time is enough (decay rate reduction of litter).
 *  Can be for phosphorus or nitrogen
 */
void reduce_decay_rates(double decay_reduction[NSOMPOOL], double net_min_pool[NSOMPOOL], const SomPoolSelection& selected, double neg_nmass_avail) {

	int neg_min_pool[NSOMPOOL] = {0};	// Keeping track on which pools that are negative

	// Add up immobilization for considered pools
	double tot_neg_min = 0.0;
	for (int p = 0; p < NSOMPOOL; p++) {
		if (selected[p] && net_min_pool[p] < 0.0) {
			tot_neg_min += net_min_pool[p];
			neg_min_pool[p] = 1;
		}
	}

	// Calculate decay reduction
	double decay_red = 0.0;

	if (tot_neg_min < neg_nmass_avail) {
		// Enough to reduce decay rates for these pools to achieve a
		// net positive mineralization
		decay_red = 1.0 - (tot_neg_min - neg_nmass_avail) / tot_neg_min;
	}
	else {
		// Need to stop these pools from decaying to be able to get
		// a net positive mineralization
		decay_red = 1.0;
	}

	// Reduce decay rate for considered pools
	for (int p = 0; p < NSOMPOOL; p++) {
		if (selected[p]) {
			decay_reduction[p] = decay_red * neg_min_pool[p];
		}
	}
}

/// Set N:C ratios for SOM pools
/** Set N:C ratios for slow, passive, humus and soil microbial pools
 *  based on mineral nitrogen pool or litter nitrogen fraction (Parton et al 1993, Fig 4)
 */
void setntoc(Soil& soil, double fac, pooltype pool, double cton_max, double cton_min,
	double fmin, double fmax) {

	if (fac <= fmin)
		soil.sompool[pool].ntoc = 1.0 / cton_max;
	else if (fac >= fmax)
		soil.sompool[pool].ntoc = 1.0 / cton_min;
	else {
		soil.sompool[pool].ntoc = 1.0 / (cton_min + (cton_max - cton_min) *
			(fmax - fac) / (fmax - fmin));
	}
}

/// Set P:C ratios for SOM pools
/** Set P:C ratios for the slow, passive and soil microbial pools from the
*  labile P pool, and for the surface microbial pool from the litter P
*  fraction. The first is Parton, Stewart and Cole (1988) Fig. 3, p. 115,
*  whose three straight lines are the (ctop_max, ctop_min) pairs the callers
*  pass and whose axis maximum is PMASS_SAT. The second has no counterpart in
*  that paper; see the PCONC_SAT comment at the top of this file. The ratio
*  set here is the P:C of the pool RECEIVING carbon, which is the paper's own
*  construction, and it is applied in transferdecomp() below.
*/
void setptoc(Soil& soil, double fac, pooltype pool, double ctop_max, double ctop_min,
	double fmin, double fmax) {

	if (fac <= fmin)
		soil.sompool[pool].ptoc = 1.0 / ctop_max;
	else if (fac >= fmax)
		soil.sompool[pool].ptoc = 1.0 / ctop_min;
	else {
		soil.sompool[pool].ptoc = 1.0 / (ctop_min + (ctop_max - ctop_min) *
			(fmax - fac) / (fmax - fmin));
	}
}


/// Temperature modifier for decomposition
/** Calculate decomposition temperature modifier (in range 0-1)
 *  [A(T_soil), Eqn A9, Comins & McMurtrie 1993; ET, Friend et al 1997; abiotic
 *  effect of soil temperature, Parton et al 1993, Fig 2)
 *
 *  \param temp_soil  Soil temperature at 25 cm depth
 */
double temperature_modifier(double temp_soil) {

	double temp_mod = temp_soil > 0.0 ? max(0.0, 0.0326 + 0.00351 * pow(temp_soil, 1.652) - pow(temp_soil / 41.748, 7.19)) : 0.0;

	// Include as an overide option when: MIN_DECOMP_TEMP < temp < 0 degC
	// This increases the respiration (from 0) between MIN_DECOMP_TEMP and 0 degC - cf Koven et al. 2011
	if (ifcarbonfreeze && temp_soil <= 0.0 && temp_soil >= MIN_DECOMP_TEMP && !iftwolayersoil) {

		double decomp_at_freezing_point = 0.0326; // temp_mod above when temp_soil = 0;
		bool linear_decrease_below_freezing = false;

		if (linear_decrease_below_freezing) {
			// Alternative: Linear approach (Koven et al. 2011)
			double slope = decomp_at_freezing_point / fabs(MIN_DECOMP_TEMP);
			temp_mod = slope * temp_soil + decomp_at_freezing_point; // i.e. a linear decrease from decomp_at_freezing_point at 0C to 0 at MIN_DECOMP_TEMP (-4C).
		} 
		else {
			// Default: Q10 relationship (Schaefer & Jafarov, 2016)
			double q10_freeze = 200.5; // i.e. average of 164 and 237 based on incubation of frozen soil samples (Mikan et al., 2002)
			temp_mod = decomp_at_freezing_point * pow(q10_freeze, temp_soil / 10.0);
		}
	}

	return temp_mod;
}


/// Water Modifier for Decomposition
/** Calculate decomposition moisture modifier (in range 0-1)
 *  Friend et al 1997, Eqn 53 (Parton et al 1993, Fig 2)
 *
 *  \param wfps  Water-filled pore space
 */
double moisture_modifier(double wfps) {

	double moist_mod = wfps < 60.0 ? exp((wfps - 60.0) * (wfps - 60.0) / -800.0) : 0.000371 * wfps * wfps - 0.0748 * wfps + 4.13;
	
	return moist_mod;
}

/// Calculates CENTURY instantaneous decay rates
/** Calculates CENTURY instantaneous decay rates given soil temperature,
 *  water content of upper soil layer
 */
void decayrates_century(Soil& soil, double temp_soil, double wcont_soil, bool tillage) {

	// Maximum exponential decay constants for each SOM pool (daily basis)
	// (Parton et al 2010, Figure 2)
	// plus Kirschbaum et al 2001 coarse woody debris decay
	// pools SURFSTRUCT,SOILSTRUCT,SOILMICRO,SURFHUMUS,SURFMICRO,SURFMETA,SURFFWD,SURFCWD,SOILMETA,SLOWSOM,PASSIVESOM
	const double K_MAX[] = {9.5e-3, 1.9e-2, 4.2e-2, 4.8e-4, 2.7e-2, 3.8e-2, 1.1e-2, 2.2e-3, 7.0e-2, 1.7e-3, 1.9e-6};

	// Modifier for effect of soil texture
	// Eqn 5, Parton et al 1993:

	// Modify decomposition below if this is a high-latitude peatland
	bool ispeatland = soil.patch.stand.is_highlatitude_peatland_stand();

	// Modify decomposition below if this is a wetland on mineral soils
	bool ismineralwetland = soil.patch.stand.is_true_wetland_stand();

	const double texture_mod = 1.0 - 0.75 * (soil.soiltype.clay_frac + soil.soiltype.silt_frac);
	const double texture_mod_peat = 1.0 - 0.75 * (soil.soiltype.clay_frac_peat + soil.soiltype.silt_frac_peat); // = 1

	// Calculate decomposition temperature modifier (in range 0-1)
	double temp_mod = temperature_modifier(temp_soil);

	// Calculate decomposition moisture modifier (in range 0-1)
	// Water Filled Pore Spaces (wfps) % water holding capacity at wilting point (wp) and saturation capacity (wsats)
	// is calculated with the help of Cosby et al 1984;
	// use Gerten equivalents here, but wfps COULD be made depth equivalent
	const double wfps = soil.wfps(0)*100.0;
	double moist_mod = moisture_modifier(wfps);
	double moist_mod_inundated_mineral = moisture_modifier(100); // 100% WFPS for wetlands on mineral soils, 0.36 approx

	// Combined moisture and temperature modifier
	
	// simple overrides for peatlands and mineral wetlands
	double moist_mod_saturated = 1.0; // no effect unless this is peatland

	if (ispeatland) {

		moist_mod = RMOIST;

		double cmass_total = 0.0;
		for (int p = 0; p < NSOMPOOL; p++) {
			cmass_total += soil.sompool[p].cmass;
		}

		double acrotelm_climit = 7.5; // kgC/m2 - Max C content in a 30cm-deep acrotelm - see Wania et al. (2009b)
		if (cmass_total > acrotelm_climit) { // take a weighted average of the aerobic and anaerobic moisture modifiers 
			moist_mod = (acrotelm_climit * RMOIST + (cmass_total-acrotelm_climit) * RMOIST_ANAEROBIC) / cmass_total;
			// moist_mod approaches a value of RMOIST_ANAEROBIC asymptotically as cmass_total increases
		}

		moist_mod_saturated = RMOIST_ANAEROBIC / moist_mod;
	}

	if (ismineralwetland)
		moist_mod = moist_mod_inundated_mineral;

	for (int p = 0; p < NSOMPOOL; p++) {

		// Calculate decay constant
		// (dC_I/dt / C_I; Parton et al 1993, Eqns 2-4)

		double k = K_MAX[p] * temp_mod * moist_mod;

		// Include effect of recalcitrance effect of lignin
		// Parton et al 1993 Eqn 2
		// Kirschbaum et al 2001 changed the exponential term
		// from 3 to 5.

		if (p == SURFSTRUCT || p == SOILSTRUCT || p == SURFFWD || p == SURFCWD) {
			k *= exp(-5.0 * soil.sompool[p].ligcfrac);
		}
		else if (p == SOILMICRO) {
			if (ispeatland)
				k *= texture_mod_peat;
			else
				k *= texture_mod;
		}

		// Increased HR for crops (tillage)
		if (tillage && (p == SURFMICRO || p == SURFHUMUS || p == SOILMICRO || p == SLOWSOM) && !ispeatland) {
			k *= TILLAGE_FACTOR;
		}

		// Reduced decomposition for the passive and slow pools in peatlands, as they are assumed to be in the catotelm
		if (p == PASSIVESOM || p == SLOWSOM) {
			k *= moist_mod_saturated; // ensures that a modifier of RMOIST_ANAEROBIC is used.
		}

		// Calculate fraction of carbon pool remaining after today's decomposition
		soil.sompool[p].fracremain = exp(-k);

		if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
			soil.sompool[p].mfracremain_mean[date.month] += soil.sompool[p].fracremain / date.ndaymonth[date.month];
		}
	}
}

/// Transfers specified fraction (frac) of today's decomposition
/** Transfers specified fraction (frac) of today's decomposition in donor pool type
  *  to receiver pool, transferring fraction respfrac of this to the accumulated CO2
 *  flux respsum (representing total microbial respiration today)
 *  Added P parameters
 */
void transferdecomp(Soil& soil, pooltype donor, pooltype receiver,
	double frac, double respfrac, double& respsum, double& nmin_actual, double& pmin_actual,
	double& nimmob, double& pimmob, double& net_min, double& net_pmin) {

	// decrement in donor carbon pool and nitrogen pools
	double cdec = soil.sompool[donor].cdec * frac;
	double ndec = soil.sompool[donor].ndec * frac;
	double pdec = soil.sompool[donor].pdec * frac;

	// associated nitrogen increment in receiver pool (Friend et al 1997, Eqn 49)
	double ninc = cdec * (1.0 - respfrac) * soil.sompool[receiver].ntoc;

	// associated phosphorus increment in receiver pool
	double pinc = cdec * (1.0 - respfrac) * soil.sompool[receiver].ptoc;

	/*if(std::isnan(pinc))
		pinc = 0.0;*/

	// if increase in receiver nitrogen greater than decrease in donor nitrogen,
	// balance must be immobilisation from mineral nitrogen pool
	// otherwise balance is nitrogen mineralisation
	if (ninc > ndec) {
		nimmob += ninc - ndec;
		net_min += ndec - ninc;
	}
	else {
		nmin_actual += ndec - ninc;
		net_min += ndec - ninc;
	}

	// if increase in receiver phosphorus greater than decrease in donor phosphorus,
	// balance must be immobilisation from mineral phosphorus pool
	// otherwise balance is phosphorus mineralisation
	if (pinc > pdec) {
		pimmob += pinc - pdec;
		net_pmin += pdec - pinc;
	}
	else {
		pmin_actual += pdec - pinc;
		net_pmin += pdec - pinc;
	}

	// "Transfer" carbon and nitrogen to receiver
	soil.sompool[receiver].delta_cmass += cdec * (1.0 - respfrac);
	soil.sompool[receiver].delta_nmass += ninc;
	soil.sompool[receiver].delta_pmass += pinc;

	// Transfer microbial respiration
	respsum += cdec * respfrac;
}

/// Fluxes between the CENTURY pools, and CO2 release to the atmosphere
/** Daily or monthly fluxes between the ten CENTURY pools, and CO2 release to the atmosphere
 *  Parton et al 1993, Fig 1; Comins & McMurtrie 1993, Appendix A
 *
 *  \param ifequilsom Whether the function is called during calculation om SOM pool equilibrium,
 *                    \see equilsom(). During this stage, somfluxes shouldn't calculate decayrates
 *                    itself, and shouldn't produce output like fluxes etc.
 */
void somfluxes(Patch& patch, bool ifequilsom, bool tillage) {

	double respsum ;
	double leachsum_cmass, leachsum_nmass, leachsum_pmass;
	double nmin_actual;	// actual (not net) nitrogen mineralisation
	double pmin_actual;	// actual (not net) phosphorus mineralisation
	double nimmob;		// nitrogen immobilisation
	double pimmob;		// phosphorus immobilisation

	const double EPS = 1.0e-16;

	Soil& soil = patch.soil;

	// mineral nitrogen mass available
	const double nmin_mass = soil.nmass_avail(NH4);// + soil.NO3_mass;
	// mineral phosphorus mass available
	const double pmin_mass = soil.pmass_labile;
	
	if (date.day == 0) {
		soil.anmin = 0.0;
		soil.animmob = 0.0;
		soil.apmin = 0.0;
		soil.apimmob = 0.0;
	}

	///////////////////////////// Balanced dynamics of P labile and P sorbed, also flux into strongly sorbed pool

	//double delta_plabile = soil.pmass_labile_delta;

	//double delta_strongly_sorbed = USORB * soil.pmass_sorbed - USSORB * soil.pmass_strongly_sorbed;

	//pmass_add(soil, delta_plabile);
	
	

	//// Silly protection against negative p values, improve.
	//if (soil.pmass_labile < 0.0) {
	//	patch.fluxes.report_flux(Fluxes::P_SOIL, soil.pmass_labile);
	//	soil.pmass_labile = 0.0;
	//}

	//if (soil.pmass_sorbed < 0.0) {
	//	patch.fluxes.report_flux(Fluxes::P_SOIL, soil.pmass_sorbed);
	//	soil.pmass_sorbed = 0.0;
	//}

	

	//soil.pmass_labile_delta = 0.0;

	////////////////////////////////////////////////////////////////

	// Warning if soil available nitrogen is negative (if happens once or so no problem, but if it propagates through time then it is)

	if (ifnlim) {
		//assert(soil.NH4_mass > -EPS);
		if(soil.NH4_mass < -EPS)
			dummy++;
		if (soil.NH4_mass < -EPS && dummy > 2)
			assert(soil.NH4_mass > -EPS);
	}

	// Warning if soil available phosphorus is negative (if happens once or so no problem, but if it propagates through time then it is)
	if (ifplim) {
		assert(soil.pmass_labile > -EPS);
	}

	// Set N:C ratios for humus, soil microbial, passive and slow pool based on estimated mineral nitrogen pool
	// (Parton et al 1993, Fig 4)

	// ForCent (Parton 2010) values
	setntoc(soil, nmin_mass, SLOWSOM, 30.0, 15.0, 0.0, NMASS_SAT);

	setntoc(soil, nmin_mass, SOILMICRO, 15.0, 6.0, 0.0, NMASS_SAT);

	setntoc(soil, nmin_mass, SURFHUMUS, 30.0, 15.0, 0.0, NMASS_SAT);

	// Set P:C ratios for the slow, passive and soil microbial pools from the
	// labile P pool. The three (ctop_max, ctop_min) pairs are Parton, Stewart
	// and Cole (1988) Fig. 3, p. 115, line for line. The saturating value of
	// the driver is that figure's axis maximum CONVERTED, because pmin_mass is
	// soil.pmass_labile, the wider Hedley-labile pool, and not the
	// resin-extractable orthophosphate the figure was drawn against. See the
	// PMASS_SAT comment at the top of this file.

	setptoc(soil, pmin_mass, SLOWSOM, 200.0, 90.0, 0.0, PMASS_SAT);

	setptoc(soil, pmin_mass, PASSIVESOM, 200.0, 20.0, 0.0, PMASS_SAT);

	// Checked against the source: 80 and 30 are Fig. 3's ACTIVE soil P line,
	// and this fork's soil microbial pool is CENTURY's active pool. The
	// paper's conclusions (p. 128) give the active range as 20 to 80 where its
	// model description (p. 115) and Fig. 3 give 30 to 80; the code follows
	// the figure.
	setptoc(soil, pmin_mass, SOILMICRO, 80.0, 30.0, 0.0, PMASS_SAT);

	// DECLARED DIVERGENCE FROM MAINLINE: surfhumus_ptoc_ramp, owner
	// WORLD-SHCP, registered in biosphere/config/somdynam.yaml. The vendored
	// CNP fork leaves this line commented out,
	//     //setptoc(soil, pmin_mass, SURFHUMUS, 200.0, 90.0, 0.0, PMASS_SAT);
	// so the surface humus pool alone among the ramped pools has no phosphorus
	// ramp, and holds the C:P soil.cpp initialises it to for the whole of a
	// run. It runs here.
	//
	// The objection this answers is that Parton, Stewart and Cole (1988) has no
	// humus pool, so 200 and 90 look like SLOWSOM's pair imported without a
	// source. They are not imported. The PAIR is Fig. 3's slow line, and the
	// IDENTIFICATION of surface humus with the slow pool is the model's own,
	// stated twice in this fork: setntoc() above gives SURFHUMUS exactly
	// SLOWSOM's (30, 15) off exactly SLOWSOM's driver, the mineral nitrogen
	// pool, and soil.cpp's own alternative P initialisation, commented out
	// beside the live one, sets SURFHUMUS and SLOWSOM to the same 1/90, which
	// is that line's phosphorus-rich end. So the nitrogen side already treats
	// this pool as the slow pool at the surface, and this is that treatment
	// carried to the other element rather than a new claim.
	//
	// The alternative was a measured surface humus C:P, and it is the wrong
	// KIND of quantity. What setptoc sets is the P:C at which a pool RECEIVES
	// carbon -- a stoichiometric target that moves with labile P, by the
	// phosphatase mechanism of McGill and Cole (1981) the paper builds on --
	// where a measured forest-floor C:P is an emergent bulk ratio. A fixed
	// number is the wrong shape for it however well sourced, which is what
	// WORLD-SHCP found the fork's 1/150 to be.
	//
	// The invariant the P immobilisation branch below keeps: A POOL WHOSE P:C
	// IS FLEXED DOWN THERE MUST BE ONE THIS BLOCK RE-DERIVES. SURFHUMUS was
	// flexed and not re-derived, so its P:C ratcheted down without bound from
	// the value soil.cpp initialises it to, and the surface humus pool
	// asymptotically received carbon carrying no phosphorus. This line is what
	// re-derives it, so it is back in that list and the two sets are the same
	// set again, which is what WORLD-16PB asked for. PASSIVESOM is the harmless
	// other direction: re-derived here, never flexed. Evidenced in
	// biosphere/notes/phosphorus-cycle-parameterisation.md.
	setptoc(soil, pmin_mass, SURFHUMUS, 200.0, 90.0, 0.0, PMASS_SAT);


	if (!ifequilsom) {

		// Calculate potential fraction remaining following decay today for all pools
		// (assumes no nitrogen limitation)
		decayrates_century(soil, soil.get_soil_temp_25(), soil.get_soil_water_upper(), tillage);

	}

	// Calculate decomposition in all pools assuming these decay rates

	// Save delta carbon and nitrogen mass

	bool net_mineralization = false;
	bool net_pmineralization = false;
	int times = 0;
	int ptimes = 0;
	double decay_reduction_np[NSOMPOOL] = { 0.0 };
	double decay_reduction_n[NSOMPOOL] = { 0.0 };
	double decay_reduction_p[NSOMPOOL] = { 0.0 };
	double init_negative_nmass, init_ntoc_reduction;
	double ntoc_reduction = 0.8;
	double init_negative_pmass, init_ptoc_reduction;
	double ptoc_reduction = 0.8;

	// If necessary, the decay rates in the pools will be reduced in groups, one group
	// is reduced after each iteration in the loop below. The groups are defined by
	// how the pools feed into each other.
	SomPoolSelection reduction_groups[4];
	reduction_groups[0].set(SURFSTRUCT).set(SURFMETA).set(SURFFWD).set(SURFCWD).set(SOILSTRUCT).set(SOILMETA);
	reduction_groups[1].set(SURFMICRO);
	reduction_groups[2].set(SURFHUMUS);
	reduction_groups[3].set(SOILMICRO).set(SLOWSOM).set(PASSIVESOM);

	// If mineralization together with soil available nitrogen is negative then decay rates are decreased
	// The SOM system have five try to get a positive result, after that all pools decay rate has been
	// affected by nitrogen limitation
	while ((!net_mineralization || !net_pmineralization) && (times < 5 || ptimes < 5)) {
	//while ((!net_mineralization && !net_pmineralization) && (times < 5 && ptimes < 5)) {

		respsum = 0.0;
		nmin_actual = 0.0;
		pmin_actual = 0.0;
		nimmob = 0.0;
		pimmob = 0.0;
		leachsum_cmass = 0.0;
		leachsum_nmass = 0.0;
		leachsum_pmass = 0.0;

		// Calculate decomposition in all pools assuming these decay rates
		for (int p = 0; p < NSOMPOOL; p++) {
			
			if(ifplim)
				decay_reduction_np[p] = max(decay_reduction_n[p], decay_reduction_p[p]);
			else
				//Mateus: Choose larger decay reduction between N and P
				decay_reduction_np[p] = decay_reduction_n[p];

			soil.sompool[p].cdec = soil.sompool[p].cmass * (1.0 - soil.sompool[p].fracremain) * (1.0 - decay_reduction_np[p]);
			soil.sompool[p].ndec = soil.sompool[p].nmass * (1.0 - soil.sompool[p].fracremain) * (1.0 - decay_reduction_np[p]);
			soil.sompool[p].pdec = soil.sompool[p].pmass * (1.0 - soil.sompool[p].fracremain) * (1.0 - decay_reduction_np[p]);

			soil.sompool[p].delta_cmass = 0.0;
			soil.sompool[p].delta_nmass = 0.0;
			soil.sompool[p].delta_pmass = 0.0;
			soil.sompool[p].delta_cmass -= soil.sompool[p].cdec;
			soil.sompool[p].delta_nmass -= soil.sompool[p].ndec;
			soil.sompool[p].delta_pmass -= soil.sompool[p].pdec;
		}

		double net_min[NSOMPOOL] = {0};
		double net_pmin[NSOMPOOL] = {0};

		// Partition potential decomposition among receiver pools

		// Donor pool SURFACE STRUCTURAL

		transferdecomp(soil, SURFSTRUCT, SURFMICRO, 1.0 - soil.sompool[SURFSTRUCT].ligcfrac,
			0.6, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFSTRUCT], net_pmin[SURFSTRUCT]);

		transferdecomp(soil, SURFSTRUCT, SURFHUMUS, soil.sompool[SURFSTRUCT].ligcfrac, 0.3,
			respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFSTRUCT], net_pmin[SURFSTRUCT]);

		// Donor pool SURFACE METABOLIC

		transferdecomp(soil, SURFMETA, SURFMICRO, 1.0, 0.6, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFMETA], net_pmin[SURFMETA]);

		// Donor pool SOIL STRUCTURAL

		transferdecomp(soil, SOILSTRUCT, SOILMICRO, 1.0 - soil.sompool[SOILSTRUCT].ligcfrac,
			0.55, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SOILSTRUCT], net_pmin[SOILSTRUCT]);

		transferdecomp(soil, SOILSTRUCT, SLOWSOM, soil.sompool[SOILSTRUCT].ligcfrac, 0.3,
			respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SOILSTRUCT], net_pmin[SOILSTRUCT]);

		// Donor pool SOIL METABOLIC

		transferdecomp(soil, SOILMETA, SOILMICRO, 1.0, 0.55, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SOILMETA], net_pmin[SOILMETA]);

		// Donor pool SURFACE FINE WOODY DEBRIS

		transferdecomp(soil, SURFFWD, SURFMICRO, 1.0 - soil.sompool[SURFFWD].ligcfrac,
			0.76, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFFWD], net_pmin[SURFFWD]);

		transferdecomp(soil, SURFFWD, SURFHUMUS, soil.sompool[SURFFWD].ligcfrac, 0.4,
			respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFFWD], net_pmin[SURFFWD]);

		// Donor pool SURFACE COARSE WOODY DEBRIS

		transferdecomp(soil, SURFCWD, SURFMICRO, 1.0 - soil.sompool[SURFCWD].ligcfrac,
			0.9, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFCWD], net_pmin[SURFCWD]);

		transferdecomp(soil, SURFCWD, SURFHUMUS, soil.sompool[SURFCWD].ligcfrac, 0.5,
			respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFCWD], net_pmin[SURFCWD]);

		// Donor pool SURFACE MICROBE

		transferdecomp(soil, SURFMICRO, SURFHUMUS, 1.0, 0.6, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFMICRO], net_pmin[SURFMICRO]);

		// Donor pool SURFACE HUMUS

		transferdecomp(soil, SURFHUMUS, SLOWSOM, 1.0, 0.6, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SURFHUMUS], net_pmin[SURFHUMUS]);

		// Donor pool SLOW SOM

		// First work out partitioning coefficients (Fig 1, Parton et al 1993)
		double csp = max(0.0, 0.003 - 0.009 * soil.get_clayfrac());
		double respfrac = 0.55;
		double csa = 1.0 - csp - respfrac;

		transferdecomp(soil, SLOWSOM, SOILMICRO, csa, 0.0, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SLOWSOM], net_pmin[SLOWSOM]);

		transferdecomp(soil, SLOWSOM, PASSIVESOM, csp, 0.0, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SLOWSOM], net_pmin[SLOWSOM]);

		// Account for respiration flux
		// Nitrogen associated with this respiration is mineralised (Parton et al 1993, p 791)
		respsum += respfrac * soil.sompool[SLOWSOM].cdec;

		if (!negligible(soil.sompool[SLOWSOM].cmass)) {
			nmin_actual += respfrac * soil.sompool[SLOWSOM].cdec * soil.sompool[SLOWSOM].nmass / soil.sompool[SLOWSOM].cmass;
			pmin_actual += respfrac * soil.sompool[SLOWSOM].cdec * soil.sompool[SLOWSOM].pmass / soil.sompool[SLOWSOM].cmass;
		}

		// Donor pool SOIL MICROBE

		// Fraction lost to  microbial respiration (F_t, Parton et al 1993 Eqn 7)
		respfrac = max(0.0, 0.85 - 0.68 * (soil.get_clayfrac() + soil.get_siltfrac()));

		// Fraction entering passive SOM pool (Parton et al 1993, Eqn 9)
		double cap = 0.003 + 0.032 * soil.get_clayfrac();

		transferdecomp(soil, SOILMICRO, PASSIVESOM, cap, 0.0, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SOILMICRO], net_pmin[SOILMICRO]);

		// Fraction entering slow SOM pool.
		//
		// A bare remainder with no sign guard, and that is the source's own
		// form: Parton et al. (1993) eqn 10, p. 788, is
		// C_AS = (1 - C_AL - C_AP - F_t). Under that paper's OWN eqn 8,
		// C_AL = (H2O_30/18)*(0.01 + 0.04*sand), the remainder is positive over
		// the whole texture simplex, bottoming out at 0.097 on pure sand. It is
		// the CENTURY 5 leaching update below -- 0.03 + 0.12*sand saturating at
		// 1.9 cm H2O per month rather than 18 -- that lets it reach -0.003.
		//
		// Reversing needs a clay-plus-silt fraction below 0.0039 AND percolation
		// at 0.98 of the leaching saturation point, both at once. No input path
		// in this model reaches that: this world's soil map tops out at 0.663
		// sand with a worst remainder of 0.263, and soilinput.h's coarse LPJ
		// soil code, the sandiest texture the model can be handed at 0.90, still
		// leaves 0.0754. So no clamp is applied, because a clamp would have to
		// decide which of the four shares absorbs an excess that cannot occur,
		// and max(0, csp) would lose mass. WORLD-T67J, bounded and closed.
		// biosphere/config/mineral_reactivity.yaml registers all four bounds and
		// biosphere/scripts/mineral_reactivity_gate.py fails on a coefficient
		// change that moves any of them.
		csp = 1.0 - respfrac - soil.orgleachfrac - cap;

		transferdecomp(soil, SOILMICRO, SLOWSOM, csp, 0.0, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[SOILMICRO], net_pmin[SOILMICRO]);

		// Account for respiration flux
		// nitrogen associated with this respiration is mineralised (Parton et al 1993, p 791)
		respsum += respfrac * soil.sompool[SOILMICRO].cdec;

		// Account for organic carbon leaching loss
		leachsum_cmass = soil.orgleachfrac * soil.sompool[SOILMICRO].cdec;

		if (!negligible(soil.sompool[SOILMICRO].cmass)) {
			nmin_actual += respfrac * soil.sompool[SOILMICRO].cdec * soil.sompool[SOILMICRO].nmass / soil.sompool[SOILMICRO].cmass;
			pmin_actual += respfrac * soil.sompool[SOILMICRO].cdec * soil.sompool[SOILMICRO].pmass / soil.sompool[SOILMICRO].cmass;

			// Account for organic nitrogen leaching loss
			leachsum_nmass = soil.orgleachfrac * soil.sompool[SOILMICRO].cdec * soil.sompool[SOILMICRO].nmass / soil.sompool[SOILMICRO].cmass;
			// Account for organic phosphorus leaching loss
			leachsum_pmass = soil.orgleachfrac * soil.sompool[SOILMICRO].cdec * soil.sompool[SOILMICRO].pmass / soil.sompool[SOILMICRO].cmass;
			// A very high sand_frac causes high organic leaching in the soil type 1, coarse inceptisoil, causing very high P limitation.
		}

		// Donor pool PASSIVE SOM

		transferdecomp(soil, PASSIVESOM, SOILMICRO, 1.0, 0.55, respsum, nmin_actual, pmin_actual, nimmob, pimmob, net_min[PASSIVESOM], net_pmin[PASSIVESOM]);

		// Total net mineralization
		double tot_net_min = nmin_actual - nimmob;

		// Estimate daily soil mineral nitrogen pool after decomposition
		// (negative value = immobilisation)
		if ((tot_net_min + nmin_mass + EPS >= 0.0) || !ifnlim) {

			net_mineralization = true;
		}
		else if (!ifnlim) {

			// Unreachable: the condition above already takes !ifnlim. Upstream's
			// form, kept because the phosphorus arm below is a copy of it and
			// the two are only legible together.
			// Not minding immobilisation higher than nmass_avail during free nitrogen years
			if (date.year > freenyears) {

				// Immobilization larger than soil available nitrogen -> reduce targeted N concentration in SOM pool with flexible N:C ratios
				if (times == 0) {
					// initial reduction
					init_negative_nmass = tot_net_min + nmin_mass;
					init_ntoc_reduction = ntoc_reduction;
				}
				else {
					// trying to match needed N:C reduction
					ntoc_reduction = min(init_ntoc_reduction, pow(init_ntoc_reduction, 1.0 / (1.0 - (tot_net_min + nmin_mass) / init_negative_nmass) + 1.0));
				}

				soil.sompool[SLOWSOM].ntoc *= ntoc_reduction;
				soil.sompool[SOILMICRO].ntoc *= ntoc_reduction;
				soil.sompool[SURFHUMUS].ntoc *= ntoc_reduction;

				net_mineralization = false;
			}
			else {
				net_mineralization = true;
			}
		}
		else {

			// Immobilization larger than soil available nitrogen -> reduce decay rates
			if (times < 4) {
				reduce_decay_rates(decay_reduction_n, net_min, reduction_groups[times], tot_net_min + nmin_mass);
			}
			net_mineralization = false;
		}
		times++;
	
		
		// Total net P mineralization
		double tot_net_pmin = pmin_actual - pimmob;

		//subtract pflux up to here to calculate total pmineralization

		// Phosphorus: reduce decay rates when immobilisation exceeds what is
		// available. (negative value = immobilisation)
		//
		// || !ifplim is this fork's, and it restores the nitrogen twin above.
		// This whole block is that branch with n substituted for p -- every
		// comment inside it still says nitrogen -- and the copy dropped the
		// short-circuit from the first condition. Without it, ifplim 0 did not
		// mean "phosphorus does not limit": it meant the P:C flexing arm ran
		// and the decay-rate arm did not, which is neither setting. With it,
		// ifplim 0 means what ifnlim 0 already means. ifplim 1 is untouched:
		// !ifplim was already false there, so the arm below has never been
		// reachable under phosphorus limitation and still is not. WORLD-16PB.
		if ((tot_net_pmin + pmin_mass + EPS >= 0.0) || !ifplim) {

			net_pmineralization = true;
		}
		else if (!ifplim) {

			// Unreachable, exactly as the nitrogen arm above it is, and kept
			// for the same reason: it is upstream's form and the two must be
			// read together. SURFHUMUS is back in the list below, because
			// WORLD-SHCP gave it the phosphorus ramp that re-derives it. The
			// set setptoc() covers and the set flexed here are one set again,
			// which is the invariant that makes flexing a within-day iteration
			// rather than a one-way ratchet. WORLD-16PB.
			if (date.year > freenyears) {

				if (ptimes == 0) {
					// initial reduction
					init_negative_pmass = tot_net_pmin + pmin_mass;
					init_ptoc_reduction = ptoc_reduction;
				}
				else {
					// trying to match needed P:C reduction
					ptoc_reduction = min(init_ptoc_reduction, pow(init_ptoc_reduction, 1.0 / (1.0 - (tot_net_pmin + pmin_mass) / init_negative_pmass) + 1.0));
				}

				soil.sompool[SLOWSOM].ptoc *= ptoc_reduction;
				soil.sompool[SOILMICRO].ptoc *= ptoc_reduction;
				soil.sompool[SURFHUMUS].ptoc *= ptoc_reduction;

				net_pmineralization = false;
			}
			else {
				net_pmineralization = true;
			}
		}
		else {

			// Immobilisation larger than soil available phosphorus -> reduce decay rates
			if (ptimes < 4) {
				reduce_decay_rates(decay_reduction_p, net_pmin, reduction_groups[ptimes], tot_net_pmin + pmin_mass);
			}
			net_pmineralization = false;
		}

		ptimes++;
 
	}

	// Most of organic leaching is retained in the ecosystem (82%, Wilcke), and mineralized (Parton 1988)
	double leachsum_pmass_retained = leachsum_pmass * 1.0;

	double leachsum_pmass_lost = leachsum_pmass * 0.0;


	// Sample the organic-leaching OPERATOR after the final nutrient-limited
	// decomposition solve. Its monthly coefficient is leached microbial C over
	// microbial C decay, not the arithmetic mean of daily hydrologic fractions:
	// a dry day with no decomposition must not carry the same weight as the day
	// on which material moved.
	if (!ifequilsom && date.year >= soil.solvesomcent_beginyr &&
		date.year <= soil.solvesomcent_endyr) {
		soil.morgleach_cmass[date.month] += leachsum_cmass;
		soil.msoilmicro_cdec[date.month] += soil.sompool[SOILMICRO].cdec;
	}

	// Update pool sizes

	for (int p = 0; p < NSOMPOOL; p++) {
		soil.sompool[p].cmass += soil.sompool[p].delta_cmass;
		soil.sompool[p].nmass += soil.sompool[p].delta_nmass;
		soil.sompool[p].pmass += soil.sompool[p].delta_pmass;
	}

	if (!ifequilsom) {

		// Updated soil fluxes. In wetlands, some portion of cflux can be emitted as CH4, so save cflux as dcflux_soil until Soil::methane() is called.  
		if (patch.stand.landcover != PEATLAND)
			patch.fluxes.report_flux(Fluxes::SOILC, respsum);
		else 
			soil.dcflux_soil=respsum; 

		// Sum annual organic nitrogen leaching

		soil.aorgNleach += leachsum_nmass;

		// Sum annual organic phosphorus leaching

		//soil.aorgPleach += leachsum_pmass;
		// Leached organic phosphorus is assumed to be mineralized (Parton 1988, pg. 117)
		soil.aorgPleach += leachsum_pmass_lost;

		// Sum annual organic carbon leaching

		soil.aorgCleach += leachsum_cmass;

		// Sum annuals
		soil.anmin += nmin_actual;
		soil.animmob += nimmob;
		soil.apmin += pmin_actual;
		soil.apimmob += pimmob;
	}

	// Fraction of microbial resp. is assumed to produce labile carbon
	soil.labile_carbon = respsum * frac_labile_carbon;

	// Adding mineral nitrogen to soil available pool
	double nmin_inc = nmin_actual - nimmob; 
		
	//double pmin_inc = pmin_actual - pimmob;
	// Leached organic phosphorus is assumed to be mineralized (Parton 1988, pg. 117)
	double pmin_inc = pmin_actual - pimmob + leachsum_pmass_retained;

	//soil.pmass_labile = max(0.0, soil.pmass_labile + pmin_inc);
	//soil.pmass_labile_delta += pmin_inc;

	pmass_add(soil, pmin_inc);

	///////////////////////////// Balanced dynamics of P labile and P sorbed, also flux into strongly sorbed pool

	//delta_plabile = soil.pmass_labile_delta;

	//pmass_add(soil, delta_plabile);

	//// Silly protection against negative p values, improve.
	//if (soil.pmass_labile < 0.0) {
	//	patch.fluxes.report_flux(Fluxes::P_SOIL, soil.pmass_labile);
	//	soil.pmass_labile = 0.0;
	//}

	//if (soil.pmass_sorbed < 0.0) {
	//	patch.fluxes.report_flux(Fluxes::P_SOIL, soil.pmass_sorbed);
	//	soil.pmass_sorbed = 0.0;
	//}

	// Wang et al. (2010) Eq. D10, and the RECEIVING SIDE OF IT, which the fork
	// left out: pmass_strongly_sorbed was declared, initialised, serialised and
	// output, and read here, but never written anywhere in the tree. With it
	// pinned at zero the back term vanished, so this was not a transfer between
	// two pools but a first-order drain of the sorbed pool that could never shut
	// off, at USORB per Earth year against a weathering input two orders of
	// magnitude smaller. It closed the phosphorus balance only because pcont()
	// excluded the destination pool and this line booked the difference as an
	// ecosystem loss.
	//
	// Assigning the pool makes both terms real. USORB equals USSORB, so the
	// strongly sorbed pool fills to the size of the sorbed pool and the net flux
	// goes to zero, which is the behaviour the cited equation describes. The
	// flux is an internal transfer and is no longer reported as a soil P loss;
	// Patch::pcont() now counts the pool instead, and double-booking it would
	// break the balance the two together close.
	//
	// This has no terminal sink beyond it. Occlusion is absent from this model
	// by decision, not by oversight: see
	// biosphere/notes/phosphorus-cycle-parameterisation.md.
	double delta_strongly_sorbed = USORB * soil.pmass_sorbed - USSORB * soil.pmass_strongly_sorbed;

	pmass_add(soil, -delta_strongly_sorbed);

	soil.pmass_strongly_sorbed += delta_strongly_sorbed;

	//soil.pmass_labile_delta = 0.0;

	////////////////////////////////////////////////////////////////

	// Estimate of N flux from soil (simple CLM-CN approach)
	if(!ifntransform) {
		double nflux = nmin_inc > 0.0 ? (nmin_inc) * 0.01 : 0.0;
		//soil.NH4_mass -= nflux;
		nmin_inc -= nflux;

	if (!ifequilsom) {
		patch.fluxes.report_flux(Fluxes::NH3_SOIL, nflux);
	}
	}
	soil.nmass_inc(nmin_inc,NH4);

	// If no nitrogen limitation or during free nitrogen years set soil
	// available nitrogen to its saturation level.
	if (!ifnlim || date.year <= freenyears) {
		if(ifntransform) {
			soil.NH4_mass = NMASS_SAT / 2.0;
			soil.NO3_mass = NMASS_SAT / 2.0;
		} 
		else {
			soil.NH4_mass = NMASS_SAT;
		}
	}

	// If no phosphorus limitation or during free phosphorus years set soil
	// available phosphorus to its saturation level.
	//
	// PMASS_SAT is the right constant for this second job as well as for the
	// ramp: "saturated" here means exactly "at the value where setptoc() stops
	// responding", so the two uses are one number and must not be split. That
	// is why WORLD-Z01O's conversion was applied to PMASS_SAT itself and this
	// pin followed. What still follows from it is that the labile P this
	// configuration reports is a CONSTANT meaning "not limiting", not a
	// simulated stock, and comparing it against a measured labile P or against
	// this fork's own P-limited run is meaningless in either direction, even
	// though the conversion has moved it to the same order as the fork's own
	// simulated pool. commonoutput.cpp writes it into PO4_mass and availp all
	// the same.
	if (!ifplim || date.year <= freenyears) {
		soil.pmass_labile = PMASS_SAT;
		//soil.pmass_labile = 0.0;
		//Value of sorbed phosphorus pool, based on labile p and soil parameters [KgP/m-2].
		//Equilibrated instantanously based on Wang 2007, 2010
		//soil.pmass_sorbed = (PMASS_SAT * soil.soiltype.spmax) / (soil.soiltype.kplab + PMASS_SAT);
		//soil.pmass_sorbed = (PMASS_SAT * soil.soiltype.spmax * soil.soiltype.kplab) / pow(soil.soiltype.kplab + PMASS_SAT, 2.0);
		soil.pmass_sorbed = soil.soiltype.spmax;
		//soil.pmass_sorbed = 0.0;
	}
}

/// Litter lignin to N ratio (for leaf and root litter)
/** Specific lignin fractions for leaf and root
 *  are specified in transfer_litter()
 *  Check if needed function for P too
 */
double lignin_to_n_ratio(double cmass_litter, double nmass_litter, double LIGCFRAC, double cton_avr) {

	if (!negligible(nmass_litter) && ifnlim) {
		return max(0.0, LIGCFRAC * cmass_litter / nmass_litter);
	}
	else {
		return max(0.0, LIGCFRAC * cton_avr / (1.0 - nrelocfrac));
	}
}

/// Metabolic litter fraction (for leaf and root litter)
/** Fm, Parton et al 1993, Eqn 1:
 *  NB: incorrect/out-of-date intercept and slope given in Eqn 1; values used in
 *  code of CENTURY 4.0 used instead (also correct in Parton et al. 1993, figure 1)
 *
 *  Check if needed function for P too
 * \param lton  Litter lignin:N ratio
 */
double metabolic_litter_fraction(double lton) {
	return max(0.0, 0.85 - lton * 0.013);
}

/// Transfers litter from growth, mortality and fire
/** Called daily to transfer last year's litter from vegetation
 *  (turnover, mortality and fire) to soil litter pools.
 *  Alternatively, with daily allocation and harvest/turnover,
 *  the litter produced a certain day.
 */
void transfer_litter(Patch& patch) {

	Soil& soil = patch.soil;

	double lat = patch.get_climate().lat;

	double EPS = -1.0e-16;

	double ligcmass_new, ligcmass_old;

	// Fire (GlobFIRM)
	double litterme[NSOMPOOL]  = {0.};
	double fireresist[NSOMPOOL]= {0.};
	if ( firemodel == GLOBFIRM ) {

		litterme[SURFSTRUCT]   = soil.sompool[SURFSTRUCT].cmass * soil.sompool[SURFSTRUCT].litterme;
		litterme[SURFMETA]     = soil.sompool[SURFMETA].cmass   * soil.sompool[SURFMETA].litterme;
		litterme[SURFFWD]      = soil.sompool[SURFFWD].cmass    * soil.sompool[SURFFWD].litterme;
		litterme[SURFCWD]      = soil.sompool[SURFCWD].cmass    * soil.sompool[SURFCWD].litterme;

		fireresist[SURFSTRUCT] = soil.sompool[SURFSTRUCT].cmass * soil.sompool[SURFSTRUCT].fireresist;
		fireresist[SURFMETA]   = soil.sompool[SURFMETA].cmass   * soil.sompool[SURFMETA].fireresist;
		fireresist[SURFFWD]    = soil.sompool[SURFFWD].cmass    * soil.sompool[SURFFWD].fireresist;
		fireresist[SURFCWD]    = soil.sompool[SURFCWD].cmass    * soil.sompool[SURFCWD].fireresist;
	}

	double leaf_litter = 0.0;
	double root_litter = 0.0;
	double wood_litter = 0.0;

	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& pft=patch.pft.getobj();

		// For stands with yearly growth, drop leaf and root litter on first month of the year for northern hemisphere
		// and first month of the second half of the year for southern hemisphere for summergreen trees. For evergreens
		// as a fraction every day and for raingreens on the drierst month from last year. 
		// For stands with daily growth, harvest and/or turnover, do this when patch.is_litter_day is true.

		double cmass_litter_leaf, nmass_litter_leaf, pmass_litter_leaf;
		double cmass_litter_root, nmass_litter_root, pmass_litter_root;

		double frac_lr = 0.0;
		//  Is a litter day, drop all leaf litter (crop)
		if (patch.is_litter_day) {
			frac_lr = 1.0;
		}
		// For summergreens drop leaf litter over all days during Jan in NH and July SH
		// and raingreens on the month with lowest phen
		else if ((pft.pft.phenology == SUMMERGREEN && ((lat >= 0.0 && date.month == 0) || (lat < 0.0 && date.month == 6))) ||
			(pft.pft.phenology == RAINGREEN && date.month == pft.driest_mth)) {
			frac_lr = 1.0 / (date.ndaymonth[date.month] - date.dayofmonth);
		}
		// Evergreens drops leaf litter every day
		else if (pft.pft.phenology == EVERGREEN || pft.pft.phenology == ANY) {
			frac_lr = 1.0 / (date.year_length() - date.day);
		}

		cmass_litter_leaf = pft.cmass_litter_leaf * frac_lr;
		nmass_litter_leaf = pft.nmass_litter_leaf * frac_lr;
		pmass_litter_leaf = pft.pmass_litter_leaf * frac_lr;
		cmass_litter_root = pft.cmass_litter_root * frac_lr;
		nmass_litter_root = pft.nmass_litter_root * frac_lr;
		pmass_litter_root = pft.pmass_litter_root * frac_lr;

		pft.cmass_litter_leaf -= cmass_litter_leaf;
		pft.nmass_litter_leaf -= nmass_litter_leaf;
		pft.pmass_litter_leaf -= pmass_litter_leaf;
		pft.cmass_litter_root -= cmass_litter_root;
		pft.nmass_litter_root -= nmass_litter_root;
		pft.pmass_litter_root -= pmass_litter_root;

		// LEAF

		//if (!negligible(cmass_litter_leaf) || !negligible(nmass_litter_leaf) || !negligible(pmass_litter_leaf)) {
		if (!negligible(cmass_litter_leaf) || !negligible(nmass_litter_leaf)) {

			// Calculate inputs to surface structural and metabolic litter

			// Leaf litter lignin:N ratio
			double leaf_lton = lignin_to_n_ratio(cmass_litter_leaf, nmass_litter_leaf, LIGCFRAC_LEAF, pft.pft.cton_leaf_avr);

			// Metabolic litter fraction for leaf litter
			double fm = metabolic_litter_fraction(leaf_lton);

			ligcmass_old = soil.sompool[SURFSTRUCT].cmass * soil.sompool[SURFSTRUCT].ligcfrac;

			// Add to pools
			soil.sompool[SURFSTRUCT].cmass += cmass_litter_leaf * (1.0 - fm);
			soil.sompool[SURFSTRUCT].nmass += nmass_litter_leaf * (1.0 - fm);
			soil.sompool[SURFSTRUCT].pmass += pmass_litter_leaf * (1.0 - fm);
			soil.sompool[SURFMETA].cmass += cmass_litter_leaf * fm;
			soil.sompool[SURFMETA].nmass += nmass_litter_leaf * fm;
			soil.sompool[SURFMETA].pmass += pmass_litter_leaf * fm;

			// Save litter input for equilsom()
			if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
				soil.litterSolveSOM.add_litter(cmass_litter_leaf * (1.0 - fm), nmass_litter_leaf * (1.0 - fm), pmass_litter_leaf * (1.0 - fm), SURFSTRUCT);
				soil.litterSolveSOM.add_litter(cmass_litter_leaf * fm, nmass_litter_leaf * fm, pmass_litter_leaf * fm, SURFMETA);

			}

			// Fire
			if (firemodel == GLOBFIRM) {

				litterme[SURFSTRUCT] += cmass_litter_leaf * (1.0 - fm) * pft.pft.litterme;
				fireresist[SURFSTRUCT] += cmass_litter_leaf * (1.0 - fm) * pft.pft.fireresist;

				litterme[SURFMETA] += cmass_litter_leaf * fm * pft.pft.litterme;
				fireresist[SURFMETA] += cmass_litter_leaf * fm * pft.pft.fireresist;
			}
			// NB: reproduction litter cannot contain nitrogen!!

			ligcmass_new = cmass_litter_leaf * (1.0 - fm) * LIGCFRAC_LEAF;

			if (negligible(soil.sompool[SURFSTRUCT].cmass)) {
				soil.sompool[SURFSTRUCT].ligcfrac = 0.0;
			}
			else {
				soil.sompool[SURFSTRUCT].ligcfrac = (ligcmass_new + ligcmass_old) /
					soil.sompool[SURFSTRUCT].cmass;
			}
		}

		// ROOT

		if (!negligible(cmass_litter_root) || !negligible(nmass_litter_root)) {

			// Calculate inputs to soil structural and metabolic litter

			// Root litter lignin:N ratio
			double root_lton = lignin_to_n_ratio(cmass_litter_root, nmass_litter_root, LIGCFRAC_ROOT, pft.pft.cton_root_avr);

			// Metabolic litter fraction for root litter
			double fm = metabolic_litter_fraction(root_lton);

			ligcmass_new = cmass_litter_root * (1.0 - fm) * LIGCFRAC_ROOT;
			ligcmass_old = soil.sompool[SOILSTRUCT].cmass * soil.sompool[SOILSTRUCT].ligcfrac;

			// Add to pools and update lignin fraction in structural pool
			soil.sompool[SOILSTRUCT].cmass += cmass_litter_root * (1.0 - fm);
			soil.sompool[SOILSTRUCT].nmass += nmass_litter_root * (1.0 - fm);
			soil.sompool[SOILSTRUCT].pmass += pmass_litter_root * (1.0 - fm);
			soil.sompool[SOILMETA].cmass += cmass_litter_root * fm;
			soil.sompool[SOILMETA].nmass += nmass_litter_root * fm;
			soil.sompool[SOILMETA].pmass += pmass_litter_root * fm;

			// Save litter input for equilsom()
			if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
				soil.litterSolveSOM.add_litter(cmass_litter_root * (1.0 - fm), nmass_litter_root * (1.0 - fm), pmass_litter_root * (1.0 - fm), SOILSTRUCT);
				soil.litterSolveSOM.add_litter(cmass_litter_root * fm, nmass_litter_root * fm, pmass_litter_root * fm, SOILMETA);

			}

			if (negligible(soil.sompool[SOILSTRUCT].cmass)) {
				soil.sompool[SOILSTRUCT].ligcfrac = 0.0;
			}
			else {
				soil.sompool[SOILSTRUCT].ligcfrac = (ligcmass_new + ligcmass_old) /
					soil.sompool[SOILSTRUCT].cmass;
			}
		}

		// WOOD

		// Woody debris enters two woody litter pools as described in
		// Kirschbaum and Paul (2002).

		// Woody litter is dropped as a fraction every day
		double frac = 1.0 / (date.year_length() - date.day);
		double cmass_litter_sap = pft.cmass_litter_sap * frac;
		double nmass_litter_sap = pft.nmass_litter_sap * frac;
		double pmass_litter_sap = pft.pmass_litter_sap * frac;
		double cmass_litter_heart = pft.cmass_litter_heart * frac;
		double nmass_litter_heart = pft.nmass_litter_heart * frac;
		double pmass_litter_heart = pft.pmass_litter_heart * frac;

		pft.cmass_litter_sap -= cmass_litter_sap;
		pft.nmass_litter_sap -= nmass_litter_sap;
		pft.pmass_litter_sap -= pmass_litter_sap;
		pft.cmass_litter_heart -= cmass_litter_heart;
		pft.nmass_litter_heart -= nmass_litter_heart;
		pft.pmass_litter_heart -= pmass_litter_heart;

		// SAP WOOD

		if (!negligible(cmass_litter_sap) || !negligible(nmass_litter_sap)) {

			// Fine woody debris

			ligcmass_new = cmass_litter_sap * LIGCFRAC_WOOD;
			ligcmass_old = soil.sompool[SURFFWD].cmass * soil.sompool[SURFFWD].ligcfrac;


			// Add to fine woody pool and update lignin fraction in pool
			soil.sompool[SURFFWD].cmass += cmass_litter_sap;
			soil.sompool[SURFFWD].nmass += nmass_litter_sap;
			soil.sompool[SURFFWD].pmass += pmass_litter_sap;

			// Save litter input for equilsom()
			if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
				soil.litterSolveSOM.add_litter(cmass_litter_sap, nmass_litter_sap, pmass_litter_sap, SURFFWD);
			}

			if (negligible(soil.sompool[SURFFWD].cmass)) {
				soil.sompool[SURFFWD].ligcfrac = 0.0;
			}
			else {
				double ligcfrac = (ligcmass_new + ligcmass_old) /
					soil.sompool[SURFFWD].cmass;
				soil.sompool[SURFFWD].ligcfrac = ligcfrac;
			}

			// Fire
			if (firemodel == GLOBFIRM) {
				litterme[SURFFWD] += cmass_litter_sap * pft.pft.litterme;
				fireresist[SURFFWD] += cmass_litter_sap * pft.pft.fireresist;
			}
		}

		// HEART WOOD

		if (!negligible(cmass_litter_heart) || !negligible(nmass_litter_heart)) {

			// Coarse woody debris

			ligcmass_new = cmass_litter_heart * LIGCFRAC_WOOD;
			ligcmass_old = soil.sompool[SURFCWD].cmass * soil.sompool[SURFCWD].ligcfrac;

			// Add to coarse woody pool and update lignin fraction in pool
			soil.sompool[SURFCWD].cmass += cmass_litter_heart;
			soil.sompool[SURFCWD].nmass += nmass_litter_heart;
			soil.sompool[SURFCWD].pmass += pmass_litter_heart;

			// Save litter input for equilsom()
			if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
				soil.litterSolveSOM.add_litter(cmass_litter_heart, nmass_litter_heart, pmass_litter_heart, SURFCWD);
			}

			if (negligible(soil.sompool[SURFCWD].cmass)) {
				soil.sompool[SURFCWD].ligcfrac = 0.0;
			}
			else {
				double ligcfrac = (ligcmass_new + ligcmass_old) /
					soil.sompool[SURFCWD].cmass;
				soil.sompool[SURFCWD].ligcfrac = ligcfrac;
			}

			// Fire
			if (firemodel == GLOBFIRM) {
				litterme[SURFCWD] += cmass_litter_heart * pft.pft.litterme;
				fireresist[SURFCWD] += cmass_litter_heart * pft.pft.fireresist;
			}

		}

		patch.pft.nextobj();
	}

	// FIRE
	if ( firemodel == GLOBFIRM ) {
		if (soil.sompool[SURFSTRUCT].cmass > 0.0) {
			soil.sompool[SURFSTRUCT].litterme   = litterme[SURFSTRUCT]   / soil.sompool[SURFSTRUCT].cmass;
			soil.sompool[SURFSTRUCT].fireresist = fireresist[SURFSTRUCT] / soil.sompool[SURFSTRUCT].cmass;
		}
		if (soil.sompool[SURFMETA].cmass > 0.0) {
			soil.sompool[SURFMETA].litterme   = litterme[SURFMETA]   / soil.sompool[SURFMETA].cmass;
			soil.sompool[SURFMETA].fireresist = fireresist[SURFMETA] / soil.sompool[SURFMETA].cmass;
		}
		if (soil.sompool[SURFFWD].cmass > 0.0) {
			soil.sompool[SURFFWD].litterme    = litterme[SURFFWD]   / soil.sompool[SURFFWD].cmass;
			soil.sompool[SURFFWD].fireresist  = fireresist[SURFFWD] / soil.sompool[SURFFWD].cmass;
		}
		if (soil.sompool[SURFCWD].cmass > 0.0) {
			soil.sompool[SURFCWD].litterme    = litterme[SURFCWD]   / soil.sompool[SURFCWD].cmass;
			soil.sompool[SURFCWD].fireresist  = fireresist[SURFCWD] / soil.sompool[SURFCWD].cmass;
		}
	}

	// Calculate total litter carbon, nitrogen and phosphorus mass for set N:C and P:C ratio of surface microbial pool
	double litter_cmass = soil.sompool[SURFSTRUCT].cmass + soil.sompool[SURFMETA].cmass +
						  soil.sompool[SURFFWD].cmass + soil.sompool[SURFCWD].cmass;
	double litter_nmass = soil.sompool[SURFSTRUCT].nmass + soil.sompool[SURFMETA].nmass +
						  soil.sompool[SURFFWD].nmass + soil.sompool[SURFCWD].nmass;
	double litter_pmass = soil.sompool[SURFSTRUCT].pmass + soil.sompool[SURFMETA].pmass +
						  soil.sompool[SURFFWD].pmass + soil.sompool[SURFCWD].pmass;

	// Set N:C ratio of surface microbial pool based on N:C ratio of litter from all PFTs
	// Parton et al 1993 Fig 4. Dry mass litter == cmass litter * 2
	//
	// The phosphorus line below has no such figure behind it. Parton, Stewart
	// and Cole (1988) carries no surface microbial pool and no C:P ramp driven
	// by a litter concentration, so both the structure and PCONC_SAT are the
	// nitrogen line's, and the pair 80 and 30 is Fig. 3's ACTIVE SOIL line
	// applied to a surface pool. See the PCONC_SAT comment at the top of this
	// file for the bound on any replacement.
	if (!negligible(litter_cmass)) {
		setntoc(soil, litter_nmass / (litter_cmass * 2.0), SURFMICRO, 20.0, 10.0, 0.0, NCONC_SAT);
		setptoc(soil, litter_pmass / (litter_cmass * 2.0), SURFMICRO, 80.0, 30.0, 0.0, PCONC_SAT);
	}

	// Add litter to solvesom array every month
	if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
		if (date.dayofmonth == 0) {
			soil.solvesom.push_back(soil.litterSolveSOM);
			soil.litterSolveSOM.clear();
		}
	}
}


/// LEACHING
/** Leaching fractions for both organic and mineral leaching
 *  Should be called every day in both daily and monthly mode
 */
void leaching(Soil& soil) {

	double minleachfrac;

	if (!negligible(soil.dperc)) {

		// Leaching from available nitrogen mineral pool
		// in proportion to amount of water drainage
		// Use Gerten equivalents here 
		minleachfrac = soil.dperc / (soil.dperc + soil.soiltype.gawc[0] * soil.get_soil_water_upper());

		// Leaching from decayed organic carbon/nitrogen
		// using Parton et al. eqn. 8; CENTURY 5 parameter update; from equation: C Leached=microbial_C*[OMLECH(1)+OMLECH(2)*sand_fraction]*[1.0f-(OMLECH(3)-water_leaching)/OMLECH(3)],
		// reorganised as: leachfrac = [water_leaching/OMLECH(3)]*[OMLECH(1)+OMLECH(2)*sand_fraction]
		// water_leaching was originally in cm/month
		double OMLECH_1 = 0.03;
		double OMLECH_2 = 0.12;
		double OMLECH_3 = 1.9;	// saturation point in leaching equation (cm H2O/month)
		double cmpermonth_to_mmperday = 10.0 * 12.0 / 365.0;

		soil.orgleachfrac = min(1.0, soil.dperc / (OMLECH_3 * cmpermonth_to_mmperday)) * (OMLECH_1 + OMLECH_2 * soil.soiltype.sand_frac);
		// A very high sand_frac causes high organic leaching in the soil type 1, coarse inceptisoil, causing very high P limitation.
	}
	else {
		minleachfrac = 0.0;
		soil.orgleachfrac = 0.0;
	}

	// Leaching of soil mineral nitrogen
	// Allowed on days with residual nitrogen following vegetation uptake
	// in proportion to amount of water drainage
	const double nmin_avail = soil.nmass_avail(NO3);
	if (nmin_avail > 0.0) {

		double leaching = nmin_avail * minleachfrac;

		soil.nmass_subtract(leaching,NO3);
		soil.aminleach += leaching;
	}

	// Leaching of soil mineral phosphorus
	// Allowed on days with residual phosphorus following vegetation uptake
	// in proportion to amount of water drainage
	const double pmin_avail = soil.pmass_labile;
	if (pmin_avail > 0.0) {

		double leaching_p = pmin_avail * minleachfrac;
		//leaching_p *= 0.18;

		/*soil.pmass_labile -= leaching_p;
		soil.pmass_labile = max(0.0, soil.pmass_labile);*/

		//soil.pmass_labile_delta -= leaching_p;
		pmass_add(soil, -leaching_p);

		soil.aminpleach += leaching_p;

	}

	if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
		soil.mminleach_survival[date.month] *= (1.0 - minleachfrac);
		soil.mminpleach_survival[date.month] *= (1.0 - minleachfrac);
	}
}


/// Nitrogen addition to the soil
/** Daily nitrogen addition to the soil
 *  from deposition and fixation
 */
void soilnadd(Patch& patch) {

	Soil& soil = patch.soil;

	if (!ifntransform) {
		double nflux = 0.01 * (soil.NH4_input + soil.NO3_input);
		soil.NH4_input -= 0.01 * soil.NH4_input;
		soil.NO3_input -= 0.01 * soil.NO3_input;

		patch.fluxes.report_flux(Fluxes::NH3_SOIL, nflux);

	// Nitrogen deposition and fertilization input to the soil (calculated in snow_ninput())
	}

	// Nitrogen fixation
	// If soil available nitrogen is above the value for minimum SOM C:N ratio, then
	// nitrogen fixation is reduced (nitrogen rich soils)
	const double nmin_avail = soil.nmass_avail(NH4);

	if (nmin_avail < NMASS_SAT) {
		// anfix_calc is kgN/m2 per EARTH year (see below), so the absolute-day
		// rate divides by the Earth year. Dividing by the simulation year would
		// deliver an Earth year of fixation every orbit.
		const double daily_nfix = soil.anfix_calc / VESPER_EARTH_YEAR_DAYS;

		if (nmin_avail + daily_nfix < NMASS_SAT) {
			soil.NH4_mass += daily_nfix;
			soil.anfix += daily_nfix;
		}
		else {
			soil.anfix += NMASS_SAT - soil.NH4_mass;
			soil.NH4_mass += NMASS_SAT - nmin_avail;
		}
	}
	// Nitrogen deposition and fertilization input to the soil (calculated in snow_ninput())
	soil.NH4_mass += soil.NH4_input;
	soil.NO3_mass += soil.NO3_input;

	// Calculate nitrogen fixation (Cleveland et al. 1999)
	// by using five year average aaet
	if (date.islastmonth && date.islastday) {

		// Add this simulation year's AET to aaet_5, which keeps the last
		// NYEARAAET simulation years. The window is a count of seasonal cycles
		// and stays as it is; what needs converting is the CONTENT, because each
		// entry is a sum over one orbit and Cleveland's regression is fitted on
		// evapotranspiration per EARTH year.
		patch.aaet_5.add(patch.aaet+patch.aevap+patch.aintercep);

		// Calculate estimated nitrogen fixation (aaet in cm per Earth year, eqn
		// gives kgN/ha per Earth year). BOTH terms needed the conversion and for
		// different reasons: the AET-proportional term because the sum it
		// multiplies spans an orbit rather than an Earth year, and the intercept
		// because it is a rate per Earth year that was being added once per
		// orbit. biosphere/notes/time-base-unit-contract.md.
		const double aaet_mm_per_earth_year =
			patch.aaet_5.mean() / VESPER_EARTH_YEARS_PER_ORBIT;
		soil.anfix_calc = max((nfix_a * aaet_mm_per_earth_year * CM_PER_MM + nfix_b) * HA_PER_M2, 0.0);

		if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
			soil.anfix_mean += soil.anfix;
		}
	}
}

/// Phosphorus addition to the soil
/** Daily phosphorus addition to the soil
*  from deposition and weathering
*/
void soilpadd(Patch& patch) {

	Soil& soil = patch.soil;

	const double pmin_avail = soil.pmass_labile;

	double wfps = soil.wfps(0)*100.0;

	double daily_pwtr, p_temp_effect;

	double R_gas_constant = 8.3144621;
	double ea = patch.get_climate().pwtr_ea;
	double soiltemp = soil.get_soil_temp_25();
	double bi = patch.get_climate().pwtr_bi;
	double pcont = patch.get_climate().pwtr_pcont;
	double shield = patch.get_climate().pwtr_shield;

	if (param["file_pwtr"].str != "") {
		// The gridded route is already a per-absolute-day calculation: runoff is
		// this day's runoff, so the result carries no year in it and needs no
		// conversion. It is inactive on this world, which supplies no file_pwtr.
		p_temp_effect = exp(-ea / R_gas_constant * (1 / (soiltemp + 273.0) - 1 / 284.15));
		daily_pwtr = bi * (pcont / 100.0) * patch.soil.runoff * p_temp_effect * shield / 1000.0;
	}
	else {
		// The texture route, and the one this world takes. Soiltype::pwtr is
		// ABSOLUTE-RATE, declared kgP/m2 per EARTH year, because rock weathering
		// does not know this world's orbit; BIO-5 emits the field against that
		// declaration. Dividing by date.year_length() delivered a whole Earth
		// year of weathered phosphorus every orbit.
		// biosphere/notes/time-base-unit-contract.md.
		//
		// It carries no temperature or moisture dependence. The gridded route
		// takes temperature through p_temp_effect and neither route reads the
		// water-filled pore space computed at the top of this function, so that
		// local stands unused. Supplying a climate dependence here is BIO-5's.
		daily_pwtr = soil.soiltype.pwtr / VESPER_EARTH_YEAR_DAYS;
	}

	//if (pmin_avail + daily_pwtr < PMASS_SAT || date.year <= freenyears) {
		// Phosphorus weathering input
		//soil.pmass_labile += daily_pwtr;
		//soil.pmass_labile_delta += max(0.0, daily_pwtr);
		pmass_add(soil, daily_pwtr);
		soil.apwtr += daily_pwtr;
	//}

	// Phosphorus fertilization and deposition input (calculated in snow_pinput())
	//soil.pmass_labile += soil.pmass_labile_input;
	//soil.pmass_labile_delta += soil.pmass_labile_input;
	pmass_add(soil, soil.pmass_labile_input);
	soil.apdep += soil.pmass_labile_input;

	if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr) {
		// Sample fluxes, not the year-to-date diagnostic accumulators. Adding
		// apwtr/apdep here forms a triangular sum and makes the accelerator's
		// input depend on which day in the year a flux occurred.
		soil.apwtr_mean += daily_pwtr;
		soil.apdep_mean += soil.pmass_labile_input;
	}
	


}


/// Vegetation nitrogen uptake
/** Daily vegetation uptake of mineral nitrogen
 *  Partitioned among individuals according to todays nitrogen demand
 *  and individuals root area
 */
void vegetation_n_uptake(Patch& patch) {

	// Daily nitrogen demand given by:
	//	 For individual:
	//     (1)  ndemand = leafndemand + rootndemand + sapndemand;
    //          where
	//          leafndemand is leaf demand based on vmax
	//			rootndemand is based on optimal leaf C:N ratio
	//          sapndemand is based on optimal leaf C:N ratio
	//
	// Actual nitrogen uptake for each day and individual given by:
	//     (2)  nuptake = ndemand * fnuptake
	//     where fnuptake is individual uptake capacity calculated in fnuptake in canexch.cpp

	double nuptake_day;

	Vegetation& vegetation=patch.vegetation;
	Soil& soil = patch.soil;

	//const double orignmass = soil.NH4_mass + soil.NO3_mass;
	const double orignmass = soil.nmass_avail();
	double ammonium_frac = orignmass ? soil.NH4_mass / orignmass : 0;

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		if (date.day == 0)
			indiv.anuptake = 0.0;

		nuptake_day           = indiv.ndemand * indiv.fnuptake;
		indiv.anuptake        += nuptake_day;
		indiv.nmass_leaf      += indiv.leaffndemand  * nuptake_day;
		indiv.nmass_root      += indiv.rootfndemand  * nuptake_day;
		indiv.nmass_sap       += indiv.sapfndemand   * nuptake_day;
		if (indiv.pft.phenology == CROPGREEN && ifnlim)
			indiv.cropindiv->nmass_agpool += indiv.storefndemand * nuptake_day;
		else
			indiv.nstore_longterm += indiv.storefndemand * nuptake_day;


		if (ifntransform) {
			double ammonium = nuptake_day * ammonium_frac;
			soil.NH4_mass -= ammonium;
			soil.NO3_mass -= nuptake_day - ammonium;
		} 
		else {
			soil.nmass_subtract(nuptake_day);
		}


		if (!negligible(indiv.phen))
			indiv.cton_leaf_aavr += min(indiv.cton_leaf(),indiv.pft.cton_leaf_max);

		vegetation.nextobj();
	}

	if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr && !negligible(orignmass)) {
		soil.fnuptake_survival[date.month] *=
			stock_survival(orignmass, soil.nmass_avail());
	}
}

/// Vegetation phosphorus uptake
/** Daily vegetation uptake of mineral phosphorus
*  Partitioned among individuals according to todays phosphorus demand
*  and individuals root area
*/
void vegetation_p_uptake(Patch& patch) {

	// Daily phosphorus demand given by:
	//	 For individual:
	//     (1)  pdemand = leafpdemand + rootpdemand + sappdemand;
	//          where
	//          leafpdemand is leaf demand based on vmax
	//			rootpdemand is based on optimal leaf C:P ratio
	//          sapnpemand is based on optimal leaf C:P ratio
	//
	// Actual phosphorus uptake for each day and individual given by:
	//     (2)  puptake = pdemand * fpuptake
	//     where fpuptake is individual uptake capacity calculated in fpuptake in canexch.cpp

	double puptake_day;

	Vegetation& vegetation = patch.vegetation;
	Soil& soil = patch.soil;

	const double origpmass = soil.pmass_labile + soil.pmass_sorbed;
	//double ammonium_frac = orignmass ? soil.NH4_mass / orignmass : 0;

	// Loop through individuals

	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv = vegetation.getobj();

		if (date.day == 0)
			indiv.apuptake = 0.0;

		puptake_day = indiv.pdemand * indiv.fpuptake;
		indiv.apuptake += puptake_day;
		indiv.pmass_leaf += indiv.leaffpdemand  * puptake_day;
		indiv.pmass_root += indiv.rootfpdemand  * puptake_day;
		indiv.pmass_sap += indiv.sapfpdemand   * puptake_day;
		if (indiv.pft.phenology == CROPGREEN && ifplim)
			indiv.cropindiv->pmass_agpool += indiv.storefpdemand * puptake_day;
		else
			indiv.pstore_longterm += indiv.storefpdemand * puptake_day;
		
		//soil.pmass_labile = max(0.0, soil.pmass_labile - puptake_day);
		//soil.pmass_labile_delta -= puptake_day;
		pmass_add(soil, -puptake_day);
				
		if (!negligible(indiv.phen))
			indiv.ctop_leaf_aavr += min(indiv.ctop_leaf(), indiv.pft.ctop_leaf_max);

		vegetation.nextobj();
	}

	if (date.year >= soil.solvesomcent_beginyr && date.year <= soil.solvesomcent_endyr && !negligible(origpmass)) {
		const double remaining = soil.pmass_labile + soil.pmass_sorbed;
		soil.fpuptake_survival[date.month] *=
			stock_survival(origpmass, remaining);
	}
}


/// Add litter to equilsom()
void add_litter(Soil& soil, int year, int pool) {

	// If LUC have occurred in between solvesomcent_beginyr and solvesomcent_endyr then array might be too smalll
	if (year >= (int)soil.solvesom.size())
		year = year%(int)soil.solvesom.size();

	// Add to litter pools
	soil.sompool[pool].cmass += soil.solvesom[year].get_clitter(pool);
	soil.sompool[pool].nmass += soil.solvesom[year].get_nlitter(pool);
	soil.sompool[pool].pmass += soil.solvesom[year].get_plitter(pool);
}

/// Iteratively solving differential flux equations for century SOM pools
/** Iteratively solving differential flux equations for century SOM pools
 *  assuming annual litter inputs, nitrogen uptake and leaching is close
 *  to long term equilibrium
 */
void equilsom(Soil& soil) {

	if (!(date.year == soil.solvesomcent_endyr && date.islastmonth && date.islastday))
		return;

	// Number of years to run SOM pools, value chosen to get cold climates to equilibrium
	const int EQUILSOM_YEARS = 40000;

	Patch& patch = soil.patch;
	const Climate& climate = soil.patch.get_climate();
	const Gridcell& gridcell = patch.stand.get_gridcell();

	// Save soil mineral nitrogen status
	double save_NH4_mass = soil.nmass_avail(NH4);
	double save_NO3_mass = 0.0;
	if (ifntransform) {
		save_NO3_mass = soil.nmass_avail(NO3);
	}

	// Save pmass_labile status
	double save_pmass_labile = soil.pmass_labile;
	double save_pmass_sorbed = soil.pmass_sorbed;
	double save_pmass_strongly_sorbed = soil.pmass_strongly_sorbed;
	double save_pmass_occluded = soil.pmass_occluded;

	// Number of years with mean input data
	int nyear = soil.solvesomcent_endyr - soil.solvesomcent_beginyr + 1;

	for (int m = 0; m < 12; m++) {

		// Monthly average decay rates
		for (int p = 0; p < NSOMPOOL; p++) {
			soil.sompool[p].mfracremain_mean[m] = pow(soil.sompool[p].mfracremain_mean[m] / nyear, (double)date.ndaymonth[m]);
		}

		// One representative month's survival is the geometric mean of the
		// sampled monthly operators. Daily survival fractions were already
		// composed while sampling.
		soil.fnuptake_survival[m] = pow(soil.fnuptake_survival[m], 1.0 / nyear);

		soil.fpuptake_survival[m] = pow(soil.fpuptake_survival[m], 1.0 / nyear);

		// Decay-weighted organic leaching. This reproduces the mass partition of
		// the daily SOILMICRO operator and keeps C, N and P on one coefficient.
		soil.morgleach_mean[m] = negligible(soil.msoilmicro_cdec[m]) ? 0.0 :
			soil.morgleach_cmass[m] / soil.msoilmicro_cdec[m];

		soil.mminleach_survival[m] = pow(soil.mminleach_survival[m], 1.0 / nyear);

		soil.mminpleach_survival[m] = pow(soil.mminpleach_survival[m], 1.0 / nyear);
	}

	// Annual average nitrogen fixation
	soil.anfix_mean /= nyear;

	// Annual phosphorus weathering
	soil.apwtr_mean /= nyear;
	soil.apdep_mean /= nyear;

	// Spin SOM pools with saved litter input, nitrogen addition and fractions of
	// nitrogen uptake and leaching for EQUILSOM_YEARS years with monthly timesteps
	for (int yr = 0; yr < EQUILSOM_YEARS; yr++) {

		// Which year in saved data set
		int savedyear = yr%nyear;

		// Monthly time steps
		for (int m = 0; m < 12; m++) {

			// Transfer yearly mean litter on first day of year
			add_litter(soil, savedyear*12+m, SURFSTRUCT);
			add_litter(soil, savedyear*12+m, SURFMETA);
			add_litter(soil, savedyear*12+m, SOILSTRUCT);
			add_litter(soil, savedyear*12+m, SOILMETA);
			add_litter(soil, savedyear*12+m, SURFFWD);
			add_litter(soil, savedyear*12+m, SURFCWD);

			// Calculate total litter carbon and nitrogen mass for set N:C ratio of surface microbial pool
			double litter_cmass = soil.sompool[SURFSTRUCT].cmass + soil.sompool[SURFMETA].cmass +
				soil.sompool[SURFFWD].cmass + soil.sompool[SURFCWD].cmass;
			double litter_nmass = soil.sompool[SURFSTRUCT].nmass + soil.sompool[SURFMETA].nmass +
				soil.sompool[SURFFWD].nmass + soil.sompool[SURFCWD].nmass;
			double litter_pmass = soil.sompool[SURFSTRUCT].pmass + soil.sompool[SURFMETA].pmass +
				soil.sompool[SURFFWD].pmass + soil.sompool[SURFCWD].pmass;

			// Set N:C ratio of surface microbial pool based on N:C ratio of litter from all PFTs
			if (!negligible(litter_cmass)) {
				setntoc(soil, litter_nmass / (litter_cmass * 2.0), SURFMICRO, 20.0, 10.0, 0.0, NCONC_SAT);
				setptoc(soil, litter_pmass / (litter_cmass * 2.0), SURFMICRO, 80.0, 30.0, 0.0, PCONC_SAT);
			}

			// Monthly nitrogen uptake
			soil.NH4_mass *= soil.fnuptake_survival[m];
			soil.NO3_mass *= soil.fnuptake_survival[m];

			// Monthly mineral nitrogen leaching
			double nmass = soil.nmass_avail(NO3);
			double nleach = nmass * (1.0 - soil.mminleach_survival[m]);
			soil.nmass_subtract(nleach, NO3);

			// Monthly phosphorus uptake
			pmass_apply_survival(soil, soil.fpuptake_survival[m]);

			// Monthly mineral phosphorus leaching
			pmass_apply_survival(soil, soil.mminpleach_survival[m]);

			// Monthly nitrogen addition to the system
			soil.nmass_inc((gridcell.aNH4dep + soil.anfix_mean) / 12.0, NH4);
			soil.nmass_inc(gridcell.aNO3dep / 12.0, NO3);

			// Monthly phosphorus addition to the system
			pmass_add(soil, (soil.apwtr_mean + soil.apdep_mean) / 12.0);

			// Monthly decomposition and fluxes between SOM pools

			// Set this months decay rates
			for (int p = 0; p < NSOMPOOL; p++) {
				soil.sompool[p].fracremain = soil.sompool[p].mfracremain_mean[m];
			}

			// Set this months organic nitrogen leaching fraction
			soil.orgleachfrac = soil.morgleach_mean[m];

			// Monthly decomposition and fluxes between SOM pools
			// and nitrogen flux from soil
			somfluxes(patch, true, false);
		}
	}

	// Reset mineral nitrogen status
	soil.NH4_mass = save_NH4_mass;
	soil.NO3_mass = save_NO3_mass;

	// Reset pmass_labile status
	soil.pmass_labile = save_pmass_labile;
	soil.pmass_sorbed = save_pmass_sorbed;
	soil.pmass_strongly_sorbed = save_pmass_strongly_sorbed;
	soil.pmass_occluded = save_pmass_occluded;

	// Reset variables for next equilsom()
	for (int m = 0; m < 12; m++) {

		for (int p = 0; p < NSOMPOOL; p++) {
			soil.sompool[p].mfracremain_mean[m] = 0.0;
		}

		soil.fnuptake_survival[m]  = 1.0;
		soil.fpuptake_survival[m]  = 1.0;
		soil.morgleach_mean[m] = 0.0;
		soil.morgleach_cmass[m] = 0.0;
		soil.msoilmicro_cdec[m] = 0.0;
		soil.morgPleach_mean[m] = 0.0;
		soil.mminleach_survival[m] = 1.0;
		soil.mminpleach_survival[m] = 1.0;
	}

	soil.anfix_mean = 0.0;
	soil.apwtr_mean = 0.0;
	soil.apdep_mean = 0.0;
	soil.solvesom.clear();
	std::vector<LitterSolveSOM>().swap(soil.solvesom); // clear array memory
}

/// SOM CENTURY DYNAMICS
/** To be called each simulation day for each modelled patch, following update
 *  of soil temperature and soil water.
 *  Transfers litter, performes nitrogen uptake and addition, leaching and decomposition.
 */
void som_dynamics_century(Patch& patch, Climate& climate, bool tillage) {

//	patch.soil.pmass_labile_delta = 0.0;

	// Transfer litter to SOM pools
	transfer_litter(patch);

	// Daily nitrogen uptake
	vegetation_n_uptake(patch);

	// Daily phosphorus uptake
	vegetation_p_uptake(patch);

	// Daily nitrogen addition to the soil
	soilnadd(patch);

	// Daily phosphorus addition to the soil
	soilpadd(patch);

	// Daily mineral and organic nitrogen leaching
	leaching(patch.soil);

	// Daily or monthly decomposition and fluxes between SOM pools
	somfluxes(patch, false, tillage);

	// Nitrogen transformation in soil
	ntransform(patch, climate);

	// Solve SOM pool sizes at end of year given by soil.solvesomcent_endyr
	equilsom(patch.soil);
}

/// Choose between CENTURY or standard LPJ SOM dynamics
/**
*/
void som_dynamics(Patch& patch, Climate& climate) {

	bool tillage = iftillage && patch.stand.landcover == CROPLAND;
	if (ifcentury) {
		som_dynamics_century(patch, climate, tillage);
	}
	else {
		som_dynamics_lpj(patch, tillage);
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Chatskikh, D., Hansen, S., Olesen, J.E. & Petersen, B.M. 2009. A simplified modelling approach
//	 for quantifying tillage effects on soil carbon stocks. Eur.J.Soil.Sci., 60:924-934.
// CENTURY Soil Organic Matter Model Version 5; Century User's Guide and Reference.
//  http://www.nrel.colostate.edu/projects/century5/reference/index.htm; accessed Oct.7, 2016.
// Comins, H. N. & McMurtrie, R. E. 1993. Long-Term Response of Nutrient-Limited
//   Forests to CO2 Enrichment - Equilibrium Behavior of Plant-Soil Models.
//   Ecological Applications, 3, 666-681.
// Cosby, B. J., Hornberger, C. M., Clapp, R. B., & Ginn, T. R. 1984 A statistical exploration
//   of the relationships of soil moisture characteristic to the physical properties of soil.
//   Water Resources Research, 20: 682-690.
// Cleveland C C (1999) Global patterns of terrestrial biological nitrogen (N2) fixation
//   in natural ecosystems. GBC 13: 623-645
// Foley J A 1995 An equilibrium model of the terrestrial carbon budget
//   Tellus (1995), 47B, 310-319
// Friend, A. D., Stevens, A. K., Knox, R. G. & Cannell, M. G. R. 1997. A
//   process-based, terrestrial biosphere model of ecosystem dynamics
//   (Hybrid v3.0). Ecological Modelling, 95, 249-287.
// Kirschbaum, M. U. F. and K. I. Paul (2002). "Modelling C and N dynamics in forest soils
//   with a modified version of the CENTURY model." Soil Biology & Biochemistry 34(3): 341-354.
// Meentemeyer, V. (1978) Macroclimate and lignin control of litter decomposition
//   rates. Ecology 59: 465-472.
// Mikan, C. J., Schimel, J. P., and Doyle, A. P. 2002. Temperature controls of microbial respiration 
//   in arctic tundra soils above and below freezing, Soil Biol.Biochem., 34, 1785–1795, 2002.
// Parton, W. J., Scurlock, J. M. O., Ojima, D. S., Gilmanov, T. G., Scholes, R. J., Schimel, D. S.,
//   Kirchner, T., Menaut, J. C., Seastedt, T., Moya, E. G., Kamnalrut, A. & Kinyamario, J. I. 1993.
//   Observations and Modeling of Biomass and Soil Organic-Matter Dynamics for the Grassland Biome
//   Worldwide. Global Biogeochemical Cycles, 7, 785-809.
// Parton, W. J., Hanson, P. J., Swanston, C., Torn, M., Trumbore, S. E., Riley, W. & Kelly, R. 2010.
//   ForCent model development and testing using the Enriched Background Isotope Study experiment.
//   Journal of Geophysical Research-Biogeosciences, 115.
// Schaefer, K. and Jafarov, E. 2016. 
//   A parameterization of respiration in frozen soils based on substrate availability
//   Biogeosciences, 13, 1991 - 2001, https ://doi.org/10.5194/bg-13-1991-2016
