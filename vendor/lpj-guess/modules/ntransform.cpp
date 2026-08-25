
///////////////////////////////////////////////////////////////////////////////////////
/// \file ntransform.cpp
/// \brief Nitrogen transformation in soil - nitrification and denitrification
///
/// \author Xu-Ri and modified for LPJ-guess by Peter Eliasson, David Wårlind and Stefan Olin.
/// $Date: 2013-10-14 14:12:00 +0100 (Mon, 10 Sep 2013) $
///
///////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////
// MODULE SOURCE CODE FILE
//
// Module:                Nitrogen transformation processes in Soil
// Header file name:      ntransform.h
// Source code file name: ntransform.cpp
// Written by:            Stefan Olin, adopted from Xu-Ri 2007-08-18.
// Version dated:         2019.
//
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
#include "guess.h"
#include "driver.h"
#include "ntransform.h"
#include <assert.h>
#include <numeric>

/// The NO and N2O fractions of the reduced NO2 flux, table 11 of Xu-Ri and
/// Prentice (2008), applied by table 9 eqns 5 and 6.
///
/// RNODN, the rate of NO production from denitrification, is given as under
/// 0.2 per cent with no mean stated, so the top of the stated bound is what is
/// taken. RN2ODN, the rate of N2O production from denitrification, is given as
/// 0.2 to 4.7 per cent with a stated mean of 2 per cent, which is what is
/// taken. The same rule settles the nitrification gas share in
/// global_soiln.ins: the paper's stated mean where it states one, the top of
/// the stated bound where it does not.
static const double R_NO_DENITRI  = 0.002;
static const double R_N2O_DENITRI = 0.02;

/// Soil NH3 volatilization
/** Daily calculation of NH3 volatilization from soil
 *
 */
void nh3_volatilization(Patch& patch, Soil& soil, double& n_budget_check){

	double f_nit_T = 0.0, nh3_inc = 0.0;
	double soil_T = soil.get_soil_temp_25();
	double wcont = soil.get_soil_water_upper();

	// The simulated soil's pH, which the soil map supplies per gridcell.
	//
	// This used to fall back to Dawson (1977)'s regression of soil pH on annual
	// precipitation, 3810 / (762 + climate.aprec_lastyear) + 3.5. The port's
	// own constant was wrong: table 5 eqn 5 of Xu-Ri and Prentice (2008) gives
	// 3810 / (762 + Precipitation_annual) + 3.8, and since pH enters
	// volatilisation as exp(2 * (pH - 10)), those 0.3 pH units are a factor of
	// 1.8 on the NH3 flux. The fallback could not work either: nothing in the
	// model ever assigns aprec_lastyear, and climate.aprec, which would feed
	// it, is reset at day 0 and never accumulated. The regression therefore
	// evaluated at zero precipitation and returned 8.5 on every gridcell for
	// every day of every run, where the paper's constant would have given 8.8. Two further things stood between it and a usable number even with a
	// live input: it is an Earth calibration, and its argument is an annual
	// precipitation sum, which on this world's shorter year is a smaller number
	// for the same precipitation rate and so reads as a drier, more alkaline
	// soil. Refuse instead of substituting a constant nobody chose: an input
	// path that supplies no pH has to say so.
	if (soil.soiltype.pH > 0.0) {
		soil.pH = soil.soiltype.pH;
	}
	else {
		fail("nh3_volatilization: ifntransform is set but the input path supplied "
		     "no soil pH for this gridcell. The soil map's ph column is what "
		     "carries it; a soil-code-only soil has none.");
	}

	// There is no nh3_max here any more. The port multiplied the whole
	// expression below by 0.001 above pH 6 and 0.00001 at or below it, and
	// table 5 of Xu-Ri and Prentice (2008) has no such factor. The paper's own
	// text says why it needs none: the NH3/NH4+ ratio "is virtually zero (about
	// 0.001) in solutions with pH below 6, but approaches unity (about 0.90) in
	// solutions with pH above 10", which is what eqn 6's fpH = exp(2 * (pH -
	// 10)) alone delivers, 3.4e-4 at pH 6 and 1.0 at pH 10. nh3_max applied
	// that ratio a second time and added a factor-100 step at pH 6 on top of a
	// smooth exponential. Over the pH range the simulated soil map carries, it
	// held NH3_SOIL three orders of magnitude below the paper's formulation on
	// the gridcells above pH 6 and five orders below it on the rest.

	// N budget check
	// Total budget before N transformations
	n_budget_check = soil.NH4_mass + soil.NO3_mass + soil.NO2_mass + soil.NO_mass + soil.N2O_mass + soil.N2_mass;

	// temperature limiting factor for volatilization, 25 degree C == 1 (table 5, eqn 7, Xu-Ri 2008)
	if (soil_T >= -40.0) {
		f_nit_T = min(1.0, exp(308.56 * (1.0 / 71.02 - 1.0 / (soil_T + 46.02))));
	}
	else {
		f_nit_T = 0.0;
	}

	// NH3 increment (table 5, eqn 2, 3, 4, 6, Xu-Ri 2008)
	nh3_inc = min(soil.NH4_mass, (min(1.0, wcont) * (1.0 - min(1.0, wcont))) * f_nit_T * f_nit_T * exp(2.0 * (soil.pH - 10.0)) * soil.NH4_mass);

	soil.NH4_mass -= nh3_inc;

	patch.fluxes.report_flux(Fluxes::NH3_SOIL, nh3_inc);

	// N budget check
	// N lost as NH3 gas
	n_budget_check -= nh3_inc;

}

/// Soil N substrate partitioning
/** Daily calculation of nitrogen substrate partition 
 *  between aerobic and anaerobic soil fractions
 */
void substrate_partition(Soil& soil){

	// An S-curve fitted to the linear increase in
	// denitrifier activity above roughly 50% WFPS with the steepest
	// change around 66% ending up at 100% around 75%.
	// Meaning that there is no nitrification at that WFPS.
	// Derived from Pilegaard 2013
	//
	// The argument is water-filled pore space, and it is what the split is
	// declared in. Xu-Ri and Prentice (2008) table 7 names the aeration
	// variable DyN allocates substrates on as "soil moisture as water-filled
	// pore space (WFPS)", and its nitrification section as the WFPS of the top
	// 50 cm, which is this layer. nitrification and denitrification below
	// already read soil.wfps(0); this function used to read
	// get_soil_water_upper(), a fraction of AVAILABLE capacity, which is a
	// different quantity with a texture-dependent offset and scale. Feeding it
	// here made the aerobic/anaerobic split a function of texture at fixed soil
	// wetness: on the current Vesper soil map, a cell at 0.50 WFPS got an
	// anaerobic share anywhere from 0.07 to 0.87 depending only on its sand and
	// clay. The curve's own midpoint and shape are the port's and are not in
	// Xu-Ri; they are registered as unsourced with their bracket in
	// biosphere/config/ntransform.yaml.
	//
	// This does NOT give the curve its full domain. wcont is bounded above by
	// field capacity, so soil.wfps(0) tops out at field capacity over
	// saturation -- a median of 0.81 and as little as 0.57 on that map. Water
	// above field capacity is the quantity the model does not carry.

	double wfps = soil.wfps(0);
	double wet = richards_curve(0.05,0.95,7.5,0.5,wfps);

	soil.NH4_mass_w = soil.NH4_mass * wet;
	soil.NO3_mass_w = soil.NO3_mass * wet;
	soil.NO2_mass_w = soil.NO2_mass * wet;
	soil.N2O_mass_w = soil.N2O_mass * wet;
	soil.NO_mass_w = soil.NO_mass * wet;
	soil.labile_carbon_w = soil.labile_carbon * wet;

	soil.NH4_mass_d = soil.NH4_mass - soil.NH4_mass_w;
	soil.NO3_mass_d = soil.NO3_mass - soil.NO3_mass_w;
	soil.NO2_mass_d = soil.NO2_mass - soil.NO2_mass_w;
	soil.N2O_mass_d = soil.N2O_mass - soil.N2O_mass_w;
	soil.NO_mass_d  = soil.NO_mass  - soil.NO_mass_w;
	soil.labile_carbon_d = soil.labile_carbon - soil.labile_carbon_w;

}

/// Soil N nitrification
/** Daily calculation of nitrification, and nitrification
 *  induced trace gas emissions 
 */
void nitrification(Patch& patch, Soil& soil) {

	double b = log(3.0) * 5.0;
	double a = exp(-b*6.0/10.0);
	double wcont = soil.get_soil_water_upper();
	double act_dry = a*exp(wcont*b);
	double act_wet = max(0.0,4.0-5.0*wcont);
	double nit_act = min(act_dry,act_wet);

	double f_nit_T, no3_inc, no_inc, n2o_inc, gross_nitrif;
	double soil_T = soil.get_soil_temp_25();

	// f_nit_T, temperature limiting factor for nitrification, 38 deg C == 1 (table 8, eqn 2, Xu-Ri 2008)
	if (soil_T < 70.0) {
		f_nit_T = min(1.0, pow((70.0 - soil_T) / (70.0 - 38.0), 12.0) * exp(12.0 * (soil_T - 38.0) / (70.0 - 38.0)));
	}
	else {
		f_nit_T = 0.0;
	}
	// gross nitrification rate, NH4_mass converted to NO3_mass (table 8, eqn 1, Xu-Ri 2008)
	no3_inc          = f_nitri_max * nit_act * f_nit_T * soil.NH4_mass_d;
	gross_nitrif     = no3_inc;
	soil.NH4_mass_d -= no3_inc;

	// The share of gross nitrification emitted as NO plus N2O. Table 8 eqns 3
	// and 4 of Xu-Ri and Prentice (2008) are NOinc = RNON * NO3inc and N2Oinc =
	// RN2ON * NO3inc, so the two together take RNON + RN2ON of the flux, and
	// f_no below splits what leaves.
	//
	// This read f_denitri_gas_max, which global_soiln.ins, parameters.h,
	// parameters.cpp and the declareitem help string all describe as the
	// maximum gaseous loss in DENITRIFICATION, while denitrification below read
	// f_nitri_gas_max, which those same four describe as the maximum gaseous
	// loss in NITRIFICATION. Four declarations agreed with the names and only
	// the two reads disagreed, and the instruction file's values corroborate
	// them: 0.33 for both steps of the reduction sequence against 0.1 and the
	// nitrification gas share. So the names were right, the reads were crossed,
	// and it is the reads that moved.
	double ngas_inc = f_nitri_gas_max * no3_inc;

	double wfps = soil.wfps(0);

	// Pilegaard 2013
	// only NO and N2O in nitrification 
	double f_no = richards_curve(1.0, 0.5, 20.0, 0.375, wfps);

	no_inc          = f_no * ngas_inc;
	no3_inc        -= no_inc;
	soil.NO_mass_d += no_inc;

	n2o_inc          = ngas_inc - no_inc;
	no3_inc         -= n2o_inc;
	soil.N2O_mass_d += n2o_inc;

	// New NO3-
	soil.NO3_mass_d += no3_inc;

	soil.NO3_mass_d += soil.NO2_mass_d;
	soil.NO2_mass_d = 0.0;

	// Report gross nitrification
	patch.fluxes.report_flux(Fluxes::GROSS_NITRIF, gross_nitrif);
}

/// Soil N denitrification
/** Daily calculation of denitrification rate, and denitrification
 *  induced trace gas emissions 
 */
void denitrification(Patch& patch,Soil& soil) {
	double soil_T = soil.get_soil_temp_25();
	double wcont = soil.get_soil_water_upper();
	// used to remove the per m3 from the constants KC and KN.
	double water_cont_m3 = wcont * soil.soiltype.gawc[0] / 1000.0;

	double wfps_upper = soil.wfps(0);

	double f_den_T, d_N_max, no2_inc, no_inc, n2o_inc, ngas_inc, gross_denitrif, n2_inc;
	if (water_cont_m3 > 0.0 && wfps_upper>0.4) {
		// temperature limiting factor for denitrification, 22 deg C == 1 (table 9, eqn 1, Xu-Ri 2008)
		//
		// Table 9 eqn 1 carries no min{1, .}, unlike table 5 eqn 7 and table 10
		// eqn 1 which both state one, and the paper's function rises past 1
		// above 22 C: 1.21 at 25 C, 1.61 at 30 C, 2.07 at 35 C, 3.15 at 45 C.
		// The port clamped it, which held denitrification in warm soil at its
		// 22 C rate. What keeps the operator conservative is not that clamp but
		// the min() against the pool on each of the two transformation lines
		// below: no2_inc cannot exceed NO3_mass_w and ngas_inc cannot exceed
		// NO2_mass_w whatever the coefficient reaches. The guard at -40 C stays
		// because the exponent diverges at -46.02 C; the function is 5e-21 there.
		if (soil_T >= -40.0) {
			f_den_T = exp(308.56 * (1.0 / 68.02 - 1.0 / (soil_T + 46.02)));
		}
		else {
			f_den_T = 0.0;
		}
		// Effect of labile carbon availability on denitrification (table 9, eqn 2, Xu-Ri 2008)
		d_N_max = soil.labile_carbon_w / (k_C * water_cont_m3 + soil.labile_carbon_w);

		// Gross denitrification ratio NO3 to NO2 (table 9, eqn 3, Xu-Ri 2008)
		no2_inc = min(soil.NO3_mass_w, soil.NO3_mass_w * f_denitri_max * d_N_max * f_den_T * soil.NO3_mass_w / (k_N * water_cont_m3 + soil.NO3_mass_w));

		gross_denitrif   = no2_inc;
		soil.NO3_mass_w -= no2_inc;
		soil.NO2_mass_w += no2_inc;

		// Gross transformation of NO2 to N2 (table 9, eqn 4, Xu-Ri 2008)

		// Denitrification rate dependence on moisture, Weier et al. 1993
		double f_den_w = min(1.0, exp(13.0360 * wfps_upper - 11.6219));

		ngas_inc = min(soil.NO2_mass_w, soil.NO2_mass_w * f_denitri_gas_max * d_N_max * f_den_w * f_den_T * soil.NO2_mass_w / (k_N * water_cont_m3 + soil.NO2_mass_w));

		soil.NO2_mass_w -= ngas_inc;

		// NO, N2O and N2 out of the reduced NO2: table 9 eqns 5, 6 and 7 of
		// Xu-Ri and Prentice (2008). NOinc = RNODN * ftemp * N2inc, N2Oinc =
		// RN2ODN * ftemp * N2inc, and N2 is what is left, with ftemp the same
		// table 9 eqn 1 response used above.
		//
		// This replaces a branch at 0.7 WFPS that produced no N2 below it and
		// no NO above it, built from three functions neither paper states: an
		// N2O:NO regression attributed to Weier et al. (1993), which measured
		// N2 and N2O by acetylene block and never measured NO; an N2O:N2 curve;
		// and a logistic in soil temperature, from a study incubated at one
		// temperature. The hard zero below 0.7 WFPS also contradicts the paper
		// the moisture response beside it is taken from: Weier tables 4 and 5
		// measured N2/N2O at 60 and 75 per cent WFPS across four soils and nine
		// carbon-by-nitrate treatments and found N2 the majority product in
		// most of those observations.
		//
		// What this does NOT settle is the size of the N2 share. At 22 C these
		// equations leave 97.8 per cent of the reduced nitrogen as N2, where
		// Weier's medians are 0.565, 0.744 and 0.796 at 60, 75 and 90 per cent
		// WFPS. The two Earth calibrations disagree, the operator runs at the
		// DyN end because every other equation in it is DyN's, and the
		// disagreement is declared as a bracket in
		// biosphere/config/ntransform.yaml.
		no_inc  = R_NO_DENITRI * f_den_T * ngas_inc;
		n2o_inc = R_N2O_DENITRI * f_den_T * ngas_inc;
		n2_inc  = ngas_inc - no_inc - n2o_inc;

		soil.NO_mass_w += no_inc;
		soil.N2O_mass_w	+= n2o_inc;
		soil.N2_mass += n2_inc;
	} 
	else {
		gross_denitrif = 0.0;
	}
	// Report fluxes of gross denitrification
	patch.fluxes.report_flux(Fluxes::GROSS_DENITRIF, gross_denitrif);
}

/// Soil N gas emissions
/** Daily calculation of soil N gas emissions. 
 */
void n_gas_emission(Patch& patch, Fluxes& fluxes, Soil& soil, double& n_budget_check) {

	//for book keeping
	soil.labile_carbon = soil.labile_carbon_w  + soil.labile_carbon_d;

	double soil_T = soil.get_soil_temp_25();
	double wcont = soil.get_soil_water_upper();
	double ftemp, no_d_flux_inc, n2o_d_flux_inc, no_w_flux_inc, n2o_w_flux_inc, n2_flux_inc;
	double net_nitrif = 0.0;
	double net_denitrif = 0.0;

	// ftemp, temperature limiting factor for gas emission, 25 degree C == 1 (table 10, eqn 1, Xu-Ri 2008)
	if (soil_T >= -40.0) {
		ftemp = min(1.0, exp(308.56 * (1.0 / 71.02 - 1.0 / (soil_T + 46.02))));
	}
	else {
		ftemp = 0.0;
	}

	// Nitrification fluxes
	// Daily NO gas released from aerobic soil to atmosphere (table 10, eqn 2, Xu-Ri 2008)
	no_d_flux_inc   = ftemp * (1.0 - min(1.0, wcont)) * soil.NO_mass_d;
	soil.NO_mass_d -= no_d_flux_inc;
	net_nitrif     += no_d_flux_inc;
	// Daily N2O gas released from aerobic soil to atmosphere (table 10, eqn 2, Xu-Ri 2008)
	n2o_d_flux_inc   = ftemp * (1.0 - min(1.0, wcont)) * soil.N2O_mass_d;
	soil.N2O_mass_d -= n2o_d_flux_inc;
	net_nitrif      += n2o_d_flux_inc;


	// Denitrification fluxes
	// Daily NO gas released from anaerobic soil to atmosphere (table 10, eqn 2, Xu-Ri 2008)
	no_w_flux_inc   = ftemp * (1.0 - min(1.0, wcont)) * soil.NO_mass_w;
	soil.NO_mass_w -= no_w_flux_inc;
	net_denitrif   += no_w_flux_inc;
	// Daily N2O gas released from anaerobic soil to atmosphere (table 10, eqn 2, Xu-Ri 2008)
	n2o_w_flux_inc   = ftemp * (1.0 - min(1.0, wcont)) * soil.N2O_mass_w;
	soil.N2O_mass_w -= n2o_w_flux_inc;
	net_denitrif    += n2o_w_flux_inc;

	// Daily N2 gas released from soil to atmosphere (table 10, eqn 2, Xu-Ri 2008)
	n2_flux_inc     = ftemp * (1.0 - min(1.0, wcont)) * soil.N2_mass;
	soil.N2_mass   -= n2_flux_inc;

	patch.fluxes.report_flux(Fluxes::NO_SOIL,     no_d_flux_inc + no_w_flux_inc);
	patch.fluxes.report_flux(Fluxes::N2O_SOIL,    n2o_d_flux_inc + n2o_w_flux_inc);

	patch.fluxes.report_flux(Fluxes::N2_SOIL,     n2_flux_inc); 
	patch.fluxes.report_flux(Fluxes::NET_NITRIF,  net_nitrif);
	patch.fluxes.report_flux(Fluxes::NET_DENITRIF,net_denitrif);

	// substrates after N gas emissions
	soil.NH4_mass = soil.NH4_mass_w  + soil.NH4_mass_d;
	soil.NO3_mass = soil.NO3_mass_w  + soil.NO3_mass_d;
	soil.NO2_mass = soil.NO2_mass_w  + soil.NO2_mass_d;
	soil.NO_mass  = soil.NO_mass_w   + soil.NO_mass_d;
	soil.N2O_mass = soil.N2O_mass_w  + soil.N2O_mass_d;
	soil.labile_carbon = soil.labile_carbon_w  + soil.labile_carbon_d;

	// N budget check
	// N lost through emissions and existing pool sizes
	n_budget_check -= (no_d_flux_inc + no_w_flux_inc + n2o_d_flux_inc + n2o_w_flux_inc + n2_flux_inc + 
	                   soil.NH4_mass + soil.NO3_mass + soil.NO2_mass + soil.NO_mass + soil.N2O_mass + soil.N2_mass);
}


/// N transformation processes in soil
/** Daily calculation of nitrification, denitrification, and trace gases emissions. 
 *  To be called each simulation day for each modelled area or patch, following update
 *  of soil organic matter dynamic submodel.
 */
void ntransform(Patch& patch, Climate& climate) {
	if (ifntransform) {
		const double EPS = 1.0e-14;
		double n_budget_check;

		// Obtain reference to Soil object
		Soil& soil = patch.soil;

		// Obtain referenc to flux object
		Fluxes& fluxes = patch.fluxes;

		// NH3 volatilization
		nh3_volatilization(patch, soil, n_budget_check);

		// N substrate partition
		substrate_partition(soil);

		// Nitrification
		nitrification(patch, soil);

		// Denitrification
		denitrification(patch, soil);

		// N gas emission
		n_gas_emission(patch, fluxes, soil, n_budget_check);

		// The operator's conservation identity: the nitrogen held in the six soil
		// pools before it runs equals what it holds after, plus the NH3, NO, N2O
		// and N2 it emitted. nh3_volatilization opens n_budget_check with the
		// pool sum and n_gas_emission closes it, so a nonzero residual is
		// nitrogen created or destroyed.
		//
		// This was an assert(), which every Release build compiles out under
		// NDEBUG, so on the binary this project builds the only check on the
		// operator did not exist. The bar is relative to the nitrogen actually
		// present, floored at the absolute EPS the assert used so an empty soil
		// is still held to something. The identity is about twenty additions and
		// subtractions of doubles, whose rounding is of order 1e-16 relative, so
		// a relative 1e-10 is four orders of slack and a residual above it is a
		// defect in the operator rather than arithmetic. The bar is set here
		// rather than after a result, because there is no result yet.
		const double REL_EPS = 1.0e-10;
		double n_total = soil.NH4_mass + soil.NO3_mass + soil.NO2_mass +
		                 soil.NO_mass + soil.N2O_mass + soil.N2_mass;
		if (fabs(n_budget_check) > max(EPS, REL_EPS * n_total)) {
			fail("ntransform: the soil nitrogen transformation operator did not "
			     "conserve mass on day %d. Residual %g kgN/m2 against %g held.",
			     date.day, n_budget_check, n_total);
		}
	}
}



//////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Xu-Ri & Prentice IC 2008 Terrestrial nitrogen cycle simulated with a
//    dynamic global vegetation model. Global Change Biology 14,1745-1764.
// Dawson, G. A. (1977). "Atmospheric Ammonia from Undisturbed Land." 
//    Transactions-American Geophysical Union 58(6): 554-554.
// Pilegaard K. 2013.  Processes regulating nitric oxide emissions from soils.
// 	  Philosophical Transactions of the Royal Society of London B: Biological Sci-
//    ences, 368(1621), 2013. ISSN 0962-8436. doi: 10.1098/rstb.2013.0126.
// K. L. Weier et al. 1993. Denitrification and the dinitrogen/nitrous oxide
//    ratio as affected by soil water, available carbon, and nitrate.
//    Soil Science Society of America Journal, 57 (1):66–72.



