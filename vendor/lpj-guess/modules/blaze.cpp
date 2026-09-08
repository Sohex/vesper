///////////////////////////////////////////////////////////////////////////////////////
/// \file blaze.cpp
/// \brief BLAZE (BLAZe induced biosphere-atmosphere flux Estimator) 
///
/// \author Lars Nieradzik
/// $Date: 2017-01-24 17:03:10 +0100 (Tue, 24 Jan 2017) $
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
#include "blaze.h"
#include "driver.h"
#include "growth.h"
#include "somdynam.h"
#include "simfire.h"
#include "guessmath.h"

/* combustion rates depending on fire-line-intensities [kW/m]
 *       <750|<3000|<7000| >=7000 (Sprouters|Seeders)
 *
 * COLUMN 3 IS UNREACHABLE. get_fire_line_intensity_index() returns 4 for every
 * intensity above 7000 and never 3, so every intense fire takes the SEEDER
 * rates -- which for rows 4 and 5, stems and branches to litter, are 0.8
 * against the sprouter 0.2. is_resprouter exists in this file but reaches only
 * survival_probability_temp_broadleaf(); it never selects a column here, and it
 * could not without restructuring, because the index is per PATCH while
 * resprouting is a per-individual trait and get_combustion_rates() writes
 * patch-level fluxes.
 *
 * Left as it is rather than repaired, because both repairs are modelling
 * decisions and not corrections: making column 3 reachable changes which rates
 * intense fires use, and deleting it discards a distinction the source
 * intended. FIRE-6 owns that choice. world-voxz,
 * biosphere/notes/fire-model-audit.md finding 11.
 */
const double TURNOVERFRACT[13][5] = {
	{ .0 ,  .0 ,  .05, .2 , .2 }, //   0 Stems       -> ATM
	{ .0 ,  .0 ,  .15, .2 , .2 }, //   1 Branches    -> ATM
	{ .03,  .13,  .25, .5 , .5 }, //   2 Bark	 -> ATM
	{ .02,  .05,  .1 , .6 , .6 }, //   3 Leaves      -> ATM
	{ .0 ,  .0 ,  .05, .2 , .8 }, //   4 Stems       -> Litter (DWD) !corrected*
	{ .0 ,  .02,  .07, .2 , .8 }, //   5 Branches    -> Litter (CWD) !corrected*
	{ .03,  .13,  .25, .5 , .5 }, //   6 Bark        -> Litter (str)
	{ .05,  .1 ,  .15, .3 , .4 }, //   7 Leaves      -> Litter (str)
	{ .0 ,  .02,  .02, .04, .04}, //   8 FDEAD roots -> ATM
	{ .5 ,  .75,  .75, .8 , .8 }, //   9 CWD         -> ATM
	{ .6 ,  .65,  .85, 1. , 1. }, //  10 Bark Litter -> ATM
	{ .6 ,  .65,  .85, 1. , 1. }, //  11 Leaf Litter -> ATM*
	{ .0 ,  .0 ,  .1 , .8 , .8 }, //  12 Deadwood    -> ATM
};

// Tuning factors for litter ready for combustion
// Boreal
const double K_LITTER_BOREAL    = 0.38;
// Temperate region
const double K_LITTER_TEMPERATE = 0.0025;
// Tropics
const double K_LITTER_TROPICS   = 0.15;
// Savanna
const double K_LITTER_SAVANNA   = 0.75 ;

// Grassy vegetation burn-rate for cohort and individual mode
const double MAX_GRASS_BURN = 0.75;

// fraction of life woody biomass that is branch
const double F_BRANCH   = 0.05;
// fraction of life woody biomass that is bark
const double F_BARK     = 0.01;
// minimum available fuel to start a fire [gC/m2]
const double MIN_FUEL   = 200.;

const double RAINFALL_AVERAGING_SPAN = 3;
	
// Internal help function for splitting up nitrogen fire fluxes into components
// Copy of report_fire_nfluxes() as used in fire() in vegdynam.cpp
void report_fire_flux_n(Patch& patch, double nflux_fire) {
	patch.fluxes.report_flux(Fluxes::NH3_FIRE, Fluxes::NH3_FIRERATIO * nflux_fire);
	patch.fluxes.report_flux(Fluxes::NOx_FIRE, Fluxes::NOx_FIRERATIO * nflux_fire);
	patch.fluxes.report_flux(Fluxes::N2O_FIRE, Fluxes::N2O_FIRERATIO * nflux_fire);
	patch.fluxes.report_flux(Fluxes::N2_FIRE,  Fluxes::N2_FIRERATIO  * nflux_fire);
}

// pixelsize() WAS HERE AND IS DELETED. It computed a gridcell's area in square
// km from the spherical-segment formula S = 2*pi*r*h, and its r was R_EARTH --
// Earth's mean radius, 6371.2213 km, from guessmath.h. config/planet.yaml
// declares this world at 1.20 Earth radii, so a fire model that started
// reporting area on it would have reported against a planet 44 per cent smaller
// in surface area, with nothing objecting: the function carried no units check
// and the error is second order in the radius.
//
// It had no callers anywhere in the vendored tree, so nothing computed on it and
// deleting it changes no number. It is deleted rather than converted because
// FIRE-6 arms this layer and nothing the port needs calls for it -- FIRE-5's
// area per fire comes from a spread rate times a duration, not from a cell's
// size. If a gridcell area is wanted later it reaches this model the way every
// other planet constant does, through the header build_vesper_header.py
// generates, and lib/nc_geometry.py owns the one expression turning
// config/planet.yaml's radius into metres.
//
// This was the only consumer of R_EARTH in the tree, and that constant goes with
// it. world-pqk7, biosphere/notes/fire-model-audit.md finding 6.

/* Get combustion rates
 * Compute the relative flux rates [frac.] between live vegetation, litter pools and
 * atmosphere given current fire-line intensity.
 */
void get_combustion_rates(Patch& patch, int fli_index, double k_tun_litter) {


	// relative fluxes from wood to atmosphere and litter pools
	patch.wood_to_atm = (1.-F_BRANCH-F_BARK) * TURNOVERFRACT[ 0][fli_index] +
		F_BRANCH * TURNOVERFRACT[ 1][fli_index] +
		F_BARK   * TURNOVERFRACT[ 2][fli_index];
	patch.wood_to_str = F_BARK               * TURNOVERFRACT[ 6][fli_index];
	patch.wood_to_fwd = F_BRANCH             * TURNOVERFRACT[ 5][fli_index];
	patch.wood_to_cwd = (1.-F_BRANCH-F_BARK) * TURNOVERFRACT[ 4][fli_index];
	
	// relative fluxes from leaf to atmosphere and litter pools
	patch.leaf_to_atm = TURNOVERFRACT[ 3][fli_index];
	patch.leaf_to_lit = TURNOVERFRACT[ 7][fli_index];

	// relative fluxes from litter pools to atmosphere
	patch.litf_to_atm = TURNOVERFRACT[11][fli_index];
	patch.lfwd_to_atm = TURNOVERFRACT[10][fli_index];
	patch.lcwd_to_atm = TURNOVERFRACT[ 9][fli_index] * k_tun_litter;
	return;
}

/* Turn fire-line intensity into an index for lookup-tables.
 * Get appropriate FLI category for look-up tables depending on 
 * computed potential FLI. The index corresponding to the entries in
 * the look-up-tables is returned.
 */
int get_fire_line_intensity_index(double fire_line_intensity) {

	// determine intensity category for combustion-lookup-tables
	int fire_line_intensity_index; 
	if ( fire_line_intensity > 7000. ) {
		fire_line_intensity_index = 4;
	}
	else if ( fire_line_intensity > 3000. ) {
		fire_line_intensity_index = 2;
	}
	else if ( fire_line_intensity > 750. ) {
		fire_line_intensity_index = 1;
	}
	else if ( fire_line_intensity > 0. ) {
		fire_line_intensity_index = 0;
	} 
	else {
		fire_line_intensity_index = -1;
	}
	return fire_line_intensity_index;
}

// Compute the amount of fuel readily available to burn
double available_fuel (Patch& patch,int fli_index, double k_tun_litter)  {
			
	get_combustion_rates(patch,fli_index,k_tun_litter);

	// Transitional-litter available to burn, 
	// i.e. last year's litter that is still to fall.
	double trans_litter_leaf  = 0.;
	double trans_litter_sap   = 0.;
	double trans_litter_heart = 0.;
	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& patchpft = patch.pft.getobj();
		trans_litter_leaf  += patchpft.cmass_litter_leaf;
		trans_litter_sap   += patchpft.cmass_litter_sap;
		trans_litter_heart += patchpft.cmass_litter_heart;
		patch.pft.nextobj();
	}

	// compute readily available fuel load for given FLI-index in PATCH
	double available_fuel = patch.litf_to_atm * (patch.soil.sompool[SURFSTRUCT].cmass + 
						  patch.soil.sompool[SURFMETA].cmass +
						  trans_litter_leaf)
	+ patch.lfwd_to_atm * ( patch.soil.sompool[SURFFWD].cmass + trans_litter_sap )
	+ patch.lcwd_to_atm * ( patch.soil.sompool[SURFCWD].cmass + trans_litter_heart) * k_tun_litter;
	
	Vegetation& vegetation=patch.vegetation;
	vegetation.firstobj();
	while (vegetation.isobj) {
		Individual& indiv=vegetation.getobj();
		if (indiv.pft.lifeform == GRASS) {
			available_fuel += 0.5 * indiv.cmass_leaf_today() ;
		}
		vegetation.nextobj();
	}
	return available_fuel;
}

/* Compute current potential fire-line intensity under given meteorological and fuel conditions.
 * Formulation following Noble 1980 derived from McArthur.
 */
void get_fireline_intensity(Patch& patch, Climate& climate) {
	
	// Energy contents of fuel [MJ/kg]. Liedloff and Cook (2007) section 2.8
	// states it as an ASSUMPTION rather than a measurement -- "the heat yield of
	// the fuel (H: assumed to be 20 MJ kg-1)" -- so the chain terminates there.
	// Bracketed at [16, 22] as fireline_intensity_heat_yield in
	// biosphere/config/fire_parameters.yaml, on the enthalpy of combustion of
	// dry plant tissue rather than on that paper.
	const double HEAT_YIELD = 20.;

	// NOT AN EMPIRICAL VALUE. This is Noble, Gill and Bary (1980)'s Mark 5
	// forest rate of spread, R = 0.0012*F*W, converted into the units this
	// function works in. Their Table 1 fixes those: R in km/h, W in tonnes/ha,
	// F dimensionless. Here fuel is g/m2 and the rate is m/s, and 1 t/ha is
	// 100 g/m2 while 1 km/h is 1/3.6 m/s, so the coefficient is
	// 0.0012 / 100 / 3.6 = 3.333333e-06.
	//
	// It was declared as 3.3333e-05 under the comment "Empirical value", which
	// is Noble's conversion with a lost decimal: 1.2/3.6 is 0.33333, and
	// 0.33333e-5 written as 3.3333e-05 drops a decade. The mantissa agreeing to
	// five digits is what identifies it as a slip rather than a different
	// calibration. At F = 20 with 2000 g/m2 of fuel it gave 53 333 kW/m against
	// Noble's 5333, which is 7.6 times above the top fire-line-intensity class,
	// so every such fire took the maximum combustion column and the maximum
	// mortality. world-c15e, biosphere/notes/fire-model-audit.md finding 10.
	//
	// Written as the conversion rather than as a corrected constant so the
	// derivation is on the line. biosphere/scripts/fire_parameter_gate.py
	// evaluates the same expression and refuses if the two disagree.
	const double NOBLE_ROS_KMH_PER_TONNE_HA = 0.0012;   // Noble et al. (1980), Mark 5 forest
	const double G_M2_PER_TONNE_HA          = 100.0;    // 1 t/ha = 100 g/m2
	const double KMH_PER_M_S                = 3.6;      // 1 m/s  = 3.6 km/h
	const double A = NOBLE_ROS_KMH_PER_TONNE_HA / G_M2_PER_TONNE_HA / KMH_PER_M_S;
	// Readily available fuel  [g/m2]
	double avail_fuel;                      
	// Rate of spread          [m/s]
	double rate_of_spread;                    
	// Fire-line intensity     [kW/m]
	double fire_line_intensity;                  
	// Fire intensity category index
	int fire_line_intensity_index = 0;

	Gridcell& gridcell = climate.gridcell;
	
	// Get available fuel for current fire-line intensity index (fli_index) and convert kg/m2 to g/m2
	avail_fuel = available_fuel(patch,fire_line_intensity_index,gridcell.k_tun_litter) * G_PER_KG;
	// Check whether there is enough fuel to ignite a fire
	if ( avail_fuel < MIN_FUEL ) { 
		fire_line_intensity  =  -1. ;
	}
	else {
		// Compute Rate-of-spread [m/s]
		rate_of_spread = A * climate.mcarthur_forest_fire_index * avail_fuel;
		
		// Fire line intensity[W/m] (Pyne, 1996 derived from Byram, 1959)
		fire_line_intensity = HEAT_YIELD * avail_fuel * rate_of_spread;
	}
	patch.fire_line_intensity = fire_line_intensity;
}
 
// Survival probability for boreal trees based on Dalziel et al. 2008
double survival_probability_boreal(double fire_line_intensity) {

	double surv_prob_boreal = exp(-fire_line_intensity/500.);
	return surv_prob_boreal;
}

// Survival probability for temperate Needleleaf trees following Kobziar 2006
double survival_probability_temp_needleleaf(double diameter_at_breast_height, double fire_line_intensity, double mass_cwd) {

	double diameter_at_breast_height_cm  = diameter_at_breast_height * CM_PER_M; // in cm
	double cwd = mass_cwd * 0.1; // in Mg/ha
	// survival probability at intensity below 750 kW/m
	double survival_probability_750;
	double survival_probability;

	if ( fire_line_intensity < 750. ) {
		survival_probability_750   = 1. - (1./(1.+ exp(-(1.0337 + 0.000151*750. 
						- .221*diameter_at_breast_height_cm + .0219*cwd))));
		survival_probability = 1. - (fire_line_intensity/750. * (1. - survival_probability_750) );
	}
	else {
		survival_probability = 1. - (1./(1.+ exp(-(1.0337 + 0.000151*fire_line_intensity
						- .221*diameter_at_breast_height_cm + .0219*cwd))));
	}

	return survival_probability;
}

// Survival probability for temperate broadleaf trees following Hickler 2004
/* Compute survival probability for Temperate Broadleaved forests
 * fire resilince parameterisation for e.g. Australian forests
 * following Hickler et al. 2004, using a generalized
 * vegetation model to simulate vegetation dynamics in NE USA.
 */
double survival_probability_temp_broadleaf(double diameter_at_breast_height, double fire_line_intensity, bool is_resprouter) {
	
	// Fire resiliance
	double resilience;
	if ( is_resprouter ) {
		resilience = 0.04;
	} else {
		resilience = 0.07;
	}

	// Compute survival probability at 3000kW/m first
	double p_surv_3000 = 0.95 - 1./(1.+ pow((diameter_at_breast_height/resilience),1.5)) ;
	double survival_probability_temp_broadleaf;
	if ( fire_line_intensity > 7000. ) {
		survival_probability_temp_broadleaf = 0.001;
	}
	else if ( fire_line_intensity > 3000 ) { 
		survival_probability_temp_broadleaf = p_surv_3000 * (1. - (fire_line_intensity-3000.)/ 4000. );
	}
	else {
		survival_probability_temp_broadleaf = exp(fire_line_intensity/3000. * log(p_surv_3000)); 
	}

	return survival_probability_temp_broadleaf;
}

// Survival probability for tropical trees following van Nieustadt 2005
double survival_probability_tropics(double diameter_at_breast_height, double fire_line_intensity) {

	// DBH in cm
	diameter_at_breast_height *= CM_PER_M;

	double survival_probability = 1.;
	// compute surv. prob. at 3000kW/m first
	double survival_probability_3000 = 1. - max( 0.82 - 0.035 * pow(diameter_at_breast_height,0.7) , 0.);
	if ( fire_line_intensity > 7000. ) {
		double scal_fac = 1. - log((fire_line_intensity/7000.)) ;
		survival_probability = scal_fac * survival_probability_3000;
	}
	else if ( fire_line_intensity > 3000. ) {
		survival_probability =  survival_probability_3000;
	}
	else {
		survival_probability = exp(fire_line_intensity/3000. * log(survival_probability_3000));
	}

	survival_probability   = max(min(1.,survival_probability), 0.001);
    	
	return survival_probability;
}

// Survival probability for savannas following Bond 2008
double survival_probability_savanna(double height, double fire_line_intensity) {

	double intensity = fire_line_intensity / 1000. ;
	double survival_probability = max(0.,1. - ( 1./(1. + exp(1.5*(height - 0.5 * intensity - 1. )))));

	survival_probability = min (1.,survival_probability);
	
	return survival_probability;
}

// Survival probability for australian savanna-resprouters Cook 2005
double survival_probability_sprouter_savanna(double height, double fire_line_intensity) {

	// Height of max survival probability [m]
	// Taller trees are vulnerable due to age
	double const MAX_PROB_HEIGHT = 8.5;

	// Fire-line intensity [MW/m]
	double intensity = fire_line_intensity / 1000.; // Conversion to MW/m 

	// Minimum height for trees to survive [m]
	double height_min = 3.7 * (1.-exp(-0.19 * intensity));

	// Survival probability [fract.]
	double survival_probability;

	// Empirically generated functions by Vanessa Haverd
	// based on observations from G. Cook.
	if (height > MAX_PROB_HEIGHT && height > height_min) {
		survival_probability = ( -.0011 * intensity - .00002) * height
			+ .0075 * intensity + 1. ;
	}
	else if (height > height_min) {
		survival_probability = ( .0178 * intensity + .0144) * height
			+ ( -.1174 * intensity + 0.9158 );
	}
	else {
		survival_probability = 0.001;
	}
	survival_probability = max(1.e-3,min(1.,survival_probability));
	return survival_probability;
}

/* Compute Individual/Cohort survival probability
 * Depending on biome and geolocation the appropriate survival_probabilities
 * will be selected.
 */
double survival_probability(Patch& patch, Individual& indiv) {
 

	const Gridcell& gridcell   = patch.stand.get_gridcell();

	double height              = indiv.height;
	double fire_line_intensity = patch.fire_line_intensity;
	double lat                 = gridcell.get_lat();

	double survival_probability = 1.;

	// Allometry function for diameter-at-breast-height as used in growth.cpp
	double dbh = pow(height * 100. / indiv.pft.k_allom2, 1.0 / indiv.pft.k_allom3) / 100.;

	int biome = gridcell.simfire_biome;

	if ( vegmode == POPULATION ) {

		// Temperate Needleleaf
		if ( biome == SF_NEEDLELEAF) { 
			survival_probability = survival_probability_temp_needleleaf(dbh, fire_line_intensity, patch.soil.sompool[SURFCWD].cmass);
		}
		// Broadleaf and mixed
		else if ( biome ==  SF_BROADLEAF || biome == SF_MIXED_FOREST ) {

			if  (lat > -30 && lat < 30 ) {
				// Tropical
				survival_probability = survival_probability_tropics(dbh,fire_line_intensity);
			} else {
				// Temperate 
				survival_probability = survival_probability_temp_broadleaf(dbh, fire_line_intensity, 0);

			}
		}
		// Savanna, shrubland, and sparsely vegetated
		else if ( biome == SF_SHRUBS || biome == SF_SAVANNA || biome == SF_BARREN) {
			survival_probability = survival_probability_savanna(height, fire_line_intensity);
		}
		// Tundra
		else if ( biome == SF_TUNDRA ) {
			survival_probability = survival_probability_boreal(fire_line_intensity);
		} 
		else {
			return 1.;
		}
	}
	else if ( vegmode == COHORT || vegmode == INDIVIDUAL ) {

		// Needleleaf
		if ( indiv.pft.leafphysiognomy == NEEDLELEAF ) {
			// Temperate Needleleaf
			if ( fabs(lat) < 50.) {
				double mass_cwd = patch.soil.sompool[SURFCWD].cmass   ;   
				survival_probability = survival_probability_temp_needleleaf(dbh, fire_line_intensity, mass_cwd);
			}
			// Tundra
			else {
				survival_probability = survival_probability_boreal(fire_line_intensity);
			}
		}
		
		// Broadleaf
		else if ( indiv.pft.leafphysiognomy == BROADLEAF ) {
			// Broadleaf, mixed Forest, and majorly needle-leaf biomes
			if ( biome == SF_CROP || biome == SF_NEEDLELEAF || biome == SF_BROADLEAF || biome == SF_MIXED_FOREST || biome == SF_TUNDRA ) {
				if  (lat > -30 && lat < 30 ) {
					// Tropical 
					survival_probability = survival_probability_tropics(dbh,fire_line_intensity);
				} else {
					// Temperate 
					survival_probability = survival_probability_temp_broadleaf(dbh, fire_line_intensity, 0);
				}
			}
			// Savanna, shrubland, and sparsely vegetated
			else if ( biome == SF_SHRUBS || biome == SF_SAVANNA || biome == SF_BARREN ) {
				survival_probability = survival_probability_savanna(height, fire_line_intensity);
			}
			else {
				dprintf("Biome %d not found in BLAZE\n",biome);
				survival_probability = 1.;
			}
		} 
		else {
			fail("Physiognomy unknown...");
		}

	}

	else {
		fail("BLAZE: case not valid");
	}

	survival_probability = min(1. ,max(survival_probability,0.0001));

	return survival_probability;
}

/// Blaze: Compute wildfire combustion and fluxes thereof
/** The combustion part of the model. Here, all fire related fluxes
 * are computed and the changes applied to the affected pools.
 * This routine handles all current available
 * settings for vegetation model and uses the century som-model
 * soil-pools.
 */
bool blaze(Patch& patch, Climate& climate) {

	Gridcell& gridcell = climate.gridcell;
	
	// Flammable area witin gridcell
	double flammable_area = gridcell.landcover.frac[PASTURE]+gridcell.landcover.frac[NATURAL];

	// Return if total burnable area too small
	if (flammable_area < 0.00001 ) {
		return false;
	}
	
	// Effective area burned as fraction of burnable area
	double area_burned = gridcell.burned_area / flammable_area;

	// Correction fractions burned earlier in the same year (vegmode = POPULATION only)
	double accumulated_fraction_burned= 1.  / (1. - gridcell.annual_burned_area / flammable_area);

	// Check whether it burns
	if (!( randfrac(patch.random_seed(STOCHASTIC_FIRE_OCCURRENCE)) <= area_burned || vegmode == POPULATION)) {
		return false;
	}
	
	// Get relative fluxes between pools
	int fli_index = get_fire_line_intensity_index(patch.fire_line_intensity);

	// A NEGATIVE INDEX MEANS NO FIRE, AND IT USED TO REACH THE LOOKUP.
	// get_fireline_intensity() sets fire_line_intensity to -1 when available
	// fuel is below MIN_FUEL -- the "not enough fuel to ignite a fire" case --
	// and get_fire_line_intensity_index() then falls through every branch and
	// returns -1. Nothing here guarded it: the only two early returns above are
	// on flammable_area and the stochastic burn draw, and neither consults
	// fuel, so a patch drawn to burn with too little fuel to carry a fire
	// reached TURNOVERFRACT[r][-1]. That array is contiguous, so [r][-1] read
	// [r-1][4] -- the >=7000 kW/m column of the row above, the largest value in
	// almost every row -- and [0][-1] read outside the array altogether, which
	// is undefined behaviour rather than merely a wrong number. The effect ran
	// the wrong way: the patch the fuel test says cannot burn combusted all of
	// its surface litter and half its leaves.
	//
	// Returning false rather than clamping the index is deliberate. -1 is the
	// signal get_fireline_intensity() already sets to mean "no fire here", so
	// the repair is to honour it, not to substitute the lowest intensity class
	// and burn the patch a little. world-voxz,
	// biosphere/notes/fire-model-audit.md finding 11.
	if (fli_index < 0) {
		return false;
	}

	// Flag patch as burned for rest of the year. Set AFTER the fuel test, not
	// before it: the flag is read for the rest of the year, so a patch that
	// turns out to have no fire must never carry it even briefly.
	patch.burned = true;

	// Adjustment factor for fluxes
	double fab = 1.0;
	if ( vegmode == POPULATION ) {
		fab = max(area_burned * accumulated_fraction_burned,1.);
	}
		
	get_combustion_rates(patch,fli_index,gridcell.k_tun_litter);

	// Compute fluxes from soil litter pools to atmosphere first
	// as they are patch-specific. The fluxes to 
	// soil litter will be added in loop over individuals below.
	double cmtb2atm = fab * patch.litf_to_atm * patch.soil.sompool[SURFMETA].cmass  ;
	double cstr2atm = fab * patch.litf_to_atm * patch.soil.sompool[SURFSTRUCT].cmass;
	double cfwd2atm = fab * patch.lfwd_to_atm * patch.soil.sompool[SURFFWD].cmass   ;
	double ccwd2atm = fab * patch.lcwd_to_atm * patch.soil.sompool[SURFCWD].cmass   ;

	// Nitrogen proportional to cmass flux [kg(C)/m2]
	double nmtb2atm = fab * patch.litf_to_atm * patch.soil.sompool[SURFMETA].nmass  ;
	double nstr2atm = fab * patch.litf_to_atm * patch.soil.sompool[SURFSTRUCT].nmass;
	double nfwd2atm = fab * patch.lfwd_to_atm * patch.soil.sompool[SURFFWD].nmass   ;
	double ncwd2atm = fab * patch.lcwd_to_atm * patch.soil.sompool[SURFCWD].nmass   ;

	// Update soil-surface-litter pools
	// Carbon
	patch.soil.sompool[SURFMETA].cmass   -= cmtb2atm;
	patch.soil.sompool[SURFSTRUCT].cmass -= cstr2atm;
	patch.soil.sompool[SURFFWD].cmass    -= cfwd2atm;
	patch.soil.sompool[SURFCWD].cmass    -= ccwd2atm;
	// Nitrogen 
	patch.soil.sompool[SURFMETA].nmass   -= nmtb2atm;
	patch.soil.sompool[SURFSTRUCT].nmass -= nstr2atm;
	patch.soil.sompool[SURFFWD].nmass    -= nfwd2atm;
	patch.soil.sompool[SURFCWD].nmass    -= ncwd2atm;

	// Report C litter -> atmosphere fluxes
	patch.fluxes.report_flux(Fluxes::FIREC, cmtb2atm + cstr2atm + cfwd2atm + ccwd2atm);

	// Report N litter -> atmosphere fluxes
	report_fire_flux_n(patch, nmtb2atm + nstr2atm + nfwd2atm + ncwd2atm );
       
	Vegetation& vegetation=patch.vegetation;
       
	// Kill trees, grass, cohorts, and population fractions
	bool killed = false;
	double frac_survive;
	if ( vegmode == POPULATION ) {

		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();

			// For this individual ...
			killed=false;

			// Apply fire by burned area
			indiv.blaze_reduce_biomass(patch, (1.-fab));

			// Remove this cohort completely if all individuals killed
			// (in individual mode: removes individual if killed)
			if (negligible(indiv.densindiv)) {
				vegetation.killobj();
				killed=true;
			}
			if (!killed) {
				vegetation.nextobj(); // ... on to next individual
			}
		}
		fab = area_burned * accumulated_fraction_burned;
	}
	else {
		/* vegmode == Individual or Cohort
		 * taken from mortality_guess (vegdynam.cpp)
		 * Impose fire in this patch with probability=Area Burned as given in the top
		 * of this routine.
		 */
		
		// Loop through individuals
		vegetation.firstobj();
		while (vegetation.isobj) {
			Individual& indiv=vegetation.getobj();
			if ( !indiv.alive ) {
				vegetation.nextobj();
				continue;
			}

			// For this individual ...
			killed=false;

			if (indiv.pft.lifeform==GRASS) {
				// Reduce individual live biomass and freshly created litter
				indiv.reduce_biomass(MAX_GRASS_BURN,MAX_GRASS_BURN);
				
				// Remove NPP and put it to fire flux
				patch.fluxes.report_flux(Fluxes::FIREC,indiv.anpp*MAX_GRASS_BURN);
				indiv.anpp *= (1. - MAX_GRASS_BURN);

				// Kill object if burn is total
				if ( MAX_GRASS_BURN == 1.0 ) {
					indiv.kill();
					vegetation.killobj();
					killed=true;
				} 
				else {
					// Update allometry
					allometry(indiv);
				}
			}
			else {
				// TREE PFT

				if (ifstochmort) {
					/* Impose stochastic mortality.
					 * Each individual in cohort dies with probability 'mort_fire'
					 * Number of individuals represented by 'indiv'
					 * (round up to be on the safe side).
					 * In individual mode densindiv is 1!
					 */
					int nindiv=(int)(indiv.densindiv*patcharea+0.5);
					int nindiv_prev=nindiv;
					for (int i=0;i<nindiv_prev;i++) {
						if (randfrac(patch.random_seed(STOCHASTIC_FIRE_MORTALITY)) > survival_probability(patch, indiv)) {
							nindiv--;
						}
					}
					
					if (nindiv_prev) {
						frac_survive = (double)nindiv / (double)nindiv_prev;
					}
					else {
						frac_survive = 0.0;
					}
				}
				else { 	
					// Deterministic mortality (cohort mode only)
					frac_survive = survival_probability(patch, indiv);
				}

				// Reduce individual biomass on patch area basis
				// to account for loss of killed individuals
				indiv.blaze_reduce_biomass(patch,frac_survive);
				
				// Remove this cohort completely if all individuals killed
				// (in individual mode: removes individual if killed)
				if (negligible(indiv.densindiv)) {
					indiv.kill();
					vegetation.killobj();
					killed=true;
				} 
				else {
					// Update allometry 
					allometry(indiv);
				}
			}
			if (!killed) vegetation.nextobj(); // ... on to next individual
		}
	}
		
	// Transitional litter (last year's litter still to fall) removed by area burned (ab).
	// Total fluxes out of patch
	double cmtb2atm_t = 0.; 
	double cstr2atm_t = 0.; 
	double cfwd2atm_t = 0.; 
	double ccwd2atm_t = 0.; 
	double nmtb2atm_t = 0.; 
	double nstr2atm_t = 0.; 
	double nfwd2atm_t = 0.;
	double ncwd2atm_t = 0.;
	
	patch.pft.firstobj();
	while (patch.pft.isobj) {
		Patchpft& patchpft = patch.pft.getobj();
		
		double lton = lignin_to_n_ratio(patchpft.cmass_litter_leaf, patchpft.nmass_litter_leaf,
						LIGCFRAC_LEAF, patchpft.pft.cton_leaf_avr);
		double fm_leaf = metabolic_litter_fraction(lton); 
		
		// Carbon transitional litter fluxes
		double cmtb2atm = fab * patch.litf_to_atm * patchpft.cmass_litter_leaf * fm_leaf;
		double cstr2atm = fab * patch.litf_to_atm * patchpft.cmass_litter_leaf * (1.-fm_leaf);
		double cfwd2atm = fab * patch.lfwd_to_atm * patchpft.cmass_litter_sap;
		double ccwd2atm = fab * patch.lcwd_to_atm * patchpft.cmass_litter_heart;

		// Nitrogen transitional litter fluxes
		double nmtb2atm = fab * patch.litf_to_atm * patchpft.nmass_litter_leaf * fm_leaf;
		double nstr2atm = fab * patch.litf_to_atm * patchpft.nmass_litter_leaf * (1.-fm_leaf);
		double nfwd2atm = fab * patch.lfwd_to_atm * patchpft.nmass_litter_sap;
		double ncwd2atm = fab * patch.lcwd_to_atm * patchpft.nmass_litter_heart;
		
		// Update transitional rest-of-year litter pools
		// Carbon
		patchpft.cmass_litter_leaf             -= (cmtb2atm + cstr2atm);
		patchpft.cmass_litter_sap              -= cfwd2atm;
		patchpft.cmass_litter_heart            -= ccwd2atm;
		// Nitrogen
		patchpft.nmass_litter_leaf       -= (nmtb2atm + nstr2atm);
		patchpft.nmass_litter_sap        -= nfwd2atm; 
		patchpft.nmass_litter_heart      -= ncwd2atm; 
		
		// Calculate total fluxes to atmosphere within patch
		// Carbon
		cmtb2atm_t += cmtb2atm; 
		cstr2atm_t += cstr2atm; 
		cfwd2atm_t += cfwd2atm; 
		ccwd2atm_t += ccwd2atm; 
		// Nitrogen
		nmtb2atm_t += nmtb2atm; 
		nstr2atm_t += nstr2atm; 
		nfwd2atm_t += nfwd2atm;
		ncwd2atm_t += ncwd2atm;
		
		patch.pft.nextobj();
	}
	
	// Report C litter -> atmosphere flux from transitional pools
	patch.fluxes.report_flux(Fluxes::FIREC, cmtb2atm_t + cstr2atm_t + cfwd2atm_t + ccwd2atm_t);

	// Report N litter -> atmosphere flux from transitional pools
	report_fire_flux_n(patch, nmtb2atm_t + nstr2atm_t + nfwd2atm_t + ncwd2atm_t );
	
	return true;
	
}  

/// Update C/N - Pools due to fire
/** Applies the fluxes computed in blaze on the class::Individual
 * level affecting the live pools, transitional
 * litter pools, and influx to CENTURY litter pools.
 * In INDIVIDUAL and COHORT mode the actual biomass killed in
 * the routine "blaze" is used to compute the fraction of the live vegetation pools
 * while in POPULATION mode the fraction equal to burned area is used.
 * Input frac_survive means fraction of surviving INDIVIDUAL/COHORT
 * in respective mode or will be (1-burned area) in case of POPULATION
 * mode.
 */
void Individual::blaze_reduce_biomass(Patch& patch, double frac_survive) {

	// When a fire doesn't provide enough heat to burn a tree
	// and it still stochastically dies, use these flux parameters
	double const DEFAULT_WOOD_TO_ATM = 0.10;
	double const DEFAULT_WOOD_TO_STR = 0.03;
	double const DEFAULT_WOOD_TO_FWD = 0.07;
	double const DEFAULT_WOOD_TO_CWD = 0.80;	
	double const DEFAULT_LEAF_TO_ATM = 0.75;
	double const DEFAULT_LEAF_TO_LIT = 0.25;

	double frac_killed = 1. - frac_survive;

	if ( negligible(frac_killed) ) return;

	// local copies of live fluxes
	double wood_to_atm = patch.wood_to_atm;
	double wood_to_str = patch.wood_to_str;
	double wood_to_fwd = patch.wood_to_fwd;
	double wood_to_cwd = patch.wood_to_cwd;
	double leaf_to_atm = patch.leaf_to_atm;
	double leaf_to_lit = patch.leaf_to_lit;

	double fab = 1.0;

	if ( vegmode == INDIVIDUAL || vegmode == COHORT ) {
		double wtotw = patch.wood_to_atm + patch.wood_to_str + patch.wood_to_fwd + patch.wood_to_cwd ;
		// Adjust relative fluxes from wood when stochastic killing has occured
		if ( wtotw > 0.0 ) {
			wood_to_atm = patch.wood_to_atm * frac_killed;
			wood_to_str = (1. - patch.wood_to_atm ) * F_BARK               * frac_killed;
			wood_to_fwd = (1. - patch.wood_to_atm ) * F_BRANCH             * frac_killed;
			wood_to_cwd = (1. - patch.wood_to_atm ) * (1.-F_BARK-F_BRANCH) * frac_killed;
		}
		else {
			wood_to_atm = frac_killed * DEFAULT_WOOD_TO_ATM;
			wood_to_str = frac_killed * DEFAULT_WOOD_TO_STR;
			wood_to_fwd = frac_killed * DEFAULT_WOOD_TO_FWD;
			wood_to_cwd = frac_killed * DEFAULT_WOOD_TO_CWD;
		}
		// Adjust relative fluxes from leaves
		double ltotw = patch.leaf_to_atm + patch.leaf_to_lit;
		if ( ltotw > 0.0 ) {
			leaf_to_atm = patch.leaf_to_atm * frac_killed ;
			leaf_to_lit = (1. - patch.leaf_to_atm) * frac_killed ;
		}			               
		else {			               
			leaf_to_atm = frac_killed * DEFAULT_LEAF_TO_ATM;
			leaf_to_lit = frac_killed * DEFAULT_LEAF_TO_LIT;
		}
		fab = 1.0;
	} 
	else if ( vegmode == POPULATION ) {
		fab = frac_killed;
		fail("BLAZE does not run in POPULATION-mode");
	}

	// Compute mass per area fluxes from live pools

	Patchpft& ppft = patchpft();
	double lton = lignin_to_n_ratio(ppft.cmass_litter_leaf, ppft.nmass_litter_leaf, LIGCFRAC_LEAF,
					ppft.pft.cton_leaf_avr);


	// Metabolic litter fraction for litter
	double fm_leaf = metabolic_litter_fraction(lton);

	// Leaves
	double cleaf2atm = fab * leaf_to_atm * cmass_leaf_today();
	double cleaf2met = fab * leaf_to_lit * cmass_leaf_today() * fm_leaf;
	double cleaf2str = fab * leaf_to_lit * cmass_leaf_today() * (1. - fm_leaf);
	double nleaf2atm = fab * leaf_to_atm * nmass_leaf;
	double nleaf2met = fab * leaf_to_lit * nmass_leaf * fm_leaf;
	double nleaf2str = fab * leaf_to_lit * nmass_leaf * (1. - fm_leaf);

	// Sap-wood
	double csapw2atm = fab * wood_to_atm * cmass_sap ; 
	double csapw2str = fab * wood_to_str * cmass_sap ;
	double csapw2fwd = fab * (wood_to_fwd + wood_to_cwd) * cmass_sap ;
	double nsapw2atm = fab * wood_to_atm * nmass_sap ; 
	double nsapw2str = fab * wood_to_str * nmass_sap ;
	double nsapw2fwd = fab * (wood_to_fwd + wood_to_cwd) * nmass_sap ;

	// Heart-wood
	double chrtw2atm = fab * wood_to_atm * cmass_heart; 
	double chrtw2str = fab * wood_to_str * cmass_heart;
	double chrtw2cwd = fab * (wood_to_fwd + wood_to_cwd) * cmass_heart;
	double nhrtw2atm = fab * wood_to_atm * nmass_heart; 
	double nhrtw2str = fab * wood_to_str * nmass_heart;
	double nhrtw2cwd = fab * (wood_to_fwd + wood_to_cwd) * nmass_heart;

	// Roots
	// Assume the same percentage of root biomass killed as for total
	// above ground woody biomass.

	double lossratio = 0.;
	if ( cmass_sap + cmass_heart > 0. ) {
		lossratio = (csapw2atm + csapw2str + csapw2fwd +
			     chrtw2atm + chrtw2str + chrtw2cwd ) / (cmass_sap + cmass_heart);
	}

	// Root litter lignin:N ratio
	lton = lignin_to_n_ratio(ppft.cmass_litter_root, ppft.nmass_litter_root, LIGCFRAC_ROOT,
					ppft.pft.cton_root_avr);

	// Metabolic litter fraction for root 
	double fm_root = metabolic_litter_fraction(lton);
	
	double croot2met = fm_root        * lossratio * cmass_root;
	double croot2str = (1. - fm_root) * lossratio * cmass_root;
	double nroot2met = fm_root        * lossratio * nmass_root;
	double nroot2str = (1. - fm_root) * lossratio * nmass_root;

	// Update pools
	
	// Live carbon
	cmass_leaf      -= (cleaf2atm + cleaf2met + cleaf2str);
	cmass_sap       -= (csapw2atm + csapw2str + csapw2fwd);
	cmass_heart     -= (chrtw2atm + chrtw2str + chrtw2cwd); 
	cmass_root      -= (croot2met + croot2str); 
	
	// Deal with c-debt before asserting to litter & atm pool
	double saploss = csapw2atm + csapw2str + csapw2fwd;
	double hrtloss = chrtw2atm + chrtw2str + chrtw2cwd;

	if ( cmass_debt > 0. ) { 
		if ( cmass_debt <= saploss + hrtloss ) {
			if ( cmass_debt <= hrtloss ) {
				double scalefac = (hrtloss - cmass_debt) / hrtloss;
				chrtw2str  *= scalefac;
				chrtw2cwd  *= scalefac; 
				chrtw2atm  *= scalefac;
				cmass_debt  = 0.;
			}
			else {
				chrtw2str  = 0.;
				chrtw2cwd  = 0.; 
				chrtw2atm  = 0.;
				double scalefac = (saploss - (cmass_debt - hrtloss)) / saploss;
				csapw2str  *= scalefac;
				csapw2fwd  *= scalefac; 
				csapw2atm  *= scalefac;
				cmass_debt  = 0.;
			}
		}
		else {
			chrtw2str  = 0.;
			chrtw2cwd  = 0.;
			chrtw2atm  = 0.;
			csapw2str  = 0.;
			csapw2fwd  = 0.;
			csapw2atm  = 0.;
			cmass_debt -= (saploss + hrtloss);
		}
	}

	// Live nitrogen
	nmass_leaf      -= (nleaf2atm + nleaf2met + nleaf2str);
	nmass_sap       -= (nsapw2atm + nsapw2str + nsapw2fwd);
	nmass_heart     -= (nhrtw2atm + nhrtw2str + nhrtw2cwd);
	nmass_root      -= (nroot2met + nroot2str); 
	
       	double d_nstore = (nstore_longterm + nstore_labile) * (1. - frac_survive);
	nstore_longterm *= frac_survive; 
	nstore_labile   *= frac_survive; 

	// Deal with accumulated NPP 
	double loss_anpp = anpp * (1. - frac_survive);
	double loss_tot  = cleaf2atm + csapw2atm + chrtw2atm + cleaf2met + 
		cleaf2str + csapw2str + chrtw2str + csapw2fwd + chrtw2cwd +
		croot2met + croot2str;
	double anpp2atm  =  0.;
	double anpp2met  =  0.;
	double anpp2str  =  0.;
	double anpp2fwd  =  0.;
	double anpp2cwd  =  0.;
	double anpp2smtb =  0.;
	double anpp2sstr =  0.;
	if ( loss_tot > 0. ) {
		anpp2atm  = loss_anpp * ( cleaf2atm + csapw2atm + chrtw2atm ) / loss_tot;
		anpp2met  = loss_anpp *                           cleaf2met   / loss_tot;
		anpp2str  = loss_anpp * ( cleaf2str + csapw2str + chrtw2str ) / loss_tot;
		anpp2fwd  = loss_anpp *                           csapw2fwd   / loss_tot;
		anpp2cwd  = loss_anpp *                           chrtw2cwd   / loss_tot;
		anpp2smtb = loss_anpp *                           croot2met   / loss_tot;
		anpp2sstr = loss_anpp *                           croot2str   / loss_tot;
		anpp     -= loss_anpp;
	}

	// Report C live -> atm flux 
	patch.fluxes.report_flux(Fluxes::FIREC, cleaf2atm + csapw2atm + chrtw2atm + anpp2atm);

	// Report N live -> atm flux 
	report_fire_flux_n(patch, nleaf2atm + nsapw2atm + nhrtw2atm + d_nstore);
 
	// Soil-surface-litter carbon
	patch.soil.sompool[SURFMETA].cmass   += anpp2met  + cleaf2met;
	patch.soil.sompool[SURFSTRUCT].cmass += anpp2str  + cleaf2str + csapw2str + chrtw2str;
	patch.soil.sompool[SURFFWD].cmass    += anpp2fwd  + csapw2fwd;
	patch.soil.sompool[SURFCWD].cmass    += anpp2cwd  + chrtw2cwd;

	// Deep soil litter carbon
	patch.soil.sompool[SOILMETA].cmass   += anpp2smtb + croot2met;
	patch.soil.sompool[SOILSTRUCT].cmass += anpp2sstr + croot2str;

	// Soil-surface-litter nitrogen 
	patch.soil.sompool[SURFMETA].nmass   += nleaf2met;
	patch.soil.sompool[SURFSTRUCT].nmass += nleaf2str + nsapw2str + nhrtw2str;
	patch.soil.sompool[SURFFWD].nmass    += nsapw2fwd;
	patch.soil.sompool[SURFCWD].nmass    += nhrtw2cwd;

	// Deep soil litter nitrogen
	patch.soil.sompool[SOILMETA].nmass   += nroot2met;
	patch.soil.sompool[SOILSTRUCT].nmass += nroot2str;

	if (pft.lifeform != GRASS) {
		densindiv *= frac_survive;

		if ( negligible(densindiv) && cmass_debt > 0.) {
			// Excess c_debt shifted from RA to NPP
			report_flux(Fluxes::NPP, cmass_debt);
			report_flux(Fluxes::RA, -cmass_debt);
			cmass_debt = 0.;
		}
	}
}

// Daily accounting of blaze relevant parameters
/* Accounting of long-term averages needed for BLAZE as well as
 * the computation of daily burned area and fire-specific parameteers like
 * the Keetch-Byram Drought-index and Forest Fire Danger Index (FFDI).
 */
void blaze_accounting_gridcell(Climate& climate) {

	Gridcell& gridcell = climate.gridcell;

	// To keep track of burned area over the year
	// reset accumulated area_burned to 0 on begining of year
	if (date.day == 0 ) {
		gridcell.annual_burned_area         = 0.0;
		gridcell.simfire_annual_burned_area = 0.0;
		climate.rainfall_cur                = 0.0;

		for (int i = 0; i < 12; i++) {
			gridcell.monthly_burned_area[i] = 0.0;
		}

		double lat = climate.gridcell.get_lat();

		// Latitude depending tuning values mortality
		if ( fabs(lat) >= 50.) {
			gridcell.k_tun_litter = K_LITTER_BOREAL;
		}
		else if ( fabs(lat) >= 30. && fabs(lat) < 50.) {
			gridcell.k_tun_litter = K_LITTER_TEMPERATE;
		}
		else if ( gridcell.simfire_biome == SF_SAVANNA ) {
			gridcell.k_tun_litter = K_LITTER_SAVANNA;
		}
		else {
			gridcell.k_tun_litter = K_LITTER_TROPICS;
		}
	}
	// Reset patch variables
	Gridcell::iterator gc_itr = gridcell.begin();
	while (gc_itr != gridcell.end()) {
		Stand& stand = *gc_itr;
		stand.firstobj();
		while (stand.isobj) {
			Patch& patch = stand.getobj();
			patch.fire_line_intensity = 0.;
			if (date.day==0) {
				patch.burned = false;
			}
			stand.nextobj();
		}
		++gc_itr;
	}
	
	// Keep track of Days-since-last-rainfall and accumulated last rainfall
	if (climate.prec > 0.01) {
		if (climate.days_since_last_rainfall > 0) {
			climate.last_rainfall = climate.prec;
		}
		else {
			climate.last_rainfall += climate.prec;
		}
		climate.days_since_last_rainfall = 0;
	}
	else {
		climate.days_since_last_rainfall++;
	}

	climate.rainfall_cur += climate.prec;

	// Update the Keetch-Byram-Drought-Index (Keetch et al. 1968)
	double v  = climate.u10    * KMH_PER_MS;       // Wind speed at 10m height [km/h] (for FFDI)
	double rh = climate.relhum * FRACT_TO_PERCENT; // relative humidity [%] (for FFDI)
	double t  = climate.tmax  ;                    // day's max temperature [deg C] (for KBDI) 

	// Gust parameterisation
	v = ( 214.7 * pow(  v + 10. ,-1.6968)  + 1. ) * v;

	double dkbdi; // change in Keetch-Byram-Drought-Index due to rainfall history
	if (climate.days_since_last_rainfall == 0) {
		if (climate.last_rainfall > 5.) {
			dkbdi = 5. - climate.last_rainfall;
		} 
		else {
			dkbdi = 0.0;
		}
	}
	else {
		// THE RAINFALL TERM IS A SITE CLIMATOLOGY IN INCHES PER EARTH YEAR.
		// Keetch and Byram (1968) use the site's MEAN ANNUAL RAINFALL as the
		// parameter setting how fast that site dries, and -.0441 was fitted
		// with it in inches per 365-day year. climate.rainfall_annual_avg is a
		// running mean of the model's own annual total, and this model's year
		// is 183 absolute days, so feeding it in directly describes every site
		// as receiving about half the rain it does.
		//
		// Because the term is exponential the resulting bias DEPENDS ON THE
		// CLIMATE rather than scaling the field -- 0.69 of the intended drought
		// accumulation at an Earth-equivalent 500 mm/yr, 0.53 at 1000 and 0.46
		// at 1500 -- so it compresses the range of simulated fire weather and
		// moves its spatial pattern. It also runs opposite to the guess: a
		// shorter year makes the world read as wetter and burn less.
		//
		// DECIDED: the argument is a climatology and not an integration, so the
		// model-year total is re-expressed per Earth year rather than the
		// coefficient being rescaled. Rescaling -.0441 would repair the unit of
		// one number inside a fitted index and leave the index reading a
		// quantity it was not fitted on. world-4iyz,
		// biosphere/notes/fire-model-audit.md finding 12.
		const double rainfall_per_earth_year =
			climate.rainfall_annual_avg / VESPER_EARTH_YEARS_PER_ORBIT;
		dkbdi = (( 800. - climate.kbdi) * (.968 * exp(.0486 * (t * 9./5. + 32.)) 
			 - 8.3) / 1000. / (1. + 10.88 * exp(-.0441 * 
			 rainfall_per_earth_year/25.4)) * .254);
	}
	climate.kbdi = max(0.0,climate.kbdi + dkbdi);

	// ...and McArthur-Drought-Factor D ... (Noble, 1980)
	double mcarthur_d = .191 * ( climate.kbdi + 104. ) * pow( climate.days_since_last_rainfall + 1.,1.5 ) / 
		( 3.52 * pow( climate.days_since_last_rainfall + 1. ,1.5 ) + climate.last_rainfall - 1. );
	mcarthur_d = max(0.0,min(10.0,mcarthur_d));
	
	// ... and finally: McArthur's Forest Fire Danger Index
	double mcarthur_fire_index = 2. * exp( -.45 + .987 * log(mcarthur_d+.001) - 
				     .03456 * rh + .0338 * t + .0234 * v );
	mcarthur_fire_index = max(0.0,mcarthur_fire_index);

	// Monthly ffdi max
	int dayx = date.day % 30;
	climate.ffdi_monthly[dayx] = mcarthur_fire_index;
	climate.mcarthur_forest_fire_index = 0.;
	for (int x=0; x<30;x++) {
		if ( climate.mcarthur_forest_fire_index < climate.ffdi_monthly[x] ) {
			climate.mcarthur_forest_fire_index = climate.ffdi_monthly[x];
		}
	}
	gridcell.effective_burned_area = 0.;

	// Get burned area from SIMFIRE
	gridcell.burned_area = simfire_burned_area(gridcell);
	
	// End of year clean-up
	if (date.islastday && date.islastmonth) {
		
		double weighting; // used to compute running average of ann rainfall
		
		// Update running mean of average annual rainfall		
		if (date.year < RAINFALL_AVERAGING_SPAN) {
			weighting = date.year + 1;
		}
		else {
			weighting = RAINFALL_AVERAGING_SPAN;
		}
		
		climate.rainfall_annual_avg = ((weighting - 1.) * climate.rainfall_annual_avg 
			+ climate.rainfall_cur ) / weighting;
		climate.rainfall_cur   = 0.0;

		// Rotate the fire-danger ring buffer so the new year's first days do
		// not land on slots still holding the old year's, out of order. The
		// buffer is written at date.day % AVERAGING_FFDI, so the shift is the
		// year's remainder modulo the window.
		//
		// This read `365 % AVERAGING_FFDI`, giving 25, on a model whose year is
		// 183 absolute days -- VESPER_YEAR_LENGTH_DAYS, which date.year_length()
		// returns. The right shift here is 30 - (183 % 30) = 27, so the rotation
		// was off by two slots at every year boundary. Taken from
		// date.year_length() rather than corrected to 27, so the buffer tracks
		// the declared year instead of carrying a second hardcoded one.
		// world-4iyz, biosphere/notes/fire-model-audit.md finding 12.
		//
		// The WINDOW is a separate question and is deliberately not touched.
		// This model's months are 15 days, so a "monthly" maximum taken over 30
		// spans two of them. Whether the window is an absolute-time memory of
		// drying, in which case 30 carries, or a fraction of the seasonal cycle,
		// in which case it does not, is FIRE-6's decision.
		const int AVERAGING_FFDI = 30;
		double tmp[AVERAGING_FFDI];
		int avg_shift = AVERAGING_FFDI - (date.year_length() % AVERAGING_FFDI);
		
		for (int i = 0; i < AVERAGING_FFDI; i++) {
			int idx = (i + avg_shift) % AVERAGING_FFDI;
			tmp[idx] = climate.ffdi_monthly[i];
		}
		
		for (int i = 0; i < AVERAGING_FFDI; i++) {
			climate.ffdi_monthly[i] = tmp[i];
		}
	}
}		     

// The driver routine for BLAZE
/* This is the driver routine for BLAZE. It retrieves potential Fire-Line-Intensity
 * and calls the blaze main routine patch-wise.
 */
void blaze_driver(Patch& patch, Climate& climate) {

	// Check whether BLAZE should be called at all
	// Has BLAZE been chosen as firemodel?
	if (firemodel != BLAZE) { 
		return;
	}

	// Do not burn before century soil has started
	if (date.year < patch.soil.solvesomcent_beginyr) {
		return;
	}
	// If fires are excluded from patch: Return
	if (!patch.has_fires()) {
		return;
	}

	// If patch has already burned this year
	if (patch.burned) {
		return;
	}

	Gridcell& gridcell = climate.gridcell;
	
	// Today's potential firelineintensity
	get_fireline_intensity(patch, climate);

	if (negligible(gridcell.burned_area) || patch.fire_line_intensity < 0.) {
		return;
	}
	
	// Call combustion model blaze
	bool caught_fire = blaze(patch, climate);
	
	// Now add Burned Area to output if patch cought fire, as the whole patch burns
	// we use the full area as burned area.
	if (caught_fire) {
		Stand& stand = patch.stand;
		double gridcell_fraction =  stand.get_gridcell_fraction() / (double)stand.npatch();
	
		gridcell.effective_burned_area           += gridcell_fraction;
		gridcell.annual_burned_area              += gridcell_fraction;
		gridcell.monthly_burned_area[date.month] += gridcell_fraction;
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// REFERENCES
//
// Dalziel, BD, Tree Mortality Following Boreal Forest Fires Reveals Scale-Dependant 
//   Interactions Between Community Structure and Fire Intensity, Ecos., 12, 2009
//   https://doi.org/10.1007/s10021-009-9272-2
// Kobziar, L, Tree mortality patterns following prescribed fires in a mixed conifer forest, 
//   Can. J. For. Res., 36, 2006, doi:10.1139/X06-183
// Hickler, T, USING A GENERALIZED VEGETATION MODEL TO SIMULATE VEGETATION DYNAMICS 
//   IN NORTHEASTERN USA, Ecology, 85, 2004, doi: 10.1890/02-0344 
// Bond, WJ, What Limits Trees in C4 Grasslands and Savannas?, 
//   Annu. Rev. Ecol. Evol. Syst. 2008. 39, doi:10.1146/annurev.ecolsys.39.110707.173411
// Cook, G, pers. comm., 2014, 
// van Nieuwstadt, MGL, Drought, fire and tree survival in a Borneo rain forest, 
//   East Kalimantan, Indonesia,J.o. Ecology, Vol.93, 1,  2005
//   https://doi.org/10.1111/j.1365-2745.2004.00954.x
// Liedloff, A, Predicting a ?tree change? in Australia?s tropical savannas: Combining different
//   types of models to understand complex ecosystem behaviour, Ecological Modelling 221, 2010
//   doi:10.1016/j.ecolmodel.2010.07.022
// Noble, IR, McArthur's fire-danger expressed as equations, Austr. J. Ecol. 5, 1980
//   doi:10.1111/j.1442-9993.1980.tb01243.x
// Keetch, JJ, A Drought Index for Forest Fire Control, Res. Pap. SE-38. Asheville, 
//   NC: U.S. Department of Agriculture, 1968
