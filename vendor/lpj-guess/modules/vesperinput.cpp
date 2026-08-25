///////////////////////////////////////////////////////////////////////////////////////
/// \file vesperinput.cpp
/// \brief Input module for Vesper, driven by an ExoPlaSim climatology
///
/// See vesperinput.h. Modelled on demoinput.cpp, which is the module the
/// LPJ-GUESS documentation points new users at, and departs from it in five
/// ways, each forced by this being another planet:
///
///  1. Insolation is supplied as NETSWRAD_TS, the net downward surface shortwave
///     ExoPlaSim already computes, rather than as percentage sunshine. The
///     climate model knows the surface albedo far better than driver.cpp's
///     global BETA constant of 0.17 does, and this project computes that albedo
///     from lithology.
///  2. Nitrogen deposition is declared per EARTH year and divided by the Earth
///     year, because deposition is an atmospheric flux in absolute time and
///     does not know this world's orbit. See
///     biosphere/notes/time-base-unit-contract.md.
///  3. The gridlist, soil codes and climate arrive in one generated binary file
///     rather than several curated text files.
///  4. The file's year length is checked against the compiled-in one.
///  5. The fourth climate array is the range of the SURFACE temperature, not
///     the near-surface air diurnal range, so it is NOT handed to climate.dtr
///     and ifbvoc is refused. The two are different variables and this world's
///     climate product carries only the first. See getclimate, and
///     biosphere/notes/ecological-forcing-field-contract.md.
///
///////////////////////////////////////////////////////////////////////////////////////

#include "config.h"
#include "vesperinput.h"
#include "soilinput.h"
#include "parallel.h"
#include "driver.h"
#include "parameters.h"
#include "guess.h"
#include <stdio.h>
#include <string.h>

REGISTER_INPUT_MODULE("vesper", VesperInput)

namespace {

/// Little-endian, and both writer and reader are x86-64. Checked via the magic.
const char DRIVER_MAGIC[8] = {'V','E','S','P','D','R','V','7'};

/// Bins per year in the driver file. These are the MODEL's months, the same
/// twelve VESPER_MONTH_LENGTHS sizes Date with, because interp_monthly_*
/// spreads each one over exactly ndaymonth days. build_lpj_driver.py remaps the
/// producer's own time bins onto them conservatively rather than assuming the
/// two partitions agree; they do not, and the last one differed by two days.
const int DRIVER_BINS = 12;

/// Length of the provenance string the driver file carries
const int PROVENANCE_BYTES = 64;

template<typename T>
void read_or_fail(FILE* in, T* target, size_t count, const char* what) {
	if (fread(target, sizeof(T), count, in) != count) {
		fail("vesperinput: driver file ended early while reading %s", what);
	}
}

} // namespace

VesperInput::VesperInput()
	: current(0), total_cells(0), nyear(1), co2(0.0), ndep(0.0),
	  have_soilmap(false), years(1), loaded_year(-1) {

	declare_parameter("nyear", &nyear, 1, 10000,
		"Number of simulation years to run after spinup");
	declare_parameter("file_driver", &file_driver, 300,
		"Path to the binary driver file built by build_lpj_driver.py");
	declare_parameter("file_soilmap", &file_soilmap, 300,
		"Optional soil map from pedology/scripts/build_soil.py. Falls back to "
		"the driver file's LPJ soil code when unset.");
}

VesperInput::~VesperInput() {
}

void VesperInput::read_driver() {

	FILE* in = fopen(file_driver, "rb");
	if (!in) {
		fail("vesperinput: could not open driver file %s", (char*)file_driver);
	}

	char magic[8];
	read_or_fail(in, magic, 8, "magic");
	if (memcmp(magic, DRIVER_MAGIC, 8) != 0) {
		fclose(in);
		char expected[9];
		memcpy(expected, DRIVER_MAGIC, 8);
		expected[8] = '\0';
		char found[9];
		memcpy(found, magic, 8);
		found[8] = '\0';
		fail("vesperinput: %s has magic '%s', expected '%s'. Rebuild it with "
		     "build_lpj_driver.py.", (char*)file_driver, found, expected);
	}

	int ncells = 0, nbins = 0, year_length = 0, nyears = 0;
	read_or_fail(in, &ncells, 1, "cell count");
	read_or_fail(in, &nbins, 1, "bin count");
	read_or_fail(in, &year_length, 1, "year length");
	read_or_fail(in, &nyears, 1, "years of climate");

	// This world's year is a function of its stellar flux. A driver file and a
	// binary built at different fluxes describe different planets, and the
	// symptom would be a slow seasonal drift rather than a crash, so refuse.
	if (year_length != Date::MAX_YEAR_LENGTH) {
		fclose(in);
		fail("vesperinput: driver file was built for a %d-day year but this "
		     "binary is compiled for %d. Regenerate framework/vesper.h with "
		     "build_vesper_header.py and rebuild, or rebuild the driver file.",
		     year_length, (int)Date::MAX_YEAR_LENGTH);
	}
	if (nbins != DRIVER_BINS) {
		fclose(in);
		fail("vesperinput: driver file has %d bins per year, expected %d",
		     nbins, DRIVER_BINS);
	}
	if (ncells < 1) {
		fclose(in);
		fail("vesperinput: driver file contains no cells");
	}
	if (nyears < 1) {
		fclose(in);
		fail("vesperinput: driver file declares %d years of climate", nyears);
	}
	years = nyears;

	read_or_fail(in, &co2, 1, "CO2");
	read_or_fail(in, &ndep, 1, "nitrogen deposition");

	char provenance_buffer[PROVENANCE_BYTES + 1];
	read_or_fail(in, provenance_buffer, PROVENANCE_BYTES, "provenance tag");
	provenance_buffer[PROVENANCE_BYTES] = '\0';
	provenance = provenance_buffer;

	cells.resize(ncells);
	for (int i = 0; i < ncells; i++) {
		Cell& cell = cells[i];
		int pad = 0;
		read_or_fail(in, &cell.lon, 1, "longitude");
		read_or_fail(in, &cell.lat, 1, "latitude");
		read_or_fail(in, &cell.soilcode, 1, "soil code");
		read_or_fail(in, &pad, 1, "padding");
		read_or_fail(in, &cell.soil_states.b, 1, "retention exponent");
		read_or_fail(in, &cell.soil_states.saturation, 1, "saturation");
		read_or_fail(in, &cell.soil_states.field_capacity, 1, "field capacity");
		read_or_fail(in, &cell.soil_states.wilting_point, 1, "wilting point");
		cell.layer_usable.resize(NSOILLAYER);
		read_or_fail(in, &cell.layer_usable[0], NSOILLAYER, "layer usable shares");
		// A cell whose driver was built without a soil map carries the sentinel
		// rather than a curve, and the soil-code path is what runs for it.
		cell.soil_states.declared = (cell.soil_states.saturation > 0.0);
		if (cell.soilcode < 0 || cell.soilcode > 9) {
			fclose(in);
			fail("vesperinput: cell %d has invalid LPJ soil code %d",
			     i, cell.soilcode);
		}
		const size_t span = (size_t)nbins * (size_t)nyears;
		cell.temp.resize(span);
		cell.prec.resize(span);
		cell.insol.resize(span);
		cell.tsrange.resize(span);
		read_or_fail(in, &cell.temp[0], span, "temperature");
		read_or_fail(in, &cell.prec[0], span, "precipitation");
		read_or_fail(in, &cell.insol[0], span, "insolation");
		read_or_fail(in, &cell.tsrange[0], span, "surface temperature range");
	}

	// Keep only this process's share of the cells.
	//
	// LPJ-GUESS's own parallel mode splits work by pre-splitting the gridlist
	// file, one per rank, which is what submit.sh does with awk and split. This
	// module has no gridlist: the whole planet arrives in one driver file, so
	// without this every rank would simulate every cell and the run would be
	// nprocs times slower than serial rather than faster.
	//
	// Strided rather than blocked, because cost per cell tracks climate and
	// climate tracks latitude, and the cells are emitted in latitude order. A
	// contiguous block would hand one rank the entire ice-free tropics and
	// another the polar desert.
	const int rank = GuessParallel::get_rank();
	const int nprocs = GuessParallel::get_num_processes();
	if (nprocs > 1) {
		std::vector<Cell> mine;
		for (size_t i = rank; i < cells.size(); i += nprocs) {
			mine.push_back(cells[i]);
		}
		total_cells = (int)cells.size();
		cells.swap(mine);
	}
	else {
		total_cells = (int)cells.size();
	}

	// A driver file with trailing content is a version mismatch we have not
	// noticed, so say so rather than silently using a prefix of it.
	char trailing;
	if (fread(&trailing, 1, 1, in) == 1) {
		fclose(in);
		fail("vesperinput: driver file has unexpected trailing data");
	}
	fclose(in);
}

void VesperInput::init() {

	read_driver();

	dprintf("Vesper driver: %s\n", (char*)file_driver);
	if (GuessParallel::get_num_processes() > 1) {
		dprintf("  %d of %d land cells on rank %d of %d, %d bins, %d-day year\n",
		        (int)cells.size(), total_cells, GuessParallel::get_rank(),
		        GuessParallel::get_num_processes(), DRIVER_BINS,
		        (int)Date::MAX_YEAR_LENGTH);
	}
	else {
		dprintf("  %d land cells, %d bins per year, %d-day year\n",
		        (int)cells.size(), DRIVER_BINS, (int)Date::MAX_YEAR_LENGTH);
	}
	if (years > 1) {
		dprintf("  %d years of climate, cycled; spin-up sees the whole cycle\n",
		        years);
	}
	else {
		dprintf("  1 year of climate, repeated: a fixed climate\n");
	}
	dprintf("  CO2 %g ppm, N deposition %g kgN/ha/yr\n", co2, ndep);
	dprintf("  built from %s\n\n", (char*)provenance);

	// A soil map from the pedology component supersedes the driver file's soil
	// code: texture weathered under this world's climate against a texture
	// inferred from parent rock alone. Which one was used is printed, because
	// the difference is not visible in the output otherwise.
	have_soilmap = (file_soilmap != "");
	if (have_soilmap) {
		soilinput.init(file_soilmap);
		dprintf("  soil from %s (pedology)\n\n", (char*)file_soilmap);
	}
	else {
		dprintf("  soil from the driver file's LPJ soil codes (parent material "
		        "only, no pedogenesis)\n\n");
	}

	// The one reader of climate.dtr is bvoc.cpp's daytime_temp, and this world's
	// forcing carries no near-surface air temperature range for it to read. See
	// getclimate below and biosphere/notes/ecological-forcing-field-contract.md.
	// Fail here rather than at the point of use: a BVOC run that starts is a run
	// whose leaf temperature came from somewhere.
	if (ifbvoc) {
		fail("vesperinput: ifbvoc 1 needs a near-surface AIR temperature range "
		     "and the Vesper forcing carries only the SURFACE temperature "
		     "range. ExoPlaSim computes the air extrema as output codes 201 and "
		     "202 and no product carries them yet. See "
		     "biosphere/notes/ecological-forcing-field-contract.md.");
	}

	landcover_input.init();
	management_input.init();

	tprogress.init();
	tmute.init();
	tprogress.settimer();
	tmute.settimer(MUTESEC);

	current = 0;
}

void VesperInput::interpolate(const Cell& cell, int year_index) {

	// The framework's own conserving interpolators, so bin means stay means and
	// bin totals stay totals. They read date.ndaymonth[], which the patched Date
	// fills with Vesper's months, so no bin length is assumed here.
	const size_t offset = (size_t)year_index * (size_t)DRIVER_BINS;
	interp_monthly_means_conserve(&cell.temp[offset], dtemp);
	interp_monthly_totals_conserve(&cell.prec[offset], dprec, 0.0);
	interp_monthly_means_conserve(&cell.insol[offset], dinsol, 0.0);
	interp_monthly_means_conserve(&cell.tsrange[offset], dtsrange, 0.0);
	loaded_year = year_index;
}

bool VesperInput::getgridcell(Gridcell& gridcell) {

	if (current >= cells.size()) {
		return false;
	}

	const Cell& cell = cells[current];

	dprintf("\nCommencing simulation for stand at (%g,%g), soil code %d\n\n",
	        cell.lon, cell.lat, cell.soilcode);

	gridcell.set_coordinates(cell.lon, cell.lat);

	// Net downward shortwave averaged over the whole timestep, which is what
	// ExoPlaSim's rss is. With this insolation type driver.cpp applies no albedo
	// correction of its own, which is the point: the surface albedo here is
	// computed from lithology and is already in the number.
	gridcell.climate.instype = NETSWRAD_TS;

	if (have_soilmap) {
		// The contract's states first: `get_soil` reads them and derives no
		// retention curve of its own. WORLD-OF6N.
		soilinput.set_contract_states(cell.soil_states);
		soilinput.get_soil(cell.lon, cell.lat, gridcell);
		apply_column_geometry(gridcell, cell.layer_usable);
	}
	else {
		// No soil map, so no contract states and no texture to derive them
		// from. The driver's single LPJ soil code is a texture class with its
		// own tabulated capacities, and it carries no column geometry either.
		soil_parameters(gridcell.soiltype, cell.soilcode);
	}

	// Year 0 of this cell, so day 0 has data before getclimate first runs.
	// getclimate re-interpolates at the start of every year after that.
	loaded_year = -1;
	interpolate(cell, 0);

	clear_all_graphs();

	return true;
}

void VesperInput::apply_column_geometry(Gridcell& gridcell,
                                        const std::vector<double>& layer_usable) {

	// LPJ-GUESS gives every gridcell the same physical profile and derives its
	// water capacity from texture alone. On a world where relief and
	// erodibility are known, that is a real omission: a plant on 0.2 m of
	// regolith over bedrock has a seventh of the water a plant on a deep
	// profile has, and no combination of sand and clay fractions can express
	// it.
	//
	// So each layer's capacity is scaled by the share of that layer the soil
	// column actually carries. THAT SHARE IS READ AND NOT DERIVED. It is the
	// land column property contract's weathered-bedrock vertical rule -- the
	// part of the layer above the regolith contact, plus the weathered-bedrock
	// water fraction of the part below it, capped so sub-bedrock material can
	// at most match the soil above -- evaluated once in
	// pedology/scripts/land_column_properties.py and carried per layer in the
	// driver file. The rule stood here and was transcribed there as well, which
	// is two copies of one rule; WORLD-OF6N left one.
	//
	// This scales capacity, not the layer geometry: the layers still exist and
	// still conduct heat, they simply hold less water. Rooting depth is
	// unchanged, which slightly overstates how much of the profile deep-rooted
	// PFTs can reach on thin soil. Recorded rather than corrected, because
	// changing rootdist would reach into PFT parameters that are Earth
	// calibrations.
	//
	// The ceiling on the bedrock share is a limitation of this model and not of
	// the world. Deeply weathered saprolite really does hold two to four times
	// what its soil does, but LPJ-GUESS ties a layer's saturation capacity to
	// its texture-derived porosity, and scaling `wsats` past that drives
	// `Frac_air` negative in Soil::update_soil_diffusivities, which is a hard
	// failure. Representing the upper half of the observed range needs a
	// genuinely separate bedrock layer with its own porosity, which is what
	// Lapides et al. (Biogeosciences 2024) added to this same model rather than
	// rescaling the existing profile.

	if ((int)layer_usable.size() != NSOILLAYER) {
		fail("vesperinput: the driver carries %d per-layer usable shares and "
		     "this binary has NSOILLAYER = %d. Rebuild the driver with "
		     "build_lpj_driver.py against this profile.",
		     (int)layer_usable.size(), (int)NSOILLAYER);
	}

	Soiltype& soil = gridcell.soiltype;

	const double upper_layer_mm = SOILDEPTH_UPPER / (double)NSOILLAYER_UPPER;
	const double lower_layer_mm = SOILDEPTH_LOWER
		/ (double)(NSOILLAYER - NSOILLAYER_UPPER);

	double kept_upper = 0.0;
	double kept_lower = 0.0;
	for (int layer = 0; layer < NSOILLAYER; layer++) {
		const double thickness = (layer < NSOILLAYER_UPPER)
			? upper_layer_mm : lower_layer_mm;
		const double usable = layer_usable[layer];

		// A hard floor is a numerical guard and not physics, and the contract
		// applies the same one where it computes the share. Checked rather than
		// re-imposed: LPJ-GUESS carries soil water as a fraction of each
		// layer's capacity, computing `wcont = Faw_layer / soiltype.awc[layer]`
		// in soilwater.cpp and canexch.cpp, so a layer at exactly zero divides
		// by zero and the NaN propagates through the nitrogen substrate and
		// zeroes the gridcell's vegetation silently, reporting zero LAI rather
		// than failing.
		if (!(usable > 0.0)) {
			fail("vesperinput: layer %d has a usable share of %g. A layer at "
			     "zero capacity divides by zero in soilwater.cpp and zeroes "
			     "the gridcell's vegetation without failing; the contract's "
			     "numerical floor should have prevented it.", layer, usable);
		}

		soil.awc[layer] *= usable;
		soil.wp[layer] *= usable;
		soil.wsats[layer] *= usable;

		if (layer < NSOILLAYER_UPPER) {
			kept_upper += usable * thickness;
		}
		else {
			kept_lower += usable * thickness;
		}
	}

	// The aggregate two-layer figures have to move with the per-layer ones or
	// the hydrology and the diagnostics disagree with each other.
	const double upper_scale = kept_upper / SOILDEPTH_UPPER;
	const double lower_scale = kept_lower / SOILDEPTH_LOWER;
	soil.gawc[0] *= upper_scale;
	soil.gawc[1] *= lower_scale;
	soil.gwp[0] *= upper_scale;
	soil.gwp[1] *= lower_scale;
	soil.gwsats[0] *= upper_scale;
	soil.gwsats[1] *= lower_scale;
	soil.wtot = soil.gawc[0] + soil.gawc[1] + soil.gwp[0] + soil.gwp[1];
}

void VesperInput::getlandcover(Gridcell& gridcell) {

	landcover_input.getlandcover(gridcell);
	landcover_input.get_land_transitions(gridcell);
}

bool VesperInput::getclimate(Gridcell& gridcell) {

	Climate& climate = gridcell.climate;
	const Cell& cell = cells[current];

	// Step through the years the driver carries, wrapping. One year repeats and
	// gives a fixed climate; several give interannual variability, which is what
	// a stellar cycle is. The wrap also means spin-up sees the whole cycle
	// rather than one arbitrary phase of it.
	const int wanted = (years > 1) ? (date.year % years) : 0;
	if (wanted != loaded_year) {
		// Whenever the loaded year is not the wanted one, and not only on day 0.
		// A run resuming from a state file at state_day >= 0 never sees day 0 of
		// the year it resumes into, and without this it would run the rest of
		// that year on whatever the arrays last held. The day-0 block below still
		// owns the seasonal landmarks and the progress report; this only owns the
		// forcing arrays, and on day 0 the two conditions coincide exactly as
		// they did when they were one test.
		interpolate(cell, wanted);
		if (date.day != 0) {
			// Day 0 sets the landmarks below. A mid-year resume has to set them
			// here or summergreen phenology spends the rest of the year on the
			// landmarks initdrivers starts from, which are not this cell's.
			climate.set_seasonal_cycle(dtemp);
		}
	}

	// Per day of ABSOLUTE time, so the Earth year and not this world's. The
	// driver file declares kgN/ha per Earth year; dividing by the simulation
	// year would deliver a whole Earth year of nitrogen every orbit, which on
	// this calendar is close to twice the declared rate. Split evenly between
	// reduced and oxidised, as demoinput does; neither the split nor the total
	// is measured on this world, and both are declared in the driver file
	// rather than assumed here.
	const double per_day = ndep / 2.0 / VESPER_EARTH_YEAR_DAYS * HA_PER_M2;
	gridcell.dNH4dep = per_day;
	gridcell.dNO3dep = per_day;

	climate.co2 = co2;
	climate.temp = dtemp[date.day];
	climate.prec = dprec[date.day];
	climate.insol = dinsol[date.day];

	// climate.dtr does NOT come from the driver file. It means the range of the
	// near-surface AIR temperature: cfinput.cpp builds it as
	// dmax_temp - dmin_temp from air temperature extrema, and its one reader,
	// bvoc.cpp's daytime_temp, reconstructs a daytime air temperature from it.
	//
	// What the driver file carries is the range of the SURFACE temperature, and
	// the two are different variables rather than two estimates of one. On the
	// bootstrap climatology maxt and mint bracket ts in every cell and fail to
	// bracket tas in 15,561 cell-bins of 24,576, by up to 28.3 K.
	//
	// The air-temperature extrema this world's climate model computes are
	// atsama and atsami, output codes 201 and 202. They are written by the model
	// and reach no product: they are absent from pyburn's ilibrary and from
	// run_exoplasim.REGULAR_CODES. Until they arrive, substituting the surface
	// range would be delivering a different variable under the right name.
	// dtsrange is carried alongside for a consumer that wants the surface range
	// for itself.
	//
	// biosphere/notes/ecological-forcing-field-contract.md carries the argument.
	//
	// Zero rather than left alone: Climate's constructor does not initialise
	// dtr, and an indeterminate read is worse than a declared one. Zero states
	// the assumption the model then runs on, that leaf temperature equals the
	// daily mean air temperature. VesperInput::init refuses ifbvoc 1 outright,
	// so nothing reaches daytime_temp on this assumption without saying so.
	climate.dtr = 0.0;

	if (date.day == 0) {

		if (date.year == nyear_spinup + nyear) {
			// This cell is finished; advance so getgridcell moves on.
			current++;
			return false;
		}

		// The whole simulation year of air temperature is already interpolated,
		// so the seasonal landmarks summergreen phenology keys on can be read off
		// it now rather than at the end of the year. Without this the first orbit
		// of every cell would run on the landmarks initdrivers starts from, which
		// assume nothing and are therefore not this cell's; with it, the forcing
		// decides them from day 0. Over a driver file carrying several years they
		// follow the cycle, since each year is handed over as it is loaded.
		climate.set_seasonal_cycle(dtemp);

		if (tmute.getprogress() >= 1.0) {
			double progress = (double)(current * (nyear_spinup + nyear) + date.year)
			                / (double)(cells.size() * (nyear_spinup + nyear));
			tprogress.setprogress(progress);
			dprintf("%3d%% complete, %s elapsed, %s remaining\n",
			        (int)(progress * 100.0),
			        tprogress.elapsed.str, tprogress.remaining.str);
			tmute.settimer(MUTESEC);
		}
	}

	return true;
}
