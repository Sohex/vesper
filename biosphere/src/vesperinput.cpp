///////////////////////////////////////////////////////////////////////////////////////
/// \file vesperinput.cpp
/// \brief Input module for Vesper, driven by an ExoPlaSim climatology
///
/// See vesperinput.h. Modelled on demoinput.cpp, which is the module the
/// LPJ-GUESS documentation points new users at, and departs from it in four
/// ways, each forced by this being another planet:
///
///  1. Insolation is supplied as NETSWRAD_TS, the net downward surface shortwave
///     ExoPlaSim already computes, rather than as percentage sunshine. The
///     climate model knows the surface albedo far better than driver.cpp's
///     global BETA constant of 0.17 does, and this project computes that albedo
///     from lithology.
///  2. Nitrogen deposition is divided by the year length rather than by 365.
///  3. The gridlist, soil codes and climate arrive in one generated binary file
///     rather than several curated text files.
///  4. The file's year length is checked against the compiled-in one.
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
const char DRIVER_MAGIC[8] = {'V','E','S','P','D','R','V','4'};

/// Bins per year in the driver file. ExoPlaSim's regular_output_bins_per_orbit.
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
		read_or_fail(in, &cell.regolith_depth_m, 1, "regolith depth");
		read_or_fail(in, &cell.bedrock_water_fraction, 1, "bedrock water fraction");
		if (cell.soilcode < 0 || cell.soilcode > 9) {
			fclose(in);
			fail("vesperinput: cell %d has invalid LPJ soil code %d",
			     i, cell.soilcode);
		}
		const size_t span = (size_t)nbins * (size_t)nyears;
		cell.temp.resize(span);
		cell.prec.resize(span);
		cell.insol.resize(span);
		cell.dtr.resize(span);
		read_or_fail(in, &cell.temp[0], span, "temperature");
		read_or_fail(in, &cell.prec[0], span, "precipitation");
		read_or_fail(in, &cell.insol[0], span, "insolation");
		read_or_fail(in, &cell.dtr[0], span, "diurnal range");
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
	interp_monthly_means_conserve(&cell.dtr[offset], ddtr, 0.0);
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
		soilinput.get_soil(cell.lon, cell.lat, gridcell);
	}
	else {
		soil_parameters(gridcell.soiltype, cell.soilcode);
	}

	apply_regolith_depth(gridcell, cell.regolith_depth_m,
	                     cell.bedrock_water_fraction);

	// Year 0 of this cell, so day 0 has data before getclimate first runs.
	// getclimate re-interpolates at the start of every year after that.
	loaded_year = -1;
	interpolate(cell, 0);

	clear_all_graphs();

	return true;
}

void VesperInput::apply_regolith_depth(Gridcell& gridcell, double depth_m,
                                       double bedrock_fraction) {

	// LPJ-GUESS gives every gridcell the same 1.5 m profile and derives its
	// water capacity from texture alone. On a world where relief and erodibility
	// are known, that is a real omission: a plant on 0.2 m of regolith over
	// bedrock has a seventh of the water a plant on a deep profile has, and no
	// combination of sand and clay fractions can express it.
	//
	// So each layer's capacity is scaled by the fraction of that layer which is
	// actually regolith rather than rock. A layer entirely above the bedrock
	// contact is untouched; one entirely below it holds nothing.
	//
	// This scales capacity, not the layer geometry: the layers still exist and
	// still conduct heat, they simply hold less water. Rooting depth is
	// unchanged, which slightly overstates how much of the profile deep-rooted
	// PFTs can reach on thin soil. Recorded rather than corrected, because
	// changing rootdist would reach into PFT parameters that are Earth
	// calibrations.

	// Material below the bedrock contact still holds plant-available water, and
	// deep-rooted woody plants demonstrably use it. `bedrock_fraction` is how
	// much, per unit volume, relative to the soil above, and it arrives per cell
	// from the pedology component as a function of weathering intensity. It can
	// legitimately exceed 1: deeply weathered saprolite holds more water than the
	// soil sitting on top of it.
	//
	// A hard floor is still needed, but only as a numerical guard and not as
	// physics. LPJ-GUESS carries soil water as a fraction of each layer's
	// capacity, computing `wcont = Faw_layer / soiltype.awc[layer]` in
	// soilwater.cpp and canexch.cpp, so a layer at exactly zero divides by zero
	// and the NaN propagates through the nitrogen substrate and kills the
	// gridcell silently, reporting zero LAI rather than failing.
	const double NUMERICAL_FLOOR = 1.0e-4;

	if (depth_m <= 0.0) {
		return;
	}

	double bedrock = bedrock_fraction;
	if (bedrock < NUMERICAL_FLOOR) {
		bedrock = NUMERICAL_FLOOR;
	}

	// Capped at 1: sub-bedrock material can at most match the soil above, never
	// exceed it. That is a limitation of this model, not of the world. Deeply
	// weathered saprolite really does hold two to four times what its soil does,
	// but LPJ-GUESS ties a layer's saturation capacity to its texture-derived
	// porosity, and scaling `wsats` past that drives `Frac_air` negative in
	// Soil::update_soil_diffusivities, which is a hard failure.
	//
	// Representing the upper half of the observed range needs a genuinely
	// separate bedrock layer with its own porosity, which is what Lapides et al.
	// (Biogeosciences 2024) added to this same model rather than rescaling the
	// existing profile. Recorded here so the ceiling is visible in results.
	if (bedrock > 1.0) {
		bedrock = 1.0;
	}

	const double profile_mm = SOILDEPTH_UPPER + SOILDEPTH_LOWER;
	double depth_mm = depth_m * 1000.0;
	if (depth_mm >= profile_mm && bedrock >= 1.0) {
		return;   // wholly regolith, or rock that holds as much; nothing to do
	}

	Soiltype& soil = gridcell.soiltype;

	const double upper_layer_mm = SOILDEPTH_UPPER / (double)NSOILLAYER_UPPER;
	const double lower_layer_mm = SOILDEPTH_LOWER
		/ (double)(NSOILLAYER - NSOILLAYER_UPPER);

	double top_mm = 0.0;
	double kept_upper = 0.0;
	double kept_lower = 0.0;
	for (int layer = 0; layer < NSOILLAYER; layer++) {
		const double thickness = (layer < NSOILLAYER_UPPER)
			? upper_layer_mm : lower_layer_mm;
		// Share of this layer above the bedrock contact. The rest is rock, and
		// holds `bedrock` times what the same volume of soil would.
		double regolith_share = (depth_mm - top_mm) / thickness;
		regolith_share = regolith_share < 0.0 ? 0.0
			: (regolith_share > 1.0 ? 1.0 : regolith_share);
		const double usable = regolith_share + (1.0 - regolith_share) * bedrock;

		soil.awc[layer] *= usable;
		soil.wp[layer] *= usable;
		soil.wsats[layer] *= usable;

		if (layer < NSOILLAYER_UPPER) {
			kept_upper += usable * thickness;
		}
		else {
			kept_lower += usable * thickness;
		}
		top_mm += thickness;
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
	if (date.day == 0 && wanted != loaded_year) {
		interpolate(cell, wanted);
	}

	// Per day of absolute time, so the year length rather than 365. Split evenly
	// between reduced and oxidised, as demoinput does; neither the split nor the
	// total is measured on this world, and both are declared in the driver file
	// rather than assumed here.
	const double per_day = ndep / 2.0 / (double)Date::MAX_YEAR_LENGTH * HA_PER_M2;
	gridcell.dNH4dep = per_day;
	gridcell.dNO3dep = per_day;

	climate.co2 = co2;
	climate.temp = dtemp[date.day];
	climate.prec = dprec[date.day];
	climate.insol = dinsol[date.day];
	climate.dtr = ddtr[date.day];

	if (date.day == 0) {

		if (date.year == nyear_spinup + nyear) {
			// This cell is finished; advance so getgridcell moves on.
			current++;
			return false;
		}

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
