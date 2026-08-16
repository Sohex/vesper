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
#include "driver.h"
#include "parameters.h"
#include "guess.h"
#include <stdio.h>
#include <string.h>

REGISTER_INPUT_MODULE("vesper", VesperInput)

namespace {

/// Little-endian, and both writer and reader are x86-64. Checked via the magic.
const char DRIVER_MAGIC[8] = {'V','E','S','P','D','R','V','1'};

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
	: current(0), nyear(1), co2(0.0), ndep(0.0) {

	declare_parameter("nyear", &nyear, 1, 10000,
		"Number of simulation years to run after spinup");
	declare_parameter("file_driver", &file_driver, 300,
		"Path to the binary driver file built by build_lpj_driver.py");
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
		fail("vesperinput: %s is not a VESPDRV1 driver file", (char*)file_driver);
	}

	int ncells = 0, nbins = 0, year_length = 0, reserved = 0;
	read_or_fail(in, &ncells, 1, "cell count");
	read_or_fail(in, &nbins, 1, "bin count");
	read_or_fail(in, &year_length, 1, "year length");
	read_or_fail(in, &reserved, 1, "reserved field");

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
		if (cell.soilcode < 0 || cell.soilcode > 9) {
			fclose(in);
			fail("vesperinput: cell %d has invalid LPJ soil code %d",
			     i, cell.soilcode);
		}
		cell.temp.resize(nbins);
		cell.prec.resize(nbins);
		cell.insol.resize(nbins);
		cell.dtr.resize(nbins);
		read_or_fail(in, &cell.temp[0], nbins, "temperature");
		read_or_fail(in, &cell.prec[0], nbins, "precipitation");
		read_or_fail(in, &cell.insol[0], nbins, "insolation");
		read_or_fail(in, &cell.dtr[0], nbins, "diurnal range");
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
	dprintf("  %d land cells, %d bins per year, %d-day year\n",
	        (int)cells.size(), DRIVER_BINS, (int)Date::MAX_YEAR_LENGTH);
	dprintf("  CO2 %g ppm, N deposition %g kgN/ha/yr\n", co2, ndep);
	dprintf("  built from %s\n\n", (char*)provenance);

	landcover_input.init();
	management_input.init();

	tprogress.init();
	tmute.init();
	tprogress.settimer();
	tmute.settimer(MUTESEC);

	current = 0;
}

void VesperInput::interpolate(const Cell& cell) {

	// The framework's own conserving interpolators, so bin means stay means and
	// bin totals stay totals. They read date.ndaymonth[], which the patched Date
	// fills with Vesper's months, so no bin length is assumed here.
	interp_monthly_means_conserve(&cell.temp[0], dtemp);
	interp_monthly_totals_conserve(&cell.prec[0], dprec, 0.0);
	interp_monthly_means_conserve(&cell.insol[0], dinsol, 0.0);
	interp_monthly_means_conserve(&cell.dtr[0], ddtr, 0.0);
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

	soil_parameters(gridcell.soiltype, cell.soilcode);

	interpolate(cell);

	clear_all_graphs();

	return true;
}

void VesperInput::getlandcover(Gridcell& gridcell) {

	landcover_input.getlandcover(gridcell);
	landcover_input.get_land_transitions(gridcell);
}

bool VesperInput::getclimate(Gridcell& gridcell) {

	Climate& climate = gridcell.climate;

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
