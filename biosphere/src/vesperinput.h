///////////////////////////////////////////////////////////////////////////////////////
/// \file vesperinput.h
/// \brief Input module for Vesper, driven by an ExoPlaSim climatology
///
/// Reads a single self-describing binary file built by
/// biosphere/scripts/build_lpj_driver.py, which carries the gridlist, the soil
/// codes derived from World Orogen lithology, and binned climate for every land
/// cell. One file rather than the several text files demoinput uses, because it
/// is generated rather than curated and its provenance travels beside it.
///
/// The file records the year length it was built for. That is checked against
/// VESPER_YEAR_LENGTH_DAYS at load and is a hard failure on mismatch: this
/// world's year is a function of its stellar flux, so a driver file and a binary
/// built at different fluxes describe different planets.
///
///////////////////////////////////////////////////////////////////////////////////////

#ifndef LPJ_GUESS_VESPERINPUT_H
#define LPJ_GUESS_VESPERINPUT_H

#include "guess.h"
#include "inputmodule.h"
#include "gutil.h"
#include "externalinput.h"
#include "soilinput.h"
#include <vector>

/// Input module driving LPJ-GUESS from an ExoPlaSim climatology of Vesper
class VesperInput : public InputModule {

public:

	VesperInput();
	~VesperInput();

	/// Reads the driver file. Called after the instruction file has been read.
	void init();

	/// See base class for documentation about this function's responsibilities
	bool getgridcell(Gridcell& gridcell);

	/// See base class for documentation about this function's responsibilities
	bool getclimate(Gridcell& gridcell);

	/// See base class for documentation about this function's responsibilities
	void getlandcover(Gridcell& gridcell);

	/// Obtains land management data for one day
	void getmanagement(Gridcell& gridcell) {management_input.getmanagement(gridcell);}

private:

	/// One land cell: where it is, what it sits on, and its binned climate
	struct Cell {
		double lon;
		double lat;
		int soilcode;
		/// regolith thickness, m, from the pedology component
		double regolith_depth_m;
		/// plant-available water below the bedrock contact, as a fraction of
		/// what the soil above holds per unit volume. From pedology, a function
		/// of weathering intensity; see pedology/config/pedogenesis.yaml.
		double bedrock_water_fraction;
		/// mean air temperature per bin, degrees C
		std::vector<double> temp;
		/// precipitation total per bin, mm
		std::vector<double> prec;
		/// net downward surface shortwave per bin, W/m2
		std::vector<double> insol;
		/// diurnal temperature range per bin, degrees C (BVOC only)
		std::vector<double> dtr;
	};

	/// Land cover input module
	LandcoverInput landcover_input;
	/// Management input module
	ManagementInput management_input;

	/// Every land cell, in the order the driver file lists them
	std::vector<Cell> cells;

	/// Index of the cell currently being simulated
	size_t current;

	/// Number of simulation years to run after spin-up
	int nyear;

	/// Path to the driver file, from the instruction file
	xtring file_driver;

	/// Optional soil map from the pedology component. When set, soil physical
	/// properties come from measured texture, organic content, pH and bulk
	/// density rather than from the driver file's single LPJ soil code.
	xtring file_soilmap;

	/// Reader for that map. LPJ-GUESS's own, so the format is theirs not ours.
	SoilInput soilinput;

	/// Whether a soil map was supplied
	bool have_soilmap;

	/// Atmospheric CO2, ppm, from the driver file
	double co2;

	/// Nitrogen deposition, kgN/ha/year, from the driver file
	double ndep;

	/// Provenance string the driver file carries, echoed to the log
	xtring provenance;

	/// Daily climate for the current cell, interpolated from the bins
	double dtemp[Date::MAX_YEAR_LENGTH];
	double dprec[Date::MAX_YEAR_LENGTH];
	double dinsol[Date::MAX_YEAR_LENGTH];
	double ddtr[Date::MAX_YEAR_LENGTH];

	/// Progress reporting, as in demoinput
	Timer tprogress, tmute;
	static const int MUTESEC = 20;

	/// Interpolates the current cell's bins onto days
	void interpolate(const Cell& cell);

	/// Scales soil water capacity by how much of each layer is really regolith
	void apply_regolith_depth(Gridcell& gridcell, double depth_m,
	                          double bedrock_fraction);

	/// Reads the driver file into `cells`, failing loudly on any mismatch
	void read_driver();
};

#endif // LPJ_GUESS_VESPERINPUT_H
