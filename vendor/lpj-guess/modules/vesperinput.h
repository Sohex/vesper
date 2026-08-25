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
/// It also carries however many years of climate it was built with, and this
/// module cycles through them. One year repeats forever, which is the right
/// forcing for a fixed climate; sixteen give a stellar cycle. Nothing here
/// assumes a number, because the earlier version did assume one, and a single
/// repeating year cannot represent a variable star at all.
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
		/// The land column property contract's retention states for this cell.
		/// Read, never derived: see SoilInput::ContractStates.
		SoilInput::ContractStates soil_states;
		/// Per physical layer, the share of that layer's capacity the column
		/// carries. The contract's weathered-bedrock vertical rule, evaluated
		/// by pedology/scripts/land_column_properties.py from regolith depth
		/// and the weathered-bedrock water fraction, so the rule lives in one
		/// place rather than here and there. One number per layer and not
		/// three, because the same share scales that layer's available water,
		/// its wilting point and its saturation: the material below the
		/// regolith contact holds less of everything.
		std::vector<double> layer_usable;
		/// Climate, flattened as [year * bins + bin]. Sized years * bins.
		/// mean air temperature, degrees C
		std::vector<double> temp;
		/// precipitation total per bin, mm
		std::vector<double> prec;
		/// net downward surface shortwave per bin, W/m2
		std::vector<double> insol;
		/// Range of the SURFACE temperature per bin, degrees C.
		/// NOT the near-surface air diurnal range. It is maxt - mint, and
		/// ExoPlaSim builds those two as extrema of dt(:,NLEP), which brackets
		/// ts and not tas. climate.dtr means an air-temperature range, so this
		/// is deliberately not assigned to it; see
		/// biosphere/notes/ecological-forcing-field-contract.md.
		std::vector<double> tsrange;
	};

	/// Land cover input module
	LandcoverInput landcover_input;
	/// Management input module
	ManagementInput management_input;

	/// Every land cell, in the order the driver file lists them
	std::vector<Cell> cells;

	/// Index of the cell currently being simulated, within this rank's share
	size_t current;

	/// Cells in the driver file before this rank took its share, for reporting
	int total_cells;

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
	double dtsrange[Date::MAX_YEAR_LENGTH];

	/// Progress reporting, as in demoinput
	Timer tprogress, tmute;
	static const int MUTESEC = 20;

	/// Years of climate the driver file carries. The forcing cycles through
	/// them, so 1 is a fixed climate and 16 is roughly one stellar cycle.
	int years;

	/// Which year of the cycle is currently interpolated, -1 for none
	int loaded_year;

	/// Interpolates one year of the current cell's bins onto days
	void interpolate(const Cell& cell, int year_index);

	/// Applies the contract's per-layer usable shares to the soil water column
	void apply_column_geometry(Gridcell& gridcell,
	                           const std::vector<double>& layer_usable);

	/// Reads the driver file into `cells`, failing loudly on any mismatch
	void read_driver();
};

#endif // LPJ_GUESS_VESPERINPUT_H
