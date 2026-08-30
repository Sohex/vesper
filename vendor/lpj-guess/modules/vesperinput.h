///////////////////////////////////////////////////////////////////////////////////////
/// \file vesperinput.h
/// \brief Input module for Vesper, driven by an ExoPlaSim climatology
///
/// Reads a single self-describing binary file built by
/// biosphere/scripts/build_lpj_driver.py, which carries the gridlist, the soil
/// codes derived from World Orogen lithology, and CHRONOLOGICAL climate for
/// every land cell. One file rather than the several text files demoinput uses,
/// because it is generated rather than curated and its provenance travels
/// beside it.
///
/// THE FILE CARRIES ITS OWN INTERVALS AND THIS MODULE INTEGRATES THEM ONTO ITS
/// DAY. Up to VESPDRV7 the format had a fixed twelve bins a year and this
/// module knew what a bin meant, so the forcing's cadence lived in the reader
/// and an interval's length lived nowhere at all. VESPDRV8 carries a table of
/// intervals with explicit start, end and duration in absolute seconds and the
/// local solar phase of each, any count, and every field as an intensive
/// quantity over its interval. The 24-hour hydrology and biogeochemistry
/// boundary is therefore this module's and not the format's: `integrate_year`
/// takes each absolute day's duration-weighted mean over the intervals
/// overlapping it. Twelve intervals a year work, one per day works, one per
/// climate-model timestep works, and none of them is a format change.
///
/// What that gives up is the smooth daily curve `interp_monthly_means_conserve`
/// manufactured between bin centres, and it had no source: the producer states
/// an interval mean and says nothing about the shape inside the interval, so
/// the smooth reconstruction was invented structure. The integration invents
/// nothing and is the identity when an interval is one absolute day.
/// biosphere/config/ecological_forcing_contract.yaml is the declaration.
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
		/// Climate, flattened as [year * nintervals + interval]. Sized
		/// years * nintervals. Every one of these is INTENSIVE over its
		/// interval, which is what lets one integration operator serve them
		/// all and what stops a depth being read as a rate.
		/// mean air temperature over the interval, degrees C
		std::vector<double> temp;
		/// mean precipitation RATE over the interval, mm per absolute day.
		/// Not a per-interval total: a total would need an interval length
		/// chosen to form it over, and the depth for a day is this times one
		/// absolute day.
		std::vector<double> prec;
		/// mean net downward surface shortwave over the interval, W/m2
		std::vector<double> insol;
		/// Mean range of the SURFACE temperature over the interval, degrees C.
		/// NOT the near-surface air diurnal range. It is maxt - mint, and
		/// ExoPlaSim builds those two as extrema of dt(:,NLEP), which brackets
		/// ts and not tas. climate.dtr means an air-temperature range, so this
		/// is deliberately not assigned to it; see
		/// biosphere/notes/ecological-forcing-field-contract.md.
		std::vector<double> tsrange;
	};

	/// One forcing interval, exactly as the driver file states it.
	///
	/// Bounds are absolute seconds from the start of the simulation year and
	/// are read, never inferred from a record index: a consumer that infers an
	/// interval from its position cannot tell a missing interval from a short
	/// one, and that distinction is what `read_driver` refuses on.
	struct Interval {
		/// Absolute seconds from the start of the year to this interval's start
		double start;
		/// Absolute seconds to its end
		double end;
		/// end - start, carried as well as the bounds so a gap is visible in
		/// one record rather than only in a difference across two
		double duration;
		/// Local solar phase at the start, as a fraction of a rotation. It
		/// advances every interval, because the absolute day and this world's
		/// rotation are not the same length. Carried for a consumer that
		/// integrates within a rotation; PCAR-1 owns that end.
		double phase;
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

	/// Manifest-recorded root of all stochastic ecological substreams.
	int root_seed;

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

	/// Daily climate for the current cell, integrated from the intervals
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

	/// Forcing intervals per simulation year, read from the file. Not a
	/// constant: the format states it and this module is indifferent to it.
	int nintervals;

	/// Samples inside each interval. Zero today, and `read_driver` refuses a
	/// non-zero count by name: EFOR-3 delivers the subdaily arm and PCAR-1 and
	/// FIRE-1 consume it, and reading a field with no operator behind it is
	/// worse than refusing to.
	int nsubdaily;

	/// The interval table, flattened as [year * nintervals + interval]. One
	/// table for the whole grid, because every cell shares one time axis.
	std::vector<Interval> intervals;

	/// Which year of the cycle is currently integrated, -1 for none
	int loaded_year;

	/// Integrates one year of the current cell's intervals onto absolute days
	void integrate_year(const Cell& cell, int year_index);

	/// Applies the contract's per-layer usable shares to the soil water column
	void apply_column_geometry(Gridcell& gridcell,
	                           const std::vector<double>& layer_usable);

	/// Reads the driver file into `cells`, failing loudly on any mismatch
	void read_driver();
};

#endif // LPJ_GUESS_VESPERINPUT_H
