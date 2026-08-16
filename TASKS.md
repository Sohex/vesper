# Tasks

Atomic, trackable work items. **Findings documents do not carry their own to-do
lists** -- an audit says what is true, this file says what to do about it, and
the two are kept apart so a finding can be read without being re-litigated and a
task can be closed without editing the argument that produced it.

Every task cites the document that justifies it. If a task has no source, it is
not yet a finding and probably needs one.

Status: `open` | `doing` | `blocked` | `done` | `wontfix` (with a reason).

## Conventions

- IDs are stable and never reused. Prefix by area: `LITH` lithology, `GRAV`
  gravity, `SURF` derived surface classes, `CLIM` climate, `HYD` hydrography,
  `BIO` biosphere, `REF` references and provenance.
- One line per task. If it needs a paragraph, it needs a findings document.
- A task that turns out to be wrong is closed `wontfix` with the reason, not
  deleted. The reasoning is the point.

---

## Open

| id | task | source | status |
| --- | --- | --- | --- |
| LITH-1 | Assign evaporite mineralogy downstream from the chemical divide, not from Orogen's single geometric `evaporite` class | `notes/audits/orogen-lithology.md` | open |
| LITH-2 | Basin-fill albedo: do NOT split evaporite by mineralogy with gypsum brighter -- band-weighted under K2.5V gypsum is 0.528 against clean halite 0.55-0.65, so that pushes the error the wrong way | `notes/audits/orogen-lithology.md` | open |
| LITH-12 | Rebuild the albedo table by solar-weighting spectra against `k25v_hr.dat`. splib07 covers the EVAPORITES (halite 11, gypsum 11, trona 6) but has zero granite/sandstone/gneiss/schist/quartzite/andesite/gabbro spectra of 3,156, so the crystalline classes still need ECOSTRESS | `notes/audits/orogen-lithology.md` | open |
| LITH-13 | Couple basin-fill albedo to wetness: measured halite crusts swing 0.64->0.24 on wetting, and `surface_water.py` already knows which basins hold water. A static 0.50 is wrong by up to 0.28 there | `notes/audits/orogen-lithology.md` | open |
| LITH-14 | Split the basalts: a shared 0.10 is the weathered state, fresh lava and tephra are below 0.05 | `notes/audits/orogen-lithology.md` | open |
| LITH-15 | `playa_clastic` 0.30 sits at the bright end of Post's 52 measured soils (mean 0.189, max 0.402); damp playa is 0.10-0.20 | `notes/audits/orogen-lithology.md` | open |
| LITH-3 | Ground the 60 uncited rock-class constants | `notes/audits/orogen-lithology.md` | doing -- erodibility and density grounded against 10 sources; albedo outstanding |
| LITH-7 | Narrow the erodibility spread to ~4x via `--lithology-strength 0.682`. MEASURED: combined with LITH-8/9 the mask is bit-identical, basins keep their ids, mean orography moves -3.0 m, but finished basin volume moves 8.1% median so hydrography and the carve verdict must be recomputed | `notes/audits/orogen-lithology.md` | open |
| LITH-16 | Rename `melange` to the subduction melange it actually models, or split blueschist out. One erodibility cannot serve a sheared block-in-matrix unit and a competent high-grade rock | `notes/audits/orogen-lithology.md` | open |
| LITH-17 | Ground or correct the evaporite share: `LITHO_SALT_CRUST_DEPTH_FRAC = 0.25` gives 10.7% of basin fill as salt against Earth's ~1.5%, about 7x. One-parameter fix, but GLiM's `ev` is a bedrock class while Orogen's is a surface crust | `notes/audits/orogen-lithology.md` | open |
| LITH-5 | Explain or correct the composition divergence from GLiM: metamorphic 2.1x Earth, evaporite ~10x, volcanic 0.5x | `notes/audits/orogen-lithology.md` | open |
| LITH-6 | Decide whether `surface_rock` should be renamed; it carries consolidated lithology, and GLiM's largest land class (unconsolidated, 24.6%) has no counterpart | `notes/audits/orogen-lithology.md` | open |
| GRAV-5 | Bracket the slope exponent | `notes/audits/orogen-gravity.md` | done -- needs NO build: gravity is a post-hoc scalar, so n=2 is the existing export rescaled by 1.1427. Land mean 0.407 -> 0.465 km, peak 4.593 -> 5.249 km, land fraction and basin identity unchanged by construction. Worth ~0.4 K by lapse rate |
| GRAV-4 | State `n = 1` explicitly at the erosion call site; the whole gravity correctness argument rests on it and nothing says so | `notes/audits/orogen-gravity.md` | open |
| GRAV-6 | Glacial erosion carries no gravity term, and Glen's law makes ice velocity go as `g^3` (2.23x here). Needs the glacial model rebuilt, not a factor bolted on; justify on its own merits | `notes/audits/orogen-gravity.md` | open |
| SURF-1 | Implement the derived-surface classifier against the two-axis design | `pedology/notes/derived-surface-classes.md` | blocked on 0.945 climatology |
| SURF-2 | Re-run `brine_paths.py` weighted by discharge rather than catchment area | `pedology/notes/derived-surface-classes.md` | blocked on 0.945 climatology |
| SURF-3 | Loess: needs a dust emission scheme and a transport path, neither of which exists | `pedology/notes/derived-surface-classes.md` | blocked |
| REF-2 | NOT HELD: Selby (1980), *A rock mass strength classification for geomorphic purposes*, Z. Geomorph. 24, 31-51, `10.1127/zfg/24/1984/31`. No full text reachable. Grounded through Bursztyn 2015, which applies it at 168 localities | `references/INDEX.md` | open |
| REF-7 | NOT HELD: Logan, Hunt, Salisbury, Balsamo (1973), *Compositional implication of Christensen frequency maximums for infrared remote sensing applications*, J. Geophys. Res. 78(23), 4983-5003. Defines the normalisation of Hunt 1982 Table 9's albedo column, which lists a value of 1.13 and so is not a plain reflectance | `references/INDEX.md` | open |
| REF-3 | NOT HELD: Zhuang et al. (2023), *Visible and near-infrared reflectance spectra of igneous rocks and their powders*, Icarus 391, 115346, `10.1016/j.icarus.2022.115346`. The best modern source for solid-rock vs powder reflectance of basalt/andesite/granite. Worth an interlibrary request | `references/INDEX.md` | open |
| REF-4 | NOT HELD: Penndorf (1956), *Luminous and Spectral Reflectance as Well as Colors of Natural Objects*, `10.21236/AD0098766`. DTIC serves a broken TLS chain | `references/INDEX.md` | open |
| REF-5 | NOT HELD: Hunt & Salisbury, *Visible and near-infrared spectra of minerals and rocks* series, Modern Geology 1970-1976. No DOIs ever deposited. Covered indirectly by ASTER and USGS libraries | `references/INDEX.md` | open |
| REF-6 | DOES NOT EXIST: no pyranometer broadband albedo for a gypsum crust, and none for an alkaline trona/natron pan. Both must be constructed by spectral integration | `notes/audits/orogen-lithology.md` | open |
| REF-1 | NOT HELD: Hardie, Eugster (1970), *The evolution of closed-basin brines*, Special Publications of the Mineralogical Society of America 3, 273-290. No DOI; society special publication. Cited through Eugster & Jones 1979 and Deocampo & Jones 2014, both of which restate it in full | `references/INDEX.md` | open |

## Done

| id | task | source | closed |
| --- | --- | --- | --- |
| LITH-4 | Fetch bare-rock and evaporite albedo literature | `notes/audits/orogen-lithology.md` | done -- 16 sources fetched and indexed, every DOI confirmed against Crossref before fetching |
| LITH-10 | Pin down what `pelagic` means | `notes/audits/orogen-lithology.md` | done -- the class is `Pelagic ooze / abyssal clay`, assigned in the ocean domain with thickness growing by seafloor age. Unconsolidated, so 3.00 is correct beside Moosdorf's 3.2 |
| LITH-8 | `carbonate` erodibility 1.30 -> 0.45 | `notes/audits/orogen-lithology.md` | done -- applied |
| LITH-9 | `schist` erodibility 1.10 -> 0.45 | `notes/audits/orogen-lithology.md` | done -- applied |
| LITH-11 | `evaporite` density 2.2 -> 2.1 | `notes/audits/orogen-lithology.md` | done -- applied |
| LITH-0 | Determine whether the missing gypsum class is an Orogen oversight | `notes/audits/orogen-lithology.md` | wontfix -- it is not. Orogen assigns evaporite geometrically and has no basis for mineralogy; the decision belongs downstream. Superseded by LITH-1 |
| GRAV-0 | Determine whether unscaled bathymetry is a bug | `notes/audits/orogen-gravity.md` | wontfix -- g cancels in the isostatic balance, so scaling ocean depth would be wrong. The rule is now stated once in `scaledHeightKm` instead of being triplicated |
| GRAV-1 | Decide whether to give Orogen's erosion physical units so gravity enters generation | `notes/audits/orogen-gravity.md` | wontfix -- the erosion law is n = 1, at which post-hoc 1/g scaling and correct-gravity erosion are the SAME operation, and the network is at grade nearly everywhere. Difference is 5.7% of local relief on single-cell headwaters and nothing elsewhere. Superseded by GRAV-4/5/6 |
| GRAV-2 | Re-examine basin statistics for the gravity texture error | `notes/audits/orogen-gravity.md` | wontfix -- rested on the superseded hillslope reasoning. The fluvial scaling is correct, so any basin-size effect must come from glacial overdeepening instead: GRAV-6 |
| GRAV-3 | Revisit the glacier rough pass for the slope error | `notes/glacier-rough-pass.md` | wontfix -- there is no slope error to inherit. The pass keys on relief amplitude, which the 1/g scaling gets right |
