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

| LITH-5 | Explain or correct the composition divergence from GLiM: metamorphic 2.1x Earth, evaporite ~10x, volcanic 0.5x | `notes/audits/orogen-lithology.md` | open |
| LITH-6 | Decide whether `surface_rock` should be renamed; it carries consolidated lithology, and GLiM's largest land class (unconsolidated, 24.6%) has no counterpart | `notes/audits/orogen-lithology.md` | open |
| GRAV-5 | Bracket the slope exponent | `notes/audits/orogen-gravity.md` | done -- needs NO build: gravity is a post-hoc scalar, so n=2 is the existing export rescaled by 1.1427. Land mean 0.407 -> 0.465 km, peak 4.593 -> 5.249 km, land fraction and basin identity unchanged by construction. Worth ~0.4 K by lapse rate |
| GRAV-4 | State `n = 1` explicitly at the erosion call site; the whole gravity correctness argument rests on it and nothing says so | `notes/audits/orogen-gravity.md` | done -- stated in `vendor/orogen/js/terrain-post.js`, including that n=1 is baked into the Braun-Willett closed form rather than being a parameter, and that the post-hoc 1/g scaling depends on it |
| GRAV-6 | Glacial erosion carries no gravity term, and Glen's law makes ice velocity go as `g^3` (2.23x here). Needs the glacial model rebuilt, not a factor bolted on; justify on its own merits | `notes/audits/orogen-gravity.md` | open |
| SURF-1 | Implement the derived-surface classifier against the two-axis design | `pedology/notes/derived-surface-classes.md` | blocked on 0.945 climatology |
| SURF-2 | Re-run `brine_paths.py` weighted by discharge rather than catchment area | `pedology/notes/derived-surface-classes.md` | blocked on 0.945 climatology |
| SURF-3 | Loess: needs a dust emission scheme and a transport path, neither of which exists | `pedology/notes/derived-surface-classes.md` | blocked |
| REF-2 | NOT HELD: Selby (1980), *A rock mass strength classification for geomorphic purposes*, Z. Geomorph. 24, 31-51, `10.1127/zfg/24/1984/31`. No full text reachable. Grounded through Bursztyn 2015, which applies it at 168 localities | `references/INDEX.md` | open |
| REF-4 | NOT HELD: Penndorf (1956), *Luminous and Spectral Reflectance as Well as Colors of Natural Objects*, `10.21236/AD0098766`. DTIC serves a broken TLS chain | `references/INDEX.md` | open |

| REF-6 | DOES NOT EXIST: no pyranometer broadband albedo for a gypsum crust, and none for an alkaline trona/natron pan. Both must be constructed by spectral integration | `notes/audits/orogen-lithology.md` | open |

## Done

| id | task | source | closed |
| --- | --- | --- | --- |
| LITH-4 | Fetch bare-rock and evaporite albedo literature | `notes/audits/orogen-lithology.md` | done -- 16 sources fetched and indexed, every DOI confirmed against Crossref before fetching |
| LITH-10 | Pin down what `pelagic` means | `notes/audits/orogen-lithology.md` | done -- the class is `Pelagic ooze / abyssal clay`, assigned in the ocean domain with thickness growing by seafloor age. Unconsolidated, so 3.00 is correct beside Moosdorf's 3.2 |
| LITH-8 | `carbonate` erodibility 1.30 -> 0.45 | `notes/audits/orogen-lithology.md` | done -- applied |
| LITH-9 | `schist` erodibility 1.10 -> 0.45 | `notes/audits/orogen-lithology.md` | done -- applied |
| LITH-11 | `evaporite` density 2.2 -> 2.1 | `notes/audits/orogen-lithology.md` | done -- applied |
| REF-3 | Zhuang et al. (2023) solid-rock vs powder reflectance | `references/INDEX.md` | done -- SUPERSEDED. Paragas 2025, via the POSEIDON surface albedo database, measures the same rocks as slab, crushed and powder, which is the same conversion. Median powder/slab 2.74 |
| LITH-2 | Basin-fill albedo: do NOT split evaporite by mineralogy with gypsum brighter -- band-weighted under K2.5V gyps | `notes/audits/orogen-lithology.md` | done -- decided NOT to split evaporite by mineralogy. Within-mineral spread (halite 0.311-0.877) exceeds the between-mineral difference, and surface state dominates both |
| LITH-3 | Ground the 60 uncited rock-class constants | `notes/audits/orogen-lithology.md` | done -- erodibility against Moosdorf/Stock/Zondervan, density against Daly 1966, albedo against Paragas slabs, Longhi rock-slab classes and Post's soils. All three grounded |
| LITH-7 | Narrow the erodibility spread to ~4x via `--lithology-strength 0.682`. MEASURED: combined with LITH-8/9 the ma | `notes/audits/orogen-lithology.md` | done -- applied, --lithology-strength 0.682 in build precarve-zoned-g1281-erod4x |
| LITH-12 | Rebuild the albedo table by solar-weighting against `k25v_hr.dat`. UNBLOCKED: POSEIDON's surface database has  | `notes/audits/orogen-lithology.md` | done -- table rebuilt by solar-weighting measured spectra against k25v_hr.dat. Basalt confirmed, felsic confirmed, quartzite 0.35->0.45, playa 0.30->0.19 |
| LITH-13 | Couple basin-fill albedo to wetness: measured halite crusts swing 0.64->0.24 on wetting, and `surface_water.py | `notes/audits/orogen-lithology.md` | done -- basin-fill albedo is now wetness-coupled: lakes composited from surface_water.nc, and salt crust derived as the zone whose lake depth is at or below local annual evaporation |
| LITH-15 | `playa_clastic` 0.30 sits at the bright end of Post's 52 measured soils (mean 0.189, max 0.402); damp playa is | `notes/audits/orogen-lithology.md` | done -- playa 0.30 -> 0.19, Post et al. (2000) pyranometer mean of 52 soils |
| LITH-17 | Ground or correct the evaporite share: `LITHO_SALT_CRUST_DEPTH_FRAC = 0.25` gives 10.7% of basin fill as salt  | `notes/audits/orogen-lithology.md` | done -- derived rather than tuned. Salt crust is the ephemeral zone from the water balance, 0.47% of land against the geometric rule's 2.85%. Fixed downstream because Orogen has no water balance; its geometric rule still serves erodibility, which runs before any climate exists |
| LITH-18 | Decide where between slab and powder this world's bare rock sits. Paragas slabs are a floor, particulate sets  | `notes/audits/orogen-lithology.md` | done -- resolved by two routes. The pedology model puts regolith at 0.11-3.15 m over 90% of land so the surface is optically regolith, and Longhi's rock-SLAB classes put quartz-rich gneisses at 26-35%, which is where the table already sat. Weathering has no consistent sign, so no single blend fraction exists |
| LITH-1 | Assign evaporite mineralogy downstream from the chemical divide, not from Orogen's single geometric `evaporite | `notes/audits/orogen-lithology.md` | wontfix -- superseded by LITH-2. The brine divide predicts which mineral a basin grows, but albedo should not be split on it |
| REF-1 | NOT HELD: Hardie, Eugster (1970), *The evolution of closed-basin brines*, Special Publications of the Mineralo | `references/INDEX.md` | done -- Hardie & Eugster (1970) obtained, MSA_SP3_273-290.pdf |
| REF-5 | NOT HELD: Hunt & Salisbury, *Visible and near-infrared spectra of minerals and rocks* series, Modern Geology 1 | `references/INDEX.md` | done -- superseded by Hunt (1982) in Carmichael's Handbook, which consolidates the whole Modern Geology series |
| LITH-14 | Split the basalts: a shared 0.10 is the weathered state, fresh lava and tephra are below 0.05 | `notes/audits/orogen-lithology.md` | done -- NO change. 0.10 is correct at 15 km cells: regolith covers 90% of land at 0.11-3.15 m, so the surface is optically regolith (the LITH-18 argument), and fresh lava is a transient state needing active resurfacing the export does not track. Bounded anyway: all 3.01% of land basalt at 0.04 is -0.0018 land albedo, +0.041 W/m2, about +0.02 K. Cannot justify a baseline |
| LITH-19 | Ground albedo and erodibility for the arc and forearc classes, never checked because their assignment rules could not fire | `notes/audits/orogen-lithology.md` | done -- albedo needed no change (arc_andesite 0.20/granite 0.30 = 0.67 against Logan's 0.69). Erodibility corrected: arc_andesite 0.90 -> 0.50, rift_bimodal 0.95 -> 0.50, arc_basalt 0.85 -> 0.65, and the andesite-above-basalt ordering was inverted against Moosdorf. melange kept at 4x. Terrain effect measured and negligible: land mean +0.3 m, -0.00004 land albedo |
| LITH-16 | Rename `melange`, or split blueschist out | `notes/audits/orogen-lithology.md` | done -- renamed to "Subduction melange (block-in-matrix)". The assignment rule produces sheared prism, not blueschist, and it now fires in the forearc where a prism belongs: 0.00% -> 6.17% of land |
| LITH-0 | Determine whether the missing gypsum class is an Orogen oversight | `notes/audits/orogen-lithology.md` | wontfix -- it is not. Orogen assigns evaporite geometrically and has no basis for mineralogy; the decision belongs downstream. Superseded by LITH-1 |
| GRAV-0 | Determine whether unscaled bathymetry is a bug | `notes/audits/orogen-gravity.md` | wontfix -- g cancels in the isostatic balance, so scaling ocean depth would be wrong. The rule is now stated once in `scaledHeightKm` instead of being triplicated |
| GRAV-1 | Decide whether to give Orogen's erosion physical units so gravity enters generation | `notes/audits/orogen-gravity.md` | wontfix -- the erosion law is n = 1, at which post-hoc 1/g scaling and correct-gravity erosion are the SAME operation, and the network is at grade nearly everywhere. Difference is 5.7% of local relief on single-cell headwaters and nothing elsewhere. Superseded by GRAV-4/5/6 |
| GRAV-2 | Re-examine basin statistics for the gravity texture error | `notes/audits/orogen-gravity.md` | wontfix -- rested on the superseded hillslope reasoning. The fluvial scaling is correct, so any basin-size effect must come from glacial overdeepening instead: GRAV-6 |
| GRAV-3 | Revisit the glacier rough pass for the slope error | `notes/glacier-rough-pass.md` | wontfix -- there is no slope error to inherit. The pass keys on relief amplitude, which the 1/g scaling gets right |
