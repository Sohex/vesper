# References

Primary sources this project leans on. Filenames are
`<firstauthor><year><suffix>-<slug>.pdf`.

**Read** means someone here has opened the paper and taken the number from it.
**Held** means the project already cites it but has only ever used it secondhand.
That distinction is the point of this file: `notes/failure-modes.md` class 9 is
"a number taken from a citation rather than from the paper", and it has cost us
twice -- once on Hartmann's phosphorus rows, once on the ozone UV weight, where
the correction was itself half wrong. Everything marked *held* is an open
exposure, not a resource.

Titles are verbatim as the publisher deposited them, verified against Crossref or
the publisher landing page. Where the published title differs from a preprint's,
the published form is given.

---

## Biosphere: photosynthesis under a non-solar spectrum

The largest single item in `analysis/error_budget.json`: bare rock versus
vegetated is worth about +4 K, so what the biosphere does under a K dwarf sets
where the planet is placed. `biosphere/notes/productivity-prediction.md`
concludes that widening the photosynthetic window from 400-700 to 400-750 nm
takes photon flux from 0.81x Earth to 0.99x -- near parity -- and that conclusion
is built entirely on the first two of these.

| file | citation | status |
| --- | --- | --- |
| `kiang2007a-photosynthesis-signatures-i.pdf` | Kiang, Siefert, Govindjee, Blankenship (2007). *Spectral Signatures of Photosynthesis. I. Review of Earth Organisms.* Astrobiology 7(1), 222-251. `10.1089/ast.2006.0105` | held |
| `kiang2007b-photosynthesis-signatures-ii.pdf` | Kiang, Segura, Tinetti, Govindjee, Blankenship, Cohen, Siefert, Crisp, Meadows (2007). *Spectral Signatures of Photosynthesis. II. Coevolution with Other Stars And The Atmosphere on Extrasolar Worlds.* Astrobiology 7(1), 252-274. `10.1089/ast.2006.0108` | held -- source of "K2V pigments peak in the red-orange" and of the anoxygenic caveat that bounds the range |
| `lehmer2021-peak-absorbance-wavelength.pdf` | Lehmer, Catling, Parenteau, Kiang, Hoehler (2021). *The Peak Absorbance Wavelength of Photosynthetic Pigments Around Other Stars From Spectral Optimization.* Front. Astron. Space Sci. 8, 689441. `10.3389/fspas.2021.689441` | held -- **the near-parity argument turns on its K2V peaks of 675/711/746 nm** |
| `marosvolgyi2010-cost-and-color-of-photosynthesis.pdf` | Marosvölgyi, van Gorkom (2010). *Cost and color of photosynthesis.* Photosynth. Res. 103(2), 105-109. `10.1007/s11120-009-9522-3` | held -- the optimisation model Lehmer applies; Lehmer is itself a secondhand use of it |

Note the published Kiang titles use periods (`. II.`); the colon form is the
arXiv preprint, and these files are the preprints.

## Biosphere: the carbon-nitrogen-phosphorus model

| file | citation | status |
| --- | --- | --- |
| `dantasdepaula2025-lpjguess-phosphorus-cycle.pdf` | Dantas de Paula, Forrest, Wårlind, Darela Filho, Fleischer, Rammig, Hickler (2025). *Including the phosphorus cycle into the LPJ-GUESS dynamic global vegetation model (v4.1, r10994) - global patterns and temporal trends of N and P primary production limitation.* Geosci. Model Dev. 18, 2249-2274. `10.5194/gmd-18-2249-2025` | read -- the fork being ported, see `biosphere/notes/cnp-fork-scoping.md` |
| `wang2010-cnp-terrestrial-biosphere.pdf` | Wang, Law, Pak (2010). *A global model of carbon, nitrogen and phosphorus cycles for the terrestrial biosphere.* Biogeosciences 7, 2261-2282. `10.5194/bg-7-2261-2010` | read |
| `prentice1993-forest-landscape-transient-model.pdf` | Prentice, Sykes, Cramer (1993). *A simulation model for the transient effects of climate change on forest landscapes.* Ecol. Model. 65(1-2), 51-70. `10.1016/0304-3800(93)90126-D` | read -- solar geometry is **equations 7-13**, not 7-14; eq 8 hardcodes e = 0.01675 and eq 10 hardcodes 23.4 degrees obliquity, both of which need replacing for this world. Not the 1992 BIOME paper. |
| `lieth1975-miami-model.pdf` | Lieth (1975). *Modeling the Primary Productivity of the World.* In Lieth & Whittaker (eds), Primary Productivity of the Biosphere, Ecological Studies 14, 237-263. `10.1007/978-3-642-80913-2_12` | held -- the Miami model, used as an empirical cross-check |

## Hydrography: evaporation over dry ground

`analysis/error_budget.json` lists "Penman over a dry land column" as the only
item whose magnitude is unknown, and its sign is one-way: over a subgrid lake in
a dry column, VPD is too high, E is overstated, and the carve verdict
under-carves. These quantify that.

| file | citation | status |
| --- | --- | --- |
| `brutsaert1979-advection-aridity.pdf` | Brutsaert, Stricker (1979). *An advection-aridity approach to estimate actual regional evapotranspiration.* Water Resour. Res. 15(2), 443-450. `10.1029/WR015i002p00443` | held |
| `morton1983a-areal-evapotranspiration.pdf` | Morton (1983). *Operational estimates of areal evapotranspiration and their significance to the science and practice of hydrology.* J. Hydrol. 66(1-4), 1-76. `10.1016/0022-1694(83)90177-4` | held |
| `morton1983b-lake-evaporation.pdf` | Morton (1983). *Operational estimates of lake evaporation.* J. Hydrol. 66(1-4), 77-100. `10.1016/0022-1694(83)90178-6` | held -- **the one that addresses open water rather than areal ET**, so the directly relevant one for playa lakes |
| `han2020-complementary-principle-review.pdf` | Han, Tian (2020). *A review of the complementary principle of evaporation: from the original linear relationship to generalized nonlinear functions.* Hydrol. Earth Syst. Sci. 24(5), 2269-2285. `10.5194/hess-24-2269-2020` | held -- modern review spanning Bouchet, CRAE and advection-aridity, from arid-basin data |
| `crago2021-comment-on-han-tian.pdf` | Crago, Szilagyi, Qualls (2021). *Comment on: "A review of the complementary principle of evaporation..." by Han and Tian (2020).* Hydrol. Earth Syst. Sci. 25(1), 63-68. `10.5194/hess-25-63-2021` | held -- read with the above, not instead of it |
| `deandrade2025-complementary-relationship-intercomparison.pdf` | Comini de Andrade, Huntington, Volk, Morton, Pearson, Albano (2025). *Multi-Model Intercomparison of the Complementary Relationship of Evaporation Across Global Environmental Settings.* Water Resour. Res. 61(9), e2024WR039740. `10.1029/2024WR039740` | held -- this copy is the ESSOAr preprint |

There is **no** dedicated review of the complementary relationship for small water
bodies in arid settings. It exists only as sections inside the above.

## Star: activity, cycles and spectrum

`config/planet.yaml` sets a bolometric cycle of 0.91 to 1.01, i.e. 10.4%
peak-to-peak, and labels it undetermined. These bound it. **The quantities are
not interchangeable** and conflating them is the trap: rotational modulation,
activity-cycle amplitude, spot filling factor and S-index are four different
things, and only the second is a brightness cycle.

| file | citation | constrains | status |
| --- | --- | --- | --- |
| `suarezmascareno2016-magnetic-cycles-late-type.pdf` | Suárez Mascareño, Rebolo, González Hernández (2016). *Magnetic cycles and rotation periods of late-type stars from photometric time series.* A&A 595, A12. `10.1051/0004-6361/201628586` | **cycle period and photometric cycle amplitude together** -- the closest to what we need | held |
| `mcquillan2014-kepler-rotation-periods.pdf` | McQuillan, Mazeh, Aigrain (2014). *Rotation Periods of 34,030 Kepler Main-Sequence Stars: The Full Autocorrelation Sample.* ApJS 211(2), 24. `10.1088/0067-0049/211/2/24` | rotational modulation amplitude, **not** a cycle | held |
| `cao2022-starspots-pleiades-m67.pdf` | Cao, Pinsonneault (2022). *Star-spots and magnetism: testing the activity paradigm in the Pleiades and M67.* MNRAS 517(2), 2165-2189. `10.1093/mnras/stac2706` | spot filling factor measured directly from spectra; a **static** ceiling, not a swing | held |
| `morris2017-hatp11-starspots.pdf` | Morris, Hebb, Davenport, Rohn, Hawley (2017). *The Starspots of HAT-P-11: Evidence for a Solar-like Dynamo.* ApJ 846(2), 99. `10.3847/1538-4357/aa8555` | filling factor on one K4 dwarf, via transit occultations | held |
| `baliunas1995-chromospheric-variations.pdf` | Baliunas et al., 27 authors (1995). *Chromospheric variations in main-sequence stars.* ApJ 438, 269-287. `10.1086/175072` | S-index cycle periods. **Scanned facsimile, no text layer.** | held |
| `borosaikia2018-chromospheric-activity-catalogue.pdf` | Boro Saikia et al., 9 authors (2018). *Chromospheric activity catalogue of 4454 cool stars.* A&A 616, A108. `10.1051/0004-6361/201629518` | S-index and log R'HK only; does **not** convert to bolometric flux without an assumed relation | held |
| `segura2003-ozone-uv-other-stars.pdf` | Segura, Krelove, Kasting, Sommerlatt, Meadows, Crisp, Cohen, Mlawer (2003). *Ozone Concentrations and Ultraviolet Fluxes on Earth-Like Planets Around Other Stars.* Astrobiology 3(4), 689-708. | read -- source of `model.ozone_scale` | read |

## Geochemistry: weathering, phosphorus and the thermostat

| file | citation | status |
| --- | --- | --- |
| `hartmann2014-weathering-phosphorus-release.pdf` | Hartmann, Moosdorf, Lauerwald, Hinderer, West (2014). *Global chemical weathering and associated P-release - The role of lithology, temperature and soil properties.* Chemical Geology 363, 145-163. | **read** -- the source of every per-lithology phosphorus number in `pedology/config/pedogenesis.yaml`, both the content row and the release row. It was read from a secondhand citation first and the wrong row was taken; the correction was then also wrong for the intended use. Both rows are now carried and this file is why. It had been sitting outside the repository entirely until 2026-08-16. |
| `walker1981-whak-thermostat.pdf` | Walker, Hays, Kasting (1981). *A negative feedback mechanism for the long-term stabilization of Earth's surface temperature.* J. Geophys. Res. 86(C10), 9776-9782. `10.1029/JC086iC10p09776` | held -- the thermostat efficiency figure came from another session and the primary was never held. **Scanned, thin text layer.** |
| `dunne1978-chemical-denudation-silicate.pdf` | Dunne, T. (1978). *Rates of chemical denudation of silicate rocks in tropical catchments.* Nature 274, 244-246. | **read** -- source of the runoff exponent in the weathering law. **Scanned, no text layer; read as an image.** |

Dunne's Fig. 1 caption:

    chemical denudation rate = 0.28 R^0.66;  n = 30(43);  S_y.x = 0.13 log units

**Do not read that 0.66 as the weathering law's exponent.** Reading this paper
also traced the chain properly, and it is not the one this project recorded.
WHAK does not cite Dunne and its runoff exponent is 1. The 0.65 used everywhere
downstream is Berner (1994) GEOCARB II, combining Dunne's Kenyan `C ~ R^-0.4`
with Peters (1984)'s `C ~ R^-0.3` into a global `C ~ R^-0.35`, hence flux
`~ R^0.65`. Dunne alone gives 0.6 on Berner's reading, because Berner fits
concentration where Dunne's caption reports flux.

The genuinely new number is **S_y.x = 0.13 log units, a factor of 1.35 scatter
about the fit** -- the real uncertainty on the runoff term, never previously
recorded. Applicability: 30 independent Kenyan catchments of 43 sampled, silicate
lithologies only, rainfall 250-3000 mm, mean annual air temperature 10-30 C.

**Still to fetch:** Berner, R. A. (1994), *GEOCARB II: A revised model of
atmospheric CO2 over Phanerozoic time*, Am. J. Sci. 294(1), 56-91,
`10.2475/ajs.294.1.56` -- the actual source of the 0.65, and currently held
secondhand. Also Peters, N. E. (1984), the US-rivers half of the synthesis.

## Aerosol optics

Assembled to settle the sign of dust radiative forcing over this world's bright
closed-basin surfaces. Resolved: cooling everywhere. See `notes/dust.md`.

| file | citation | status |
| --- | --- | --- |
| `dibiagio2019-dust-shortwave-refractive-index.pdf` | Di Biagio et al., 15 authors (2019). *Complex refractive indices and single-scattering albedo of global dust aerosols in the shortwave spectrum and relationship to size and iron content.* Atmos. Chem. Phys. 19, 15503-15531. `10.5194/acp-19-15503-2019` | read -- band-1 indices, Table 4 |
| `rochalima2018-fennec-saharan-dust.pdf` | Rocha-Lima et al., 11 authors (2018). *A detailed characterization of the Saharan dust collected during the Fennec campaign in 2011: in situ ground-based and laboratory measurements.* Atmos. Chem. Phys. 18, 1023-1043. `10.5194/acp-18-1023-2018` | read -- band-2 k, **digitised from Fig. 10**; the title does not mention refractive indices but the spectral k is inside |
| `balkanski2007-mineral-dust-radiative-forcing.pdf` | Balkanski, Schulz, Claquin, Guibert (2007). *Reevaluation of Mineral aerosol radiative forcings suggests a better agreement with satellite and AERONET data.* Atmos. Chem. Phys. 7, 81-95. `10.5194/acp-7-81-2007` | read -- source of the size distribution used; tabulates no indices |
| `hess1998-opac.pdf` | Hess, Koepke, Schult (1998). *Optical Properties of Aerosols and Clouds: The Software Package OPAC.* Bull. Amer. Meteor. Soc. 79(5), 831-844. `10.1175/1520-0477(1998)079<0831:OPOAAC>2.0.CO;2` | read -- validation target for `exoplasim/scripts/mie_dust.py` |
| `dibiagio2017-dust-longwave-refractive-index.pdf` | Di Biagio et al., 16 authors (2017). *Global scale variability of the mineral dust long-wave refractive index: a new dataset of in situ measurements for climate modeling and remote sensing.* Atmos. Chem. Phys. 17, 1901-1929. `10.5194/acp-17-1901-2017` | read -- **its k is pinned constant below 6 um by the authors' own statement; use the n, not the k** |
| `dibiagio2014-african-dust-ir-refractive-index.pdf` | Di Biagio, Boucher, Caquineau, Chevaillier, Cuesta, Formenti (2014). *Variability of the infrared complex refractive index of African mineral dust: experimental estimation and implications for radiative transfer and satellite remote sensing.* Atmos. Chem. Phys. 14, 11093-11116. `10.5194/acp-14-11093-2014` | held -- 2.5-25 um, k genuinely retrieved, but pellet method |
| `deschutter2022-gobi-dust-optical-properties.pdf` | Deschutter (2022). *Propriétés optiques des poussières désertiques de Gobi et de ses composés purs.* Thèse, Université de Lille. NNT 2022ULILR018, HAL `tel-03917567` | read -- **dead end for 0.95-4 um**: no tabulated indices anywhere, and that band is literature-seeded |
| `volz1973-ir-optical-constants.pdf` | Volz (1973). *Infrared Optical Constants of Ammonium Sulfate, Sahara Dust, Volcanic Pumice, and Flyash.* Appl. Opt. 12(3), 564-568. `10.1364/AO.12.000564` | held -- the longwave basis of OPAC "mineral transported" |
| `kong2024-dust-refractive-index-uncertainty.pdf` | Kong, Wang, Bi (2024). *Uncertainties in laboratory-measured shortwave refractive indices of mineral dust aerosols and derived optical properties: a theoretical assessment.* Atmos. Chem. Phys. 24, 6911-6935. `10.5194/acp-24-6911-2024` | held -- quantifies the error the sphere assumption injects, which is what `mie_dust.py` assumes |

## Hydrology: the Earth comparator

| file | citation | status |
| --- | --- | --- |
| `wang2018-endorheic-water-storage.pdf` | Wang et al. (2018). *Recent global decline in endorheic basin water storages.* Nature Geoscience 11, 926-932. `10.1038/s41561-018-0265-7` | read -- source of "Earth is one-fifth endorheic", which corrected the 13% this project had been using. 31.8 million km2 across 48,813 landlocked watersheds in 15 arcsec HydroSHEDS. Measures catchment area draining internally, the same quantity as our own endorheic share. |

---

# Derived surface classes

Fetched 2026-08-16 for the derived-surface-class work: deriving diatomite,
duricrust, deflation armour and loess as products of surface process history
rather than as rock types Orogen could know about. See `notes/dust.md` and
`pedology/analysis/phosphorus_budget.json`.

## Bodélé Depression and biogenic lacustrine phosphorus

The exemplar the whole aeolian phosphorus return leg is modelled on.

| file | citation | status |
| --- | --- | --- |
| `hudsonedwards2014-bodele-phosphorus-speciation.pdf` | Hudson-Edwards, Bristow, Cibin, Mason, Peacock (2014). *Solid-phase phosphorus speciation in Saharan Bodélé Depression dusts and source sediments.* Chemical Geology 384, 16-26. `10.1016/j.chemgeo.2014.06.014` | **read** -- Table 1 is the load-bearing one: diatomite 600-610 ppm P, emitted dust 590-930, **aeolian sand 40**. Only 2-4% water soluble. The authigenic apatite is fish bone and scale, "the first-ever report of fish material in aeolian dust" |
| `bristow2010-fertilizing-amazon-dust.pdf` | Bristow, Hudson-Edwards, Chappell (2010). *Fertilizing the Amazon and equatorial Atlantic with West African dust.* Geophys. Res. Lett. 37(14), L14807. `10.1029/2010GL043486` | held -- 0.12 Tg P/yr exported |
| `bristow2009-deflation-bodele.pdf` | Bristow, Drake, Armitage (2009). *Deflation in the dustiest place on Earth: The Bodélé Depression, Chad.* Geomorphology 105(1-2), 50-58. `10.1016/j.geomorph.2007.12.014` | held |
| `bristow2018-auto-abrasion-diatomite.pdf` | Bristow, Moller (2018). *Testing the auto-abrasion hypothesis for dust production using diatomite dune sediments from the Bodélé Depression in Chad.* Sedimentology 65(4), 1322-1330. `10.1111/sed.12423` | held -- the mechanism for WHY diatomite deflates so readily: saltating grains shattering each other |
| `washington2005-bodele-low-level-jet.pdf` | Washington, Todd (2005). *Atmospheric controls on mineral dust emission from the Bodélé Depression, Chad: The role of the low level jet.* Geophys. Res. Lett. 32(17), L17701. `10.1029/2005GL023597` | held -- the low-level jet is the reason this particular basin exports so much; the mechanism is meteorological, not just geological |
| `washington2003-dust-storm-source-areas.pdf` | Washington, Todd, Middleton, Goudie (2003). *Dust-Storm Source Areas Determined by the Total Ozone Monitoring Spectrometer and Surface Observations.* Annals AAG 93(2), 297-313. `10.1111/1467-8306.9302003` | held -- the printed title says "Monitoring"; the instrument is Total Ozone MAPPING Spectrometer. The error is the publisher's, reproduce it |
| `yu2020-disproving-bodele-amazon.pdf` | Yu et al. (2020). *Disproving the Bodélé Depression as the Primary Source of Dust Fertilizing the Amazon Rainforest.* Geophys. Res. Lett. 47(13), e2020GL088020. `10.1029/2020GL088020` | held -- **the premise is contested**; adopt or reject it deliberately rather than by inheritance |
| `chadwick1999-changing-nutrient-sources.pdf` | Chadwick, Derry, Vitousek, Huebert, Hedin (1999). *Changing sources of nutrients during four million years of ecosystem development.* Nature 397, 491-497. `10.1038/17276` | held -- the empirical demonstration that aeolian P sustains an ecosystem once bedrock P is exhausted |

**The reframe this set forces.** Diatomite at ~600 ppm sits essentially at this
world's land mean of 628. The Bodélé is not P-rich *material*; it is an enormous,
exceptionally deflatable source of ordinary-P material. Its export is a mass-flux
result, not a concentration result, and the 15x contrast is against aeolian sand.

## Dust emission physics and supply limitation

| file | citation | status |
| --- | --- | --- |
| `kok2012-physics-windblown-sand-dust.pdf` | Kok, Parteli, Michaels, Bou Karam (2012). *The physics of wind-blown sand and dust.* Rep. Prog. Phys. 75(10), 106901. `10.1088/0034-4885/75/10/106901` | held -- **covers Earth AND Mars**, so the only saltation-threshold review treating non-Earth gravity. Directly relevant at 12.81 m/s2 |
| `marticorena1995-dust-emission-scheme.pdf` | Marticorena, Bergametti (1995). *Modeling the atmospheric dust cycle: 1. Design of a soil-derived dust emission scheme.* JGR 100(D8), 16415-16430. `10.1029/95JD00690` | held -- the standard scheme |
| `shao2000-threshold-friction-velocity.pdf` | Shao, Lu (2000). *A simple expression for wind erosion threshold friction velocity.* JGR 105(D17), 22437-22443. `10.1029/2000JD900304` | held |
| `kok2014-improved-dust-emission-model.pdf` | Kok et al. (2014). *An improved dust emission model - Part 1: Model description and comparison against measurements.* Atmos. Chem. Phys. 14(23), 13023-13041. `10.5194/acp-14-13023-2014` | held |
| `bullard2011-preferential-dust-sources.pdf` | Bullard et al. (2011). *Preferential dust sources: A geomorphological classification designed for use in global dust-cycle models.* JGR 116(F4), F04034. `10.1029/2011JF002061` | held -- classifies sources by supply and availability regime, written for global models |
| `macpherson2008-supply-limited-dust-emission.pdf` | Macpherson, Nickling, Gillies, Etyemezian (2008). *Dust emissions from undisturbed and disturbed supply-limited desert surfaces.* JGR Earth Surf. 113(F2), F02S04. `10.1029/2007JF000800` | held -- **not "Macpherson and King"**, a citation that does not exist |
| `kocurek1999-aeolian-sediment-state.pdf` | Kocurek, Lancaster (1999). *Aeolian system sediment state: theory and Mojave Desert Kelso dune field example.* Sedimentology 46(3), 505-515. `10.1046/j.1365-3091.1999.00227.x` | held -- where supply / availability / transport-capacity is formalised |

## Desert pavement

| file | citation | status |
| --- | --- | --- |
| `wells1995-cosmogenic-stone-pavements.pdf` | Wells, McFadden, Poths, Olinger (1995). *Cosmogenic 3He surface-exposure dating of stone pavements: Implications for landscape evolution in deserts.* Geology 23(7), 613-616. `10.1130/0091-7613(1995)023<0613:CHSEDO>2.3.CO;2` | **read** -- and it **inverts the intuitive model**. Cosmogenic ages show the clasts never moved: *"stone pavements are born at the surface."* Deflation lag, water winnowing and shrink-swell are all rejected. Dust is trapped BENEATH the clast mosaic into a cumulic Av horizon, so the surface is a net dust SINK under a non-erodible cover. Citing this as evidence that deflation armours a surface inverts the paper |
| `mcfadden1987-desert-pavement-origin.pdf` | McFadden, Wells, Jercinovich (1987). *Influences of eolian and pedogenic processes on the origin and evolution of desert pavements.* Geology 15(6), 504-508. `10.1130/0091-7613(1987)15<504:IOEAPP>2.0.CO;2` | held -- proposes the model Wells 1995 tests. Title is "Influences" plural; a widely-copied citing reference list prints the singular |
| `pelletier2007-desert-pavement-dynamics.pdf` | Pelletier, Cline, DeLong (2007). *Desert pavement dynamics: numerical modeling and field-based calibration.* Earth Surf. Process. Landforms 32(13), 1913-1927. `10.1002/esp.1500` | held -- a numerical model rather than a concept |
| `haff1996-pavement-healing-disturbance.pdf` | Haff, Werner (1996). *Dynamical Processes on Desert Pavements and the Healing of Surficial Disturbances.* Quaternary Research 45(1), 38-46. `10.1006/qres.1996.0004` | held -- re-formation timescale after disturbance |

## Duricrusts

| file | citation | status |
| --- | --- | --- |
| `watson1983-gypsum-crusts-deserts-vol1.pdf` + `watson1983-gypsum-crusts-deserts-vol2.pdf` | Watson (1983). *The origin, nature and distribution of gypsum crusts in deserts.* D.Phil. thesis, Oxford. `10.5287/ora-wv9z40k84` | **read** -- **gypcrete forms below ~250 mm/yr with potential evaporation exceeding precipitation every month**, and about 22% of Earth's land lies within those bounds. Explicitly: gypsum crusts occupy the driest zones while calcretes take the less arid parts |
| `alonsozarza2003-palustrine-carbonates-calcretes.pdf` | Alonso-Zarza (2003). *Palaeoenvironmental significance of palustrine carbonates and calcretes in the geological record.* Earth-Sci. Rev. 60(3-4), 261-298. `10.1016/S0012-8252(02)00106-X` | **read** -- **calcrete favoured below 500-600 mm/yr, optimum 100-500, upper bound arguably 1000; lower bound as low as 50.** Mineralogy-keyed bands in section 3.6.1 |
| `machette1985-calcic-soils-sw-usa.pdf` | Machette (1985). *Calcic soils of the southwestern United States.* GSA Special Paper 203, 1-22. `10.1130/SPE203-p1` | held -- Stages I-VI of carbonate morphology |
| `bachman1977-calcic-soils-calcretes-sw-usa.pdf` | Bachman, Machette (1977). *Calcic soils and calcretes in the southwestern United States.* USGS Open-File Report 77-794. `10.3133/ofr77794` | **read** -- carbonate accumulation **0.22 to 0.51 g/cm2/kyr**, and the two-sided constraint: too dry and solutions never infiltrate, too wet and carbonate leaches out |
| `nash2011-desert-crusts-rock-coatings.pdf` | Nash (2011). *Desert Crusts and Rock Coatings*, ch. 8, 131-180, in Thomas (ed.), Arid Zone Geomorphology, 3rd edn. `10.1002/9780470710777.ch8` | held -- the only chapter treating calcrete, gypcrete and silcrete together |
| `ullyott2016-pedogenic-nonpedogenic-silcretes.pdf` | Ullyott, Nash (2016). *Distinguishing pedogenic and non-pedogenic silcretes in the landscape and geological record.* Proc. Geol. Assoc. 127(3), 311-319. `10.1016/j.pgeola.2016.03.001` | held |
| `fenske2025-duricrust-water-table.html` | Fenske, Braun, Guillocheau, Robin (2025). *A numerical model for duricrust formation by water table fluctuations.* Earth Surf. Dynam. 13(1), 119-146. `10.5194/esurf-13-119-2025` | held -- HTML, the publisher's TLS chain was broken. States the ordering compactly and, usefully, that silcrete cannot be climatically calibrated |

**Silcrete has no climatic window and should be dropped from any classifier**
rather than given a fabricated threshold. Nash and Ullyott call identifying its
environmental parameters "highly problematic"; Fenske et al. say climatic
calibration "will not provide duricrust formation boundaries" for it.

## Aeolian phosphorus deposition, and loess

| file | citation | status |
| --- | --- | --- |
| `okin2004-dust-phosphorus-terrestrial.pdf` | Okin, Mahowald, Chadwick, Artaxo (2004). *Impact of desert dust on the biogeochemistry of phosphorus in terrestrial ecosystems.* Global Biogeochem. Cycles 18(2), GB2005. `10.1029/2003GB002145` | held -- deposition to LAND, which is our return leg specifically rather than the ocean |
| `mahowald2008-global-phosphorus-deposition.pdf` | Mahowald et al. (2008). *Global distribution of atmospheric phosphorus sources, concentrations and deposition rates, and anthropogenic impacts.* Global Biogeochem. Cycles 22(4), GB4026. `10.1029/2008GB003240` | held -- mineral aerosols are 82% of total P globally |
| `myriokefalitakis2016-bioavailable-phosphorus-ocean.pdf` | Myriokefalitakis, Nenes, Baker, Mihalopoulos, Kanakidou (2016). *Bioavailable atmospheric phosphorous supply to the global ocean: a 3-D global modeling study.* Biogeosciences 13(24), 6519-6543. `10.5194/bg-13-6519-2016` | held -- 1.300 Tg-P/yr total deposited, 0.455 dissolved, a 35% soluble fraction |
| `nenes2011-atmospheric-acidification-phosphorus.pdf` | Nenes, Krom, Mihalopoulos, Van Cappellen, Shi, Bougiatioti, Zarmpas, Herut (2011). *Atmospheric acidification of mineral aerosols: a source of bioavailable phosphorus for the oceans.* Atmos. Chem. Phys. 11(13), 6265-6272. `10.5194/acp-11-6265-2011` | held -- the acid-processing mechanism itself. Cited by Hudson-Edwards for the finding that acid treatment releases 81-96% of dust P against 2-4% in circum-neutral water, which is why the Bodélé soluble fraction is a floor |
| `stockdale2016-acid-processing-mineral-dusts.pdf` | Stockdale, Krom, Mortimer, Benning, Carslaw, Herbert, Shi, Myriokefalitakis, Kanakidou, Nenes (2016). *Understanding the nature of atmospheric acid processing of mineral dusts in supplying bioavailable phosphorus to the oceans.* PNAS 113(51), 14639-14644. `10.1073/pnas.1608136113` | held -- apatite is 79-96% of dust P, which bounds how much of a dust-borne P subsidy can ever be bioavailable |
| `herbert2018-acid-processing-bioavailable-p.pdf` | Herbert et al. (2018). *The Effect of Atmospheric Acid Processing on the Global Deposition of Bioavailable Phosphorus From Dust.* Global Biogeochem. Cycles 32(9), 1367-1385. `10.1029/2018GB005880` | held -- bioavailable fraction rises from ~10% labile to a 22% ocean mean after acid processing |
| `bettis2003-last-glacial-loess-usa.pdf` | Bettis, Muhs, Roberts, Wintle (2003). *Last Glacial loess in the conterminous USA.* Quat. Sci. Rev. 22(18-19), 1907-1946. `10.1016/S0277-3791(03)00169-0` | held -- mass accumulation rates to 17,500 g/m2/yr, many near-source above 1500. **Last-glacial mid-continent North America, the high end of the global range, not a global typical** |
| `muhs2013-geologic-records-dust-quaternary.pdf` | Muhs (2013). *The geologic records of dust in the Quaternary.* Aeolian Research 9, 3-48. `10.1016/j.aeolia.2012.08.001` | held -- the broad Quaternary dust-deposition-rate compilation |

## Major-ion weathering by lithology

Needed because **Hartmann et al. (2014) does not supply it.** Section 2.1 defines
the weathering rate as the fluvial export of *total* Ca + Mg + Na + K + SiO2 plus
carbonate CO3 -- the cations are summed before the model runs and never separated
again -- and evaporite dissolution beyond carbonate is explicitly excluded. Its
Table A1-2 looks like the answer and is a trap: bulk rock composition in weight
percent, not release.

| file | citation | status |
| --- | --- | --- |
| `meybeck1987-global-chemical-weathering.pdf` | Meybeck (1987). *Global chemical weathering of surficial rocks estimated from river dissolved loads.* Am. J. Sci. 287(5), 401-428. `10.2475/ajs.287.5.401` | **read** -- **the per-lithology major-ion source.** Table 2C gives concentrations and Table 5 release fractions for all eight species across 10 rock types. Two limits: it is a **scanned image with no text layer**, and Table 5 covers **exorheic continents only**, awkward for a 76%-endorheic world |
| `hartmann2011-japan-silicate-weathering-phosphorus.pdf` | Hartmann, Moosdorf (2011). *Chemical weathering rates of silicate-dominated lithological classes and associated liberation rates of phosphorus on the Japanese Archipelago - Implications for global scale analysis.* Chem. Geol. 287(3-4), 125-157. `10.1016/j.chemgeo.2010.12.004` | held -- the upstream calibration for Hartmann 2014, and the one paper that might carry cations ON GLiM CLASSES, which would save the Meybeck mapping. Worth opening first |
| `gaillardet1999-silicate-weathering-co2.pdf` | Gaillardet, Dupré, Louvat, Allègre (1999). *Global silicate weathering and CO2 consumption rates deduced from the chemistry of large rivers.* Chem. Geol. 159(1-4), 3-30. `10.1016/S0009-2541(99)00031-5` | held -- river-basin endmembers, not a lithological yield table |
| `bluth1994-lithologic-climatologic-river-chemistry.pdf` | Bluth, Kump (1994). *Lithologic and climatologic controls of river chemistry.* Geochim. Cosmochim. Acta 58(10), 2341-2359. `10.1016/0016-7037(94)90015-9` | held -- 101 monolithologic rivers, but **bicarbonate and silica only**: two of the eight species |
| `hartmann2012-terrestrial-surface-composition.pdf` | Hartmann, Dürr, Moosdorf, Meybeck, Kempe (2012). *The geochemical composition of the terrestrial surface (without soils) and comparison with the upper continental crust.* Int. J. Earth Sci. 101(1), 365-376. `10.1007/s00531-010-0635-x` | held -- the source of Hartmann 2014's Table A1-2 compositions. Composition, still not release |
| `hartmannmoosdorf2012-glim-lithological-map.pdf` | Hartmann, Moosdorf (2012). *The new global lithological map database GLiM: A representation of rock properties at the Earth surface.* Geochem. Geophys. Geosyst. 13(12), 2012GC004370. `10.1029/2012GC004370` | held -- the source of our lithology classes. A map, not a chemistry dataset |
