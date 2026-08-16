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

## Geochemistry: the weathering thermostat

| file | citation | status |
| --- | --- | --- |
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
