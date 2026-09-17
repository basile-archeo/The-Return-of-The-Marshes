# The return of the marshes — Lucca plain (1st–6th c. AD)
Reproducibility materials for:

Basile, S. (2026). The return of the marshes: spatial modelling of settlement
dynamics and hydraulic equilibria in the Lucca plain (Italy, 1st–6th centuries AD).
*Journal of Mediterranean Archaeology* 39(2).

This repository contains the R code, the derived susceptibility model and the
spatial layers needed to reproduce the spatial-statistical results reported in
the paper and in its supplementary note: Monte Carlo simulation of mean flood
susceptibility at site locations, inhomogeneous Poisson point-process (PPM)
models, and the sampling-effort sensitivity analysis.

## Contents

```
.
├── R_marshes_supplementary.R   # analysis script (Monte Carlo, PPM, effort-weighted null)
├── data/
│   ├── modello_paludi.tif       # flood-susceptibility raster, values 1–10 (EPSG:3003)
│   ├── maschera.*               # study-area coverage mask
│   ├── mura_lucca.*             # Lucca city walls (closed line, urban/rural split)
│   ├── interventi_archeologici.* # all 427 archaeological interventions
│   ├── I d.C. piana_valle.*     # rural sites, phase I  (1st c. AD)
│   ├── III d.C. piana_valle.*   # rural sites, phase II–III (2nd–3rd c. AD)
│   ├── IV.*                     # rural sites, phase IV (4th c. AD)
│   └── V_VI.*                   # rural sites, phase V–VI (5th–6th c. AD)
├── CITATION.cff
└── README.md
```

* Each shapefile is a set of sidecar files sharing one name
(`.shp`, `.shx`, `.dbf`, `.prj`, `.cpg`) — keep them together.

## How to reproduce

1. Install R (≥ 4.2) and the three dependencies:

   ```r
   install.packages(c("sf", "terra", "spatstat"))
   ```

2. Put all the files from `data/` in one working folder, then in the script set
   the working directory to that folder and leave `shp_dir <- "."`. All layers
   are read relative to that single folder.

3. Run `R_marshes_supplementary.R`. It prints a results table and the summary
   lines used in the paper. Runtime is a few minutes (1,000 Monte Carlo
   simulations per phase, uniform and effort-weighted).

The random seed is fixed (`set.seed(42)`). Monte Carlo p-values carry
simulation noise of a few thousandths and may vary slightly with a different
seed; every phase-level conclusion is stable across runs.

### Expected output

| phase            |  n  | PPM coef. | PPM p | MC p (uniform) | MC p (effort-weighted) |
|------------------|:---:|:---------:|:-----:|:--------------:|:----------------------:|
| I (1st c. AD)    | 111 |   +0.359  | <0.05 |     1.000      |         1.000          |
| II–III (2–3 c.)  |  42 |   +0.011  | n.s.  |     ~0.52      |         ~0.30          |
| IV (4th c. AD)   |  28 |   −0.217  | <0.05 |     ~0.037     |         ~0.007         |
| V–VI (5–6 c.)    |  18 |   −0.354  | <0.05 |     ~0.002     |         ~0.000         |

A positive PPM coefficient (and a Monte Carlo p near 1) indicates preference for
flood-prone ground; a negative coefficient (and p near 0) indicates avoidance.

## Data sources and attribution

The `modello_paludi` susceptibility raster, the `maschera` mask and the phase
site layers are derived products created by the author. The composite model was
built from open regional base data (DEM from LiDAR, soil and geological maps)
published by **Regione Toscana** via the Geoscopio geoportal
(https://www502.regione.toscana.it/geoscopio/), used under their open-data
terms. The `interventi_archeologici` layer derives from the MAPPA/MAGOH
database (University of Pisa); the wider dataset is documented at
https://digitallib.unipi.it/it/raccolta/Lucca-e-lager-lucensis-dalleta-tardo-repubblicana-al-tardoantico-le-trasformazioni-di-una-citta-e-del-suo-territorio/

Please cite the article (and this repository's DOI, once minted) when reusing
these materials.

## Licence

- **Code** (`R_marshes_supplementary.R`), Data and derived layers** (rasters, shapefiles): 
Creative Commons Attribution 4.0 International (CC BY 4.0),
  https://creativecommons.org/licenses/by/4.0/

Author: Salvatore Basile, MAPPA Laboratory, Dipartimento di Civiltà e Forme del
Sapere, University of Pisa.

