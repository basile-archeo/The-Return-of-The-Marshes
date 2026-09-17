
#  Supplementary R code
#  "The return of the marshes: spatial modelling of settlement dynamics and
#   hydraulic equilibria in the Lucca plain (1st-6th c. AD)"
#  Journal of Mediterranean Archaeology, 39.2 2026
#  Author: Salvatore Basile, MAPPA Laboratory, Dep. of Civilization and Forms
#  of Knowledge, 
#  University of Pisa (Italy)
#
#  Spatial statistical analysis:
#    (1) Monte Carlo simulation of mean flood susceptibility at site locations
#    (2) Inhomogeneous Poisson point process models (PPM)
#    (3) Sampling-effort sensitivity analysis (effort-weighted Monte Carlo)
#
#  Plotting code (ggplot2 figures) is intentionally omitted: this script
#  produces only the numerical results reported in the paper. 
#  Figures are generated separately.
#
#  Dependencies: sf, terra, spatstat

library(sf)
library(terra)
library(spatstat)      # loads spatstat.geom / .explore / .model

set.seed(42)
N_SIM <- 1000


# 1.  INPUT DATA  (all layers in EPSG:3003)
# All layers are provided in the data/ folder of this repository. The wider
# MAPPA/MAGOH dataset from which the intervention layer derives is documented at:
# https://digitallib.unipi.it/it/raccolta/Lucca-e-lager-lucensis-dalleta-tardo-repubblicana-al-tardoantico-le-trasformazioni-di-una-citta-e-del-suo-territorio/
 

# Set the working directory to the data/ folder (all layers sit together there),
# e.g. setwd("path/to/lucca-marshes-jma/data"); then shp_dir <- "." resolves every file.
setwd("C:/Users/Mappa21/University of Pisa/OneDrive - University of Pisa/Pubblicazioni/pubblicazioni/articolo paludi/revised_new/lucca-marshes-jma/lucca-marshes-jma/data")
shp_dir <- "." 

#upload all the files 

susc     <- rast("modello_paludi.tif")            # susceptibility raster, 1-10
maschera <- st_read("maschera.shp", quiet = TRUE) # study-area coverage mask
walls    <- st_read("mura_lucca.shp", quiet = TRUE)     # Roman city walls of Lucca (line, closed below)
interv   <- st_read("interventi_archeologici.shp", quiet = TRUE)  # all 427 interventions (geometry + Id)
phase_I     <- st_read(file.path(shp_dir, "I d.C. piana_valle.shp"),   quiet = TRUE)
phase_II_III<- st_read(file.path(shp_dir, "III d.C. piana_valle.shp"), quiet = TRUE)
phase_IV    <- st_read(file.path(shp_dir, "IV.shp"),                   quiet = TRUE)
phase_V_VI  <- st_read(file.path(shp_dir, "V_VI.shp"),                 quiet = TRUE)


# phase site layers (rural sites of the plain, one shapefile per phase)
phases <- list(
  "I (1st c. AD)"      = phase_I,
  "II-III (2-3 c. AD)" = phase_II_III,
  "IV (4th c. AD)"     = phase_IV,
  "V-VI (5-6 c. AD)"   = phase_V_VI
)



# 2.  DERIVED OBJECTS: masked susceptibility, spatstat window, cell table

# Clip the raster to the study-area mask. Cells outside the mask become NA;
# the non-NA cells (n = 932847) are the spatial universe for every test.
susc_masked <- mask(susc, vect(maschera))

# spatstat covariate image and observation window.
# Building the window from the non-NA pixels of the masked image avoids any
# polygon-orientation pitfalls and keeps window and covariate perfectly
# co-registered (this also removes the raster-flip risk of manual matrix
# construction).
susc_df <- as.data.frame(susc_masked, xy = TRUE, na.rm = FALSE)
names(susc_df) <- c("x", "y", "z")
susc_im <- as.im(susc_df)          
W       <- as.owin(susc_im)        

# Table of in-mask cells: coordinates + susceptibility value.
cell_id  <- which(!is.na(values(susc_masked)))
cell_xy  <- xyFromCell(susc_masked, cell_id)
cell_val <- as.numeric(values(susc_masked)[cell_id])
cat(sprintf("In-mask cells: %d | mean susceptibility: %.3f\n",
            length(cell_val), mean(cell_val)))


# 3.  susceptibility at a set of site points, dropping points that
#     fall outside raster coverage. 
#     This is the same handling used throughout the analysis, so that n per phase
#     is 111 / 42 / 28 / 18.

site_susc <- function(pts) {
  v <- terra::extract(susc, vect(pts))[, 2]
  v[is.finite(v)]
}


# 4.  MONTE CARLO SIMULATION
#     Compare observed mean susceptibility at site locations against a null
#     distribution of the same number of points redistributed within the
#     study area. p = proportion of simulated means <- observed mean
#     (a p-value near 0 -> sites occupy significantly LOWER susceptibility
#     than chance; a p-value near 1 => significantly HIGHER).
#     `weights` = NULL gives the uniform null used in the paper; a vector of
#     per-cell weights gives the effort-weighted null (section 6).

run_mc <- function(obs_vals, weights = NULL, n_sim = N_SIM) {
  n        <- length(obs_vals)
  obs_mean <- mean(obs_vals)
  sims     <- numeric(n_sim)
  for (s in seq_len(n_sim)) {
    idx      <- sample.int(length(cell_val), n, replace = FALSE, prob = weights)
    sims[s]  <- mean(cell_val[idx])
  }
  list(n = n, obs_mean = obs_mean, p = mean(sims <= obs_mean),
       null_mean = mean(sims))
}


# 5.  INHOMOGENEOUS POISSON POINT PROCESS MODEL
#     First-order model with susceptibility as a continuous covariate.
#     A negative slope => site intensity falls as susceptibility rises
#     (avoidance); a positive slope => preference for susceptible ground.

run_ppm <- function(pts) {
  xy <- st_coordinates(pts)
  ok <- is.finite(terra::extract(susc, vect(pts))[, 2])   # keep in-coverage pts
  pp <- ppp(xy[ok, 1], xy[ok, 2], window = W, checkdup = FALSE)
  fit <- ppm(pp ~ susc, covariates = list(susc = susc_im)) ## susceptibility entered as a continuous first-order covariate
  co  <- summary(fit)$coefs.SE.CI
  est <- co["susc", "Estimate"]; se <- co["susc", "S.E."]
  z   <- est / se
  list(coef = est, se = se,
       ci_lo = est - 1.96 * se, ci_hi = est + 1.96 * se,
       p = 2 * pnorm(-abs(z)))
}


# 6.  SAMPLING-EFFORT SENSITIVITY: effort-weighted Monte Carlo null
#     The uniform null assumes fieldwork was spatially even. To test this,
#     build a kernel-density surface of ACTUAL excavation effort from the
#     rural interventions and draw the null points in proportion to it.
#
#     Rural vs urban: an intervention is urban if its centroid lies inside
#     the (closed) line of the Lucca city walls, rural otherwise. All four
#     phase site-sets lie outside the walls, so the effort surface is built
#     from rural interventions only.

# close the city-walls line into a polygon (endpoints ~7 m apart)
walls_line  <- st_line_merge(st_combine(st_geometry(walls)))
wc          <- st_coordinates(walls_line)[, c("X", "Y")]
if (!all(wc[1, ] == wc[nrow(wc), ])) wc <- rbind(wc, wc[1, ])
walls_poly  <- st_sfc(st_polygon(list(wc)), crs = st_crs(walls))

interv_cent <- st_centroid(st_geometry(interv))
is_urban    <- lengths(st_within(interv_cent, walls_poly)) > 0
rural_cent  <- interv_cent[!is_urban]
rural_in    <- rural_cent[lengths(st_within(rural_cent, st_geometry(maschera))) > 0]
cat(sprintf("Interventions: %d urban (in walls), %d rural; %d rural within mask\n",
            sum(is_urban), sum(!is_urban), length(rural_in)))

# kernel-density effort surface (Scott's-rule bandwidth) and per-cell weights
rxy        <- st_coordinates(rural_in)
rural_ppp  <- ppp(rxy[, 1], rxy[, 2], window = W, checkdup = FALSE)
effort_im  <- density.ppp(rural_ppp, sigma = bw.scott(rural_ppp), #Scott's rule chosen because intervention density is sufficiently high and objective bandwidth selection avoids subjective smoothing.
                          positive = TRUE)
w_cell     <- interp.im(effort_im, cell_xy[, 1], cell_xy[, 2])
w_cell[!is.finite(w_cell) | w_cell < 0] <- 0
w_cell     <- w_cell / sum(w_cell)
cat(sprintf("Effort-weighted mean susceptibility: %.3f (uniform: %.3f)\n",
            sum(w_cell * cell_val), mean(cell_val)))

# ---------------------------------------------------------------------------
# 7.  PRINT THE RESULTS TABLE
# ---------------------------------------------------------------------------

# Could take a few minutes

res <- data.frame()
for (nm in names(phases)) {
  sv   <- site_susc(phases[[nm]])
  mc_u <- run_mc(sv, weights = NULL)     # uniform null (published)
  mc_e <- run_mc(sv, weights = w_cell)   # effort-weighted null (sensitivity)
  pm   <- run_ppm(phases[[nm]])
  res  <- rbind(res, data.frame(
    phase       = nm,
    n           = mc_u$n,
    ppm_coef    = round(pm$coef, 3),
    ppm_ci_lo   = round(pm$ci_lo, 3),
    ppm_ci_hi   = round(pm$ci_hi, 3),
    ppm_p       = signif(pm$p, 3),
    mc_p_unif   = round(mc_u$p, 3),
    mc_p_effort = round(mc_e$p, 3)
  ))
}
print(res, row.names = FALSE)

# For reproducibility, results should closely match those reported in the manuscript.


#  EXPECTED OUTPUT (values reported in the paper and supplementary note)
#
#  phase                 n   ppm_coef   ppm_p    mc_p_unif   mc_p_effort
#  I (1st c. AD)       111     +0.359   <0.05        1.000         1.000
#  II-III (2-3 c. AD)   42     +0.011    n.s.        ~0.52         ~0.30
#  IV (4th c. AD)       28     -0.217   <0.05        ~0.037        ~0.007
#  V-VI (5-6 c. AD)     18     -0.354   <0.05        ~0.002        ~0.000
#
#  HOW TO READ THE MONTE CARLO p-VALUES
#  p is the proportion of simulated means that fall at or below the observed
#  mean, so the two tails have opposite meanings:
#    p close to 1  => sites occupy SIGNIFICANTLY HIGHER susceptibility than
#                     expected by chance (preference for flood-prone ground);
#    p close to 0  => sites occupy SIGNIFICANTLY LOWER susceptibility than
#                     expected by chance (avoidance of flood-prone ground);
#    p around 0.5  => no departure from chance.
#  The PPM coefficient carries the same information independently: positive
#  values indicate preference, negative values avoidance.
#
#  Monte Carlo p-values carry simulation noise of a few thousandths and vary
#  slightly between runs with different seeds; the values above were stable
#  across repeated runs. The effort-weighted column leaves every phase-level
#  conclusion unchanged and, in the two significant late phases, lowers the
#  p-value further: excavation effort is itself very weakly biased towards
#  higher-susceptibility ground, so the effort-weighted null is the more
#  conservative test of late-antique avoidance.
# ---------------------------------------------------------------------------