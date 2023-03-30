#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(lubridate)

## Dependencies
# ICC volume and ScaleFactors
icc_scales      <- here("data/derivatives/adni_icc_scale.csv") |> fread()
icc_scales[, SCANDATE := ymd(SCANDATE)]

# Curated HC CSF volumes
f_volumes       <- here("data/rds/adni-bl_volumes_hcvc.rds")
if (file.exists(f_volumes)) {
  volumes       <- read_rds(f_volumes)
  rm(f_volumes)
} else {
  s_volumes     <- here("code/qc_segmentations_adni-bl.R")
  source(s_volumes)
  rm(f_volumes, s_volumes, volumes_hcvcag)
}

# Merge ICC and volumes
volumes     <- icc_scales[volumes, on = .(PTID, SCANDATE)]

# Bring back to native scale
# To convert stx volumes to native space DIVIDE by SCALEFACTOR
volumes[, `:=`(HC_l_nat     = LHC         / SCALEFACTOR,
               HC_r_nat     = RHC         / SCALEFACTOR,
               CSF_l_nat    = LCSF        / SCALEFACTOR,
               CSF_r_nat    = RCSF        / SCALEFACTOR)]

# Calculate HC_norm and HVR
volumes[, `:=`(HC_l_norm    = HC_l_nat    / ICC,
               HC_r_norm    = HC_r_nat    / ICC,
               HVR_l        = HC_l_nat    / (HC_l_nat + CSF_l_nat),
               HVR_r        = HC_r_nat    / (HC_r_nat + CSF_r_nat))]

# Average two sides
volumes[, `:=`(HC_mean      = (HC_l_nat   + HC_l_nat  ) / 2,
               HC_stx_mean  = (LHC        + RHC       ) / 2,
               HC_norm_mean = (HC_l_norm  + HC_r_norm ) / 2,
               HVR_mean     = (HVR_l      + HVR_r     ) / 2)]

# Export RDS
volumes     <- volumes[, .(PTID, SCANDATE, ICC, SCALEFACTOR,
                           METHOD       = METHOD,
                           HC_r         = HC_r_nat,
                           HC_l         = HC_l_nat,
                           HC_mean,
                           HC_stx_r     = RHC,
                           HC_stx_l     = LHC,
                           HC_stx_mean,
                           HC_norm_r    = HC_r_norm,
                           HC_norm_l    = HC_l_norm,
                           HC_norm_mean,
                           HVR_l        = HVR_l,
                           HVR_r        = HVR_l,
                           HVR_mean
                           )]

volumes |>
  write_rds(here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds"))
