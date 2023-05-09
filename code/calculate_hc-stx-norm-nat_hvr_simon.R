#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(stringr)

## Dependencies
# ICC volume and ScaleFactors
f_icc_scales    <- here("data/derivatives/simon_icc_scale.csv")
icc_scales      <- fread(f_icc_scales)
rm(f_icc_scales)

# Extracted HC CSF volumes
f_volumes       <- here("data/derivatives/simon_volumes_hcvc_cnn.csv")
volumes         <- fread(f_volumes)
rm(f_volumes)

# Parse SCANs from IDs
volumes[, SCAN  := str_extract(ID, "(?<=SIMON_).*(?=_t1)")]
volumes[, ID    := NULL]

## Remove Failed QC segmentation:
# stx2_SIMON_ses-048_acq-transverse_run-3_t1_hcvc-cnn.mnc
fail_seg        <- "ses-048_acq-transverse_run-3"
volumes         <- volumes[!SCAN == fail_seg]

## Merge ICC and volumes
volumes         <- icc_scales[volumes, on = "SCAN"]

## Bring back to native scale
# To convert stx volumes to native space DIVIDE by SCALEFACTOR
volumes[, `:=`(HC_l_nat     = LHC         / SCALEFACTOR,
               HC_r_nat     = RHC         / SCALEFACTOR,
               CSF_l_nat    = LCSF        / SCALEFACTOR,
               CSF_r_nat    = RCSF        / SCALEFACTOR
               )]

## Calculate HC_norm and HVR
volumes[, `:=`(HC_l_norm    = HC_l_nat    / ICC,
               HC_r_norm    = HC_r_nat    / ICC,
               HVR_l        = HC_l_nat    / (HC_l_nat + CSF_l_nat),
               HVR_r        = HC_r_nat    / (HC_r_nat + CSF_r_nat)
               )]

## Average two sides
volumes[, `:=`(HC_mean      = (HC_l_nat   + HC_l_nat  ) / 2,
               HC_stx_mean  = (LHC        + RHC       ) / 2,
               HC_norm_mean = (HC_l_norm  + HC_r_norm ) / 2,
               HVR_mean     = (HVR_l      + HVR_r     ) / 2
               )]

## Export RDS
volumes         <- volumes[, .(SCAN, ICC, SCALEFACTOR,
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
                               HVR_r        = HVR_r,
                               HVR_mean
                               )]

volumes |>
  write_rds(here("data/rds/simon_volumes_hc-stx-norm-nat_hvr.rds"))
