#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(progress)
library(ggplot2)
library(ggridges)
library(ggsignif)
library(ggtext)
library(gridExtra)

## Calculate and compare linear regressions of HVR and HCv
## ADNI data CN|MCI|AD

# Read RDS objects
adnimerge     <- here("data/rds/adnimerge_baseline.rds") |> read_rds()
volumes       <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds") |>
                read_rds()

# Merge
DT            <- volumes[adnimerge, on = "PTID",
                         .(PTID, METHOD, DX, AGE, PTEDUCAT, PTGENDER, ADAS13,
                           #RAVLT_immediate, RAVLT_learning, RAVLT_forgetting,
                           RAVLT_learning = as.numeric(RAVLT_learning),
                           #HC = HC_mean,     # Native space
                           HCv = HC_stx_mean, # Head-size normalized
                           HVR = HVR_mean)
                         ][DX != ""       &
                           !is.na(HVR)    &
                           !is.na(ADAS13) &
                           !is.na(RAVLT_learning)]



