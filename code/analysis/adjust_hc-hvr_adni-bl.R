#!/usr/bin/env Rscript

library(here)
library(readr)
library(glue)
library(data.table)
library(lubridate)

## Dependencies
# ICC volume and ScaleFactors
fpath           <- here("data/derivatives/adni_icc_scale.csv")
if (!file.exists(fpath)) glue("File: {fpath} ",
                              "is required but could not be found.") |> stop()
icc_scales      <- fread(fpath)
icc_scales[, SCANDATE := ymd(SCANDATE)]

# Curated HC CSF volumes
fpath       <- here("data/rds/adni-bl_volumes_hcvc.rds")
if (file.exists(fpath)) {
  volumes       <- read_rds(fpath)
} else {
  here("code/data_parsing/qc_segmentations_adni-bl.R") |> source()
  rm(volumes_hcvcag) # Unused
}

fpath    <- here("data/rds/adni-bl_volumes_freesurfer.rds")
if (file.exists(fpath)) {
  fs_vols       <- read_rds(fpath)
} else {
  here('code/data_parsing/parse_freesurfer-vols.R') |> source()
}

## Controls
fpath    <- here("data/rds/adnimerge_baseline.rds")
if (file.exists(fpath)) {
  adnimerge     <- read_rds(fpath)
} else {
  here('code/data_parsing/parse_adnimerge-bl.R') |> source()
}

controls    <- adnimerge[DX == "CH", PTID]
rm(adnimerge, fpath)

# Remove failed segmentations
fs_vols[!is.na(FS_house), QC := "Pass"]
fs_vols[is.na(QC), QC := "Fail"]

# Merge ICC and volumes
volumes     <- rbindlist(list(volumes,
                              fs_vols[, .(LHC, RHC, HC, LCSF, RCSF, CSF, QC,
                                          #ICC_fs = BRAIN, # use this??
                                          PTID, SCANDATE = ymd(SCANDATE),
                                          METHOD = "fs6")]),
                         fill = TRUE)
rm(fs_vols)

volumes     <- icc_scales[volumes, on = .(PTID, SCANDATE)]
rm(icc_scales)
volumes[, ICC := ICC / 1000]
icc_mean_cn <- volumes[QC == "Pass"
                       ][PTID %in% controls, .(ICC_cn = mean(ICC)), METHOD]

# Bring back to native scale
# To convert stx volumes to native space DIVIDE by SCALEFACTOR
# Scale everything to cm^3
volumes[QC == "Pass" & METHOD == 'fs6',
        `:=`(HC_l_nat     = LHC               / 1000,
             HC_r_nat     = RHC               / 1000,
             CSF_l_nat    = LCSF              / 1000,
             CSF_r_nat    = RCSF              / 1000,
             HC_l_stx     = LHC * SCALEFACTOR / 1000,
             HC_r_stx     = RHC * SCALEFACTOR / 1000)]

volumes[QC == "Pass" & METHOD != 'fs6',
        `:=`(HC_l_nat     = LHC   / (SCALEFACTOR * 1000),
             HC_r_nat     = RHC   / (SCALEFACTOR * 1000),
             CSF_l_nat    = LCSF  / (SCALEFACTOR * 1000),
             CSF_r_nat    = RCSF  / (SCALEFACTOR * 1000),
             HC_l_stx     = LHC   /                1000,
             HC_r_stx     = RHC   /                1000)]

## Regression slopes for PCP & Residual normalizations
# Use only Controls for the models
volumes_lng <- volumes[QC == "Pass" & PTID %in% controls,
                       .(METHOD, ICC,
                         HC_l_nat, HC_r_nat,
                         CSF_l_nat, CSF_r_nat)] |>
  melt(id.vars = c("METHOD", "ICC"),
       variable.name = "ROI",
       value.name = "VAL")

rm(controls)

volumes_lng[, ROI := stringr::str_remove(ROI, "_nat")]
#volumes_lng[, ROI := stringr::str_to_lower(ROI)]

# Power-corrected proportion:
# VOL_adj = VOL / ICC ** b
# b: slope of log(VOL) ~ log(ICC)
b_pcp       <- volumes_lng[, summary(lm(log(VAL) ~ log(ICC)))$coefficients[2],
                           .(ROI, METHOD)] |>
              dcast(METHOD ~ ROI, value.var = "V1")
setnames(b_pcp, names(b_pcp)[-1], paste0(names(b_pcp)[-1], "_b_pcp"))

# Residuals
# Remove the residuals from VOL ~ ICC regression
# VOL_adj = VOL - b(ICC - ICC_cn)
b_res       <- volumes_lng[, summary(lm(VAL ~ ICC))$coefficients[2],
                           .(ROI, METHOD)] |>
              dcast(METHOD ~ ROI, value.var = "V1")
setnames(b_res, names(b_res)[-1], paste0(names(b_res)[-1], "_b_res"))

# Add slopes and ICC_cn to volumes
volumes     <- volumes[b_pcp, on = "METHOD"
                       ][b_res, on = "METHOD"
                       ][icc_mean_cn, on = "METHOD"]
rm(volumes_lng, icc_mean_cn, b_pcp, b_res)

## Apply adjustment methods
# HC & CSF
volumes[,
        `:=`(HC_l_prop    = HC_l_nat  / ICC,
             HC_r_prop    = HC_r_nat  / ICC,
             HC_l_pcp     = HC_l_nat  / ICC ** HC_l_b_pcp,
             HC_r_pcp     = HC_r_nat  / ICC ** HC_r_b_pcp,
             CSF_l_pcp    = CSF_l_nat / ICC ** CSF_l_b_pcp,
             CSF_r_pcp    = CSF_r_nat / ICC ** CSF_r_b_pcp,
             HC_l_res     = HC_l_nat  - HC_l_b_res  * (ICC - ICC_cn),
             HC_r_res     = HC_r_nat  - HC_r_b_res  * (ICC - ICC_cn),
             CSF_l_res    = CSF_l_nat - CSF_l_b_res * (ICC - ICC_cn),
             CSF_r_res    = CSF_r_nat - CSF_r_b_res * (ICC - ICC_cn))]

# HVR
volumes[,
        `:=`(HVR_l        = HC_l_nat  / (HC_l_nat + CSF_l_nat),
             HVR_r        = HC_r_nat  / (HC_r_nat + CSF_r_nat),
             HVR_l_pcp    = HC_l_pcp  / (HC_l_pcp + CSF_l_pcp),
             HVR_r_pcp    = HC_r_pcp  / (HC_r_pcp + CSF_r_pcp),
             HVR_l_res    = HC_l_res  / (HC_l_res + CSF_l_res),
             HVR_r_res    = HC_r_res  / (HC_r_res + CSF_r_res))]

## Average two sides
volumes[, `:=`(HC_mean      = (HC_l_nat   + HC_r_nat  ) / 2,
               HC_stx_mean  = (HC_l_stx   + HC_r_stx  ) / 2,
               HC_prop_mean = (HC_l_prop  + HC_r_prop ) / 2,
               HC_pcp_mean  = (HC_l_pcp   + HC_r_pcp  ) / 2,
               HC_res_mean  = (HC_l_res   + HC_r_res  ) / 2,
               HVR_mean     = (HVR_l      + HVR_r     ) / 2,
               HVR_pcp_mean = (HVR_l_pcp  + HVR_r_pcp ) / 2,
               HVR_res_mean = (HVR_l_res  + HVR_l_res ) / 2)]

## Export RDS
volumes     <- volumes[, .(PTID, SCANDATE, ICC, SCALEFACTOR, METHOD,
                           HC_r         = HC_r_nat,
                           HC_l         = HC_l_nat,
                           HC_mean,
                           HC_stx_r     = HC_r_stx,
                           HC_stx_l     = HC_l_stx,
                           HC_stx_mean,
                           HC_prop_r    = HC_r_prop,
                           HC_prop_l    = HC_l_prop,
                           HC_prop_mean,
                           HC_pcp_r    = HC_r_pcp,
                           HC_pcp_l    = HC_l_pcp,
                           HC_pcp_mean,
                           HC_res_r    = HC_r_res,
                           HC_res_l    = HC_l_res,
                           HC_res_mean,
                           HVR_l        = HVR_l,
                           HVR_r        = HVR_r,
                           HVR_mean,
                           HVR_pcp_l    = HVR_l_pcp,
                           HVR_pcp_r    = HVR_r_pcp,
                           HVR_pcp_mean,
                           HVR_res_l    = HVR_l_res,
                           HVR_res_r    = HVR_r_res,
                           HVR_res_mean)]

write_rds(volumes, here("data/rds/adni-bl_volumes_icv-adjusted.rds"))
