#!/usr/bin/env Rscript

library(here)
library(data.table)
library(stringr)
library(lubridate)
library(readr)

### ADNIMERGE
f_adnimerge     <- here("data/rds/adnimerge_baseline.rds")
if (file.exists(f_adnimerge)) {
  adnimerge     <- read_rds(f_adnimerge)
  rm(f_adnimerge)
} else {
  s_adnimerge   <- here("code/parse_adnimerge-bl.R")
  source(s_adnimerge)
  rm(f_adnimerge, s_adnimerge)
}

### ADNI QC Fails
acq_fails       <- fread(here('lists/adni_acquisition_failures.lst'),
                         header = FALSE, col.names = "FAILS")

acq_fails[, `:=`(PTID = str_extract(FAILS, "\\d{3}_S_\\d{4}"),
                 SCANDATE = ymd(str_extract(FAILS, "(?<=\\d{4}_)S*\\d+")))]

### MALF Segmentation QC
malf_qc1        <- fread(here('lists/qrater_malf_2022-12-20.csv'),
                         header = FALSE, col.names = c("ID", "QC"))

malf_qc2        <- fread(here('lists/qrater_malf_reg_fails_2022-12-20.csv'),
                         header = FALSE, col.names = c("ID", "QC"))

malf_qc         <- rbindlist(list(malf_qc1, malf_qc2))

malf_qc_fails   <- malf_qc[QC == "Fail"]
malf_qc_fails[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                     SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)S*\\d+")))]
rm(malf_qc, malf_qc1, malf_qc2)

## HC/VC volumes ## 1746
f_malf_vols     <- here('data/derivatives/adni-bl_volumes_hcvc_malf.csv')
malf_vols       <- fread(f_malf_vols)
rm(f_malf_vols)

malf_vols[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                 SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)\\d+")))]

## Anti-joins
# Dropped by ADNI -> 70
dropped         <- malf_vols[!adnimerge, on = "PTID", PTID]
malf_vols       <- malf_vols[!PTID %in% dropped]    # 1676

firstpass       <- malf_vols[!acq_fails, on = "PTID"]  # Acquisition -> 1641
malf_vols_final <- firstpass[!malf_qc_fails, on = "PTID"] # Segm -> 1638

# Remove useless col and add method
malf_vols_final[, `:=`(ID = NULL, METHOD = "malf")]
rm(malf_vols, dropped, firstpass, malf_qc_fails)

### Non-local Patch-based Segmentation QC
nlpb_qc1        <- fread(here('lists/qrater_nlpb_2022-12-20.csv'),
                         header = FALSE, col.names = c("ID", "QC"))

nlpb_qc2        <- fread(here('lists/qrater_nlpb_reg_fails_2022-12-20.csv'),
                         header = FALSE, col.names = c("ID", "QC"))

nlpb_qc         <- rbindlist(list(nlpb_qc1, nlpb_qc2))

nlpb_qc_fails   <- nlpb_qc[QC == "Fail"]
nlpb_qc_fails[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                     SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)S*\\d+")))]
rm(nlpb_qc, nlpb_qc1, nlpb_qc2)

## HC/VC volumes
f_nlpb_vols     <- here('data/derivatives/adni-bl_volumes_hcvc_nlpb.csv')
nlpb_vols       <- fread(f_nlpb_vols) ## 1746
rm(f_nlpb_vols)

nlpb_vols[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                 SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)\\d+")))]

## Anti-joins
dropped         <- nlpb_vols[!adnimerge, on = "PTID", PTID]
nlpb_vols       <- nlpb_vols[!PTID %in% dropped]  # 1676

firstpass       <- nlpb_vols[!acq_fails, on = "PTID"]       # Acquisition -> 1641
nlpb_vols_final <- firstpass[!nlpb_qc_fails, on = "PTID"] # Segm -> 1638

# Remove useless cols and add method
nlpb_vols_final[, `:=`(ID = NULL, METHOD = "nlpb")]
rm(nlpb_vols, dropped, firstpass, nlpb_qc_fails)

### CNN Segmentations QC
## Only QCed HV/VC segmentations
cnn_qc1         <- fread(here('lists/qrater_cnn1_2022-12-20.csv'),
                         header = FALSE, col.names = c("ID", "QC"))

cnn_qc2         <- fread(here('lists/qrater_cnn1_reg_fails_2022-12-20.csv'),
                         header = FALSE, col.names = c("ID", "QC"))

cnn_qc          <- rbindlist(list(cnn_qc1, cnn_qc2))

cnn_qc_fails    <- cnn_qc[QC == "Fail"]
cnn_qc_fails[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                    SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)S*\\d+")))]
rm(cnn_qc, cnn_qc1, cnn_qc2)

## HC/VC volumes
cnn_vols        <- fread(here('data/derivatives/adni-bl_volumes_hcvc_cnn.csv')) ## 1746

cnn_vols[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)\\d+")))]

## Anti-joins
dropped         <- cnn_vols[!adnimerge, on = "PTID", PTID]    # Dropped by ADNI -> 70
cnn_vols        <- cnn_vols[!PTID %in% dropped]  # 1676

firstpass       <- cnn_vols[!acq_fails, on = "PTID"]       # Acquisition -> 1641
cnn_vols_final  <- firstpass[!cnn_qc_fails, on = "PTID"] # Segm -> 1639

# Remove useless cols and add method
cnn_vols_final[, `:=`(ID = NULL, METHOD = "cnn")]
rm(cnn_vols, dropped, firstpass)

## Merge data.tables
volumes         <- rbindlist(list(malf_vols_final,
                                  nlpb_vols_final,
                                  cnn_vols_final),
                             use.names = TRUE)

rm(malf_vols_final, nlpb_vols_final, cnn_vols_final)


### HC/VC/AMY volumes
cnn_vols2       <- fread(here('data/derivatives/adni-bl_volumes_hcvc-ag_cnn.csv')) ## 1746

cnn_vols2[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                 SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)\\d+")))]

## Anti-joins
dropped         <- cnn_vols2[!adnimerge, on = "PTID", PTID]    # Dropped by ADNI -> 70
cnn_vols2       <- cnn_vols2[!PTID %in% dropped]  # 1676

firstpass       <- cnn_vols2[!acq_fails, on = "PTID"]   # Acquisition -> 1641
volumes_hcvcag  <- firstpass[!cnn_qc_fails, on = "PTID"] # Segm -> 1639
rm(cnn_vols2, dropped, firstpass, cnn_qc_fails)
rm(acq_fails, adnimerge)

# Remove useless cols
volumes_hcvcag[, ID := NULL]

## Export
write_rds(acq_fails, here("data/rds/ptid_qc_discarded.rds"))
write_rds(volumes, here("data/rds/adni-bl_volumes_hcvc.rds"))
write_rds(volumes_hcvcag, here("data/rds/adni-bl_volumes_hcvc-ag.rds"))
