#!/usr/bin/env Rscript
library(here)
library(data.table)
library(stringr)
library(lubridate)
library(readr)

### ADNIMERGE
fpath         <- here("data/rds/adnimerge_baseline.rds")
if (file.exists(fpath)) {
  adnimerge   <- read_rds(fpath)
} else {
  here("code/parse_adnimerge-bl.R") |> source()
}

### ADNI QC Fails
fpath         <- here('lists/adni_acquisition_failures.lst')
if (!file.exists(fpath)) str_glue("File: {fpath} ",
                                  "is required but could not be found.") |>
                          stop()
acq_fails     <- fread(fpath, header = FALSE, col.names = "FAILS")

acq_fails[, `:=`(PTID = str_extract(FAILS, "\\d{3}_S_\\d{4}"),
                 SCANDATE = ymd(str_extract(FAILS, "(?<=\\d{4}_)S*\\d+")))]

### MALF Segmentation QC
fpath         <- here('lists/qrater_malf_2022-12-20.csv')
if (!file.exists(fpath)) str_glue("File: {fpath} ",
                                  "is required but could not be found.") |>
                          stop()
malf_qc1      <- fread(fpath, header = FALSE, col.names = c("ID", "QC"))

fpath         <- here('lists/qrater_malf_reg_fails_2022-12-20.csv')
if (!file.exists(fpath)) str_glue("File: {fpath} ",
                                  "is required but could not be found.") |>
                          stop()
malf_qc2      <- fread(fpath, header = FALSE, col.names = c("ID", "QC"))

malf_qc       <- rbindlist(list(malf_qc1, malf_qc2))

malf_qc_fails <- malf_qc[QC == "Fail"]
malf_qc_fails[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                     SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)S*\\d+")))]
rm(malf_qc, malf_qc1, malf_qc2)

## HC/VC volumes ## 1746
fpath         <- here('data/derivatives/adni-bl_volumes_hcvc_malf.csv')
if (!file.exists(fpath)) str_glue("File: {fpath} ",
                                  "is required but could not be found.") |>
                          stop()
malf_vols     <- fread(fpath)
malf_vols[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                 SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)\\d+")))]

## Anti-joins
# Dropped by ADNI -> 70
dropped       <- malf_vols[!adnimerge, on = "PTID", PTID]
malf_vols     <- malf_vols[!PTID %in% dropped]    # 1676

firstpass     <- malf_vols[!acq_fails, on = "PTID"]  # Acquisition -> 1641
malf_vols_final <- firstpass[!malf_qc_fails, on = "PTID"] # Segm -> 1638

# Remove useless col and add method
malf_vols_final[, `:=`(ID = NULL, METHOD = "malf")]
rm(malf_vols, dropped, firstpass, malf_qc_fails)

### Non-local Patch-based Segmentation QC
fpath         <- here('lists/qrater_nlpb_2022-12-20.csv')
if (!file.exists(fpath)) str_glue("File: {fpath} ",
                                  "is required but could not be found.") |>
                          stop()
nlpb_qc1      <- fread(fpath, header = FALSE, col.names = c("ID", "QC"))

fpath         <- here('lists/qrater_nlpb_reg_fails_2022-12-20.csv')
if (!file.exists(fpath)) str_glue("File: {fpath} ",
                                  "is required but could not be found.") |>
                          stop()
nlpb_qc2      <- fread(fpath, header = FALSE, col.names = c("ID", "QC"))

nlpb_qc       <- rbindlist(list(nlpb_qc1, nlpb_qc2))

nlpb_qc_fails <- nlpb_qc[QC == "Fail"]
nlpb_qc_fails[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                     SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)S*\\d+")))]
rm(nlpb_qc, nlpb_qc1, nlpb_qc2)

## HC/VC volumes
fpath         <- here('data/derivatives/adni-bl_volumes_hcvc_nlpb.csv')
if (!file.exists(fpath)) str_glue("File: {fpath} ",
                                  "is required but could not be found.") |>
                          stop()
nlpb_vols     <- fread(fpath) ## 1746
nlpb_vols[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                 SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)\\d+")))]

## Anti-joins
dropped       <- nlpb_vols[!adnimerge, on = "PTID", PTID]
nlpb_vols     <- nlpb_vols[!PTID %in% dropped]  # 1676

firstpass     <- nlpb_vols[!acq_fails, on = "PTID"]       # Acquisition -> 1641
nlpb_vols_final <- firstpass[!nlpb_qc_fails, on = "PTID"] # Segm -> 1638

# Remove useless cols and add method
nlpb_vols_final[, `:=`(ID = NULL, METHOD = "nlpb")]
rm(nlpb_vols, dropped, firstpass, nlpb_qc_fails)

### CNN Segmentations QC
## Only QCed HV/VC segmentations
fpath         <- here('lists/qrater_cnn1_2022-12-20.csv')
if (!file.exists(fpath)) str_glue("File: {fpath} ",
                                  "is required but could not be found.") |>
                          stop()
cnn_qc1       <- fread(fpath, header = FALSE, col.names = c("ID", "QC"))

fpath         <- here('lists/qrater_cnn1_reg_fails_2022-12-20.csv')
if (!file.exists(fpath)) str_glue("File: {fpath} ",
                                  "is required but could not be found.") |>
                          stop()
cnn_qc2       <- fread(fpath, header = FALSE, col.names = c("ID", "QC"))

cnn_qc        <- rbindlist(list(cnn_qc1, cnn_qc2))

cnn_qc_fails  <- cnn_qc[QC == "Fail"]
cnn_qc_fails[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                    SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)S*\\d+")))]
rm(cnn_qc, cnn_qc1, cnn_qc2)

## HC/VC volumes
fpath         <- here('data/derivatives/adni-bl_volumes_hcvc_cnn.csv')
if (!file.exists(fpath)) str_glue("File: {fpath} ",
                                  "is required but could not be found.") |>
                          stop()
cnn_vols      <- fread(fpath) ## 1746

cnn_vols[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)\\d+")))]

## Anti-joins
dropped       <- cnn_vols[!adnimerge, on = "PTID", PTID]    # Dropped by ADNI -> 70
cnn_vols      <- cnn_vols[!PTID %in% dropped]  # 1676

firstpass     <- cnn_vols[!acq_fails, on = "PTID"]       # Acquisition -> 1641
cnn_vols_final  <- firstpass[!cnn_qc_fails, on = "PTID"] # Segm -> 1639

# Remove useless cols and add method
cnn_vols_final[, `:=`(ID = NULL, METHOD = "cnn")]
rm(cnn_vols, dropped, firstpass)

## Merge data.tables
volumes       <- rbindlist(list(malf_vols_final,
                                nlpb_vols_final,
                                cnn_vols_final),
                           use.names = TRUE)

rm(malf_vols_final, nlpb_vols_final, cnn_vols_final)


### HC/VC/AMY volumes
fpath         <- here('data/derivatives/adni-bl_volumes_hcvc-ag_cnn.csv')
if (!file.exists(fpath)) str_glue("File: {fpath} ",
                                  "is required but could not be found.") |>
                          stop()
cnn_vols2     <- fread(fpath) ## 1746
cnn_vols2[, `:=`(PTID = str_extract(ID, "\\d{3}_S_\\d{4}"),
                 SCANDATE = ymd(str_extract(ID, "(?<=\\d{4}_)\\d+")))]

## Anti-joins
dropped       <- cnn_vols2[!adnimerge, on = "PTID", PTID]    # Dropped by ADNI -> 70
cnn_vols2     <- cnn_vols2[!PTID %in% dropped]  # 1676

firstpass     <- cnn_vols2[!acq_fails, on = "PTID"]   # Acquisition -> 1641
volumes_hcvcag <- firstpass[!cnn_qc_fails, on = "PTID"] # Segm -> 1639
rm(cnn_vols2, dropped, firstpass, cnn_qc_fails, adnimerge)

# Remove useless cols
volumes_hcvcag[, ID := NULL]

## Save RDS
outdir        <- here("data/rds")
if (!file.exists(outdir)) dir.create(outdir, recursive = TRUE)
write_rds(acq_fails, here(outdir, "ptid_qc_discarded.rds"))
write_rds(volumes, here(outdir, "adni-bl_volumes_hcvc.rds"))
write_rds(volumes_hcvcag, here(outdir, "adni-bl_volumes_hcvc-ag.rds"))
rm(fpath, outdir)
