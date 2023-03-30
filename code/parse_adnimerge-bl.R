#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lubridate)

## List of images
bl_list           <- here("lists/adni_baseline.lst") |>
                    fread(header = FALSE, col.names = c("PTID", "SCANDATE"))
bl_list[, SCANDATE := ymd(SCANDATE)]

## "Baseline" segmented images are not necessarily "Baseline" in ADNIMERGE
## Need to identify visit with SCANDATE of MRILIST

mrilist           <- here("data/MRILIST.csv") |>
                    fread(select = c("SUBJECT", "VISIT", "SCANDATE")) |>
                    unique()

mrilist[, SCANDATE := ymd(SCANDATE)]

bl_mri            <- mrilist[bl_list, on = .(SUBJECT = PTID, SCANDATE)]

## Pick the closest ADNIMERGE EXAMDATE for each PTID's SCANDATE
adnimerge_dates   <- here("data/ADNIMERGE.csv") |>
                    fread(select = c("PTID", "VISCODE", "EXAMDATE"))

bl_mri            <- adnimerge_dates[bl_mri,
                                     on = .(PTID = SUBJECT)
                                     ][!is.na(EXAMDATE)]
bl_mri[, DATEDIFF := abs(ymd(EXAMDATE) - ymd(SCANDATE))]
bl_mri            <- bl_mri[bl_mri[order(DATEDIFF), .I[1], PTID]$V1]


## ADNIMERGE
adnimerge         <- fread(here("data/ADNIMERGE.csv"),
                           select = c("PTID", "RID",
                                      "DX_bl", "DX",
                                      "AGE", "PTGENDER", "PTEDUCAT",
                                      "PTETHCAT", "PTRACCAT",
                                      "ADAS13", "CDRSB", "MMSE",
                                      "RAVLT_immediate",
                                      "RAVLT_learning",
                                      "RAVLT_forgetting",
                                      "RAVLT_perc_forgetting",
                                      "EXAMDATE", "Month",
                                      "ABETA", "PIB", "AV45",
                                      "APOE4", "TAU", "PTAU",
                                      "Ventricles", "Hippocampus", "ICV",
                                      "FSVERSION"))

## ApoE4
#apoe4             <- fread(here("data/APOERES.csv"),
                           #select = c("RID", "APGEN1", "APGEN2"))

#adnimerge_apoe4   <- apoe4[adnimerge, on = "RID"]
#adnimerge_apoe4[, RID := NULL]

# Select only PTIDs in bl_list
adnimerge         <- adnimerge[PTID %in% bl_mri[, PTID]]

# Latest recorded session
latest_visit      <- adnimerge[adnimerge[order(-Month), .I[1], PTID]$V1,
                               .(PTID, LATEST_VISIT = Month)]

adnimerge         <- latest_visit[adnimerge, on = "PTID"]

# Parse conversion
# Need to double check given not all SCANDATES are baseline
conversion        <- adnimerge[DX_bl != "AD" & DX == "Dementia",
                               .(PTID, CONV_MONTHS = Month)]

# Pick earliest date of conversion
conversion        <- conversion[conversion[order(CONV_MONTHS), .I[1], PTID]$V1]

# Merge All data
adnimerge         <- conversion[adnimerge, on = "PTID"
                                ][bl_mri[, .(PTID, EXAMDATE, SCANDATE)],
                                  on = .(PTID, EXAMDATE)]

# Shift latest visits and Conversion according to SCANDATE
adnimerge[, CONV_MONTHS := CONV_MONTHS - Month]
adnimerge[, LATEST_VISIT := LATEST_VISIT - Month]

## Conversion on 2y 3y 5y
# Converters
adnimerge[CONV_MONTHS <= Month, CONV_MONTHS := NA]
adnimerge[CONV_MONTHS <= 24, CONV_2Y := 1]
adnimerge[CONV_MONTHS <= 36, CONV_3Y := 1]
adnimerge[CONV_MONTHS <= 72, CONV_5Y := 1]

# Stables; Check that there is known DX for that time with LATEST_VISIT
adnimerge[is.na(CONV_2Y) & LATEST_VISIT >= 24, CONV_2Y := 0]
adnimerge[is.na(CONV_3Y) & LATEST_VISIT >= 36, CONV_3Y := 0]
adnimerge[is.na(CONV_5Y) & LATEST_VISIT >= 72, CONV_5Y := 0]

adnimerge[, `:=`(CONV_MONTHS = NULL, LATEST_VISIT = NULL)]

# Save RDS
readr::write_rds(adnimerge, here("data/rds/adnimerge_baseline.rds"))
