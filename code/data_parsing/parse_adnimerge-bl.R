#!/usr/bin/env Rscript

library(here)
library(data.table)
library(glue)
library(lubridate)


## List of images
fpath             <- here("lists/adni_baseline.lst")
if (!file.exists(fpath)) glue("File: {fpath} ",
                              "is required but could not be found.") |> stop()

bl_list           <- fpath |>
                    fread(header = FALSE, col.names = c("PTID", "SCANDATE"))
bl_list[, SCANDATE := ymd(SCANDATE)]

## "Baseline" segmented images are not necessarily "Baseline" in ADNIMERGE
## Need to identify visit with SCANDATE of MRILIST
fpath             <- here("data/MRILIST.csv")
if (!file.exists(fpath)) glue("File: {fpath} ",
                              "is required but could not be found.") |> stop()

mrilist           <- fpath |>
                    fread(select = c("SUBJECT", "VISIT", "SCANDATE")) |>
                    unique()

mrilist[, SCANDATE := ymd(SCANDATE)]

bl_mri            <- mrilist[bl_list, on = .(SUBJECT = PTID, SCANDATE)]

## Pick the closest ADNIMERGE EXAMDATE for each PTID's SCANDATE
fpath               <- here("data/ADNIMERGE.csv")
if (!file.exists(fpath)) glue("File: {fpath} ",
                              "is required but could not be found.") |> stop()

adnimerge_dates   <- fpath |> fread(select = c("PTID", "VISCODE", "EXAMDATE"))

bl_mri            <- adnimerge_dates[bl_mri,
                                     on = .(PTID = SUBJECT)
                                     ][!is.na(EXAMDATE)]
bl_mri[, DATEDIFF := abs(ymd(EXAMDATE) - ymd(SCANDATE))]
bl_mri            <- bl_mri[bl_mri[order(DATEDIFF), .I[1], PTID]$V1]

## ADNIMERGE
adnimerge         <- fpath |>
                     fread(select = c("PTID", "RID",
                                      "DX_bl", "DX",
                                      "AGE", "PTGENDER", "PTEDUCAT",
                                      "PTETHCAT", "PTRACCAT",
                                      "ADAS13", "CDRSB", "MMSE",
                                      "RAVLT_immediate",
                                      "RAVLT_learning",
                                      "RAVLT_forgetting",
                                      "RAVLT_perc_forgetting",
                                      "EXAMDATE", "Month",
                                      "ABETA", "PIB", "AV45", "FBB",
                                      "APOE4", "TAU", "PTAU",
                                      "Ventricles", "Hippocampus", "ICV",
                                      "FSVERSION"))

adnimerge[, DX := factor(DX, levels = c("CN", "MCI", "Dementia"),
                         labels = c("CH", "MCI", "AD"))]

# Select only Baseline data
adnimerge         <- adnimerge[bl_mri[, .(PTID, EXAMDATE, SCANDATE)],
                               on = .(PTID, EXAMDATE)]

# Use DX_bl for missing DX
adnimerge[is.na(DX), DX := DX_bl]
adnimerge[DX_bl == "CN", DX := "CH"]
adnimerge[DX_bl %like% "MCI", DX := "MCI"]
adnimerge[, DX := factor(DX)]

# Save RDS
outdir            <- here("data/rds")
if (!file.exists(outdir)) dir.create(outdir, recursive = TRUE)
readr::write_rds(adnimerge, here(outdir, "adnimerge_baseline.rds"))

# Clean workspace
rm(fpath, adnimerge_dates, bl_list, bl_mri, mrilist, outdir)
