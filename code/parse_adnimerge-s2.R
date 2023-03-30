#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lubridate)

## List of images
segm_list         <- here("lists/adni_hcvc-segm.lst") |>
                    fread(header = FALSE, select = c(1:4),
                          col.names = c("PTID", "SCANDATE",
                                        "EXAMDATE", "SESSION"))
segm_list[, `:=`(SCANDATE = ymd(SCANDATE),
                 EXAMDATE = ymd(EXAMDATE))]

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


## Select only PTIDs who have 2 sessions
subs_2sess        <- segm_list[PTID %in% segm_list[SESSION == 2, PTID],
                               .(PTID, EXAMDATE, SCANDATE)]

adnimerge         <- adnimerge[subs_2sess, on = .(PTID, EXAMDATE)]

## Save RDS
readr::write_rds(adnimerge, here("data/rds/adnimerge_s2.rds"))
