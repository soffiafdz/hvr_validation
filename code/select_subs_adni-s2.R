#!/usr/bin/env Rscript

library(here)
library(data.table)
library(lubridate)
library(stringr)

## List of preprocessed ADNI t1w
preprocessed    <- here("lists/adni_preproc.lst") |>
                  fread(col.names = c("PTID", "SCANDATE", "TESLA", "PATH"))
# Parse dates
preprocessed[, SCANDATE := ymd(SCANDATE)]

# Clean TESLA
preprocessed[TESLA == 2.9 | TESLA == 30000.0, TESLA := 3.0]
preprocessed[TESLA == 15000.0, TESLA := 1.5]

## Already segmented images
seg_list        <- here("lists/adni_baseline.lst") |>
                  fread(header = FALSE, col.names = c("PTID", "SCANDATE"))
# Parse dates
seg_list[, SCANDATE := ymd(SCANDATE)]

# Identify segmented
seg_list[, SEGMENTED := 1]
preprocessed    <- seg_list[preprocessed, on = .(PTID, SCANDATE)]
#rm(seg_list)

## MRILIST
# SCANDATES do not match ADNIMERGE's EXAMDATEs
# Must use MRILIST's SCANDATE
mrilist         <- here("data/MRILIST.csv") |>
                  fread(select = c("SUBJECT", "VISIT", "SCANDATE")) |>
                  unique()

# Parse dates
mrilist[, SCANDATE := ymd(SCANDATE)]

# Merge
preproc <- copy(preprocessed)
preprocessed    <- mrilist[preprocessed,
                           on = .(SUBJECT = PTID, SCANDATE),
                           .(PTID = SUBJECT,
                             TESLA, SCANDATE, VISIT, SEGMENTED, PATH)]
rm(mrilist)

## ADNIMERGE
adnimerge       <- here("data/ADNIMERGE.csv") |>
                  fread(select = c("PTID", "DX", "VISCODE", "EXAMDATE"))

# SCANDATE == EXAMDATE
match           <- adnimerge[preprocessed,
                             on = .(PTID, EXAMDATE = SCANDATE)
                             ][!is.na(VISCODE),
                             .(PTID, VISCODE, SCANDATE = EXAMDATE, TESLA)]

# Find closest date to session by PTID
unmatch         <- preprocessed[!match, on = .(PTID, SCANDATE)]
merged_dates    <- adnimerge[unmatch,
                             on = .(PTID),
                             .(PTID, VISCODE, SCANDATE, TESLA,
                               DIFF = abs(ymd(EXAMDATE) - ymd(SCANDATE))),
                             allow.cartesian = TRUE]
matched         <- merged_dates[merged_dates[order(DIFF),
                                             .I[1],
                                             .(PTID, SCANDATE)]$V1]

matched[is.na(VISCODE), .N] #70 PTID not in ADNIMERGE

# Final merging
viscodes        <- rbindlist(list(match,
                                  matched[!is.na(VISCODE),
                                          .(PTID, VISCODE,
                                            SCANDATE, TESLA)]))
preprocessed    <- preprocessed[viscodes, on = .(PTID, SCANDATE)]
preprocessed    <- adnimerge[preprocessed, on = .(PTID, VISCODE)]

#rm(match, unmatch, matched, merged_dates, viscodes)

##
# Get SCANDATEs for Segmentations
segm_dates      <- preprocessed[SEGMENTED == 1,
                                  .(PTID, SEGM_DATE = SCANDATE,
                                    DX_bl = DX, TESLA_bl = TESLA)]

# Remove subjects that don't have a first segmentation
segm_volumes    <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds") |>
                  read_rds()
unseg_dates     <- segm_dates[PTID %in% segm_volumes[METHOD == "cnn", PTID]
                              ][preprocessed[is.na(SEGMENTED)],
                              on = "PTID",
                              .(PTID, DX_bl, DX, TESLA, TESLA_bl, SCANDATE,
                                DIFF = abs(ymd(SEGM_DATE) - ymd(SCANDATE))),
                              allow.cartesian = TRUE]

# Between 300 and 500 days ~ 1 year
unseg_dates     <- unseg_dates[DIFF > 300 & DIFF < 500]

# Pick closer SCANDATE to baseline by subject
# and same Teslage as baseline
unseg_dates     <- unseg_dates[TESLA == TESLA_bl]
unseg_dates     <- unseg_dates[unseg_dates[order(DIFF), .I[1], PTID]$V1]
unseg_dates     <- unseg_dates[!is.na(PTID)]

# Extract 100 random sessions:
# Stable CN, MCI, AD and MCI converters

set.seed(1618)
cn_stable       <- unseg_dates[DX_bl == DX][DX == "CN"][sample(.N, 100)]
mci_stable      <- unseg_dates[DX_bl == DX][DX == "MCI"][sample(.N, 100)]
#mci_conv        <- unseg_dates[DX_bl != DX][DX == "Dementia"][sample(.N, 80)]
ad_stable       <- unseg_dates[DX_bl == DX][DX == "Dementia"][sample(.N, 100)]

to_be_segmented <- rbindlist(list(cn_stable, mci_stable, ad_stable))
#rm(segm_dates, unseg_dates, cn_stable, mci_stable, mci_conv, ad_stable)


# Export list
longit_list     <- to_be_segmented[,
                                   .(PTID, SCANDATE, TO_segment = 2)
                                   ][preprocessed, on = .(PTID, SCANDATE)]
longit_list[!is.na(TO_segment), SEGMENTED := 2]
longit_list     <- longit_list[!is.na(SEGMENTED),
                               .(PTID, SCANDATE, EXAMDATE,
                                 SEGMENTATION = SEGMENTED, PATH)]
setorder(longit_list, SCANDATE)
fwrite(longit_list, here("lists/adni_hcvc-segm.lst"), col.names = FALSE)
