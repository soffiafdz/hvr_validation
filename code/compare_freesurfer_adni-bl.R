#!/usr/bin/env Rscript
library(here)
library(readr)
library(data.table)
library(ggplot2)
library(GGally)

## Parse volumes
# ADNI FSvols
adnimerge     <- here("data/rds/adnimerge_baseline.rds") |> read_rds()
adni_vols     <- adnimerge[, .(PTID, DX, SCANDATE, Hippocampus, FSVERSION)]
rm(adnimerge)

# House FSvols
fs6_vols      <- here("data/ADNI_FS_hc.csv") |>
                fread(select = c(1,3,5,6),
                      col.names = c("PTID", "DATE", "LHC", "RHC"))
fs6_vols[, FS_house := LHC + RHC]

# Merge
fs_vols       <- fs6_vols[adni_vols[!is.na(Hippocampus)],
                          on = .(PTID, DATE = SCANDATE),
                          .(PTID, DX, LHC, RHC,
                            FS_house, FSVERSION, FS_adni = Hippocampus)]
rm(adni_vols, fs6_vols)

fs_vols[, DX := factor(DX,
                       levels = c("CN", "MCI", "Dementia"),
                       label = c("CN", "MCI", "AD"))]

# Remove discarded from QC
discarded     <- here("data/rds/ptid_qc_discarded.rds") |> read_rds()
fs_vols       <- fs_vols[!discarded, on = "PTID"]

# Versions
fs_vols[, .N, FSVERSION]
## 5.1: 712
## 4.3: 545

# Write fs_vols rds
write_rds(fs_vols, here("data/rds/adni-bl_volumes_freesurfer.rds"))

## Plot
# Palette
cbPalette     <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

g <- ggpairs(fs_vols[!is.na(DX) & !is.na(FS_house),
                     .(PTID, DX, FS_house, FS_adni)],
             columns = 3:4,
             aes(colour = DX, alpha = .5)) +
  theme_linedraw(base_size = 12) +
  theme(text = element_text(size = 12)) +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette) +
  labs(title = "Similarity between ADNI FS and In-house FS6")

png(here("plots/adni-bl_similarity_freesurfer.png"),
    width = 10, height = 5, units = "in", res = 600)
print(g)
dev.off()

tiff(here("plots/adni-bl_similarity_freesurfer.tiff"),
     width = 10, height = 5, units = "in", res = 600)
print(g)
dev.off()
