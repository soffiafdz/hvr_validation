#!/usr/bin/env Rscript
library(here)
library(readr)
library(data.table)
library(ggplot2)
library(GGally)

## Parse volumes
# ADNI FSvols
adnimerge     <- here("data/rds/adnimerge_baseline.rds") |> read_rds()
adni_vols     <- adnimerge[, .(PTID, DX, SCANDATE, Hippocampus)]
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
                            FS_house, FS_adni = Hippocampus)]
rm(adni_vols, fs6_vols)

fs_vols[, DX := factor(DX,
                       levels = c("CN", "MCI", "Dementia"),
                       label = c("CN", "MCI", "AD"))]

# Write fs_vols rds
write_rds(fs_vols, here("data/rds/adni-bl_volumes_freesurfer.rds"))

## Plot
# Palette
cbPalette     <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

png(here("plots/adni-bl_similarity_freesurfer.png"),
    width = 20, height = 10, units = "in", res = 300)

g <- ggpairs(fs_vols[DX != "" & !is.na(FS_house),
                     .(PTID, DX, FS_house, FS_adni)],
             columns = 3:4,
             aes(colour = DX, alpha = .5)) +
  theme_linedraw(base_size = 24) +
  theme(text = element_text(size = 24)) +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette) +
  labs(title = "Similarity between ADNI FS and In-house FS6")
print(g)
dev.off()
