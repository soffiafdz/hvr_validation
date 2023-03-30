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

# Write fs_vols rds
write_rds(fs_vols, here("data/rds/adni-bl_volumes_freesurfer.rds"))

## Plot
png(here("plots/adni-bl_similarity_freesurfer.png"),
    width = 20, height = 10, units = "in", res = 300)
g <- ggpairs(fs_vols[DX != "" & !is.na(FS_house),
                     .(PTID, DX, FS_house, FS_adni)],
             columns = 3:4,
             aes(colour = DX, alpha = .5)) +
  theme_linedraw(base_size = 24) +
  theme(text = element_text(size = 24)) +
  labs(title = "Similarity between ADNI FS and In-house FS6")
print(g)
dev.off()
