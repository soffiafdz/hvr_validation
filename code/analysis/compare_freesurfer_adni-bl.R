#!/usr/bin/env Rscript
library(here)
library(readr)
library(data.table)
library(ggplot2)
library(GGally)

## INPUT
fpath       <- here("data/rds/adni-bl_volumes_freesurfer.rds")
if (file.exists(fpath)) {
  fs_vols   <- read_rds(fpath)
} else {
  here("code/data_parsing/parse_freesurfer-vols.R") |> source()
}

## Plot
# Palette
cbPalette   <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

g <- ggpairs(fs_vols[!is.na(FS_house) & !is.na(FS_ucsf),
                     .(PTID, DX, v6.0 = FS_house, v4.3_v5.1 = FS_ucsf)],
             columns = 3:4,
             aes(colour = DX, alpha = .7),
             upper = list(continuous = wrap("cor", method = "spearman"))) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 14)) +
  scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
  scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
  labs(title = "Similarity between UCSF FS and In-house FS6")

png(here("plots/adni-bl_similarity_freesurfer.png"),
    width = 10, height = 5, units = "in", res = 600)
print(g)
dev.off()

tiff(here("plots/adni-bl_similarity_freesurfer.tiff"),
     width = 10, height = 5, units = "in", res = 600)
print(g)
dev.off()
