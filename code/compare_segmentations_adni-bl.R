#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(ggplot2)
library(GGally)

## Read RDS objects
adnimerge     <- here("data/rds/adnimerge_baseline.rds") |>
                read_rds()
fs_volumes    <- here("data/rds/adni-bl_volumes_freesurfer.rds") |> read_rds()
seg_volumes   <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds") |>
                read_rds()


## Merge
adni          <- adnimerge[, .(PTID, PTGENDER)][fs_volumes, on = "PTID"]
volumes       <- seg_volumes[adni,
                            on = "PTID",
                            .(PTID, METHOD, DX, PTGENDER, FS_house, FS_adni,
                              HC = HC_l + HC_r)]
#rm(adnimerge, fs_volumes, seg_volumes, adni)

volumes       <- dcast(volumes[!is.na(METHOD)],
                       ... ~ METHOD,
                       value.var = "HC")
setnames(volumes,
         c("cnn", "malf", "nlpb"),
         c("CNN", "MALF", "NLPB"))

# By Sex
# Palette
cbPalette     <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

png(here("plots/adni-bl_similarity_sex.png"),
    width = 20, height = 10, units = "in", res = 300)

g <- ggpairs(volumes[DX != ""],
             columns = 4:8,
             ggplot2::aes(colour = PTGENDER, alpha = 0.7)) +
  theme_linedraw(base_size = 24) +
  theme(text = element_text(size = 24)) +
  scale_fill_manual(values = cbPalette[-1]) +
  scale_colour_manual(values = cbPalette[-1]) +
  labs(title = "Similarity between segmentations — Sex")
print(g)
dev.off()


# By Diagnosis
png(here("plots/adni-bl_similarity_dx.png"),
    width = 20, height = 10, units = "in", res = 300)
g <- ggpairs(volumes[DX != ""],
             columns = 4:8,
             ggplot2::aes(colour = DX, alpha = 0.7)) +
  theme_linedraw(base_size = 24) +
  theme(text = element_text(size = 24)) +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette) +
  labs(title = "Similarity between segmentations — Diagnosis")
print(g)
dev.off()
