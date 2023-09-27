#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(bootES)
library(ggplot2)
library(GGally)
library(ggtext)
library(rlang)
library(gtsummary)
library(dunn.test)

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
                            .(PTID, METHOD, DX, PTGENDER,
                              FS_ADNI = FS_adni * SCALEFACTOR / 2000,
                              #HC = HC_stx_l + HC_stx_r,
                              HC = HC_stx_mean,
                              HVR = HVR_mean)]
#rm(adnimerge, fs_volumes, seg_volumes, adni)

## HC volume
hcv.dt        <- dcast(volumes[!is.na(METHOD), -"HVR"],
                       ... ~ METHOD, value.var = "HC")

setnames(hcv.dt,
         c("cnn", "malf", "nlpb", "fs6"),
         c("CNN", "MALF", "NLPB", "FS_V6"))

setcolorder(hcv.dt, c(1:4, 6))

hcv.dt.long   <- melt(hcv.dt, id.vars = c(1:3),
                      variable.name = "METHOD", value.name = "HCV")
hcv.dt.long   <- hcv.dt.long[!is.na(HCV)]

# Table
hcv.dt[, -c("PTID", "PTGENDER")] |>
  tbl_summary(by = DX,
              label = list(FS_ADNI ~ "FreeSurfer (v4.3 & v5.1)",
                           FS_V6 ~ "FreeSurfer (v6.0)"),
              statistic = all_continuous() ~ "{mean} ({sd})",
              missing_text = "Failures") |>
  modify_header(label ~ "**Automatic Method**") |>
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Clinical Label**") |>
  add_n() |>
  as_flex_table() |>
  flextable::save_as_docx(path = "data/derivatives/adni-bl_hcv-table.docx")

# Post-hoc differences between Methods by DX
#hcv_dunn.dt <- hcv.dt.long[, dunn.test(HCV, METHOD,
                                       #method = "bonferroni", list = TRUE), DX]

#hcv_dunn.dt[order(DX, comparisons)] |>
  #flextable() |>
  #separate_header() |>
  #colformat_double(digits = 3) |>
  #labelizor(part = "header", labels = c("chi2" = "Chi-squared",
                                        #"P" = "p-value",
                                        #"P.adjusted" = "Adjusted p-value",
                                        #"comparisons" = "Comparisons")) |>
  #bold(part = "header") |>
  ##fix_border_issues() |>
  ##bold(~ P < 0.05, j = 4) |>
  ##bold(~ P.adjusted < 0.05, j = 5) |>
  #autofit() |>
  #save_as_docx(path = "data/derivatives/adni-bl_hcv-posthoc.docx")





## HVR
hvr.dt        <- dcast(volumes[!is.na(METHOD), -"HC"],
                       ... ~ METHOD, value.var = "HVR")

hvr.dt[, FS_ADNI := NULL]
setnames(hvr.dt,
         c("cnn", "malf", "nlpb", "fs6"),
         c("CNN", "MALF", "NLPB", "FS_V6"))

setcolorder(hvr.dt, c(1:3, 5))

# Table
hvr.dt[, -c("PTID", "PTGENDER")] |>
  tbl_summary(by = DX,
              label = FS_V6 ~ "FreeSurfer (v6.0)",
              statistic = all_continuous() ~ "{mean} ({sd})",
              missing = "no") |>
  modify_header(label ~ "**Automatic Method**") |>
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Clinical Label**") |>
  add_n() |>
  as_flex_table() |>
  flextable::save_as_docx(path = "data/derivatives/adni-bl_hvr-table.docx")


## Effect sizes
# HC volume (sum of sides)
mtds  <- hcv.dt.long[, levels(METHOD)]
dxs   <- hcv.dt.long[, levels(DX)][-2] # Focus on CN-AD difference

# Sex
# Control by DX
effects <- bounds_l <- bounds_h <- vector()
for (mtd in mtds) {
  effect <- bootES(hcv.dt.long[METHOD == mtd],
                   data.col = "HCV",
                   group.col = "PTGENDER",
                   block.col = "DX",
                   contrast = c("Male", "Female"),
                   effect.type = "cohens.d")
  effects   <- c(effects, effect$t0)
  bounds_l  <- c(bounds_l, effect$bounds[1])
  bounds_h  <- c(bounds_h, effect$bounds[2])
}

effvals_hcv_sex <- data.table(METHOD    = mtds,
                              EFFECT    = round(effects, 2),
                              BOUNDS_l  = round(bounds_l, 2),
                              BOUNDS_h  = round(bounds_h, 2))

effvals_hcv_sex[, `:=`(PTGENDER = NA, DX = NA,
                       LABEL = paste0("d = ", EFFECT,
                                      " [", BOUNDS_l,
                                      ", ", BOUNDS_h, "]"))]

# CN vs AD
# Glass' delta (CN sd only)
effects <- bounds_l <- bounds_h <- vector()
for (mtd in mtds) {
  effect <- bootES(hcv.dt.long[METHOD == mtd & DX %in% dxs],
                   data.col       = "HCV",
                   group.col      = "DX",
                   contrast       = c("CN", "AD"),
                   effect.type    = "cohens.d",
                   glass.control  = "CN")
  effects   <- c(effects, effect$t0)
  bounds_l  <- c(bounds_l, effect$bounds[1])
  bounds_h  <- c(bounds_h, effect$bounds[2])
}

effvals_hcv_dx <- data.table(METHOD    = mtds,
                             EFFECT    = round(effects, 2),
                             BOUNDS_l  = round(bounds_l, 2),
                             BOUNDS_h  = round(bounds_h, 2))

effvals_hcv_dx[, `:=`(PTGENDER = NA, DX = NA,
                      LABEL = paste0("&Delta; = ", EFFECT,
                                     " [", BOUNDS_l,
                                     ", ", BOUNDS_h, "]"))]

# HVR (average of sides)
hvr.dt.long   <- melt(hvr.dt, id.vars = c(1:3),
                      variable.name = "METHOD", value.name = "HVR")
hvr.dt.long   <- hvr.dt.long[!is.na(HVR)]

mtds  <- hvr.dt.long[, levels(METHOD)]
dxs   <- hvr.dt.long[, levels(DX)][-2] # Focus on CN-AD difference

# Sex
# Control by DX
effects <- bounds_l <- bounds_h <- vector()
for (mtd in mtds) {
  effect <- bootES(hvr.dt.long[METHOD == mtd],
                   data.col = "HVR",
                   group.col = "PTGENDER",
                   block.col = "DX",
                   contrast = c("Male", "Female"),
                   effect.type = "cohens.d")
  effects   <- c(effects, effect$t0)
  bounds_l  <- c(bounds_l, effect$bounds[1])
  bounds_h  <- c(bounds_h, effect$bounds[2])
}

effvals_hvr_sex <- data.table(METHOD    = mtds,
                              EFFECT    = round(effects, 2),
                              BOUNDS_l  = round(bounds_l, 2),
                              BOUNDS_h  = round(bounds_h, 2))

effvals_hvr_sex[, `:=`(PTGENDER = NA, DX = NA,
                       LABEL = paste0("d = ", EFFECT,
                                      " [", BOUNDS_l,
                                      ", ", BOUNDS_h, "]"))]

# CN vs AD
# Glass' delta (CN sd only)
effects <- bounds_l <- bounds_h <- vector()
for (mtd in mtds) {
  effect <- bootES(hvr.dt.long[METHOD == mtd & DX %in% dxs],
                   data.col       = "HVR",
                   group.col      = "DX",
                   contrast       = c("CN", "AD"),
                   effect.type    = "cohens.d",
                   glass.control  = "CN")
  effects   <- c(effects, effect$t0)
  bounds_l  <- c(bounds_l, effect$bounds[1])
  bounds_h  <- c(bounds_h, effect$bounds[2])
}

effvals_hvr_dx <- data.table(METHOD    = mtds,
                             EFFECT    = round(effects, 2),
                             BOUNDS_l  = round(bounds_l, 2),
                             BOUNDS_h  = round(bounds_h, 2))

effvals_hvr_dx[, `:=`(PTGENDER = NA, DX = NA,
                      LABEL = paste0("&Delta; = ", EFFECT,
                                     " [", BOUNDS_l,
                                     ", ", BOUNDS_h, "]"))]

## Plots
# Palette
cbPalette     <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Diagonal plots: add mean vertical lines and effect sizes
diag_fun  <- function(data, mapping, var, labels.dt,...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(alpha = 1) +
    stat_summary(aes(xintercept = ..x.., y = 0), fun = mean,
                 geom = "vline", orientation = "y",
                 linetype = "dashed", alpha = 1) +
    geom_richtext(data = labels.dt[labels.dt$METHOD == as_label(mapping$x)],
                  aes(label = LABEL), inherit.aes = FALSE,
                  colour = "Black", fill = "White",
                  size = 3.5, x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.25)
}

# HC By Sex
g <- ggpairs(hcv.dt, columns = 4:8,
             aes(colour = PTGENDER, alpha = 0.7),
             diag = list(continuous = wrap(diag_fun,
                                           labels.dt = effvals_hcv_sex))) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 12)) +
  scale_fill_manual(values = cbPalette[-1]) +
  scale_colour_manual(values = cbPalette[-1]) +
  labs(title = "Similarity of HCvol between segmentations — Sex")

png(here("plots/adni-bl_similarity_hcv_sex.png"),
    width = 13, height = 7, units = "in", res = 600)
print(g)
dev.off()

tiff(here("plots/adni-bl_similarity_hcv_sex.tiff"),
    width = 13, height = 7, units = "in", res = 600)
print(g)
dev.off()

# HVR By Sex
g <- ggpairs(hvr.dt, columns = 4:7,
             aes(colour = PTGENDER, alpha = 0.7),
             diag = list(continuous = wrap(diag_fun,
                                           labels.dt = effvals_hvr_sex))) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 12)) +
  scale_fill_manual(values = cbPalette[-1]) +
  scale_colour_manual(values = cbPalette[-1]) +
  labs(title = "Similarity of HVR between segmentations — Sex")

png(here("plots/adni-bl_similarity_hvr_sex.png"),
    width = 13, height = 7, units = "in", res = 600)
print(g)
dev.off()

tiff(here("plots/adni-bl_similarity_hvr_sex.tiff"),
    width = 13, height = 7, units = "in", res = 600)
print(g)
dev.off()

# HC By DX
g <- ggpairs(hcv.dt, columns = 4:8,
             aes(colour = DX, alpha = 0.7),
             diag = list(continuous = wrap(diag_fun,
                                           labels.dt = effvals_hcv_dx))) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 12)) +
  scale_fill_manual(values = cbPalette[-1]) +
  scale_colour_manual(values = cbPalette[-1])
  #labs(title = "Similarity of HCvol between segmentations — DX")

png(here("plots/adni-bl_similarity_hcv_dx.png"),
    width = 13, height = 7, units = "in", res = 600)
print(g)
dev.off()

tiff(here("plots/adni-bl_similarity_hcv_dx.tiff"),
    width = 13, height = 7, units = "in", res = 600)
print(g)
dev.off()

# HVR By DX
g <- ggpairs(hvr.dt, columns = 4:7,
             aes(colour = DX, alpha = 0.7),
             diag = list(continuous = wrap(diag_fun,
                                           labels.dt = effvals_hvr_dx))) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 12)) +
  scale_fill_manual(values = cbPalette[-1]) +
  scale_colour_manual(values = cbPalette[-1])
  #labs(title = "Similarity of HVR between segmentations — DX")

png(here("plots/adni-bl_similarity_hvr_dx.png"),
    width = 13, height = 7, units = "in", res = 600)
print(g)
dev.off()

tiff(here("plots/adni-bl_similarity_hvr_dx.tiff"),
    width = 13, height = 7, units = "in", res = 600)
print(g)
dev.off()
