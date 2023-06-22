#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(stringr)
library(progress)
library(bootES)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(gridExtra)

## Read RDS objects
adnimerge     <- here("data/rds/adnimerge_baseline.rds") |> read_rds()
seg_volumes   <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds") |>
                read_rds()

## Merge
volumes       <- adnimerge[seg_volumes[METHOD == "cnn"],
                           on = "PTID",
                           .(PTID, DX, PTGENDER, AGE, ADAS13,
                             RAVLT_learning, RAVLT_immediate, RAVLT_forgetting,
                             HC_vol = (HC_l       + HC_r      ) / 2,
                             #HC_stx = (HC_stx_l   + HC_stx_r  ) / 2,
                             HC_icv = (HC_norm_l  + HC_norm_r ) / 2,
                             HVR    = (HVR_l      + HVR_r     ) / 2)]

volumes[DX == "Dementia", DX := "AD"]

# Scaling
volumes[, HC_vol_z := scale(HC_vol,
                            center  = volumes[DX == "CN", mean(HC_vol)],
                            scale   = volumes[DX == "CN", sd(HC_vol)])]
volumes[, HC_icv_z := scale(HC_icv,
                            center  = volumes[DX == "CN", mean(HC_icv)],
                            scale   = volumes[DX == "CN", sd(HC_icv)])]
volumes[, HVR_z := scale(HVR,
                         center  = volumes[DX == "CN", mean(HVR)],
                         scale   = volumes[DX == "CN", sd(HVR)])]

volumes[, c("HC_vol", "HC_icv", "HVR") := NULL]

volumes[, DX := factor(DX, c("CN", "MCI", "AD"))]

volumes_w     <- volumes[DX != ""] |>
                melt(measure        = patterns("^H"),
                     variable.name  = "HC_measure",
                     value.name     = "HC_value")

volumes_w[, HC_measure := factor(HC_measure,
                                 levels = c("HC_vol_z", "HC_icv_z", "HVR_z"),
                                 labels = c("HCvol", "HCvol/ICV", "HVR"))]

cbPalette     <- c("#999999A0", "#E69F00A0", "#56B4E9A0", "#009E73A0",
                   "#F0E442A0", "#0072B2A0", "#D55E00A0", "#CC79A7A0")

# AGE & HC by Sex
p1  <- ggplot(volumes_w[PTGENDER == "Male"],
              aes(AGE, HC_value, colour = DX)) +
  theme_classic(base_size = 18) +
  theme(text = element_text(size = 18), legend.position = "none") +
  scale_colour_manual(values = cbPalette[2:4]) +
  geom_point(shape = 21, alpha = .3) + geom_smooth(method = lm) +
  stat_cor(label.y.npc = "bottom") +
  facet_grid(cols = vars(HC_measure), scales = "free") +
  labs(title = "Correlation between Age and HC measures", subtitle = "Males",
       x = "Age", y = "HC Z-scores")

p2  <- ggplot(volumes_w[PTGENDER == "Female"],
              aes(AGE, HC_value, colour = DX)) +
  theme_classic(base_size = 18) +
  theme(text = element_text(size = 18), legend.position = "bottom") +
  scale_colour_manual(values = cbPalette[2:4]) +
  geom_point(shape = 21, alpha = .3) + geom_smooth(method = lm) +
  stat_cor(label.y.npc = "bottom") +
  facet_grid(cols = vars(HC_measure), scales = "free") +
  labs(subtitle = "Females", x = "Age", y = "HC Z-scores", colour = "Diagnosis")

p   <- grid.arrange(p1, p2, nrow = 2)

ggsave(here("plots/aaic_24_adni-bl_corr_age2.png"), p,
       width = 12, height = 10, units = "in", dpi = "retina", bg = "white")

# ADAS13 & HC
p1  <- ggplot(volumes_w[PTGENDER == "Male"],
              aes(ADAS13, HC_value, colour = DX)) +
  theme_classic(base_size = 18) +
  theme(text = element_text(size = 18), legend.position = "none") +
  scale_colour_manual(values = cbPalette[2:4]) +
  geom_point(shape = 21, alpha = .3) + geom_smooth(method = lm) +
  stat_cor(label.x.npc = .5) +
  facet_grid(cols = vars(HC_measure)) +
  labs(title = "Correlation between ADAS13 and HC measures", subtitle = "Males",
       x = "ADAS13 scores", y = "HC Z-scores")

p2  <- ggplot(volumes_w[PTGENDER == "Female"],
              aes(ADAS13, HC_value, colour = DX)) +
  theme_classic(base_size = 18) +
  theme(text = element_text(size = 18), legend.position = "bottom") +
  scale_colour_manual(values = cbPalette[2:4]) +
  geom_point(shape = 21, alpha = .3) + geom_smooth(method = lm) +
  stat_cor(label.x.npc = .5) +
  facet_grid(cols = vars(HC_measure)) +
  labs(subtitle = "Females", x = "ADAS13 scores", y = "HC Z-scores",
       colour = "Diagnosis")

p   <- grid.arrange(p1, p2, nrow = 2)

ggsave(here("plots/aaic_24_adni-bl_corr_cog2.png"), p,
       width = 12, height = 10, units = "in", dpi = "retina", bg = "white")

# RAVLT & HC
volumes_w2     <- volumes_w |>
                melt(measure        = patterns("^RAVLT_"),
                     variable.name  = "RAVLT_test",
                     value.name     = "RAVLT_score")

ggplot(volumes_w2[PTGENDER == "Male"],
       aes(RAVLT_score, HC_value, colour = DX)) +
  theme_classic(base_size = 18) +
  theme(text = element_text(size = 18), legend.position = "bottom") +
  scale_colour_manual(values = cbPalette[2:4]) +
  geom_point(shape = 21, alpha = .3) + geom_smooth(method = lm) + stat_cor() +
  facet_grid(rows = vars(HC_measure), cols = vars(RAVLT_test), scales = "free") +
  labs(title = "Correlation between RAVLT subtests and HC measures",
       subtitle = "Males",
       x = "RAVLT scores", y = "HC Z-scores", colour = "Diagnosis")

ggsave(here("plots/aaic_24_adni-bl_corr_mem1.png"),
       width = 12, height = 12, units = "in", dpi = "retina", bg = "white")

ggplot(volumes_w2[PTGENDER == "Female"],
       aes(RAVLT_score, HC_value, colour = DX)) +
  theme_classic(base_size = 18) +
  theme(text = element_text(size = 18), legend.position = "bottom") +
  scale_colour_manual(values = cbPalette[2:4]) +
  geom_point(shape = 21, alpha = .3) + geom_smooth(method = lm) + stat_cor() +
  facet_grid(rows = vars(HC_measure), cols = vars(RAVLT_test), scales = "free") +
  labs(title = "Correlation between RAVLT subtests and HC measures",
       subtitle = "Females",
       x = "RAVLT scores", y = "HC Z-scores", colour = "Diagnosis")

ggsave(here("plots/aaic_24_adni-bl_corr_mem2.png"),
       width = 12, height = 12, units = "in", dpi = "retina", bg = "white")

## Correlation comparisons
f_rcomps    <- here('data/rds/adni-bl_r-comps_aaic.rds')
if (file.exists(f_rcomps)) {
#if (F) {
  rcomps    <- read_rds(f_rcomps)
  rm(f_rcomps)
} else {
  covars    <- c("AGE", "ADAS13", "RAVLT_learning",
                 "RAVLT_immediate", "RAVLT_forgetting")
  dxs       <- volumes_w[, levels(DX)]
  measures  <- volumes_w[, levels(HC_measure)]
  comps     <- list(c(1,2), c(1,3), c(2,3))
  r_diffs <- bounds_l <- bounds_h <- vector()
  pb <- progress_bar$new(format = "R diffs | :what [:bar] :current/:total",
                         total = length(covars) * length(dxs) * length(comps),
                         clear = FALSE, width = 75)
  for (covar in covars) {
    for (dx in dxs) {
      for (comp in seq_along(comps)) {
        cols  <- c(covar, "HC_measure", "HC_value")
        pb$tick(tokens = list(what = paste(covar, dx, comp, sep = ":")))
        estim <- bootES(volumes_w[HC_measure %in% measures[comps[[comp]]] &
                                  DX == dx &
                                  !is.na(get(covar)),
                                  ..cols],
                        group.col = "HC_measure",
                        R = 3000)

        r_diffs   <- c(r_diffs, estim$t0)
        bounds_l  <- c(bounds_l, estim$bounds[1])
        bounds_h  <- c(bounds_h, estim$bounds[2])
      }
    }
  }
  rcomps  <- data.table(COVAR = rep(covars, each = 3 * 3),
                        DX    = rep(dxs, times = 5, each = 3),
                        COMP  = rep(c("HCvol-HCicv", "HCvol-HVR", "HCicv-HVR"),
                                    times = 3 * 3),
                        DIFF  = round(r_diffs, 3),
                        B_low = round(bounds_l, 3),
                        B_hi  = round(bounds_h, 3))
  write_rds(rcomps, f_rcomps)
  fwrite(rcomps, here('data/derivatives/adni-bl_rdiffs_aaic.csv'))
  rm(f_rcomps)
}

