#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(stringr)
library(bootES)
library(glue)
library(ggplot2)
library(ggridges)
library(ggrepel)
library(ggtext)

## Read RDS objects
adnimerge     <- here("data/rds/adnimerge_s2.rds") |> read_rds()
segmentations <- here("lists/adni_hcvc-segm.lst") |>
                fread(select = 1:3,
                      col.names = c("PTID", "SCANDATE", "EXAMDATE"))
seg1_volumes  <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds") |>
               read_rds()
seg2_volumes  <- here("data/rds/adni-s2_volumes_hc-stx-norm-nat_hvr.rds") |>
               read_rds()

## Merge
seg1_volumes  <- seg1_volumes[PTID %in% seg2_volumes[, PTID] & METHOD == "cnn"]
seg1_volumes[, `:=`(METHOD = NULL, SESSION = "Baseline")]
seg2_volumes[, SESSION := "Follow_up"]

seg_volumes   <- rbindlist(list(seg1_volumes, seg2_volumes))

seg_volumes   <- segmentations[seg_volumes, on = .(PTID, SCANDATE)]

hc_dt         <- adnimerge[seg_volumes, on = .(PTID, SCANDATE),
                           .(PTID, SESSION, DX, PTGENDER,
                             HC_vol         = (HC_l + HC_r) / 2,
                             HC_icv         = (HC_norm_l + HC_norm_r) / 2,
                             HC_stx         = (HC_stx_l + HC_stx_r) / 2,
                             HVR            = (HVR_l + HVR_r) / 2)]
rm(adnimerge, segmentations, seg1_volumes, seg2_volumes, seg_volumes)

## Wide -> Long
hc_dt         <- melt(hc_dt,
                      measure = patterns("^H"),
                      variable.name = "HC_measure",
                      value.name = "HC_value")

hc_dt[, `:=`(SESSION    = factor(SESSION),
             DX         = factor(DX,
                                 levels  = c("CN", "MCI", "Dementia"),
                                 labels  = c("CN", "MCI", "AD")),
             HC_measure = factor(HC_measure,
                                 levels = c("HC_vol", "HC_icv", "HC_stx", "HVR")))]

hc_dt_wide    <- dcast(hc_dt, ... ~ SESSION, value.var = "HC_value")
hc_dt_wide[, DIFF := Baseline - Follow_up]

## Effect sizes — DX
f_effvals_s2  <- here("data/rds/adni-s2_effsizes.rds")
#if (file.exists(f_effvals_s2)) {
if (FALSE) {
  effvals_s2  <- read_rds(f_effvals_s2)
} else {
  #sides       <- hc_dt[, levels(SIDE)]
  dxs         <- hc_dt[, levels(DX)]
  hc_measures <- hc_dt[, levels(HC_measure)]
  #contrast    <- hc_dt[, levels(SESSION)]
  effects <- bounds_l <- bounds_h <- vector()
  for (dx in dxs) {
    for (hc_measure in hc_measures) {
      sprintf("Dx: %s; HC: %s", dx, hc_measure)
      effect   <- bootES(hc_dt_wide[DX == dx & HC_measure  == hc_measure],
                         data.col    = "DIFF",
                         effect.type = "cohens.d")
      effect
      #effects  <- c(effects, mean(effect$t))
      effects  <- c(effects, effect$t0)
      bounds_l <- c(bounds_l, effect$bounds[1])
      bounds_h <- c(bounds_h, effect$bounds[2])
    }
  }
  effvals_s2  <- data.table(
                            DX         = rep(dxs, each = 4),
                            HC_measure = rep(hc_measures, times = 3),
                            EFFECT     = round(effects, 3),
                            BOUNDS_l   = round(bounds_l, 3),
                            BOUNDS_h   = round(bounds_h, 3))
  write_rds(effvals_s2, f_effvals_s2)
  rm(effects, bounds_l, bounds_h, f_effvals_s2)
}

## GGridges
f_plot        <- here("plots/adni-s2_effsizes")
fp_png        <- paste0(f_plot, ".png")
fp_tiff       <- paste0(f_plot, ".tiff")
if(!file.exists(fp_png) || !file.exists(fp_tiff)) {
#if (TRUE) {
  # Palette
  cbPalette   <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  # Labels
  efflabs     <- effvals_s2[, .(HC_measure, DX,
                                LABEL = paste0("d = ", EFFECT,
                                               " [", BOUNDS_l,
                                               ", ", BOUNDS_h, "]"))]

  # Plot
  ggplot(hc_dt, aes(x = HC_value)) +
  theme_linedraw(base_size = 12) +
  theme(text = element_text(size = 12), legend.position = "bottom") +
  geom_density(aes(fill = SESSION, colour = SESSION), alpha = .2) +
  #geom_density_ridges(aes(fill = SESSION, colour = SESSION),
                      #rel_min_height = .01, alpha = .2) +
  geom_vline(data = hc_dt[, mean(HC_value), .(DX, SESSION, HC_measure)],
             aes(xintercept = V1, colour = SESSION), linetype = "dashed") +
  geom_richtext(data = efflabs, aes(label = LABEL), size = 2.5,
                x = -Inf, y = -Inf, hjust = -0.05, vjust = -0.1) +
  facet_wrap(factor(DX, levels = c("CN", "MCI", "AD")) ~
             factor(HC_measure, levels = c("HC_vol", "HC_icv", "HC_stx", "HVR")),
             scales = "free") +
  scale_fill_manual(values = cbPalette[-1]) +
  scale_colour_manual(values = cbPalette[-1], guide = "none") +
  labs(title = "1Y Follow-up effect sizes of CN vs MCI vs AD using HC volume (raw and normalized) and HVR",
       x = "Measure", y = NULL, fill = "Session")

  if(!file.exists(fp_png)){
    ggsave(fp_png, width = 13, height = 7, units = "in", dpi = 600)
  }

  if(!file.exists(fp_tiff)){
    ggsave(fp_tiff, width = 13, height = 7, units = "in",
           device = "tiff", dpi = 600)
  }
}
