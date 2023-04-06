#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(stringr)
library(bootES)
library(ggplot2)
library(ggridges)
library(ggtext)

## Read RDS objects
adnimerge   <- here("data/rds/adnimerge_baseline.rds") |> read_rds()
seg_volumes <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds") |>
              read_rds()

## Merge
volumes     <- adnimerge[seg_volumes[METHOD == "cnn"], on = "PTID",
                         .(PTID, DX, PTGENDER, ABETA, PIB, AV45,
                           HC       = HC_mean,
                           HC_stx   = HC_stx_mean,
                           HC_norm  = HC_norm_mean,
                           HVR      = HVR_mean)]

## ABeta positivity:
## 666 subjects
## 131 CN; 348 CMI; 185 AD; 2 NA
## 367 Male / 299 Female

# SUVR       >  1.11 AV45 PET
# SUVR       >  1.2  Pitts compound-B PET
# CSF AB1-42 <= 980pg/ml
volumes[ABETA == ">1700", ABETA := "1700"][, ABETA := as.numeric(ABETA)]
volumes[ABETA != ""  | !is.na(PIB) | !is.na(AV45), ABETA_pos := "Negative"]
volumes[ABETA <= 980 | PIB > 1.11  | AV45 > 1.2  , ABETA_pos := "Positive"]

volumes     <- volumes[, `:=`(ABETA = NULL, PIB = NULL, AV45 = NULL)]

amy_posit   <- volumes[DX != "" & ABETA_pos == "Positive"] |>
              melt(measure        = patterns("^H"),
                   variable.name  = "HC_measure",
                   value.name     = "HC_value")

amy_posit[, `:=`(ABETA_pos  = NULL,
                 DX         = factor(DX,
                                     levels  = c("CN", "MCI", "Dementia"),
                                     labels  = c("CN", "MCI", "AD")),
                 PTGENDER   = factor(PTGENDER),
                 HC_measure = factor(HC_measure))]

## Effect sizes — sex
f_effvals_sex   <- here("data/rds/adni-bl_effect-sizes_sex.rds")
if (file.exists(f_effvals_sex)) {
  effvals_sex   <- read_rds(f_effvals_sex)
} else {
  dxs           <- amy_posit[, levels(DX)]
  hc_measures   <- amy_posit[, levels(HC_measure)]
  effects <- bounds_l <- bounds_h <- vector()
  for (dx in dxs) {
    for (hc_measure in hc_measures) {
      effect    <- bootES(amy_posit[DX == dx & HC_measure == hc_measure],
                          R           = 5000,
                          data.col    = "HC_value",
                          group.col   = "PTGENDER",
                          contrast    = c("Male", "Female"),
                          effect.type = "cohens.d")
      effect
      effects   <- c(effects, mean(effect$t))
      bounds_l  <- c(bounds_l, effect$bounds[1])
      bounds_h  <- c(bounds_h, effect$bounds[2])
    }
  }
  effvals_sex   <- data.table(HC_measure  = rep(hc_measures, times = 3),
                              DX          = rep(dxs, each = 4),
                              EFFECT      = round(effects, 3),
                              BOUNDS_l    = round(bounds_l, 3),
                              BOUNDS_h    = round(bounds_h, 3))
  write_rds(effvals_sex, f_effvals_sex)
  rm(effects, bounds_l, bounds_h, f_effvals_sex)
}

## GGridges
f_plot_sex      <- here("plots/adni-bl_effect-sizes_sex.png")
#if (!file.exists(f_plot_sex)) {
if (TRUE) {
  # Palette
   cbPalette    <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  # Labels
  efflabs_sex   <- effvals_sex[, .(HC_measure, DX,
                                   LABEL = paste0(EFFECT, "<br>[", BOUNDS_l,
                                                  ", ", BOUNDS_h, "]"))]
  # Plot
  ggplot(amy_posit, aes(x = HC_value)) +
  theme_linedraw(base_size = 24) +
  theme(text = element_text(size = 24), legend.position = "bottom",
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_density(aes(fill = PTGENDER, colour = PTGENDER), alpha = .2) +
  geom_vline(data = amy_posit[, mean(HC_value), .(DX, PTGENDER, HC_measure)],
             aes(xintercept = V1, colour = PTGENDER),
             linetype = "dashed") +
  geom_richtext(data = efflabs_sex, aes(label = LABEL), size = 6,
                x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2) +
                #fill = NA, label.color = NA,
                #label.padding = grid::unit(rep(0, 4), "pt")) +
  facet_wrap(factor(DX, levels = c("CN", "MCI", "AD")) ~ HC_measure,
             scales = "free") +
  scale_fill_manual(values = cbPalette[-1]) +
  scale_colour_manual(values = cbPalette[-1], guide = "none") +
  labs(title = "Effect sizes of sex by Dx using HC volume (raw and normalized) and HVR",
       x = "Measure", y = NULL, fill = "Sex")

  ggsave(f_plot_sex, width = 30, height = 13, units = "in", dpi = "retina")
  #rm(efflabs_sex, f_plot_sex)
}

## Effect sizes — Dx by Sex
f_effvals_dx    <- here("data/rds/adni-bl_effect-sizes_dx.rds")
if (file.exists(f_effvals_dx)) {
  effvals_dx    <- read_rds(f_effvals_dx)
} else {
  contrasts     <- list(c("AD", "CN"), c("MCI", "CN"), c("AD", "MCI"))
  sexes         <- amy_posit[, levels(PTGENDER)]
  hc_measures   <- amy_posit[, levels(HC_measure)]
  effects <- bounds_l <- bounds_h <- vector()
  for (contrast in contrasts) {
    for (hc_measure in hc_measures) {
      for (sex in sexes) {
        effect    <- bootES(amy_posit[!is.na(DX) &
                                      HC_measure == hc_measure &
                                      PTGENDER == sex],
                            R           = 5000,
                            data.col    = "HC_value",
                            group.col   = "DX",
                            contrast    = contrast,
                            effect.type = "cohens.d")
        effect
        effects   <- c(effects, mean(effect$t))
        bounds_l  <- c(bounds_l, effect$bounds[1])
        bounds_h  <- c(bounds_h, effect$bounds[2])
      }
    }
  }
  effvals_dx    <- data.table(CONTRAST    = rep(c("AD—CN", "MCI—CN", "AD—MCI"),
                                                each = 8),
                              HC_measure  = rep(hc_measures,
                                                times = 3, each = 2),
                              PTGENDER    = rep(sexes, times = 12),
                              EFFECT      = round(effects, 3),
                              BOUNDS_l    = round(bounds_l, 3),
                              BOUNDS_h    = round(bounds_h, 3))
  write_rds(effvals_dx, f_effvals_dx)
  rm(effects, bounds_l, bounds_h, f_effvals_dx)
}

## GGridges
f_plot_dx       <- here("plots/adni-bl_effect-sizes_dx.png")
#if (!file.exists(f_plot_dx)) {
if (TRUE) {
  # Palette
   cbPalette    <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  # Labels
  cn_label      <- paste0("<span style='color:", cbPalette[1], "'>CN</span>")
  mci_label     <- paste0("<span style='color:", cbPalette[2], "'>MCI</span>")
  ad_label      <- paste0("<span style='color:", cbPalette[3], "'>AD</span>")
  rtext         <- data.table(
                    CONTRAST  = c("AD—CN", "MCI—CN", "AD—MCI"),
                    PREFIX    = c(paste0(ad_label, " vs ", cn_label, ": "),
                                  paste0(ad_label, " vs ", mci_label, ": "),
                                  paste0(mci_label, " vs ", cn_label, ": ")))
  efflabs_dx    <- rtext[effvals_dx, on = "CONTRAST",
                         .(CONTRAST, HC_measure, PTGENDER,
                           LABEL = paste0(PREFIX, "<br>", EFFECT,
                                          " [", BOUNDS_l, ",", BOUNDS_h, "]"))]
  efflabs_dx    <- efflabs_dx[,
                              .(LABEL = paste(LABEL[CONTRAST == "AD—CN"],
                                              LABEL[CONTRAST ==  "MCI—CN"],
                                              LABEL[CONTRAST ==  "AD—MCI"],
                                              sep = "<br>")),
                              .(HC_measure, PTGENDER)]

  # Plot
  ggplot(amy_posit, aes(x = HC_value)) +
  theme_linedraw(base_size = 24) +
  theme(text = element_text(size = 24), legend.position = "bottom",
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_density(aes(fill = DX, colour = DX), alpha = .2) +
  geom_vline(data = amy_posit[, mean(HC_value), .(DX, PTGENDER, HC_measure)],
             aes(xintercept = V1, colour = DX), linetype = "dashed") +
  geom_richtext(data = efflabs_dx, aes(label = LABEL), size = 5,
                x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1) +
  facet_wrap(PTGENDER ~ HC_measure, nrow = 2,  scales = "free") +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette, guide = "none") +
  labs(title = "Effect sizes of Dx by sex using HC volume (raw and normalized) and HVR",
       x = "Measure", y = NULL, fill = "Dx")

  ggsave(f_plot_dx, width = 30, height = 10, units = "in", dpi = "retina")
  #rm(efflabs_dx, f_plot_dx)
}
