#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(stringr)
library(progress)
library(bootES)
library(ggplot2)
#library(ggridges)
library(ggtext)

## Read RDS objects
adnimerge     <- here("data/rds/adnimerge_baseline.rds") |> read_rds()
seg_volumes   <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds") |>
                read_rds()

## Merge
volumes       <- adnimerge[seg_volumes, on = "PTID",
                           .(PTID, DX, PTGENDER, METHOD,
                             HC_vol = (HC_l       + HC_r      ) / 2,
                             #HC_stx = (HC_stx_l   + HC_stx_r  ) / 2,
                             HC_icv = (HC_norm_l  + HC_norm_r ) / 2,
                             HVR    = (HVR_l      + HVR_r     ) / 2)]

volumes       <- volumes[DX != ""] |>
                melt(measure        = patterns("^H"),
                     variable.name  = "HC_measure",
                     value.name     = "HC_value")

volumes[DX == "Dementia", DX := "AD"]

## Effect sizes — Dx by Sex
f_effvals_dx  <- here("data/rds/adni-bl_effect-sizes_aaic.rds")
if (file.exists(f_effvals_dx)) {
effvals_dx    <- read_rds(f_effvals_dx)
} else {
  contrasts   <- list(c("AD", "CN"), c("MCI", "CN"), c("AD", "MCI"))
  methods     <- volumes[METHOD != "fs6", unique(METHOD)]
  hc_measures <- volumes[, levels(HC_measure)]
  eff_vals <- bounds_l <- bounds_h <- vector()
  pb <- progress_bar$new(format = "Effect Sizes | :what [:bar] :current/:total",
                         total = length(contrasts) *
                           length(methods) *
                           length(hc_measures),
                         clear = FALSE, width = 75)
  for (contrast in contrasts) {
    for (method in methods) {
      for (hc_measure in hc_measures) {
        pb$tick(tokens = list(what = paste(contrast, method, hc_measure,
                                           sep = ":")))
        eff       <- bootES(volumes[HC_measure == hc_measure &
                                    METHOD == method],
                            R           = 5000,
                            data.col    = "HC_value",
                            group.col   = "DX",
                            #block.col   = "PTGENDER",
                            contrast    = contrast,
                            effect.type = "cohens.d")
        eff
        eff_vals  <- c(eff_vals, mean(eff$t))
        bounds_l  <- c(bounds_l, eff$bounds[1])
        bounds_h  <- c(bounds_h, eff$bounds[2])
      }
    }
  }
  effvals_dx  <- data.table(CONTRAST    = rep(c("AD—CN", "MCI—CN", "AD—MCI"),
                                              each = 3 * 3),
                            METHOD      = rep(methods, times = 3, each = 3),
                            HC_measure  = rep(hc_measures, times = 3 * 3),
                            EFFECT      = round(eff_vals, 3),
                            BOUNDS_l    = round(bounds_l, 3),
                            BOUNDS_h    = round(bounds_h, 3))
  write_rds(effvals_dx, f_effvals_dx)
  #rm(effects, bounds_l, bounds_h, f_effvals_dx)
}

## GGridges
 f_plot_dx    <- here("plots/aaic23_adni-bl_effect-sizes_dx.png")
#if (!file.exists(f_plot_dx)) {
if (T) {
  # Palette
  cbPalette   <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  plot.dt     <- volumes[METHOD != "fs6"]
  plot.dt[, METHOD := factor(METHOD, levels = c("malf", "nlpb", "cnn"),
                             labels = c("MALF", "NLPB", "CNN"))]
  plot.dt[, HC_measure := factor(HC_measure, c("HC_vol", "HC_icv", "HVR"),
                                 c("HCvol", "HCvol/ICV", "HVR"))]
  plot.dt[, DX := factor(DX, c("CN", "MCI", "AD"))]

  # Labels
  cn_label    <- paste0("<span style='color:", cbPalette[2], "'>CN</span>")
  mci_label   <- paste0("<span style='color:", cbPalette[3], "'>MCI</span>")
  ad_label    <- paste0("<span style='color:", cbPalette[4], "'>AD</span>")
  rtext       <- data.table(CONTRAST  = c("AD—CN", "MCI—CN", "AD—MCI"),
                            PREFIX    = c(paste0(ad_label, " vs ",
                                                 cn_label, ": "),
                                          paste0(ad_label, " vs ",
                                                 mci_label, ": "),
                                          paste0(mci_label, " vs ",
                                                 cn_label, ": ")))
  efflabs_dx  <- rtext[effvals_dx, on = "CONTRAST",
                       .(CONTRAST, HC_measure, METHOD,
                         LABEL = paste0(PREFIX, "<br>", round(EFFECT, 2),
                                        " [", round(BOUNDS_l, 2), ",",
                                        round(BOUNDS_h, 2), "]"))]
  efflabs_dx  <- efflabs_dx[, .(LABEL = paste(LABEL[CONTRAST == "AD—CN"],
                                              LABEL[CONTRAST ==  "MCI—CN"],
                                              LABEL[CONTRAST ==  "AD—MCI"],
                                              sep = "<br>")),
                            .(HC_measure, METHOD)]

  efflabs_dx[, METHOD := factor(METHOD, levels = c("malf", "nlpb", "cnn"),
                             labels = c("MALF", "NLPB", "CNN"))]
  efflabs_dx[, HC_measure := factor(HC_measure, c("HC_vol", "HC_icv", "HVR"),
                                    c("HCvol", "HCvol/ICV", "HVR"))]

  # Plot
  ggplot(plot.dt, aes(x = HC_value)) +
    theme_classic(base_size = 24) +
    theme(text = element_text(size = 24), legend.position = "bottom",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) +
    geom_histogram(aes(colour = DX), fill = "transparent", bins = 45) +
    geom_vline(data = plot.dt[, mean(HC_value), .(DX, HC_measure, METHOD)],
               aes(xintercept = V1, colour = DX),
               linetype = "dashed", linewidth = 1.2, alpha = .8) +
    geom_richtext(data = efflabs_dx, aes(label = LABEL),
                  size = 7, x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3) +
    facet_grid(METHOD ~ HC_measure, scales = "free") +
    scale_colour_manual(values = cbPalette[-1], guide = "none") +
    labs(title = "Effect sizes of Dx using HC volume (raw and normalized) and HVR",
         x = "Measure", y = NULL, fill = "Dx")

  ggsave(f_plot_dx, width = 23, height = 17, units = "in", dpi = "retina")
  rm(efflabs_dx, f_plot_dx)
}
