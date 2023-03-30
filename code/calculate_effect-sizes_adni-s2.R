#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(stringr)
library(bootES)
library(glue)
library(ggplot2)
library(ggridges)
library(ggtext)

## Read RDS objects
adnimerge     <- here("data/rds/adnimerge_s2.rds") |> read_rds()
seg1_volumes  <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds") |>
               read_rds()
seg2_volumes  <- here("data/rds/adni-s2_volumes_hc-stx-norm-nat_hvr.rds") |>
               read_rds()

## Merge
seg1_volumes  <- seg1_volumes[PTID %in% seg2_volumes[, PTID] & METHOD == "cnn"]
seg1_volumes[, METHOD := NULL]
seg_volumes   <- rbindlist(list(seg1_volumes, seg2_volumes))



hc_dt         <- adnimerge[fs_volumes, on = "PTID"]
hc_dt         <- hc_dt[seg_volumes, on = "PTID",
                       .(PTID, METHOD, DX, PTGENDER, ABETA, PIB, AV45,
                         HC           = HC_l + HC_r,
                         HC_stx       = HC_stx_l + HC_stx_r,
                         HC_norm      = HC_norm_l + HC_norm_r,
                         HVR          = HVR_l + HVR_r,
                         HC_fs        = FS_house,
                         HC_stx_fs    = FS_house * SCALEFACTOR,
                         HC_norm_fs   = FS_house / ICC)]
rm(adnimerge, volumes)

## ABeta positivity:
# SUVR       >  1.11 AV45 PET
# SUVR       >  1.2  Pitts compound-B PET
# CSF AB1-42 <= 980pg/ml
hc_dt[ABETA == ">1700", ABETA := "1700"][, ABETA := as.numeric(ABETA)]
hc_dt[ABETA != ""  | !is.na(PIB) | !is.na(AV45), ABETA_pos := "Negative"]
hc_dt[ABETA <= 980 | PIB > 1.11  | AV45 > 1.2  , ABETA_pos := "Positive"]

hc_dt         <- hc_dt[, `:=`(ABETA = NULL, PIB = NULL, AV45 = NULL)]

## Long -> Wide
hc_dt         <- dcast(hc_dt[!is.na(METHOD) & DX != ""],
                       ... ~ METHOD,
                       value.var = c("HC", "HC_stx", "HC_norm", "HVR"))

## Wide -> Long again
hc_dt         <- melt(hc_dt,
                      measure = patterns("^H"),
                      variable.name = "MEASURE",
                      value.name = "HC_value")

hc_dt[, `:=`(HC_measure = str_extract(MEASURE, "(HC_stx|HC_norm|HC|HVR)"),
             METHOD     = str_extract(MEASURE, "(fs|malf|nlpb|cnn)")
             )][,
        `:=`(MEASURE    = NULL,
             DX         = factor(DX,
                                 levels  = c("Dementia", "MCI", "CN"),
                                 labels  = c("AD", "MCI", "CN")),
             HC_measure = factor(HC_measure),
             METHOD     = factor(METHOD,
                                 levels  = c("fs", "malf", "nlpb", "cnn"),
                                 labels  = c("FSv6", "MALF", "NLPB", "CNN"))
             )]

## Effect sizes — sex
f_effvals_sex <- here("data/rds/effect_size_sex.rds")
if (file.exists(f_effvals_sex)) {
  effvals_sex <- read_rds(f_effvals_sex)
} else {
  dxs               <- hc_dt[, levels(DX)]
  hc_measures       <- hc_dt[, levels(HC_measure)]
  methods           <- hc_dt[, levels(METHOD)]
  effects <- bounds_l <- bounds_h <- vector()
  for (dx in dxs) {
    for (hc_measure in hc_measures) {
      for (method in methods) {
        if (method == "ADNI" && hc_measure == "HVR" ) {
          effects   <- c(effects, NA)
          bounds_l  <- c(bounds_l, NA)
          bounds_h  <- c(bounds_h, NA)
        } else {
          sprintf("Dx: %s; HC: %s; dt: %s", dx, hc_measure, dt)
          effect    <- bootES(hc_dt[DX          == dx         &
                                    METHOD      == method     &
                                    HC_measure  == hc_measure
                                    ][!is.na(HC_value)],
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
    }
  }

  effvals_sex <- data.table(METHOD        = rep(methods, times = 12),
                            HC_measure    = rep(hc_measures,
                                                times = 3, each = 4),
                            DX            = rep(dxs, each = 16),
                            EFFECT        = round(effects, 2),
                            BOUNDS_l      = round(bounds_l, 2),
                            BOUNDS_h      = round(bounds_h, 2))
  write_rds(effvals_sex[!is.na(EFFECT)], f_effvals_sex)
  rm(effects, bounds_l, bounds_h, f_effvals_sex)
}

## GGridges
f_plot_sex    <- here("plots/effect_sizes/effect_sizes_sex.png")
if (!file.exists(f_plot_sex)) {
  # Labels
  efflabs_sex <- effvals_sex[, .(METHOD,
                                 HC_measure,
                                 DX,
                                 LABEL = paste0(DX, ": ", EFFECT,
                                                " [", BOUNDS_l, ", ",
                                                BOUNDS_h, "]"))]
  efflabs_sex[DX == "MCI", LABEL := paste0("<br>", LABEL)]
  efflabs_sex[DX == "AD", LABEL := paste0("<br><br>", LABEL)]

  # Plot
  ggplot(hc_dt, aes(x = HC_value, y = DX, fill = PTGENDER)) +
  theme_linedraw(base_size = 24) +
  theme(text = element_text(size = 24), legend.position = "bottom") +
  geom_density_ridges(scale = 3, rel_min_height = 0.01, alpha = .2) +
  geom_richtext(data = efflabs_sex, aes(label = LABEL), size = 5,
                x = -Inf, y = Inf, hjust = -0.05, vjust = 2,
                fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  facet_grid(METHOD ~ HC_measure, scales = "free") +
  scale_fill_viridis_d() +
  labs(title = "Sex Effect sizes of CN vs MCI vs AD using HC volume (raw and normalized) and HVR",
       x = "Measure", y = "Diagnosis", fill = "Sex")

  ggsave(f_plot_sex, width = 30, height = 15, units = "in", dpi = "retina")
  #rm(efflabs_sex, f_plot_sex)
}

## Effect sizes — Amy+
f_effvals_amy <- here("data/rds/effect_size_amy.rds")
if (file.exists(f_effvals_amy)) {
  effvals_amy <- read_rds(f_effvals_amy)
} else {
  dxs               <- hc_dt[, levels(DX)]
  hc_measures       <- hc_dt[, levels(HC_measure)]
  methods           <- hc_dt[, levels(METHOD)]
  effects <- bounds_l <- bounds_h <- vector()
  for (dx in dxs) {
    for (hc_measure in hc_measures) {
      for (method in methods) {
        if (method == "ADNI" && hc_measure == "HVR" ) {
          effects   <- c(effects, NA)
          bounds_l  <- c(bounds_l, NA)
          bounds_h  <- c(bounds_h, NA)
        } else {
          sprintf("Dx: %s; HC: %s; dt: %s", dx, hc_measure, dt)
          effect    <- bootES(hc_dt[DX          == dx         &
                                    METHOD      == method     &
                                    HC_measure  == hc_measure
                                    ][!is.na(HC_value)],
                              R           = 5000,
                              data.col    = "HC_value",
                              group.col   = "ABETA_pos",
                              contrast    = c("Positive", "Negative"),
                              effect.type = "cohens.d")
          effect
          effects   <- c(effects, mean(effect$t))
          bounds_l  <- c(bounds_l, effect$bounds[1])
          bounds_h  <- c(bounds_h, effect$bounds[2])
        }
      }
    }
  }

  effvals_amy <- data.table(METHOD        = rep(methods, times = 12),
                            HC_measure    = rep(hc_measures,
                                                times = 3, each = 4),
                            DX            = rep(dxs, each = 16),
                            EFFECT        = round(effects, 2),
                            BOUNDS_l      = round(bounds_l, 2),
                            BOUNDS_h      = round(bounds_h, 2)
                            )
  rm(effects, bounds_l, bounds_h)
  write_rds(effvals_amy[!is.na(EFFECT)], f_effvals_amy)
}

## GGridges
f_plot_amy    <- here("plots/effect_sizes/effect_sizes_amy.png")
if (!file.exists(f_plot_amy)) {
  # Labels
  efflabs_amy <- effvals_amy[, .(METHOD,
                                 HC_measure,
                                 DX,
                                 LABEL = paste0(DX, ": ", EFFECT,
                                                " [", BOUNDS_l, ", ",
                                                BOUNDS_h, "]"))]
  efflabs_amy[DX == "MCI", LABEL := paste0("<br>", LABEL)]
  efflabs_amy[DX == "AD", LABEL := paste0("<br><br>", LABEL)]

  # Plot
  ggplot(hc_dt[!is.na(ABETA_pos)], aes(x = HC_value, y = DX, fill = ABETA_pos)) +
  theme_linedraw(base_size = 24) +
  theme(text = element_text(size = 24), legend.position = "bottom") +
  geom_density_ridges(scale = 3, rel_min_height = 0.01, alpha = .2) +
  geom_richtext(data = efflabs_amy, aes(label = LABEL), size = 5,
                x = -Inf, y = Inf, hjust = -0.05, vjust = 2,
                fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  facet_grid(METHOD ~ HC_measure, scales = "free") +
  scale_fill_viridis_d() +
  labs(title = "Amyloid Effect sizes of CN vs MCI vs AD using HC volume (raw and normalized) and HVR",
       x = "Measure", y = "Diagnosis", fill = "Amyloid+")

  ggsave(f_plot_amy, width = 30, height = 15, units = "in", dpi = "retina")
  #rm(efflabs_amy, f_plot_amy)
}

# Effect sizes by DX — ABeta+ only
# CN: 145; MCI: 391; AD: 210
f_effvals_dx  <- here("data/rds/effect_size_dx_abpos.rds")
if (file.exists(f_effvals_dx)) {
  effvals_dx  <- read_rds(f_effvals_dx)
} else {
  contrasts     <- list(c("AD", "CN"), c("MCI", "CN"), c("AD", "MCI"))
  hc_measures   <- hc_dt[, levels(HC_measure)]
  methods       <- hc_dt[, levels(METHOD)]
  effects <- bounds_l <- bounds_h <- vector()
  for (contrast in contrasts) {
    for (hc_measure in hc_measures) {
      for (method in methods) {
        if (method == "ADNI" && hc_measure == "HVR" ) {
          effects   <- c(effects, NA)
          bounds_l  <- c(bounds_l, NA)
          bounds_h  <- c(bounds_h, NA)
        } else {
          sprintf("Contrast: %s; HC: %s; method: %s",
                  contrast, hc_measure, method)
          effect    <- bootES(hc_dt[ABETA_pos  == "Positive"  &
                                    METHOD     == method      &
                                    HC_measure == hc_measure
                                    ][!is.na(HC_value)],
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
  }

  effvals_dx <- data.table(CONTRAST       = rep(c("AD-CN", "MCI-CN", "AD-MCI"),
                                                each = 16),
                           HC_measure     = rep(hc_measures,
                                               times = 3, each = 4),
                           METHOD         = rep(methods, times = 12),
                           EFFECT         = round(effects, 2),
                           BOUNDS_l       = round(bounds_l, 2),
                           BOUNDS_h       = round(bounds_h, 2))
  write_rds(effvals_dx[!is.na(EFFECT)], f_effvals_dx)
  rm(effects, bounds_l, bounds_h, f_effvals_dx)
}

## GGridges
f_plot_dx     <- here("plots/effect_sizes/effect_sizes_dx_abpos.png")
if (!file.exists(f_plot_dx)) {
  # Labels
  cn_label    <- paste0("<span style='color:",
                        viridisLite::viridis(3)[3],
                        "'>CN</span>")
  mci_label   <- paste0("<span style='color:",
                        viridisLite::viridis(3)[2],
                        "'>MCI</span>")
  ad_label    <- paste0("<span style='color:",
                        viridisLite::viridis(3)[1],
                        "'>AD</span>")

  rtext       <- c(glue("{ad_label} vs {cn_label}: "),
                   glue("<br>{mci_label} vs {cn_label}: "),
                   glue("<br><br>{ad_label} vs {mci_label}: "))
  rm(cn_label, mci_label, ad_label)

  efflabs_dx  <- effvals_dx[, .(METHOD, HC_measure,
                                suffix = paste0(": ", EFFECT,
                                                " [", BOUNDS_l, ", ",
                                                BOUNDS_h, "]"))]
  efflabs_dx[, rtext := rep(rtext, each = 15)]
  efflabs_dx[, LABEL := paste0(rtext, suffix)]
  efflabs_dx[, `:=`(rtext = NULL, suffix = NULL)]

  # Plot
  ggplot(hc_dt[ABETA_pos == "Positive"],
         aes(x = HC_value, y = DX, fill = DX)) +
  theme_linedraw(base_size = 24) +
  theme(text = element_text(size = 24), legend.position = "none") +
  geom_density_ridges(scale = 8, rel_min_height = 0.01, alpha = .2) +
  geom_richtext(data = efflabs_dx, aes(label = LABEL), size = 4,
                x = -Inf, y = Inf, hjust = -0.05, vjust = 2,
                fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, 4), "pt")) +
  facet_grid(METHOD ~ HC_measure, scales = "free") +
  scale_fill_viridis_d() +
  labs(title = "Effect sizes of Dx in Amy+ using HC volume (raw and normalized) and HVR",
       x = "Measure", y = "Diagnosis")

  ggsave(f_plot_dx, width = 30, height = 15, units = "in", dpi = "retina")
  rm(efflabs_dx, f_plot_dx)
}
