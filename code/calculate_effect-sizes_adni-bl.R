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

f_volumes   <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds")
if (file.exists(f_volumes)) {
  seg_volumes <- read_rds(f_volumes)
} else {
  here("code/adjust_hc-hvr_adni-bl.R") |> source()
}
rm(f_volumes)

## Merge
volumes     <- adnimerge[seg_volumes[METHOD == "cnn"], on = "PTID",
                         .(PTID, DX, PTGENDER, ABETA, PIB, AV45,
                           HC_raw   = HC_mean,
                           HC_stx   = HC_stx_mean,
                           HC_prop  = HC_prop_mean,
                           HC_pcp   = HC_pcp_mean,
                           HC_res   = HC_res_mean,
                           HVR      = HVR_mean,
                           HVR_pcp  = HVR_pcp_mean,
                           HVR_res  = HVR_res_mean)]


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

amy_posit[, ABETA_pos := NULL]
#amy_posit[, `:=`(ABETA_pos  = NULL,
                 #SIDE       = str_extract(HC_measure, "Left|Right"),
                 #HC_measure = str_extract(HC_measure, "(?<=_).*"))]

amy_posit[, `:=`(DX         = factor(DX,
                                     levels  = c("CN", "MCI", "Dementia"),
                                     labels  = c("CN", "MCI", "AD")),
                 PTGENDER   = factor(PTGENDER),
                 HC_measure = factor(HC_measure))]

## Effect sizes — sex
f_effvals_sex   <- here("data/rds/adni-bl_effect-sizes_sex_factorial_ext.rds")
if (file.exists(f_effvals_sex)) {
#if (FALSE) {
  effvals_sex   <- read_rds(f_effvals_sex)
} else {
  #sides         <- amy_posit[, levels(SIDE)]
  dxs           <- amy_posit[, levels(DX)]
  hc_measures   <- amy_posit[, levels(HC_measure)]
  effects <- bounds_l <- bounds_h <- vector()
  for (hc_measure in hc_measures) {
    effect    <- bootES(amy_posit[HC_measure == hc_measure],
                        data.col    = "HC_value",
                        group.col   = "PTGENDER",
                        block.col   = "DX",
                        contrast    = c("Male", "Female"),
                        effect.type = "cohens.d")
    effect
    #effects   <- c(effects, mean(effect$t))
    effects   <- c(effects, effect$t0)
    bounds_l  <- c(bounds_l, effect$bounds[1])
    bounds_h  <- c(bounds_h, effect$bounds[2])
  }
  effvals_sex   <- data.table(HC_measure  = hc_measures,
                              EFFECT      = round(effects, 3),
                              BOUNDS_l    = round(bounds_l, 3),
                              BOUNDS_h    = round(bounds_h, 3))
  write_rds(effvals_sex, f_effvals_sex)
  rm(effects, bounds_l, bounds_h, f_effvals_sex)
}

## GGridges
f_plot_sex      <- here("plots/adni-bl_effect-sizes_sex")
fp1_png         <- paste0(f_plot_sex, ".png")
fp1_tiff        <- paste0(f_plot_sex, ".tiff")
if(!file.exists(fp1_png) || !file.exists(fp1_tiff)) {
#if (TRUE) {
  # Palette
   cbPalette    <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  # Labels
  efflabs_sex   <- effvals_sex[, .(HC_measure,
                                   LABEL = paste0("d = ", EFFECT,
                                                  " [", BOUNDS_l,
                                                  ", ", BOUNDS_h, "]"))]
  # Plot
  ggplot(amy_posit, aes(x = HC_value)) +
  theme_linedraw(base_size = 12) +
  theme(text = element_text(size = 12), legend.position = "bottom",
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_density(aes(fill = PTGENDER, colour = PTGENDER), alpha = .2) +
  #geom_density_ridges(aes(fill = PTGENDER, colour = PTGENDER),
                      #scale = 3,
                      #rel_min_height = 0.01, alpha = .2) +
  geom_vline(data = amy_posit[, mean(HC_value),
                              .(DX, PTGENDER, HC_measure)],
             aes(xintercept = V1, colour = PTGENDER), linetype = "dashed") +
  geom_richtext(data = efflabs_sex, aes(label = LABEL), size = 2.5,
                x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.1) +
                #fill = NA, label.color = NA,
                #label.padding = grid::unit(rep(0, 4), "pt")) +
  #facet_wrap(DX ~ factor(HC_measure,
                         #levels = c("HC_raw", "HC_prop", "HC_stx",
                                    #"HC_pcp", "HC_res", "HVR", "HVR_pcp",
                                    #"HVR_res")),
  facet_wrap(DX ~ HC_measure, scales = "free") +
  scale_fill_manual(values = cbPalette[-1]) +
  scale_colour_manual(values = cbPalette[-1], guide = "none") +
  labs(title = "Effect sizes of sex by Dx using HC volume (unnormalized and normalized)",
       x = "Measure", y = NULL, fill = "Sex")

  if(!file.exists(fp1_png)){
    ggsave(fp1_png, width = 13, height = 13, units = "in", dpi = 600)
  }

  if(!file.exists(fp1_tiff)){
    ggsave(fp1_tiff, width = 13, height = 13, units = "in",
           device = "tiff", dpi = 600)
  }
  #rm(efflabs_sex, f_plot_sex)
}

## Effect sizes — Dx by Sex
f_effvals_dx    <- here("data/rds/adni-bl_effect-sizes_dx_glass_ext.rds")
if (file.exists(f_effvals_dx)) {
#if (FALSE) {
  effvals_dx    <- read_rds(f_effvals_dx)
} else {
  #contrasts     <- list(c("AD", "CN"), c("MCI", "CN"), c("AD", "MCI"))
  sexes         <- amy_posit[, levels(PTGENDER)]
  hc_measures   <- amy_posit[, levels(HC_measure)]
  effects <- bounds_l <- bounds_h <- vector()
  for (hc_measure in hc_measures) {
    effect    <- bootES(amy_posit[HC_measure == hc_measure],
                        #R           = 5000,
                        data.col    = "HC_value",
                        group.col   = "DX",
                        contrast    = c(CN=3, MCI=-1, AD=-2),
                        effect.type = "cohens.d",
                        glass.control = "CN")
    effect
    #effects   <- c(effects, mean(effect$t))
    effects   <- c(effects, effect$t0)
    bounds_l  <- c(bounds_l, effect$bounds[1])
    bounds_h  <- c(bounds_h, effect$bounds[2])
  }
  effvals_dx    <- data.table(HC_measure  = hc_measures,
                              EFFECT      = round(effects, 3),
                              BOUNDS_l    = round(bounds_l, 3),
                              BOUNDS_h    = round(bounds_h, 3))
  write_rds(effvals_dx, f_effvals_dx)
  rm(effects, bounds_l, bounds_h, f_effvals_dx)
}

# GGridges
f_plot_dx       <- here("plots/adni-bl_effect-sizes_dx")
fp2_png         <- paste0(f_plot_dx, ".png")
fp2_tiff        <- paste0(f_plot_dx, ".tiff")

if(!file.exists(fp2_png) || !file.exists(fp2_tiff)) {
#if (TRUE) {
  # Palette
   cbPalette    <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  # Labels

  efflabs_dx    <- effvals_dx[, .(HC_measure,
                                  LABEL = paste0("&Delta; = ", EFFECT,
                                                 " [", BOUNDS_l,
                                                 ", ", BOUNDS_h, "]"))]

  ggplot(amy_posit, aes(x = HC_value)) +
  theme_linedraw(base_size = 12) +
  theme(text = element_text(size = 12), legend.position = "bottom",
        axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_density(aes(fill = DX, colour = DX), alpha = .2) +
  geom_vline(data = amy_posit[, mean(HC_value), .(DX, HC_measure)],
             aes(xintercept = V1, colour = DX), linetype = "dashed") +
  geom_richtext(data = efflabs_dx, aes(label = LABEL), size = 5,
                x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.1) +
  facet_wrap(~ HC_measure, scales = "free") +
  #facet_wrap(~ factor(HC_measure,
                         #levels = c("HC_vol", "HC_icv", "HC_stx", "HVR")),
             #nrow = 1, scales = "free") +
  scale_fill_manual(values = cbPalette) +
  scale_colour_manual(values = cbPalette, guide = "none") +
  labs(title = "Effect sizes of Dx using HC volume (unnormalized and normalized)",
       x = "Measure", y = NULL, fill = "Dx")

  if(!file.exists(fp2_png)) {
    ggsave(fp2_png, width = 13, height = 13, units = "in", dpi = 600)
  }

  if(!file.exists(fp2_tiff)) {
    ggsave(fp2_tiff, width = 13, height = 13, units = "in",
           device = "tiff", dpi = 600)
  }
}
