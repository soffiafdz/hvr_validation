#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(stringr)
library(lme4)
library(glue)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(ggpubr)
library(Hmisc)
library(gridExtra)

## Read volume data
f_volumes     <- here("data/rds/simon_volumes_hc-stx-norm-nat_hvr.rds")
if (file.exists(f_volumes)) {
  volumes     <- read_rds(f_volumes)
  rm(f_volumes)
} else {
  s_volumes   <- here("code/calculate_hc-stx-norm-nat_hvr_simon.R")
  source(s_volumes)
  rm(f_volumes, s_volumes)
}

# Parse SESS, ACQ, RUN
volumes[, `:=`(SESS      = str_extract(SCAN, "ses-\\d{3}"),
               ACQ       = str_extract(SCAN, "(?<=_).*(?=_)"),
               RUN       = str_extract(SCAN, "run-\\d+"))]

# Z-scoring
volumes[, `:=`(HC_z       = scale(HC_mean),
               HC_norm_z  = scale(HC_norm_mean),
               HC_stx_z   = scale(HC_stx_mean),
               HVR_z      = scale(HVR_mean))]

## Merge Covariate data
f_covars      <- here("data/SIMON_pheno.csv")
covars        <- fread(f_covars)
rm(f_covars)

# Clean
covars[, `:=`(SESS       = sprintf("ses-%03d", Session),
              date       = lubridate::mdy(Acquisition_date))]

# Merge
simon_z       <- covars[volumes, on = "SESS",
                        .(SESS          = factor(SESS),
                          SCANNER_man   = factor(manufacturer),
                          SCANNER_model = factor(man_model_name),
                          TIME_shift_bl = as.numeric((date - min(date))/365.25),
                          HC_z, HC_norm_z, HC_stx_z, HVR_z)]

simon_side    <- covars[volumes, on = "SESS",
                        .(SESS          = factor(SESS),
                          SCANNER_man   = factor(manufacturer),
                          SCANNER_model = factor(man_model_name),
                          TIME_shift_bl = as.numeric((date - min(date))/365.25),
                          HC_l, HC_r, HC_norm_l, HC_norm_r,
                          HC_stx_l, HC_stx_r, HVR_l, HVR_r)]

#simon_side    <- simon_side |> melt(measure.vars = patterns("(_l|_r)$"),
                                    #variable.name = "MEASURE",
                                    #value.name = "VALUE")
#simon_side[MEASURE %like% "_l", SIDE := "Left"]
#simon_side[MEASURE %like% "_r", SIDE := "Right"]
#simon_side    <- simon_side[, MEASURE := str_extract(MEASURE, ".*(?=_l|_r)")]

rm(volumes, covars)

# Rename
colnames     <- c("HC", "HC_stx", "HC_norm", "HVR")
setnames(simon_z, paste(colnames, "z", sep = "_"), colnames)

# Mixed Effects models
mod_hc        <- lmer(HC ~ TIME_shift_bl + (1|SCANNER_man) + (1|SESS),
                      data = simon_z)
mod_hc_stx    <- lmer(HC_stx ~ TIME_shift_bl + (1|SCANNER_man) + (1|SESS),
                      data = simon_z)
mod_hc_norm   <- lmer(HC_norm ~ TIME_shift_bl + (1|SCANNER_man) + (1|SESS),
                      data = simon_z)
mod_hvr       <- lmer(HVR ~ TIME_shift_bl + (1|SCANNER_man) + (1|SESS),
                      data = simon_z)

mod_hc2       <- lmer(TIME_shift_bl ~ HC_l + HC_r + (1 | SCANNER_man),
                      data = simon_side)
mod_hc_stx2   <- lmer(TIME_shift_bl ~ HC_stx_l + HC_stx_r + (1 | SCANNER_man),
                      data = simon_side)
mod_hc_norm2  <- lmer(TIME_shift_bl ~ HC_norm_l + HC_norm_r + (1 | SCANNER_man),
                      data = simon_side)
mod_hvr2      <- lmer(TIME_shift_bl ~ HVR_l + HVR_r + (1 | SCANNER_man),
                      data = simon_side)

## Residuals comparisons
resids_dt     <- data.table(ID      = 1:93,
                            HC_vol  = residuals(mod_hc),
                            HC_icv  = residuals(mod_hc_norm),
                            HC_stx  = residuals(mod_hc_stx),
                            HVR     = residuals(mod_hvr))

resids_dt     <- melt(resids_dt, id.vars = "ID",
                    variable.name = "HC_measure", value.name = "RESIDUAL")

resids_dt[, HC_measure := factor(HC_measure, levels = c("HC_vol", "HC_icv",
                                                        "HC_stx", "HVR"))]

resids_dt[, `:=`(
                 MEDIAN   = median(RESIDUAL),
                 MEAN     = mean(RESIDUAL),
                 SD       = sd(RESIDUAL),
                 MED_abs  = median(abs(RESIDUAL)),
                 MEAN_abs = mean(abs(RESIDUAL)),
                 SD_abs   = sd(abs(RESIDUAL))),
          HC_measure]

## AIC | BIC
aic_bic_dt    <- data.table(HC_metric = colnames,
                            AIC(mod_hc2, mod_hc_stx2, mod_hc_norm2, mod_hvr2),
                            BIC(mod_hc2, mod_hc_stx2, mod_hc_norm2, mod_hvr2))

fwrite(aic_bic_dt[, .(MEASURE = HC_metric, df, AIC, BIC)],
       here("data/derivatives/simon_models_aic_bic.csv"))

### Plots
## Palette
#cbPalette     <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   #"#F0E442", "#0072B2", "#D55E00", "#CC79A7")
## Residuals plots
## HC
#p1 <- plot(mod_hc,
          #main = "Residual plot: Hippocampus (Z-scored)",
          #xlim = c(-6, 3),
          #ylim = c(-3, 3))

## HC-stx
#p2 <- plot(mod_hc_stx,
          #main = "Residual plot: Hippocampus-stx (Z-scored)",
          #xlim = c(-6, 3),
          #ylim = c(-3, 3))

## HC-norm
#p3 <- plot(mod_hc_norm,
          #main = "Residual plot: Hippocampus-icv (Z-scored)",
          #xlim = c(-6, 3),
          #ylim = c(-3, 3))

## HVR
#p4 <- plot(mod_hvr,
          #main = "Residual plot: HVR (Z-scored)",
          #xlim = c(-6, 3),
          #ylim = c(-3, 3))

#resids_plot <- arrangeGrob(p1, p2, p3, p4, nrow = 2)

#ggsave(here("plots/simon_residuals.png"), plot = resids_plot,
       #width = 10, height = 10, units = "in", dpi = 600)

#ggsave(here("plots/simon_residuals.tiff"), plot = resids_plot,
       #width = 10, height = 10, units = "in", device = "tiff", dpi = 600)

## Measures comparison of residuals
#ggplot(resids_dt, aes(x = HC_measure, y = RESIDUAL)) +
  #theme_linedraw(base_size = 12) +
  #theme(text = element_text(size = 12)) +
  #geom_jitter(height = 0, width = .05, shape = 21) +
  #geom_violin(trim = FALSE, alpha = .1) +
  #stat_summary(fun.data = "median_hilow", geom = "pointrange",
               #alpha = .7, colour = "red") +
  #geom_label_repel(data = resids_dt[order(RESIDUAL), .SD[.N/2], HC_measure],
                   #aes(y = MEDIAN,
                       #label = glue("{signif(MEAN_abs, 3)}\n({round(SD_abs, 3)})")),
                #alpha = .8, size = 3, nudge_x = 0.25) +
  #ylab("Model residuals") +
  #xlab("Hippocampal measures") +
  #ggtitle("Residuals:\n<HC> ~ TIME_shift + (1|SCANNER_manuf) + (1|Session)")

#ggsave(here("plots/simon_residuals_comparison.png"),
       #width = 7, height = 5, units = "in", dpi = 600)

#ggsave(here("plots/simon_residuals_comparison.tiff"),
       #width = 7, height = 5, units = "in", device = "tiff", dpi = 600)

## Point plot
#setnames(simon_z, c("HC", "HC_norm"), c("HC_vol", "HC_icv"))
#simon_wide <- melt(simon_z, measure.vars = patterns("^H"))
#ggplot(simon_wide, aes(x = TIME_shift_bl, y = value, colour = SCANNER_man)) +
  #theme_linedraw(base_size = 12) +
  #theme(text = element_text(size = 12),
        #legend.position = "bottom") +
  #geom_point(size = 2, alpha = .3, shape = 21) +
  #geom_smooth(method = "lm", se = FALSE) +
  ##stat_cor(size = 5) +
  #facet_grid(cols = vars(variable)) +
  ##scale_shape_manual(values = c(21, 22)) +
  #scale_colour_manual(values = cbPalette) +
  #labs(title = "HC measures through time",
       #x = "Time (years)",
       #y = "HC measure (Z-scored)",
       #colour = "Scanner manufacturer")

#ggsave(here("plots/simon_hc-measures_time.png"),
       #width = 7, height = 5, units = "in", dpi = 600)

#ggsave(here("plots/simon_hc-measures_time.tiff"),
       #width = 7, height = 5, units = "in", device = "tiff", dpi = 600)
