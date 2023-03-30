#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(lme4)
library(glue)
library(ggplot2)
library(ggrepel)
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
volumes[, `:=`(SESS      = stringr::str_extract(SCAN, "ses-\\d{3}"),
               ACQ       = stringr::str_extract(SCAN, "(?<=_).*(?=_)"),
               RUN       = stringr::str_extract(SCAN, "run-\\d+"))]

# Z-scoring
volumes[, `:=`(HC_z      = (HC_mean - mean(HC_mean)) / sd(HC_mean),
               HC_stx_z  = (HC_stx_mean - mean(HC_stx_mean)) / sd(HC_stx_mean),
               HC_norm_z = (HC_norm_mean - mean(HC_norm_mean)) / sd(HC_norm_mean),
               HVR_z     = (HVR_mean - mean(HVR_mean)) / sd(HVR_mean))]


## Merge Covariate data
f_covars      <- here("data/SIMON_pheno.csv")
covars        <- fread(f_covars)
rm(f_covars)

# Clean
covars[, `:=`(SESS       = sprintf("ses-%03d", Session),
              date       = lubridate::mdy(Acquisition_date))]

# Merge
simon_dt      <- covars[volumes, on = "SESS",
                        .(SCANNER_man   = factor(manufacturer),
                          SCANNER_model = factor(man_model_name),
                          TIME_shift_bl = as.numeric((date - min(date))/365.25),
                          HC            = HC_z,
                          HC_stx        = HC_stx_z,
                          HC_norm       = HC_norm_z,
                          HVR           = HVR_z)]
#rm(volumes, covars)

# Remove one outlier

## Mixed Effects models
mod_hc        <- lmer(HC ~ TIME_shift_bl + (1 | SCANNER_man),
                      data = simon_dt)
mod_hc_stx    <- lmer(HC_stx ~ TIME_shift_bl + (1 | SCANNER_man),
                      data = simon_dt)
mod_hc_norm   <- lmer(HC_norm ~ TIME_shift_bl + (1 | SCANNER_man),
                      data = simon_dt)
mod_hvr       <- lmer(HVR ~ TIME_shift_bl + (1 | SCANNER_man),
                      data = simon_dt)

mod_hc2       <- lmer(TIME_shift_bl ~ HC + (1 | SCANNER_man),
                      data = simon_dt)
mod_hc_stx2   <- lmer(TIME_shift_bl ~ HC_stx + (1 | SCANNER_man),
                      data = simon_dt)
mod_hc_norm2  <- lmer(TIME_shift_bl ~ HC_norm + (1 | SCANNER_man),
                      data = simon_dt)
mod_hvr2      <- lmer(TIME_shift_bl ~ HVR + (1 | SCANNER_man),
                      data = simon_dt)

## Residuals comparisons
resids_dt     <- data.table(ID      = 1:93,
                            HC      = residuals(mod_hc),
                            HC_stx  = residuals(mod_hc_stx),
                            HC_norm = residuals(mod_hc_norm),
                            HVR     = residuals(mod_hvr))

resids_dt     <- melt(resids_dt, id.vars = "ID",
                    variable.name = "HC_measure", value.name = "RESIDUAL")

resids_dt[, `:=`(MEDIAN = median(RESIDUAL),
                 MEAN   = mean(RESIDUAL),
                 SD     = sd(RESIDUAL)),
          HC_measure]

## AIC | BIC
aic_bic_dt    <- data.table(c("HC", "HC_stx", "HC_norm", "HVR"),
                            AIC(mod_hc2, mod_hc_stx2, mod_hc_norm2, mod_hvr2),
                            BIC(mod_hc2, mod_hc_stx2, mod_hc_norm2, mod_hvr2))

fwrite(aic_bic_dt[, .(MEASURE = V1, df, AIC, BIC)],
       here("data/derivatives/simon_models_aic_bic.csv"))

## Plots
# Residuals plots
# HC
png(here("plots/simon_residuals_hc.png"),
    width = 10, height = 5, units = "in", res = 300)
p <- plot(mod_hc, main = "Residual plot: Hippocampus (Z-scored)")
print(p)
dev.off()

# HC-stx
png(here("plots/simon_residuals_hc_stx.png"),
    width = 10, height = 5, units = "in", res = 300)
p <- plot(mod_hc_stx, main = "Residual plot: Hippocampus-stx (Z-scored)")
print(p)
dev.off()

# HC-norm
png(here("plots/simon_residuals_hc_norm.png"),
    width = 10, height = 5, units = "in", res = 300)
p <- plot(mod_hc_norm, main = "Residual plot: Hippocampus-norm (Z-scored)")
print(p)
dev.off()

# HVR
png(here("plots/simon_residuals_hvr.png"),
    width = 10, height = 5, units = "in", res = 300)
p <- plot(mod_hvr, main = "Residual plot: HVR (Z-scored)")
print(p)
dev.off()

# Measures comparison of residuals
ggplot(resids_dt, aes(x = HC_measure, y = RESIDUAL)) +
  theme_linedraw(base_size = 24) +
  theme(text = element_text(size = 24)) +
  geom_jitter(height = 0, width = .05, shape = 21) +
  geom_violin(trim = FALSE, alpha = .1) +
  stat_summary(fun.data = "median_hilow", geom = "pointrange",
               alpha = .7, colour = "red") +
  geom_label_repel(data = resids_dt[, .SD[.N/2], HC_measure],
                   aes(y = MEDIAN,
                       label = glue("{signif(MEAN, 4)} ({round(SD, 3)})")),
                   alpha = .8, size = 5, nudge_x = 0.21) +
  ylab("Model residuals") +
  xlab("Hippocampal measures") +
  ggtitle("Residuals: <HC> ~ TIME_shift + (1 | SCANNER_manuf)")

ggsave(here("plots/simon_residuals_comparison.png"),
       width = 15, height = 10, units = "in", dpi = "retina")

# Point plot
ggplot(melt(simon_dt, measure.vars = patterns("^H")),
       aes(x = TIME_shift_bl, y = value, colour = SCANNER_man)) +
  theme_linedraw(base_size = 24) +
  theme(text = element_text(size = 24),
        legend.position = "bottom") +
  geom_point(size = 2, shape = 21) +
  geom_smooth(method = "lm", se = FALSE) +
  #stat_cor(size = 5) +
  facet_grid(cols = vars(variable)) +
  labs(title = "SIMON dataset. HC measures through time.",
       x = "Time (years)",
       y = "HC measure (Z-scored)",
       colour = "Scanner manufacturer")

ggsave(here("plots/simon_hc-measures_time.png"),
       width = 15, height = 7, units = "in", dpi = "retina")
