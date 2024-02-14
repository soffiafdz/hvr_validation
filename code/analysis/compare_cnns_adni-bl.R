#!/usr/bin/env Rscript

library(here)
library(data.table)
library(ggplot2)
library(glue)
library(ggrepel)

### Read comparison file of CNN simple vs CNN simplified and plot
fpath         <- here('data/derivatives/adni-bl_kappa_hcvc_cnn.csv')
if (!file.exists(fpath)) glue("File: {fpath} ",
                              "is required but could not be found.") |> stop()

overlap <- fread(fpath)
setnames(overlap,
         c("xcorr", "lhc", "lcsf", "rhc", "rcsf"),
         c("XCorrelation", "KAPPA_HC-left", "KAPPA_VC-left",
           "KAPPA_HC_right", "KAPPA_VC_right"))

overlap_long  <- melt(overlap, id.vars = "id",
                      variable.name = "OVERLAP", value.name = "VALUE")

overlap_long[, `:=`(MEDIAN = median(VALUE), SD = sd(VALUE)), by = OVERLAP]
rm(overlap)

### Violinplots
ggplot(overlap_long, aes(x = OVERLAP, y = VALUE)) +
  theme_linedraw(base_size = 24) +
  theme(text = element_text(size = 24),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.2)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  #geom_jitter(height = 0, width = .05) +
  stat_summary(fun.data = "median_hilow", geom = "pointrange",
               alpha = 0.7, colour = "red") +
  geom_label_repel(data = unique(overlap_long[VALUE == MEDIAN,
                                 .(OVERLAP, VALUE, MEDIAN, SD)]),
                   aes(y = VALUE,
                       label = glue("{round(MEDIAN, 3)} ({round(SD, 3)})")),
                   size = 5, nudge_x = .15) +
  ylab(NULL) +
  xlab('Overlap Measure') +
  coord_cartesian(ylim = c(0.875, 1)) +
  ggtitle("ADNI baseline — Comparison CNN (Whole) and CNN (Summed-portions)")

ggsave(here("plots/adni-bl_comparison_cnns.png"),
    width = 25, height = 15, units = "in", dpi = "retina")
