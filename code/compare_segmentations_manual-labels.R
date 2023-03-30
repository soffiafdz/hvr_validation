#!/usr/env/bin Rscript

library(here)
library(readr)
library(stringr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(glue)
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(gridExtra))

### Read KAPPA DTs (rds) and data cleaning

## MALF
# Full and reduced labelsets
malf_f <- read_rds(here("data/rds/manual-labels_dice_malf-full.rds"))
malf_r <- read_rds(here("data/rds/manual-labels_dice_malf-reduced.rds"))

# Merge DTs
malf <- rbindlist(list(malf_f, malf_r), use.names = TRUE)

# Pad id col
malf[, new_id := sprintf("%03d", sub)]

# Keep only useful columns
keep_cols <- c("new_id", "dice", "roi", "side", "portion")
malf <- malf[, ..keep_cols]
setnames(malf, "new_id", "id")

# Add method column
malf[, method := "malf"]

# Clean memory
rm(malf_f, malf_r, keep_cols)

## Non-local patch based
# Full and reduced labelset
nlpb_f <- read_rds(here("data/rds/manual-labels_dice_nlpb-full.rds"))
nlpb_r <- read_rds(here("data/rds/manual-labels_dice_nlpb-reduced.rds"))

# Merge DTs
nlpb <- rbindlist(list(nlpb_f, nlpb_r), use.names = TRUE)
rm(nlpb_f, nlpb_r)

# Pad id column
nlpb[, new_id := sprintf("%03d", id)]

# Remove label column and old id
nlpb[, `:=`(label = NULL, id = NULL)]

# Rename columns for consistency
# new_id -> id; kappa -> dice (volume_similarity's kappa == minccmp's dice)
setnames(nlpb, c("new_id", "kappa"), c("id", "dice"))

# Add method column
nlpb[, method := "nlpb"]


## nnUNet
cnn_f <- read_rds(here("data/rds/manual-labels_dice_cnn-full.rds"))
cnn_r <- read_rds(here("data/rds/manual-labels_dice_cnn-reduced.rds"))

# Clean
cnn_f[, `:=`(roi = str_to_lower(str_extract(cls, "HC|VC|AG")),
             portion = str_to_lower(str_extract(cls, "TAIL|BODY|HEAD")),
             side = str_extract(cls, "(?<=_)L|R"))][, cls := NULL]

cnn_r[, `:=`(roi = str_to_lower(str_extract(cls, "HC|CSF")),
             portion = "whole",
             side = str_extract(cls, "L|R"))][, cls := NULL]

# Merge DTs
cnn <- rbindlist(list(cnn_f, cnn_r), use.names = TRUE)

# Consistency
cnn[roi == "vc", roi := "csf"]

cnn[roi == "ag", portion := "whole"]

cnn[side == "L", side := "left"]
cnn[side == "R", side := "right"]

setnames(cnn, "kappa", "dice")

# Add method column
cnn[, method := "cnn"]

### Merge all three
scores <- rbindlist(list(malf, nlpb, cnn), use.names = TRUE, fill = TRUE)
rm(malf, nlpb, cnn)

scores[is.na(fold), fold := 0]

scores[, `:=`(roi = str_to_upper(roi),
              side = str_to_title(side),
              method = str_to_upper(method), portion = str_to_title(portion))]

scores[, roi2 := paste(roi, portion, sep = "_")]

scores[, `:=`(roi = factor(roi, levels = c("HC", "CSF", "AG")),
              roi2 = factor(roi2,
                            levels = c("HC_Head", "HC_Body", "HC_Tail",
                                       "CSF_Head", "CSF_Body", "CSF_Tail",
                                       "HC_Whole", "CSF_Whole", "AG_Whole"),
                            labels = c("HC_head", "HC_body", "HC_tail",
                                       "CSF_head", "CSF_body", "CSF_tail",
                                       "HC", "CSF", "AG")))]

scores[, `:=`(mean = mean(dice),
              median = median(dice),
              sd = sd(dice)),
       .(method, roi2, side)]

setorder(scores, dice)

### Boxplot: HC+VC
f_plot1 <- here("plots/manual-labels_segm-dice_hcvc.png")
#if(!file.exists(f_plot1)) {
if(TRUE) {
  ggplot(scores[portion == "Whole" & roi != "AG" & fold == 0],
         aes(x = method, y = dice)) +
      theme_linedraw(base_size = 24) +
      theme(
         text = element_text(size = 24),
         axis.text.y = element_text(size = 20),
         axis.text.x = element_text(size = 20, angle = 90,
                                    hjust = 0.95, vjust = 0.2)) +
      geom_jitter(height = 0, width = .05, shape = 21, size = .5) +
      geom_violin(trim = FALSE, alpha = .1) +
      stat_summary(fun.data = "median_hilow",
                   geom = "pointrange",
                   alpha = 0.7, colour = "red") +
      geom_label_repel(data = scores[portion == "Whole" & roi != "AG",
                                     .SD[.N/2], .(method, roi, side)],
                       aes(y = dice,
                           label = glue("{round(mean, 3)} ({round(sd, 3)})")),
                       alpha = 0.7, size = 3) +
      facet_grid(factor(roi, levels = c("HC", "CSF")) ~ side,
                 scales = "free") +
      ylab('Overlap Kappa') +
      xlab('Segmentation Method') +
      #ylim(0.4,1.0) +
      ggtitle("Automatic segmentation: HC and CSF")

  ggsave(f_plot1, width = 10, height = 10, units = "in", dpi = "retina")
}

### Boxplot:
f_plot2 <- here("plots/manual-labels_segm-dice_hcvc-ag.png")
if(!file.exists(f_plot2)) {
#if(TRUE) {
  ggplot(scores[!roi2 %in% c("CSF", "HC")], aes(x = method, y = dice)) +
      theme_linedraw(base_size = 24) +
      theme(
         text = element_text(size = 24),
         axis.text.y = element_text(size = 20),
         axis.text.x = element_text(size = 20, angle = 90,
                                    hjust = 0.95, vjust = 0.2)) +
      geom_jitter(height = 0, width = .05, shape = 21, size = .5) +
      geom_violin(trim = FALSE, alpha = .1) +
      stat_summary(fun.data = "median_hilow",
                   geom = "pointrange",
                   alpha = 0.7, colour = "red") +
      geom_label_repel(data = scores[!roi2 %in% c("CSF", "HC"),
                                     .SD[.N/2], .(method, roi2, side)],
                       aes(y = dice,
                           label = glue("{round(mean, 3)} ({round(sd, 3)})")),
                       alpha = 0.7, size = 4) +
      facet_grid(side ~ roi2, scales = "free") +
      ylab('Overlap Kappa') +
      xlab('Segmentation Method') +
      #ylim(0.4,1.0) +
      ggtitle("Automatic segmentation: AG and subportions of HC and CSF")

  ggsave(f_plot2, width = 30, height = 10, units = "in", dpi = "retina")
}


## Volumes Correlation
# Read Volume CSVs
manual  <- here("data/derivatives/man-seg_volumes_hcvc_manual.csv") |> fread()
malf    <- here("data/derivatives/man-seg_volumes_hcvc_malf.csv") |> fread()
nlpb    <- here("data/derivatives/man-seg_volumes_hcvc_nlpb.csv") |> fread()
cnn     <- here("data/derivatives/man-seg_volumes_hcvc_cnn.csv") |> fread()

vols    <- rbindlist(list(manual[, METHOD := "manual"],
                          malf[, METHOD := "malf"],
                          nlpb[, METHOD := "nlpb"],
                          cnn[, METHOD := "cnn"]))
rm(manual, malf, nlpb, cnn)

# Clean
vols[, `:=`(ID = str_extract(ID, "\\d{3}"),
            HC = (LHC + RHC) / 2, CSF = (LCSF + RCSF) / 2)]

vols    <- melt(vols[, .(ID, METHOD, HC, CSF)],
                measure.vars = c("HC", "CSF"),
                value.name = "VOL",
                variable.name = "ROI")

manual  <- vols[METHOD == "manual"]
segs    <- dcast(vols[METHOD != "manual"], ... ~ METHOD, value.var = "VOL")

vols    <- manual[segs, on = .(ID, ROI),
                  .(ID, ROI, MANUAL = VOL, malf, nlpb, cnn)]
vols    <- melt(vols,
                measure.vars = c("malf", "nlpb", "cnn"),
                variable.name = "SEGMENTATION",
                value.name = "VOLUME")
rm(manual, segs)

vols[, SEGMENTATION := factor(SEGMENTATION,
                              levels = c("cnn", "malf", "nlpb"),
                              labels = c("CNN", "MALF", "NLPB"))]

f_plot3 <- here("plots/manual-labels_segm-corr_hcvc.png")
#if(!file.exists(f_plot3)) {
if(TRUE) {
  ggplot(vols, aes(x = MANUAL, y = VOLUME, colour = SEGMENTATION)) +
    theme_linedraw(base_size = 24) +
    theme(text = element_text(size = 24), legend.position = "bottom") +
    geom_point(size = 2, shape = 21) +
    geom_smooth(method = "lm") +
    stat_cor(size = 5) +
    facet_wrap(vars(ROI), scales = "free") +
    #scale_color_viridis_d(option = "H") +
    labs(x = "Mean computed volume", y = "Mean manual volume",
         colour = "Segmentation method")

    ggsave(f_plot3, width = 20, height = 10, units = "in", dpi = "retina")
}

# Bland-Altman plot
vols[, `:=`(AVG = (MANUAL + VOLUME) / 2,
            DIFF = MANUAL - VOLUME)]
vols[, `:=`(MEAN_DIFF = mean(DIFF),
            LOW_CI = mean(DIFF) - 1.96 * sd(DIFF),
            HI_CI = mean(DIFF) + 1.96 * sd(DIFF)),
     .(ROI, SEGMENTATION)]



f_plot4 <- here("plots/manual-labels_segm-bland-altman_hcvc.png")
#if(!file.exists(f_plot4)) {
if(TRUE) {
  ggplot(vols, aes(x = AVG, y = DIFF)) +
    theme_linedraw(base_size = 24) +
    theme(text = element_text(size = 24), legend.position = "bottom") +
    geom_point(shape = 21) +
    geom_smooth(method = "lm", se = FALSE) +
    geom_hline(data = vols[, .SD[1], .(ROI, SEGMENTATION)],
               aes(yintercept = MEAN_DIFF),
               colour = "red", linetype = "dashed", alpha = .7) +
    geom_text(data = vols[, .SD[1], .(ROI, SEGMENTATION)],
               aes(x = Inf, y = MEAN_DIFF, label = round(MEAN_DIFF, 3)),
               hjust = "inward", nudge_y = 70) +
    geom_hline(data = vols[, .SD[1], .(ROI, SEGMENTATION)],
               aes(yintercept = LOW_CI),
               colour = "red", linetype = "dashed", alpha = .7) +
    geom_text(data = vols[, .SD[1], .(ROI, SEGMENTATION)],
               aes(x = Inf, y = LOW_CI, label = round(LOW_CI, 3)),
               hjust = "inward", nudge_y = 70) +
    geom_hline(data = vols[, .SD[1], .(ROI, SEGMENTATION)],
               aes(yintercept = HI_CI),
               colour = "red", linetype = "dashed", alpha = .7) +
    geom_text(data = vols[, .SD[1], .(ROI, SEGMENTATION)],
               aes(x = Inf, y = HI_CI, label = round(HI_CI, 3)),
               hjust = "inward", nudge_y = 70) +
    facet_grid(rows = vars(ROI), cols = vars(SEGMENTATION), scales = "free") +
    labs(title = "Bland-Altman plots",
         x = "Mean manual and computed volumes",
         y = "Manual - computed volumes")

  ggsave(f_plot4, width = 20, height = 10, units = "in", dpi = "retina")
}
