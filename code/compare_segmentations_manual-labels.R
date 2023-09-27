#!/usr/env/bin Rscript

library(here)
library(readr)
library(stringr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(glue)

### Read KAPPA DTs (rds) and data cleaning

## MALF
## Full and reduced labelsets
#malf_f <- read_rds(here("data/rds/manual-labels_dice_malf-full.rds"))
#malf_r <- read_rds(here("data/rds/manual-labels_dice_malf-reduced.rds"))
malf_r      <- here("data/derivatives/man-seg_kappa_hcvc_malf.csv") |> fread()

## Merge DTs
#malf <- rbindlist(list(malf_f, malf_r), use.names = TRUE)
malf        <- malf_r

# Pad id col
malf[, new_id := sprintf("%03d", sub)]

# Keep only useful columns
keep_cols   <- c("new_id", "kappa", "roi", "side")
malf        <- malf[, ..keep_cols]
setnames(malf, "new_id", "id")

# Add method column
malf[, method := "malf"]

# Clean memory
#rm(malf_f, malf_r, keep_cols)
rm(malf_r, keep_cols)

## Non-local patch based
## Full and reduced labelset
#nlpb_f <- read_rds(here("data/rds/manual-labels_dice_nlpb-full.rds"))
#nlpb_r <- read_rds(here("data/rds/manual-labels_dice_nlpb-reduced.rds"))
nlpb_r      <- here("data/derivatives/man-seg_kappa_hcvc_nlpb.csv") |> fread()

## Merge DTs
#nlpb <- rbindlist(list(nlpb_f, nlpb_r), use.names = TRUE)
nlpb        <- nlpb_r

# Add info
# Side
nlpb[label %in% 1:2, side := "left"]
nlpb[label %in% 3:4, side := "right"]
# ROI
nlpb[label %in% c(1, 3), roi := "hc"]
nlpb[label %in% c(2, 4), roi := "csf"]

# Keep only useful columns
keep_cols   <- c("id", "kappa", "roi", "side")
nlpb        <- nlpb[, ..keep_cols]

# Add method column
nlpb[, method := "nlpb"]

# Clean memory
rm(nlpb_r, keep_cols)

## nnUNet
#cnn_f <- read_rds(here("data/rds/manual-labels_dice_cnn-full.rds"))
#cnn_r <- read_rds(here("data/rds/manual-labels_dice_cnn-reduced.rds"))
cnn_r       <- here("data/derivatives/man-seg_kappa_hcvc_cnn.csv") |> fread()

## Clean
#cnn_f[, `:=`(roi = str_to_lower(str_extract(cls, "HC|VC|AG")),
             #portion = str_to_lower(str_extract(cls, "TAIL|BODY|HEAD")),
             #side = str_extract(cls, "(?<=_)L|R"))][, cls := NULL]

#cnn_r[, `:=`(roi = str_to_lower(str_extract(cls, "HC|CSF")),
             ##portion = "whole",
             #side = str_extract(cls, "L|R"))][, cls := NULL]

## Merge DTs
#cnn <- rbindlist(list(cnn_f, cnn_r), use.names = TRUE)
cnn         <- cnn_r

## Consistency
#cnn[roi == "vc", roi := "csf"]

##cnn[roi == "ag", portion := "whole"]

#cnn[side == "L", side := "left"]
#cnn[side == "R", side := "right"]

# Add info
# Side
cnn[label %in% 1:2, side := "left"]
cnn[label %in% 3:4, side := "right"]
# ROI
cnn[label %in% c(1, 3), roi := "hc"]
cnn[label %in% c(2, 4), roi := "csf"]

# Keep only useful columns
keep_cols   <- c("id", "kappa", "roi", "side")
cnn         <- cnn[, ..keep_cols]

# Add method column
cnn[, method := "cnn"]

# Clean memory
rm(cnn_r)

### Merge all three
scores <- rbindlist(list(malf, nlpb, cnn), use.names = TRUE, fill = TRUE)

#scores[is.na(fold), fold := 0]

scores[, `:=`(roi = str_to_upper(roi),
              side = str_to_title(side),
              #portion = str_to_title(portion),
              method = str_to_upper(method))]

scores[, roi := factor(roi, levels = c("HC", "CSF"))]

scores[, `:=`(mean = mean(kappa),
              median = median(kappa),
              sd = sd(kappa)),
       .(method, roi, side)]

setorder(scores, kappa)

### Boxplot Kappas HCVC
f_plot1  <- here("plots/manual-labels_segm-dice_hcvc")
fp1_png  <- paste0(f_plot1, ".png")
fp1_tiff <- paste0(f_plot1, ".tiff")
if(!file.exists(fp1_png) || !file.exists(fp1_tiff)) {
  ggplot(scores, aes(x = method, y = kappa)) +
    theme_classic(base_size = 12) +
    theme(
       text = element_text(size = 12),
       axis.text.y = element_text(size = 10),
       axis.text.x = element_text(size = 10, angle = 90,
                                  hjust = 0.95, vjust = 0.2)) +
    #geom_jitter(height = 0, width = .05, shape = 21, size = .4) +
    geom_violin(trim = FALSE, linewidth = 0.3) +
    stat_summary(fun.data = "median_hilow", geom = "pointrange",
                 colour = "darkred", size = 0.01, linewidth = 0.2) +
    geom_label_repel(data = scores[order(kappa), .SD[.N/2], .(method, roi, side)],
                     aes(y = kappa,
                         label = glue("{round(mean, 2)}\n({round(sd, 2)})")),
                     size = 3, box.padding = 1.5, alpha = 0.75) +
    facet_grid(factor(roi, levels = c("HC", "CSF")) ~ side,
               scales = "free") +
    labs(y = "Overlap Kappa", x = "Segmentation Method")

  if(!file.exists(fp1_png)){
    ggsave(fp1_png, width = 5, height = 5, units = "in", dpi = 600)
  }

  if(!file.exists(fp1_tiff)){
    ggsave(fp1_tiff, width = 5, height = 5, units = "in",
           device = "tiff", dpi = 600)
  }
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
vols[, ID := str_extract(ID, "\\d{3}")]
#vols[, `:=`(ID = str_extract(ID, "\\d{3}"),
            #HC = (HC) / 2, CSF = (CSF) / 2)]

vols    <- melt(vols[, .(ID, METHOD, LHC, LCSF, RHC, RCSF)],
                measure.vars = patterns("(HC|CSF)$"),
                value.name = "VOL",
                variable.name = "ROI")

vols[, SIDE := factor(str_extract(ROI, "L|R"),
                      levels = c("L", "R"),
                      labels = c("Left", "Right"))]
vols[ROI %like% "HC", ROI := "HC"]
vols[ROI %like% "CSF", ROI := "CSF"]

manual  <- vols[METHOD == "manual"]
segs    <- dcast(vols[METHOD != "manual"], ... ~ METHOD, value.var = "VOL")

vols    <- manual[segs, on = .(ID, SIDE, ROI),
                  .(ID, SIDE, ROI, MANUAL = VOL, malf, nlpb, cnn)]
vols    <- melt(vols,
                measure.vars = c("malf", "nlpb", "cnn"),
                variable.name = "SEGMENTATION",
                value.name = "VOLUME")
rm(manual, segs)

vols[, SEGMENTATION := factor(SEGMENTATION,
                              levels = c("cnn", "malf", "nlpb"),
                              labels = c("CNN", "MALF", "NLPB"))]

f_plot3  <- here("plots/manual-labels_segm-corr_hcvc")
fp3_png  <- paste0(f_plot3, ".png")
fp3_tiff <- paste0(f_plot3, ".tiff")
if(!file.exists(fp3_png) || !file.exists(fp3_tiff)) {
  # Palette
   cbPalette    <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  ggplot(vols, aes(x = VOLUME, y = MANUAL, colour = SEGMENTATION)) +
    theme_classic(base_size = 12) +
    theme(text = element_text(size = 12), legend.position = "bottom") +
    geom_point(size = 2, shape = 21) +
    geom_abline(intercept = 0, slope = 1,
                colour = cbPalette[1], linetype = "dashed") +
    geom_smooth(method = "lm", alpha = .2) +
    stat_cor(size = 2.7, label.x.npc = "right", label.y.npc = "bottom",
             hjust = "inward") +
    facet_grid(rows = vars(ROI), cols = vars(SIDE), scales = "free") +
    scale_colour_manual(values = cbPalette[-1]) +
    labs(x = "Computed volume", y = "Manual volume",
         colour = "Segmentation method")

  if(!file.exists(fp3_png)){
    ggsave(fp3_png, width = 8, height = 5, units = "in", dpi = 600)
  }

  if(!file.exists(fp3_tiff)){
    ggsave(fp3_tiff, width = 8, height = 5, units = "in",
           device = "tiff", dpi = 600)
  }
}

# Bland-Altman plot
vols[, `:=`(AVG = (MANUAL + VOLUME) / 2,
            DIFF = MANUAL - VOLUME)]
vols[, `:=`(MEAN_DIFF = mean(DIFF),
            LOW_CI = mean(DIFF) - 1.96 * sd(DIFF),
            HI_CI = mean(DIFF) + 1.96 * sd(DIFF)),
     .(ROI, SEGMENTATION)]



f_plot4  <- here("plots/manual-labels_segm-bland-altman_hcvc")
fp4_png  <- paste0(f_plot4, ".png")
fp4_tiff <- paste0(f_plot4, ".tiff")
if(!file.exists(fp4_png) || !file.exists(fp4_tiff)) {
  # Palette
   cbPalette    <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  ggplot(vols, aes(x = AVG, y = DIFF)) +
    theme_classic(base_size = 12) +
    theme(text = element_text(size = 12), legend.position = "bottom") +
    geom_point(shape = 21, colour = cbPalette[1]) +
    geom_smooth(method = "lm", alpha = .1, colour = cbPalette[3]) +
    geom_hline(data = vols[, .SD[1], .(ROI, SEGMENTATION)],
               aes(yintercept = MEAN_DIFF),
               colour = cbPalette[2], linetype = "dashed", alpha = .7) +
    geom_text(data = vols[ROI == "HC", .SD[1], SEGMENTATION],
              aes(x = 1200, y = MEAN_DIFF, label = round(MEAN_DIFF, 2)),
              size = 2.5) +
    geom_text(data = vols[ROI == "CSF", .SD[1], SEGMENTATION],
              aes(x = 5500, y = MEAN_DIFF, label = round(MEAN_DIFF, 2)),
              size = 2.5) +
    geom_hline(data = vols[, .SD[1], .(ROI, SEGMENTATION)],
               aes(yintercept = LOW_CI),
               colour = cbPalette[2], linetype = "dashed", alpha = .7) +
    geom_text(data = vols[ROI == "HC", .SD[1], SEGMENTATION],
              aes(x = 1200, y = LOW_CI, label = round(LOW_CI, 2)),
              size = 2.5) +
    geom_text(data = vols[ROI == "CSF", .SD[1], SEGMENTATION],
              aes(x = 5500, y = LOW_CI, label = round(LOW_CI, 2)),
              size = 2.5) +
    geom_hline(data = vols[, .SD[1], .(ROI, SEGMENTATION)],
               aes(yintercept = HI_CI),
               colour = cbPalette[2], linetype = "dashed", alpha = .7) +
    geom_text(data = vols[ROI == "HC", .SD[1], SEGMENTATION],
              aes(x = 1200, y = HI_CI, label = round(HI_CI, 2)),
              size = 2.5) +
    geom_text(data = vols[ROI == "CSF", .SD[1], SEGMENTATION],
              aes(x = 5500, y = HI_CI, label = round(HI_CI, 2)),
              size = 2.5) +
    xlim(950, 6200) +
    facet_grid(rows = vars(ROI), cols = vars(SEGMENTATION), scales = "free") +
    labs(x = "Mean manual and computed volumes",
         y = "Manual - computed volumes")

  if(!file.exists(fp4_png)){
    ggsave(fp4_png, width = 10, height = 5, units = "in", dpi = 600)
  }

  if(!file.exists(fp4_tiff)){
    ggsave(fp4_tiff, width = 10, height = 5, units = "in",
           device = "tiff", dpi = 600)
  }
}
