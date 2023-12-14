#!/usr/bin/env Rscript
library(here)
library(readr)
library(data.table)
library(ggplot2)
library(GGally)

## Parse volumes
# ADNI FSvols
adnimerge     <- here("data/rds/adnimerge_baseline.rds") |> read_rds()
adni_vols     <- adnimerge[DX != "",
                           .(PTID, RID, DX, SCANDATE, Hippocampus, FSVERSION)]
#rm(adnimerge)

## UCSF HCvols
fsv4            <- "Cross-Sectional FreeSurfer (FreeSurfer Version 4.3)"
f_ucsf_v4       <- here("data/UCSFFSX_11_02_15_20Nov2023.csv")
fs4_vols        <- fread(f_ucsf_v4, select = c("RID", "OVERALLQC",
                                               "ST29SV", # LeftHC
                                               "ST88SV"))# RightHC
fs4_vols        <- fs4_vols[, .(RID, OVERALLQC, LHC = ST29SV, RHC = ST88SV,
                                Hippocampus = ST29SV + ST88SV)
                            ][adni_vols[FSVERSION == fsv4],
                              on = .(RID, Hippocampus)]
fs4_vols[, UCSFFS := 4.3]


f_ucsf_v5       <- here("data/UCSFFSX51_11_08_19_20Nov2023.csv")
fsv5            <- "Cross-Sectional FreeSurfer (5.1)"
fs5_vols        <- fread(f_ucsf_v5, select = c("RID", "EXAMDATE", "OVERALLQC",
                                               "LHIPQC", "RHIPQC",
                                               "ST29SV", # LeftHC
                                               "ST88SV"))# RightHC
fs5_vols        <- fs5_vols[, .(RID, OVERALLQC, LHIPQC, RHIPQC,
                                LHC = ST29SV, RHC = ST88SV,
                                Hippocampus = ST29SV + ST88SV)
                            ][adni_vols[FSVERSION == fsv5],
                              on = .(RID, Hippocampus)]
fs5_vols[, UCSFFS := 5.1]
rm(f_ucsf_v4, fsv4, f_ucsf_v5, fsv5)

ucsf_vols     <- rbindlist(list(fs4_vols[!is.na(Hippocampus),
                                .(PTID, SCANDATE, DX, UCSFFS, Hippocampus)],
                                fs5_vols[!is.na(Hippocampus),
                                .(PTID, SCANDATE, DX, UCSFFS, Hippocampus)]))

# House FSvols
fs6_vols1     <- here("data/ADNI_FS_hc.csv") |>
                fread(select = c(1, 3, 5:7),
                      col.names = c("PTID", "DATE", "LHC", "RHC", "BRAIN"))
fs6_vols2     <- here("data/ADNI_FS_hc_vc.csv") |>
                fread(select = c(1, 3, 8:9),
                      col.names = c("PTID", "DATE", "LCSF", "RCSF"))

fs6_vols      <- fs6_vols1[fs6_vols2, on = .(PTID, DATE)]
fs6_vols[, FS_house := LHC + RHC]

# Merge
fs_vols       <- fs6_vols[ucsf_vols,
                          on = .(PTID, DATE = SCANDATE),
                          .(PTID, DX, SCANDATE = DATE,
                            LHC, RHC, HC = LHC + RHC,
                            LCSF, RCSF, CSF = LCSF + RCSF,
                            BRAIN, UCSFFS, FS_house, FS_ucsf = Hippocampus)]
#rm(adni_vols, fs6_vols)

fs_vols[, DX := factor(DX,
                       levels = c("CN", "MCI", "Dementia"),
                       label = c("CN", "MCI", "AD"))]

# Remove discarded from QC
discarded     <- here("data/rds/ptid_qc_discarded.rds") |> read_rds()
fs_vols       <- unique(fs_vols[!discarded, on = "PTID"])

# Versions
fs_vols[, .N, UCSFFS]
## 5.1: 701
## 4.3: 545

# Write fs_vols rds
write_rds(fs_vols, here("data/rds/adni-bl_volumes_freesurfer.rds"))

## Plot
# Palette
cbPalette     <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                   "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

g <- ggpairs(fs_vols[!is.na(FS_house) & !is.na(FS_ucsf),
                     .(PTID, DX, v6.0 = FS_house, v4.3_v5.1 = FS_ucsf)],
             columns = 3:4,
             aes(colour = DX, alpha = .7),
             upper = list(continuous = wrap("cor", method = "spearman"))) +
  theme_classic(base_size = 12) +
  theme(text = element_text(size = 14)) +
  scale_fill_manual(values = cbPalette[c(2:3, 8)]) +
  scale_colour_manual(values = cbPalette[c(2:3, 8)]) +
  labs(title = "Similarity between UCSF FS and In-house FS6")

png(here("plots/adni-bl_similarity_freesurfer.png"),
    width = 10, height = 5, units = "in", res = 600)
print(g)
dev.off()

tiff(here("plots/adni-bl_similarity_freesurfer.tiff"),
     width = 10, height = 5, units = "in", res = 600)
print(g)
dev.off()
