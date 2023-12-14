#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(gtsummary)
library(dunn.test)

## Read RDS objects
adnimerge     <- here("data/rds/adnimerge_baseline.rds") |>
                read_rds()
fs_volumes    <- here("data/rds/adni-bl_volumes_freesurfer.rds") |> read_rds()
seg_volumes   <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds") |>
                read_rds()


# Merge
adni          <- adnimerge[fs_volumes, on = "PTID"]
DT            <- seg_volumes[adni, on = "PTID",
                            .(PTID, METHOD, DX, PTGENDER, AGE, PTEDUCAT, ADAS13,
                              RAVLT_learning = as.numeric(RAVLT_learning),
                              adni = FS_adni * SCALEFACTOR / 2000,
                              HC = HC_stx_mean,
                              HVR = HVR_mean)]

DT_hvr        <- DT[, .(PTID, METHOD, HVR)]
DT_hc         <- dcast(DT[, -"HVR"], ... ~ METHOD, value.var = "HC") |>
                melt(id.vars = 1:7, variable.name = "METHOD", value.name = "HC")

DT            <- DT_hvr[DT_hc, on = .(PTID, METHOD)]

demog.dt      <- unique(DT[, -c("PTID", "METHOD", "HC", "HVR")])
demog.dt[, DX := factor(DX,
                        levels = c("CN", "MCI", "Dementia"),
                        labels = c("CN", "MCI", "AD"))]

hc.dt         <- DT[, .(METHOD, HC, HVR)]

# N
demog.dt[, .N, DX]

# Females (N; percentage)
demog.dt[, .(N = sum(PTGENDER == "Female"),
             Perc = sum(PTGENDER == "Female") / .N * 100),
         DX]

# Education
demog.dt[, .(M = mean(PTEDUCAT), SD = sd(PTEDUCAT)), DX]
# Age
demog.dt[, .(M = mean(AGE), SD = sd(AGE)), DX]

# Cognition
demog.dt[!is.na(ADAS13),
         .(M = mean(ADAS13),
           SD = sd(ADAS13)), DX]

# Memory
demog.dt[!is.na(RAVLT_learning),
         .(M = mean(RAVLT_learning),
           SD = sd(RAVLT_learning)), DX]

# Demographics table
demog.dt |>
  tbl_summary(by = DX,
              label = list(PTGENDER ~ "Sex",
                           AGE ~ "Age (years)",
                           PTEDUCAT ~ "Education (years)",
                           RAVLT_learning ~ "RAVLT (learning subtest)"),
              statistic = all_continuous() ~ "{mean} ({sd})",
              missing_text = "Missing") |>
  modify_header(label ~ "**Variable**") |>
  modify_spanning_header(c("stat_1", "stat_2", "stat_3") ~ "**Clinical Label**") |>
  add_n() |>
  add_p() |>
  as_flex_table() |>
  flextable::save_as_docx(path = "data/derivatives/adni-bl_demog-table.docx")


## Post-hoc analyses
# Chi2: DX v Sex
demog.dt[, chisq.posthoc.test::chisq.posthoc.test(table(DX, PTGENDER))]

# Dunn.tests
demog.dt[, dunn.test(AGE, DX, method = "bonferroni")]
demog.dt[, dunn.test(PTEDUCAT, DX, method = "bonferroni")]
demog.dt[, dunn.test(ADAS13, DX, method = "bonferroni")]
demog.dt[, dunn.test(RAVLT_learning, DX, method = "bonferroni")]
