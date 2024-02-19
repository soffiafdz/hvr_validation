#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(gtsummary)
library(dunn.test)

## Read RDS objects
fpath         <- here("data/rds/adnimerge_baseline.rds")
if (file.exists(fpath)) {
  adnimerge   <- fpath |> read_rds()
} else {
  here("code/data_parsing/parse_adnimerge-bl.R") |> source()
}

fpath         <- here("data/rds/adni-bl_volumes_hcvc.rds")
if (file.exists(fpath)) {
  volumes     <- fpath |> read_rds()
} else {
  here("code/data_parsing/qc_segmentations_adni-bl.R") |> source()
}


# Merge
demog.dt      <- adnimerge[volumes[METHOD == "cnn"], on = "PTID",
                         .(PTID, DX, PTGENDER, AGE, PTEDUCAT, ADAS13,
                           RAVLT_learning = as.numeric(RAVLT_learning))]

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
demog.dt[, -1] |>
  tbl_summary(by = DX,
              label = list(PTGENDER ~ "Sex",
                           AGE ~ "Age (years)",
                           PTEDUCAT ~ "Education (years)",
                           RAVLT_learning ~ "RAVLT (learning)"),
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
