#!/usr/bin/env Rscript

library(here)
library(data.table)

adnimerge     <- fread(here("data/ADNIMERGE_bl.csv"))
scanner       <- adnimerge[, .(PTID, TESLA = substr(FLDSTRENG, 1, 1))]

scanner[TESLA == "1", TESLA := "1.5"]

adnilst       <- fread(here("lists/adni_baseline.lst"),
                   col.names = c("PTID", "SESS"))

adni_scanner  <- scanner[adnilst, on = "PTID"]
setcolorder(adni_scanner, c("PTID", "SESS", "TESLA"))

fwrite(adni_scanner,
       file = here("lists/adni-bl_scanner.lst"),
       col.names = FALSE)
