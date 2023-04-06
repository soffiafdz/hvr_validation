#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(caret)
library(pROC)
library(purrr)
library(ggrepel)
library(glue)

# Read RDS objects
adnimerge <- here("data/rds/adnimerge_baseline.rds") |> read_rds()
volumes   <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds") |>
            read_rds()

# Merge
DT        <- volumes[METHOD == "cnn"
                   ][adnimerge, on = "PTID",
                     .(PTID, DX, AGE, PTGENDER, PTEDUCAT, CONV_3Y,
                       ABETA, PIB, AV45, APOE4, ADAS13, MMSE, CDRSB,
                       RAVLT_immediate, RAVLT_learning, RAVLT_forgetting,
                       HC      = HC_mean,
                       HC_stx  = HC_stx_mean,
                       HC_norm = HC_norm_mean,
                       HVR     = HVR_mean)
                   ][!is.na(HC)]

## ApoE4: binarized
DT[APOE4 == 0, APOE4_bin := 0]
DT[APOE4 != 0, APOE4_bin := 1]

## ABeta positivity:
## SUVR       >  1.11 AV45 PET
## SUVR       >  1.2  Pitts compound-B PET
## CSF AB1-42 <= 980pg/ml
DT[ABETA == ">1700", ABETA := "1700"][, ABETA := as.numeric(ABETA)]
DT[ABETA != ""  | !is.na(PIB) | !is.na(AV45), ABETA_pos := "Negative"]
DT[ABETA <= 980 | PIB > 1.11  | AV45 > 1.2  , ABETA_pos := "Positive"]

# Remove NAs
amy_posit <- DT[!is.na(CONV_3Y) &
                !is.na(ADAS13) &
                DX %in% c("CN", "MCI") &
                ABETA_pos == "Positive"]

amy_posit[, CONV_3Y := factor(CONV_3Y, levels = c(0, 1),
                              labels = c("Stable", "Progressor"))]

## 70% train/test split
set.seed(1618)
train_idx <- createDataPartition(amy_posit[!is.na(CONV_3Y), CONV_3Y],
                                 p = .8, list = FALSE)

dt_train  <- amy_posit[!is.na(CONV_3Y)][ train_idx]
dt_test   <- amy_posit[!is.na(CONV_3Y)][-train_idx]
rm(train_idx)

## PCA — Cognition
cog_vars  <- c("ADAS13", "CDRSB", "MMSE")
pca_cog   <- prcomp(dt_train[, ..cog_vars], scale = T)
dt_train[, COG_PC1 := pca_cog$x[, "PC1"]]
dt_test[, COG_PC1 := predict(pca_cog, dt_test[, ..cog_vars])[, "PC1"]]

## PCA — Memory
mem_vars  <- paste("RAVLT", c("immediate", "learning", "forgetting"),
                   sep = "_")
pca_mem   <- prcomp(dt_train[, ..mem_vars], scale = T)
dt_train[, MEM_PC1 := pca_mem$x[, "PC1"]]
dt_test[, MEM_PC1 := predict(pca_mem, dt_test[, ..mem_vars])[, "PC1"]]

## Covars
covars    <- c("AGE", "PTGENDER", "PTEDUCAT", "APOE4", cog_vars, mem_vars)

fit_ctl   <- trainControl(method = "repeatedcv", number = 5, repeats = 3,
                          classProbs = TRUE, summaryFunction = twoClassSummary)

### Adaboost
tune_grid <- expand.grid(mfinal = seq(5, 50, 15), maxdepth = c(1, 3, 5))
res_tst   <- vector("list", length = 5)

# Base
f_aboost_base   <- here("data/rds/adni-bl_model-conv3_adaboost_base.rds")
if (file.exists(f_aboost_base)) {
#if (FALSE) {
  aboost_base   <- read_rds(f_aboost_base)
} else {
  cols          <- c("CONV_3Y", covars)
  aboost_base   <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "AdaBag",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(aboost_base, f_aboost_base)
}

preds_base_dt   <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(aboost_base, dt_test, type = "prob"),
                              pred  = predict(aboost_base, dt_test))

res_tst[[1]]    <- twoClassSummary(preds_base_dt,
                                   lev = levels(preds_base_dt$obs))
preds_base      <- fifelse(preds_base_dt$Progressor > .5, 1, 0)

# HC
f_aboost_hc     <- here("data/rds/adni-bl_model-conv3_adaboost_hc.rds")
if (file.exists(f_aboost_hc)) {
#if (FALSE) {
  aboost_hc     <- read_rds(f_aboost_hc)
} else {
  cols          <- c("CONV_3Y", "HC", covars)
  aboost_hc     <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "AdaBag",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(aboost_hc, f_aboost_hc)
}

preds_hc_dt     <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(aboost_hc, dt_test, type = "prob"),
                              pred  = predict(aboost_hc, dt_test))

res_tst[[2]]    <- twoClassSummary(preds_hc_dt, lev = levels(preds_hc_dt$obs))
preds_hc        <- fifelse(preds_hc_dt$Progressor > .5, 1, 0)

# HC_stx
f_aboost_hcstx  <- here("data/rds/adni-bl_model-conv3_adaboost_hc_stx.rds")
if (file.exists(f_aboost_hcstx)) {
#if (FALSE) {
  aboost_hcstx  <- read_rds(f_aboost_hcstx)
} else {
  cols          <- c("CONV_3Y", "HC_stx", covars)
  aboost_hcstx  <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "AdaBag",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(aboost_hcstx, f_aboost_hcstx)
}

preds_hcstx_dt  <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(aboost_hcstx, dt_test, type = "prob"),
                              pred  = predict(aboost_hcstx, dt_test))

res_tst[[3]]    <- twoClassSummary(preds_hcstx_dt,
                                   lev = levels(preds_hcstx_dt$obs))
preds_hcstx     <- fifelse(preds_hcstx_dt$Progressor > .5, 1, 0)

# HC_norm
f_aboost_hcnorm <- here("data/rds/adni-bl_model-conv3_adaboost_hc_norm.rds")
if (file.exists(f_aboost_hcnorm)) {
#if (FALSE) {
  aboost_hcnorm <- read_rds(f_aboost_hcnorm)
} else {
  cols          <- c("CONV_3Y", "HC_norm", covars)
  aboost_hcnorm <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "AdaBag",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(aboost_hcnorm, f_aboost_hcnorm)
}

preds_hcnorm_dt <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(aboost_hcnorm, dt_test, type = "prob"),
                              pred  = predict(aboost_hcnorm, dt_test))

res_tst[[4]]    <- twoClassSummary(preds_hcnorm_dt,
                                   lev = levels(preds_hcnorm_dt$obs))
preds_hcnorm    <- fifelse(preds_hcnorm_dt$Progressor > .5, 1, 0)

# HVR
f_aboost_hvr    <- here("data/rds/adni-bl_model-conv3_adaboost_hvr.rds")
if (file.exists(f_aboost_hvr)) {
#if(FALSE) {
  aboost_hvr    <- read_rds(f_aboost_hvr)
} else {
  cols          <- c("CONV_3Y", "HVR", covars)
  aboost_hvr    <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "AdaBag",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(aboost_hvr, f_aboost_hvr)
}

preds_hvr_dt    <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(aboost_hvr, dt_test, type = "prob"),
                              pred  = predict(aboost_hvr, dt_test))

res_tst[[5]]    <- twoClassSummary(preds_hvr_dt,
                                   lev = levels(preds_hvr_dt$obs))
preds_hvr       <- fifelse(preds_hvr_dt$Progressor > .5, 1, 0)

# ROC
preds_roc       <- data.table(OBS     = dt_test[, CONV_3Y],
                              BASE    = preds_base,
                              HC      = preds_hc,
                              HC_stx  = preds_hcstx,
                              HC_norm = preds_hcnorm,
                              HVR     = preds_hvr)
aboost_roc      <- roc(OBS ~ BASE + HC + HC_stx + HC_norm + HVR,
                       data = preds_roc)
rm(preds_base, preds_hc, preds_hcstx, preds_hcnorm, preds_hvr, preds_roc)


## Table of results
results_aboost  <- rbindlist(lapply(res_tst, data.table))
results_aboost[, `:=`(metric = rep(c("ROC", "Sens", "Spec"), times = 5),
                      var    = rep(c("Base", "HC", "HC_stx", "HC_norm", "HVR"),
                                   each = 3))]

results_aboost  <- dcast(results_aboost, ... ~ metric, value.var = "V1")
rm(res_tst)
fwrite(results_aboost,
       here("data/derivatives/adni-bl_conv3_adaboost_performance.csv"))

## ROCplots
auc_aboost      <- data.table(x     = rep(1, 5),
                              y     = rep(0, 5),
                              name  = names(aboost_roc),
                              auc   = map_dbl(aboost_roc, pluck("auc")))

png(here("plots/adni-bl_conv3_adaboost_roc.png"),
    width = 15, height = 4, units = "in", res = 300)
g <- ggroc(aboost_roc, aes = "group", legacy.axes = TRUE) +
  theme_linedraw(base_size = 15) +
  facet_grid(~name) +
  geom_label(data = auc_aboost, hjust = "inward", vjust = "inward", size = 6,
             aes(x = x, y = y, label = paste("AUC:", round(auc, 3)))) +
  ggtitle("Adaboost | ROC:\nConv-3Y ~ <HC> + Age + Sex + Education + Cog.PC1 + Mem.PC1 + APOE4")
print(g)
dev.off()


### SVM Linear with class weights
tune_grid       <- expand.grid(cost   = c(5, 10, 25, 50),
                               weight = seq(3, 7, 2))
res_tst         <- vector("list", length = 5)

# Base
f_svm_base      <- here("data/rds/adni-bl_model-conv3_svm_base.rds")
if (file.exists(f_svm_base)) {
#if (FALSE) {
  svm_base      <- read_rds(f_svm_base)
} else {
  cols          <- c("CONV_3Y", covars)
  svm_base      <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "svmLinearWeights",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(svm_base, f_svm_base)
}

preds_base_dt   <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(svm_base, dt_test, type = "prob"),
                              pred  = predict(svm_base, dt_test))

res_tst[[1]]    <- twoClassSummary(preds_base_dt,
                                   lev = levels(preds_hc_dt$obs))
preds_base      <- fifelse(preds_base_dt$Progressor > .5, 1, 0)

# HC
f_svm_hc        <- here("data/rds/adni-bl_model-conv3_svm_hc.rds")
if (file.exists(f_svm_hc)) {
#if (FALSE) {
  svm_hc        <- read_rds(f_svm_hc)
} else {
  cols          <- c("CONV_3Y", "HC", covars)
  svm_hc        <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "svmLinearWeights",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(svm_hc, f_svm_hc)
}

preds_hc_dt     <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(svm_hc, dt_test, type = "prob"),
                              pred  = predict(svm_hc, dt_test))

res_tst[[2]]    <- twoClassSummary(preds_hc_dt, lev = levels(preds_hc_dt$obs))
preds_hc        <- fifelse(preds_hc_dt$Progressor > .5, 1, 0)

# HC_stx
f_svm_hcstx     <- here("data/rds/adni-bl_model-conv3_svm_hc_stx.rds")
if (file.exists(f_svm_hcstx)) {
#if (FALSE) {
  svm_hcstx     <- read_rds(f_svm_hcstx)
} else {
  cols          <- c("CONV_3Y", "HC_stx", covars)
  svm_hcstx     <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "svmLinearWeights",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(svm_hcstx, f_svm_hcstx)
}

preds_hcstx_dt  <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(svm_hcstx, dt_test, type = "prob"),
                              pred  = predict(svm_hcstx, dt_test))

res_tst[[3]]    <- twoClassSummary(preds_hcstx_dt,
                                   lev = levels(preds_hcstx_dt$obs))
preds_hcstx     <- fifelse(preds_hcstx_dt$Progressor > .5, 1, 0)

# HC_norm
f_svm_hcnorm <- here("data/rds/adni-bl_model-conv3_svm_hc_norm.rds")
if (file.exists(f_svm_hcnorm)) {
#if (FALSE) {
  svm_hcnorm <- read_rds(f_svm_hcnorm)
} else {
  cols          <- c("CONV_3Y", "HC_norm", covars)
  svm_hcnorm <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "svmLinearWeights",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(svm_hcnorm, f_svm_hcnorm)
}

preds_hcnorm_dt <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(svm_hcnorm, dt_test, type = "prob"),
                              pred  = predict(svm_hcnorm, dt_test))

res_tst[[4]]    <- twoClassSummary(preds_hcnorm_dt,
                                   lev = levels(preds_hcnorm_dt$obs))
preds_hcnorm    <- fifelse(preds_hcnorm_dt$Progressor > .5, 1, 0)

# HVR
f_svm_hvr    <- here("data/rds/adni-bl_model-conv3_svm_hvr.rds")
if (file.exists(f_svm_hvr)) {
#if(FALSE) {
  svm_hvr    <- read_rds(f_svm_hvr)
} else {
  cols          <- c("CONV_3Y", "HVR", covars)
  svm_hvr    <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "svmLinearWeights",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(svm_hvr, f_svm_hvr)
}

preds_hvr_dt    <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(svm_hvr, dt_test, type = "prob"),
                              pred  = predict(svm_hvr, dt_test))

res_tst[[5]]    <- twoClassSummary(preds_hvr_dt,
                                   lev = levels(preds_hvr_dt$obs))
preds_hvr       <- fifelse(preds_hvr_dt$Progressor > .5, 1, 0)

# ROC
preds_roc       <- data.table(OBS     = dt_test[, CONV_3Y],
                              BASE    = preds_base,
                              HC      = preds_hc,
                              HC_stx  = preds_hcstx,
                              HC_norm = preds_hcnorm,
                              HVR     = preds_hvr)
svm_roc         <- roc(OBS ~ BASE + HC + HC_stx + HC_norm + HVR,
                       data = preds_roc)
rm(preds_base, preds_hc, preds_hcstx, preds_hcnorm, preds_hvr, preds_roc)


## Table of results
results_svm     <- rbindlist(lapply(res_tst, data.table))
results_svm[, `:=`(metric = rep(c("ROC", "Sens", "Spec"), times = 5),
                   var    = rep(c("Base", "HC", "HC_stx", "HC_norm", "HVR"),
                                each = 3))]

results_svm     <- dcast(results_svm, ... ~ metric, value.var = "V1")
rm(res_tst)
fwrite(results_svm,
       here("data/derivatives/adni-bl_conv3_svm_performance.csv"))

## ROCplots
auc_svm         <- data.table(x     = rep(1, 5),
                              y     = rep(0, 5),
                              name  = names(svm_roc),
                              auc   = map_dbl(svm_roc, pluck("auc")))

png(here("plots/adni-bl_conv3_svm_roc.png"),
    width = 15, height = 4, units = "in", res = 300)
g <- ggroc(svm_roc, aes = "group", legacy.axes = TRUE) +
  theme_linedraw(base_size = 15) +
  facet_grid(~name) +
  geom_label(data = auc_svm, hjust = "inward", vjust = "inward", size = 6,
             aes(x = x, y = y, label = paste("AUC:", round(auc, 3)))) +
  ggtitle("SVM | ROC:\nConv-3Y ~ <HC> + Age + Sex + Education + Cog.PC1 + Mem.PC1 + APOE4")
print(g)
dev.off()

### Naive Bayes
tune_grid       <- expand.grid(usekernel  = c(TRUE, FALSE),
                               laplace    = seq(0, 1, .5),
                               adjust     = seq(.75, 1.5, .25))
res_tst         <- vector("list", length = 5)

# Base
f_nbayes_base      <- here("data/rds/adni-bl_model-conv3_nbayes_base.rds")
if (file.exists(f_nbayes_base)) {
#if (FALSE) {
  nbayes_base      <- read_rds(f_nbayes_base)
} else {
  cols          <- c("CONV_3Y", covars)
  nbayes_base      <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "naive_bayes",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(nbayes_base, f_nbayes_base)
}

preds_base_dt   <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(nbayes_base, dt_test, type = "prob"),
                              pred  = predict(nbayes_base, dt_test))

res_tst[[1]]    <- twoClassSummary(preds_base_dt,
                                   lev = levels(preds_hc_dt$obs))
preds_base      <- fifelse(preds_base_dt$Progressor > .5, 1, 0)

# HC
f_nbayes_hc        <- here("data/rds/adni-bl_model-conv3_nbayes_hc.rds")
if (file.exists(f_nbayes_hc)) {
#if (FALSE) {
  nbayes_hc        <- read_rds(f_nbayes_hc)
} else {
  cols          <- c("CONV_3Y", "HC", covars)
  nbayes_hc        <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "naive_bayes",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(nbayes_hc, f_nbayes_hc)
}

preds_hc_dt     <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(nbayes_hc, dt_test, type = "prob"),
                              pred  = predict(nbayes_hc, dt_test))

res_tst[[2]]    <- twoClassSummary(preds_hc_dt, lev = levels(preds_hc_dt$obs))
preds_hc        <- fifelse(preds_hc_dt$Progressor > .5, 1, 0)

# HC_stx
f_nbayes_hcstx     <- here("data/rds/adni-bl_model-conv3_nbayes_hc_stx.rds")
if (file.exists(f_nbayes_hcstx)) {
#if (FALSE) {
  nbayes_hcstx     <- read_rds(f_nbayes_hcstx)
} else {
  cols          <- c("CONV_3Y", "HC_stx", covars)
  nbayes_hcstx     <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "naive_bayes",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(nbayes_hcstx, f_nbayes_hcstx)
}

preds_hcstx_dt  <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(nbayes_hcstx, dt_test, type = "prob"),
                              pred  = predict(nbayes_hcstx, dt_test))

res_tst[[3]]    <- twoClassSummary(preds_hcstx_dt,
                                   lev = levels(preds_hcstx_dt$obs))
preds_hcstx     <- fifelse(preds_hcstx_dt$Progressor > .5, 1, 0)

# HC_norm
f_nbayes_hcnorm <- here("data/rds/adni-bl_model-conv3_nbayes_hc_norm.rds")
if (file.exists(f_nbayes_hcnorm)) {
#if (FALSE) {
  nbayes_hcnorm <- read_rds(f_nbayes_hcnorm)
} else {
  cols          <- c("CONV_3Y", "HC_norm", covars)
  nbayes_hcnorm <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "naive_bayes",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(nbayes_hcnorm, f_nbayes_hcnorm)
}

preds_hcnorm_dt <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(nbayes_hcnorm, dt_test, type = "prob"),
                              pred  = predict(nbayes_hcnorm, dt_test))

res_tst[[4]]    <- twoClassSummary(preds_hcnorm_dt,
                                   lev = levels(preds_hcnorm_dt$obs))
preds_hcnorm    <- fifelse(preds_hcnorm_dt$Progressor > .5, 1, 0)

# HVR
f_nbayes_hvr    <- here("data/rds/adni-bl_model-conv3_nbayes_hvr.rds")
if (file.exists(f_nbayes_hvr)) {
#if(FALSE) {
  nbayes_hvr    <- read_rds(f_nbayes_hvr)
} else {
  cols          <- c("CONV_3Y", "HVR", covars)
  nbayes_hvr    <- train(CONV_3Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "naive_bayes",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(nbayes_hvr, f_nbayes_hvr)
}

preds_hvr_dt    <- data.frame(obs   = dt_test[, CONV_3Y],
                              predict(nbayes_hvr, dt_test, type = "prob"),
                              pred  = predict(nbayes_hvr, dt_test))

res_tst[[5]]    <- twoClassSummary(preds_hvr_dt,
                                   lev = levels(preds_hvr_dt$obs))
preds_hvr       <- fifelse(preds_hvr_dt$Progressor > .5, 1, 0)

# ROC
preds_roc       <- data.table(OBS     = dt_test[, CONV_3Y],
                              BASE    = preds_base,
                              HC      = preds_hc,
                              HC_stx  = preds_hcstx,
                              HC_norm = preds_hcnorm,
                              HVR     = preds_hvr)
nbayes_roc      <- roc(OBS ~ BASE + HC + HC_stx + HC_norm + HVR,
                       data = preds_roc)
rm(preds_base, preds_hc, preds_hcstx, preds_hcnorm, preds_hvr, preds_roc)


## Table of results
results_nbayes  <- rbindlist(lapply(res_tst, data.table))
results_nbayes[, `:=`(metric = rep(c("ROC", "Sens", "Spec"), times = 5),
                   var    = rep(c("Base", "HC", "HC_stx", "HC_norm", "HVR"),
                                each = 3))]

results_nbayes  <- dcast(results_nbayes, ... ~ metric, value.var = "V1")
rm(res_tst)
fwrite(results_nbayes,
       here("data/derivatives/adni-bl_conv3_nbayes_performance.csv"))

## ROCplots
auc_nbayes      <- data.table(x     = rep(1, 5),
                              y     = rep(0, 5),
                              name  = names(nbayes_roc),
                              auc   = map_dbl(nbayes_roc, pluck("auc")))

png(here("plots/adni-bl_conv3_nbayes_roc.png"),
    width = 15, height = 4, units = "in", res = 300)
g <- ggroc(nbayes_roc, aes = "group", legacy.axes = TRUE) +
  theme_linedraw(base_size = 15) +
  facet_grid(~name) +
  geom_label(data = auc_nbayes, hjust = "inward", vjust = "inward", size = 6,
             aes(x = x, y = y, label = paste("AUC:", round(auc, 3)))) +
  ggtitle("Naive Bayes | ROC:\nConv-3Y ~ <HC> + Age + Sex + Education + Cog.PC1 + Mem.PC1 + APOE4")
print(g)
dev.off()
