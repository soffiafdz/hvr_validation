#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(caret)
library(pROC)
library(purrr)

# Read RDS objects
adnimerge <- here("data/rds/adnimerge_baseline.rds") |> read_rds()
volumes   <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds") |>
            read_rds()

# Merge
DT        <- volumes[METHOD == "cnn"
                   ][adnimerge, on = "PTID",
                     .(PTID, DX, AGE, PTGENDER, PTEDUCAT, CONV_5Y,
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
amy_posit <- DT[!is.na(CONV_5Y) &
                !is.na(ADAS13) &
                DX %in% c("CN", "MCI") &
                ABETA_pos == "Positive"]

amy_posit[, CONV_5Y := factor(CONV_5Y, levels = c(0, 1),
                              labels = c("Stable", "Progressor"))]

## 70% train/test split
set.seed(1618)
train_idx <- createDataPartition(amy_posit[!is.na(CONV_5Y), CONV_5Y],
                                 p = .8, list = FALSE)

dt_train  <- amy_posit[!is.na(CONV_5Y)][ train_idx]
dt_test   <- amy_posit[!is.na(CONV_5Y)][-train_idx]
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
f_aboost_base   <- here("data/rds/adni-bl_model-conv5_adaboost_base.rds")
if (file.exists(f_aboost_base)) {
#if (FALSE) {
  aboost_base   <- read_rds(f_aboost_base)
} else {
  cols          <- c("CONV_5Y", covars)
  aboost_base   <- train(CONV_5Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "AdaBag",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(aboost_base, f_aboost_base)
}

preds_base_dt   <- data.frame(obs   = dt_test[, CONV_5Y],
                              predict(aboost_base, dt_test, type = "prob"),
                              pred  = predict(aboost_base, dt_test))

res_tst[[1]]    <- twoClassSummary(preds_base_dt,
                                   lev = levels(preds_base_dt$obs))
preds_base      <- fifelse(preds_base_dt$Progressor > .5, 1, 0)

# HC
f_aboost_hc     <- here("data/rds/adni-bl_model-conv5_adaboost_hc.rds")
if (file.exists(f_aboost_hc)) {
#if (FALSE) {
  aboost_hc     <- read_rds(f_aboost_hc)
} else {
  cols          <- c("CONV_5Y", "HC", covars)
  aboost_hc     <- train(CONV_5Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "AdaBag",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(aboost_hc, f_aboost_hc)
}

preds_hc_dt     <- data.frame(obs   = dt_test[, CONV_5Y],
                              predict(aboost_hc, dt_test, type = "prob"),
                              pred  = predict(aboost_hc, dt_test))

res_tst[[2]]    <- twoClassSummary(preds_hc_dt, lev = levels(preds_hc_dt$obs))
preds_hc        <- fifelse(preds_hc_dt$Progressor > .5, 1, 0)

# HC_stx
f_aboost_hcstx  <- here("data/rds/adni-bl_model-conv5_adaboost_hc_stx.rds")
if (file.exists(f_aboost_hcstx)) {
#if (FALSE) {
  aboost_hcstx  <- read_rds(f_aboost_hcstx)
} else {
  cols          <- c("CONV_5Y", "HC_stx", covars)
  aboost_hcstx  <- train(CONV_5Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "AdaBag",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(aboost_hcstx, f_aboost_hcstx)
}

preds_hcstx_dt  <- data.frame(obs   = dt_test[, CONV_5Y],
                              predict(aboost_hcstx, dt_test, type = "prob"),
                              pred  = predict(aboost_hcstx, dt_test))

res_tst[[3]]    <- twoClassSummary(preds_hcstx_dt,
                                   lev = levels(preds_hcstx_dt$obs))
preds_hcstx     <- fifelse(preds_hcstx_dt$Progressor > .5, 1, 0)

# HC_norm
f_aboost_hcnorm <- here("data/rds/adni-bl_model-conv5_adaboost_hc_norm.rds")
if (file.exists(f_aboost_hcnorm)) {
#if (FALSE) {
  aboost_hcnorm <- read_rds(f_aboost_hcnorm)
} else {
  cols          <- c("CONV_5Y", "HC_norm", covars)
  aboost_hcnorm <- train(CONV_5Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "AdaBag",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(aboost_hcnorm, f_aboost_hcnorm)
}

preds_hcnorm_dt <- data.frame(obs   = dt_test[, CONV_5Y],
                              predict(aboost_hcnorm, dt_test, type = "prob"),
                              pred  = predict(aboost_hcnorm, dt_test))

res_tst[[4]]    <- twoClassSummary(preds_hcnorm_dt,
                                   lev = levels(preds_hcnorm_dt$obs))
preds_hcnorm    <- fifelse(preds_hcnorm_dt$Progressor > .5, 1, 0)

# HVR
f_aboost_hvr    <- here("data/rds/adni-bl_model-conv5_adaboost_hvr.rds")
if (file.exists(f_aboost_hvr)) {
#if(FALSE) {
  aboost_hvr    <- read_rds(f_aboost_hvr)
} else {
  cols          <- c("CONV_5Y", "HVR", covars)
  aboost_hvr    <- train(CONV_5Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "AdaBag",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(aboost_hvr, f_aboost_hvr)
}

preds_hvr_dt    <- data.frame(obs   = dt_test[, CONV_5Y],
                              predict(aboost_hvr, dt_test, type = "prob"),
                              pred  = predict(aboost_hvr, dt_test))

res_tst[[5]]    <- twoClassSummary(preds_hvr_dt,
                                   lev = levels(preds_hvr_dt$obs))
preds_hvr       <- fifelse(preds_hvr_dt$Progressor > .5, 1, 0)

# ROC
preds_roc       <- data.table(OBS     = dt_test[, CONV_5Y],
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
#results_aboost[, `:=`(covars = rep(c("Demog", "Extended"), each = 12),
                      #metric = rep(c("ROC", "Sens", "Spec"), times = 8),
                      #var    = rep(c("HC", "HC_stx", "HC_norm", "HVR"),
                                   #times = 2, each = 3))]

#results_aboost  <- dcast(results_aboost, ... ~ metric, value.var = "V1")
#rm(res_tst)
#fwrite(results_aboost,
       #here("data/derivatives/adni-bl_adaboost_performance.csv"))

## ROCplots
#auc_3y_base     <- data.table(x     = rep(0, 4),
                              #y     = rep(1, 4),
                              #name  = names(conv3y_base_roc),
                              #auc   = map_dbl(conv3y_base_roc, pluck("auc")))

#png(here("plots/adni-bl_conv_3y_base_roc.png"),
    #width = 20, height = 10, units = "in", res = 300)
#g <- ggroc(conv3y_base_roc, aes = "group", legacy.axes = TRUE) +
  #theme_linedraw(base_size = 24) +
  #facet_grid(~name) +
  #geom_label(data = auc_3y_base, hjust = "inward", size = 6,
             #aes(x = x, y = y, label = paste("AUC:", round(auc, 3)))) +
  #ggtitle("Adaboost | ROC: Conv-3Y ~ <HC> + Age + Sex + Education")
#print(g)
#dev.off()


### SVM Linear with class weights
tune_grid       <- expand.grid(cost   = c(5, 10, 25, 50),
                               weight = seq(3, 7, 2))
res_tst         <- vector("list", length = 5)

# Base
f_svm_base      <- here("data/rds/adni-bl_model-conv5_svm_base.rds")
if (file.exists(f_svm_base)) {
#if (FALSE) {
  svm_base      <- read_rds(f_svm_base)
} else {
  cols          <- c("CONV_5Y", covars)
  svm_base      <- train(CONV_5Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "svmLinearWeights",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(svm_base, f_svm_base)
}

preds_base_dt   <- data.frame(obs   = dt_test[, CONV_5Y],
                              predict(svm_base, dt_test, type = "prob"),
                              pred  = predict(svm_base, dt_test))

res_tst[[1]]    <- twoClassSummary(preds_base_dt,
                                   lev = levels(preds_hc_dt$obs))
preds_base      <- fifelse(preds_base_dt$Progressor > .5, 1, 0)

# HC
f_svm_hc        <- here("data/rds/adni-bl_model-conv5_svm_hc.rds")
if (file.exists(f_svm_hc)) {
#if (FALSE) {
  svm_hc        <- read_rds(f_svm_hc)
} else {
  cols          <- c("CONV_5Y", "HC", covars)
  svm_hc        <- train(CONV_5Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "svmLinearWeights",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(svm_hc, f_svm_hc)
}

preds_hc_dt     <- data.frame(obs   = dt_test[, CONV_5Y],
                              predict(svm_hc, dt_test, type = "prob"),
                              pred  = predict(svm_hc, dt_test))

res_tst[[2]]    <- twoClassSummary(preds_hc_dt, lev = levels(preds_hc_dt$obs))
preds_hc        <- fifelse(preds_hc_dt$Progressor > .5, 1, 0)

# HC_stx
f_svm_hcstx     <- here("data/rds/adni-bl_model-conv5_svm_hc_stx.rds")
if (file.exists(f_svm_hcstx)) {
#if (FALSE) {
  svm_hcstx     <- read_rds(f_svm_hcstx)
} else {
  cols          <- c("CONV_5Y", "HC_stx", covars)
  svm_hcstx     <- train(CONV_5Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "svmLinearWeights",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(svm_hcstx, f_svm_hcstx)
}

preds_hcstx_dt  <- data.frame(obs   = dt_test[, CONV_5Y],
                              predict(svm_hcstx, dt_test, type = "prob"),
                              pred  = predict(svm_hcstx, dt_test))

res_tst[[3]]    <- twoClassSummary(preds_hcstx_dt,
                                   lev = levels(preds_hcstx_dt$obs))
preds_hcstx     <- fifelse(preds_hcstx_dt$Progressor > .5, 1, 0)

# HC_norm
f_svm_hcnorm <- here("data/rds/adni-bl_model-conv5_svm_hc_norm.rds")
if (file.exists(f_svm_hcnorm)) {
#if (FALSE) {
  svm_hcnorm <- read_rds(f_svm_hcnorm)
} else {
  cols          <- c("CONV_5Y", "HC_norm", covars)
  svm_hcnorm <- train(CONV_5Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "svmLinearWeights",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(svm_hcnorm, f_svm_hcnorm)
}

preds_hcnorm_dt <- data.frame(obs   = dt_test[, CONV_5Y],
                              predict(svm_hcnorm, dt_test, type = "prob"),
                              pred  = predict(svm_hcnorm, dt_test))

res_tst[[4]]    <- twoClassSummary(preds_hcnorm_dt,
                                   lev = levels(preds_hcnorm_dt$obs))
preds_hcnorm    <- fifelse(preds_hcnorm_dt$Progressor > .5, 1, 0)

# HVR
f_svm_hvr    <- here("data/rds/adni-bl_model-conv5_svm_hvr.rds")
if (file.exists(f_svm_hvr)) {
#if(FALSE) {
  svm_hvr    <- read_rds(f_svm_hvr)
} else {
  cols          <- c("CONV_5Y", "HVR", covars)
  svm_hvr    <- train(CONV_5Y ~ .,
                         data       = dt_train[, ..cols],
                         method     = "svmLinearWeights",
                         trControl  = fit_ctl,
                         tuneGrid   = tune_grid,
                         metric     = "ROC")
  write_rds(svm_hvr, f_svm_hvr)
}

preds_hvr_dt    <- data.frame(obs   = dt_test[, CONV_5Y],
                              predict(svm_hvr, dt_test, type = "prob"),
                              pred  = predict(svm_hvr, dt_test))

res_tst[[5]]    <- twoClassSummary(preds_hvr_dt,
                                   lev = levels(preds_hvr_dt$obs))
preds_hvr       <- fifelse(preds_hvr_dt$Progressor > .5, 1, 0)

# ROC
preds_roc       <- data.table(OBS     = dt_test[, CONV_5Y],
                              BASE    = preds_base,
                              HC      = preds_hc,
                              HC_stx  = preds_hcstx,
                              HC_norm = preds_hcnorm,
                              HVR     = preds_hvr)
svm_roc      <- roc(OBS ~ BASE + HC + HC_stx + HC_norm + HVR,
                       data = preds_roc)
rm(preds_base, preds_hc, preds_hcstx, preds_hcnorm, preds_hvr, preds_roc)


## Table of results
results_svm  <- rbindlist(lapply(res_tst, data.table))

### 3 Years Conversion
## 70% train/test split
#set.seed(1618)
#train_idx         <- createDataPartition(DT[!is.na(CONV_3Y), CONV_3Y],
                                         #p = .8, list = FALSE)

#dt_3y_train       <- DT[!is.na(CONV_3Y)][ train_idx]
#dt_3y_test        <- DT[!is.na(CONV_3Y)][-train_idx]
#rm(train_idx)

#res_tst           <- vector("list", length = 8)

### Demographic covariates
## HC
#f_svm_3y_hc1      <- here("data/rds/adni-bl_model-conv3_svm_hc_base.rds")
##if (file.exists(f_svm_3y_hc1)) {
#if (FALSE) {
  #svm_3y_hc1      <- read_rds(f_svm_3y_hc1)
#} else {
  #cols            <- c("CONV_3Y", "HC", covars_base)
  #svm_3y_hc1      <- train(CONV_3Y ~ .,
                           #data       = dt_3y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_3y_hc1, f_svm_3y_hc1)
#}

#preds_hc_dt       <- data.frame(obs   = dt_3y_test[, CONV_3Y],
                                #predict(svm_3y_hc1, dt_3y_test, type = "prob"),
                                #pred  = predict(svm_3y_hc1, dt_3y_test))

#res_tst[[1]]      <- twoClassSummary(preds_hc_dt,
                                     #lev = levels(preds_hc_dt$obs))
#preds_hc          <- fifelse(preds_hc_dt$Progressor > .5, 1, 0)

## HC_stx
#f_svm_3y_hc_stx1  <- here("data/rds/adni-bl_model-conv3_svm_hc_stx_base.rds")
#if (file.exists(f_svm_3y_hc_stx1)) {
##if (FALSE) {
  #svm_3y_hc_stx1  <- read_rds(f_svm_3y_hc_stx1)
#} else {
  #cols            <- c("CONV_3Y", "HC_stx", covars_base)
  #svm_3y_hc_stx1  <- train(CONV_3Y ~ .,
                           #data       = dt_3y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_3y_hc_stx1, f_svm_3y_hc_stx1)
#}

#preds_hc_stx_dt   <- data.frame(obs   = dt_3y_test[, CONV_3Y],
                                #predict(svm_3y_hc_stx1, dt_3y_test, type = "prob"),
                                #pred  = predict(svm_3y_hc_stx1, dt_3y_test))

#res_tst[[2]]      <- twoClassSummary(preds_hc_stx_dt,
                                     #lev = levels(preds_hc_stx_dt$obs))
#preds_hc_stx      <- fifelse(preds_hc_stx_dt$Progressor > .5, 1, 0)

## HC_norm
#f_svm_3y_hc_norm1 <- here("data/rds/adni-bl_model-conv3_svm_hc_norm_base.rds")
#if (file.exists(f_svm_3y_hc_norm1)) {
##if (FALSE) {
  #svm_3y_hc_norm1 <- read_rds(f_svm_3y_hc_norm1)
#} else {
  #cols            <- c("CONV_3Y", "HC_norm", covars_base)
  #svm_3y_hc_norm1 <- train(CONV_3Y ~ .,
                           #data       = dt_3y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_3y_hc_norm1, f_svm_3y_hc_norm1)
#}

#preds_hc_norm_dt  <- data.frame(obs   = dt_3y_test[, CONV_3Y],
                                #predict(svm_3y_hc_norm1, dt_3y_test, type = "prob"),
                                #pred  = predict(svm_3y_hc_norm1, dt_3y_test))

#res_tst[[3]]      <- twoClassSummary(preds_hc_norm_dt,
                                     #lev = levels(preds_hc_norm_dt$obs))
#preds_hc_norm     <- fifelse(preds_hc_norm_dt$Progressor > .5, 1, 0)

## HVR
#f_svm_3y_hvr1     <- here("data/rds/adni-bl_model-conv3_svm_hvr_base.rds")
#if (file.exists(f_svm_3y_hvr1)) {
##if(FALSE) {
  #svm_3y_hvr1     <- read_rds(f_svm_3y_hvr1)
#} else {
  #cols            <- c("CONV_3Y", "HVR", covars_base)
  #svm_3y_hvr1     <- train(CONV_3Y ~ .,
                           #data       = dt_3y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_3y_hvr1, f_svm_3y_hvr1)
#}

#preds_hvr_dt      <- data.frame(obs   = dt_3y_test[, CONV_3Y],
                                #predict(svm_3y_hvr1, dt_3y_test, type = "prob"),
                                #pred  = predict(svm_3y_hvr1, dt_3y_test))

#res_tst[[4]]      <- twoClassSummary(preds_hvr_dt,
                                     #lev = levels(preds_hvr_dt$obs))
#preds_hvr         <- fifelse(preds_hvr_dt$Progressor > .5, 1, 0)

## ROC
#preds_roc         <- data.table(OBS     = dt_3y_test[, CONV_3Y],
                                #HC      = preds_hc,
                                #HC_stx  = preds_hc_stx,
                                #HC_norm = preds_hc_norm,
                                #HVR     = preds_hvr)
#conv3y_base_roc   <- roc(OBS ~ HC + HC_stx + HC_norm + HVR, data = preds_roc)
#rm(preds_hc, preds_hc_stx, preds_hc_norm, preds_hvr, preds_roc)

### Extended covariates
## HC
#f_svm_3y_hc2      <- here("data/rds/adni-bl_model-conv3_svm_hc_ext.rds")
#if (file.exists(f_svm_3y_hc2)) {
##if (FALSE) {
  #svm_3y_hc2      <- read_rds(f_svm_3y_hc2)
#} else {
  #cols            <- c("CONV_3Y", "HC", covars_ext)
  #svm_3y_hc2      <- train(CONV_3Y ~ .,
                           #data       = dt_3y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_3y_hc2, f_svm_3y_hc2)
#}

#preds_hc_dt       <- data.frame(obs   = dt_3y_test[, CONV_3Y],
                                #predict(svm_3y_hc2, dt_3y_test, type = "prob"),
                                #pred  = predict(svm_3y_hc2, dt_3y_test))

#res_tst[[5]]      <- twoClassSummary(preds_hc_dt,
                                     #lev = levels(preds_hc_dt$obs))
#preds_hc          <- fifelse(preds_hc_dt$Progressor > .5, 1, 0)

## HC_stx
#f_svm_3y_hc_stx2  <- here("data/rds/adni-bl_model-conv3_svm_hc_stx_ext.rds")
#if (file.exists(f_svm_3y_hc_stx2)) {
##if (FALSE) {
  #svm_3y_hc_stx2  <- read_rds(f_svm_3y_hc_stx2)
#} else {
  #cols            <- c("CONV_3Y", "HC_stx", covars_ext)
  #svm_3y_hc_stx2  <- train(CONV_3Y ~ .,
                           #data       = dt_3y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_3y_hc_stx2, f_svm_3y_hc_stx2)
#}

#preds_hc_stx_dt   <- data.frame(obs   = dt_3y_test[, CONV_3Y],
                                #predict(svm_3y_hc_stx2, dt_3y_test, type = "prob"),
                                #pred  = predict(svm_3y_hc_stx2, dt_3y_test))

#res_tst[[6]]      <- twoClassSummary(preds_hc_stx_dt,
                                     #lev = levels(preds_hc_stx_dt$obs))
#preds_hc_stx      <- fifelse(preds_hc_stx_dt$Progressor > .5, 1, 0)

## HC_norm
#f_svm_3y_hc_norm2 <- here("data/rds/adni-bl_model-conv3_svm_hc_norm_ext.rds")
#if (file.exists(f_svm_3y_hc_norm2)) {
##if (FALSE) {
  #svm_3y_hc_norm2 <- read_rds(f_svm_3y_hc_norm2)
#} else {
  #cols            <- c("CONV_3Y", "HC_norm", covars_ext)
  #svm_3y_hc_norm2 <- train(CONV_3Y ~ .,
                           #data       = dt_3y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_3y_hc_norm2, f_svm_3y_hc_norm2)
#}

#preds_hc_norm_dt  <- data.frame(obs   = dt_3y_test[, CONV_3Y],
                                #predict(svm_3y_hc_norm2, dt_3y_test, type = "prob"),
                                #pred  = predict(svm_3y_hc_norm2, dt_3y_test))

#res_tst[[7]]      <- twoClassSummary(preds_hc_norm_dt,
                                     #lev = levels(preds_hc_norm_dt$obs))
#preds_hc_norm     <- fifelse(preds_hc_norm_dt$Progressor > .5, 1, 0)

## HVR
#f_svm_3y_hvr2     <- here("data/rds/adni-bl_model-conv3_svm_hvr_ext.rds")
#if (file.exists(f_svm_3y_hvr2)) {
##if(FALSE) {
  #svm_3y_hvr2     <- read_rds(f_svm_3y_hvr2)
#} else {
  #cols            <- c("CONV_3Y", "HVR", covars_ext)
  #svm_3y_hvr2     <- train(CONV_3Y ~ .,
                           #data       = dt_3y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_3y_hvr2, f_svm_3y_hvr2)
#}

#preds_hvr_dt      <- data.frame(obs   = dt_3y_test[, CONV_3Y],
                                #predict(svm_3y_hvr2, dt_3y_test, type = "prob"),
                                #pred  = predict(svm_3y_hvr2, dt_3y_test))

#res_tst[[8]]      <- twoClassSummary(preds_hvr_dt,
                                     #lev = levels(preds_hvr_dt$obs))
#preds_hvr         <- fifelse(preds_hvr_dt$Progressor > .5, 1, 0)

## ROC
#preds_roc         <- data.table(OBS     = dt_3y_test[, CONV_3Y],
                                #HC      = preds_hc,
                                #HC_stx  = preds_hc_stx,
                                #HC_norm = preds_hc_norm,
                                #HVR     = preds_hvr)
#conv3y_ext_roc    <- roc(OBS ~ HC + HC_stx + HC_norm + HVR, data = preds_roc)
#rm(preds_hc, preds_hc_stx, preds_hc_norm, preds_hvr, preds_roc)

### Table of results
#results_svm_3y    <- rbindlist(lapply(res_tst, data.table))
#results_svm_3y[, `:=`(covars = rep(c("Demog", "Extended"), each = 12),
                      #metric = rep(c("ROC", "Sens", "Spec"), times = 8),
                      #var    = rep(c("HC", "HC_stx", "HC_norm", "HVR"),
                                   #times = 2, each = 3))]

#results_svm_3y    <- dcast(results_svm_3y, ... ~ metric, value.var = "V1")
#rm(res_tst)
#fwrite(results_svm_3y,
       #here("data/derivatives/adni-bl_svm_conv-3y_performance.csv"))

### 5 Years Conversion
## 70% train/test split
#set.seed(1618)
#train_idx         <- createDataPartition(DT[!is.na(CONV_5Y), CONV_5Y],
                                         #p = .8, list = FALSE)

#dt_5y_train       <- DT[!is.na(CONV_5Y)][ train_idx]
#dt_5y_test        <- DT[!is.na(CONV_5Y)][-train_idx]
#rm(train_idx)

#res_tst           <- vector("list", length = 8)

### Demographic covariates
## HC
#f_svm_5y_hc1      <- here("data/rds/adni-bl_model-conv5_svm_hc_base.rds")
#if (file.exists(f_svm_5y_hc1)) {
##if (FALSE) {
  #svm_5y_hc1      <- read_rds(f_svm_5y_hc1)
#} else {
  #cols            <- c("CONV_5Y", "HC", covars_base)
  #svm_5y_hc1      <- train(CONV_5Y ~ .,
                           #data       = dt_5y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_5y_hc1, f_svm_5y_hc1)
#}

#preds_hc_dt       <- data.frame(obs   = dt_5y_test[, CONV_3Y],
                                #predict(svm_5y_hc1, dt_5y_test, type = "prob"),
                                #pred  = predict(svm_5y_hc1, dt_5y_test))

#res_tst[[1]]      <- twoClassSummary(preds_hc_dt,
                                     #lev = levels(preds_hc_dt$obs))
#preds_hc          <- fifelse(preds_hc_dt$Progressor > .5, 1, 0)

## HC_stx
#f_svm_5y_hc_stx1  <- here("data/rds/adni-bl_model-conv5_svm_hc_stx_base.rds")
#if (file.exists(f_svm_5y_hc_stx1)) {
##if (FALSE) {
  #svm_5y_hc_stx1  <- read_rds(f_svm_5y_hc_stx1)
#} else {
  #cols            <- c("CONV_5Y", "HC_stx", covars_base)
  #svm_5y_hc_stx1  <- train(CONV_5Y ~ .,
                           #data       = dt_5y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_5y_hc_stx1, f_svm_5y_hc_stx1)
#}

#preds_hc_stx_dt   <- data.frame(obs   = dt_5y_test[, CONV_3Y],
                                #predict(svm_5y_hc_stx1, dt_5y_test, type = "prob"),
                                #pred  = predict(svm_5y_hc_stx1, dt_5y_test))

#res_tst[[2]]      <- twoClassSummary(preds_hc_stx_dt,
                                     #lev = levels(preds_hc_stx_dt$obs))
#preds_hc_stx      <- fifelse(preds_hc_stx_dt$Progressor > .5, 1, 0)

## HC_norm
#f_svm_5y_hc_norm1 <- here("data/rds/adni-bl_model-conv5_svm_hc_norm_base.rds")
#if (file.exists(f_svm_5y_hc_norm1)) {
##if (FALSE) {
  #svm_5y_hc_norm1 <- read_rds(f_svm_5y_hc_norm1)
#} else {
  #cols            <- c("CONV_3Y", "HC_norm", covars_base)
  #svm_5y_hc_norm1 <- train(CONV_3Y ~ .,
                           #data       = dt_5y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_5y_hc_norm1, f_svm_5y_hc_norm1)
#}

#preds_hc_norm_dt  <- data.frame(obs   = dt_5y_test[, CONV_5Y],
                                #predict(svm_5y_hc_norm1, dt_5y_test, type = "prob"),
                                #pred  = predict(svm_5y_hc_norm1, dt_5y_test))

#res_tst[[3]]      <- twoClassSummary(preds_hc_norm_dt,
                                     #lev = levels(preds_hc_norm_dt$obs))
#preds_hc_norm     <- fifelse(preds_hc_norm_dt$Progressor > .5, 1, 0)

## HVR
#f_svm_5y_hvr1     <- here("data/rds/adni-bl_model-conv5_svm_hvr_base.rds")
#if (file.exists(f_svm_5y_hvr1)) {
##if(FALSE) {
  #svm_5y_hvr1     <- read_rds(f_svm_5y_hvr1)
#} else {
  #cols            <- c("CONV_5Y", "HVR", covars_base)
  #svm_5y_hvr1     <- train(CONV_5Y ~ .,
                           #data       = dt_5y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_5y_hvr1, f_svm_5y_hvr1)
#}

#preds_hvr_dt      <- data.frame(obs   = dt_5y_test[, CONV_5Y],
                                #predict(svm_5y_hvr1, dt_5y_test, type = "prob"),
                                #pred  = predict(svm_5y_hvr1, dt_5y_test))

#res_tst[[4]]      <- twoClassSummary(preds_hvr_dt,
                                     #lev = levels(preds_hvr_dt$obs))
#preds_hvr         <- fifelse(preds_hvr_dt$Progressor > .5, 1, 0)

## ROC
#preds_roc         <- data.table(OBS     = dt_5y_test[, CONV_5Y],
                                #HC      = preds_hc,
                                #HC_stx  = preds_hc_stx,
                                #HC_norm = preds_hc_norm,
                                #HVR     = preds_hvr)
#conv5y_base_roc   <- roc(OBS ~ HC + HC_stx + HC_norm + HVR, data = preds_roc)
#rm(preds_hc, preds_hc_stx, preds_hc_norm, preds_hvr, preds_roc)

### Extended covariates
## HC
#f_svm_5y_hc2      <- here("data/rds/adni-bl_model-conv5_svm_hc_ext.rds")
#if (file.exists(f_svm_5y_hc2)) {
##if (FALSE) {
  #svm_5y_hc2      <- read_rds(f_svm_5y_hc2)
#} else {
  #cols            <- c("CONV_5Y", "HC", covars_ext)
  #svm_5y_hc2      <- train(CONV_5Y ~ .,
                           #data       = dt_5y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_5y_hc2, f_svm_5y_hc2)
#}

#preds_hc_dt       <- data.frame(obs   = dt_5y_test[, CONV_5Y],
                                #predict(svm_5y_hc2, dt_5y_test, type = "prob"),
                                #pred  = predict(svm_5y_hc2, dt_5y_test))

#res_tst[[5]]      <- twoClassSummary(preds_hc_dt,
                                     #lev = levels(preds_hc_dt$obs))
#preds_hc          <- fifelse(preds_hc_dt$Progressor > .5, 1, 0)

## HC_stx
#f_svm_5y_hc_stx2  <- here("data/rds/adni-bl_model-conv5_svm_hc_stx_ext.rds")
#if (file.exists(f_svm_5y_hc_stx2)) {
##if (FALSE) {
  #svm_5y_hc_stx2  <- read_rds(f_svm_5y_hc_stx2)
#} else {
  #cols            <- c("CONV_5Y", "HC_stx", covars_ext)
  #svm_5y_hc_stx2  <- train(CONV_5Y ~ .,
                           #data       = dt_5y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_5y_hc_stx2, f_svm_5y_hc_stx2)
#}

#preds_hc_stx_dt   <- data.frame(obs   = dt_5y_test[, CONV_5Y],
                                #predict(svm_5y_hc_stx2, dt_5y_test, type = "prob"),
                                #pred  = predict(svm_5y_hc_stx2, dt_5y_test))

#res_tst[[6]]      <- twoClassSummary(preds_hc_stx_dt,
                                     #lev = levels(preds_hc_stx_dt$obs))
#preds_hc_stx      <- fifelse(preds_hc_stx_dt$Progressor > .5, 1, 0)

## HC_norm
#f_svm_5y_hc_norm2 <- here("data/rds/adni-bl_model-conv5_svm_hc_norm_ext.rds")
#if (file.exists(f_svm_5y_hc_norm2)) {
##if (FALSE) {
  #svm_5y_hc_norm2 <- read_rds(f_svm_5y_hc_norm2)
#} else {
  #cols            <- c("CONV_5Y", "HC_norm", covars_ext)
  #svm_5y_hc_norm2 <- train(CONV_5Y ~ .,
                           #data       = dt_5y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_5y_hc_norm2, f_svm_5y_hc_norm2)
#}

#preds_hc_norm_dt  <- data.frame(obs   = dt_5y_test[, CONV_5Y],
                                #predict(svm_5y_hc_norm2, dt_5y_test, type = "prob"),
                                #pred  = predict(svm_5y_hc_norm2, dt_5y_test))

#res_tst[[7]]      <- twoClassSummary(preds_hc_norm_dt,
                                     #lev = levels(preds_hc_norm_dt$obs))
#preds_hc_norm     <- fifelse(preds_hc_norm_dt$Progressor > .5, 1, 0)

## HVR
#f_svm_5y_hvr2     <- here("data/rds/adni-bl_model-conv5_svm_hvr_ext.rds")
#if (file.exists(f_svm_5y_hvr2)) {
##if(FALSE) {
  #svm_5y_hvr2     <- read_rds(f_svm_5y_hvr2)
#} else {
  #cols            <- c("CONV_5Y", "HVR", covars_ext)
  #svm_5y_hvr2     <- train(CONV_5Y ~ .,
                           #data       = dt_5y_train[, ..cols],
                           #method     = "svmLinearWeights",
                           #trControl  = fit_control,
                           #tuneGrid   = tune_grid,
                           #metric     = "ROC")
  #write_rds(svm_5y_hvr2, f_svm_5y_hvr2)
#}

#preds_hvr_dt      <- data.frame(obs   = dt_5y_test[, CONV_5Y],
                                #predict(svm_5y_hvr2, dt_5y_test, type = "prob"),
                                #pred  = predict(svm_5y_hvr2, dt_5y_test))

#res_tst[[8]]      <- twoClassSummary(preds_hvr_dt,
                                     #lev = levels(preds_hvr_dt$obs))
#preds_hvr         <- fifelse(preds_hvr_dt$Progressor > .5, 1, 0)

## ROC
#preds_roc         <- data.table(OBS     = dt_5y_test[, CONV_5Y],
                                #HC      = preds_hc,
                                #HC_stx  = preds_hc_stx,
                                #HC_norm = preds_hc_norm,
                                #HVR     = preds_hvr)
#conv5y_ext_roc   <- roc(OBS ~ HC + HC_stx + HC_norm + HVR, data = preds_roc)
#rm(preds_hc, preds_hc_stx, preds_hc_norm, preds_hvr, preds_roc)

### Table of results
#results_svm_5y    <- rbindlist(lapply(res_tst, data.table))
#results_svm_5y[, `:=`(covars = rep(c("Demog", "Extended"), each = 12),
                      #metric = rep(c("ROC", "Sens", "Spec"), times = 8),
                      #var    = rep(c("HC", "HC_stx", "HC_norm", "HVR"),
                                   #times = 2, each = 3))]

#results_svm_5y    <- dcast(results_svm_5y, ... ~ metric, value.var = "V1")
#rm(res_tst)
#fwrite(results_svm_5y,
       #here("data/derivatives/adni-bl_svm_conv-5y_performance.csv"))

#### ROCplots
### SVM_3y_base
##auc_3y_base       <- data.table(x     = rep(0, 4),
                                ##y     = rep(1, 4),
                                ##name  = names(conv3y_base_roc),
                                ##auc   = map_dbl(conv3y_base_roc, pluck("auc")))

##png(here("plots/adni-bl_conv_3y_base_roc.png"),
    ##width = 20, height = 10, units = "in", res = 300)
##g <- ggroc(conv3y_base_roc, aes = "group", legacy.axes = TRUE) +
  ##theme_linedraw(base_size = 24) +
  ##facet_grid(~name) +
  ##geom_label(data = auc_3y_base, hjust = "inward", size = 6,
             ##aes(x = x, y = y, label = paste("AUC:", round(auc, 3)))) +
  ##ggtitle("SVM | ROC: Conv-3Y ~ <HC> + Age + Sex + Education")
##print(g)
##dev.off()

### SVM_3y_ext
##auc_3y_ext        <- data.table(x     = rep(0, 4),
                                ##y     = rep(1, 4),
                                ##name  = names(conv3y_ext_roc),
                                ##auc   = map_dbl(conv3y_ext_roc, pluck("auc")))

##png(here("plots/adni-bl_conv_3y_ext_roc.png"),
    ##width = 20, height = 10, units = "in", res = 300)
##g <- ggroc(conv3y_ext_roc, aes = "group", legacy.axes = TRUE) +
  ##theme_linedraw(base_size = 24) +
  ##facet_grid(~name) +
  ##geom_label(data = auc_3y_ext, hjust = "inward", size = 6,
             ##aes(x = x, y = y, label = paste("AUC:", round(auc, 3)))) +
  ##ggtitle("SVM | ROC: Conv-3Y ~ <HC> + Age + Sex + Education + AMY + APOE4")
##print(g)
##dev.off()

### SVM_5y_base
##auc_5y_base       <- data.table(x     = rep(0, 4),
                                ##y     = rep(1, 4),
                                ##name  = names(conv5y_base_roc),
                                ##auc   = map_dbl(conv5y_base_roc, pluck("auc")))

##png(here("plots/adni-bl_conv_5y_base_roc.png"),
    ##width = 20, height = 10, units = "in", res = 300)
##g <- ggroc(conv5y_base_roc, aes = "group", legacy.axes = TRUE) +
  ##theme_linedraw(base_size = 24) +
  ##facet_grid(~name) +
  ##geom_label(data = auc_5y_base, hjust = "inward", size = 6,
             ##aes(x = x, y = y, label = paste("AUC:", round(auc, 3)))) +
  ##ggtitle("SVM | ROC: Conv-5Y ~ <HC> + Age + Sex + Education")
##print(g)
##dev.off()

### SVM_5y_ext
##auc_5y_ext        <- data.table(x     = rep(0, 4),
                                ##y     = rep(1, 4),
                                ##name  = names(conv5y_ext_roc),
                                ##auc   = map_dbl(conv5y_ext_roc, pluck("auc")))

##png(here("plots/adni-bl_conv_5y_ext_roc.png"),
    ##width = 20, height = 10, units = "in", res = 300)
##g <- ggroc(conv5y_ext_roc, aes = "group", legacy.axes = TRUE) +
  ##theme_linedraw(base_size = 24) +
  ##facet_grid(~name) +
  ##geom_label(data = auc_5y_ext, hjust = "inward", size = 6,
             ##aes(x = x, y = y, label = paste("AUC:", round(auc, 3)))) +
  ##ggtitle("SVM | ROC: Conv-5Y ~ <HC> + Age + Sex + Education + AMY + APOE4")
##print(g)
##dev.off()
