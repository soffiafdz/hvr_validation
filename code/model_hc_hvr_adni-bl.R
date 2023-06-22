#!/usr/bin/env Rscript

library(here)
library(readr)
library(data.table)
library(caret)


# Read RDS objects
adnimerge     <- here("data/rds/adnimerge_baseline.rds") |> read_rds()
volumes       <- here("data/rds/adni-bl_volumes_hc-stx-norm-nat_hvr.rds") |>
                read_rds()

# Merge
hc_dt         <- volumes[METHOD == "cnn"
                         ][adnimerge, on = "PTID",
                         .(PTID, DX, AGE, PTGENDER, PTEDUCAT,
                           ABETA, PIB, AV45, APOE4, ADAS13, MMSE, CDRSB,
                           RAVLT_immediate, RAVLT_learning, RAVLT_forgetting,
                           Left_HC        = HC_l,
                           Left_HC_stx    = HC_stx_l,
                           Left_HC_norm   = HC_norm_l,
                           Left_HVR       = HVR_l,
                           Right_HC       = HC_r,
                           Right_HC_stx   = HC_stx_r,
                           Right_HC_norm  = HC_norm_r,
                           Right_HVR      = HVR_r)
                         ][!is.na(Left_HC)]

## ApoE4: binarized
hc_dt[APOE4 == 0, APOE4_bin := 0]
hc_dt[APOE4 != 0, APOE4_bin := 1]

## ABeta positivity:
## SUVR       >  1.11 AV45 PET
## SUVR       >  1.2  Pitts compound-B PET
## CSF AB1-42 <= 980pg/ml
hc_dt[ABETA == ">1700", ABETA := "1700"][, ABETA := as.numeric(ABETA)]
hc_dt[ABETA != ""  | !is.na(PIB) | !is.na(AV45), ABETA_pos := "Negative"]
hc_dt[ABETA <= 980 | PIB > 1.11  | AV45 > 1.2  , ABETA_pos := "Positive"]

# Keep just Amy+ and APOE4 data subs
amy_positive    <- hc_dt[ABETA_pos == "Positive" & DX != "" &
                         !is.na(APOE4_bin) &
                         !is.na(ADAS13) &
                         !is.na(CDRSB) &
                         !is.na(MMSE) &
                         !is.na(RAVLT_immediate) &
                         !is.na(RAVLT_learning) &
                         !is.na(RAVLT_forgetting)]
#rm(hc_dt)

# Wide -> long
amy_positive    <- melt(amy_positive, measure.vars = patterns("Right|Left"))
amy_positive[variable %like% "Left", SIDE := "Left"]
amy_positive[variable %like% "Right", SIDE := "Right"]
amy_positive[, variable := stringr::str_extract(variable, "(?<=_).*")]
amy_positive    <- dcast(amy_positive, ... ~ variable, value.var = "value")

## CNN: HC, HC_stx, HC_norm, HVR
# Train/test split
set.seed(1618)
train_perc      <- .8
train_idx       <- amy_positive[SIDE == "Right",
                                .SD[sample(1:.N, .N * train_perc)],
                                "DX"][, PTID]
#train_idx       <- createDataPartition(amy_positive[, HC],
                                       #p = .8, list = FALSE)
train_dt        <- amy_positive[PTID %in% train_idx]
test_dt         <- amy_positive[!PTID %in% train_idx]
rm(train_idx)

# Z-scoring using training data
train_dt[, `:=`(HC_z      = scale(HC),
                HC_stx_z  = scale(HC_stx),
                HC_norm_z = scale(HC_norm),
                HVR_z     = scale(HVR))]

hc_mean         <- train_dt[, mean(HC)]
hc_std          <- train_dt[, sd(HC)]
hc_stx_mean     <- train_dt[, mean(HC_stx)]
hc_stx_std      <- train_dt[, sd(HC_stx)]
hc_norm_mean    <- train_dt[, mean(HC_norm)]
hc_norm_std     <- train_dt[, sd(HC_norm)]
hvr_mean        <- train_dt[, mean(HVR)]
hvr_std         <- train_dt[, sd(HVR)]

test_dt[, `:=`(HC_z      = (HC - hc_mean) / hc_std,
               HC_stx_z  = (HC_stx - hc_stx_mean) / hc_stx_std,
               HC_norm_z = (HC_norm - hc_norm_mean) / hc_norm_std,
               HVR_z     = (HVR - hvr_mean) / hvr_std)]


# PCA Cognition
cog_vars      <- c("ADAS13", "CDRSB", "MMSE")
pca_cog       <- prcomp(train_dt[, ..cog_vars], scale = T)
train_dt[, COG_PC1 := pca_cog$x[, "PC1"]]
test_dt[, COG_PC1 := predict(pca_cog, test_dt[, ..cog_vars])[, "PC1"]]

# PCA Memory
mem_vars      <- paste("RAVLT",
                       c("immediate", "learning", "forgetting"),
                       sep = "_")
pca_mem       <- prcomp(train_dt[, ..mem_vars], scale = T)
train_dt[, MEM_PC1 := pca_mem$x[, "PC1"]]
test_dt[, MEM_PC1 := predict(pca_mem, test_dt[, ..mem_vars])[, "PC1"]]

## HC, HC_norm, HC_stx, HVR X lm, gam = 8
res_trn       <- vector("list", length = 8)
res_tst       <- vector("list", length = 8)
fit_control   <- trainControl(method = "repeatedcv", number = 5, repeats = 3)
tune_grid     <- expand.grid(span = seq(0.1, 0.6, len = 10), degree = 1:2)
covars <- c("SIDE", "DX", "AGE", "PTGENDER", "PTEDUCAT",
            "COG_PC1", "MEM_PC1", "APOE4_bin")

# Linear Models
set.seed(987)
cols          <- c("HC_z", covars)
lm_hc         <- train(HC_z ~ ., data = train_dt[, ..cols],
                       method = "lm", trControl = fit_control)
res_trn[[1]]  <- lm_hc$results[2:4]
res_tst[[1]]  <- postResample(predict(lm_hc, test_dt), test_dt[, HC_z])

cols          <- c("HC_stx_z", covars)
lm_hc_stx     <- train(HC_stx_z ~ ., data = train_dt[, ..cols],
                       method = "lm", trControl = fit_control)
res_trn[[2]]  <- lm_hc_stx$results[2:4]
res_tst[[2]]  <- postResample(predict(lm_hc_stx, test_dt),
                              test_dt[, HC_stx_z])

cols          <- c("HC_norm_z", covars)
lm_hc_norm    <- train(HC_norm_z ~ ., data = train_dt[, ..cols],
                       method = "lm", trControl = fit_control)
res_trn[[3]]  <- lm_hc_norm$results[2:4]
res_tst[[3]]  <- postResample(predict(lm_hc_norm, test_dt),
                              test_dt[, HC_norm_z])

cols          <- c("HVR_z", covars)
lm_hvr        <- train(HVR_z ~ ., data = train_dt[, ..cols],
                       method = "lm", trControl = fit_control)
res_trn[[4]]  <- lm_hvr$results[2:4]
res_tst[[4]]  <- postResample(predict(lm_hvr, test_dt),
                              test_dt[, HVR_z])

# General Additive Models — LOESS
set.seed(987)
cols          <- c("HC_z", covars)
glm_hc        <- train(HC_z ~ ., data = train_dt[, ..cols],
                       method = "gam", trControl = fit_control)
res_trn[[5]]  <- data.table(glm_hc$results)[select == glm_hc$bestTune[[1]],
                                            .(RMSE, Rsquared, MAE)]
res_tst[[5]]  <- postResample(predict(glm_hc, test_dt), test_dt[, HC_z])

cols          <- c("HC_stx_z", covars)
glm_hc_stx    <- train(HC_stx_z ~ ., data = train_dt[, ..cols],
                       method = "gam", trControl = fit_control)
res_trn[[6]]  <- data.table(glm_hc_stx$results)[select == glm_hc_stx$bestTune[[1]],
                                                .(RMSE, Rsquared, MAE)]
res_tst[[6]]  <- postResample(predict(glm_hc_stx, test_dt),
                              test_dt[, HC_stx_z])

cols          <- c("HC_norm_z", covars)
glm_hc_norm   <- train(HC_norm_z ~ ., data = train_dt[, ..cols],
                       method = "gam", trControl = fit_control)
res_trn[[7]]  <- data.table(glm_hc_norm$results)[select == glm_hc_norm$bestTune[[1]],
                                                 .(RMSE, Rsquared, MAE)]
res_tst[[7]]  <- postResample(predict(glm_hc_norm, test_dt),
                              test_dt[, HC_norm_z])

cols          <- c("HVR_z", covars)
glm_hvr       <- train(HVR_z ~ ., data = train_dt[, ..cols],
                       method = "gam", trControl = fit_control)
res_trn[[8]]  <- data.table(glm_hvr$results)[select == glm_hvr$bestTune[[1]],
                                             .(RMSE, Rsquared, MAE)]
res_tst[[8]]  <- postResample(predict(glm_hvr, test_dt),
                              test_dt[, HVR_z])
rm(cols)

## Concatenate Result metrics
# CV
res_trn_dt    <- rbindlist(res_trn)
setnames(res_trn_dt, c("rmse", "r2", "mae"))
res_trn_dt[, `:=`(data    = rep("train", each = 8),
                  #side    = rep(c("left", "right"), each = 8),
                  model   = rep(c("lm", "gam"), each = 4),
                  var     = rep(c("hc", "hc_stx", "hc_norm", "hvr"),
                                times = 2))]
#setcolorder(res_trn_dt, c(4:6, 1:3))

res_tst_dt    <- rbindlist(lapply(res_tst, data.table))
res_tst_dt[, `:=`(data    = rep("test", each = 24),
                  #side    = rep(c("left", "right"), each = 24),
                  model   = rep(c("lm", "gam"), each = 12),
                  var     = rep(c("hc", "hc_stx", "hc_norm", "hvr"),
                                each = 3, times = 2),
                  metric  = rep(c("rmse", "r2", "mae"), times = 8))]
res_tst_dt    <- dcast(res_tst_dt, ... ~ metric, value.var = "V1")

results_dt    <- rbindlist(list(res_trn_dt, res_tst_dt), use.names = TRUE)
rm(res_trn, res_trn_dt, res_tst, res_tst_dt)
fwrite(results_dt, here("data/derivatives/adni-bl_model_fits.csv"))
