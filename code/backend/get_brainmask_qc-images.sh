#!/usr/bin/env bash

### Obtain brain-mask QC images from the preprocessed ADNI directory

set -xue
BASE_DIR="${HVR_ADNI_DIR:?path of the hvr_adni project}"
ADNI_DIR="${ADNI_PREPROC_DIR:?path of the preprocessed ADNI data}"

LIST=${BASE_DIR}/lists/adni_baseline.lst

OUT_DIR=${BASE_DIR}/plots/qc_adni/skull_masks
[[ -d $OUT_DIR ]] || mkdir $OUT_DIR

MISSING_SUBS=${BASE_DIR}/lists/missing_skull_masks_qc.lst

mapfile -t IDS < $LIST

for id in ${IDS[@]}
do
	sub=$(echo $id | cut -d, -f1)
	sess=$(echo $id | cut -d, -f2)

	qc_img=${ADNI_DIR}/${sub}/qc/qc_stx2_mask_${sub}_${sess}.jpg
	cp -u $qc_img $OUT_DIR

	[[ -f ${OUT_DIR}/$(basename $qc_img) ]] \
		|| printf "%s,%s\n" $sub $sess >> $MISSING_SUBS
done
