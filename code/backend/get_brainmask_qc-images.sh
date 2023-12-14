#!/usr/bin/env bash

### Obtain Brain Mask QC images from Mahsa's LP2013 directory

set -xue
BASE_DIR=/ipl/ipl27/sfernandez/hvr_adni
ADNI_DIR=/ipl/scratch15/Mahsa/ADNI/LP_2013

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
