#!/usr/bin/env bash

## Compare overlap similarity between segmentations:
## CNN and manual labels
## XCorrelation

TMPDIR=$(mktemp -d --tmpdir)
trap "rm -rf $TMPDIR" 0 1 2 15

set -ux

# Directories
HERE=/ipl/ipl27/sfernandez/hvr_validation
MANUAL_LABS=${HERE}/data/labels_dorothee/labels_reduced
CNN_LABS=${HERE}/data/orig_segmentations/cnn/hcvc
KAPPAS=${HERE}/data/derivatives/man-seg_kappa_hcvc_cnn.csv

# Recode from 1-4 for comparison
# $1 In $2 Out
recode() {
	itk_resample \
		--labels \
		--clobber \
		--lut-string "11 1; 12 2; 21 3; 22 4" \
		$1 $2
}

# Calculate kappas & dice
compare() {
	minccmp \
		-similarity \
		$1 $2
}

printf "id,label,dice,kappa\n" > $KAPPAS

for segmented in $CNN_LABS/*
do
	bname=$(basename $segmented)
	id=${bname#tal_as_}
	id=${id%%_*}

	segmented_recoded=$TMPDIR/${id}_cnn.mnc
	recode $segmented $segmented_recoded

	manual=$MANUAL_LABS/tal_${id}_hp_csf_ld.mnc
	manual_recoded=$TMPDIR/${id}_manual.mnc
	recode $manual $manual_recoded

	compare \
		$segmented_recoded \
		$manual_recoded | awk -v id=$id \
		'$1 ~ /[1234]/ {print id "," $1 "," $2 "," $6}' >> $KAPPAS
done
