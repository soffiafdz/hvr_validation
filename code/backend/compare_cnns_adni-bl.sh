#!/usr/bin/env bash

## Compare overlap similarity between segmentations:
## Simple CNN and Simplified CNN
## XCorrelation

TMPDIR=$(mktemp -d --tmpdir)
trap "rm -rf $TMPDIR" 0 1 2 15

set -ux

printf "id,xcorr,lhc,lcsf,rhc,rcsf\n" > $OUT_FILE

compare_cnns() {
	local BASE_DIR=/ipl/ipl27/sfernandez/hvr_validation
	local OUT_FILE=${BASE_DIR}/data/derivatives/adni-bl_dice_cnn_simple_simplified.csv
	local SIMPLE=${BASE_DIR}/data/derivatives/adni-bl_cnn_hcvc
	local SIMPLIFIED=${BASE_DIR}/data/derivatives/adni-bl_cnn_hcvc_simplified
	local img=$1
	local id=${img%_hcvc-cnn.mnc}
	local img1=${SIMPLE}/$img
	local img2=${SIMPLIFIED}/${id}_hcvc-ag-cnn_simplified.mnc

	# Xcorr
	local xcorr=$(minccmp $img1 $img2 | awk '/xcorr/ {print $2}')

	# Labels
	local kappas=$TMPDIR/${id}_dice.txt
	volume_similarity $img1 $img2 > $kappas
	local lhc=$(awk '/11,Kappa/ {print $2}' $kappas)
	local lcsf=$(awk '/12,Kappa/ {print $2}' $kappas)
	[ -z $lhc ] && lhc=0
	[ -z $rhc ] && rhc=0

	local rhc=$(awk '/21,Kappa/ {print $2}' $kappas)
	local rcsf=$(awk '/22,Kappa/ {print $2}' $kappas)
	[ -z $lcsf ] && lcsf=0
	[ -z $rcsf ] && rcsf=0

	printf "%s,%.3f,%.3f,%.3f,%.3f,%.3f\n" \
		$id $xcorr $lhc $lcsf $rhc $rcsf >> $OUT_FILE
}
export -f compare_cnns

IMGS=($(ls ${BASE_DIR}/data/derivatives/adni-bl_cnn_hcvc))
parallel -j +0 compare_cnns ::: "${IMGS[@]}"
