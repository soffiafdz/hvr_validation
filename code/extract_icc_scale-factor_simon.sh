#!/usr/bin/env bash

## Shell script for extracting the first session of all subjects of ADNI
## and exporting them into a list

set -xu

HERE=/ipl/ipl27/sfernandez/hvr_validation
SIMON=/ipl/scratch20/Mahsa/Simon/LP/SIMON
VOLUMES=${HERE}/data/derivatives/simon_icc_scale.csv

echo "SCAN,ICC,SCALEFACTOR" > $VOLUMES

for dir in ${SIMON}/ses-*
do
	scan=$(basename $dir)
	mask=${dir}/stx2/stx2_SIMON_${scan}_mask.mnc
	xfm=${mask/mask.mnc/t1.xfm}

	# SCALEFACTOR from STX2 xfm
	scale=$(xfm2param $xfm |
		awk '/-scale/{print $2*$3*$4}')

	# ICC (native space)
	icc=$(print_all_labels $mask |
		awk -v scale=$scale '{printf "%.10f", $NF / scale}')

	# Print to file
	printf "%s,%f,%f\n" \
		$scan $icc $scale >> $VOLUMES
done
