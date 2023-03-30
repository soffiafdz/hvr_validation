#!/usr/bin/env bash

## Shell script for extracting the ICC and SCALE_factor
## from all preprocessed subjects of ADNI
## and exporting them into a list

set -xu

HERE=/ipl/ipl27/sfernandez/hvr_validation
LIST=${HERE}/lists/adni_hcvc-segm.lst
VOLUMES=${HERE}/data/derivatives/adni-s2_icc_scale.csv

echo "PTID,SCANDATE,ICC,SCALEFACTOR" > $VOLUMES

mapfile -t IDS < $LIST

for id in ${IDS[@]}
do
	segm=$(printf $id | cut -d, -f4)
	if [[ $segm -eq 2 ]]
	then
		sub=$(printf $id | cut -d, -f1)
		visit=$(printf $id | cut -d, -f2)
		visit=${visit//-/}

		stx=$(printf $id | cut -d, -f5)
		mask=${stx/t1/mask}
		xfm=${stx/mnc/xfm}

		# SCALEFACTOR from STX2 xfm
		scale=$(xfm2param $xfm |
			awk '/-scale/{print $2*$3*$4}')

		# ICC (native space)
		icc=$(print_all_labels $mask |
			awk -v scale=$scale '{printf "%.10f", $NF / scale}')

		printf "%s,%s,%f,%f\n" \
			$sub $visit $icc $scale >> $VOLUMES
	fi
done < $LIST
