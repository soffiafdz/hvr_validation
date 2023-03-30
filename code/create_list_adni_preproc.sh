#!/usr/bin/env bash

## Shell script for extracting all sessions of all subjects of ADNI
## and exporting them into a list
## Include MAG_strength

set -ux

ADNI=/data/ipl/scratch15/Mahsa/ADNI/LP_2013
HERE=/ipl/ipl27/sfernandez/hvr_validation
LIST=${HERE}/lists/adni_preproc.lst

[[ -f $LIST ]] && rm $LIST

for sub_path in $ADNI/[0-9]*
do
	sub=$(basename $sub_path)
	for sess_path in $sub_path/20*
	do
		session=$(basename $sess_path)

		# TESLA
		clp=${ADNI}/${sub}/${session}/clp/clp_${sub}_${session}_t1.mnc
		tesla=$(mincheader $clp | awk '/field_value/ {print $3}')

		# STX img
		stx2=${ADNI}/${sub}/${session}/stx2/stx2_${sub}_${session}_t1.mnc
		[[ -e $stx2 ]] \
			&& printf "%s,%s,%.1f,%s\n" $sub $session $tesla $stx2 >> $LIST
	done
done
