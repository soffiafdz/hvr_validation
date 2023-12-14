#!/usr/bin/env bash

## Shell script for extracting the first session of all subjects of ADNI
## and exporting them into a list

set -ux

ADNI=/data/ipl/scratch15/Mahsa/ADNI/LP_2013
BASE_DIR=/ipl/ipl27/sfernandez/hvr_validation
LIST=${BASE_DIR}/lists/adni_baseline.lst

[[ -f $LIST ]] && rm $LIST

for full_path in $ADNI/[0-9]*; do
	sub=$(basename $full_path)
	session=$(ls --group-directories-first $full_path | head -n 1)
	[[ -e ${ADNI}/$sub/$session/stx2/stx2_${sub}_${session}_t1.mnc ]] \
		&& echo "$sub,$session" >> $LIST
done
