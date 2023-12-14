#! /bin/bash

## Run MALF single ADNI subject from list
## Made for GNU Parallel

BASE_DIR=/ipl/ipl27/sfernandez/hvr_validation

# Read subject list into an Array
mapfile -t IDS < ${BASE_DIR}/lists/adni_baseline.lst
#mapfile -t IDS < ${BASE_DIR}/lists/adni_baseline_rev.lst
#IDS=( $(shuf -e "${IDS[@]}") )

# Main
id="${IDS[$1]}"
sub=$(printf $id | cut -d, -f1)
session=$(printf $id | cut -d, -f2)

source ${BASE_DIR}/code/run_segmentations/adni-bl_malf_hcvc.sh \
	$sub \
	$session
