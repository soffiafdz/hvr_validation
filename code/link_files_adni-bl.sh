#! /usr/bin/env bash

## Link adni stx2 files listed in subject list
WORK_PATH=/ipl/ipl27/sfernandez/hvr_validation

# Default list vs input
# ID,SESS
if [ $# -eq 1 ]
then
	LIST=$1
else
	LIST=${WORK_PATH}/lists/adni_baseline.lst
fi

mapfile -t IDS < $LIST
[ ${#IDS[@]} -eq 0 ] && echo "$LIST is empty" && exit 1

ADNI_PATH=/ipl/scratch15/Mahsa/ADNI/LP_2013
DESTINATION_PATH=${WORK_PATH}/data/adni_baselines

for id in ${IDS[@]}
do
	sub=$(echo $id | cut -d, -f1)
	sess=$(echo $id | cut -d, -f2)
	outdir=${DESTINATION_PATH}/$sub
	[[ -d $outdir ]] || mkdir -p $outdir
	ln -sfv \
		${ADNI_PATH}/${sub}/${sess}/stx2/stx2_${sub}_${sess}_mask.mnc $outdir
	ln -sfv \
		${ADNI_PATH}/${sub}/${sess}/stx2/stx2_${sub}_${sess}_t1.mnc $outdir
done
