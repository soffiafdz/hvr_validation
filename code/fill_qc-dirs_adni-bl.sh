#!/usr/bin/env bash

## Create directories of qc images to be loaded on QRATER
## Fill them with subjects from lists

set -xu
BASE_DIR=/ipl/ipl27/sfernandez/hvr_validation

# Read PTIDS lists into Arrays
mapfile -t REGIS_FAILS < ${BASE_DIR}/lists/adni-bl_qc_lin-reg_ids.lst
mapfile -t PASSES < ${BASE_DIR}/lists/adni-bl_to-qc_ids.lst

for method in malf nlpb cnn/hcvc cnn/hcvc-ag
do
	qc_dir=${BASE_DIR}/plots/qc_adni-bl/${method}

	# Linear registration
	reg_dir=${qc_dir}/reg_fails
	[[ -d $reg_dir ]] || mkdir -p $reg_dir
	for id in "${REGIS_FAILS[@]}"
	do
		ln -f ${qc_dir}/all/adni_${id}*.jpg ${reg_dir}
	done

	# Curated
	out_dir=${qc_dir}/curated
	[[ -d $out_dir ]] || mkdir -p $out_dir
	for id in "${PASSES[@]}"
	do
		ln -f ${qc_dir}/all/adni_${id}*.jpg ${out_dir}
	done
done
