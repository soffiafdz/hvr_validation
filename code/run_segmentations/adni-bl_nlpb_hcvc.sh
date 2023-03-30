#! /bin/bash

## Run SNIPE ADNI subjects from list

# Directories
BASE_DIR=/ipl/ipl27/sfernandez/hvr_validation
LIB_DIR=${BASE_DIR}/libraries/snipe_adni
TMP_DIR=${BASE_DIR}/tmp/nlpb_hcvc
QC_DIR=${BASE_DIR}/plots/qc_adni-bl/nlpb/all
OUT_DIR=${BASE_DIR}/data/derivatives/adni-bl_nlpb_hcvc

# Read subject list into an Array
mapfile -t IDS < ${BASE_DIR}/lists/adni_bl_scanner.lst

# Main
for id in "${IDS[@]}"
do
	sub=$(printf $id | cut -d, -f1)
	session=$(printf $id | cut -d, -f2)
	tesla=$(printf $id | cut -d, -f3)

	[[ $tesla == 3. ]] \
		&& argument_3_tesla="--3t" \
		|| argument_3_tesla=

	# Input
	in_mri=${BASE_DIR}/data/adni_baselines/${sub}/stx2_${sub}_${session}_t1.mnc
	bname=adni_${sub}_${session}_hcvc-nlpb

	# SNIPE
	echo "SNIPE $sub"

	# Scratch
	tmpSubdir=${TMP_DIR}/${sub}
	[[ ! -d "$tmpSubdir" ]] && mkdir "$tmpSubdir"

	output=${tmpSubdir}/stx_${sub}_${session}_hcvc.mnc
	changed=
	if [[ ! -e $output ]]
	then
		date
		${BASE_DIR}/code/nlpb/snipe_minipipe_reduc.pl \
			$in_mri $tmpSubdir \
			--model-dir $LIB_DIR \
			--subject ${sub}_${session} \
			--nuc \
			$argument_3_tesla
		changed=Y
	else
		echo "$sub is already processed"
	fi

	outfile=${OUT_DIR}/${bname}.mnc
	if [[ ! -e $outfile || -n $changed ]]
	then
		echo "Saving output of $sub"
		minclookup \
			-clobber \
			-discrete \
			-int \
			-lut '11 11; 12 12; 21 21; 22 22' \
			$output $outfile
		changed=Y
	fi

	# QC images
	qcImg=${QC_DIR}/${bname}.jpg
	if [[ ! -e $qcImg || -n $changed ]]
	then
		echo "Making QC images for HC"
		date
		mkdir -p $QC_DIR
		if [ -e  $qcImg ]
		then
			echo "Image already exists. Overwriting it."
			rm $qcImg
		fi
		${BASE_DIR}/code/qc/qc_plot_reduc.pl $in_mri $outfile -title $bname $qcImg
	fi
done
