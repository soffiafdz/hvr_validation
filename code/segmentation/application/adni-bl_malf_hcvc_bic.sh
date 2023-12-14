#! /bin/bash

## Run MALF ADNI subjects from list

# Directories
BASE_DIR=/ipl/ipl27/sfernandez/hvr_validation
LIB_DIR=${BASE_DIR}/libraries/malf_dorothee
TMP_DIR=${BASE_DIR}/tmp/malf_hcvc
QC_DIR=${BASE_DIR}/plots/qc_adni-bl/malf/all
OUT_DIR=${BASE_DIR}/data/derivatives/adni-bl_malf_hcvc

# Read subject list into an Array
mapfile -t IDS < ${BASE_DIR}/lists/adni_baseline.lst

# Main
for id in "${IDS[@]}"
do
	sub=$(printf $id | cut -d, -f1)
	session=$(printf $id | cut -d, -f2)

	# Input
	in_mri=${BASE_DIR}/data/adni_baselines/${sub}/stx2_${sub}_${session}_t1.mnc
	in_mask=${BASE_DIR}/data/adni_baselines/${sub}/stx2_${sub}_${session}_mask.mnc
	in_xfm=${BASE_DIR}/data/identity.xfm
	bname=adni_${sub}_${session}_hcvc-malf
	out=${OUT_DIR}/${bname}.mnc

	## Run on BIC
	if [[ ! -e $out ]]
	then
		echo "Submitting job for $out"
		qsub -q ipl.q -V -cwd -j y -N MALF_$sub -o proc/$bname.log <<END
#! /bin/sh
#MALF
	echo "MALF $sub"

	# Scratch
	tmpSubdir=${TMP_DIR}/${sub}
	[[ ! -d "$tmpSubdir" ]] && mkdir "$tmpSubdir"

	date
	${BASE_DIR}/code/malf/seg_hippo.pl \
		-tmpdir $tmpSubdir \
		$in_mri $in_mask $sub $in_xfm \
		$bname $out

	 QC images
	echo "Making QC images for HC"
	date
	qcImg=${QC_DIR}/${bname}.jpg
	mkdir -p $QC_DIR
	if [ -e  $qcImg ]
	then
		echo "Image already exists. Overwriting it."
		rm $qcImg
	fi
	${BASE_DIR}/code/qc/qc_plot_reduc.pl $in_mri $out -title $bname $qcImg
END
	fi
done
