#! /bin/bash

## Run MALF with MCCV ##

# Directories
LIB_DIR=/ipl/scratch22/sfernandez/HVR_MALF/lib_dorothee2
TMP_DIR=/ipl/scratch22/sfernandez/HVR_MALF/tmp/mccv8
#QC_DIR=/ipl/scratch22/sfernandez/HVR_MALF/qc_mccv6
OUT_DIR=/ipl/scratch22/sfernandez/HVR_MALF/proc/mccv8/original

# Read subject list into an Array
mapfile -t IDS < subjects.lst
toRun=()

# If argument(s) given check if in IDS, else run IDS
if [ $# -gt 0 ]
then
	for arg in "$@"
	do
		[[ "${IDS[@]}" =~ "$arg" ]] && toRun+=( "$arg" )
	done
else
	toRun=( "${IDS[@]}" )
fi

# Main
for id in "${toRun[@]}"
do
	# Input
	in_mri=${LIB_DIR}/orig_labels/mri/${id}_t1.mnc
	in_mask=${LIB_DIR}/orig_labels/mri/${id}_mask.mnc
	in_xfm=${LIB_DIR}/orig_labels/identity.xfm
	bname=malf_mccv_${id}
	out=${OUT_DIR}/${bname}.mnc

	sublibDir=mccv_submodels2/$id
	#submodelLABS=${sublibDir}/labels
	submodelLABS=${sublibDir}/labels_reduced
	submodelMOD=${sublibDir}/models
	submodelMRI=${submodelMOD}_mri
	submodelSEG=${submodelMOD}_seg
	submodelXFM=${submodelMOD}_xfms

	# Monte-Carlo Sub-libraries with 20 random subjects' MRI/Labels

	MCCV=64

	tmpSubdir=$TMP_DIR/$id
	[[ ! -d "$tmpSubdir" ]] && mkdir "$tmpSubdir"

	if [[ ! -d "$sublibDir" ]]
	then
		echo "Creating skeleton sublibrary"
		date
		idsUnbiased=( ${IDS[@]/$id/} ) # IDS without input ID
		libIDs=()
		i=0
		while [ ${#libIDs[@]} -lt $MCCV ]
		do
			random_id=${idsUnbiased[$RANDOM % ${#idsUnbiased[@]}]}
			if [[ ! ${libIDs[@]} =~ "$random_id" ]]
			then
				(( i++ ))
				echo "Random subject #${i}: $random_id"
				libIDs+=( $random_id )
			fi
		done

		# Create submodel directory for specific image
		mkdir -p "$sublibDir"
		ln -s ${LIB_DIR}/models $sublibDir
		mkdir -p ${sublibDir}/{labels{,_reduced},models_{mri,seg,xfms}}

		# Fill it
		for libID in ${libIDs[@]}
		do
			# Label
			ln -s ${LIB_DIR}/orig_labels/labels/*${libID}*.mnc $submodelLABS
			# MRI
			ln -s ${LIB_DIR}/models_mri/${libID}* $submodelMRI
			# SEG
			ln -s ${LIB_DIR}/models_seg/${libID}* $submodelSEG
			# XFM
			ln -s ${LIB_DIR}/models_xfms/${libID}* $submodelXFM
		done
	fi

	## Run on BIC
	#qsub -q ipl.q -V -cwd -j y -N MALF_$id -o proc/$bname.log <<END
#! /bin/sh
# MALF
	echo "MALF $id"
	date
	./bin/seg_hippo.pl \
		-tmpdir $tmpSubdir -keeptmp \
		-modelDIR $submodelMOD \
		-modelSEG $submodelSEG \
		-modelMRI $submodelMRI \
		-modelLAB $submodelLABS \
		-modelXFM $submodelXFM \
		$in_mri $in_mask $id $in_xfm \
		$bname $out

	# QC images
	#echo "Making QC images for HC"
	#date
	#qcImg=${QC_DIR}/${bname}_mccv64.jpg
	#mkdir -p $QC_DIR
	#if [ -e  $qcImg ]
	#then
		#echo "Image already exists. Overwriting it."
		#rm $qcImg
	#fi
	#./bin/make_hc_verif_jpg $in_mri $out -title $bname $qcImg
	#date
#END
done
