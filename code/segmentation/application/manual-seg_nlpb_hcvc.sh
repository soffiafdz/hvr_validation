#! /bin/bash

## RUN SNIPE using MCCV64 with Dorothee's original labels

# Directories
HERE=/ipl/ipl27/sfernandez/snipe_hvr
SCRATCH=${HERE}/proc_4

for img in ${HERE}/original_labels/t1w/*
do
	id=$(basename $img _m00_t1w.mnc)
	id=${id##*_}
	outdir=${SCRATCH}/${id}

	[[ -d $outdir ]] || mkdir $outdir

	${HERE}/snipe_minipipe_mccv_reduc.pl \
		$img $outdir \
		--model-dir $HERE \
		--subject $id \
		--nuc \
		--3t \
		--mccv 64
done
