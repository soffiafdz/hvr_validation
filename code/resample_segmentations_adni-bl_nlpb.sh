#!/usr/bin/env bash

set -xu

HERE=/ipl/ipl27/sfernandez/hvr_validation
OUTDIRS=${HERE}/tmp/nlpb_hcvc

for dir in $OUTDIRS/*
do
	id=$(basename $dir)

	labels=$(ls $dir/stx_${id}_*_HVR.mnc)
	cleaned=${labels/HVR/hcvc}
	minclookup \
		-discrete \
		-float \
		-clobber \
		-lut "11 11; 12 12; 21 21; 22 22" \
		$labels $cleaned


	xfm=$(ls $dir/stx_${id}_*_bbox.xfm)
	xfminv=${xfm/bbox/bbox_inv}
	xfminvert -clobber $xfm $xfminv


	resampled=${cleaned/.mnc/_resampled.mnc}
	itk_resample \
		--labels \
		--transform $xfminv \
		--clobber \
		$cleaned $resampled


	bname=$(basename $resampled _resampled.mnc)
	outfile=$HERE/data/derivatives/adni-bl_nlpb_hcvc/${bname/stx/adni}-nlpb.mnc
	ln -sfv $resampled $outfile
done

