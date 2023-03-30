#! /usr/bin/env bash

## Relabel CNN HVR-AG labels to simplified HC & VC

set -xue

HERE=/ipl/ipl27/sfernandez/hvr_validation
OUTDIR=${HERE}/data/derivatives/adni-bl_cnn_hcvc_simplified
[[ -d $OUTDIR ]] || mkdir $OUTDIR

IN=$1
OUT=${OUTDIR}/$(basename $IN .mnc)_simplified.mnc

itk_resample $IN $OUT \
	--clobber \
	--labels \
	--byte \
	--lut-string \
	'111 11; 112 11; 113 11; 121 12; 122 12; 123 12; 211 21; 212 21; 213 21; 221 22; 222 22; 223 22'
