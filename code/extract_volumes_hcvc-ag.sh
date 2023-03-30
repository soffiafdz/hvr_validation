#!/usr/bin/env bash

### Parse HCVC segmentations (detailed + AMY)
## Extract voxel num for each label and calculate complete HC/CSF
## Labels:
## 111 - L_HC_Tail
## 112 - L_HC_Body
## 113 - L_HC_Head
## 121 - L_CSF_Tail
## 122 - L_CSF_Body
## 123 - L_CSF_Head
## 130 - L_AG
## 211 - R_HC_Tail
## 212 - R_HC_Body
## 213 - R_HC_Head
## 221 - R_CSF_Tail
## 222 - R_CSF_Body
## 223 - R_CSF_Head
## 230 - R_AG

set -ex

IN_DIR=$1
OUT_FILE=$2

[ $# -ne 2 ] || [ ! -d $IN_DIR ] && exit 1

[ -f $OUT_FILE ] && rm $OUT_FILE
printf "ID,L_HC_T,L_HC_B,L_HC_H,L_VC_T,L_VC_B,L_VC_H,L_AMY," > $OUT_FILE
printf "R_HC_T,R_HC_B,R_HC_H,R_VC_T,R_VC_B,R_VC_H,R_AMY\n" >> $OUT_FILE

for img in $IN_DIR/*
do
	bname=$(basename $img .mnc)

	# Left side
	lhct=$(print_all_labels $img | awk '/Label: 111/ {print $3}')
	lhcb=$(print_all_labels $img | awk '/Label: 112/ {print $3}')
	lhch=$(print_all_labels $img | awk '/Label: 113/ {print $3}')
	lvct=$(print_all_labels $img | awk '/Label: 121/ {print $3}')
	lvcb=$(print_all_labels $img | awk '/Label: 122/ {print $3}')
	lvch=$(print_all_labels $img | awk '/Label: 123/ {print $3}')
	lamy=$(print_all_labels $img | awk '/Label: 130/ {print $3}')
	printf "%s,%i,%i,%i,%i,%i,%i,%i," \
		$bname $lhct $lhcb $lhch $lvct $lvcb $lvch $lamy >> $OUT_FILE

	# Right side
	rhct=$(print_all_labels $img | awk '/Label: 211/ {print $3}')
	rhcb=$(print_all_labels $img | awk '/Label: 212/ {print $3}')
	rhch=$(print_all_labels $img | awk '/Label: 213/ {print $3}')
	rvct=$(print_all_labels $img | awk '/Label: 221/ {print $3}')
	rvcb=$(print_all_labels $img | awk '/Label: 222/ {print $3}')
	rvch=$(print_all_labels $img | awk '/Label: 223/ {print $3}')
	ramy=$(print_all_labels $img | awk '/Label: 230/ {print $3}')
	printf "%i,%i,%i,%i,%i,%i,%i\n" \
		$rhct $rhcb $rhch $rvct $rvcb $rvch $ramy >> $OUT_FILE
done
