#! /usr/bin/env sh

HERE=/ipl/ipl27/sfernandez/hvr_validation
orig_list=${HERE}/lists/adni_baseline.lst
new_list=${HERE}/lists/adni-bl_scanner.lst

[[ -e $new_list ]] && rm $new_list

while read -r line
do
	sub=$(printf $line | cut -d, -f1)
	sess=$(printf $line | cut -d, -f2)
	img=/ipl/scratch15/Mahsa/ADNI/LP_2013/$sub/$sess/clp/clp_${sub}_${sess}_t1.mnc
	tesla=$(mincheader $img | awk '/field_value/ {print $3}')
	if [[ -z $tesla ]]
	then
		printf "%s: %s\n" $sub $sess
	else
		printf "%s,%s,%s\n" $sub $sess $tesla >> $new_list
	fi
	tesla=
done < $orig_list

## Add manually:
# Philips 3T
#012_S_5121,20130315,3.
