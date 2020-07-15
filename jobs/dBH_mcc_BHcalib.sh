#!/bin/bash

ml R

cd ../R/
while read LINE;
do 
    ng=$(echo $LINE | cut -d ' ' -f 2)
    nr=$(echo $LINE | cut -d ' ' -f 3)
    pi1=$(echo $LINE | cut -d ' ' -f 4)
    mutype=$(echo $LINE | cut -d ' ' -f 5)
    side=$(echo $LINE | cut -d ' ' -f 6)
    filename="../data/BH_calib_mcc_ng${ng}_nr${nr}_pi1${pi1}_mutype${mutype}_side${side}.RData"
    if [ ! -f $filename ];
    then
	Rscript dBH_mcc_expr_BHcalib.R --ng "${ng}" --nr "${nr}" --pi1 "${pi1}" --mutype "${mutype}" --side "${side}"
    fi
done < "$1"
