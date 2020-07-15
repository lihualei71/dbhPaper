#!/bin/bash

ml R

cd ../R/
while read LINE;
do 
    n=$(echo $LINE | cut -d ' ' -f 2)
    pi1=$(echo $LINE | cut -d ' ' -f 3)
    mutype=$(echo $LINE | cut -d ' ' -f 4)
    side=$(echo $LINE | cut -d ' ' -f 5)
    filename="../data/BH_calib_mvgauss_n${n}_pi1${pi1}_mutype${mutype}_side${side}.RData"
    if [ ! -f $filename ];
    then
	Rscript dBH_mvgauss_expr_BHcalib.R --n "${n}" --pi1 "${pi1}" --mutype "${mutype}" --side "${side}"
    fi
done < "$1"
