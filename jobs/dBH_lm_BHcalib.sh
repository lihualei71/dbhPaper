#!/bin/bash

ml R

cd ../R/
while read LINE;
do 
    Xseed=$(echo $LINE | cut -d ' ' -f 2)
    n=$(echo $LINE | cut -d ' ' -f 3)
    p=$(echo $LINE | cut -d ' ' -f 4)
    pi1=$(echo $LINE | cut -d ' ' -f 5)
    mutype=$(echo $LINE | cut -d ' ' -f 6)
    side=$(echo $LINE | cut -d ' ' -f 7)
    filename="../data/BH_calib_lm_n${n}_p${p}_pi1${pi1}_mutype${mutype}_side${side}_Xseed${Xseed}.RData"
    if [ ! -f $filename ];
    then
	Rscript dBH_lm_expr_BHcalib.R --n "${n}" --p "${p}" --pi1 "${pi1}" --mutype "${mutype}" --side "${side}" --Xseed "${Xseed}"
    fi
done < "$1"
