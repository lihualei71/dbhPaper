#!/bin/bash

ml R

cd ../R/
while read LINE;
do 
    n=$(echo $LINE | cut -d ' ' -f 2)
    df=$(echo $LINE | cut -d ' ' -f 3)
    pi1=$(echo $LINE | cut -d ' ' -f 4)
    mutype=$(echo $LINE | cut -d ' ' -f 5)
    side=$(echo $LINE | cut -d ' ' -f 6)
    filename="../data/BH_calib_mvt_n${n}_df${df}_pi1${pi1}_mutype${mutype}_side${side}.RData"
    if [ ! -f $filename ];
    then
	Rscript dBH_mvt_expr_BHcalib.R --n "${n}" --df "${df}" --pi1 "${pi1}" --mutype "${mutype}" --side "${side}"
    fi
done < "$1"
