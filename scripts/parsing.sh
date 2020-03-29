#!/bin/bash

natoms=108
POS_DAT=./results/positions.dat
DATA_ANLYSIS_DIR=./data_analysis

for Ar in `seq -f '%03g' 1 ${natoms}`; do
    pattern=Ar${Ar}
    awk -v pat="$pattern" -F ":" '$0~pat {print $0}' ${POS_DAT} \
    | sed 's/   /*/g' | sed 's/  /*/g' | sed 's/ /*/g' | sed 's/*/ /g' \
    | awk ' {print $2, $3, $4}' > ${DATA_ANLYSIS_DIR}/Ar${Ar}.dat
done

