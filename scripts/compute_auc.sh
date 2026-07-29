#!/bin/bash

input_tsv=$1
target_tsv=$2
lookup=$3
out_roc=$4

script="${BASH_SOURCE[0]}"

gawk 'BEGIN {OFS="\t"} ARGIND==1 {a[$1]=$2;next;} ARGIND==2 {b[$1]=$2;next;} ARGIND==3 {print $0, b[$1], b[$2], a[$2]}' ${target_tsv} ${lookup} ${input_tsv} > ${input_tsv}.m8

bash ${script}/get_tp_cnt.sh ${lookup} ${input_tsv}.m8 ${out_roc}

rm ${input_tsv}.m8
