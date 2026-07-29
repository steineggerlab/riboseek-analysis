#!/bin/bash

lookup=$1
aln=$2
output=$3

script="${BASH_SOURCE[0]}"

${script}/roc1.awk $lookup $aln | sort -k4,4rn > $output
awk '{sum+=$4} END {print sum/NR, NR}' $output
