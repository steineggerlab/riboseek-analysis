#!/usr/bin/awk

in_file=$1
out_tsv=$2

awk 'BEGIN {OFS="\t"} (substr($1, 1, 1) != "#" && !($3$1 in a)) {print $3, $1, $(NF-3);a[$3$1]}' ${in_file} > ${out_tsv}
