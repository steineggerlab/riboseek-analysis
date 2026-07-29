#!/usr/bin/awk

in_file=$1
out_tsv=$2

awk 'BEGIN {OFS="\t"} (substr($1, 1, 1) != "#" && !($1$2 in a)) {print $1, $2, $(NF-1);a[$1$2]}' ${in_file} > ${out_tsv}
