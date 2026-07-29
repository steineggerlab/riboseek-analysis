#!/usr/bin/mawk -f
BEGIN {OFS="\t"} 
NR==FNR {
	fam[$1]=$2;
	if (length($3) > 0) {
	    clan[$1]=$3;
	} else {
	    clan[$1]="";
	}
	tq[$1]=0;
	cnt[$2]+=1;
	next;} 
	
!($1 in fp) {
	if (fam[$1] == fam[$2]) {
		tp[$1]+=1
	} else if (clan[$1] != "" && clan[$1] == clan[$2]) {
            next
        } else {
	fp[$1]
}} 

END {for (qry in tp) {
	print qry, tp[qry], cnt[fam[qry]], tp[qry]/cnt[fam[qry]]}
} 
