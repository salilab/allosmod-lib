#!/bin/bash
#inserts gap in alignment file from LRES to URES
#IGAP is the ith sequence in alignment file used to define the residue index
#insert 4 residue gap for sugars

FIL=$1 #alignment file
IGAP=$2
LRES=$3
URES=$4

iLRES=`awk 'BEGIN{FS="";iseq=0}{head=0}($1==">"){head=1;iseq+=1;ires=0;ctr=0}($1$2$3$4$5=="struc"||$1$2$3$4$5=="seque"){head=1}(NF==0){head=1}\
(head==0){for(a=1;a<=NF;a++){if($a!="-"&&$a!="/"){ires+=1}{ctr+=1};\
if(ires=='${LRES}'&&iseq=='${IGAP}'){print ctr;exit}}}' $FIL`

awk 'BEGIN{FS="";iseq=0}{head=0}($1==">"){print $0;head=1;iseq+=1;ires=0;ctr=0}($1$2$3$4$5=="struc"||$1$2$3$4$5=="seque"){print $0;head=1}(NF==0){print $0;head=1}\
(head==0){for(a=1;a<=NF;a++){if($a!="-"&&$a!="/"){ires+=1}{ctr+=1};\
if(ctr=='${iLRES}'&&iseq!='${IGAP}'){for(b='${LRES}';b<='${URES}';b++){printf "-"}}\
{printf $a}\
if(ires=='${URES}'&&iseq=='${IGAP}'){for(b='${LRES}';b<='${URES}';b++){printf "-"}}\
if(a==NF){printf "\n"}}}' $FIL

