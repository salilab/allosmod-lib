#!/bin/bash

MINL=10 #min loop length to output adjacent residues (all loops > MINL will be considered)
ADJ=4 #number of residues adjacent to loop that will be outputted
SCL=2.0 #scale interactions by factor of SCL

S1=`grep -n pm.pdb align.ali | awk 'BEGIN{FS=":"}(NR==1){print $1}'`
R_S1=(`grep -n ">P" align.ali | awk 'BEGIN{FS=":"}{print $1}'`)

awk 'BEGIN{FS=""}(NR>=('${S1}'+2)){for(a=1;a<=NF;a++){if($a=="*"){exit};print $a}}' align.ali | tr [a-z.] - >tempglar001 #turn all non amino acids to -
awk 'BEGIN{a=0}{if($1!="-"&&$1!="/"){a+=1; print a}else{print "-"}}' tempglar001 >tempglar000

ctr=1
for i in ${R_S1[@]}; do
    if test `echo "${S1}!=${i}" |bc -l` -eq 1; then
	ctr=$((${ctr} + 1))
	awk 'BEGIN{FS=""}(NR>=('${i}'+2)){for(a=1;a<=NF;a++){if($a=="*"){exit};print $a}}' align.ali >tempglar00${ctr}
    fi
done

echo $ctr | awk '{printf "paste tempglar000 "}{for(a=1;a<=$1;a++){printf "tempglar00"a" "}}' | awk '{print $0}' |sh |\
   awk 'BEGIN{ctrLOOP=0;START=1}
        {for(a=1;a<=NF;a++){if(a==1){II=$a};if(a==2){PM=$a};\
           if(PM!="-"&&PM!="/"){\
             if(a==1){LOOP=1}\
             if(a>2&&$a!="-"&&$a!="/"){LOOP=0};\
             if(a==NF){\
               if(LOOP==1){ctrLOOP+=1};\
               if(ctrLOOP>='${MINL}'&&LOOP==1&&START==0){print II-'${MINL}'-'${ADJ}'+1","II-'${MINL}',"S"; START=1};\
               if(ctrLOOP>='${MINL}'&&LOOP==0){print II","II+'${ADJ}'-1};\
               if(LOOP==0){ctrLOOP=0; START=0};\
             }\
           }\
           if(PM=="/"){ctrLOOP=0;START=1}  }}' |\
   awk 'BEGIN{FS=","}{for(a=$1;a<=$2;a++){print a,"'${SCL}'"}}' >>break.dat


rm tempglar00*
