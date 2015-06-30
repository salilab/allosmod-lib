#!/bin/bash

MODEL_FILE=$1
R_START=(`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $0}' ${MODEL_FILE} | awk 'BEGIN{FS="";lc=""}(lc!=$22){li=0}\
          (li!=$23$24$25$26&&$1$2$3=="HET"){print NR-1}{li=$23$24$25$26;lc=$22}'`)
R_INDEX=(`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $0}' ${MODEL_FILE} | awk 'BEGIN{FS="";lc=""}(lc!=$22){li=0}\
          (li!=$23$24$25$26&&$1$2$3=="HET"){print $23$24$25$26}{li=$23$24$25$26;lc=$22}'`)
R_GTYPE=(`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $0}' ${MODEL_FILE} | awk 'BEGIN{FS="";lc=""}(lc!=$22){li=0}\
          (li!=$23$24$25$26&&$1$2$3=="HET"){print $18$19$20}{li=$23$24$25$26;lc=$22}'`)
R_CHAIN=(`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $0}' ${MODEL_FILE} | awk 'BEGIN{FS="";lc=""}(lc!=$22){li=0}\
          (li!=$23$24$25$26&&$1$2$3=="HET"){print $22}{li=$23$24$25$26;lc=$22}'`)

R_CONECT_IND=(`awk '{print $1}' get_rest.in`)
R_CONECT_ATM=(`awk '{print $2}' get_rest.in`)

if test ${#R_CONECT_IND[@]} != ${#R_START[@]}; then
    echo "Number of sugars in get_rest.in: "${#R_CONECT_IND[@]}
    echo "Number of sugars in "$MODEL_FILE": "${#R_START[@]}
    echo Connections for all sugars are not defined
    exit
fi

ctr=-1
for i in ${R_START[@]}; do
    ctr=$((${ctr}+1))
    index=${R_INDEX[$ctr]}
    gtype=${R_GTYPE[$ctr]}
    chain=${R_CHAIN[$ctr]}
    con_ind=${R_CONECT_IND[$ctr]}
    con_atm=${R_CONECT_ATM[$ctr]}

    if test $gtype == "NAN"; then
    R_SHIFT=(`awk 'BEGIN{FS=""}{print $14$15$16,$22,$23$24$25$26}' $MODEL_FILE | awk '($3=="'${index}'"&&$2=="'${chain}'"){print $1}' |\
       awk '($1=="C1"){print "1 "NR}($1=="C2"){print "2 "NR}($1=="C3"){print "3 "NR}($1=="C4"){print "4 "NR}($1=="C5"){print "5 "NR}($1=="O6"){print "6 "NR}\
            ($1=="C7"){print "7 "NR}($1=="H32"){print "8 "NR}($1=="O4"){print "9 "NR}($1=="C6"){print "10 "NR}($1=="N"){print "11 "NR}\
            ($1=="H6"){print "12 "NR}($1=="H31"){print "13 "NR}($1=="H4"){print "14 "NR}($1=="H5"){print "15 "NR}' |\
       sort -nk1 | awk '{print $2+'${i}'}'`)
    else
    R_SHIFT=(`awk 'BEGIN{FS=""}{print $14$15$16,$22,$23$24$25$26}' $MODEL_FILE | awk '($3=="'${index}'"&&$2=="'${chain}'"){print $1}' |\
       awk '($1=="C1"){print "1 "NR}($1=="C2"){print "2 "NR}($1=="C3"){print "3 "NR}($1=="C4"){print "4 "NR}($1=="C5"){print "5 "NR}($1=="O5"){print "6 "NR}\
            ($1=="O2"||$1=="N"){print "7 "NR}($1=="O3"){print "8 "NR}($1=="O4"){print "9 "NR}($1=="C6"){print "10 "NR}($1=="H1"){print "11 "NR}\
            ($1=="H2"){print "12 "NR}($1=="H3"){print "13 "NR}($1=="H4"){print "14 "NR}($1=="H5"){print "15 "NR}' |\
       sort -nk1 | awk '{print $2+'${i}'}'`)
    fi

    connect=`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $0}' ${MODEL_FILE} | awk 'BEGIN{FS=""}{print $14$15$16" "$23$24$25$26}' |\
	     awk '($1=="'${con_atm}'"&&$2=='${con_ind}'){print NR}'`

    cat ${gtype}.rsr | sed "s/AAAAA/${R_SHIFT[0]}/g" | sed "s/BBBBB/${R_SHIFT[1]}/g" | sed "s/CCCCC/${R_SHIFT[2]}/g" |\
                       sed "s/DDDDD/${R_SHIFT[3]}/g" | sed "s/EEEEE/${R_SHIFT[4]}/g" | sed "s/FFFFF/${R_SHIFT[5]}/g" |\
                       sed "s/GGGGG/${R_SHIFT[6]}/g" | sed "s/HHHHH/${R_SHIFT[7]}/g" | sed "s/IIIII/${R_SHIFT[8]}/g" |\
                       sed "s/JJJJJ/${R_SHIFT[9]}/g" | sed "s/KKKKK/${R_SHIFT[10]}/g"| sed "s/LLLLL/${R_SHIFT[11]}/g" |\
                       sed "s/MMMMM/${R_SHIFT[12]}/g"| sed "s/NNNNN/${R_SHIFT[13]}/g"| sed "s/OOOOO/${R_SHIFT[14]}/g" |\
                       sed "s/XXXXX/${connect}/g" |\
                       awk '{if(NF==12){printf "R    3   1   1  27   2   2   1%6.0f%6.0f   %10.4f%10.4f\n",$9,$10,$11,$12}\
                             else{printf "R    3   1   1  27   2   2   1%6.0f%6.0f%6.0f%6.0f   %10.4f%10.4f\n",$9,$10,$11,$12,$13,$14}}'

done
