#!/bin/bash
#averages two files... must have exact same #atoms

# Absolute path containing this and other scripts
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

FIL1=$1 #if hetatm's, then taken from this file first (not randomized)
FIL2=$2
ID1=$3 #input structures
ID2=$4
AFIL=align.ali

NATOM1=`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"){print $0}' $FIL1 | awk 'END{print NR}'`
NATOM2=`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"){print $0}' $FIL2 | awk 'END{print NR}'`
if test $NATOM1 != $NATOM2; then
  echo "ERROR pm files not same length" >&2
  exit 1
fi

allosmod setchain $FIL1 A > temp9941
allosmod setchain $FIL2 A > temp9942

echo >>run.log
#structurally align 2 allosteric states
echo "getavgpdb: Salign 2 allosteric states (only consider similar residues) to generate random interpolation structure" >>run.log
NRES=`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"){print $0}' temp9941 |\
      awk 'BEGIN{FS="";lc=""}(lc!=$22){li=0}(li!=$23$24$25$26){a+=1}{li=$23$24$25$26;lc=$22}END{print a}'`
rcut=10.0
allosmod get_qiavg_ca temp9941 $rcut temp9942 >>run.log
QI=(`awk '{printf $2" "}' qi_1.dat`)
if test -e tempfit.pdb; then rm tempfit.pdb; fi
for i in `seq ${NRES}`; do
    if test `echo "${QI[$i-1]}>0.8" |bc -l` -eq 1; then
	$SCRIPT_DIR/getsegpdb temp9942 A $i $i >>tempfit.pdb
    fi
done

$SCRIPT_DIR/salign_sub.sh temp9941 tempfit.pdb temp9942

echo >>run.log

#if salign fails, set avgpdb to either FIL1 or FIL2, randomly 
if (test ! -e 1temp88_fit.pdb -o ! -e tempfit_fit.pdb); then 
    date --rfc-3339=ns | awk -F. '{print $2}' | awk -F- '{print $1}' | awk 'BEGIN{FS=""}{print $NF%2}' |\
         awk '{if($1==1){print "cp '${FIL1}' avgpdb.pdb"}else{print "cp '${FIL2}' avgpdb.pdb"}}' |sh
    echo "Salign of two input structures failed, using only one to define starting structure" >>run.log
    echo "Salign of two input structures failed, using only one to define starting structure"
    exit
fi

awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"){print $0}' 1temp88_fit.pdb >temp9941
awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"){print $0}' tempfit_fit.pdb >temp9942

#make random numbers between 0 and 1
if test -e tempN; then rm tempN; fi
NN=`awk 'END{print NR}' temp9941`
for s in `seq $NN`; do
    date --rfc-3339=ns | awk -F. '{print $2}' | awk -F- '{print $1/1000000000}' >>tempN
done
NMIN=`$SCRIPT_DIR/getmin tempN 1`
NMAX=`$SCRIPT_DIR/getmax tempN 1`
perl $SCRIPT_DIR/randomize_list_2.pl tempN | awk '{a=($1-"'${NMIN}'")/("'${NMAX}'"-"'${NMIN}'")}{print a,1-a}' >tempNN

#interpolate between two structures
awk 'BEGIN{FS=""}{for(a=31;a<=38;a++){printf $a}}{printf " "}\
{for(a=39;a<=46;a++){printf $a}}{printf " "}\
{for(a=47;a<=54;a++){printf $a}}{printf "\n"}' temp9941 >temp9991
awk 'BEGIN{FS=""}{for(a=31;a<=38;a++){printf $a}}{printf " "}\
{for(a=39;a<=46;a++){printf $a}}{printf " "}\
{for(a=47;a<=54;a++){printf $a}}{printf " "$18$19$20}{printf "\n"}' temp9942 >temp9992
#paste temp9991 temp9992 | awk '{if($7=="PHE"||$7=="HIS"||$7=="TYR"||$7=="TRP"||$7=="HID"||$7=="HIE"||$7=="HIP"||$7=="HSD"||$7=="HSE"||$7=="HSP"){print $1,$1,$2,$2,$3,$3}\
#                               else{print $4,$1,$5,$2,$6,$3}}' >tempdd
#paste tempdd tempNN | awk '{printf "%8.3f%8.3f%8.3f\n",$7*$1+$8*$2,$7*$3+$8*$4,$7*$5+$8*$6}' >tempd
paste temp9991 temp9992 | awk '{print $4,$1,$5,$2,$6,$3,$7}' >tempdd
paste tempdd tempNN | awk 'BEGIN{lastcyc=0;last8=0}{if($7=="PHE"||$7=="HIS"||$7=="TYR"||$7=="TRP"||$7=="HID"||$7=="HIE"||$7=="HIP"||$7=="HSD"||$7=="HSE"||$7=="HSP"){\
                              if((lastcyc==0&&$8<0.5)||(lastcyc==1&&last8<0.5)){printf "%8.3f%8.3f%8.3f\n",$1,$3,$5}\
                              else{printf "%8.3f%8.3f%8.3f\n",$2,$4,$6}if(lastcyc==0){last8=$8}{lastcyc=1}}\
                            else{printf "%8.3f%8.3f%8.3f\n",$8*$1+$9*$2,$8*$3+$9*$4,$8*$5+$9*$6; lastcyc=0}}' >tempd

#reassemble PDB
awk 'BEGIN{FS=""}{for (a=1; a<=30 ;a++) {printf $a}}{printf "\n"}' temp9941 >tempa
awk 'BEGIN{FS=""}{for (a=55; a<=76 ;a++) {printf $a}}{printf "\n"}' temp9941 >tempe
paste -d"\0" tempa tempd tempe >avgpdb.pdb

#get hetatm's for avgpdb.pdb
NSTRT=`grep -n pm.pdb ${AFIL} | awk 'BEGIN{FS=":"}{print $1;exit}'`
R_iHET0=(`awk 'BEGIN{FS="";ires=0}{for(a=1;a<=NF;a++){if(NR>'${NSTRT}'+1&&$a!="-"){ires+=1}if(NR>'${NSTRT}'+1&&($a=="."||$a=="h")){print ires}\
          if(NR>'${NSTRT}'+1&&$a=="*"){exit}}}' ${AFIL}`)
NSTRT=`grep -n ${ID1} ${AFIL} | awk 'BEGIN{FS=":"}{print $1;exit}'`
R_iHET1=(`awk 'BEGIN{FS="";ires=0}{for(a=1;a<=NF;a++){if(NR>'${NSTRT}'+1&&$a!="-"){ires+=1}if(NR>'${NSTRT}'+1&&($a=="."||$a=="h")){print ires}\
          if(NR>'${NSTRT}'+1&&$a=="*"){exit}}}' ${AFIL}`)
NSTRT=`grep -n ${ID2} ${AFIL} | awk 'BEGIN{FS=":"}{print $1;exit}'`
R_iHET2=(`awk 'BEGIN{FS="";ires=0}{for(a=1;a<=NF;a++){if(NR>'${NSTRT}'+1&&$a!="-"){ires+=1}if(NR>'${NSTRT}'+1&&($a=="."||$a=="h")){print ires}\
          if(NR>'${NSTRT}'+1&&$a=="*"){exit}}}' ${AFIL}`)

for het0 in ${R_iHET0[@]}; do
    match=0
    for het1 in ${R_iHET1[@]}; do
	if test `echo "${het0}==${het1}" |bc -l` -eq 1; then
            awk 'BEGIN{FS=""}{print $23$24$25$26$27" "$0}' $FIL1 |\
            awk '($1=='${het1}'){print $0}' |\
            awk 'BEGIN{FS=""}{for(a=7;a<=NF;a++){printf $a;if(a==NF){printf "\n"}}}' >>avgpdb.pdb
	    match=1
	fi
    done
    if test `echo "${match}==0" |bc -l` -eq 1; then
	for het2 in ${R_iHET2[@]}; do
	    if test `echo "${het0}==${het2}" |bc -l` -eq 1; then
                awk 'BEGIN{FS=""}{print $23$24$25$26$27" "$0}' $FIL2 |\
                awk '($1=='${het2}'){print $0}' |\
                awk 'BEGIN{FS=""}{for(a=7;a<=NF;a++){printf $a;if(a==NF){printf "\n"}}}' >>avgpdb.pdb
		match=1
	    fi
	done
    fi
    if test `echo "${match}==0" |bc -l` -eq 1; then
      echo "ERROR getavgpdb: can't locate hetatm" >&2
      exit 1
    fi
done

#put missing hetatm's back into pm files
for het0 in ${R_iHET0[@]}; do
    match1=0
    for het1 in ${R_iHET1[@]}; do
	if test `echo "${het0}==${het1}" |bc -l` -eq 1; then
	    match1=1
	fi
    done
    match2=0
    for het2 in ${R_iHET2[@]}; do
	if test `echo "${het0}==${het2}" |bc -l` -eq 1; then
	    match2=1
	fi
    done

    if test `echo "${match1}==1 && ${match2}==0" |bc -l` -eq 1; then
	awk 'BEGIN{FS=""}{print $23$24$25$26$27" "$0}' $FIL1 |\
            awk '($1=='${het1}'){print $0}' |\
            awk 'BEGIN{FS=""}{for(a=7;a<=NF;a++){printf $a;if(a==NF){printf "\n"}}}' >tempgap2441
	allosmod translatepdb -- tempgap2441 200 0 0 >>$FIL2
	rm tempgap2441
    fi
    if test `echo "${match1}==0 && ${match2}==1" |bc -l` -eq 1; then
	awk 'BEGIN{FS=""}{print $23$24$25$26$27" "$0}' $FIL2 |\
	    awk '($1=='${het2}'){print $0}' |\
	    awk 'BEGIN{FS=""}{for(a=7;a<=NF;a++){printf $a;if(a==NF){printf "\n"}}}' >tempgap2441
	allosmod translatepdb -- tempgap2441 200 0 0 >>$FIL1
	rm tempgap2441
    fi
done

rm tempfit.pdb 1temp88.pdb 2temp88.pdb 3temp88.pdb 1dsio.ali temp9991 temp9992 modeller.in
rm temp5773.ali 2temp88_fit.pdb modeller.in.log modeller2.in modeller2.in.log tempfit_sub.pdb qi_1.dat qi_avg.dat
rm temp[ade] temp994[12] tempN tempNN tempdd 1temp88_fit.pdb tempfit_fit.pdb
