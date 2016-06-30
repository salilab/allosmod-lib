#!/bin/bash
# print restraint between protein and sugars

PM=$1 #pm.pdb file
ALLOSPY=$2 #allosmod.py with patches, to be converted to harmonic restraints

handle_missing_residue() {
  local resnum="$1"
  echo "Could not parse glycosylation patch in $ALLOSPY;" >&2
  echo "residue $resnum not found in $PM." >&2
  echo "Note that residue numbering in the model may not be the same" >&2
  echo "as that in your input structure(s), particularly if those" >&2
  echo "inputs have missing residues." >&2
}

R_CA=( `grep -A1 GLB ${ALLOSPY}| grep residues | sed "s/'//g" | awk 'BEGIN{FS="["}{print $2}' | awk 'BEGIN{FS=":"}(NR%2==1){print $1}'` \
       `grep -A1 GPA ${ALLOSPY}| grep residues | sed "s/'//g" | awk 'BEGIN{FS="["}{print $2}' | awk 'BEGIN{FS=":"}(NR%2==1){print $1}'` \ 
       `grep -A1 GPB ${ALLOSPY}| grep residues | sed "s/'//g" | awk 'BEGIN{FS="["}{print $2}' | awk 'BEGIN{FS=":"}(NR%2==1){print $1}'` )

R_C1=( `grep -A1 GLB ${ALLOSPY}| grep residues | sed "s/'//g" | awk 'BEGIN{FS="["}{print $2}' | awk 'BEGIN{FS=":"}(NR%2==0){print $1}'` \
       `grep -A1 GPA ${ALLOSPY}| grep residues | sed "s/'//g" | awk 'BEGIN{FS="["}{print $2}' | awk 'BEGIN{FS=":"}(NR%2==0){print $1}'` \ 
       `grep -A1 GPB ${ALLOSPY}| grep residues | sed "s/'//g" | awk 'BEGIN{FS="["}{print $2}' | awk 'BEGIN{FS=":"}(NR%2==0){print $1}'` )

if test `echo "${#R_CA[@]}!=${#R_C1[@]}" |bc -l` -eq 1; then exit; fi

ctr=-1
for i in ${R_CA[@]}; do
    ctr=$((${ctr} + 1))
    j=${R_C1[${ctr}]}

    iatm1=`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $23$24$25$26$27,$14$15}' $PM | awk '($1=='${i}'&&$2=="CA"){print NR}'`
    iatm2=`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $23$24$25$26$27,$14$15}' $PM | awk '($1=='${j}'&&$2=="C1"){print NR}'`
    if [ "$iatm1" = "" ]; then
      handle_missing_residue $i
      exit 1
    fi
    if [ "$iatm2" = "" ]; then
      handle_missing_residue $j
      exit 1
    fi

    echo $iatm1 $iatm2 | awk '{printf "R    3   1   1   1   2   2   1    %5.0f %5.0f    5.0000    0.0350\n",$1,$2}'

done

grep -v special_patches ${ALLOSPY} | grep -v self.residues >tempggr2041
mv tempggr2041 ${ALLOSPY}
