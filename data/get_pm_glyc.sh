#!/bin/bash

# Absolute path containing this and other scripts
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

TARG_SEQ=$1 #pm.pdb in alignment
LIST_KNOWNS=$2
RAND=$3 #random seed
REP_OPT=$4 #option to allow repeat optimization
ATT_GAP=`echo $5 | tr [A-Z] [a-z]` #inserts gaps to allow flexibility at glycosylation sites
GLYCPM=$6

get_array_index()
{
    #first arg is search_query, followed by expanded array. prints index of search_query in array, or nothing if not in array
    R_INP=(`echo $@`)
    search_query=${R_INP[0]}
    ctr=-2
    for i in ${R_INP[@]}; do
	ctr=$((${ctr} + 1))
	if test `echo "${ctr}>-1" |bc -l` -eq 1; then
	    if test $search_query == $i; then echo $ctr; break; fi
	fi
    done
}

# get pdb file names/chains
NFIL=`cat ${LIST_KNOWNS} |awk 'END{print NR-1}'`
F=(`cat ${LIST_KNOWNS}`)

#script to create initial model with sugar
cat << EOF > model_ini.py
from modeller import *
from modeller.scripts import complete_pdb
import allosmod

env =environ(rand_seed=${RAND}, restyp_lib_file='$SCRIPT_DIR/restyp.dat', copy=None)

# Read in HETATM records from template PDBs
env.io.atom_files_directory = ['.', '../atom_files']
env.libs.topology.read(file='$SCRIPT_DIR/top_all_glyco.lib')
env.libs.parameters.read(file='$SCRIPT_DIR/par_all_glyco.lib')
env.io.hetatm = True
env.io.hydrogen = True

a = allosmod.AllosModel(env, deviation=None, alnfile='align2.ali', 
EOF

echo -n "knowns=('${F[0]}'" >>model_ini.py
awk '(NR>1&&NF>0){printf ",@"$1"@"}' ${LIST_KNOWNS} | sed "s/@/'/g" >>model_ini.py
echo "), sequence='${TARG_SEQ}')" >>model_ini.py
echo >>model_ini.py

cat <<EOF >>model_ini.py
# Very thorough VTFM optimization:
a.library_schedule = allosmod.MDnone
a.max_var_iterations = 500

# MD optimization:
a.md_level = allosmod.none

# Repeat the whole cycle 1 time and do not stop unless obj.func. > 1E9
a.repeat_optimization = 1
a.max_molpdf = 1e9

a.starting_model = 1
a.ending_model = 1
a.make()
EOF

#script to create initial model with protein only
cat << EOF > model_ini0.py
from modeller import *
from modeller.scripts import complete_pdb
import allosmod

env =environ(rand_seed=${RAND}, restyp_lib_file='$SCRIPT_DIR/restyp.dat', copy=None)

# Read in HETATM records from template PDBs
env.io.atom_files_directory = ['.', '../atom_files']
env.libs.topology.read(file='$SCRIPT_DIR/top_all_glyco.lib')
env.libs.parameters.read(file='$SCRIPT_DIR/par_all_glyco.lib')
env.io.hetatm = True
env.io.hydrogen = True

a = allosmod.AllosModel(env, deviation=None, alnfile='align.ali', 
EOF

echo -n "knowns=('${F[0]}'" >>model_ini0.py
awk '(NR>1&&NF>0){printf ",@"$1"@"}' ${LIST_KNOWNS} | sed "s/@/'/g" >>model_ini0.py
echo "), sequence='${TARG_SEQ}')" >>model_ini0.py
echo >>model_ini0.py

cat <<EOF >>model_ini0.py
# Very thorough VTFM optimization:
a.library_schedule = allosmod.MDnone
a.max_var_iterations = 500

# MD optimization:
a.md_level = allosmod.none

# Repeat the whole cycle 1 time and do not stop unless obj.func. > 1E9
a.repeat_optimization = 1
a.max_molpdf = 1e9

a.starting_model = 1
a.ending_model = 1
a.make()
EOF

#script to generate glycosylated models
cat <<EOF >model_glyc.py
from modeller import *
from modeller.scripts import complete_pdb
import allosmod

env =environ(rand_seed=${RAND}, restyp_lib_file='$SCRIPT_DIR/restyp.dat', copy=None)

# Read in HETATM records from template PDBs
env.io.atom_files_directory = ['.', '../atom_files']
env.libs.topology.read(file='$SCRIPT_DIR/top_all_glyco.lib')
env.libs.parameters.read(file='$SCRIPT_DIR/par_all_glyco.lib')
env.io.hetatm = True
env.io.hydrogen = True

a = allosmod.AllosModel(env, csrfile='converted.rsr', deviation=None, alnfile='align2.ali', 
EOF

echo -n "knowns=('${F[0]}'" >>model_glyc.py
awk '(NR>1&&NF>0){printf ",@"$1"@"}' ${LIST_KNOWNS} | sed "s/@/'/g" >>model_glyc.py
echo "), sequence='${TARG_SEQ}')" >>model_glyc.py
echo >>model_glyc.py

cat <<EOF >>model_glyc.py
# Very thorough VTFM optimization:
a.library_schedule = allosmod.MDopt
a.max_var_iterations = 500

# MD optimization:
a.md_level = allosmod.moderate

# Repeat the whole cycle 1 time and do not stop unless obj.func. > 1E9
a.repeat_optimization = ${REP_OPT}
a.max_molpdf = 1e9

a.starting_model = 1
a.ending_model = 1
a.make()
EOF

if test $GLYCPM == "script"; then exit; fi

#add glycosidic bonds into allosmod.py
S1=`grep -n ${TARG_SEQ} align.ali | grep P1 | awk 'BEGIN{FS=":"}{print $1+2}'`
S2=`grep -n P1 align.ali | awk 'BEGIN{FS=":"}($1>'${S1}'){print $1-1;exit}'`
LAST_RES=`cat -v align.ali | sed "s/\^M//g" | awk 'BEGIN{FS=""}(NR>='${S1}'&&NR<='${S2}'){for(a=1;a<=NF;a++){print $a}}' | \
          awk '($1!="-"&&$1!="/"&&$1!="*"&&$1!="."){a+=1}($1=="*"){print a;exit}'`
ATTACH_PNTS=(`awk '($2=="NGLA"||$2=="NGLB"||$2=="SGPA"||$2=="SGPB"||$2=="TGPA"||$2=="TGPB"){print $3}' glyc.dat`)
ATTACH_IND=( `awk '($2=="NGLA"||$2=="NGLB"||$2=="SGPA"||$2=="SGPB"||$2=="TGPA"||$2=="TGPB"){print NR}END{print NR+1}' glyc.dat`)

  #relist chains from A,B,... for protein and sugars
R_CHAIN=( A B C D E F G H I J K L M N O P Q R S T U V W X Y Z ) #used to specify 1 chain for all sugars, if >26 protein chains, then error
R_C_PDB=(`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $22}' $GLYCPM | sort -u`)
iOFFSET=${#R_C_PDB[@]}

ctr=-1
LI=$LAST_RES
GLYC_STRING="?/"
if test -e get_rest.in; then rm get_rest.in; fi
echo "# define bonds between sugar monomers" >>allosmod.py
echo "    def special_patches(self, aln):" >>allosmod.py
for ap in ${ATTACH_PNTS[@]}; do
    ctr=$((${ctr} + 1))

    INP1=${ATTACH_IND[$ctr]}
    INP2=$((${ATTACH_IND[$(($ctr + 1))]} - 1))
#    INP_CHAIN=`awk '{if(NR=='${INP1}'){if(NF==4){print $4}else{print ""}}}' glyc.dat`
#    PROT_CHAIN=${R_CHAIN[`get_array_index $INP_CHAIN ${R_C_PDB[@]}`]}
    PROT_CHAIN=`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"){print $23$24$25$26$27,$22}' ${GLYCPM} | awk '($1=='${ap}'){print $2; exit}'`
    GCHAIN=${R_CHAIN[$((${iOFFSET}))]} #${R_CHAIN[$((${ctr} + ${iOFFSET}))]} #use same chain for all sugars, this allows >26 chains

    R_TYPE=(`awk '(NR>='${INP1}'&&NR<='${INP2}'&&NF>0){print $1}' glyc.dat | \
             awk '($1=="NAG"){print 1}($1=="MAN"){print 2}($1=="BMA"){print 3}($1=="GLB"){print 4}($1=="FUC"){print 5}($1=="NAN"){print 8}($1=="NGA"){print 9}'`)
    R_BOND=(`awk '(NR>='${INP1}'&&NR<='${INP2}'&&NF>0){print $2}' glyc.dat`)
    R_CONN=(`awk '(NR>='${INP1}'&&NR<='${INP2}'&&NF>0){print $3+'${LI}'}' glyc.dat`)
    if test ${#R_TYPE[@]} != ${#R_BOND[@]}; then echo "INCORRECT SUGAR TYPE SPECIFIED"; exit; fi

    #add bonds for all sugar monomers i
    i=$LI
    ctr2=-1
    for j in ${R_CONN[@]}; do
	ctr2=$((${ctr2} + 1))
	BOND=${R_BOND[$ctr2]}
	connect_atom=`echo $BOND | awk '($1=="NGLA"||$1=="NGLB"){print "ND2"}($1=="SGPA"||$1=="SGPB"){print "OG"}($1=="TGPA"||$1=="TGPB"){print "OG1"}'`
	i=$((${i} + 1))

	if test `echo "${ctr2}==0" |bc -l` -eq 1; then
	    echo "        self.patch(residue_type='"${BOND}"', residues=(self.residues['"${ap}":"${PROT_CHAIN}"']," >>allosmod.py
	    echo "                                                  self.residues['"${i}":"${GCHAIN}"']))" >>allosmod.py
	    echo ${ap} $connect_atom ${BOND} >>get_rest.in
	else
	    echo "        self.patch(residue_type='"${BOND}"', residues=(self.residues['"${j}":"${GCHAIN}"']," >>allosmod.py
	    echo "                                                  self.residues['"${i}":"${GCHAIN}"']))" >>allosmod.py
	     echo ${j} `echo $BOND | awk 'BEGIN{FS=""}{if($1==1){print "O"$2}else{print "O"$4}}'` ${BOND} >>get_rest.in
	fi
    done
    echo >>allosmod.py
#    GLYC_STRING=(${GLYC_STRING[@]} "?/" ${R_TYPE[@]})
    GLYC_STRING=(${GLYC_STRING[@]}${R_TYPE[@]})
    LI=$i
done

#make align2.ali: alignment file with sugar entries
GLYC_STRING=`echo ${GLYC_STRING[@]} | sed "s/ //g" | sed -e 's/?/\\\/g'`
if test -e align2.ali; then rm align2.ali; fi
R_SEQ=(${TARG_SEQ} `cat ${LIST_KNOWNS}`)
for seq in ${R_SEQ[@]}; do
    S1=`grep -n ${seq} align.ali | grep P1 | awk 'BEGIN{FS=":"}{print $1}'`
    S2=`grep -n P1 align.ali | awk 'BEGIN{FS=":"}($1>'${S1}'){print $1-1;exit}' | awk 'END{if(NF>0){print $0}else{print 99999}}'`
    
    awk '(NR>='${S1}'&&NR<='${S1}'+1){print $0}' align.ali >>align2.ali
    awk '(NR>='${S1}'+2&&NR<='${S2}'){print $0}' align.ali | sed "s/\*/${GLYC_STRING}\*/g" >>align2.ali
    GLYC_STRING=`echo ${GLYC_STRING} | sed -e 's/\\\//g' | awk 'BEGIN{FS=""}{for(a=1;a<=NF;a++){if($a=="/"){printf "?/"}else{printf "-"}}}' | sed -e 's/?/\\\/g'`
done

LINE2REP=`grep -n ${TARG_SEQ} align2.ali | grep P1 | awk 'BEGIN{FS=":"}{print $1+1}'`
NUMGLYC=`awk 'BEGIN{a=0}(NF>1){a+=1}END{print a}' glyc.dat`
IND_PM2=`awk 'BEGIN{FS=":"}(NR=='${LINE2REP}'){print $5+'${NUMGLYC}'}' align2.ali`
awk '(NR<'${LINE2REP}'){print $0}' align2.ali >tempali0
awk 'BEGIN{FS=":"}(NR=='${LINE2REP}'){print $1":"$2":"$3":"$4":"'${IND_PM2}'":"$6":"$7":"$8":"$9":"$10}' align2.ali >>tempali0
awk '(NR>'${LINE2REP}'){print $0}' align2.ali >>tempali0
mv tempali0 align2.ali

#add gaps near insertion sites
if test $ATT_GAP == "true"; then
    R_LB=(`awk '($2=="NGLA"||$2=="NGLB"||$2=="SGPA"||$2=="SGPB"||$2=="TGPA"||$2=="TGPB"){print $3,$3-2}' glyc.dat | sort -nk1 | awk '{print $2}'`)
    R_UB=(`awk '($2=="NGLA"||$2=="NGLB"||$2=="SGPA"||$2=="SGPB"||$2=="TGPA"||$2=="TGPB"){print $3,$3+2}' glyc.dat | sort -nk1 | awk '{print $2}'`)
    #test for LB and UB overlaps
    for i in `echo 1 2 3 4 5`; do
	ctr=-1; last_ub=-9999; NGLYC=${#R_LB[@]}
	OK=1
	for lb in ${R_LB[@]}; do
	    ctr=$((${ctr} + 1))
	    ub=${R_UB[$ctr]}
	    if test `echo "${lb}-(${last_ub})<=0" |bc -l` -eq 1; then
		R_LB[${ctr}]=$((${lb} + 1))
		R_UB[$((${ctr} - 1))]=$((${last_ub} - 1))
		OK=0
	    fi
	    last_ub=$ub
            if test `echo "${R_LB[${ctr}]}>${R_UB[${ctr}]}" |bc -l` -eq 1; then OK=0; fi
	done
	if test `echo "${OK}==1" |bc -l` -eq 1; then break; fi
    done

    #check for correct alignment insert positions
    if test `echo "${OK}==0" |bc -l` -eq 1; then
	echo "ERROR inserting loops into the alignment file, check input" >>error.log
	echo "ERROR inserting loops into the alignment file, check input"
    else
        #insert gaps
	ctr=-1
	for lb in ${R_LB[@]}; do
	    ctr=$((${ctr} + 1))
	    ub=${R_UB[$ctr]}
	    $SCRIPT_DIR/ali_insertgap.sh align2.ali 1 $lb $ub >tempgpg11
	    mv tempgpg11 align2.ali
	done
    fi
fi

