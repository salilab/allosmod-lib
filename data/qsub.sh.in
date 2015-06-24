
RUNDIR=`pwd`
OUTDIR=$RUNDIR/DDD/XASPDB_${jobname}

mkdir -p $OUTDIR
echo $OUTDIR

# Get the 'allosmod' binary in the path
module load allosmod

# Make a temporary directory on the local scratch disk,
# specific to the user and SGE job.
TMPDIR="@LOCAL_SCRATCH@/XASPDB_${jobname}_$JOB_ID/"
mkdir -p $TMPDIR
echo $TMPDIR

# Copy input files to $TMPDIR here...
awk '{print "cp "$1" '${TMPDIR}'"}' list |sh
cp align.ali $TMPDIR/align.ali
cp list $TMPDIR/
cp lig.pdb $TMPDIR/
# If optional files exit, copy without newline problems
if test -e glyc.dat; then awk '{print $0}' glyc.dat | sed 's/\o015//g' > $TMPDIR/glyc.dat; fi
if test -e allosmod.py; then awk '{print $0}' allosmod.py | sed 's/\o015//g' > $TMPDIR/allosmod.py; fi
if test -e break.dat; then awk '{print $0}' break.dat | sed 's/\o015//g' > $TMPDIR/break.dat; fi
if test -e model_run.py; then awk '{print $0}' model_run.py | sed 's/\o015//g' > $TMPDIR/model_run.py; fi

########## do the actual job 
cd $TMPDIR
echo $TMPDIR >tempdir
hostname
date >>run.log
date

#########
#check alignment file and pdb files
#########
#if no chains specified, set to A
CHAIN_PM=`grep pm.pdb align.ali | grep structureX | awk 'BEGIN{FS=":"}{print $4}' | head -n1`
if test -z $CHAIN_PM; then CHAIN_PM=`grep pm.pdb align.ali | grep sequence | awk 'BEGIN{FS=":"}{print $4}' | head -n1`; fi
if test -z $CHAIN_PM; then
    CHAIN_PM=A
    LINE2REP=`grep -n pm.pdb align.ali | grep structureX | awk 'BEGIN{FS=":"}{print $1}' | head -n1`
    if test -z $LINE2REP; then LINE2REP=`grep -n pm.pdb align.ali | grep sequence | awk 'BEGIN{FS=":"}{print $1}' | head -n1`; fi
    if test -z $LINE2REP; then echo entry for pm.pdb has error >>${OUTDIR}/error.log; mv * $OUTDIR; exit ; fi
    awk '(NR<'${LINE2REP}'){print $0}' align.ali >tempali0
    awk 'BEGIN{FS=":"}(NR=='${LINE2REP}'){print $1":"$2":"$3":A"$4":"$5":A"$6":"$7":"$8":"$9":"$10}' align.ali >>tempali0
    awk '(NR>'${LINE2REP}'){print $0}' align.ali >>tempali0
    mv tempali0 align.ali
fi
#set from 1 to NRES
IND_PM=1
LINE2REP=`grep -n pm.pdb align.ali | grep structureX | awk 'BEGIN{FS=":"}{print $1}' | head -n1`
if test -z $LINE2REP; then LINE2REP=`grep -n pm.pdb align.ali | grep sequence | awk 'BEGIN{FS=":"}{print $1}' | head -n1`; fi
if test -z $LINE2REP; then echo entry for pm.pdb has error >>${OUTDIR}/error.log; mv * $OUTDIR; exit ; fi
IND_PM2=`cat -v align.ali | sed "s/\^M//g" | awk 'BEGIN{FS="";ctr=0}(NR>'${LINE2REP}'){for(a=1;a<=NF;a++){if($a!="/"&&$a!="-"&&$a!="*"&&$a!=" "){ctr+=1}\
             if($a=="*"){exit}}}END{print ctr}'`
awk '(NR<'${LINE2REP}'){print $0}' align.ali >tempali0
awk 'BEGIN{FS=":"}(NR=='${LINE2REP}'){print $1":"$2":"'${IND_PM}'":"$4":"'${IND_PM2}'":"$6":"$7":"$8":"$9":"$10}' align.ali >>tempali0
awk '(NR>'${LINE2REP}'){print $0}' align.ali >>tempali0
mv tempali0 align.ali
#test for hetatms in pm.pdb entry, if not filter out all hetatm entries... right now looking for block and heme entries only
isHET=`awk 'BEGIN{FS="";isHET=0}{for(a=1;a<=NF;a++){if($a==">"){tit=1}if(tit>2&&($a=="."||$a=="h")){isHET=1}if(a==NF){tit+=1}}}END{print isHET}' align.ali`
if test `echo "${isHET}==0" |bc -l` -eq 1; then
    for s in `cat list`; do
	awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"){print $0}' $s > tempiq4; mv tempiq4 $s
    done
    awk 'BEGIN{FS=""}{for(a=1;a<=NF;a++){if($a==">"){tit=1}if(tit<3||$a!="."){printf $a}\
         if(a==NF){tit+=1;printf "\n"}}}(NF==0){printf "\n"}' align.ali > tempiq4; mv tempiq4 align.ali
fi
#modify residues that Modeller cannot recognize
for s in `cat list`; do
    allosmod pdb_fix_res --in-place $s
done

#set up random numbers for structure generation
RR=`date --rfc-3339=ns | awk -F. '{print $2}' | awk -F- '{print $1}'`
RAND_NUM=`echo "scale=0; -1*(155421*${JOB_ID}+${RR})%40000-2" |bc -l | awk '{printf "%9.0f",$1}' | awk '{print $1}'`
RR1=`echo "scale=0; (155421*${JOB_ID}+${RR})%360" |bc -l | awk '{printf "%9.0f",$1}' | awk '{print $1}'`
RR2=`echo "scale=0; (5433*${JOB_ID}+2*${RR})%360" |bc -l | awk '{printf "%9.0f",$1}' | awk '{print $1}'`
RR3=`echo "scale=0; (11*${JOB_ID}+3*${RR})%360" |bc -l | awk '{printf "%9.0f",$1}' | awk '{print $1}'`

##########
### if AllosMod landscape, set up landscape here ###
##########
if test XSAMP == "simulation" -o XSAMP == "moderate_am"; then

##########
#make initial structure: 1) randomized, interpolation of pdbs if NTOT>1 or ligand binding site is modeled 2) single model wtih rotations/translations if NTOT=1
##########
NTOT=`allosmod count_alignments align.ali list pm.pdb`
isLIGMOD=`awk 'BEGIN{a=1}($3=="XX"&&NR==1){a=0}END{print a}' lig.pdb`
if test `echo "${NTOT}>1" |bc -l` -eq 1 -o `echo "${isLIGMOD}==1" |bc -l` -eq 1; then
    #generate pm's for input structures
    ctr=0
    for s in `cat list`; do
	echo $s >listin
echo making input structure $s
	allosmod get_pm_initialstruct --target pm.pdb --keep-alignment \
                                      align.ali listin ./ 1 slow
	awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $0}'  pred_${s}/pm.pdb.B99990001.pdb >pm_${s}
	if (test ! -s pm_${s}); then
	    echo "MODELLER has failed to create an initial model of the following structure: "${s} >>${OUTDIR}/error.log
	    echo "Perhaps there is an issue with the alignment file, check model_ini0.log in dir: "pred_${s} >>${OUTDIR}/error.log
	    echo "The alignment headers and sequence should look something like: " >>${OUTDIR}/error.log
	    allosmod pdb2ali $s >>${OUTDIR}/error.log
	    echo "Also, double check that the alignments themselves make sense." >>${OUTDIR}/error.log
	    break
	fi
echo completed input structure $s
	#align onto ligpdb so can be visualed with ligand
	allosmod salign0 XLPDB pm_${s}
	PMFIT=`echo pm_${s} | awk 'BEGIN{FS=""}{if($(NF-3)$(NF-2)$(NF-1)$NF==".pdb"){for(a=1;a<=NF-4;a++){printf $a}}else{printf $0}}END{print "_fit.pdb"}'`
	if (test -s $PMFIT); then
	    awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $0}' ${PMFIT} >pm_${s}
	else
	    echo pm_${s} not aligned, continue anyways
	fi

	ctr=$((${ctr} + 1))
    done

    #exit here if initial structures not generated
    if (test -e ${OUTDIR}/error.log); then
	echo "" >>${OUTDIR}/error.log
	allosmod get_auto_align align.ali pm.pdb list align_suggested.ali
	echo "See align_suggested.ali for a possible alignment file.  " >>${OUTDIR}/error.log
	echo "***WARNING*** Small errors in the alignment can cause big errors during a simulation due to energy conservation problems. " >>${OUTDIR}/error.log
	echo "Make sure there are no misalignments in which adjacent residues are aligned far apart \
	    (alignment programs often do this at the beginning or end of chains)." >>${OUTDIR}/error.log
	mv * $OUTDIR; exit
    fi

    @SCRIPT_DIR@/getavgpdb2.sh pm_XASPDB pm_XOTHPDB XASPDB XOTHPDB
    cp list list4contacts

else #NTOT=1
    echo generate avgpdb.pdb using all structures in list, which are translated/rotated
    echo generate avgpdb.pdb using all structures in list, which are translated/rotated >>run.log

    #translate and rotate input pdbs
    XVECT=( 0 1 0 0 -1 0 0 )
    YVECT=( 0 0 1 0 0 -1 0 )
    ZVECT=( 0 0 0 1 0 0 -1 )
    ctr=-1; iVECT=0
    for s in `cat list`; do
	ctr=$((${ctr} + 1))
	RR=`date --rfc-3339=ns | awk -F. '{print $2}' | awk -F- '{print $1}'`
	RAND_NUM=`echo "scale=0; -1*(155421*${JOB_ID}+${RR})%40000-2" |bc -l | awk '{printf "%9.0f",$1}' | awk '{print $1}'`
	RR1=`echo "scale=0; ((${ctr}+1)*${JOB_ID}+${RR})%360" |bc -l | awk '{printf "%9.0f",$1}' | awk '{print $1}'`
	RR2=`echo "scale=0; ((${ctr}+13)*${JOB_ID}+2*${RR})%360" |bc -l | awk '{printf "%9.0f",$1}' | awk '{print $1}'`
	RR3=`echo "scale=0; ((${ctr}+23)*${JOB_ID}+3*${RR})%360" |bc -l | awk '{printf "%9.0f",$1}' | awk '{print $1}'`
	#rotate
echo $RAND_NUM $RR1 $RR2 $RR3
	allosmod rotatepdb -- $s $RR1 $RR2 $RR3 >random.ini
	#translate
	RG[$ctr]=`allosmod getrofg random.ini`
	COFM=(`allosmod getcofm random.ini`)
	DIST=`echo "1.3*(${RG[0]}+${RG[${ctr}]})" | bc -l`
	iD=`echo $ctr | awk '{if($1<7){print $1%7}else{print ($1-7)%6+1}}'`
	if test `echo "${iD}==1" |bc -l` -eq 1; then iVECT=$((${iVECT} + 1)); fi
echo $iD $iVECT $DX $DY $DZ
	DX=`echo "-1*(${COFM[0]})+(${XVECT[${iD}]}*${iVECT}*${DIST})" | bc -l`
	DY=`echo "-1*(${COFM[1]})+(${YVECT[${iD}]}*${iVECT}*${DIST})" | bc -l`
	DZ=`echo "-1*(${COFM[2]})+(${ZVECT[${iD}]}*${iVECT}*${DIST})" | bc -l`
	allosmod translatepdb -- random.ini $DX $DY $DZ > $s
	rm random.ini
	
	#reset lig.pdb to c of m
	awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"){{for(a=31;a<=38;a++){printf $a}}{printf " "}\
            {for(a=39;a<=46;a++){printf $a}}{printf " "}\
            {for(a=47;a<=54;a++){printf $a}}{printf " \n"}}' XLPDB |\
            awk 'BEGIN{a=0;b=0;c=0}{a+=$1;b+=$2;c+=$3}END\
            {printf "ATOM      1  XX  ALA A   1    %8.3f%8.3f%8.3f  1.00 99.99           C\n",a/NR,b/NR,c/NR}' >lig.pdb
    done
    #generate input structure
    cp list listin
    s=`head -n1 list`
    allosmod get_pm_initialstruct --target pm.pdb --keep-alignment \
                                  align.ali listin ./ 1 slow
    awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $0}'  pred_${s}/pm.pdb.B99990001.pdb >avgpdb.pdb
    if (test ! -s avgpdb.pdb); then
	echo "MODELLER has failed to create an initial model of the following structure: avgpdb.pdb" >>${OUTDIR}/error.log
	echo "Perhaps there is an issue with the alignment file, check model_ini0.log in dir: "pred_${s} >>${OUTDIR}/error.log
	mv * $OUTDIR; exit
    fi

    cp avgpdb.pdb pm_XASPDB
    cp avgpdb.pdb pm_XOTHPDB
    echo XASPDB >list4contacts
fi
echo input structure completed

############
#testfirst: if first job setup initial files, else copy from directory
############
if test `echo "${jobname}==0" |bc -l` -eq 1; then
#generate allosteric site
allosmod get_allosteric_site --output_pdb allostericsite.pdb \
               --atom_list atomlistASRS XLPDB lig.pdb pm_XLPDB XrAS

############
#get restraints
############
TEST_NUC=`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"){print $18$19$20}' pm_XASPDB | sort -u | awk 'BEGIN{test=""}\
    ($1=="ADE"||$1=="CYT"||$1=="GUA"||$1=="THY"||$1=="DA"||$1=="DC"||$1=="DG"||$1=="DT"){test="--nucleic-acids"}END{print test}'` #use longer sc-sc dist for nucleotides
if test "XrAS" != "1000"; then #use ASPDB restraints if AS defined
    echo XASPDB >listin
    allosmod make_mod_inputs ${TEST_NUC} -- pm.pdb listin \
                             -3333 3 3 3 XDEV > model_ini.log
    mv pm.pdb.rsr listAS.rsr
fi
awk '(NF>0){print $0}' list >listin
allosmod make_mod_inputs ${TEST_NUC} -- pm.pdb listin \
                         -3333 3 3 3 XDEV > model_ini.log
mv pm.pdb.rsr listOTH.rsr
if test "XrAS" == "1000"; then cp listOTH.rsr listAS.rsr; fi #use all restraints if AS not defined

#if delEmax was unspecified in input.dat, calculate here (ignore DNA/RNA because handled differently) using hypothetical E landscape
delEmax=XdE
if test "XdE" == "CALC"; then
  echo "Determining delEmax" >>run.log
  awk 'BEGIN{FS=""}{a=$18$19$20}(a!="ADE"&&a!="  A"&&a!=" DA"&&a!="THY"&&a!="  T"&&a!=" DT"&&a!="URA"&&a!="  U"&&\
       a!=" DU"&&a!="GUA"&&a!="  G"&&a!=" DG"&&a!="CYT"&&a!="  C"&&a!=" DC"&&($1$2$3$4=="ATOM"||$1$2$3$4=="HETA")){print $0}' pm_XASPDB >tempiq7781
  awk 'BEGIN{FS=""}{a=$18$19$20}(a!="ADE"&&a!="  A"&&a!=" DA"&&a!="THY"&&a!="  T"&&a!=" DT"&&a!="URA"&&a!="  U"&&\
       a!=" DU"&&a!="GUA"&&a!="  G"&&a!=" DG"&&a!="CYT"&&a!="  C"&&a!=" DC"&&($1$2$3$4=="ATOM"||$1$2$3$4=="HETA")){print $0}' pm_XOTHPDB >tempiq7782

  if (test -s tempiq7781); then
    NATOM1=`awk 'END{print NR}' pm_XASPDB`
    NATOM2=`awk 'END{print NR}' tempiq7781`
    if test `echo "${NATOM1}==${NATOM2}" |bc -l` -eq 1; then 
	NTOT=`allosmod count_alignments align.ali list pm.pdb`
	allosmod edit_restraints --sigma_AS=2.0 --sigma_RS=2.0 \
                 --sigma_inter=2.0 --cutoff=11.0 --ntotal=$NTOT --delEmax=0.1 \
                 listOTH.rsr listAS.rsr tempiq7781 list4contacts \
                 atomlistASRS > edited.rsr
    else #redo steps without nucleotides
	echo redo restraints without nucleotides
        allosmod get_allosteric_site --atom_list atomlistASRS2 \
                                  XLPDB lig.pdb tempiq7781 XrAS
	cp align.ali align.ali.bak
	allosmod pdb2ali tempiq7781 >>align.ali
	echo tempiq7781 >listin
	allosmod make_mod_inputs -- tempiq7781 listin \
                                 -3333 3 3 3 XDEV > model_ini.log
	mv tempiq7781.rsr listAS2.rsr
	allosmod pdb2ali tempiq7782 >>align.ali
	echo tempiq7782 >listin
	allosmod make_mod_inputs -- tempiq7782 listin \
                                 -3333 3 3 3 XDEV > model_ini.log
	mv tempiq7782.rsr listOTH2.rsr #will this make sense in all cases?

	NTOT=`allosmod count_alignments align.ali list pm.pdb`
	mv align.ali.bak align.ali

	allosmod edit_restraints --sigma_AS=2.0 --sigma_RS=2.0 \
                 --sigma_inter=2.0 --cutoff=11.0 --ntotal=$NTOT --delEmax=0.1 \
                 listOTH2.rsr listAS2.rsr tempiq7781 list4contacts \
                 atomlistASRS2 > edited.rsr

	rm listOTH2.rsr listAS2.rsr atomlistASRS2 tempiq778[12].ini
    fi
    NHET=`awk 'BEGIN{FS=""}($1$2$3$4=="HETA"){print $18$19$20,$22,$23$24$25$26}' tempiq7781 |\
          awk 'BEGIN{nhet=0;li="";lc="";lt==""}{if($1!=lt||$2!=lc||$3!=li){nhet+=1;lt=$1;lc=$2;li=$3}}END{print nhet}'`
    NRES=`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"){print $23$24$25$26}' tempiq7781 | awk 'END{print $1}'`
    #if no distance interactions, delEmax defaults to 0.1
    delEmax=`awk '($2==50){a+=1}END{if(a>0){printf "%9.3f",3.6*('${NRES}'+'${NHET}')/a}else{print 0.1}}' edited.rsr | awk '{print $1}'`
    NDISTCONT=`awk '($2==50){a+=1}END{print a}' edited.rsr`
    TF=`echo "367*${delEmax}*${NDISTCONT}/(4*(${NRES}+${NHET}))" | bc -l | awk '{printf "%9.0f",$1}' | awk '{print $1}'`
  else
    delEmax=999
    NDISTCONT=NA
    TF=NA
  fi

  if test -z $delEmax; then echo "delEmax not correctly calculated" >>${OUTDIR}/error.log; mv * $OUTDIR; exit; fi
  echo "delEmax set to: "$delEmax >>run.log
  echo "Number of distance contacts: "$NDISTCONT >>run.log
  echo "Estimated Tf: "$TF >>run.log

  rm edited.rsr tempiq778[12]
fi

#handle break.dat options
if test "XBREAK" == "true"; then
    awk '(NF>0){print $0}' list >listin
    #redo here because nucleotides change max_sc_sc_distance
    allosmod make_mod_inputs -- pm.pdb listin \
                       -3333 3 3 3 XDEV > model_ini.log
    if [ "XCHEMFR" = "cdensity" ]; then
      allosmod contpres --cdensity_cutoff XZCUTOFF pm.pdb.rsr pm_XASPDB XSCLBREAK
    else
      allosmod contpres pm.pdb.rsr pm_XASPDB XSCLBREAK
    fi
    echo "Chemical frustration is implemented: contacts with buried acidic/basic residues are scaled by XSCLBREAK, z-score cutoff is XZCUTOFF" >>run.log
    echo "Chemical frustration type is XCHEMFR" >>run.log
fi
if test "XLOCRIGID" == "true"; then
    @SCRIPT_DIR@/get_loopadjres.sh #strengthen residues adjacent to loops
fi
if (test -e break.dat); then
    break="break.dat"
else
    break=""
fi

########
#merge restraint files and add additional restraints
########
if (test ! -e ${OUTDIR}/error.log); then 
    NTOT=`allosmod count_alignments align.ali list pm.pdb`
    allosmod edit_restraints --sigma_AS=2.0 --sigma_RS=2.0 \
             --sigma_inter=2.0 --cutoff=11.0 --ntotal=$NTOT \
             --delEmax=${delEmax} XCOARSE XLOCRIGID ${break} \
             listOTH.rsr listAS.rsr pm_XASPDB list4contacts \
             atomlistASRS > edited.rsr
    #add restraints between protein and sugar
    if test `echo "XGLYC2==1" |bc -l` -eq 1; then 
	@SCRIPT_DIR@/get_glyc_restraint.sh pm_XASPDB allosmod.py >>edited.rsr
    fi
    #add restraints for bonds and distance upper/lower bounds
    allosmod get_add_restraint ${RUNDIR}/input.dat pm_XASPDB HARM >tempaddrestr
    allosmod get_add_restraint ${RUNDIR}/input.dat pm_XASPDB UPBD >>tempaddrestr
    allosmod get_add_restraint ${RUNDIR}/input.dat pm_XASPDB LOBD >>tempaddrestr
    if (test -s tempaddrestr); then
	cat tempaddrestr >>edited.rsr
	echo "" >>run.log; echo "additional restraints: " >>run.log; cat tempaddrestr >>run.log; echo "" >>run.log
	#redo initial structure with additional restraints
	for dir in `ls -1d pred_*`; do
	    cd $dir
	    awk '($2==2||$2==3){print $0}' ../tempaddrestr >> pm.pdb.rsr
	    awk '($2==1){print $0}' ../tempaddrestr | awk 'BEGIN{FS=""}{for(a=1;a<=45;a++){printf $a;if(a==45){printf "   20.0000    0.0100\n"}}}' >> pm.pdb.rsr
	    mv pm.pdb.rsr edit.rsr
	    cat model_ini0.py | sed "s/'align.ali',/'align.ali', csrfile = 'edit.rsr',/g" >tempiqmi5; mv tempiqmi5 model_ini0.py
	    /salilab/diva1/home/modeller/modSVN model_ini0.py
	    cd ..
	done
	if test `echo "${NTOT}>1" |bc -l` -eq 1 -o `echo "${isLIGMOD}==1" |bc -l` -eq 1; then
	    for s in `cat list`; do
		awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $0}'  pred_${s}/pm.pdb.B99990001.pdb >pm_${s}
	    done
	    @SCRIPT_DIR@/getavgpdb2.sh pm_XASPDB pm_XOTHPDB XASPDB XOTHPDB
	else
	    s=`head -n1 list`
	    awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $0}'  pred_${s}/pm.pdb.B99990001.pdb >avgpdb.pdb
	    cp avgpdb.pdb pm_XASPDB
	    cp avgpdb.pdb pm_XOTHPDB
	fi
    fi
    rm tempaddrestr
fi

else #notfirst: use previously calculated atomlistASRS and restraints
delEmax=XdE
echo searching for converted.rsr
date
for c in `@SCRIPT_DIR@/count.pl 1 60`; do
    echo $c
    FIRSTDIR=${RUNDIR}/pred_dEXdErASXrAS/XASPDB_0
    echo $FIRSTDIR
    ls $FIRSTDIR
    TESTF=`ls -1 ${FIRSTDIR}/converted.rsr | awk '(NR==1){print "basename "$1}' |sh`
    echo "test: "$TESTF
    if (test ! -z $TESTF); then
	cp ${FIRSTDIR}/atomlistASRS .
	cp ${FIRSTDIR}/converted.rsr .
	cp ${FIRSTDIR}/allostericsite.pdb .
	cp ${FIRSTDIR}/contacts.dat .
	if (test -e ${FIRSTDIR}/allosmod.py); then cp ${FIRSTDIR}/allosmod.py .; fi
	echo cp done
	break
    else
	echo sleep
	sleep 2m
    fi
    echo continue
    if test `echo "${c}==60" |bc -l` -eq 1; then echo restraints not found >> ${OUTDIR}/error.log; mv * $OUTDIR; exit; fi
done
fi #testfirst

echo edited.rsr file complete

#catch errors
if (test -e ${OUTDIR}/error.log); then mv * $OUTDIR; exit; fi

#initialize starting structure
allosmod make_mod_inputs -- pm.pdb list $RAND_NUM \
                         $RR1 $RR2 $RR3 XDEV > model_ini.log
allosmod setchain random.ini $CHAIN_PM | awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $0}' >tempini
mv tempini random.ini
if (test ! -s random.ini); then echo structure not initialized >> ${OUTDIR}/error.log; mv * $OUTDIR; exit; fi

#convert restraints to splines
if test `echo "${jobname}==0" |bc -l` -eq 1; then
    if test `echo "XGLYC1==0" |bc -l` -eq 1; then #skip if hydrogens are needed
	echo convert to splines 
	allosmod spline random.ini edited.rsr converted.rsr
	cp converted.rsr atomlistASRS allostericsite.pdb contacts.dat ${RUNDIR}/pred_dEXdErASXrAS/XASPDB_0
	if (test -e allosmod.py); then cp allosmod.py ${RUNDIR}/pred_dEXdErASXrAS/XASPDB_0; fi
    fi
    echo converted.rsr file complete, copied to ${RUNDIR}/pred_dEXdErASXrAS/XASPDB_0
    ls ${RUNDIR}/pred_dEXdErASXrAS/XASPDB_0
    date
fi
if (test ! -s converted.rsr); then echo restraints failed to be converted >> ${OUTDIR}/error.log; mv * $OUTDIR; exit; fi

fi ### end setup AllosMod landscape ###

#get AllosModel subclasses from allosmod-lib
if (test ! -e allosmod.py); then 
    ALLOSMOD_PY=`python -c 'import allosmod.modeller, os; print(os.path.splitext(allosmod.modeller.__file__)[0])'`
    cp ${ALLOSMOD_PY}.py allosmod.py
fi
MDTEMP=`echo XMDTEMP | tr [A-Z] [a-z] | awk '{if($1=="scan"){print ('${jobname}'%5)*50+300.0}else{print $1}}'`

#generate scripts and structures
if test `echo "XGLYC1==0" |bc -l` -eq 1; then
    if test "XSAMP" == "simulation"; then
	echo "randomized numbers: $RAND_NUM $RR1 $RR2 $RR3" >>run.log
	@SCRIPT_DIR@/get_pm2.sh pm.pdb list $RAND_NUM XDEV script $MDTEMP
	cp @SCRIPT_DIR@/README_user $RUNDIR/README
	if test "XSCRAPP" == "true"; then
	    python model_run.py > model_run.log
	    JOBID=`echo $RUNDIR | awk 'BEGIN{FS="/"}{print $(NF-1)}'`
	    OUTDIR=@GLOBAL_SCRATCH@/${JOBID}/DDD/XASPDB_${jobname}
	    mkdir -p $OUTDIR
	    echo @GLOBAL_SCRATCH@/${JOBID}/DDD > $RUNDIR/DDD/XASPDB_${jobname}/OUTPUT_IS_HERE
	fi
    elif test "XSAMP" == "moderate_cm" -o "XSAMP" == "moderate_am" -o "XSAMP" == "fast_cm"; then
	@SCRIPT_DIR@/get_pm2.sh pm.pdb list $RAND_NUM XDEV XSAMP $MDTEMP
	echo running modeller
	python model_run.py > model_run.log
	echo done with modeller
#    elif test "XSAMP" == "moderate_cm_simulation"; then
#	@SCRIPT_DIR@/get_pm2.sh pm.pdb list run $RAND_NUM $RR1 $RR2 $RR3 XDEV XSAMP
#	python model_run.py > model_run.log
#	cp pm.pdb.B99990001.pdb random.ini
#	@SCRIPT_DIR@/get_pm2.sh pm.pdb list run $RAND_NUM $RR1 $RR2 $RR3 0.0 script $MDTEMP
    else
	echo "sampling type not properly defined" >>${OUTDIR}/error.log; mv * $OUTDIR; exit
    fi
else #glycosylation
    ls -lrt
    cp @SCRIPT_DIR@/[A-Z][A-Z][A-Z].rsr .
    @SCRIPT_DIR@/get_pm_glyc.sh pm.pdb list $RAND_NUM XREP_OPT XATT_GAP script $MDTEMP
    #generate protein model with correct chains
    python model_ini0.py > model_ini0.log
    CHAIN=`awk 'BEGIN{FS=""}($1$2$3$4=="ATOM"||$1$2$3$4=="HETA"){print $22$23}' pm.pdb.B99990001.pdb | awk '{print $1}' | sort -u | head -n1`
    if test -z $CHAIN; then
	NC=A #`grep pm.pdb align.ali | grep -v P1 | awk 'BEGIN{FS=":"}{print $4}' | awk 'BEGIN{a="A"}(NR==1){a=$1}{print a}'`
	allosmod setchain pm.pdb.B99990001.pdb --in-place $NC
    fi
    #generate initial glycosylation model and restraints
    @SCRIPT_DIR@/get_pm_glyc.sh pm.pdb list $RAND_NUM XREP_OPT XATT_GAP pm.pdb.B99990001.pdb 
    python model_ini.py > model_ini.log
    @SCRIPT_DIR@/get_rest.sh pm.pdb.B99990001.pdb >> pm.pdb.rsr
    mv pm.pdb.rsr converted.rsr
    #sample glycosylated structures
    python model_glyc.py > model_glyc.log
#    if test "XSAMP" == "moderate_cm_simulation"; then #to run simulations with sugar
#	cp pm.pdb.B99990001.pdb random.ini
#	cp converted.rsr tempsav.rsr
#	cp @SCRIPT_DIR@/allosmod.py ./allosmod2.py
#	#convert protein restraints (with hydrogens)
#	python model_ini0.py > model_ini0.log
#	MAX_PROT_AIND=`awk '($6>1){print $9"\n"$10}' pm.pdb.rsr | sort -nk1 | tail -n1`
#	grep ATOM pm.pdb.B99990001.pdb >pm_rand.pdb
#	awk '{print NR" AS"}' pm_rand.pdb >atomlistASRS
#	@SCRIPT_DIR@/editrestraints2.sh pm.pdb.rsr pm.pdb.rsr pm_rand.pdb pm_rand.pdb atomlistASRS 2.0,2.0,2.0 11.0 XNTOT ${delEmax} >>run.log
#	allosmod spline pm.pdb edited.rsr converted.rsr
#	#get glyc restraints
#	awk '($6==2&&($9>'${MAX_PROT_AIND}'||$10>'${MAX_PROT_AIND}')){print $0}' tempsav.rsr >>converted.rsr
#	awk '($6==3&&($9>'${MAX_PROT_AIND}'||$10>'${MAX_PROT_AIND}'||$11>'${MAX_PROT_AIND}')){print $0}' tempsav.rsr >>converted.rsr
#	awk '($6==4&&($9>'${MAX_PROT_AIND}'||$10>'${MAX_PROT_AIND}'||$11>'${MAX_PROT_AIND}'||$12>'${MAX_PROT_AIND}')){print $0}' tempsav.rsr >>converted.rsr
#	#generate simulation script
#	@SCRIPT_DIR@/get_pm2.sh pm.pdb list run $RAND_NUM $RR1 $RR2 $RR3 0.0 script_glyc $MDTEMP
#	rm pm_rand.pdb atomlistASRS model_ini0.py allosmod2.py tempsav.rsr
#    fi
    rm [A-Z][A-Z][A-Z].rsr get_rest.in
fi

date >>run.log
date

# Copy back output files from $TMPDIR here...
rm get_pm.params
if test -e targlist; then rm targlist; fi
if test -e model_ini.log; then rm model_ini.log; fi
rm pm.pdb.rsr crap
rm *fit.pdb listin list4contacts
rm edited.rsr
rm pm.pdb.B0*.pdb model_ini.py model_ini0.py 
rm list*rsr avgpdb.pdb

#zip -r out.zip *
#mv out.zip $OUTDIR

#cd $OUTDIR
#sleep 10s

#rm -rf $TMPDIR
#unzip out.zip
#sleep 10s
#rm out.zip
#sleep 10s



mv * $OUTDIR

cd $OUTDIR
sleep 10s

rm -rf $TMPDIR
sleep 10s
