#!/bin/bash 
# 1) FILE1 2) FILE2_subset 3) FILE2_whole
# saligns FILE2_subset onto FILE1 and then superimposes FILE2_whole onto FILE2_subset

PDIR=./
cp ${PDIR}/$1 ./1temp88.pdb
cp ${PDIR}/$2 ./2temp88.pdb
cp ${PDIR}/$3 ./3temp88.pdb

FF1=1temp88.pdb
#CC1=`awk '(NR==1){print $5}' $FF1`
#CC1=`echo $1 | awk -F_ '{print $2}' | awk 'BEGIN{FS=""}{print $5}'`
FF2=2temp88.pdb
#CC2=`awk '(NR==1){print $5}' $FF2`
#CC2=`echo $2 | awk -F_ '{print $2}' | awk 'BEGIN{FS=""}{print $5}'`
#SNAME2=`echo $2 | awk -F_ '{print $2}' | awk 'BEGIN{FS=""}{print $1$2$3$4}'`
SFILE2=`echo $2 | awk 'BEGIN{FS="."}{print $1}'`
#SRES2=`echo $2 | awk -F_ '{print $3}'`
#ERES2=`echo $2 | awk -F_ '{print $4}' |awk -F. '{print $1}'`
#cp $FF2 3temp88.pdb
#awk 'BEGIN{FS=""}($22!=0){print $0}' ../allcont/all_${SNAME2}${CC2}_${SRES2}_${ERES2}.ent >>3temp88.pdb #no cryst contacts
FF3=3temp88.pdb

NCA=`awk 'BEGIN{FS="";a=0}($14$15=="CA"){a+=1}END{print a}' $FF1`
NP=`awk 'BEGIN{FS="";a=0}($14$15=="P "){a+=1}END{print a}' $FF1`
if test `echo "${NCA}>=${NP}" |bc -l` -eq 1; then
    OPT=CA
else
    OPT=P
fi

#get alignment file
if test -e list; then cp list list.bak; fi
ls -1 $FF2 $FF3 >list
allosmod get_auto_align align.ali pm.pdb list 1dsio.ali >>run.log
if test -e list.bak; then mv list.bak list; fi

cat <<EOF >modeller.in

from modeller import *

log.verbose()
env = environ()
env.io.atom_files_directory = ['.', '../atom_files']

aln = alignment(env)
#for (code, chain) in (('$FF1', '$CC1'), ('$FF2', '$CC2')):
#    mdl = model(env, file=code, model_segment=('FIRST:'+chain, 'LAST:'+chain)) 
for (code) in (('$FF1'), ('$FF2')):
    mdl = model(env, file=code)
    aln.append_model(mdl, atom_files=code, align_codes=code)

for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
               rr_file='\$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50),
               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
               dendrogram_file='temp5773.tree', fit_atoms='${OPT}',
               alignment_type='tree', # If 'progresive', the tree is not
                                      # computed and all structues will be
                                      # aligned sequentially to the first
               #ext_tree_file='1is3A_exmat.mtx', # Tree building can be avoided
                                                 # if the tree is input
               feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
               improve_alignment=True, fit=True, write_fit=write_fit,
               write_whole_pdb=whole, output='ALIGNMENT QUALITY')

# aln.write(file='2g75.pap', alignment_format='PAP')
aln.write(file='temp5773.ali', alignment_format='PIR')

# The number of equivalent positions at different RMS_CUTOFF values can be
# computed by changing the RMS value and keeping all feature weights = 0
#aln.salign(rms_cutoff=1.0,
#           normalize_pp_scores=False, rr_file='\$(LIB)/as1.sim.mat', overhang=30,
#           gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3), fit_atoms='${OPT}',
#           gap_gap_score=0, gap_residue_score=0, dendrogram_file='temp5773.tree',
#           alignment_type='progressive', feature_weights=[0]*6,
#           improve_alignment=False, fit=False, write_fit=True,
#           write_whole_pdb=False, output='QUALITY')
EOF

echo "run salign of FILE2_subset onto FILE1" >>run.log
python modeller.in >& modeller.in.log

if (test ! -e 2temp88_fit.pdb); then exit; fi

SCOR=`grep percentage modeller.in.log  |awk '(NR==4){print $4}'`
#SI=`./get_percid.sh temp5773.ali` 
echo "quality percentage: " $SCOR >>run.log
#echo "QS, %SI: " $FF1 $FF2 $SCOR $SI   

cat <<EOF >modeller2.in
from modeller import *

log.verbose()
env = environ()    
env.io.atom_files_directory = ['.', '../atom_files']
# Read in HETATM records from template PDBs
env.io.hetatm = True

aln = alignment(env)
mdl  = model(env, file='2temp88_fit.pdb')
mdl2 = model(env, file='$FF3')
aln = alignment(env, file='1dsio.ali', align_codes=('$FF2', '$FF3'))

atmsel = selection(mdl).only_atom_types('${OPT}')
r = atmsel.superpose(mdl2, aln,superpose_refine=True,rms_cutoff=6.0)

mdl2.write(file='${SFILE2}_fit.pdb')
EOF
echo "run superimpose FILE2 onto FILE2_subset" >>run.log
python modeller2.in >& modeller2.in.log

#save only interacting residues
NLINE=`awk 'END{print NR}' $FF2`
awk '(NR>'${NLINE}'+1){print $0}' ${SFILE2}_fit.pdb >${SFILE2}_sub.pdb

#rm 1temp88.pdb 2temp88.pdb 1temp88_fit.pdb 2temp88_fit.pdb
#rm modeller.in modeller.in.log
#rm temp5773.ali
