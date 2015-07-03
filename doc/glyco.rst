Modeling protein glycosylation
******************************

AllosMod can also be used to model protein glycosylation. For full details,
see `Guttman et al. 2013 Structure <http://www.ncbi.nlm.nih.gov/pubmed/23473666>`_.
This tool is meant to sample the range of conformations accessible to sugars.
There are a few ways to run the protocol depending on the desired
conformational flexibility of the protein and/or sugars:

 #. Input: a protein without sugars. Output: a semi-rigid protein (motions
    will be limited to loops and surface side chains) with flexible sugars.

 #. Input: a protein with sugars. Output: a flexible protein with rigid sugars.

 #. Input: a protein without sugars. Output: a flexible protein with flexible
    sugars.

For option 1, input the protein sequence only and add a file describing sugar
connectivity (``glyc.dat``). See below for more information.

For option 2, run option 1 and then use a single structure from the output as
input for another run. In the input sequence, specify the protein sequence
and a "." for each sugar. The "." specifies a block residue in MODELLER
and will remain rigid. Such rigid bodies can cause problems during dynamics
so care should be given here. Also include the ``allosmod.py`` file that was
generated in the first step.

For option 3, first run simulations of the protein only using AllosMod.
Then take the resulting structures and do option 1 (one for each structure)
using AllosMod. Take all the glycosylated protein structures as input into
the FoXS server, if desired.

Glycosylation input files
=========================

Glycosylation is controlled by adding a file ``glyc.dat`` to the files used
in the regular AllosMod protocol. For help making the ``glyc.dat`` file from
a PDB file that includes sugars, see the GlycoSciences
`webpage <http://www.glycosciences.de/tools/pdbcare/>`_ and set the
"Assign connections by atom distances" to "HETATOM only." This tool will
also assess the structural quality of the sugars. 

Parameter file
--------------

A number of options in the AllosMod parameter file (``input.dat``) pertain
specifically to glycosylation. They are all optional.

ATTACH_GAPS = True/False (default True)
    Gaps will be inserted into the alignment file at the sugar attachment
    sites. This allows more flexibility of the protein to accommodate the
    sugars.

SAMPLING = X
    defaults to ``moderate_cm``: Sampling is performed using simulated
    annealing at moderate temperatures. The energy function for the protein
    is composed of restraints from MODELLER while the energy function for
    sugar molecules is a combination of dihedral terms from CHARMM and
    harmonic restraints as described in the published paper
    (Guttman et al. 2013 Structure). 

REPEAT_OPTIMIZATION = X
    is the number of optimization steps to be performed. Default is one. 

Glycosylation file
------------------

The glycosylation file (which must be called ``glyc.dat``) contains one line
per sugar monomer grouped by chains of monomers. The line with the
protein-bonded monomer specifies the beginning of a chain and is followed
by all monomers within the chain. Three columns define the sugar types
and connectivity: 1) monomer name, 2) O1 bond type, and
3) O1 attachment residue index:

 * monomer name:
 
   * NAG - b-N-Acetyl-D-Glucosamine
   * NGA - b-N-Acetyl-D-Galactosamine
   * GLB - b-Galactose
   * FUC - a-Fucose
   * MAN - a-Mannose
   * BMA - b-Mannose
   * NAN - a-Neuraminic acid

 * O1 bond type:
 
   * NGLA or NGLB - axial or equatorial bond to ASN*
   * SGPA or SGPB - bond between alpha or beta position and SER*
   * TGPA or TGPB - bond between alpha or beta position and THR*
   * 16ab - (i) 1->6 (i-1) axial at C1 and equatorial at C6
   * 16fu - (i) 1->6 (i-1) axial at C1 and equatorial at C6
   * 14bb - (i) 1->4 (i-1) equatorial at C1 and equatorial at C4
   * 13ab - (i) 1->3 (i-1) axial at C1 and equatorial at C3
   * 13bb - (i) 1->4 (i-1) equatorial at C1 and equatorial at C3
   * 12aa - (i) 1->2 (i-1) axial at C1 and axial at C2
   * 12ba - (i) 1->2 (i-1) equatorial at C1 and axial at C2
   * sa23 - og sialic acid alpha 2->3 equatorial
   * sa26 - og sialic acid alpha 2->6 equatorial

        \*entry specifies the beginning of a new chain 

 * O1 attachment residue index:
 
    For every new sugar chain, specify the amino acid index of the
    residue bound to the first O1 atom in the sugar. The amino acid
    index must correspond to the input sequence in the alignment,
    i.e. the first residue in the input sequence is 1 and the
    Nth is N. For sugars within the chain, specify the index of
    the sugar monomer bound to the O1 atom defined as such: the
    first listed monomer in a chain has index 1, the second listed
    monomer in a chain has index 2,...

Output
======

On output a directory is created for each run, each containing a glycosylated
structure named ``pm.pdb.B99990001.pdb``.
