Basic usage
***********

Once installed, the AllosMod protocol can be run by means of a command line
tool (``allosmod``). Each component of the protocol is also a Python package,
which can be called directly from other Python software
(via ``import allosmod``).

Overview
========

Running the basic protocol consists of three steps:

#. Create a set of input files specifying the sequence to model, any
   structures to use as templates, and AllosMod parameters.

#. Run ``allosmod setup`` to verify the inputs and generate a script file.

#. Run the script file (typically on a Linux cluster) to set up the AllosMod
   energy landscape and generate Modeller input files to sample it. In some
   cases this sampling is then carried out; in others you will need to run
   Modeller on the input files.

Input files
===========

Several input files are needed to run the protocol:

PDB files
---------

All structures used to create the energy landscape for the simulation should
be provided in PDB format.

Alignment file
--------------

This should be named ``align.ali`` and should contain one entry for each
PDB file above (the align code, i.e. the part after the ``>P1;`` header,
should match the
PDB filename) and another entry (named ``pm.pdb``) with the sequence to
be simulated. This alignment file should be generated after an alignment
procedure, as this alignment will be used to generate restraints for
the simulation. Multiple chains can be specified by using a "/" as a
separator, the same specifications used in MODELLER. There are many ways
to create an alignment file including `MODELLER <https://salilab.org/modeller/tutorial/basic.html>`_ and `ClustalW <https://www.ebi.ac.uk/Tools/msa/clustalw2/>`_.

*WARNING* Small errors in the alignment can cause big errors during a
simulation due to energy conservation problems. Make sure there are no
misalignments in which adjacent residues are aligned far apart in
sequence (alignment programs often do this at the beginning or end
of chains).

Structure list
--------------

All the PDB files used to create the energy landscape for the simulation
should be listed in a file called ``list``, one per line.
Refer to the LIGPDB and ASPDB options in ``input.dat`` to
define interactions in the allosteric site.

Ligand file
-----------

If desired, a ligand file (called ``lig.pdb``) can be provided; this contains
the structure of the ligand extracted from a ligand bound PDB file
(defined by LIGPDB in ``input.dat``). A radius (rAS) around the ligand
is used to define the allosteric site. If ``lig.pdb`` is excluded,
AllosMod will set up a landscape with as many energy minima as are
described by structures in the list file.

AllosMod parameter file
-----------------------

This should be named ``input.dat`` and should contain
one line per parameter as follows. All parameters are optional except for
NRUNS.

NRUNS = X
    is the number of independent simulations to run.

rAS = X
    is the radius (in Ångstroms) around the coordinates of the ligand that
    will specify the allosteric site. If the file ``lig.pdb`` is included,
    the allosteric site will be calculated using rAS and the coordinates
    in LIGPDB. Therefore, ``lig.pdb`` must be extracted from LIGPDB.

SAMPLING = X
    can be one of:

    ``simulation`` (Default)
        A simulation is set up to be run later.

    ``moderate_am``
        Sampling is performed using a quick, unequilibrated simulation.
        This quick sampling will give a representation of the types of
        conformations that are consistent with the modeled energy landscape.
        Set "SAMPLING = simulation" to predict the relative populations
        of the conformations at equilibrium.

delEmax = X
    is the maximum energy for each pairwise atomic distance contact,
    typically between 0.09 and 0.12 kcal/mol. If not given or set to the
    special value "CALC", the value will be assigned according to
    3.6*(number of residues/number of distance interactions).
    See paper for more details.

LIGPDB = X
    is the PDB file used to define the allosteric site. AllosMod defines
    the allosteric site using the distance (rAS) from the effector
    (``lig.pdb``) with respect to the LIGPDB coordinates.

ASPDB = X
    is the PDB file used to define the contacts in the allosteric site,
    i.e. the pairwise atomic distances in ASPDB are used to determine the
    nonbonded distance energy. As an example, to run an effector unbound
    simulation: 1) include the effector bound and unbound PDB files in
    ``align.ali`` and ``list``, 2) set ASPDB to the effector unbound
    PDB file, and 3) set LIGPDB to the effector bound PDB file.

DEVIATION = X
    is the distance (in Ångstroms) that the atoms will be randomized
    when creating the initial structure (default is 1-10 Å depending on
    simulation type).

MDTEMP = X
    is the temperature (in degrees Kelvin) for the simulation (default
    is 300 K). Alternatively, set MDTEMP to "scan" and the simulation
    temperature will alternate between 300 K, 350 K, 400 K, 450 K, and
    500 K. Therefore, directory 0 will have a 300 K simulation, directory 1
    will have a 350 K simulation, and so on until directory 5 that will
    restart the sequence with a 300 K simulation.

BREAK = True/False
    is an option to include chemical frustration (Weinkam et al. 2009
    Biochemistry, p2394-2402). Chemical frustration is modeled by breaking
    all interactions involving buried, charged residues. Regions with
    many buried, charged residues will have high conformational variability.

SCLBREAK = X
    if BREAK=True, this number is used to scale the contacts with residues
    that cause chemical frustration.

CHEMFR = cdensity/charge
    if BREAK=True, this selects the type of chemical frustration to use.
    If set to 'cdensity' (the default) then a distribution of charged contacts per
    residue is calculated; all residues with a z-score above ``ZCUTOFF``
    (see below) are predicted to cause chemical frustration. If set to 'charge'
    then all residues with a certain number of charged contacts are used.
    
ZCUTOFF = X
    if BREAK=True and CHEMFR=cdensity, this number is used to select which
    residues cause chemical frustration. ZCUTOFF is the z-score cutoff of the
    distribution involving the number of charged contacts per residue; residues with
    a z-score above this threshold are predicted to cause chemical frustration.

LOCALRIGID = True/False
    if set to True, secondary structure, corresponding to the input PDB
    files, will have increased stability in the simulation. Increased
    stability is maintained by increasing the energy by a factor of 10
    for all Cα-Cα contacts between 2 and 5 residues apart.

COARSE = True/False
    is an option to coarse grain the energy landscape by restricting the
    nonbonded distance energy to include Cα and Cβ atoms only.
    This allows very large proteins to be simulated without overwhelming
    the computer's memory. This option is automatically set to True for
    proteins over 1500 residues.

{ADDITIONAL_RESTRAINT} {DISTANCE} {STANDARD_DEVIATION} {INDICES}
    is used to add additional restraints between residues.
    ADDITIONAL_RESTRAINT can be HARM, LOBD, or UPBD corresponding to
    distance restraints that are harmonic, lower bounded only, or upper
    bounded only, respectively. DISTANCE and STANDARD_DEVIATION corresponds
    to the distance (in Ångstroms) between two atoms in the residues
    specified in INDICES. If residue index is an amino acid, atom type
    will be CA, otherwise atom type will be the first present: N, P, C,
    or O. INDICES is a list of residue indices separated by commas.
    Restraints are added between each successive pair of indices,
    i.e. between i1 and i2, between i3 and i4, ... The residue index
    corresponds to the position in the input alignment file. Therefore,
    if there are multiple chains, the index for the first residue in the
    second chain will be one more than the index for the last residue in
    the first chain (refer to any output PDB for simplicity).

Alter residue contact energies
------------------------------

If desired, a file ``break.dat`` can be provided, which contains
a list of residues whose pairwise contact energies (delEmax) will be
scaled by a specified value. Each line contains one residue index
(corresponding to simulated sequence) in the first column and one
scaling factor in the second column. For example, to reduce all
contact energies for residue 30 by 90 %, ``break.dat`` would have one
line with "30 0.1". ``break.dat`` is created automatically by setting
BREAK=True, however, the user may specify any desired residues and
scaling factors by including ``break.dat`` in a batch run.

Set up AllosMod protocol
========================

Once all the input files are prepared, run ``allosmod setup`` in the directory
containing them. The ``allosmod`` command line tool provides many subfunctions
(use ``allosmod help`` to list them all). ``allosmod setup`` will check the
input files for problems, and if they all look OK, it will generate a
script file called ``qsub.sh``. This script can be run on any Linux machine,
although it is intended to be run on an SGE cluster using something like
``qsub -S /bin/sh -l arch=linux-x64 -cwd -t 1-N qsub.sh``,
where ``N`` is the value of NRUNS in ``input.dat``.

This script file will set up the AllosMod landscape. If SAMPLING in
``input.dat`` is set to 'simulation' (the default) MODELLER input files are
generated. These can then be run to perform the simulation. Otherwise, the
sampling is performed by ``qsub.sh`` itself.
