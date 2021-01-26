"""Make input files for Modeller."""

from __future__ import print_function, absolute_import
import allosmod.util


def make_mod_inputs(target, templates, rand_seed, rand_ang, deviation,
                    nucleic_acids):
    import modeller.automodel
    import modeller.scripts
    env = modeller.Environ(rand_seed=rand_seed)
    env.io.atom_files_directory = ['../atom_files']
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    # Read in HETATM records from template PDBs
    env.io.hetatm = True

    mdl = modeller.scripts.complete_pdb(env, 'avgpdb.pdb')

    # Act on all atoms in the model
    sel = modeller.Selection(mdl)

    # Change all existing X,Y,Z for +- dev:
    sel.rotate_mass_center([0, 0, 1], rand_ang[0])
    sel.rotate_mass_center([0, 1, 0], rand_ang[1])
    sel.rotate_mass_center([1, 0, 0], rand_ang[2])
    sel.randomize_xyz(deviation=deviation)
    mdl.write(file='random.ini')
    a = modeller.automodel.AutoModel(env, deviation=None, alnfile='align.ali',
                                     knowns=templates, sequence=target)
    if nucleic_acids:
        # use longer restraints for nucleic acids
        a.max_sc_sc_distance = 14.0

    # Would perhaps make more sense to set this to False, so we can edit the
    # restraint later on (and restraints get splined in the last step anyway)
    # but the original code has it True
    a.spline_on_site = True
    a.make(exit_stage=1)


def parse_args():
    usage = """%prog [opts] <target> <template file>
                            <rand seed> <dz> <dy> <dx> <deviation>

Make input files for Modeller. Given an alignment file (align.ali) and
a structure file (avgpdb.pdb), it will generate a randomized initial
model (random.ini) and a restraints file.

<target> is the target sequence.
<template file> is a file listing the templates to use, one per line.
<rand seed> is a number to use as the Modeller random number seed.
<dz>, <dy>, <dx> are angles (in degrees) to rotate the initial model by
                 about each axis (z, y, x).
<deviation> is an amount (in angstroms) to randomize the model by.
"""
    parser = allosmod.util.ModellerOptionParser(usage)
    parser.add_option("--nucleic-acids", action="store_true", dest="nuc",
                      help="use longer restraints for nucleic acids")

    opts, args = parser.parse_args()
    if len(args) != 7:
        parser.error("incorrect number of arguments")
    return (args[0], args[1], int(args[2]),
            (float(args[3]), float(args[4]), float(args[5])),
            float(args[6]), opts)


def main():
    target, template_file, rand_seed, rand_ang, deviation, opts = parse_args()
    make_mod_inputs(target, allosmod.util.read_templates(template_file),
                    rand_seed, rand_ang, deviation, opts.nuc)


if __name__ == '__main__':
    main()
