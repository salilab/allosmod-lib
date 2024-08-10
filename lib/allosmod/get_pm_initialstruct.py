"""Make initial perturbation model (PM) from sequence."""

import shutil
import os
import modeller
import modeller.automodel
import modeller.scripts
import random
import allosmod.util


def get_target(e, target, aln_file):
    """Get the target sequence (either the provided value, or the first entry
       in the alignment if none was provided)"""
    if target:
        return target
    else:
        aln = modeller.Alignment(e, file=aln_file)
        return aln[0].code


def get_pm_initialstruct(aln_file, templates, pdb_dir, nmodel,
                         refine_level, opts):
    # Note the assumption is made in this code that the align code and the
    # PDB file are the same. This is not necessarily the case in Modeller.
    env = modeller.Environ(rand_seed=random.randint(-40000, -2))
    target = get_target(env, opts.target, aln_file)
    dirname = 'pred_%s' % templates[0]
    if opts.chdir:
        if os.path.exists(dirname):
            print("%s exists, overwriting" % dirname)
        else:
            os.mkdir(dirname)
        for template in templates:
            shutil.copy(os.path.join(pdb_dir, template), dirname)
        shutil.copy(aln_file, dirname)
        os.chdir(dirname)

    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    # Read in HETATM records from template PDBs
    env.io.hetatm = True

    class MyModel(modeller.automodel.AutoModel):
        def special_patches(self, aln):
            # If only one chain, set a suitable chain ID (Modeller default
            # is to leave it blank)
            if len(self.chains) == 1:
                chain = aln[self.sequence].range[0].split(':')[1]
                self.chains[0].name = chain if chain else 'A'

    a = MyModel(env, deviation=None, alnfile=aln_file, knowns=templates,
                sequence=target, csrfile=opts.csrfile)
    a.library_schedule = modeller.automodel.autosched.normal
    a.md_level = getattr(modeller.automodel.refine, refine_level)
    a.repeat_optimization = 1
    a.max_molpdf = 1e9

    a.starting_model = 1
    a.ending_model = nmodel
    if not opts.keepaln:
        a.auto_align()  # get an automatic alignment
    a.make(exit_stage=1 if opts.rsronly else 0)


def parse_args():
    usage = """%prog [opts] <alignment file> <template file>
                            <pdb dir> <nmodel> <refine level>

Makes an initial perturbation model (PM) from a sequence.

<alignment file> is a Modeller alignment file.
<template file> is a file listing the templates to use, one per line.
<pdb dir> is the directory containing the template PDB files.
<nmodel> is the number of models to generate.
<refine level> is a Modeller refinement level, e.g. "slow".
"""
    parser = allosmod.util.ModellerOptionParser(usage)
    parser.add_option("--restraints-only", action="store_true", dest="rsronly",
                      help="exit upon making the restraint file "
                           "(make no models)")
    parser.add_option("--target", default=None,
                      help="select the target sequence from the alignment "
                           "file by identifier (by default, the first "
                           "sequence is used)")
    parser.add_option("--keep-alignment", action="store_true", dest="keepaln",
                      help="use the alignment for modeling (by default, "
                           "Modeller will automatically align the target "
                           "with the template)")
    parser.add_option("--no-chdir", action="store_false", dest="chdir",
                      default=True,
                      help="run in the current directory; don't make and "
                           "change into a temporary one")
    parser.add_option("--csrfile", default=None,
                      help="use the given Modeller restraints file")

    opts, args = parser.parse_args()
    if len(args) != 5:
        parser.error("incorrect number of arguments")
    return args[0], args[1], args[2], int(args[3]), args[4], opts


def main():
    aln_file, template_file, pdb_dir, nmodel, refine_level, opts = parse_args()
    get_pm_initialstruct(aln_file, allosmod.util.read_templates(template_file),
                         pdb_dir, nmodel, refine_level, opts)


if __name__ == '__main__':
    main()
