"""Structurally align two PDBs using SALIGN."""

from __future__ import print_function, absolute_import
import optparse

def determine_fit_atoms(mdl):
    n_ca = n_p = 0
    for atom in mdl.atoms:
        if atom.name == 'CA':
            n_ca += 1
        elif atom.name == 'P':
            n_p += 1
    if n_ca >= n_p:
        return 'CA'
    else:
        return 'P'

def salign0(ff1, ff2):
    import modeller
    modeller.log.verbose()
    env = modeller.environ()
    env.io.atom_files_directory = ['.', '../atom_files']

    # Read in HETATM records from template PDBs
    env.io.hetatm = True

    aln = modeller.alignment(env)
    code = ff1
    mdl = modeller.model(env, file=code, model_segment=('FIRST:@', 'END:'))
    fit_atoms = determine_fit_atoms(mdl)
    aln.append_model(mdl, atom_files=code, align_codes=code)
    code = ff2
    mdl = modeller.model(env, file=code, model_segment=('FIRST:@', 'END:'))
    aln.append_model(mdl, atom_files=code, align_codes=code)

    for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                        ((1.,0.5, 1., 1., 1., 0.), False, True),
                                        ((1.,1., 1., 1., 1., 0.), True, False)):
        aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                   rr_file='$(LIB)/as1.sim.mat', overhang=30,
                   gap_penalties_1d=(-450, -50),
                   gap_penalties_3d=(0, 3), gap_gap_score=0,
                   gap_residue_score=0, fit_atoms=fit_atoms,
                   alignment_type='tree', feature_weights=weights,
                   improve_alignment=True, fit=True, write_fit=write_fit,
                   write_whole_pdb=whole, output='ALIGNMENT QUALITY')
    return aln

def parse_args():
    usage = """%prog <PDB file> <PDB file> <alignment file>

Structurally align two PDBs using SALIGN and write out an alignment file
and fit PDBs.

The fit is done using either P or CA atoms depending on which are found more
in the first PDB file.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")
    return args[0], args[1], args[2]

def main():
    ff1, ff2, aln_file = parse_args()
    aln = salign0(ff1, ff2)
    aln.write(file=aln_file)

if __name__ == '__main__':
    main()
