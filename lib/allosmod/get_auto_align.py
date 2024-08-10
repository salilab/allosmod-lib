"""Suggest an alignment of target with templates."""

import argparse
import os
import allosmod.util
from allosmod.pdb2ali import pdb2ali


def get_auto_align(in_aln_file, target, templates, out_aln_file):
    import modeller
    modeller.log.none()
    env = modeller.Environ()
    env.io.atom_files_directory = ['.']
    aln = modeller.Alignment(env, file=in_aln_file, align_codes=target)

    with allosmod.util.temporary_directory() as tempd:
        temp_aln = os.path.join(tempd, "templates.ali")
        with open(temp_aln, 'w') as fh:
            for template in templates:
                pdb2ali(template, fh=fh)
        aln.append(file=temp_aln)
        aln.salign(overhang=30, gap_penalties_1d=(-450, -50),
                   alignment_type='tree', output='ALIGNMENT')
    aln.write(file=out_aln_file)


def parse_args():
    parser = argparse.ArgumentParser(description="""
Suggest an alignment of target with templates.

The sequence for <target> is read from <input align file>, is
combined with the sequences for each template in <templates file>,
and the combination is automatically aligned with SALIGN and written
out as <output align file>.
""")
    parser.add_argument("in_align", help="Input alignment file")
    parser.add_argument("target", help="Target alignment code")
    parser.add_argument("templates",
                        help="File containing a list of templates")
    parser.add_argument("out_align", help="Output alignment file")
    return parser.parse_args()


def main():
    args = parse_args()
    get_auto_align(args.in_align, args.target,
                   allosmod.util.read_templates(args.templates),
                   args.out_align)


if __name__ == '__main__':
    main()
