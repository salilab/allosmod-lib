"""Suggest an alignment of target with templates."""

from __future__ import print_function, absolute_import
import optparse
import os
import allosmod.util
from allosmod.pdb2ali import pdb2ali


def get_auto_align(in_aln_file, target, templates, out_aln_file):
    import modeller
    modeller.log.none()
    env = modeller.environ()
    env.io.atom_files_directory = ['.']
    aln = modeller.alignment(env, file=in_aln_file, align_codes=target)

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
    usage = """%prog <input align file> <target> <templates file>
                  <output align file>

Suggest an alignment of target with templates.

The sequence for <target> is read from <input align file>, is
combined with the sequences for each template in <templates file>,
and the combination is automatically aligned with SALIGN and written
out as <output align file>.
"""
    parser = optparse.OptionParser(usage)
    options, args = parser.parse_args()
    if len(args) != 4:
        parser.error("incorrect number of arguments")
    return args


def main():
    in_aln_file, target, templates_file, out_aln_file = parse_args()
    get_auto_align(in_aln_file, target,
                   allosmod.util.read_templates(templates_file), out_aln_file)


if __name__ == '__main__':
    main()
