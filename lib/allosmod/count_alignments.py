"""Count the number of alignments per residue."""

from __future__ import print_function, absolute_import, division
import optparse
import allosmod.util

def count_alignments(aln_file, target, templates):
    import modeller
    modeller.log.none()
    env = modeller.environ()
    aln = modeller.alignment(env, file=aln_file)
    target = aln[target]
    templates = [aln[t] for t in templates]
    num_align = 0
    for r in target.residues:
        for template in templates:
            if r.get_aligned_residue(template) is not None:
                num_align += 1
    return num_align, len(target.residues)

def parse_args():
    usage = """%prog [opts] <alignment file>
                                 <template file> <target>

Count the number of alignments per residue. For every residue in <target>,
count the number of residues it aligns with in <alignment file> in each
template listed in <template file>. Return the total number of such
alignments divided by the number of residues in <target>, rounded to the
nearest whole number.
"""
    parser = optparse.OptionParser(usage)

    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")
    return (args[0], args[1], args[2])

def main():
    aln_file, template_file, target = parse_args()
    num_align, num_res = count_alignments(aln_file, target,
                                allosmod.util.read_templates(template_file))
    print(int(num_align / num_res + 0.5))

if __name__ == '__main__':
    main()
