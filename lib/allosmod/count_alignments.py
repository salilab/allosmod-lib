"""Count the number of alignments per residue."""

import argparse
import allosmod.util


def count_alignments(aln_file, target, templates):
    import modeller
    modeller.log.none()
    env = modeller.Environ()
    aln = modeller.Alignment(env, file=aln_file)
    target = aln[target]
    templates = [aln[t] for t in templates]
    num_align = 0
    for r in target.residues:
        for template in templates:
            if r.get_aligned_residue(template) is not None:
                num_align += 1
    return num_align, len(target.residues)


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
Count the number of alignments per residue. For every residue in <target>,
count the number of residues it aligns with in <alignment file> in each
template listed in <template file>. Return the total number of such
alignments divided by the number of residues in <target>, rounded to the
nearest whole number.""")
    parser.add_argument("alignment", help="alignment file")
    parser.add_argument("template", help="template file")
    parser.add_argument("target", help="target align code")

    return parser.parse_args()


def main():
    args = parse_args()
    num_align, num_res = count_alignments(
        args.alignment, args.target,
        allosmod.util.read_templates(args.template))
    print(int(num_align / num_res + 0.5))


if __name__ == '__main__':
    main()
