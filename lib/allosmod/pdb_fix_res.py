"""Modify PDB residues that Modeller cannot recognize."""

import optparse
import fileinput


def pdb_fix_res(pdb_file, inplace):
    for line in fileinput.input(pdb_file, inplace):
        if line.startswith('ATOM') or line.startswith('HETATM'):
            if line[17:20] in ('HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP'):
                line = 'ATOM  ' + line[6:17] + 'HIS' + line[20:]
            if line[17:20] == 'MSE':
                line = 'ATOM  ' + line[6:17] + 'MET' + line[20:]
                if line[12:14] == 'SE':
                    line = line[:12] + ' SD' + line[15:]
            print(line, end='')


def parse_args():
    usage = """%prog <PDB file>

Modify PDB residues that Modeller cannot recognize, and write the new PDB file
to standard output (or modify it inplace if --in-place is specified).

The following substitutions are made:
- Everything except ATOM and HETATM records is stripped from the file.
- MSE residues (either ATOM or HETATM) are mapped to HIS ATOM records,
  and any SE atoms in those residues are mapped to SD.
"""
    parser = optparse.OptionParser(usage)
    parser.add_option("--in-place", action="store_true", dest="inplace",
                      help="modify the file in place")

    opts, args = parser.parse_args()
    if len(args) != 1:
        parser.error("incorrect number of arguments")
    return args[0], opts


def main():
    pdb_file, opts = parse_args()
    pdb_fix_res(pdb_file, opts.inplace)


if __name__ == '__main__':
    main()
