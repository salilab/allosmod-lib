"""Generate Modeller script to model with glycosylation"""

from __future__ import print_function, absolute_import
import optparse
import subprocess
import allosmod.util

def get_pm_glyc(target, templates, rand, rep_opt, att_gap, glycpm):
    # No Python implementation yet - call the original script:
    script = allosmod.util.get_data_file('get_pm_glyc.sh')
    subprocess.check_call([script, target, templates, rand, rep_opt,
                           att_gap, glycpm])

def parse_args():
    usage = """%%prog [opts] <target> <templates> <rand> <rep_opt>
                 <att_gap> <glycpm>

Generate Modeller script to model with glycosylation.

<target> alignment code of the target (usually pm.pdb)
<templates> file containing a list of all templates
<rand> random seed
<rep_opt> option to allow repeat optimization
<att_gap> if set, inserts gaps to allow flexibility at glycosylation sites
<glycpm> either "script" or the name of a PDB file
"""
    parser = optparse.OptionParser(usage)

    opts, args = parser.parse_args()
    if len(args) != 6:
        parser.error("incorrect number of arguments")
    return args

def main():
    target, templates, rand, rep_opt, att_gap, glycpm = parse_args()
    get_pm_glyc(target, templates, rand, rep_opt, att_gap, glycpm)

if __name__ == '__main__':
    main()
