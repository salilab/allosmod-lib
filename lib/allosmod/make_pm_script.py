"""Make a script to generate perturbation models (PM)."""

from __future__ import print_function, absolute_import
import optparse
import allosmod.util

def make_pm_script(target, templates, rand_seed, deviation, typ, mdtemp):
    template = {'script':'pm_script.py.in',
                'moderate_cm':'pm_moderate_cm.py.in',
                'moderate_cm_simulation':'pm_moderate_cm.py.in',
                'fast_cm':'pm_fast_cm.py.in',
                'fast_cm_simulation':'pm_fast_cm.py.in',
                'moderate_am':'pm_moderate_am.py.in'}[typ]
    subs = {'RAND':str(rand_seed), 'DEVIATION':str(deviation),
            'MDTEMP':str(mdtemp), 'KNOWNS':repr(templates),
            'TARG_SEQ':repr(target)}
    allosmod.util.subst_file(open(allosmod.util.get_data_file(template)),
                             open('model_run.py', 'w'), subs)

def parse_args():
    script_types = ['script', 'moderate_cm', 'moderate_cm_simulation',
                    'fast_cm', 'fast_cm_simulation', 'moderate_am']
    usage = """%%prog [opts] <target> <template file>
                            <rand seed> <deviation> <type> <mdtemp>

Make a script file to generate perturbation models (PM).

<target> is the target sequence.
<template file> is a file listing the templates to use, one per line.
<rand seed> is a number to use as the Modeller random number seed.
<deviation> is an amount (in angstroms) to randomize the model by.
<type> is the type of script to generate; can be %s.
<mdtemp> is the temperature for MD simulations, in Kelvin.
""" % ", ".join(script_types)
    parser = optparse.OptionParser(usage)

    opts, args = parser.parse_args()
    if len(args) != 6:
        parser.error("incorrect number of arguments")
    if args[4].lower() not in script_types:
        parser.error("type should be one of %s" % ", ".join(script_types))
    return (args[0], args[1], args[2], float(args[3]), args[4].lower(),
            float(args[5]))

def main():
    target, template_file, rand_seed, deviation, typ, mdtemp = parse_args()
    make_pm_script(target, allosmod.util.read_templates(template_file),
                   rand_seed, deviation, typ, mdtemp)

if __name__ == '__main__':
    main()
