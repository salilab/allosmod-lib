"""Convert restraints into splines."""

from __future__ import print_function, absolute_import
import allosmod.util


def convert_restraints(rsr):
    from modeller import forms, physical, features
    from allosmod.modeller.forms import TruncatedGaussian

    for group in (physical.ca_distance, physical.n_o_distance,
                  physical.sd_mn_distance, physical.sd_sd_distance,
                  physical.xy_distance, physical.disulfide_distance):
        for form in (TruncatedGaussian, forms.multi_gaussian):
            rsr.spline(form, features.distance, group,
                       spline_dx=0.05, spline_range=4.0,
                       spline_min_points=5, edat=None)


def spline(pdb_file, in_restraints, out_restraints):
    import modeller
    # Needed to keep our custom form alive for restraints.read()
    from allosmod.modeller.forms import TruncatedGaussian  # noqa: F401

    e = modeller.environ()
    m = modeller.model(e, file=pdb_file)
    m.restraints.read(file=in_restraints)
    convert_restraints(m.restraints)
    m.restraints.write(file=out_restraints)


def parse_args():
    usage = """%prog [opts] <pdb file> <restraints in> <restraints out>

Convert restraints into splines.
Selected restraints from the Modeller restraints file <restraints in>
(which apply to <pdb file>) are converted into cubic splines and written
out to <restraints out>.

Currently only the regular Modeller multi-Gaussian and the AllosMod-specific
TruncatedGaussian restraints are converted to splines. Splines are much faster
to evaluate than these restraint types, and it also means that the end user
of an AllosMod protocol need not have the TruncatedGaussian implementation.
"""
    parser = allosmod.util.ModellerOptionParser(usage)

    opts, args = parser.parse_args()
    if len(args) != 3:
        parser.error("incorrect number of arguments")
    return args


def main():
    pdb_file, restraints_in, restraints_out = parse_args()
    spline(pdb_file, restraints_in, restraints_out)


if __name__ == '__main__':
    main()
