from __future__ import print_function, absolute_import

import modeller
from allosmod.modeller import _truncated_gaussian


class TruncatedGaussian(modeller.forms.restraint_form):
    """AllosMod truncated Gaussian restraint.
       This is implemented as a C extension (_truncated_gaussian.so)
       to Modeller."""

    _builtin_index = _truncated_gaussian.truncated_gaussian_create()

    def __init__(self, group, feature, dele_max, slope, scl_delx, weights,
                 means, stdevs):
        lv = -1
        for var in (weights, means, stdevs):
            if (lv >= 0 and lv != len(var)) \
               or not isinstance(var, (tuple, list)):
                raise TypeError("weights, means and stdevs should all be "
                                "sequences of the same length")
            lv = len(var)
        modeller.forms.restraint_form.__init__(
            self, group, feature, len(weights),
            (dele_max, slope, scl_delx) + tuple(weights) + tuple(means)
            + tuple(stdevs))
