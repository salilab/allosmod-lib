# Dummy config; should not be installed. This is here simply so that the
# test suite works.

import os

datadir = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                       '..', '..', '..', 'data')
local_scratch = '/tmp'
global_scratch = '/scrapp'
