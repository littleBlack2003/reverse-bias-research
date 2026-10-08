#!/usr/bin/env python3
"""NEW formula/import smoke check only. Does not calibrate or solve candidate JV."""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'code'))
import numpy as np
from transport import density_factor, field_factor
from density_device import DensityDevice
from sclc import SCLC
from extend_device import parameters
assert density_factor(0., .074, 300) == 1.
assert field_factor(0., .074, 300, (2.4e20)**(-1/3)) == 1.
assert density_factor(.015, .074, 225) > density_factor(.001, .074, 225) > 1.
try:
    density_factor(.015, .074, 100)
except ValueError:
    pass
else:
    raise AssertionError('100 K candidate guard failed')
d = DensityDevice(parameters(300)[0], N=21, mode='density_only')
o = d.evaluate(d.initial(0.), 0.)
assert np.all(np.isfinite(o['mobility_factor_n']))
assert np.all(np.isfinite(o['mobility_factor_p']))
print('PASS: candidate imports, formula limits, 100 K rejection, finite edge factors')
print('NOT ADOPTED; no candidate device solve or experimental validation performed')
