"""wfn: wind farm noise.

Sub-packages are independent and only share the input file format and
``xarray`` data structures:

* :mod:`wfn.source`: wind turbine source model (TEN, TIN)
* :mod:`wfn.flow`: flow models (interface + simple analytical flows)
* :mod:`wfn.propagation`: WAPE propagation (setup, local runs, deltaL assembly)
* :mod:`wfn.receivers`: receiver grids
* :mod:`wfn.units`: dB helpers
* :mod:`wfn.spl`: source + propagation coupling (:func:`wfn.spl.combine`)
"""
import logging

from . import receivers, spl, units

logging.getLogger(__name__).addHandler(logging.NullHandler())

__version__ = "0.1.0"
__all__ = ["receivers", "spl", "units"]
