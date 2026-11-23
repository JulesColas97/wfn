"""WAPE propagation model (next step).

Planned content:

* ``config.py``: ``propagation`` section of the input file
* ``inputs.py``: write the Fortran namelist ``input.nml`` for each (source height, direction)
* ``launcher.py``: run the WAPE executable locally or through a scheduler
* ``results.py``: read WAPE outputs into an ``xr.Dataset`` of ``deltaL`` with dims
  ``(*grid_dims, freq, source_height)``

Coupling with the source is then a function adding ``deltaL`` interpolated at
the blade segment heights ``source_ds.z_seg`` to the source SPL.
"""
