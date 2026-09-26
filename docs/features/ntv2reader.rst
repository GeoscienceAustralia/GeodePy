.. _features/ntv2reader:

NTV2 Reader
===============

Tihs module provides functionality to read and utilize NTV2 grid files for coordinate transformations.
It has been adapted from Jaimie Dodd's ntv2reader.py

Interpolation near subgrid boundaries
------------------------------------

Bicubic interpolation requires a 4-by-4 neighbourhood of nodes entirely within
the selected subgrid. It raises ``ValueError`` where this neighbourhood is
unavailable, including cells along the subgrid edges and subgrids with fewer
than four rows or columns. This prevents reads from neighbouring rows,
subgrid headers, or beyond the file.

For these locations, explicitly select ``method="bilinear"`` when calling
``interpolate_ntv2`` or ``geodepy.transform.ntv2_2d``. The interpolation method
is not changed automatically. Subgrid selection and the existing exclusion
of coordinates on the northern and western extents are unchanged.

Classes
-------

.. autoclass:: geodepy.ntv2reader.NTv2Grid
   :members:

.. autoclass:: geodepy.ntv2reader.SubGrid
   :members:

Functions
---------

.. autofunction:: geodepy.ntv2reader.read_ntv2_file
.. autofunction:: geodepy.ntv2reader.interpolate_ntv2
.. autofunction:: geodepy.ntv2reader.read_node
.. autofunction:: geodepy.ntv2reader.bilinear_interpolation
.. autofunction:: geodepy.ntv2reader.bicubic_interpolation

