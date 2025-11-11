.. _features/ntv2reader:

NTV2 Reader
===============

Tihs module provides functionality to read and utilize NTV2 grid files for coordinate transformations.
It has been adapted from Jaimie Dodd's ntv2reader.py

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

