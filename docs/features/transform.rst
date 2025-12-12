.. _features/transform:

Transforming
=================

This module provides functions for transforming coordinates between different geodetic datums. It includes methods for performing Helmert transformations, functions for various 
common geodetic transformations and support for NTv2 2D grid-based transformations.

Helmert Transformations
-----------------------

.. autofunction:: geodepy.transform.conform7

.. autofunction:: geodepy.transform.conform14

Common Geodetic Transformations
--------------------------------

.. autofunction:: geodepy.transform.transform_mga94_to_mga2020
.. autofunction:: geodepy.transform.transform_mga2020_to_mga94
.. autofunction:: geodepy.transform.transform_atrf2014_to_gda2020
.. autofunction:: geodepy.transform.transform_gda2020_to_atrf2014

NTv2 Transformations
----------------------

.. autofunction:: geodepy.transform.ntv2_2d
