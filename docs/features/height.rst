.. _features/height:  

Vertical Datums
=================

This module includes functions for dealing with height datums and geoids.

.. note:: GDAL is needed for this module. Please consult the GDAL pypi for instructions on installing GDAL

Converting Between Height Datums
--------------------------------

The following functions all convert between different height datums with a particular focus on Australian datums.

.. autofunction:: geodepy.height.GPS_to_AHD
.. autofunction:: geodepy.height.AHD_to_GPS
.. autofunction:: geodepy.height.GPS_to_AVWS
.. autofunction:: geodepy.height.AVWS_to_GPS
.. autofunction:: geodepy.height.AHD_to_AVWS
.. autofunction:: geodepy.height.AVWS_to_AHD
.. autofunction:: geodepy.height.GPS_to_AUSGeoid98
.. autofunction:: geodepy.height.AUSGeoid98_to_GPS
.. autofunction:: geodepy.height.GPS_to_AUSGeoid09
.. autofunction:: geodepy.height.AUSGeoid09_to_GPS
.. autofunction:: geodepy.height.DOV
.. autofunction:: geodepy.height.DOV_09
.. autofunction:: geodepy.height.DOV_98

Calculating Gravity
-------------------

The functions below can be used to calculate different components of gravity.

.. autofunction:: geodepy.height.mean_surface_grav
.. autofunction:: geodepy.height.interp_grav
.. autofunction:: geodepy.height.normal_grav
.. autofunction:: geodepy.height.mean_normal_grav
.. autofunction:: geodepy.height.normal_correction
.. autofunction:: geodepy.height.normal_orthometric_correction

Auxilary Function
-----------------

Functions used to enable other functions in the module.

.. autofunction:: geodepy.height.interp_file