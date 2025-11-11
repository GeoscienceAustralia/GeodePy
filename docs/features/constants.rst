.. _features/constants:

Constants
=========

This module contains constants commonly used in geodetic calculations, including ellipsoids, projections and transformations.

Classes
--------

GeodePy provides four classes for handling commonly used geodetic constants. These include a class for ellipsoids, projections, tranformations and tranformations sigmas.
These classes can be used to create objects that store the relevant parameters for each type of constant, and provide methods for accessing and manipulating these parameters.


.. autoclass:: geodepy.constants.Ellipsoid
    :members:

.. autoclass:: geodepy.constants.Projection
    :members:   

.. autoclass:: geodepy.constants.Transformation
    :members:

    .. automethod:: __neg__
    .. automethod:: __add__

.. autoclass:: geodepy.constants.TransformationSD   
    :members:

Predefined Constants
----------------------

GeodePy also includes a set of predefined constants for commonly used ellipsoids, projections and transformations. These can be seen below.

Ellipsoid
^^^^^^^^^

For :class:`~geodepy.constants.Ellipsoid`:

+-----------------------------+---------------------+-------------------------------------------------------------+
| Name                        | Type                | Description                                                 |
+=============================+=====================+=============================================================+
| grs80                       | Ellipsoid           | Geodetic Reference System 1980                              |
+-----------------------------+---------------------+-------------------------------------------------------------+
| wgs84                       | Ellipsoid           | World Geodetic System 1984                                  |
+-----------------------------+---------------------+-------------------------------------------------------------+
| ans                         | Ellipsoid           | Australian National Spheroid                                |
+-----------------------------+---------------------+-------------------------------------------------------------+
| intl24                      | Ellipsoid           | International (Hayford) 1924                                |
+-----------------------------+---------------------+-------------------------------------------------------------+

Projection
^^^^^^^^^^

For :class:`~geodepy.constants.Projection`:

+-----------------------------+---------------------+-------------------------------------------------------------+
| Name                        | Type                | Description                                                 |
+=============================+=====================+=============================================================+
| utm                         | Projection          | Universal Transverse Mercator projection                    |
+-----------------------------+---------------------+-------------------------------------------------------------+
| isg                         | Projection          | Integrated Survey Grid (NSW, AGD66)                         |
+-----------------------------+---------------------+-------------------------------------------------------------+

Transformation
^^^^^^^^^^^^^^

For :class:`~geodepy.constants.Transformation`:

+-----------------------------+---------------------+-------------------------------------------------------------+
| Name                        | Type                | Description                                                 |
+=============================+=====================+=============================================================+
| gda94_to_gda2020            | Transformation      | GDA94 to GDA2020 (national)                                 |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2014_to_gda2020         | Transformation      | ITRF2014 to GDA2020                                         |
+-----------------------------+---------------------+-------------------------------------------------------------+
| atrf2014_to_gda2020         | Transformation      | ATRF2014 to GDA2020                                         |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2008_to_gda94           | Transformation      | ITRF2008 to GDA94                                           |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2005_to_gda94           | Transformation      | ITRF2005 to GDA94                                           |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2000_to_gda94           | Transformation      | ITRF2000 to GDA94                                           |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf97_to_gda94             | Transformation      | ITRF97 to GDA94                                             |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf96_to_gda94             | Transformation      | ITRF96 to GDA94                                             |
+-----------------------------+---------------------+-------------------------------------------------------------+
| agd84_to_gda94              | Transformation      | AGD84 to GDA94                                              |
+-----------------------------+---------------------+-------------------------------------------------------------+
| agd66_to_gda94              | Transformation      | AGD66 to GDA94                                              |
+-----------------------------+---------------------+-------------------------------------------------------------+
| agd66_to_gda94_act          | Transformation      | AGD66 to GDA94 (ACT region)                                 |
+-----------------------------+---------------------+-------------------------------------------------------------+
| agd66_to_gda94_tas          | Transformation      | AGD66 to GDA94 (Tasmania region)                            |
+-----------------------------+---------------------+-------------------------------------------------------------+
| agd66_to_gda94_vicnsw       | Transformation      | AGD66 to GDA94 (Vic/NSW region)                             |
+-----------------------------+---------------------+-------------------------------------------------------------+
| agd66_to_gda94_nt           | Transformation      | AGD66 to GDA94 (NT region)                                  |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf2014        | Transformation      | ITRF2020 to ITRF2014                                        |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf2014_vel    | Transformation      | ITRF2020 to ITRF2014 (velocity only)                        |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf2008        | Transformation      | ITRF2020 to ITRF2008                                        |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf2005        | Transformation      | ITRF2020 to ITRF2005                                        |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf2000        | Transformation      | ITRF2020 to ITRF2000                                        |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf97          | Transformation      | ITRF2020 to ITRF97                                          |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf96          | Transformation      | ITRF2020 to ITRF96                                          |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf94          | Transformation      | ITRF2020 to ITRF94                                          |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf93          | Transformation      | ITRF2020 to ITRF93                                          |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf92          | Transformation      | ITRF2020 to ITRF92                                          |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf91          | Transformation      | ITRF2020 to ITRF91                                          |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf90          | Transformation      | ITRF2020 to ITRF90                                          |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf89          | Transformation      | ITRF2020 to ITRF89                                          |
+-----------------------------+---------------------+-------------------------------------------------------------+
| itrf2020_to_itrf88          | Transformation      | ITRF2020 to ITRF88                                          |
+-----------------------------+---------------------+-------------------------------------------------------------+


All other combinations of ITRF transofmrations are available.

IERS to GeodePy Transformation
------------------------------

GeeodePy also includes a function for converting from IERS transfomration parameters to a GeodePy Transformation object.

.. autofunction:: geodepy.constants.iers2trans

Height
---------

The module also contains locations for files commonly used in height converstions that will be used in the :ref:`height <features/height>` module.