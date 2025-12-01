.. _tutorials/coord:

Coordinate Classes and Conversions
=============================================

GeodePy has 3 coordinate classes that can be used to store coordinates and convert between coordinate types. 

The three different classes are:

- :ref:`CoordCart - Cartesian Coordinates (x, y, z) <tut/cart>`
- :ref:`CoordGeo - Geographic Coordinates (lat, long, H) <tut/geo>`
- :ref:`CoordTM - Transverse Mercator Coordinates (e, n, H) <tut/tm>`

To learn more about these corrdinate types refer to the `GDA2020 technical manual <https://www.anzlic.gov.au/sites/default/files/files/GDA2020%20Technical%20Manual%20V1.8_published.pdf>`_.

Classes
--------

.. _tut/cart:

Cartesian Coordinates
^^^^^^^^^^^^^^^^^^^^^^
Cartesian coordinates represent points in three dimensions (X, Y, Z), typically in an Earth-Centered, Earth-Fixed (ECEF) system. 
In this class an n value can also be added representing seperation between ellipsoid and geiod.

- **Description:** Defines a point by its distance along three perpendicular axes.
- **Format:** ``(X, Y, Z)`` in meters.
- **Example:** ``( -4052051.0, 4212831.0, -2545100.0 )``

To initalise a cartesian coordinate class:

.. code:: python

    coord1 = geodepy.coord.CoordCart(x, y, z, n=None)

.. _tut/geo:

Geographic Coordinates
^^^^^^^^^^^^^^^^^^^^^^
Geographic coordinates express positions on the Earth's surface using latitude, longitude, and optionally height.

- **Description:** Latitude and longitude define angular position relative to the equator and prime meridian.
- **Format:** ``(latitude, longitude, height)``
- **Example:** ``(-33.8650°, 151.2094°, 58)``

To initalise a geographic coordinate class:

.. code:: python

    coord1 = geodepy.coord.CoordGeo(lat, long, ell_ht=None)

.. _tut/tm:

Transverse Mercator Coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A projected coordinate system that maps the curved Earth onto a flat plane using the Transverse Mercator projection.

- **Description:** Represents positions as Easting and Northing values in meters.
- **Format:** ``(Zone, Easting, Northing, Height)``
- **Example:** ``(55, 334567.89, 6254321.12, 58.2)``

To initalise a transverse mercator coordinate class:

.. code:: python

    coord1 = geodepy.coord.CoordTM(zone, east, north, ell_ht=None)

Converting Between Classes
--------------------------

First import GeodePy

.. code:: python

    import geodepy.coord
    import geodepy.geodesy

We can now create a coordinate obect. For this example we will use a transverse mercator coordinate.

.. code:: python

    coord1 = geodepy.coord.CoordTM(55, 696053.337, 6086610.13)
    print(coord1)

    >>CoordTM: Zone: 55 East: 696053.337 North: 6086610.13 Ell_Ht: None Orth_Ht: None Hemisphere: South

This object can now be transformed into the other classes using the inbuilt methods.

.. code:: python

    print(coord1.cart())
    print(coord1.geo())

    >>CoordCart: X: -4471828.926838844 Y: 2670252.9985762094 Z: -3669113.8962611817 NVal: None
    >>CoordGeo: Lat: -35.3445551951 Lon: 149.15740394128 Ell_Ht: None Orth_Ht: None

Individual variables from a coordinate class can be used within different functions.

.. code:: python

    coord1Geo = coord1.geo()
    print(geodepy.geodesy.rho(coord1Geo.lat)) # to calculate radius of curvature of ellipsoid

    >>6356788.983764104

Using Convert Functions
-----------------------

Instead of using classes, function can also be used to convert between coordinate types. The two main conversions function are:

- geo2grid - Converts from geographic (lat, long, h) to grid (E, N, u)
- xyz2llh - Converts from cartesian (x, y, z) to geographic (lat, long, h)

Both of these functions can be reversed to convert the other way.

To convert using function first geodepy needs to be imported

.. code:: python

    import geodepy.geodesy
    import geodepy.convert

Now either of the functions can be used

.. code:: python

    coordGeo = geodepy.convert.grid2geo(55, 696053.337, 6086610.13)
    print(coordGeo)

    >>(-35.34455523, 149.15740394, 1.0000737, 1.2484390010290551)

    coordllh = geodepy.convert.xyz2llh(-4471828.926838844, 2670252.9985762094, -3669113.8962611817)
    print(coordllh)

    >>(-35.34455523, 149.15740394, 0)

These function can be used together to convert between grid and cartesian.

.. code:: python

    coordGeo = geodepy.convert.grid2geo(55, 696053.337, 6086610.13)
    coordCart = geodepy.convert.llh2xyz(coordGeo[0],coordGeo[1])
    print(coordCart)

    >>(-4471828.924896575, 2670252.9973158445, -3669113.8995236363)