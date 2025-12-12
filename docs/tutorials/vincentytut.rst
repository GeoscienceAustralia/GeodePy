.. _tutorials/vincenty:

Connecting Two Points
=====================

GeodePy can calculate the relationship between two points in two different ways, either on a flat plane or curved plane.
To learn more about geodetic formulas refer to the `GDA2020 technical manual <https://www.anzlic.gov.au/sites/default/files/files/GDA2020%20Technical%20Manual%20V1.8_published.pdf>`_.

Flat plane
----------

To calculate the distance and bearing between two points the following command can be used. Here the coordinates of two 
points are entered.

.. code:: python

    import geodepy.survey
    import geodepy.angles

    connection = geodepy.survey.joins(
        696053.337, #easting of A
        6086610.13, #northing of A
        696770.781, #easting of B
        6086089.772 #northing of B
    )
    bearing = connection[1]
    bearing = geodepy.angles.dec2dms(bearing)

    print(f"The connection is: {connection[0]:.4f} @ {bearing}")

    >>The connection is: 886.2834 @ 125 57 11.37834570

Here the first number is the horizontal distance while the second is the bearing, here converted to DMS.

You can also determine the coordinate of a second point given the first coordinate and a bearing and distance. 

.. code:: python

    coord = geodepy.survey.radiations(
        696053.337, #easting of A
        6086610.13, #northing of A
        bearing.dec(), #bearing A to B
        886.2834 #distance A to B
    )
    print(f"The point 2 coord is: {coord[0]:.3f} {coord[1]:.3f}")

    >>The point 2 coord is: 696770.781 6086089.772

As can be seen this value matches the coordinate of the second point entered in the joins function.

Curved Plane
------------

To calculate the joins between two points on a curved plane, vincenty's formula can be used. Within GeodePy there 
are two types of Vincenty's that either use geographic coordinates or grid coordinates. In this example we will 
use the utm vincenty functions that uses grid coordinates but this same process can be undertaken using latitude and longitude.

.. code:: python 

    import geodepy.geodesy
    import geodepy.angles

    connection2 = geodepy.geodesy.vincinv_utm(
        55, #zone of A
        696053.337, #easting of A
        6086610.13, #northing of A
        55, #zone of B
        696770.781, #easting of B
        6086089.772 #northing of B
    )
    bearing2 = connection2[1]
    bearing2 = geodepy.angles.dec2dms(bearing2)

    print(f"The connection is: {connection2[0]:.4f} @ {bearing2}")

    >>The connection is: 886.2839 @ 125 57 11.1186109

As seen here the connection is slightly different between the flat and curved plane.

The coordinate of a second point can also be calculated.

.. code:: python

    coord2 = geodepy.geodesy.vincdir_utm(
        55, #zone of A
        696053.337, #easting of A
        6086610.13, #northing of A
        bearing2.dec(), #bearing A to B
        886.2839 #distance A to B
    )
    print(f"The point 2 coord is: {coord2[1]:.3f} {coord2[2]:.3f}")

    >>The point 2 coord is: 696770.781 6086089.772

Again it can be seen that this coordinate matches what was entered earlier.