.. _tutorials/vincenty:

Connecting Two Points
=====================

GeodePy has the ability to calculate the relationship between two points in two different ways, either on a flat plane or curved plane.

Flat plane
----------

To calculate the distance and bearing between two points the following command can be used.

.. code:: python

    import geodepy.survey
    import geodepy.angles

    connection = geodepy.survey.joins(696053.337, 6086610.13, 696770.781, 6086089.772)
    bearing = connection[1]
    bearing = geodepy.angles.dec2dms(bearing)

    print(f"The connection is: {connection[0]:.4f} @ {bearing}")

    >>The connection is: 886.2834 @ 125 57 11.37834570

Here the first number is distance while the second is the bearing. It might look clearer when converting the bearng to DMS

You can also determine the coordinate of a second point given first coordinate and a bearing and distance. 

.. code:: python

    coord = geodepy.survey.radiations(696053.337, 6086610.13, bearing.dec(), 886.2834)
    print(f"The point 2 coord is: {coord[0]:.3f} {coord[1]:.3f}")

    >>The point 2 coord is: 696770.781 6086089.772

As can be seen this value matches the coordinate of the second point entered in the joins function.

Curved Plane
------------

To calculate the joins between two points on a curved plan vincenty's formula can be used. In this example we will 
use the utm vincenty functions but tihs same process can be undertaken using the normal formulas and lat and longs.

.. code:: python 

    import geodepy.geodesy
    import geodepy.angles

    connection2 = geodepy.geodesy.vincinv_utm(55, 696053.337, 6086610.13, 55, 696770.781, 6086089.772)
    bearing2 = connection2[1]
    bearing2 = geodepy.angles.dec2dms(bearing2)

    print(f"The connection is: {connection2[0]:.4f} @ {bearing2}")

    >>The connection is: 886.2839 @ 125 57 11.1186109

As seen here the connection is slightly different between the flat and curved plane.

The coordinate of a second point can also be calculated.

.. code:: python

    coord2 = geodepy.geodesy.vincdir_utm(55, 696053.337, 6086610.13, bearing2.dec(), 886.2839)
    print(f"The point 2 coord is: {coord2[1]:.3f} {coord2[2]:.3f}")

    >>The point 2 coord is: 696770.781 6086089.772

Again it can be seen that this coordinate matches what was entered earlier.