.. _tutorials/vert:

Vertical Datums
===============

This tutorial will work through converting between different vertical datums and working with gravity.
To learn more about vertical datums refer to the `GDA2020 technical manual <https://www.anzlic.gov.au/sites/default/files/files/GDA2020%20Technical%20Manual%20V1.8_published.pdf>`_.

Converting Between Vertical Datums
----------------------------------

Geodepy allows you to convert between different vertical datums within Australia. Below we will convert from 
an ellipsoidal height to AHD and then to AVWS.

.. code:: python

    import geodepy.height

First we will convert from ellipsoidal height to AHD.

.. code:: python

    rl = geodepy.height.GPS_to_AHD(-35.34405212, 149.15847673, 594.495)

    print(rl[0])

    >>[575.176]

This is the AHD height of this point. Now we can convert to AVWS.

.. code:: python

    rl = geodepy.height.AHD_to_AVWS(-35.34405212, 149.15847673,rl[0])

    print(rl[0])

    >>[575.303]

Finding Gravity Values
----------------------

Geodepy can also find gravity values at any point around Australia. This can be seen below.

.. code:: python

    import geodepy.height

To find the gravity at any point the geodepy.height.mean_surface_grav function can be used.
Normally this function is used to find the mean surface gravty between two points however
if the same lat and long are used for both point A and B then the gravity at a set location can be found.

.. code:: python

    grav = geodepy.height.mean_surface_grav(-35.34405212, 149.15847673, 575.176, -35.34405212, 149.15847673, 575.176,)

    print(grav)

    >>[9.79603591]

To calculate the average gravity between the ellipsoid and a certain height the following function can be used.

.. code:: python

    grav = geodepy.height.mean_normal_grav(-35.34405212, 594.495)

    print(grav)

    >>9.79671297

