.. _tutorials/vert:

Vertical Datums
===============

This tutorial will work through converting between different vertical datums and finding gravity values.
To learn more about vertical datums refer to the 
`GDA2020 technical manual <https://www.anzlic.gov.au/sites/default/files/files/GDA2020%20Technical%20Manual%20V1.8_published.pdf>`_.

Converting Between Vertical Datums
----------------------------------

GeodePy allows you to convert between different vertical datums within Australia. Below we will convert from 
an ellipsoidal height to AHD and then to AVWS.

.. code:: python

    import geodepy.height

First we will convert from ellipsoidal height to AHD.

.. code:: python

    rl = geodepy.height.GPS_to_AHD(
        -35.34405212, #lat
        149.15847673, #long
        594.495 #ellipsoidal height
    )

    print(rl[0])

    >>[575.176]

This is the AHD height. Now we can convert to AVWS.

.. code:: python

    rl = geodepy.height.AHD_to_AVWS(
        -35.34405212, #lat
        149.15847673, #long
        rl[0] #AHD
    )

    print(rl[0])

    >>[575.303]

Finding Gravity Values
----------------------

GeodePy can also find gravity values at any point around Australia. This can be seen below.

.. code:: python

    import geodepy.height

To find the gravity at any point the mean_surface_grav function can be used.
Normally this function is used to find the mean surface gravty between two points, however
if the same lat and long are used for both point A and B then the gravity at a set location can be found.

.. code:: python

    grav = geodepy.height.mean_surface_grav(
        -35.34405212, #lat of A
        149.15847673, #long of A
        575.176, #AHD
        -35.34405212, #lat of B
        149.15847673, #lat of B
        575.176 #AHD
    )

    print(grav)

    >>[9.79603591]

To calculate the average gravity between the ellipsoid and a certain height the following function can be used.

.. code:: python

    grav = geodepy.height.mean_normal_grav(
        -35.34405212, #lat
        594.495 #ellipsoid height
    )

    print(grav)

    >>9.79671297

