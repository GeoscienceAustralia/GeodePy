Quickstart
==========

Ready to get started? This page gives a good introduction to using GeodePy.

First, ensure that:

* GeodePy is :ref:`installed <installation>`
* GeodePy is :ref:`up-to-date <updating>`

Let's get started with a simple example.

Transforming from MGAA94 to MGA2020
-----------------------------------

Begin by importing GeodePy:

.. code:: python

    import geodepy.transform as tr   

Next, define some coordinate values in MGA94:

.. code:: python

    zone_94 = 55
    east_94 = 696053.3373
    north_94 = 6086610.1338

Now, transform this coordinate to MGA2020:

.. code:: python

    (zone_20, east_20, north_20, _, _) = tr.transform_mga94_to_mga2020(
        zone_94,
        east_94,
        north_94
    )

Finally, print the transformed coordinate:

.. code:: python

    print(zone_20, east_20, north_20)

This will output the latitude, longitude, and height of the coordinate in MGA2020.

.. code::

    >> 55 696053.872 6086611.549

For more detailed examples and tutorials, please refer to the :ref:`features <features/index>` and :ref:`tutorial <tutorials/index>` section.