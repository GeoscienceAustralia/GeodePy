Quickstart
==========

Ready to get started? This page gives a good introduction to using GeodePy.

First, ensure that:

* GeodePy is :ref:`installed <installation>`
* GeodePy is :ref:`up-to-date <updating>`

Let's get started with a simple example.

Transforming from GDA94 to GDA2020
----------------------------------

Begin by importing GeodePy:

.. code:: python

    import geodepy as gp   

Next, define a coordinate in GDA94:

.. code:: python

    coord_gda94 = gp.Coordinate(
        latitude=-33.865143,
        longitude=151.209900,
        height=30,
        datum=gp.datums.GDA94
    )

Now, transform this coordinate to GDA2020:

.. code:: python

    coord_gda2020 = gp.transform_datum(
        coord_gda94,
        target_datum=gp.datums.GDA2020
    )

Finally, print the transformed coordinate:

.. code:: python

    print(coord_gda2020)

This will output the latitude, longitude, and height of the coordinate in GDA2020.
For more detailed examples and tutorials, please refer to the :ref:`tutorials section <tutorials/index>`.