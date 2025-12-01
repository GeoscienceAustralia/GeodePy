.. _tutorials/timetrans:

Time Dependant Transformations
==============================

In this tutorial we will discuss time dependant transformations. 
To see transformations between static datums view the :ref:`datum transformation <tutorials/transform>` tutorial. 
To learn more about time dependant transformations refer to the `GDA2020 technical manual <https://www.anzlic.gov.au/sites/default/files/files/GDA2020%20Technical%20Manual%20V1.8_published.pdf>`_.

Time dependant transformation are much more complex then transformation between static datums. First we will complete a simply example.

Common Example
--------------

Here we will convert between ATRF2014 and GDA2020. This is transforming from a dynamic datum to a static datum.

In this example we will transform a coordinate in ATRF2014 at epoch 1/1/2011.

First we will import geodepy

.. code:: python

    import geodepy.tranform
    from datetime import date

Now the function for transforming from ATRF2014 to GDA2020 can be used.

.. code:: python

    x_20, y_20, z_20, vcv = geodepy.transform.transform_atrf2014_to_gda2020(-4050762.770917, 4220880.800229, -2533400.199554, date(2011,1,1))
    print(x_20, y_20, z_20)

    >>-4050763.124034 4220880.753100 -2533399.713463

This is the GDA2020 coordinate.

Transforming between Dynamic and Static Datums
----------------------------------------------

The simpliest time dependant transformations go between static and dynamic datums. This is because often the 
transformation paramters are already present. Below we will complete an example that doesnt include a dedicated function.

Here we will transform from GDA94 to ITRF2008 at 1/1/2007.

.. code:: python

    import geodepy.tranform
    import geodepy.constants
    from datetime import date

To go from GDA94 to ITRF2008 we need to investigate what transformations are present in GeodePy. This can be found 
in :ref:`this <features/constants/transform>` table. Here we can see that there is a direct transformation from GDA94 to ITRF2008.
We will use this to complete our transformation.

.. code:: python

    x, y, z, vcv = geodepy.transform.conform14(-4050763.124034, 4220880.753100, -2533399.713463, date(2007,1,1), geodepy.constants.gda94_to_itrf2008)

    print(x, y, z)

    >>-4050763.585602 4220880.607715 -2533398.980318

This is the ITRF2008 coordinate on 1/1/2007

Transforming Between Two Dynamic Datums
---------------------------------------

Transforming between two dynamic datums is more complex, requiring a few more  steps and considerations. 
For tihs example we will transform from ITRF2008 at 1/1/2007 to ITRF2020 at 1/1/2030.

.. code:: python

    import geodepy.constant
    import geodepy.transform

First the ITRF2008 coordinates need to be converted to ITRF2014. This is completed so that the plate motion 
between 2007 and 2030 can be applied. The plate motion can only be applied to coordinates in ITRF2014 or 
ATRF2014. When completeing this transformation the date of the ITRF2008 epoch is entered. This means the 
resulting ITRF2020 cooridnate will be at the IRTF2008 epoch.

.. code:: python

    x, y, z, vcv = geodepy.transform.conform14(-4050762.612530, 4220880.821783, -2533400.416214, date(2007, 1, 1), geodepy.constants.itrf2008_to_itrf2014)

    print(x, y, z)

    >>-4050762.614575 4220880.820347 -2533400.419192

Now we have an ITRF2014 coordinate at 1/1/2007. Now this needs to be moved to the 1/1/2030. This can be 
done using the ITRF2014 to GDA2020 transformation which approximates plate motion in Australia. To complete 
this transformation on another plate a different plate motion model should be used. 
Some carfeul math needs to be completed here. To go from 2007 to 2030, 23 years 
of plate motion needs to be added. The reference epoch of the ITRF2014 to GDA2020 transformation is 2020. 
As such 23 needs to be subtracted from 2020 to get the desired motion. This means the epoch 1/1/1993 should be entered.

.. caution::
    Transformations using the plate motion model of itrf2014_to_gda2020 should only be completed for epochs between
    2005 - 2035. For transformations outside of this range refer to the :ref:`next <tutorials/transold>` section.


.. code:: python

    x, y, z, vcv = geodepy.transform.conform14(x, y, z, date(1997, 1, 1), geodepy.constants.itrf2014_to_gda2020)

    print(x, y, z)

    >>-4050763.516973 4220880.699906 -2533399.176976

This is now the ITRF2014 corrdinate at 1/1/2030. Now we can convert this ITRF2014 cooridnate to ITRF2020.

.. code:: python

    x, y, z, vcv = geodepy.transform.conform14(x, y, z, date(2030,1,1), geodepy.constants.itrf2014_to_itrf2020)

    print(x, y, z)

    >>-4050763.517274 4220880.704079 -2533399.182440

This is the final cooridnate in ITRF2020 at 1/1/2030.

.. _tutorials/transold:

Transforming Between Older Dynamic Datums
-----------------------------------------

The Australia plate motion model should only be used between the years of 2005 to 2035. If a datum older then 
this needs to be transformed a different method should be used. For this example we will tranform from ITRF88 
at 1/1/1988 to ITRF2014 at 1/1/2030.

.. caution:: This method only works for coordinates within Australia.

.. code:: python

    import geodepy.constant
    import geodepy.transform

To go from an old dynamic datum to a more current datum, transformations to static datums should be completed first.
In this case transforming to GDA2020 is most accurate, before then transforming back to ITRF2014. To get to GDA2020
the ITRF88 coordinate first needs to be transformed to ITRF2014 at the 1/1/1988.

.. code:: python

    x, y, z, vcv = geodepy.transform.conform14(-4050763.124645, 4220880.752269, -2533399.717044, date(1988,1,1), geodepy.constants.itrf88_to_itrf2014)

    print(x, y, z)

    >>-4050763.116490 4220880.700494 -2533399.614981

This is the ITRF2014 coordinate at 1/1/1988. Now this needs to be transformed into GDA2020

.. code:: python

    x, y, z, vcv = geodepy.transform.conform14(x, y, z, date(1988,1,1), geodepy.constants.itrf2014_to_gda2020)

    print(x, y, z)

    >>-4050764.372112 4220880.532909 -2533397.886526

Now this corrinate can be changed to ITRF2014 at 1/1/2023.

.. code:: python

    x, y, z, vcv = geodepy.transform.conform14(x, y, z, date(2030,1,1), geodepy.constants.gda2020_to_itrf2014)

    print(x, y, z)

    >>-4050764.764547 4220880.480531 -2533397.346309

This is now the ITRF2014 coordinate at 1/1/2030.