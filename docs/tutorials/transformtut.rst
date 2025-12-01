.. _tutorials/transform:

Datum Transformation
=====================

GeodePy has the ability to tranform between datums. Here we will discuss how to 
change between datums without changing the reference epoch while in the 
:ref:`time dependant <tutorials/timetrans>` tutorial we will discuss changing epochs.
To learn more about transformation refer to the `GDA2020 technical manual <https://www.anzlic.gov.au/sites/default/files/files/GDA2020%20Technical%20Manual%20V1.8_published.pdf>`_.

Common Example
--------------

The most common example for datum transformation is converting from MGA94 to MGA2020. 
This is handled by a function in the transformation module that converts the input to 
grid input to caresian coordinate (xyz), and then runs a 7 paramter helmert transformation 
using the GDA94_to_gda2020 transformation constant. This process will be explored more in 
the :ref:`latter <tutorials/transformfunc>` part of this tutorial but first lets use this function.

Begin by importing GeodePy:

.. code:: python

    import geodepy.transform   

Next, define some coordinate values in MGA94:

.. code:: python

    zone_94 = 55
    east_94 = 696053.3373
    north_94 = 6086610.1338

Now, transform this coordinate to MGA2020:

.. code:: python

    (zone_20, east_20, north_20, _, _) = geodepy.transform.transform_mga94_to_mga2020(
        zone_94,
        east_94,
        north_94
    )

    print(zone_20, east_20, north_20)

    >>55 696053.872 6086611.549

This it the MGA2020 coordinates.

To complete this transformation a function from the transformation module was used. This 
will not always be the case. In the next section this will be explored more.

.. _tutorials/transformfunc:

Constructing a new Transformation Function
------------------------------------------

Here we will explore how to complete a transformation between datums without a dedicated function.
We will then create a new function for this new transformation. Here we are going to transform between AGD84 to GDA2020

First import GeodePy.

.. code:: python

    import geodepy.transform 
    import geodepy.constants
    import geodepy.angles

Now we need some starting coordinates

.. code:: python

    lat = geodepy.angles.DMSAngle(-23,33,25.21)
    long = geodepy.angles.DMSAngle(133,49,13.87)
    height = 427.863 

    print(f"The AGD84 position is {lat}, {long}, {height}"

    >>The AGD84 position is -23 33 25.21, 133 49 13.87, 427.863

All transformations in GeodePy need to be completed with corrdinates in cartesian (xyz) form. Lets transform to xyz.

.. code:: python

    x, y, z = geodepy.transform.llh2xyz(lat, long, height)
    print(x, y, z)

    >>-4050634.051819 4220935.13646 -2533555.369303

Now we need some transformation parameters. Within Geodepy there are many transformaton parameters already 
present within the constants module. A table of these can be seen :ref:`here <features/constants/transform>`. If the 
transformation needed isn't currently in Geodepy, new transformations can be added using the :ref:`transformation class <transclass>`. 
Here the agd66_to_gda94 transformation will be used.

.. code:: python

    print(geodepy.constants.agd84_to_gda94)

    >>Transformation: From 'AGD84' to 'GDA94'
    >>Reference Epoch: 0
    >>tx: -117.763m + 0.0m/yr
    >>ty: -51.51m + 0.0m/yr
    >>tz: 139.061m + 0.0m/yr
    >>sc: -0.191ppm + 0.0ppm/yr
    >>rx: -0.292" + 0.0"/yr
    >>ry: -0.443" + 0.0"/yr
    >>rz: -0.277" + 0.0"/yr

Now this will be used to complete a 7 paramter helmert transformation.

.. code:: python

    x_94, y_94, z_94, _ = geodepy.transform.conform7(x, y, z, geodepy.constants.agd84_to_gda94)
    print(x_94, y_94, z_94)

    >>-4050762.150962 4220880.96717 -2533401.14935

.. tip:: 
    The "_" in the conform7 variables is for a vcv matrix. If you dont input a vcv matrix 
    then the resulting variable will be None. However as the variable needs to be assigned 
    for the function to work using "_" meets to requirement but doesn't store the variable 
    in a meaningful way.

Now we need to transform from GDA94 to GDA2020 using the same method but with the gda94_to_gda2020 
transformation class.

.. code:: python

    x_20, y_20, z_20, _ = geodepy.transform.conform7(x_94, y_94, z_94, geodepy.constants.gda94_to_gda2020)

    print(x_20, y_20, z_20)

    >>-4050763.124034 4220880.75310 -2533399.713463

Now the cartesian coordinates can be converted back to geographic coordinates (lat, long)

.. code:: python

    lat_20, long_20, height_20 = geodepy.transform.xyz2llh(x_20, y_20, z_20)

    print(f"The GDA2020 position is {geodepy.angles.dec2dms(lat_20)}, {geodepy.angles.dec2dms(long_20)}, {height_20}")

    >>The GDA2020 position is -23 33 19.921108332, 133 49 18.481110119, 411.61055134

All of this can now be combined into one function.

.. code:: python

    def transform_agd84_to_gda2020(lat, long, height):

        x, y, z = geodepy.transform.llh2xyz(lat, long, height)
        x_94, y_94, z_94, _ = geodepy.transform.conform7(x, y, z, geodepy.constants.agd84_to_gda94)
        x_20, y_20, z_20, _ = geodepy.transform.conform7(x_94, y_94, z_94, geodepy.constants.gda94_to_gda2020)

        return geodepy.transform.xyz2llh(x_20, y_20, z_20)

    lat_new, long_new, height_new = transform_agd84_to_gda2020(lat, long, height)

    print(f"The GDA2020 position is {geodepy.angles.dec2dms(lat_new)}, {geodepy.angles.dec2dms(long_new)}, {height_new}")

    >>The GDA2020 position is -23 33 19.921108332, 133 49 18.481110119, 411.61055134

