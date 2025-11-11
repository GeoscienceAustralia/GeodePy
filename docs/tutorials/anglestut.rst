.. _tutorials/angles:

Angle Classes and Converstions
=================================

GeodePy has 5 main angle classes to represent angles in different formats. These will be explored here along with how to convert between these types.

Creating Angle Class
--------------------

First import GeodePy

.. code:: python

    import geodepy.angles as angles  

In this example a DMS angle will be created. This object can be initalised by including the degrees, minutes and seconds as arguments.

.. code:: python

    angle1 = angles.DMSAngle(30, 5, 42)
    print(angle1)

    >> 30 5 42

The methods within the class can be used to convert the angle into different types. This is seen below:

.. code:: python

    print(angle1.ddm())
    print(angle1.dec())
    print(angle1.gona())
    print(angle1.hpa())
    print(angle1.rad())

    >>30 5.7
    >>30.095
    >>33.4388888889
    >>30.0542
    >>0.5252568383876934

This can be completed with any of the 5 angle classes within GeodePy. To complete math on the angle classes the following can be completed.


