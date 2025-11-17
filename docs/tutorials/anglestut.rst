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

First a new angle needs to be defined. This will be done using the DDM Angle class

.. code:: python

    angle2 = angles.DDMAngle(40, 10.52)

Now this new anlge class can be added or subtracted from the first class

.. code:: python

    angle3 = angle1 + angle2
    angle4 = angle2 - angle1

    print(angle3)
    print(type(angle3))
    print(angle4)
    print(type(angle4))

    >>70 16 13.2
    >><class 'geodepy.angles.DMSAngle'>
    >>10 4.82
    >><class 'geodepy.angles.DDMAngle'>

As can be seen here, simple math can be completed on these classes. It should be noted that 
the result will have the class of the first variable in the calculation. 

.. note:: Operations can only be preformed when both angles are an object.

The following operators can be preformed on angle objects:

+----------------------+------------------+
| Operation            | Method           |
+======================+==================+
| Addition             |        '+'       |
+----------------------+------------------+
| Subtraction          |        '-'       |
+----------------------+------------------+
| Multiplication       |        '*'       |
+----------------------+------------------+
| Division             |        '/'       |
+----------------------+------------------+
| Equality             |        '=='      |
+----------------------+------------------+
| Not Equal            |        '!='      |
+----------------------+------------------+
| Less Than            |        '<'       |
+----------------------+------------------+
| Greater Than         |        '>'       |
+----------------------+------------------+
| Modulo               |        '%'       |
+----------------------+------------------+
| Round                |     round()      |
+----------------------+------------------+



