.. _tutorials/angles:

Angle Classes and Converstions
=================================

GeodePy has 5 main angle classes to represent angles in different formats. These will be explored here along with how to convert between these types.

The 5 classes are:

* :ref:`Degrees, Minutes and Second (dms) <tut/dms>`
* :ref:`Degrees and Decimal Minutes (ddm) <tut/ddm>`
* :ref:`Decimal Degrees (dec) <tut/dd>`
* :ref:`HP notation (hpa) <tut/hp>`
* :ref:`Gradians (gona) <tut/grad>`

Classes
-------

.. _tut/dms:

Degrees, Minutes and Seconds (dms)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Angles expressed in degrees, minutes, and seconds.

- **Format:** ``ddd° mm' ss.s"``
- **Conversion:** 1° = 60′, 1′ = 60″
- **Example:** ``123° 34' 56.2"``

To initalise a dms class:

.. code:: python

    angle1 = geodepy.angles.DMSAngle(d, m, s)

.. _tut/ddm:

Degrees and Decimal Minutes (ddm)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Angles expressed in degrees and minutes, with minutes shown as a decimal fraction.

- **Format:** ``ddd° mm.mm'``
- **Example:** ``123° 34.933'``

To initalise a ddm class:

.. code:: python

    angle1 = geodepy.angles.DDMAngle(d, dm)

.. _tut/dd:

Decimal Degrees (dec)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Angles expressed entirely in decimal degrees.

- **Format:** ``ddd.ddd°``
- **Example:** ``123.5823°``

To initalise a dec class:

.. code:: python

    angle1 = geodepy.angles.DECAngle(d)

.. _tut/hp:

HP Notation (hpa)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Hemispheric Positive notation expresses latitude and longitude as positive values with hemisphere indicators.

- **Format:** ``ddd.mmssss``
- **Example:** ``123.231524°``

To initalise a hpa class:

.. code:: python

    angle1 = geodepy.angles.HPAngle(hp)

.. _tut/grad:

Gradians (gona)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A metric-based angle unit where a full circle equals 400 gradians.

- **Format** ``ggg.ggg``
- **Conversion:** 1 grad = 0.9°
- **Example:** ``137.5``

To initalise a gona class:

.. code:: python

    angle1 = geodepy.angles.GONAngle(gon)


Using Angle Classes
--------------------

First import GeodePy

.. code:: python

    import geodepy.angles as angles  

In this example a DMS angle will be created. This object can be initalised by including the degrees, minutes and seconds as arguments.

.. code:: python

    angle1 = angles.DMSAngle(30, 5, 42)
    print(angle1)

    >>30 5 42

Using this class we can get individual variables for degree minute and seconds componets seperately.

.. code:: python

    print(angle1.degree)
    print(angle1.minute)
    print(angle1.second)

    >>30
    >>5
    >>42

The methods within the class can also be used to convert the angle into different types. This is seen below:

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

.. caution:: Basic arthimitc should not be completed on HPA class. These should be converted to decimal degree first.

