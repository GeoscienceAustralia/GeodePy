.. _features/angles:
 
Angles
======

Angle Types
-----------

GeodePy supports Angular Notation in 9 different formats

+------------------------+--------------------------------------------------------------+
| ABRV                   | FORMAT (type)                                                |
+========================+==============================================================+
| rad                    | Radians (stored as float)                                    |
+------------------------+--------------------------------------------------------------+
| dec                    | Decimal Degrees (stored as float)                            |
+------------------------+--------------------------------------------------------------+
| :ref:`dec <DECAngle>`  | Decimal Degrees (via DECAngle class)                         |
+------------------------+--------------------------------------------------------------+
| hp                     | Hewlett Packard (HP) Notation (stored as float)              |
+------------------------+--------------------------------------------------------------+
| :ref:`hpa <HPAngle>`   | Hewlett Packard (HP) Notation (via HPAngle class)            |
+------------------------+--------------------------------------------------------------+
| gon                    | Gradians (stored as float)                                   |
+------------------------+--------------------------------------------------------------+
| :ref:`gona <GONAngle>` | Gradians (via GONAngle class)                                |
+------------------------+--------------------------------------------------------------+
| :ref:`dms <DMSAngle>`  | Degrees, Minutes and Seconds Notation (via DMSAngle class)   |
+------------------------+--------------------------------------------------------------+
| :ref:`ddm <DDMAngle>`  | Degrees and Decimal Minutes Notation (via DDMAngle class)    |
+------------------------+--------------------------------------------------------------+

Conversion between all formats is supported as shown below:

* Radians to/from Decimal Degrees via builtin math.radians and math.degrees

* Formats as floats to all other types via functions in the form abrv2abrv
  
  e.g. gon2hpa()

* DECAngle, HPAngle, GONAngle, DMSAngle and DDMAngle class objects via methods in
  the form CLASS.abrv()
  
  e.g. HPAngle(value).dec()

All **angle classes** can be seen below:

.. _DECAngle:

.. autoclass:: geodepy.angles.DECAngle
    :members:

.. _HPAngle:

.. autoclass:: geodepy.angles.HPAngle
    :members:

.. _GONAngle:

.. autoclass:: geodepy.angles.GONAngle
    :members:

.. _DMSAngle:

.. autoclass:: geodepy.angles.DMSAngle
    :members:

.. _DDMAngle:

.. autoclass:: geodepy.angles.DDMAngle
    :members:

.. _converstions:

-------------------

Angle Conversions
-----------------

All conversion functions can be seen below.

For converting **decimal degree** to other formats:

.. autofunction:: geodepy.angles.dec2hp

.. autofunction:: geodepy.angles.dec2hpa

.. autofunction:: geodepy.angles.dec2gon

.. autofunction:: geodepy.angles.dec2gona

.. autofunction:: geodepy.angles.dec2dms

.. autofunction:: geodepy.angles.dec2ddm

-------------------

For converting **Hewlett Packard** to other formats:


.. autofunction:: geodepy.angles.hp2dec

.. autofunction:: geodepy.angles.hp2deca

.. autofunction:: geodepy.angles.hp2rad

.. autofunction:: geodepy.angles.hp2gon

.. autofunction:: geodepy.angles.hp2gona

.. autofunction:: geodepy.angles.hp2dms

.. autofunction:: geodepy.angles.hp2ddm

-------------------

For converting **Gradians** to other formats:

.. autofunction:: geodepy.angles.gon2dec

.. autofunction:: geodepy.angles.gon2deca

.. autofunction:: geodepy.angles.gon2hp

.. autofunction:: geodepy.angles.gon2hpa

.. autofunction:: geodepy.angles.gon2rad

.. autofunction:: geodepy.angles.gon2dms

.. autofunction:: geodepy.angles.gon2ddm

-------------------

Other Converstion Functions:

.. autofunction:: geodepy.angles.dd2sec

..

  .. autofunction:: geodepy.angles.dec2hp_v

  .. autofunction:: geodepy.angles.hp2dec_v

Angle Type Checking
-------------------

.. autofunction:: geodepy.angles.angular_typecheck
