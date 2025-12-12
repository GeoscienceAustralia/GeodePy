.. _installation:

Installation
============

This section contains instructions for installing GeodePy.

Installing via pip
--------------------
The recommended way to install GeodePy is via pip. You can do this by running the following command in terminal:

.. code:: bash
    
    pip install geodepy


Installing from source
----------------------
If you prefer to install GeodePy from source, you can clone the repository from GitHub and use the python files:

.. code:: bash

   git clone https://github.com/GeoscienceAustralia/GeodePy.git

.. _updating:

Requirements
------------

GeodePy has some requirements that also need to be installed. Not all packages will be used for all modules.
The requirements can be seen below.

.. code::

    NumPy
    SciPy
    GDAL #only needed for height module

GDAL is only used for the heights module and sometimes can be difficult to install. For more information on 
installing GDAL see the `GDAL pypi <https://pypi.org/project/GDAL/>`_ page.

Updating GeodePy
-----------------
To update GeodePy to the latest version, you can use pip with the upgrade flag:

.. code:: bash

    pip install --upgrade geodepy  

.. _testing:

Testing GeodePy
----------------- 

GeodePy can be tested using the a unit test. This can be run as seen below.

.. caution:: The directory may need to be changed to point to the correct test folder

.. code:: python

    python3 -m unittest discover GeodePy/geodepy/tests --verbose

All tests should complete successfully. If not you may need to install some required 
packages.
