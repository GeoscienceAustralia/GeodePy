.. GeodePy documentation master file, created by
   sphinx-quickstart on Wed Nov  5 03:43:35 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

GeodePy: Geodesy in Python
======================================================

Release v0.6.0

.. image:: https://static.pepy.tech/personalized-badge/geodepy?period=total&units=INTERNATIONAL_SYSTEM&left_color=GREY&right_color=BLUE&left_text=downloads
    :target: https://pepy.tech/project/geodepy
    :alt: Total Downloads Badge

.. image:: https://img.shields.io/badge/license-Apache_2.0-green
    :target: https://opensource.org/licenses/Apache-2.0
    :alt: License Badge

.. image:: https://img.shields.io/badge/PyPi-GeodePy-yellow
    :target: https://pypi.org/project/geodepy/
    :alt: Pypi Badge

.. image:: https://img.shields.io/badge/Paper-GDA2020_Technical_Manual-556472
    :target: https://www.anzlic.gov.au/ICSM
    :alt: Paper Badge

**GeodePy** is a package of tools for manipulating geospatial datasets using Python.

-------------------

Features
--------

GeodePy includes a variety of :ref:`features <features/index>` for geodesy and geospatial data manipulation, including:

* :ref:`Converting <features/convert>` between coordinate types
* :ref:`Transforming <features/transform>` between datums
* Calculating :ref:`geodetic <features/geodesy>` distances and bearings
* Working with :ref:`geoid <features/height>` models
* :ref:`Surveying <features/survey>` calculations
* Various classes for :ref:`angles <features/angles>` , :ref:`corrdinates <features/coord>`, and :ref:`datums <features/constants>`
* :ref:`Statistics <features/statistics>`
* And more!

User Guide
--------------

This section contains the user guide for GeodePy, including installation, instructions, features and tutorials.

.. toctree::
   :maxdepth: 3
   :caption: User Guide

   installation
   quickstart
   features/index
   tutorials/index

Community Guide
---------------

This section contains information for contributors to the GeodePy project.

.. toctree::
   :maxdepth: 3
   :caption: Community Guide

   community/contributing
   community/authors
