.. _tutorials/edm:

EDM Corrections
================

Within GeodePy, EDM corrections can be completed. To complete an EDM correction the following data is needed.

- Carrier wavelength of instrument
- Modulation frequency and unit length or reference refractive index
- Distance measured
- Zentith angle
- Temperature
- Pressure
- Relative humidity
- (Optional) Co2 ppm

Using these we can calculate the first velocity correction and reduce to a horizontal distance.

First the variables need to be defined. Below is the values for a Lecia Viva.

.. code:: python

    import geodepy.survey

    wavelength = 0.658 #micrometers
    modFreq = 9.9902213*10e6 #hz
    unitLength = 1.5 #m
    dist = 145.265 #m
    zAngle = 91.15678 #dec
    temp = 26 #°c
    pressure = 1010.8 #hPa
    relHum = 37 #%
    co2ppm = 345 #ppm

First Velocity Parameters
-------------------------

Now we can calculate the first velocity parameters

.. code:: python

    params = first_vel_params(wavelength, modFreq, None, unitLength)

    print(params)

    >>(286.3433 80.6752)

First Velocity Correction
-------------------------

Using these parameters we can calculate the first velocity correction.

.. code:: python

    correction = geodepy.survey.first_vel_corrn(dist, params, temp, pressure, relHum)

    print(correction)

    >>0.0020656

Now this an be applied to the distance to get a corrected distance.

.. code:: python

    corrDist = dist + correction

    print(corrDist)

    >>145.2670656

The first velocity correction can also be calculated using the Co2 ppm. 
This gives a more accurate result but requires the co2 ppm and the wavelength.

.. code:: python

    correction2 = geodepy.survey.first_vel_corrn(dist, params, temp, pressure, relHum, 
                                                 wavelength=wavelength, CO2_ppm=co2ppm)

    print(correction2)

    print(dist + correction2)

    >>0.0020756
    >>145.2670756

Now the horizontal distance can be found using the zenith angle and corrected distance.

.. code:: python

    horzDist = geodepy.survey.va_conv(zAngle, corrDist)
    
    print(horzDist)

    >>(-1.15678, 145.2670656, 145.2374597, -2.9326875)

Where the third value in the tuple is the horizontal distance and the last value is the change in height.
