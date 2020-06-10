
__all__ = ["grid2geoio", "geo2gridio", "gdatrans7"]

import geodepy.convert
import pandas as pd
import tkinter as ttk
import geodepy.transform as tf
import geodepy.convert as cv
import geodepy.constants as cs

'''
The inputoutput module acts as the backend for the GUI and manages the calls
to the GeodePy functions. The functions are rewritten  from the transform
module but use pandas to improve efficiency. Functions placed here to not
overcrowd the transform module.
'''


def grid2geoio(fn, fn_out, easting, northing, utmzone, geotypeout):
    """
        Input:
         The CSV data must have headers.

         Inputs passed the the function:
         :param fn: input file path
         :param fn_out: file output path
         :param easting: is the column name in the csv where the eastings are stored
         :param northing: is the column name in the csv where the northings are stored
         :param utmzone: is the column name in the csv where the UTM Zone is stored
         :param geotypeout: format of latitude and longitude output e.g. DD or DMS

        Output:
        Data in the file is output to the fn_out location as a CSV.
        """

    # Check whether geotypeout value is from GUI or if function called within other code, then check that a valid
    # geotypeout value was provided
    if isinstance(geotypeout, ttk.StringVar):
        geotypeout = geotypeout.get()
    geotypeout_options = ["DD", "DMS"]
    if geotypeout not in geotypeout_options:
        raise ValueError("Invalid geotypeout. Expected one of: %s" % geotypeout_options)

    # Opens the CSV as a pandas dataframe
    csvdf = pd.read_csv(fn, low_memory=False)
    # Sends data to tf.grid2geo
    returned = csvdf.apply(lambda x: geodepy.convert.grid2geo(x[utmzone], x[easting], x[northing]), axis=1)
    # Convert tuple returned to DafaFrame
    resultdf = pd.DataFrame(list(returned), columns=['Latitude', 'Longitude', 'Point Scale Factor', 'Grid Convergence'])
    if geotypeout == 'DMS':
        resultdf['Latitude'] = resultdf['Latitude'].apply(cv.dec2hp)
        resultdf['Longitude'] = resultdf['Longitude'].apply(cv.dec2hp)
    # Adds the results to the original dataframe from the csv
    csvdf = pd.concat([csvdf, resultdf], axis=1)
    # Writes the dataframe to a csv
    csvdf.to_csv(fn_out, index=False)


def geo2gridio(fn, fn_out, latitude, longitude, geotypein):
    """
        Input:
        The CSV data must have headers.

         Inputs passed the the function:
         :param fn: input file path
         :param fn_out: file output path
         :param latitude: is the column name in the csv where the latitudes are stored
         :param longitude: is the column name in the csv where the longitudes are stored
         :param geotypein: format of latitude and longitude e.g. DD or DMS

        Output:
        Data in the file is output to the fn_out location as a CSV.
        """

    # Check whether geotypein value is from GUI or if function called within other code, then check that a valid
    # geotypein value was provided
    if isinstance(geotypein, ttk.StringVar):
        geotypein = geotypein.get()
    geotypein_options = ["DD", "DMS"]
    if geotypein not in geotypein_options:
        raise ValueError("Invalid geotypein. Expected one of: %s" % geotypein_options)

    # Opens the CSV as a pandas dataframe
    csvdf = pd.read_csv(fn, low_memory=False)

    # Converts DMS lat and long to Decimal Degrees if required
    if geotypein == 'DMS':
        csvdf['LatitudeDD'] = csvdf[latitude].apply(cv.hp2dec)
        csvdf['LongitudeDD'] = csvdf[longitude].apply(cv.hp2dec)
        latitude = 'LatitudeDD'
        longitude = 'LongitudeDD'

    # Sends data to tf.geo2grid
    returned = csvdf.apply(lambda x: geodepy.convert.geo2grid(x[latitude], x[longitude]), axis=1)
    # Convert tuple returned from tf.gepo2grid into DafaFrame
    resultdf = pd.DataFrame(list(returned), columns=['Hemisphere', 'UTMZone', 'Easting', 'Northing',
                                                     'Point Scale Factor', 'Grid Convergence'])
    # Adds the results to the original dataframe from the csv
    csvdf = pd.concat([csvdf, resultdf], axis=1)
    # Writes the dataframe to a csv
    csvdf.to_csv(fn_out, index=False)


def gdatrans7(fn, fn_out, latitude, longitude, ellht, gdageotypein, direction):

    """
    Input:
     The CSV data must have headers.

     Inputs passed the the function:
     :param fn: input file path
     :param fn_out: file output path
     :param latitude: is the column name in the csv where the latitudes are stored
     :param longitude: is the column name in the csv where the longitudes are stored
     :param ellht: is the column name in the csv where the ellipsoidal heights are stored
     :param gdageotypein: format of latitude and longitude e.g. DD or DMS
     :param direction: either "94to2020" or "2020to94". Specifies the datum to transform from and to.

     trans: Transformation parameter from the geodepy.constants module. The trans parameters are hard coded below.
     Parameters were sourced from GDA2020 Technical Manual. Need to work on a way of entering and passing different
     transformation parameters to the function.

    Output:
    Data in the file is output to the fn_out location as a CSV.
    """

    # Check whether direction value is from GUI or if function called within other code, then check that a valid
    # direction value was provided
    if isinstance(direction, ttk.StringVar):
        direction = direction.get()
    direction_options = ["94to2020", "2020to94"]
    if direction not in direction_options:
        raise ValueError("Invalid direction. Expected one of: %s" % direction_options)

    # Check whether gdageotypein value is from GUI or if function called within other code, then check that a valid
    # gdageotypein value was provided
    if isinstance(gdageotypein, ttk.StringVar):
        gdageotypein = gdageotypein.get()
    gdageotypein_options = ["DD", "DMS"]
    if gdageotypein not in gdageotypein_options:
        raise ValueError("Invalid gdageotypein. Expected one of: %s" % gdageotypein_options)

    # Sets up the direction of the transformation and sets up the output column naming
    if direction == "94to2020":
        trans = cs.Transformation('MGA94', 'MGA2020', '1994', 0.06155, -0.01087, -0.04019, -0.009994, -0.0394924,
                                  -0.0327221, -0.0328979)
        translat = "GDA2020Latitude"
        translong = "GDA2020Longitude"
        transellht = "GDA2020EllHt"
    else:
        trans = cs.Transformation('MGA94', 'MGA2020', '1994', 0.06155, -0.01087, -0.04019, -0.009994, -0.0394924,
                                  -0.0327221, -0.0328979).__neg__()
        translat = "GDA1994Latitude"
        translong = "GDA1994Longitude"
        transellht = "GDA1994EllHt"

    # Opens the CSV as a pandas dataframe
    csvdf = pd.read_csv(fn, low_memory=False)

    # Converts DMS lat and long to Decimal Degrees if required
    if gdageotypein == 'DMS':
        csvdf['LatitudeDMS'] = csvdf[latitude]
        csvdf['LongitudeDMS'] = csvdf[longitude]
        csvdf['Latitude'] = csvdf[latitude].apply(cv.hp2dec)
        csvdf['Longitude'] = csvdf[longitude].apply(cv.hp2dec)

    # Converts Lat, Long & ElipHt XYZ
    returned = csvdf.apply(lambda x: geodepy.convert.llh2xyz(x[latitude], x[longitude], x[ellht]), axis=1)
    # Converts the results from llh2xyz into a dataframe
    xyzresultdf = pd.DataFrame(list(returned), columns=['CartesianX', 'CartesianY', 'CartesianZ'])

    # Combines the results with the original dataframe
    csvdf = pd.concat([csvdf, xyzresultdf], axis=1)

    # Transforms Lat, Long & ElipHt XYZ using a 7 parameter transformation
    transformed = csvdf.apply(lambda x: tf.conform7(x.CartesianX, x.CartesianY, x.CartesianZ, trans), axis=1)
    # Converts the results from conform7 into a dataframe
    transresultdf = pd.DataFrame(list(transformed), columns=['TransCartesianX', 'TransCartesianY', 'TransCartesianZ'])

    # Combines the results with the original dataframe
    csvdf = pd.concat([csvdf, transresultdf], axis=1)

    # Converts the the transformed XYZ coordinates back to Lats, Longs and Ellipsoidal Heights
    returned = csvdf.apply(lambda x: geodepy.convert.xyz2llh(x.TransCartesianX, x.TransCartesianY, x.TransCartesianZ), axis=1)
    # Converts the results from llh2xyz into a dataframe
    llhresultdf = pd.DataFrame(list(returned), columns=[translat, translong, transellht])

    # Combines the results with the original dataframe
    csvdf = pd.concat([csvdf, llhresultdf], axis=1)
    csvdf = csvdf.drop(['CartesianX', 'CartesianY', 'CartesianZ', 'TransCartesianX', 'TransCartesianY',
                        'TransCartesianZ'], axis=1)

    # Writes the csv dataframe to a csv file
    csvdf.to_csv(fn_out, index=False)
