# -*- coding: utf-8 -*-

# example curl commands:
# curl 'http://www.ga.gov.au/bin/geodesy/run/sunrisenset' --data 'end=end&Location=&loc=yes&lathemi=south&LatDeg=-37&LatMin=57&LatSec=0&longhemi=east&LongDeg=144&LongMin=25&LongSec=0&TimeZonename=+&austzone=%2B10&TimeZone=%2B10&Event=7&Event=4&ZenithDeg=&ZenithMin=&ZenithSec=0&height=0&end=end&Date=17%2F07%2F2018&end=end'
# curl 'http://www.ga.gov.au/bin/geodesy/run/gda_vincenty?inverse=1&lat_degrees1=-37&lat_minutes1=57&lat_seconds1=03.72030&NamePoint1=FlindersPeak&lon_degrees1=144&lon_minutes1=25&lon_seconds1=29.52440&lat_degrees2=-37&lat_minutes2=39&lat_seconds2=10.15610&NamePoint2=Buninyong&lon_degrees2=143&lon_minutes2=55&lon_seconds2=35.38390&lat_deg1=-37+deg&lat_min1=57+min&lat_sec1=3.7203+sec&lon_deg1=144+deg&lon_min1=25+min&lon_sec1=29.5244+sec&lat_deg2=-37+deg&lat_min2=39+min&lat_sec2=10.1561+sec&lon_deg2=143+deg&lon_min2=55+min&lon_sec2=35.3839+sec&Submit=Submit+Data'
# curl 'http://www.ga.gov.au/bin/geodesy/run/gda_vincenty?inverse=0&lat_degrees1=-37&lat_minutes1=57&lat_seconds1=03.72030&NamePoint1=Flinders+Peak&lon_degrees1=144&lon_minutes1=25&lon_seconds1=29.52440&forward_azimuth_deg=306&forward_azimuth_min=52&forward_azimuth_sec=05.37&NamePoint2=Buninyong&ellipsoidal_dist=54972.217&lat_deg1=-37+deg&lat_min1=57+min&lat_sec1=3.7203+sec&lon_deg1=144+deg&lon_min1=25+min&lon_sec1=29.5244+sec&f_az_deg=306+deg&f_az_min=52+min&f_az_sec=5.37+sec&Submit=Submit+Data'

from urllib.request import urlopen
import os
import numpy as np
from bs4 import BeautifulSoup
import re
from geodepy.geodesy import vincdir, vincinv
from geodepy.transform import dd2dms, dms2dd


def parse_coord(coord):
    return re.split('[°\'"]+', str(coord))


def join_dms(degrees, minutes, seconds):
    seconds, subseconds = seconds.split('.')
    return degrees + '.' + minutes.zfill(2) + seconds.zfill(2) + subseconds


def split_dms(dms):
    degrees, rest = str(dms).split('.')
    return degrees, rest[:2], rest[2:4] + '.' + rest[4:]


def query_vincenty_inverse(lat1_deg, lat1_min, lat1_sec,
                           lon1_deg, lon1_min, lon1_sec,
                           lat2_deg, lat2_min, lat2_sec,
                           lon2_deg, lon2_min, lon2_sec):
    url = "http://www.ga.gov.au/bin/geodesy/run/gda_vincenty?inverse=1&" \
          "lat_degrees1={}&" \
          "lat_minutes1={}&" \
          "lat_seconds1={}&" \
          "NamePoint1=&" \
          "lon_degrees1={}&" \
          "lon_minutes1={}&" \
          "lon_seconds1={}&" \
          "lat_degrees2={}&" \
          "lat_minutes2={}&" \
          "lat_seconds2={}&" \
          "NamePoint2=&" \
          "lon_degrees2={}&" \
          "lon_minutes2={}&" \
          "lon_seconds2={}&" \
          "Submit=Submit+Data".format(lat1_deg, lat1_min, lat1_sec,
                                      lon1_deg, lon1_min, lon1_sec,
                                      lat2_deg, lat2_min, lat2_sec,
                                      lon2_deg, lon2_min, lon2_sec)

    data = urlopen(url).read()
    # with open('inverse.html', 'wb') as f:
    #     f.write(data)
    soup = BeautifulSoup(data, 'html.parser')
    # abs_path = os.path.abspath(os.path.dirname(__file__))
    # soup = BeautifulSoup(open(abs_path + '/inverse.html').read(), 'html.parser')

    fa_deg, fa_min, fa_sec = parse_coord(soup.find('td', text='Forward Azimuth:')
                                         .nextSibling
                                         .nextSibling
                                         .text
                                         .replace(' ', '')
                                         .replace('\'\'', ''))

    ra_deg, ra_min, ra_sec = parse_coord(soup.find('td', text='Reverse Azimuth:')
                                         .nextSibling
                                         .nextSibling
                                         .text
                                         .replace(' ', '')
                                         .replace('\'\'', ''))

    ell_dist = (soup.find('td', text='Ellipsoidal Distance:')
                .nextSibling
                .nextSibling
                .text
                .replace(' ', '')
                .replace('meters', '')
                .replace('\'\'', ''))

    return (fa_deg, fa_min, fa_sec,
            ra_deg, ra_min, ra_sec,
            ell_dist)


def query_vincenty_direct(lat1_deg, lat1_min, lat1_sec,
                          lon1_deg, lon1_min, lon1_sec,
                          fa_deg, fa_min, fa_sec,
                          ell_dist):
    url = ("http://www.ga.gov.au/bin/geodesy/run/gda_vincenty?inverse=0&"
           "lat_degrees1={}&"
           "lat_minutes1={}&"
           "lat_seconds1={}&"
           "NamePoint1=&"
           "lon_degrees1={}&"
           "lon_minutes1={}&"
           "lon_seconds1={}&"
           "forward_azimuth_deg={}&"
           "forward_azimuth_min={}&"
           "forward_azimuth_sec={}&"
           "NamePoint2=&"
           "ellipsoidal_dist={}&"
           "lat_deg1={}+deg&"
           "lat_min1={}+min&"
           "lat_sec1={}+sec&"
           "lon_deg1={}+deg&"
           "lon_min1={}+min&"
           "lon_sec1={}+sec&"
           "f_az_deg={}+deg&"
           "f_az_min={}+min&"
           "f_az_sec={}+sec&"
           "Submit=Submit+Data").format(lat1_deg, lat1_min, lat1_sec,
                                        lon1_deg, lon1_min, lon1_sec,
                                        fa_deg, fa_min, fa_sec,
                                        ell_dist,
                                        lat1_deg, lat1_min, lat1_sec,
                                        lon1_deg, lon1_min, lon1_sec,
                                        fa_deg, fa_min, fa_sec)

    data = urlopen(url).read()
    with open('direct2.html', 'wb') as f:
        f.write(data)
    soup = BeautifulSoup(data, 'html.parser')
    # abs_path = os.path.abspath(os.path.dirname(__file__))
    # soup = BeautifulSoup(open(abs_path + "/direct.html").read(), 'html.parser')

    lat_row = soup.find('td', text='Latitude:').parent
    lat2_deg, lat2_min, lat2_sec = tuple(parse_coord(lat.text
                                                     .replace(' ', '')
                                                     .replace('\'\'', ''))
                                         for lat in lat_row.td.next_siblings)[1]

    lon_row = soup.find('td', text='Longitude:').parent
    lon2_deg, lon2_min, lon2_sec = tuple(parse_coord(lon.text
                                                     .replace(' ', '')
                                                     .replace('\'\'', ''))
                                         for lon in lon_row.td.next_siblings)[1]

    ra_deg, ra_min, ra_sec = parse_coord(soup.find('td', text='Reverse Azimuth:')
                                         .nextSibling
                                         .nextSibling
                                         .text
                                         .replace(' ', '')
                                         .replace('\'\'', ''))

    return (lat2_deg, lat2_min, lat2_sec,
            lon2_deg, lon2_min, lon2_sec,
            ra_deg, ra_min, ra_sec)


def vincenty_from_website():
    np.set_printoptions(suppress=True)
    abs_path = os.path.abspath(os.path.dirname(__file__))
    coords = np.genfromtxt(os.path.join(abs_path, 'resources/natadjust_rvs_example.dat'),
                           usecols=[5, 6])

    # pair (1st, last), (2nd, 1st), ..., (2nd last, 3rd last), (last, 2nd last)
    coord_pairs = np.column_stack([coords, np.roll(coords, 2)])

    # coord_pairs = np.array([[-37.570372030, 144.252952440, -37.391015610, 143.553538390]])

    for row_number, (lat1, lon1, lat2, lon2) in enumerate(coord_pairs):
        # inverse vincentys
        print(("Testing vincentys_inverse({}, {}, {}, {})"
               " on row {}").format(lat1, lon1, lat2, lon2, str(row_number)))

        (web_inv_fa_deg, web_inv_fa_min, web_inv_fa_sec,
         web_inv_ra_deg, web_inv_ra_min, web_inv_ra_sec,
         web_inv_ell_dist) = query_vincenty_inverse(*split_dms(lat1), *split_dms(lon1),
                                                    *split_dms(lat2), *split_dms(lon2))
        web_inv_az1to2, web_inv_az2to1 = (join_dms(web_inv_fa_deg, web_inv_fa_min, web_inv_fa_sec),
                                          join_dms(web_inv_ra_deg, web_inv_ra_min, web_inv_ra_sec))

        inv_ell_dist, inv_az1to2, inv_az2to1 = vincinv(dms2dd(lat1), dms2dd(lon1), dms2dd(lat2), dms2dd(lon2))

        np.testing.assert_almost_equal(inv_ell_dist, float(web_inv_ell_dist),
                                       decimal=3,
                                       err_msg=("vincentys_inverse({}, {}, {}, {})"
                                                " on row {}").format(lat1, lon1, lat2, lon2, str(row_number)))
        np.testing.assert_almost_equal((dd2dms(inv_az1to2), dd2dms(inv_az2to1)),
                                       (float(web_inv_az1to2), float(web_inv_az2to1)),
                                       decimal=6,
                                       err_msg=("vincentys_inverse({}, {}, {}, {})"
                                                " on row {}").format(lat1, lon1, lat2, lon2, str(row_number)))

        # direct vincentys
        print("Testing vincentys_direct({}, {}, {}, {}) on row {}".format(lat1,
                                                                          lon1,
                                                                          dd2dms(inv_az1to2),
                                                                          dd2dms(inv_ell_dist),
                                                                          str(row_number)))
        (lat2_deg, lat2_min, lat2_sec,
         lon2_deg, lon2_min, lon2_sec,
         web_dir_ra_deg, web_dir_ra_min, web_dir_ra_sec) = query_vincenty_direct(*split_dms(lat1), *split_dms(lon1),
                                                                                 web_inv_fa_deg, web_inv_fa_min,
                                                                                 web_inv_fa_sec,
                                                                                 web_inv_ell_dist)
        web_dir_lat2, web_dir_lon2, web_dir_az2to1 = (join_dms(lat2_deg, lat2_min, lat2_sec),
                                                      join_dms(lon2_deg, lon2_min, lon2_sec),
                                                      join_dms(web_dir_ra_deg, web_dir_ra_min, web_dir_ra_sec))

        dir_lat2, dir_lon2, dir_az2to1 = vincdir(dms2dd(lat1), dms2dd(lon1), inv_az1to2, inv_ell_dist)

        np.testing.assert_almost_equal((float(web_dir_lat2), float(web_dir_lon2)),
                                       (dd2dms(dir_lat2), dd2dms(dir_lon2)),
                                       decimal=8,
                                       err_msg="vincentys_direct({}, {}, {}, {})"
                                               "on row {}".format(lat1, lon1, dd2dms(inv_az1to2),
                                                                  inv_ell_dist, str(row_number)))
        np.testing.assert_almost_equal(float(web_dir_az2to1),
                                       dd2dms(dir_az2to1),
                                       decimal=6,
                                       err_msg="vincentys_direct({}, {}, {}, {})"
                                               "on row {}".format(lat1, lon1, dd2dms(inv_az1to2),
                                                                  inv_ell_dist, str(row_number)))
        print("--------")


if __name__ == "__main__":
    vincenty_from_website()
