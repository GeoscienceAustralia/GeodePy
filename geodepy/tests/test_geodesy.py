import unittest
import os.path
import numpy as np
import numpy.lib.recfunctions as rfn
from geodepy.convert import (
    hp2dec,
    dec2hp,
    rect2polar,
    polar2rect,
    grid2geo,
    llh2xyz,
    DMSAngle,
)
from geodepy.geodesy import vincinv, vincdir, vincinv_utm, vincdir_utm, enu2xyz, xyz2enu


class TestGeodesy(unittest.TestCase):
    def test_enu2xyz(self):
        MOBS_MGA2020 = (55, 321820.085, 5811181.510, 40.570)
        MOBS_MGA1994 = (55, 321819.594, 5811180.038, 40.659)

        # Convert UTM Projection Coordinates to Geographic Coordinates
        MOBS_GDA2020 = grid2geo(MOBS_MGA2020[0], MOBS_MGA2020[1], MOBS_MGA2020[2])
        MOBS_GDA1994 = grid2geo(MOBS_MGA1994[0], MOBS_MGA1994[1], MOBS_MGA1994[2])

        # Convert Geographic Coordinates to Cartesian XYZ Coordinates
        MOBS_GDA2020_XYZ = llh2xyz(MOBS_GDA2020[0], MOBS_GDA2020[1], MOBS_MGA2020[3])
        MOBS_GDA1994_XYZ = llh2xyz(MOBS_GDA1994[0], MOBS_GDA1994[1], MOBS_MGA1994[3])

        # Generate Vector Between UTM Projection Coordinates
        mga_vector = [
            MOBS_MGA2020[1] - MOBS_MGA1994[1],
            MOBS_MGA2020[2] - MOBS_MGA1994[2],
            MOBS_MGA2020[3] - MOBS_MGA1994[3],
        ]

        # Generate Vector Between Cartesian XYZ Coordinates
        xyz_vector = (
            MOBS_GDA2020_XYZ[0] - MOBS_GDA1994_XYZ[0],
            MOBS_GDA2020_XYZ[1] - MOBS_GDA1994_XYZ[1],
            MOBS_GDA2020_XYZ[2] - MOBS_GDA1994_XYZ[2],
        )

        # Rotate UTM Projection Vector by Grid Convergence
        grid_dist, grid_brg = rect2polar(mga_vector[0], mga_vector[1])
        local_east, local_north = polar2rect(grid_dist, grid_brg - MOBS_GDA2020[3])
        local_vector = (local_east, local_north, mga_vector[2])

        # Calculate XYZ Vector using Local Vector Components
        x, y, z = enu2xyz(MOBS_GDA2020[0], MOBS_GDA2020[1], *local_vector)
        self.assertAlmostEqual(x, xyz_vector[0], 4)
        self.assertAlmostEqual(y, xyz_vector[1], 4)
        self.assertAlmostEqual(z, xyz_vector[2], 4)

        # Calculate Local Vector using XYZ Vector Components
        e, n, u = xyz2enu(MOBS_GDA2020[0], MOBS_GDA2020[1], *xyz_vector)
        self.assertAlmostEqual(e, local_vector[0], 4)
        self.assertAlmostEqual(n, local_vector[1], 4)
        self.assertAlmostEqual(u, local_vector[2], 4)

    def test_vincinv(self):
        # Flinders Peak
        lat1 = hp2dec(-37.57037203)
        lon1 = hp2dec(144.25295244)
        lat1_DMS = DMSAngle(-37, 57, 3.7203)
        lon1_DMS = DMSAngle(144, 25, 29.5244)

        # Buninyong
        lat2 = hp2dec(-37.39101561)
        lon2 = hp2dec(143.55353839)
        lat2_DMS = DMSAngle(-37, 39, 10.1561)
        lon2_DMS = DMSAngle(143, 55, 35.3839)

        # Test Decimal Degrees Input
        ell_dist, azimuth1to2, azimuth2to1 = vincinv(lat1, lon1, lat2, lon2)
        self.assertEqual(round(ell_dist, 3), 54972.271)
        self.assertEqual(round(dec2hp(azimuth1to2), 6), 306.520537)
        self.assertEqual(round(dec2hp(azimuth2to1), 6), 127.102507)

        # additional test case:
        pl1 = (-29.85, 140.71666666666667)
        pl2 = (-29.85, 140.76666666666667)
        ell_dist, azimuth1to2, azimuth2to1 = vincinv(pl1[0], pl1[1], pl2[0], pl2[1])
        self.assertEqual(round(ell_dist, 3), 4831.553)
        self.assertEqual(round(dec2hp(azimuth1to2), 6), 90.004480)
        self.assertEqual(round(dec2hp(azimuth2to1), 6), 269.591520)

        test2 = vincinv(lat1, lon1, lat1, lon1)
        self.assertEqual(test2, (0, 0, 0))

        # Test coincident coordinates (within float precision)
        pl1 = (-30.645230675, 152.996285475)
        pl2 = (-30.645230674999997, 152.996285475)
        test3 = vincinv(pl1[0], pl1[1], pl2[0], pl2[1])
        self.assertEqual(test3, (0, 0, 0))

        # Test DMSAngle Input
        ell_dist, azimuth1to2, azimuth2to1 = vincinv(
            lat1_DMS, lon1_DMS, lat2_DMS, lon2_DMS
        )
        self.assertEqual(round(ell_dist, 3), 54972.271)
        self.assertEqual(round(dec2hp(azimuth1to2), 6), 306.520537)
        self.assertEqual(round(dec2hp(azimuth2to1), 6), 127.102507)

        test2 = vincinv(lat1_DMS, lon1_DMS, lat1_DMS, lon1_DMS)
        self.assertEqual(test2, (0, 0, 0))

        # Test DDMAngle Input
        (ell_dist, azimuth1to2, azimuth2to1) = vincinv(
            lat1_DMS.ddm(), lon1_DMS.ddm(), lat2_DMS.ddm(), lon2_DMS.ddm()
        )
        self.assertEqual(round(ell_dist, 3), 54972.271)
        self.assertEqual(round(dec2hp(azimuth1to2), 6), 306.520537)
        self.assertEqual(round(dec2hp(azimuth2to1), 6), 127.102507)

        test2 = vincinv(lat1_DMS.ddm(), lon1_DMS.ddm(), lat1_DMS.ddm(), lon1_DMS.ddm())
        self.assertEqual(test2, (0, 0, 0))

    def test_vincdir(self):
        # Flinders Peak
        lat1 = hp2dec(-37.57037203)
        lon1 = hp2dec(144.25295244)
        lat1_DMS = DMSAngle(-37, 57, 3.7203)
        lon1_DMS = DMSAngle(144, 25, 29.5244)

        # To Buninyong
        azimuth1to2 = hp2dec(306.520537)
        azimuth1to2_DMS = DMSAngle(306, 52, 5.37)
        ell_dist = 54972.271

        # Test Decimal Degrees Input
        lat2, lon2, azimuth2to1 = vincdir(lat1, lon1, azimuth1to2, ell_dist)
        self.assertEqual(round(dec2hp(lat2), 8), -37.39101561)
        self.assertEqual(round(dec2hp(lon2), 8), 143.55353839)
        self.assertEqual(round(dec2hp(azimuth2to1), 6), 127.102507)

        # Test DMSAngle Input
        lat2, long2, azimuth2to1 = vincdir(
            lat1_DMS, lon1_DMS, azimuth1to2_DMS, ell_dist
        )
        self.assertEqual(round(dec2hp(lat2), 8), -37.39101561)
        self.assertEqual(round(dec2hp(long2), 8), 143.55353839)
        self.assertEqual(round(dec2hp(azimuth2to1), 6), 127.102507)

        # Test DDMAngle Input
        lat2, long2, azimuth2to1 = vincdir(
            lat1_DMS.ddm(), lon1_DMS.ddm(), azimuth1to2_DMS.ddm(), ell_dist
        )
        self.assertEqual(round(dec2hp(lat2), 8), -37.39101561)
        self.assertEqual(round(dec2hp(long2), 8), 143.55353839)
        self.assertEqual(round(dec2hp(azimuth2to1), 6), 127.102507)

    def test_vincinv_utm(self):
        # Flinders Peak (UTM 55)
        zone1 = 55
        east1 = 273741.2966
        north1 = 5796489.7769
        # Buninyong (UTM 55)
        zone2 = 55
        east2 = 228854.0513
        north2 = 5828259.0384
        # Buninyong (UTM 54)
        zone3 = 54
        east3 = 758173.7973
        north3 = 5828674.3402

        # Test Coordinates in Zone 55 only
        grid_dist, grid1to2, grid2to1, lsf = vincinv_utm(
            zone1, east1, north1, zone2, east2, north2
        )
        self.assertAlmostEqual(lsf, 1.00036397, 8)
        self.assertAlmostEqual(grid_dist, 54992.279, 3)
        self.assertAlmostEqual(dec2hp(grid1to2), 305.17017259, 7)
        self.assertAlmostEqual(dec2hp(grid2to1), 125.17418518, 7)

        # Test Coordinates in Different Zones (55 and 54)
        # (Point 2 Grid Bearing Different (Zone 54 Grid Bearing))
        grid_dist, grid1to2, grid2to1, lsf = vincinv_utm(
            zone1, east1, north1, zone3, east3, north3
        )
        self.assertAlmostEqual(lsf, 1.00036397, 8)
        self.assertAlmostEqual(grid_dist, 54992.279, 3)
        self.assertAlmostEqual(dec2hp(grid1to2), 305.17017259, 7)
        self.assertAlmostEqual(dec2hp(grid2to1), 128.57444307, 7)

    def test_vincdir_utm(self):
        # Flinders Peak (UTM 55)
        zone1 = 55
        east1 = 273741.2966
        north1 = 5796489.7769
        # Grid Dimensions to Point 2 (Buninyong)
        grid_dist = 54992.279
        grid1to2 = hp2dec(305.17017259)
        grid1to2_DMS = DMSAngle(305, 17, 1.7259)

        # Test Decimal Degrees Input
        (zone2, east2, north2, grid2to1, lsf) = vincdir_utm(
            zone1, east1, north1, grid1to2, grid_dist
        )
        self.assertEqual(zone2, zone1)
        self.assertAlmostEqual(east2, 228854.0513, 3)
        self.assertAlmostEqual(north2, 5828259.0384, 3)
        self.assertAlmostEqual(dec2hp(grid2to1), 125.17418518, 7)
        self.assertAlmostEqual(lsf, 1.00036397, 8)

        # Test DMSAngle Input
        (zone2, east2, north2, grid2to1, lsf) = vincdir_utm(
            zone1, east1, north1, grid1to2_DMS, grid_dist
        )
        self.assertEqual(zone2, zone1)
        self.assertAlmostEqual(east2, 228854.0513, 3)
        self.assertAlmostEqual(north2, 5828259.0384, 3)
        self.assertAlmostEqual(dec2hp(grid2to1), 125.17418518, 7)
        self.assertAlmostEqual(lsf, 1.00036397, 8)

    def test_equality_vincentys(self):
        # Test multiple point-to-point vincinv calculations
        abs_path = os.path.abspath(os.path.dirname(__file__))

        test_geo_coords = np.genfromtxt(
            os.path.join(abs_path, "resources/Test_Conversion_Geo.csv"),
            delimiter=",",
            dtype="S4,f8,f8",
            names=["site", "lat1", "long1"],
            usecols=("lat1", "long1"),
        )

        test_geo_coord2 = np.genfromtxt(
            os.path.join(abs_path, "resources/Test_Conversion_Geo.csv"),
            delimiter=",",
            dtype="S4,f8,f8",
            names=["site", "lat2", "long2"],
            usecols=("lat2", "long2"),
        )

        # Form array with point pairs from test file
        test_pairs = rfn.merge_arrays(
            [test_geo_coords, np.roll(test_geo_coord2, 1)], flatten=True
        )

        # Calculate Vincenty's Inverse Result using Lat, Long Pairs
        vincinv_result = np.array(
            list(vincinv(*x) for x in test_pairs[["lat1", "long1", "lat2", "long2"]])
        )

        # Calculate Vincenty's Direct Result using Results from Inverse Function
        vincdir_input = rfn.merge_arrays(
            [test_geo_coords, vincinv_result[:, 1], vincinv_result[:, 0]], flatten=True
        )
        vincdir_input.dtype.names = ["lat1", "long1", "az1to2", "ell_dist"]
        vincdir_result = np.array(
            list(
                vincdir(*x)
                for x in vincdir_input[["lat1", "long1", "az1to2", "ell_dist"]]
            )
        )

        np.testing.assert_almost_equal(
            test_pairs["lat2"], vincdir_result[:, 0], decimal=8
        )
        np.testing.assert_almost_equal(
            test_pairs["long2"], vincdir_result[:, 1], decimal=8
        )
        np.testing.assert_almost_equal(vincinv_result[:, 2], vincdir_result[:, 2])

    def test_vincinv_edgecases(self):
        lat1 = -32.153892
        lon1 = -15.394827
        lat2 = -31.587369
        lon2 = -13.487739
        gdist, az12, az21 = vincinv(lat1, lon1, lat2, lon2)
        lon1 = lon1 + 14
        lon2 = lon2 + 14
        gdist_2, az12_2, az21_2 = vincinv(lat1, lon1, lat2, lon2)
        self.assertEqual(gdist, gdist_2)
        self.assertEqual(az12, az12_2)
        self.assertEqual(az21, az21_2)


if __name__ == "__main__":
    unittest.main()
