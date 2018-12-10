import unittest
from api.app import app


class TestAPI(unittest.TestCase):
    def test_vincinv(self):
        query = {
            'lat1': -37.57037203,
            'lon1': 144.25295244,
            'lat2': -37.39101561,
            'lon2': 143.5535383,
            'from_coord': 'dms',
            'to_coord': 'dms'
        }

        expected_response = b'(54972.17204, 306.52053231124, 127.10250207968)'
        response = app.test_client().get('/vincinv', query_string=query)
        self.assertEqual(response.data, expected_response)

    def test_vincdir(self):
        query = {
            'lat1': -37.57037203,
            'lon1': 144.25295244,
            'azimuth1to2': 306.520537,
            'ell_dist': 54972.271,
            'from_coord': 'dms',
            'to_coord': 'dms'
        }

        expected_response = b'(-37.3910156124268, 143.5535383883988, 127.10250671432)'
        response = app.test_client().get('/vincdir', query_string=query)
        self.assertEqual(response.data, expected_response)


