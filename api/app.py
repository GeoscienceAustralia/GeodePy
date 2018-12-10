from flask import Flask
from flask import request
from flask import url_for
from geodepy.geodesy import vincdir
from geodepy.geodesy import vincinv
from geodepy.convert import dms2dd
from geodepy.convert import dd2dms

app = Flask(__name__)

coord2dd = {
    'dd': lambda x: x,
    'dms': dms2dd,
}

dd2coord = {
    'dd': lambda x: x,
    'dms': dd2dms,
}


@app.route('/')
def list_routes():
    return str(tuple(url_for(rule.endpoint) for rule in app.url_map.iter_rules() if rule.endpoint != 'static'))


@app.route('/vincinv')
def handle_vincinv():
    from_coord = request.args.get('from_coord')
    to_coord = request.args.get('to_coord')
    lat1 = request.args.get('lat1', type=float)
    lon1 = request.args.get('lon1', type=float)
    lat2 = request.args.get('lat2', type=float)
    lon2 = request.args.get('lon2', type=float)

    dd = coord2dd[from_coord]
    lat1_dd, lon1_dd, lat2_dd, lon2_dd = dd(lat1), dd(lon1), dd(lat2), dd(lon2)

    ell_dist_dd, azimuth1to2_dd, azimuth2to1_dd = vincinv(lat1_dd, lon1_dd, lat2_dd, lon2_dd)

    coord = dd2coord[to_coord]
    ell_dist, azimuth1to2, azimuth2to1 = coord(ell_dist_dd), coord(azimuth1to2_dd), coord(azimuth2to1_dd)

    return str((ell_dist, azimuth1to2, azimuth2to1)), 200


@app.route('/vincdir')
def handle_vincdir():
    from_coord = request.args.get('from_coord')
    to_coord = request.args.get('to_coord')
    lat1 = request.args.get('lat1', type=float)
    lon1 = request.args.get('lon1', type=float)
    azimuth1to2 = request.args.get('azimuth1to2', type=float)
    ell_dist = request.args.get('ell_dist', type=float)

    dd = coord2dd[from_coord]
    lat1_dd, lon1_dd, azimuth1to2_dd = dd(lat1), dd(lon1), dd(azimuth1to2)

    lat2_dd, lon2_dd, azimuth2to1_dd = vincdir(lat1_dd, lon1_dd, azimuth1to2_dd, ell_dist)

    coord = dd2coord[to_coord]
    lat2, lon2, azimuth2to1 = coord(lat2_dd), coord(lon2_dd), coord(azimuth2to1_dd)

    return str((lat2, lon2, azimuth2to1)), 200


if __name__ == '__main__':
    app.run()
