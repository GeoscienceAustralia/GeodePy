from flask import Flask
from flask import jsonify
from flask import request
from flask import url_for
from geodepy.geodesy import vincdir
from geodepy.geodesy import vincinv
from geodepy.convert import hp2dec
from geodepy.convert import dec2hp


app = Flask(__name__)

angle_type_to_dd = {
    'dd': lambda x: x,
    'dms': hp2dec,
}

dd_to_angle_type = {
    'dd': lambda x: x,
    'dms': dec2hp,
}


@app.route('/')
def list_routes():
    return str(tuple(url_for(rule.endpoint) for rule in
                     app.url_map.iter_rules() if rule.endpoint != 'static'))


@app.route('/vincinv')
def handle_vincinv():
    from_angle_type = request.args.get('from_angle_type', default='dd')
    to_angle_type = request.args.get('to_angle_type', default='dd')
    lat1 = request.args.get('lat1', type=float)
    lon1 = request.args.get('lon1', type=float)
    lat2 = request.args.get('lat2', type=float)
    lon2 = request.args.get('lon2', type=float)

    dd = angle_type_to_dd[from_angle_type]
    lat1_dd, lon1_dd, lat2_dd, lon2_dd = dd(lat1), dd(lon1), dd(lat2), dd(lon2)

    ell_dist, azimuth1to2_dd, azimuth2to1_dd = vincinv(lat1_dd, lon1_dd,
                                                       lat2_dd, lon2_dd)

    angle = dd_to_angle_type[to_angle_type]
    azimuth1to2, azimuth2to1 = angle(azimuth1to2_dd), angle(azimuth2to1_dd)

    return jsonify({
        'ell_dist': ell_dist,
        'azimuth1to2': azimuth1to2,
        'azimuth2to1': azimuth2to1
    }), 200


@app.route('/vincdir')
def handle_vincdir():
    from_angle_type = request.args.get('from_angle_type', default='dd')
    to_angle_type = request.args.get('to_angle_type', default='dd')
    lat1 = request.args.get('lat1', type=float)
    lon1 = request.args.get('lon1', type=float)
    azimuth1to2 = request.args.get('azimuth1to2', type=float)
    ell_dist = request.args.get('ell_dist', type=float)

    dd = angle_type_to_dd[from_angle_type]
    lat1_dd, lon1_dd, azimuth1to2_dd = dd(lat1), dd(lon1), dd(azimuth1to2)

    lat2_dd, lon2_dd, azimuth2to1_dd = vincdir(lat1_dd, lon1_dd,
                                               azimuth1to2_dd, ell_dist)

    angle_type = dd_to_angle_type[to_angle_type]
    lat2, lon2, azimuth2to1 = angle_type(lat2_dd), angle_type(lon2_dd),\
        angle_type(azimuth2to1_dd)

    return jsonify({
        'lat2': lat2,
        'lon2': lon2,
        'azimuth2to1': azimuth2to1
    }), 200


if __name__ == '__main__':
    app.run()
