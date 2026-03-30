import json
from datetime import date
import geodepy.transform as tr
import geodepy.constants as gc
import geodepy.geodesy as gg
import geodepy.convert as con
import numpy as np
import datetime
import geodepy.point_in_polygon as pp

# Setup for transformation paths
with open("other_files/transformation_routes_v3.json", "r", encoding="utf-8") as f:
    ROUTING = json.load(f)

ROUTES = ROUTING["routes"]
DIRECT = ROUTING["direct_pairs"]

def mga_parse(ref):
    """
    If referenc frame is MGA, change to equivalent GDA.
    MGA conversions are handled individually in functions.

    :param ref: Reference frame
    :type ref: str
    :return: Reference frame changed if needed
    :rtype: str
    """
    mga_names = {"MGA94": "GDA94",
                "MGA2020": "GDA2020"}
    
    if ref in mga_names:
        ref = mga_names[ref]

    return ref

def resolve_path(src, dst):
    """
    Find path between source and destination reference paths from json file.

    :param src: Source reference frame
    :type src: str
    :param dst: Destination reference frame
    :type dst: str
    :return: Return path from source to destination as array
    :rtype: list of str
    """
    if src == dst:
        return [src]

    # Change MGA to equivalent GDA
    src = mga_parse(src)
    dst = mga_parse(dst)

    key = f"{src}|{dst}"
    if key in ROUTES:
        return ROUTES[key]

    rev_key = f"{dst}|{src}"
    if rev_key in ROUTES:
        return list(reversed(ROUTES[rev_key]))

    if key in DIRECT:
        return [src, dst]

    raise KeyError(f"No route found for {src} -> {dst}")

def ref_frame_parser(ref_name):
    """
    Return GeodePy compliant refernece frame name from generic name.

    :param ref_name: Reference frame name
    :type ref_name: str
    :return: Return GeodePy compliant reference frame name
    :rtype: str
    """

    # Dictionary of reference frame names
    ref_names = {"wgs84 (g2296)": "wgs84g2296",
                 "wgs84 (g2139)": "wgs84g2139",
                 "wgs84 (g1762)": "wgs84g1762",
                 "wgs84 (g1674)": "wgs84g1674",
                 "wgs84 (g1150)": "wgs84g1150",
                 "wgs84 (g873)": "wgs84g873",
                 "wgs84 (g730)": "wgs84g730",
                 "wgs84 (transit)": "wgs84trans",
                 "wgs84 ensemble": "wgs84ensemble"}

    ref_name = ref_name.lower()

    if ref_name in ref_names:
        ref_name = ref_names[ref_name]

    return ref_name

def transformation_type(path):
    """
    Find the transformation types. The options are:
    1. static_to_static
    2. static_to_dynamic
    3. dynamic_to_static
    4. dynamic_to_dynamic

    :param path: The path between reference frames
    :type path: list of str
    :return: Returns string with one of four options
    :rtype: str
    """

    # Reference frame types dictionary
    trans_types = {  "GDA2020": "static",
    "GDA94": "static",
    "AGD66": "static",
    "AGD84": "static",
    "ATRF2014": "dynamic",
    "ITRF88": "dynamic",
    "ITRF89": "dynamic",
    "ITRF90": "dynamic",
    "ITRF91": "dynamic",
    "ITRF92": "dynamic",
    "ITRF93": "dynamic",
    "ITRF94": "dynamic",
    "ITRF96": "dynamic",
    "ITRF97": "dynamic",
    "ITRF2000": "dynamic",
    "ITRF2005": "dynamic",
    "ITRF2008": "dynamic",
    "ITRF2014": "dynamic",
    "ITRF2020": "dynamic",
    "WGS84 (G2296)": "dynamic",
    "WGS84 (G2139)": "dynamic",
    "WGS84 (G1762)": "dynamic",
    "WGS84 (G1674)": "dynamic",
    "WGS84 (G1150)": "dynamic",
    "WGS84 (G873)": "dynamic",
    "WGS84 (G730)": "dynamic",
    "WGS84 (Transit)": "dynamic",
    "WGS84 Ensemble": "dynamic",
    }

    start_type = trans_types[path[0]]
    end_type = trans_types[path[-1]]

    return f"{start_type}_to_{end_type}"

def static_to_static_trans(x, y, z, path, vcv):
    """
    Completes static to static transformation

    :param x: Cartesian X (m)
    :type x: float
    :param y: Cartesian Y (m)
    :type y: float
    :param z: Cartesian Z (m)
    :type z: float
    :param path: Array with path between reference frames
    :type path: list of str
    :param vcv: Optional 3*3 numpy array in Cartesian units to propagate tf uncertainty
    :type vcv: numpy.ndarray
    :return: Cartesian X, Y, Z co-ordinates and vcv matrix transformed
    :rtype: tuple
    """
    
    i = 1

    # Completes each step of transformation path
    while i < len(path):

        from_frame = ref_frame_parser(path[i-1])
        to_frame = ref_frame_parser(path[i])
        
        # Find name of correct transformation
        transform = f"{from_frame}_to_{to_frame}"
        print(transform)

        # Get transformation from geodepy.constants
        trans = getattr(gc, transform)

        # Complete 7 paramater transformation
        x, y, z, vcv = tr.conform7(x, y, z, trans, vcv)

        i+=1


    print("")
    print("Completed Transformation")
    return round(x,4), round(y, 4), round(z, 4), vcv

def static_to_dynamic_trans(x, y, z, to_epoch, path, vcv):
    """
    Completes static to dynamic transformation

    :param x: Cartesian X (m)
    :type x: float
    :param y: Cartesian Y (m)
    :type y: float
    :param z: Cartesian Z (m)
    :type z: float
    :param to_epoch: The target epoch for the dynamic reference frame (datetime.date Object)
    :type to_epoch: datetime.date
    :param path: Array with path between reference frames
    :type path: list of str
    :param vcv: Optional 3*3 numpy array in Cartesian units to propagate tf uncertainty
    :type vcv: numpy.ndarray
    :return: Cartesian X, Y, Z co-ordinates and vcv matrix transformed
    :rtype: tuple
    """

    i = 1

    # Completes each step of transformation path
    while i < len(path):

        from_frame = ref_frame_parser(path[i-1])
        to_frame = ref_frame_parser(path[i])

        # Find name of correct transformation
        transform = f"{from_frame}_to_{to_frame}"
        print(transform)

        # Get transformation from geodepy.constants
        trans = getattr(gc, transform)

        # If transformation is for 7 parameter complete conform7 otherwise conform14
        if trans.ref_epoch == 0:

            x, y, z, vcv = tr.conform7(x, y, z, trans, vcv)

        else:
            x, y, z, vcv = tr.conform14(x, y, z, to_epoch, trans, vcv)

        i+=1


    print("")
    print("Completed Transformation")
    return round(x,4), round(y, 4), round(z, 4), vcv

def dynamic_to_static_trans(x, y, z, from_epoch, path, vcv):
    """
    Completes dynamic to static transformation

    :param x: Cartesian X (m)
    :type x: float
    :param y: Cartesian Y (m)
    :type y: float
    :param z: Cartesian Z (m)
    :type z: float
    :param from_epoch: The source epoch for the dynamic reference frame (datetime.date Object)
    :type from_epoch: datetime.date
    :param path: Array with path between reference frames
    :type path: list of str
    :param vcv: Optional 3*3 numpy array in Cartesian units to propagate tf uncertainty
    :type vcv: numpy.ndarray
    :return: Cartesian X, Y, Z co-ordinates and vcv matrix transformed
    :rtype: tuple
    """

    i = 1

    # Completes each step of transformation path
    while i < len(path):

        from_frame = ref_frame_parser(path[i-1])
        to_frame = ref_frame_parser(path[i])

        # Find name of correct transformation
        transform = f"{from_frame}_to_{to_frame}"
        print(transform)

        # Get transformation from geodepy.constants
        trans = getattr(gc, transform)

        # If transformation is for 7 parameter complete conform7 otherwise conform14
        if trans.ref_epoch == 0:

            x, y, z, vcv = tr.conform7(x, y, z, trans, vcv)

        else:
            x, y, z, vcv = tr.conform14(x, y, z, from_epoch, trans, vcv)

        i+=1


    print("")
    print("Completed Transformation")
    return round(x,4), round(y, 4), round(z, 4), vcv

def dynamic_to_dynamic_trans(x, y, z, from_ref, to_ref, from_epoch, to_epoch, plate_motion, vcv):
    """
    Completes dynamic to dynamic transformation

    :param x: Cartesian X (m)
    :type x: float
    :param y: Cartesian Y (m)
    :type y: float
    :param z: Cartesian Z (m)
    :type z: float
    :param from_ref: Source reference frame
    :type from_ref: str
    :param to_ref: Target reference frame
    :type to_ref: str
    :param from_epoch: The source epoch for the dynamic reference frame (datetime.date Object)
    :type from_epoch: datetime.date
    :param to_epoch: The target epoch for the dynamic reference frame (datetime.date Object)
    :type to_epoch: datetime.date
    :param plate_motion: The plate motion model to use
    :type plate_motion: str
    :param vcv: Optional 3*3 numpy array in Cartesian units to propagate tf uncertainty
    :type vcv: numpy.ndarray
    :return: Cartesian X, Y, Z co-ordinates and vcv matrix transformed
    :rtype: tuple
    """

    # If from and to epoch are the same then direct transformation without plate motion 
    # can be completed
    if from_epoch == to_epoch:
        path = resolve_path(from_ref, to_ref)

        i = 1

        # Completes each step of transformation path
        while i < len(path):

            from_frame = ref_frame_parser(path[i-1])
            to_frame = ref_frame_parser(path[i])

            # Find name of correct transformation
            transform = f"{from_frame}_to_{to_frame}"
            print(transform)

            # Get transformation from geodepy.constants
            trans = getattr(gc, transform)

            # If transformation is for 7 parameter complete conform7 otherwise conform14
            if trans.ref_epoch == 0:

                x, y, z, vcv = tr.conform7(x, y, z, trans, vcv)

            else:
                x, y, z, vcv = tr.conform14(x, y, z, from_epoch, trans, vcv)

            i+=1


        print("")
        print("Completed Transformation")
        return round(x,4), round(y, 4), round(z, 4), vcv

    # Otherwise if epochs are not the same and plate motion is required

    # Create paths to ITRF2014 and from ITRF2014, allows for plate motion on ITRF2014
    path_pre = resolve_path(from_ref, "ITRF2014")
    path_post = resolve_path("ITRF2014", to_ref)

    i = 1

    # Convert to ITRF2014 coord, completing each step
    while i < len(path_pre):

        from_frame = ref_frame_parser(path_pre[i-1])
        to_frame = ref_frame_parser(path_pre[i])

        # Find name of correct transformation
        transform = f"{from_frame}_to_{to_frame}"
        print(transform)

        # Get transformation from geodepy.constants
        trans = getattr(gc, transform)

        # If transformation is for 7 parameter complete conform7 otherwise conform14
        if trans.ref_epoch == 0:

            x, y, z, vcv = tr.conform7(x, y, z, trans, vcv)

        else:
            x, y, z, vcv = tr.conform14(x, y, z, from_epoch, trans, vcv)

        i+=1

    # Decide which plate motion to use
    if plate_motion == "auto":
        x, y, z, vcv = pp.universal_plate_motion_transformation(x, y, z, from_epoch, to_epoch, vcv)

    elif plate_motion == "aus":
        x, y, z, vcv = tr.plate_motion_transformation(x, y, z, from_epoch, to_epoch, gc.itrf2014_to_gda2020, vcv)

    else:
        raise ValueError("plate_motion must be either auto or aus")

    # Convert to final ref frame

    i = 1

    # convert from ITRF2014 coord, completing each step
    while i < len(path_post):

        from_frame = ref_frame_parser(path_post[i-1])
        to_frame = ref_frame_parser(path_post[i])

        # Find name of correct transformation
        transform = f"{from_frame}_to_{to_frame}"
        print(transform)

        # Get transformation from geodepy.constants
        trans = getattr(gc, transform)

        # If transformation is for 7 parameter complete conform7 otherwise conform14
        if trans.ref_epoch == 0:

            x, y, z, vcv = tr.conform7(x, y, z, trans, vcv)

        else:
            x, y, z, vcv = tr.conform14(x, y, z, to_epoch, trans, vcv)

        i+=1

    print("")
    print("Completed Transformation")
    return round(x,4), round(y, 4), round(z, 4), vcv

def universal_transform(x, y, z, from_ref, to_ref, from_epoch=None, to_epoch=None, plate_motion="auto", vcv=None, return_type="xyz", ignore_errors=False):
    """
    Completes a transformation from any reference frame to any other

    :param x: Cartesian X (m)
    :type x: float
    :param y: Cartesian Y (m)
    :type y: float
    :param z: Cartesian Z (m)
    :type z: float
    :param from_ref: Source reference frame
    :type from_ref: str
    :param to_ref: Target reference frame
    :type to_ref: str
    :param from_epoch: The source epoch for the dynamic reference frame (datetime.date Object)
    :type from_epoch: datetime.date
    :param to_epoch: The target epoch for the dynamic reference frame (datetime.date Object)
    :type to_epoch: datetime.date
    :param plate_motion: Set what plate motion to use: "auto" Automatically finds plate motion, "aus" use australian plate motion model.
    :type plate_motion: str
    :param vcv: Optional 3*3 numpy array in Cartesian units to propagate tf uncertainty
    :type vcv: numpy.ndarray
    :param return_type: The type of coordinates to return, either "xyz", "llh" or "enu"
    :type return_type: str
    :param ignore_errors: If True, will not raise errors about point not being on Australia's plate
    :type ignore_errors: bool
    :return: Co-ordinates and vcv matrix transformed and in the return type given
    :rtype: dict
    """
    # MGA can not be enters as from_ref for xyz
    if from_ref.lower() in ("mga94", "mga2020"):
        raise ValueError("Source reference frame can not be MGA when using cartesian coordinates (xyz).")

    # Find path for transformation and check frames given are real
    path = resolve_path(from_ref, to_ref)

    # If path exists find transformation type
    trans_type = transformation_type(path)

    # Input validation of XYZ
    if type(x) != float:
        raise TypeError("X value not float")
    if type(y) != float:
        raise TypeError("Y value not float")
    if type(z) != float:
        raise TypeError("Z value not float")

    # Input validation of from_epoch and to_epoch
    if trans_type == "static_to_dynamic" and type(to_epoch) != datetime.date:
        raise ValueError("to_epoch is required for static to dynamic transformation and must be a datetime.date object")
    if trans_type == "dynamic_to_static" and type(from_epoch) != datetime.date:
        raise ValueError("from_epoch is required for dynamic to static transformation and must be a datetime.date object") 
    if trans_type == "dynamic_to_dynamic" and (type(to_epoch) != datetime.date or type(from_epoch) != datetime.date):
        raise ValueError("to_epoch and from_epoch are required for dynamic to dynamic transformation and must be datetime.date objects")       
    
    # Input validation of plate_Motion
    if plate_motion not in ("auto", "aus"):
        raise ValueError("plate_motion must be either 'auto' or 'aus'")
    
    if vcv is not None:
        if not isinstance(vcv, np.ndarray):
            raise TypeError("VCV must be a numpy array")
        if vcv.shape != (3, 3):
            raise ValueError(f"VCV must be shape (3,3), got {vcv.shape}")
    
    # If using Australian plate motion ("aus"), check point is on australian plate
    plates = pp.build_plate_index("other_files/MORVEL56_plates.dig")
    plate_id = pp.plate_from_xyz(x, y, z, plates)

    if plate_motion == "aus" and plate_id != "AU" and not ignore_errors:
        raise ValueError("Point is not on the australia plate. Can not use Australian plate motion.")
    
    # If using any australian reference frame, check point is within Australia's EEZ
    aus_frames = ("GDA94", "GDA2020", "AGD66", "AGD84", "ATRF2020")
    
    if (from_ref in aus_frames) or (to_ref in aus_frames):
        plates = pp.build_plate_index("other_files/EEZ_australia_approx.dig")
        plate_id = pp.plate_from_xyz(x, y, z, plates)

        if plate_id is None and not ignore_errors:
            raise ValueError("Point is not within Australia's EEZ while using Austalian refernce frame.")

    wgs84_frames = ["WGS84 (G2296)",
                    "WGS84 (G2139)",
                    "WGS84 (G1762)",
                    "WGS84 (G1674)",
                    "WGS84 (G1150)",
                    "WGS84 (G873)",
                    "WGS84 (G730)",
                    "WGS84 (Transit)",
                    "WGS84 Ensemble"]

    # Check for WGS84 and change date to middle of year
    if from_ref in wgs84_frames:

        from_epoch = date(from_epoch.year, 6, 30)

    if to_ref in wgs84_frames:

        to_epoch = date(to_epoch.year, 6, 30)

    # Check return_type is an option and correct for to_ref
    if return_type.lower() not in ("xyz", "llh", "enu"):
        raise ValueError("return_type must be either \"xyz\", \"llh\", \"enu\"")
    
    if (return_type.lower() == "enu") and (to_ref.lower() not in ("mga94", "mga2020")):
        raise ValueError("Can only give an enu result if to_ref is an MGA reference frame.")
    
    if (return_type.lower() != "enu") and (to_ref.lower() in ("mga94", "mga2020")):
        raise ValueError("Can only return mga coordinates in enu. return_type must be \"enu\"")

    #Finished input validation
    print()
    print(f"Completing a {trans_type} transformation")
    print("")
    print("Steps: ")

    # Choose the correct function for the transformation type
    if trans_type == "static_to_static":
        x, y, z, vcv = static_to_static_trans(x, y, z, path, vcv)

    elif trans_type == "static_to_dynamic":
        x, y, z, vcv = static_to_dynamic_trans(x, y, z, to_epoch, path, vcv)

    elif trans_type == "dynamic_to_static":
        x, y, z, vcv = dynamic_to_static_trans(x, y, z, from_epoch, path, vcv)

    elif trans_type == "dynamic_to_dynamic":
        x, y, z, vcv = dynamic_to_dynamic_trans(x, y, z, from_ref, to_ref, from_epoch, to_epoch, plate_motion, vcv)

    else:
        raise(ValueError("Transformation type not valid"))

    # Assign output to dictionary

    output = {
        "type": return_type,
        "vcv": vcv,
    }

    if return_type.lower() == "xyz":
        output["coords"] = {
            "x": round(x, 4),
            "y": round(y, 4),
            "z": round(z, 4)
        }
    elif return_type.lower() == "llh":
        lat, lon, el_height = con.xyz2llh(x, y, z)
        output["coords"] = {
            "lat": round(lat,8),
            "lon": round(lon,8),
            "el_height": round(el_height,4)
        }
    elif return_type.lower() == "enu":
        lat, lon, el_height = con.xyz2llh(x, y, z)        
        hem, zone, east, north, psf, converge = con.geo2grid(lat, lon)
        output["coords"] = {
            "east": east,
            "north": north,
            "height": el_height,
            "zone": zone
        }
    
    return output

def universal_transform_llh(lat, lon, el_height, from_ref, to_ref, from_epoch=None, to_epoch=None, plate_motion="auto", vcv=None, return_type="llh", ignore_errors=False):
    """
    Completes a transformation from any reference frame to any other

    :param lat: Latitude (degrees)
    :type lat: float
    :param lon: Longitude (degrees)
    :type lon: float
    :param el_height: Ellipsoidal height (m)
    :type el_height: float
    :param from_ref: Source reference frame
    :type from_ref: str
    :param to_ref: Target reference frame
    :type to_ref: str
    :param from_epoch: The source epoch for the dynamic reference frame (datetime.date Object)
    :type from_epoch: datetime.date
    :param to_epoch: The target epoch for the dynamic reference frame (datetime.date Object)
    :type to_epoch: datetime.date
    :param plate_motion: Set what plate motion to use: "auto" Automatically finds plate motion, "aus" use australian plate motion model.
    :type plate_motion: str
    :param vcv: Optional 3*3 numpy array in Cartesian units to propagate tf uncertainty
    :type vcv: numpy.ndarray
    :param return_type: The type of coordinates to return, either "xyz", "llh" or "enu"
    :type return_type: str
    :param ignore_errors: If True, will not raise errors about point not being on Australia's plate
    :type ignore_errors: bool
    :return: Co-ordinates and vcv matrix transformed and in the return type given
    :rtype: dict
    """
    # Can not have mga as from_ref for llh input
    if from_ref.lower() in ("mga94", "mga2020"):
        raise ValueError("Source reference frame can not be MGA when using geographic coordinates (llh).")

    # Input validation of lat, lon, el_height
    if type(lat) != float:
        raise TypeError("lat value not float")
    if type(lon) != float:
        raise TypeError("lon value not float")
    if type(el_height) != float:
        raise TypeError("el_height value not float") 

    # Check return_type is an option and correct for to_ref
    if return_type.lower() not in ("xyz", "llh", "enu"):
        raise ValueError("return_type must be either \"xyz\", \"llh\", \"enu\"")
    
    if (return_type.lower() == "enu") and (to_ref.lower() not in ("mga94", "mga2020")):
        raise ValueError("Can only give an enu result if to_ref is an MGA reference frame.")
    
    if (return_type.lower() != "enu") and (to_ref.lower() in ("mga94", "mga2020")):
        raise ValueError("Can only return mga coordinates in enu. return_type must be \"enu\"")

    # Finished input validation

    # Convert to xyz and run xyz function
    x, y, z = con.llh2xyz(lat, lon, el_height)

    # Change to_ref to GDA equivalent
    to_ref = mga_parse(to_ref)

    result = universal_transform(x, y, z, from_ref, to_ref, from_epoch, to_epoch, plate_motion, vcv, ignore_errors=ignore_errors)

    x = result["coords"]["x"]
    y = result["coords"]["y"]
    z = result["coords"]["z"]
    vcv = result["vcv"]

    # Assign output to dictionary

    output = {
        "type": return_type,
        "vcv": vcv,
    }

    if return_type.lower() == "xyz":
        output["coords"] = {
            "x": round(x, 4),
            "y": round(y, 4),
            "z": round(z, 4)
        }
    elif return_type.lower() == "llh":
        lat, lon, el_height = con.xyz2llh(x, y, z)
        output["coords"] = {
            "lat": round(lat,8),
            "lon": round(lon,8),
            "el_height": round(el_height,4)
        }
    elif return_type.lower() == "enu":
        lat, lon, el_height = con.xyz2llh(x, y, z)        
        hem, zone, east, north, psf, converge = con.geo2grid(lat, lon)
        output["coords"] = {
            "east": east,
            "north": north,
            "height": el_height,
            "zone": zone
        }
    
    return output
    
def universal_transform_enu(east, north, height, zone, from_ref, to_ref, from_epoch=None, to_epoch=None, plate_motion="auto", vcv=None, return_type="enu", ignore_errors=False):
    """
    Completes a transformation from any reference frame to any other

    :param east: Easting (m)
    :type east: float
    :param north: Northing (m)
    :type north: float
    :param height: Ellipsoidal height (m)
    :type height: float
    :param zone: UTM zone
    :type zone: int
    :param from_ref: Source reference frame
    :type from_ref: str
    :param to_ref: Target reference frame
    :type to_ref: str
    :param from_epoch: The source epoch for the dynamic reference frame (datetime.date Object)
    :type from_epoch: datetime.date
    :param to_epoch: The target epoch for the dynamic reference frame (datetime.date Object)
    :type to_epoch: datetime.date
    :param plate_motion: Set what plate motion to use: "auto" Automatically finds plate motion, "aus" use australian plate motion model.
    :type plate_motion: str
    :param vcv: Optional 3*3 numpy array in Cartesian units to propagate tf uncertainty
    :type vcv: numpy.ndarray
    :param return_type: The type of coordinates to return, either "xyz", "llh" or "enu"
    :type return_type: str
    :param ignore_errors: If True, will not raise errors about point not being on Australia's plate
    :type ignore_errors: bool
    :return: Co-ordinates and vcv matrix transformed and in the return type given
    :rtype: dict
    """
    # Can not have mga as from_ref for llh input
    if from_ref.lower() not in ("mga94", "mga2020"):
        raise ValueError("Source reference frame can only be MGA when using projected coordinates (enu). from_ref must be \"MGA94\" or \"MGA2020\".")

    # Input validation of lat, lon, el_height
    if type(east) != float:
        raise TypeError("east value not float")
    if type(north) != float:
        raise TypeError("north value not float")
    if type(height) != float:
        raise TypeError("height value not float") 
    if type(zone) != int:
        raise TypeError("zone value not int")

    # Check return_type is an option and correct for to_ref
    if return_type.lower() not in ("xyz", "llh", "enu"):
        raise ValueError("return_type must be either \"xyz\", \"llh\", \"enu\"")
    
    if (return_type.lower() == "enu") and (to_ref.lower() not in ("mga94", "mga2020")):
        raise ValueError("Can only give an enu result if to_ref is an MGA reference frame.")
    
    if (return_type.lower() != "enu") and (to_ref.lower() in ("mga94", "mga2020")):
        raise ValueError("Can only return mga coordinates in enu. return_type must be \"enu\"")

    # Finished input validation

    # Convert to enu to xyz and run xyz function
    lat, lon, psf, converge = con.grid2geo(zone, east, north)
    x, y, z = con.llh2xyz(lat, lon, height)

    # Change from_ref and to_ref to GDA equivalent if needed
    to_ref = mga_parse(to_ref)
    from_ref = mga_parse(from_ref)

    result = universal_transform(x, y, z, from_ref, to_ref, from_epoch, to_epoch, plate_motion, vcv, ignore_errors=ignore_errors)

    x = result["coords"]["x"]
    y = result["coords"]["y"]
    z = result["coords"]["z"]
    vcv = result["vcv"]

    # Assign output to dictionary

    output = {
        "type": return_type,
        "vcv": vcv,
    }

    if return_type.lower() == "xyz":
        output["coords"] = {
            "x": round(x, 4),
            "y": round(y, 4),
            "z": round(z, 4)
        }
    elif return_type.lower() == "llh":
        lat, lon, el_height = con.xyz2llh(x, y, z)
        output["coords"] = {
            "lat": round(lat,8),
            "lon": round(lon,8),
            "el_height": round(el_height,4)
        }
    elif return_type.lower() == "enu":
        lat, lon, el_height = con.xyz2llh(x, y, z)        
        hem, zone, east, north, psf, converge = con.geo2grid(lat, lon)
        output["coords"] = {
            "east": east,
            "north": north,
            "height": el_height,
            "zone": zone
        }
    
    return output
