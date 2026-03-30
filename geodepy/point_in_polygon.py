from __future__ import annotations

import re
import math
from pathlib import Path
from dataclasses import dataclass
from datetime import date

from shapely.geometry import Polygon, Point
from shapely.ops import unary_union
from shapely.prepared import prep

from pyproj import CRS, Transformer

import geodepy.convert
from geodepy.constants import Transformation, TransformationSD
import geodepy.transform

PLATE_ID_RE = re.compile(r"^[A-Z]{2}$")

def lon_to_180(lon):
    """
    Normalize longitude to between -180 and 180
    
    :param lon: Longitude
    :type lon: float
    :return: Longitude between -180 to 180
    :rtype: float
    """
    return ((lon + 180.0) % 360.0) - 180.0


def circular_mean_lon(lons_deg):
    """
    Circular mean longitude in degrees, returned in [-180,180).

    :param lons_deg: Longitdue in degrees
    :type lons_deg: list of floats
    :return: Returns the mean longitude
    :rtype: float
    """
    angles = [math.radians(l) for l in lons_deg]
    sx = sum(math.cos(a) for a in angles)
    sy = sum(math.sin(a) for a in angles)
    return lon_to_180(math.degrees(math.atan2(sy, sx)))


def unwrap_lons_to_reference(coords, ref_lon):
    """
    Shift each lon by +/-360 so it stays within 180° of ref_lon.
    This avoids seam jumps when projecting.

    :param coords: Array with lat and lon
    :type coords: list of tuples
    :param ref_lon: Mean Longitude
    :type ref_lon: float
    :return: Unwrapped Longitude
    :rtype: list of tuples
    """
    out = []
    for lon, lat in coords:
        lon = lon_to_180(lon)
        while lon - ref_lon > 180.0:
            lon -= 360.0
        while lon - ref_lon < -180.0:
            lon += 360.0
        out.append((lon, lat))
    return out

def read_morvel_plate_file(path):
    """
    Read MORVEL file and place into dictionary
    
    :param path: Path to morvel plate file
    :type path: str
    :return: Dictionary of plates with their boundaries
    :rtype: dict
    """
    plates= {}
    current_plate = None
    current_coords = []

    with open(path, "r", encoding="utf-8") as f:
        
        # Iterate through each line
        for raw in f:
            line = raw.strip()
            if not line:
                continue

            # If line is end of plate coords
            if "end of line segment" in line:

                # If coords have been collected, append to plates
                if current_plate and current_coords:
                    plates.setdefault(current_plate, []).append(current_coords)
                
                # Reset variables for new plate
                current_plate = None
                current_coords = []
                continue

            # Skip comments
            if line.startswith("*"):
                continue

            # If line is a plate ID
            if PLATE_ID_RE.fullmatch(line):

                # Catch in case previous plate isnt ended properly
                if current_plate and current_coords:
                    plates.setdefault(current_plate, []).append(current_coords)
                
                current_plate = line
                current_coords = []
                continue

            # If a current plate exists, map line to lon, lat and append
            if current_plate:
                lon, lat = map(float, line.replace(",", " ").split()[:2])
                lon = lon_to_180(lon)
                current_coords.append((lon, lat))
    
    # Catch at end if file not ended properly
    if current_plate and current_coords:
        plates.setdefault(current_plate, []).append(current_coords)

    return plates

def build_plate_index(dig_path):
    """
    Builds a projected polygon for each plate in its own local LAEA projection.
    This avoids the pole + seam problems of planar lon/lat polygons.

    Returns:
      index[pid] = {
        "prep": prepared geometry,
        "transformer": Transformer (EPSG:4326 -> local),
        "bounds": (minx, miny, maxx, maxy) in projected meters
      }
    
    :param dig_path: Path to morvel plate file
    :type dig_path: str
    :return: Dictionary of plates
    :rtype: dict
    """
    plates = read_morvel_plate_file(dig_path)
    index = {}

    for pid, segments in plates.items():
        all_lons = [lon for seg in segments for lon, lat in seg]
        all_lats = [lat for seg in segments for lon, lat in seg]
        if len(all_lons) < 4:
            continue

        # Calculate centre lat lon so projection can have least distortion
        lon0 = circular_mean_lon(all_lons)
        lat0 = sum(all_lats) / len(all_lats)

        # Create CRS with centre latl on of plate
        local_crs = CRS.from_proj4(
            f"+proj=laea +lat_0={lat0} +lon_0={lon0} +datum=WGS84 +units=m +no_defs"
        )

        # Create Transformer object from lat lon to Easting Northing
        transformer = Transformer.from_crs("EPSG:4326", local_crs, always_xy=True)

        polys = []
        for seg in segments:
            
            if len(seg) < 4:
                continue

            # Close ring if needed for big plates
            if seg[0] != seg[-1]:
                seg = seg + [seg[0]]

            # Unwrap to reference BEFORE projecting
            seg_adj = unwrap_lons_to_reference(seg, lon0)

            # Project coordinates, then build polygon in easting Northing
            xs, ys = transformer.transform(
                [lon for lon, lat in seg_adj],
                [lat for lon, lat in seg_adj]
            )

            ring_xy = list(zip(xs, ys))

            poly = Polygon(ring_xy)

            # Repair in planar coordinates
            if not poly.is_valid:
                poly = poly.buffer(0)

            if not poly.is_empty:
                polys.append(poly)

        if not polys:
            continue

        geom = unary_union(polys)
        if geom.is_empty:
            continue

        index[pid] = {
            "prep": prep(geom),
            "transformer": transformer,
            "bounds": geom.bounds,
        }

    return index

def plate_from_ll(lat, lon, plate_index):
    """
    Find plate for geographic lat/lon (degrees).

    :param lat: Latitude of point
    :type lat: float
    :param lon: Longitude of point
    :type lon: float
    :param plate_index: Dictionary of plates
    :type plate_index: dict
    :return: Plate ID
    :rtype: Array of strings or None
    """
    lon = lon_to_180(lon)

    hits = []
    for pid, rec in plate_index.items():
        # Project point into this plate's local CRS
        x, y = rec["transformer"].transform(lon, lat)

        # Quick bounds reject
        minx, miny, maxx, maxy = rec["bounds"]
        if x < minx or x > maxx or y < miny or y > maxy:
            continue
        
        # Test if point is in plate
        if rec["prep"].covers(Point(x, y)):
            hits.append(pid)

    if not hits:
        return None

    # deterministic (boundary cases)
    return sorted(hits)[0]

def plate_from_xyz(x, y, z, plate_index):
    """
    Converts XYZ to LLH and runs plate_from_ll to find plate

    :param x: Cartesian X (m)
    :type x: float
    :param y: Cartesian Y (m)
    :type y: float
    :param z: Cartesian Z (m)
    :type z: float
    :param plate_index: Dictionary of plates
    :type plate_index: dict
    :return: plate ID
    :rtype: Array of strings or None
    """
    lat, lon, h = geodepy.convert.xyz2llh(x, y, z)
    return plate_from_ll(lat, lon, plate_index)

@dataclass(frozen=True)
class PlatePole:
    """
    Class for storing plate rotation information
    """
    code: str
    lat_deg: float
    lon_deg: float
    rot_deg_ma: float
    rms_mm_yr: float = None

def load_poles(path):
    """
    Read poles file with plate motion information.

    Skips blank lines and comment lines beginning with '*'.

    Expected columns:
      ID lat lon rot [name...]
    where rot is deg/Ma anti-clockwise.

    :param path: Path to plate motion file
    :type path: str
    :return: Array of PlatePole objects containing each plate
    :rtype: list of PlatePole
    """
    path = Path(path)
    poles = []

    # Iterate through each line
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        
        # Skips commented lines
        if not line or line.startswith("*"):
            continue

        parts = line.split()
        if len(parts) < 4:
            continue

        # Store each part of line to different variables
        code = parts[0].upper()
        lat = float(parts[1])
        lon = float(parts[2])
        rot = float(parts[3])

        rms = None

        # If rms is present assign value
        if len(parts) >=5:
            try:
                rms = float(parts[4])
            except ValueError:
                rms = None
                pass

        # Add PlatePole object to poles array
        poles.append(PlatePole(code=code, lat_deg=lat, lon_deg=lon, rot_deg_ma=rot, rms_mm_yr=rms))

    return poles

def euler_to_drx_dry_drz(lat_deg, lon_deg, rot_deg_ma):
    """
    Convert Euler pole (lat, lon, rot deg/Ma) to GeodePy rotation rates d_rx,d_ry,d_rz in arcsec/yr.

    :param lat_deg: Euler pole Latitude in degrees
    :type lat_deg: float
    :param lon_deg: Euler pole Longitude in degrees
    :type lon_deg: float
    :param rot_deg_ma: Euler pole rotation in deg/Ma
    :type rot_deg_ma: float
    :return: Rotation paramters in rx, ry, rz
    :rtype: tuple of floats (d_rx, d_ry, d_rz)
    """
    phi = math.radians(lat_deg)
    lam = math.radians(lon_deg)
    omega_arcsec_yr = rot_deg_ma * 3600.0 / 1_000_000.0

    d_rx = round(omega_arcsec_yr * math.cos(phi) * math.cos(lam), 8)
    d_ry = round(omega_arcsec_yr * math.cos(phi) * math.sin(lam), 8)
    d_rz = round(omega_arcsec_yr * math.sin(phi), 8)
    return d_rx, d_ry, d_rz


def sd_drot_from_rms_and_pole(rms_mm_yr, lat_deg, lon_deg, radius_m=6378137.0):
    """
    Compute sd_d_rx, sd_d_ry, sd_d_rz (arcsec/yr) from:
      - RMS speed (1σ) in mm/yr (Table RMS)
      - Euler pole latitude/longitude (deg)

    Assumes uncertainty is magnitude-only along pole direction.
    
    :param rms_mm_yr: RMS of pole in mm/yr
    :type rms_mm_yr: float
    :param lat_deg: Euler pole Latitude in degrees
    :type lat_deg: float
    :param lon_deg: Euler pole Longitude in degrees
    :type lon_deg: float
    :return: RMS values for rotation paramters
    :rtype: tuple of floats (sd_d_rx, sd_d_ry, sd_d_rz)
    """
    # RMS speed -> m/yr
    sigma_v = rms_mm_yr / 1000.0

    # rad/yr
    sigma_Omega_rad_yr = sigma_v / radius_m

    # arcsec/yr
    sigma_Omega_arcsec_yr = sigma_Omega_rad_yr * (180.0 / math.pi) * 3600.0

    # unit vector to pole
    phi = math.radians(lat_deg)
    lam = math.radians(lon_deg)
    ux = math.cos(phi) * math.cos(lam)
    uy = math.cos(phi) * math.sin(lam)
    uz = math.sin(phi)

    # Component SDs
    sd_d_rx = abs(ux) * sigma_Omega_arcsec_yr
    sd_d_ry = abs(uy) * sigma_Omega_arcsec_yr
    sd_d_rz = abs(uz) * sigma_Omega_arcsec_yr

    return sd_d_rx, sd_d_ry, sd_d_rz

def plate_transformation(
    plate_code,
    poles_file,
    ref_epoch = date(2020, 1, 1),
    from_datum = "",
    to_datum = "plate motion",
):
    """
    Build a GeodePy Transformation object for a plate code.

    The Transformation is defined with:
      tx/ty/tz/sc/rx/ry/rz = 0
      d_rx/d_ry/d_rz = Euler pole converted to arcsec/yr

    :param plate_code: Plate code eg "AU"
    :type plate_code: str
    :param poles_file: Path to plate motion file
    :type poles_file: str
    :param ref_epoch: Reference epoch of transformation
    :type ref_epoch: datetime.date
    :param from_datum: Source datum for transformation
    :type from_datum: str
    :param to_datum: Destination datum for transformation
    :type to_datum: str
    :return: Geodepy.constants.Transformation Object
    """
    poles = load_poles(poles_file)

    wanted = plate_code.strip().upper()
    pole_map = {p.code: p for p in poles}

    if wanted not in pole_map:
        raise KeyError(f"Plate '{plate_code}' not found in {Path(poles_file).name}")

    pole = pole_map[wanted]
    d_rx, d_ry, d_rz = euler_to_drx_dry_drz(pole.lat_deg, pole.lon_deg, pole.rot_deg_ma)

    tf_sd = None

    # Create TransformationSD object
    if pole.rms_mm_yr is not None:
        sd_d_rx, sd_d_ry, sd_d_rz = sd_drot_from_rms_and_pole(pole.rms_mm_yr, pole.lat_deg, pole.lon_deg)
        tf_sd = TransformationSD(
            sd_tx=0.0, sd_ty=0.0, sd_tz=0.0,
            sd_sc=0.0,
            sd_rx=0.0, sd_ry=0.0, sd_rz=0.0,
            sd_d_tx=0.0, sd_d_ty=0.0, sd_d_tz=0.0,
            sd_d_sc=0.0,
            sd_d_rx=sd_d_rx, sd_d_ry=sd_d_ry, sd_d_rz=sd_d_rz
        )
    # Return Transformation object
    return Transformation(
        from_datum=f"{from_datum}{pole.code}",
        to_datum=to_datum,
        ref_epoch=ref_epoch,
        tx=0.0, ty=0.0, tz=0.0,
        sc=0.0,
        rx=0.0, ry=0.0, rz=0.0,
        d_tx=0.0, d_ty=0.0, d_tz=0.0,
        d_sc=0.0,
        d_rx=d_rx, d_ry=d_ry, d_rz=d_rz,
        tf_sd=tf_sd
    )

def universal_plate_motion_transformation(x ,y, z, 
                                          from_epoch, to_epoch, vcv=None,
                                          plate_file = "other_files/MORVEL56_plates.dig", 
                                          poles_file = "other_files/NNR-MORVEL56_poles.txt", 
                                          ref_epoch = date(2020, 1, 1)):
    """
    Given ECEF XYZ coordinates, find the plate, find the plate motion and transform point.

    :param x: Cartesian X (m)
    :type x: float
    :param y: Cartesian Y (m)
    :type y: float
    :param z: Cartesian Z (m)
    :type z: float
    :param from_epoch: Source epoch for plate motion
    :type from_epoch: datetime.date
    :param to_epoch: Destination epoch for plate motion
    :type to_epoch: datetime.date
    :param vcv: Optional 3*3 numpy array in Cartesian units to propagate tf uncertainty
    :type vcv: numpy.ndarray
    :param plate_file: Path to morvel plate file
    :type plate_file: str
    :param poles_file: Path to plate motion file
    :type poles_file: str
    :param ref_epoch: Reference epoch of transformation
    :type ref_epoch: datetime.date
    :return: Transformed XYZ and VCV
    :rtype: tuple (xtrans, ytrans, ztrans, vcv)
    """

    # Create Index of plates
    plates = build_plate_index(plate_file)
    
    # Find which plate point is on
    plate_id = plate_from_xyz(x, y, z, plates)
    if not plate_id:
        raise ValueError("No plate found for given coordinates")
    
    print(f"Plate motion on {plate_id} plate")

    # Create GeodePy Transformation object for plate
    transformation = plate_transformation(plate_id, poles_file, ref_epoch=ref_epoch)

    # Complete transformation between given dates
    xtrans, ytrans, ztrans, vcv = geodepy.transform.plate_motion_transformation(x, y, z, from_epoch, to_epoch, transformation, vcv)

    return xtrans, ytrans, ztrans, vcv
