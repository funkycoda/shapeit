import os
import sys
import argparse
import csv
import logging
from functools import partial
import json

from shapely.geometry import mapping, shape, Point, Polygon, MultiPoint
from shapely.wkt import dumps, loads
from shapely.ops import transform as stran

import fiona
from fiona.crs import from_epsg, from_string

from pyproj import Proj, transform


#setup logging
logging.basicConfig(filename='shape_it.log', format='%(asctime)s %(message)s', level=logging.DEBUG)

flog = logging.getLogger('Fiona')
flog.setLevel(logging.ERROR)


# setup the variables
code_3857 = '+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs'
proj_3857 = from_string(code_3857)
proj_4326 = from_epsg(4326)
p_in = Proj(proj_4326)
p_out = Proj(proj_3857)
# from_epsg(3857)

project = partial(transform, p_out, p_in)
output = "output"

IDX_PROGRESS = 0


class Within(object):
    def __init__(self, o):
        self.o = o
    def __lt__(self, other):
        return self.o['polygon'].within(other.o['polygon'])


def print_exit_message(msg):
    print("\nUnable to process the file. \n Messages:\n")
    print("=========================")
    # print("   Shapefile generator   ")
    # print("=========================")
    print(msg)
    print("=========================")
    print("\n")
    sys.exit(1)


def check_output_directory():
    if os.path.exists(output):
        if os.path.isfile(output):
            print_exit_message("Output folder already exists but is a file. Please delete it before proceeding.")
        else:
            print("Output folder already exists. Please rename or delete it before proceeding.")
            return True
    else:
        os.mkdir(output)
        return True

    return False


def process_file_with_species_location(csvfile, spcol, lngcol, latcol):
    global IDX_PROGRESS

    logging.info("Processing file: %s" % csvfile)
    print("Processing file: %s" % csvfile)

    snlist = {}

    with open(csvfile) as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        for row in reader:
            sciname = row[spcol]
            coords = [float(row[lngcol]), float(row[latcol])]

            if not snlist.has_key(sciname):
                snlist[sciname] = []

            snlist[sciname].append(coords)

    IDX_PROGRESS = 0
    for key in snlist:

        IDX_PROGRESS += 1
        update_progress(IDX_PROGRESS)

        coords = snlist[key]
        logging.debug("---------------------")

        # Sort the coordinates so it's easier to match them
        scoords = sorted(coords, key=lambda crd: crd[1])
        scoords = sorted(scoords, key=lambda crd: crd[0])

        logging.info('Processing species: %s', key)
        if len(coords) > 1:
            output_geoms = get_multi_polygons(scoords)
        else:
            output_geoms = get_single_polygon(scoords)

        write_shapefile(key.replace(' ', '_'), output_geoms)

        logging.debug("---------------------\n")

    # clean up the progress bar
    sys.stdout.write("\n")
    print("Completed")


def process_file_with_latlon(csvfile, lngcol, latcol):
    global IDX_PROGRESS

    logging.info("Processing file: %s", csvfile)
    print("Processing file: %s" % csvfile)

    lyrname = os.path.basename(csvfile)
    outfilename = os.path.splitext(lyrname)[0]
    lyrname = outfilename

    with open(csvfile) as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        rows = list(reader)

        coords = []
        points = []
        for i, row in enumerate(rows):
            IDX_PROGRESS += 1
            update_progress(IDX_PROGRESS)
            coords.append([float(row[lngcol]), float(row[latcol])])
            points.append(Point([float(row[lngcol]), float(row[latcol])]))

        # Sort the coordinates so it's easier to match them
        scoords = sorted(coords, key=lambda crd: crd[1])
        scoords = sorted(scoords, key=lambda crd: crd[0])

        output_geoms = get_multi_polygons(scoords)

        write_shapefile(lyrname, output_geoms)

    # clean up the progress bar
    sys.stdout.write("\n")
    print("Completed.\nCheck shape_it.log for additional details")


def get_single_polygon(coords):
    global IDX_PROGRESS

    logging.info("Processing %d points", len(coords))

    polygons = [{
        'polygon': MultiPoint([coords[0]]).convex_hull,
        'coords': [coords[0]]
    }]
    return polygons


def get_multi_polygons(coords):
    global IDX_PROGRESS

    logging.info("Processing %d points", len(coords))

    tcoords = coords

    polygons = [{
        'polygon': None,
        'coords': []
    }]
    poly_coords = []
    pts_coords = []

    # 1. Take the first 2 set of coords and
    # create a polygon if close <80km
    # or create 2 polygons
    c1 = tcoords[0]
    c2 = tcoords[1]

    if not is_feature_far(Point(c1), Point(c2)):  # 0.720
        polygons[0]['polygon'] = MultiPoint([c1, c2]).convex_hull
        polygons[0]['coords'].append(c1)
        polygons[0]['coords'].append(c2)
    else:
        polygons[0]['polygon'] = MultiPoint([c1]).convex_hull
        polygons[0]['coords'].append(c1)
        polygons.append({'polygon': MultiPoint([c2]).convex_hull, 'coords': [c2]})

    # 2. Starting with the 3rd coord, check which polygon it falls in and add
    # that point
    for crd in tcoords[2:]:
        IDX_PROGRESS += 1
        update_progress(IDX_PROGRESS)
        coord_added = False
        for i, poly in enumerate(polygons):
            pp = Point(crd[0], crd[1])
            d = poly['polygon'].distance(pp) * 111.111
            is_far = is_feature_far(poly['polygon'], pp)
            if not is_far:
                ecoords = poly['coords']
                ecoords.append(crd)
                poly['polygon'] = MultiPoint(ecoords).convex_hull
                coord_added = True
                break

        if not coord_added:
            polygons.append({'polygon': MultiPoint([crd]).convex_hull, 'coords': [crd]})

    # 3. Check and merge any polygons that are <80km
    polygons = check_and_merge_polygons(polygons)

    logging.info("Generated buffered point")
    return polygons

def check_and_merge_polygons(polygons):
    global IDX_PROGRESS

    logging.info('Check and sort...')
    # if it's a single polygon, we are good
    if len(polygons) == 1:
        return polygons

    spolys = sorted(polygons, key=Within)
    out_polygons = [spolys[0]]
    for poly in spolys[1:]:
        IDX_PROGRESS += 1
        update_progress(IDX_PROGRESS)
        poly_merged = False
        for idx, opoly in enumerate(out_polygons):
            if not is_feature_far(opoly['polygon'], poly['polygon']):
                opoly['coords'] = opoly['coords'] + poly['coords']
                opoly['polygon'] = MultiPoint(opoly['coords'] + poly['coords']).convex_hull
                poly_merged = True

        if not poly_merged:
            out_polygons.append(poly)

    return out_polygons


def write_shapefile(lyrname, features):
    logging.info("Writing layer: %s ", lyrname)

    output_path = "%s/%s/" % (output, lyrname)
    output_file = "%s/%s.shp" % (output_path, lyrname)
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    schema = {'geometry': 'Polygon', 'properties': {}}
    with fiona.open(output_file, 'w', crs=proj_4326, driver='ESRI Shapefile', schema=schema) as out:

        for feature in features:
            feature['polygon'] = feature['polygon'].buffer(0.20)
            out.write({
                'properties': {},
                'geometry': mapping(feature['polygon'])
            })

    logging.info("Generated polygon shapefile")


def is_feature_far(f1, f2, dist=80):
    return (f1.distance(f2) * 111.111) > dist


def update_progress(iprg):
    sys.stdout.write('\r[{0}]'.format('#'*(iprg/2)))
    sys.stdout.flush()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate shapefiles')

    # mandatory arguments
    parser.add_argument("input", metavar="input", help="Input CSV file")

    # optional arguments
    parser.add_argument("-g", "--group", dest="group", default=None,
                        help="Group column, e.g. scientific name. Default: None.")
    parser.add_argument("-lon", "--longitude", dest="lon", default='longitude',
                        help="Longitude column name, e.g. longitude. Default: longitude.")
    parser.add_argument("-lat", "--latitude", dest="lat", default='latitude',
                        help="Latitude column name, e.g. latitude. Default: latitude.")
    parser.add_argument("-o", "--output", dest="output", default=output,
                        help="Output directory. (default: ./output/)")

    args = parser.parse_args()
    output = args.output
    if check_output_directory():
        if args.group:
            process_file_with_species_location(args.input, args.group, args.lon, args.lat)
        else:
            process_file_with_latlon(args.input, args.lon, args.lat)
