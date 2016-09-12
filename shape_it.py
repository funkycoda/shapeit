import os
import sys
import argparse
import csv
import logging
from functools import partial
import json

from shapely.geometry import mapping, shape, Point, Polygon
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


def process_file(csvfile, colno):
    logging.info("Processing file: %s" % csvfile)
    print("Processing file: %s" % csvfile)

    snlist = {}

    with open(csvfile) as f:
        reader = csv.DictReader(f, skipinitialspace=True)
        # print(reader.fieldnames)
        for row in reader:
            sciname = row['MasterTax_SciName']
            coords = [float(row['lat_dec']), float(row['lon_dec'])]

            if not snlist.has_key(sciname):
                snlist[sciname] = []

            snlist[sciname].append(coords)

    iprg = 0
    for key in snlist:

        iprg += 1
        update_progress(iprg)

        coords = snlist[key]
        logging.debug("---------------------")

        if len(coords) > 2:
            generate_shapefile_polygon(key, coords)
        else:
            generate_shapefile_multipoint(key, coords)

        logging.debug("---------------------\n")

    # clean up the progress bar
    sys.stdout.write("\n")
    print("Completed")


def generate_shapefile_polygon(sciname, coords):
    logging.info("Processing species: %s ", sciname)

    output_path = "%s/%s/" % (output, sciname)
    output_file = "%s/%s.shp" % (output_path, sciname)
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    tcoords = []
    for coord in coords:
        x, y = transform(p_in, p_out, float(coord[1]), float(coord[0]))
        tcoords.append([x, y])

    polygon = Polygon(tcoords)
    buffered = polygon.buffer(20000.0)
    polygon_out = stran(project, buffered)

    schema = { 'geometry': 'Polygon', 'properties': { 'name': 'str' } }
    with fiona.open(output_file, 'w', crs=proj_4326, driver='ESRI Shapefile', schema=schema) as out:

        out.write({
            'properties': {
                'name': sciname
            },
            'geometry': mapping(polygon_out)
        })

    logging.info("Generated polygon shapefile")


def generate_shapefile_multipoint(sciname, coords):
    logging.info("Processing species: %s ", sciname)

    output_path = "%s/%s/" % (output, sciname)
    output_file = "%s/%s.shp" % (output_path, sciname)
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    tcoords = []
    for coord in coords:
        x, y = transform(p_in, p_out, float(coord[1]), float(coord[0]))
        point = Point(x, y)
        buffered = point.buffer(20000.0)
        point_out = stran(project, buffered)
        tcoords.append(point_out)

    schema = { 'geometry': 'Polygon', 'properties': { 'name': 'str' } }
    with fiona.open(output_file, 'w', crs=proj_4326, driver='ESRI Shapefile', schema=schema) as out:

        for poly in tcoords:        
            out.write({
                'properties': {
                    'name': sciname
                },
                'geometry': mapping(poly)
            })

    logging.info("Generated point shapefile")


def generate_shapefile_point(sciname, lat, lng):
    logging.info("Processing species: %s ", sciname)

    output_path = "%s/%s/" % (output, sciname)
    output_file = "%s/%s.shp" % (output_path, sciname)
    if not os.path.exists(output_path):
        os.mkdir(output_path)

    x, y = transform(p_in, p_out, lat, lng)
    point = Point(x, y)
    buffered = point.buffer(20000.0)
    point_out = stran(project, buffered)

    schema = { 'geometry': 'Polygon', 'properties': { 'name': 'str' } }
    with fiona.open(output_file, 'w', crs=proj_4326, driver='ESRI Shapefile', schema=schema) as out:

        out.write({
            'properties': {
                'name': sciname
            },
            'geometry': mapping(point_out)
        })

    logging.info("Generated point shapefile")


def update_progress(iprg):
    sys.stdout.write('\r[{0}]'.format('#'*(iprg/2)))
    sys.stdout.flush()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate shapefiles')

    # mandatory arguments
    parser.add_argument("input", metavar="input", help="Input CSV file")

    # optional arguments
    parser.add_argument("-g", "--group", metavar="group", dest="group", type=int, default=2, help="Group column, e.g. scientific name. Default: 3.")
    parser.add_argument("-o", "--output", metavar="output", dest="output", default=output, help="Output directory. (default: ./output/)")

    args = parser.parse_args()
    output = args.output
    folder_checks_out = check_output_directory()
    if folder_checks_out:
        process_file(args.input, args.group)
