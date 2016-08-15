import numpy as np
from sedpy import ds9region


__all__ = ["read_brown_apertures", "read_brown_coordinates",
           "read_sings_apertures",
           # "process_aperture_line", "process_coordinate_line",
           "make_ds9_region"]


def make_ds9_region(info, **extras):
    """Given a dictionary of region info (ra, dec, PA and width, height or a,
    b), build and return a ds9 region. If 'width' is a key of the dictionary,
    returns a Polygon object describing this rectangle as a serties of
    vertices, otherwise returns an Ellipse object.
    """
    if 'width' in info:
        # get coordinates of vertices
        ra, dec = dumb_corners(**info)
        # make a list of ["ra","dec",....] for the vertices
        defstring = [str(v) for pair in zip(ra, dec) for v in pair]
        # Join the list of string coordinates
        defstring = ','.join(defstring)
        # Now instantiate the region object
        reg = ds9region.Polygon(defstring)
        return reg
    elif 'a' in info:
        # deg, deg, arcsec, arcsec, deg
        defstring = '{ra},{dec},{a},{b},{PA}'.format(**info)
        return ds9region.Ellipse(defstring)


def read_sings_apertures():
    # 235749.9 −323525
    cat = {}
    cat['ngc7793'] = {'ra': 359.4579166666667, 'dec': -32.590277,
                      'a': 716 / 2., 'b': 526 / 2., 'PA': 98, }
    return cat


def read_brown_apertures(cat={}, aperturefile='data/brown_apertures.txt'):
    """Read the Brown et al. aperture data from a text file and return a
    dictionary with the aperture data dictionaries (height, width, PA, and
    reference), keyed by galaxy name.
    """
    with open(aperturefile) as f:
        lines = f.readlines()

    hdr_line = len(lines)
    for i, line in enumerate(lines):
        if line.startswith('(arcsec)'):
            hdr_line = i
        if i <= hdr_line:
            # We haven't hit the end-of-hdr indicator, go to the next i
            continue
        n, h, w, p, s = process_aperture_line(line)
        aperture = {'height': h, 'width': w, 'PA': p, 'specref': s}
        if n in cat:
            cat[n].update(aperture)
        else:
            cat[n] = aperture
    return cat


def process_aperture_line(line):
    """Split a line of the txt file containing aperture data into the actual
    info.
    """
    cols = line.split('\t')
    name = cols[0].lower().replace(' ', '')
    height, width = [float(c.strip()) for c in cols[1].split('x')]
    pa = float(cols[2])
    specsource = cols[3]
    return name, height, width, pa, specsource


def read_brown_coordinates(cat={}, coordfile='data/brown_coordinates.txt'):
    """Read the Brown et al. aperture central coordinates form a txt file.
    Returns a dictionary of coordinate dictionaries, keyed by galaxy name.
    """
    with open(coordfile) as f:
        lines = f.readlines()
    hdr_line = len(lines)
    for i, line in enumerate(lines):
        if line.startswith('Class'):
            hdr_line = i
        if i <= hdr_line:
            # We haven't hit the end-of-hdr indicator, go to the next i
            continue
        if line.startswith('\n'):
            break
        n, ra, dec = process_coordinate_line(line)
        if n in cat:
            cat[n].update({'ra': ra, 'dec': dec})
        else:
            cat[n] = {'ra': ra, 'dec': dec}
    return cat


def process_coordinate_line(line):
    cols = line.split('\t')
    name = cols[0].lower().replace(' ', '')
    coords = [float(c) for c in cols[1].split()]
    ra = 15 * (coords[0] + coords[1] / 60. + coords[2] / 3600.)
    dec = (coords[3] + coords[4] / 60. + coords[5] / 3600.)
    if '-' in cols[1] and dec > 0:
        dec *= -1

    return name, ra, dec


def dumb_corners(ra=None, dec=None, PA=None,
                 height=None, width=None, **extras):
    """Given an ra,dec center, dimensions, and PA, return a list of rectangular
    vertices for use in constructing ds9 polygon apertures. Approximate using
    2-D plane.
    """
    size = np.array([width / 3600., height / 3600.])
    theta = np.deg2rad(PA)
    T = np.array([[np.cos(theta), np.sin(theta)],
                  [-np.sin(theta), np.cos(theta)]])
    corners = np.array([[-0.5, -0.5],
                        [-0.5, 0.5],
                        [0.5, 0.5],
                        [0.5, -0.5]])
    newcorners = np.dot(corners * size[None, :], T)
    outra = ra + newcorners[:, 0] / np.cos(np.deg2rad(dec))
    outdec = dec + newcorners[:, 1]
    return outra, outdec


def smart_corners(ra=None, dec=None, PA=None,
                  height=None, width=None):
    # Convert rotation axis and corners to cartesian coordinates

    # Build rotation matrix

    pass
