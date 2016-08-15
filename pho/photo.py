import glob
import numpy as np
from sedpy import photometer, ds9region
import astropy.io.fits as pyfits
import astropy.wcs as pywcs


__all__ = ["photometer"]


def photometer(image, header, regions, pad=[0, 0], mef=False):
    """Given an image (result of pyfits.getdata), a FITS header (result of
    pyfits.getheader), and a set of regions, produce total fluxes within those
    regions.  Returns a flux, a platescale, and a flux unit if available in the
    header.

    Note that this is not particularly precise in the case where the region
    size is comparable to or only slightly larger than the pixel size.

    :param image:
        Result of pyfits.getdata

    :param header:
        Result of pyfits.getheader

    :param regions:
        List of ds9region.Region() objects.

    :param mef:
        Set to True if the header is from a multi-extension FITS file

    :returns fluxes:
        An ndarray of the same length as `regions`, giving the summed image
        values for pixels whose centers fall within the regions.

    :returns ps:
        An ndarray of shape (2,) that gives the size of each side of a pixel,
        in arcsec.

    :returns bunit:
        The BUNIT field of the header, if avaialable.  Otherwise `None`.
    """
    # Get a WCS object
    wcs = pywcs.WCS(header)
    # Get the pixel transformation matrix
    try:
        cd = wcs.wcs.cd
    except:
        # Sometimes full matrix isn't available...
        cd = np.diag(wcs.wcs.cdelt[0:2])
    # get the pixel dimensions, in arcsec
    ps = np.hypot(*cd*3600.)
    # get a huge list of x and y coordinates for each pixel
    yy, xx = np.indices(image.shape)
    # convert these to ra, dec and put in an ndarray of shape (npix, 2)
    if mef:
        ra, dec, _ = wcs.wcs_pix2world(xx.flatten(), yy.flatten(), np.array([0]), 0)
    else:
        ra, dec = wcs.wcs_pix2world(xx.flatten(), yy.flatten(), 0)
    ra = ra.flatten()
    dec = dec.flatten()
    points = np.vstack((ra, dec)).T
    flatim = image.flatten()
    flux = []
    # Loop over regions
    for region in regions:
        # Find which pixels are in the region...
        if region.shape == 'polygon':
            # Use RA, Dec, because figuring out the conventions is a nightmare
            sel = region.contains(points=points, wcs=wcs, fast=False, pad=pad)
        else:
            # Use x, y because equation for ellipse in shperical is a disaster.
            sel = region.contains(x=xx.flatten(), y=yy.flatten(), wcs=None)
        # sum the image in these pixels and attach that sum to the output list
        flux.append(np.nansum(flatim[sel]))

    return np.array(flux), ps, header.get('BUNIT', None)
