import numpy as np
import glob, os

from pho.apertures import *
from pho.photo import *

import astropy.io.fits as pyfits

    
def load_uvot_images(galaxy, band, imdir=''):
    """ugh:
    :returns hdr:

    :returns images:
        List of images  [cps, exposure, sensitivity]

    :returns masks:
        Mask images.  [bkg, reg, reg20]
    """
    imroot = os.path.join(imdir, galaxy, galaxy)
    info = {'root': imroot, 'band': band}
    
    names = ['{root}_t{band}.fits', '{root}_t{band}_ex.fits', '{root}_t{band}_lss.fits']
    hdr = pyfits.getheader(names[0].format(**info), 1)
    images = [pyfits.getdata(n.format(**info)) for n in names]

    names = ['{root}_bkg.fits', '{root}_reg.fits', '{root}_reg20.fits']
    masks = [pyfits.getdata(n.format(**info)) for n in names]

    return hdr, images, masks


def brown_coi_correction(cpspa):
    """
    :param cpspa:
        counts/sec/square arcsec.  Note typical UVOT pixels are 1 square arcsec.

    :returns ccor:
        Counts corrected for coincidence loss according to Brown et al.
    """
    ft, df = 0.011088, 0.0155844
    cprime = 80. * cpspa * ft
    ccor = -np.log(1. - cprime - 2 * cprime**2) / (80. * ft * (1-df))
    return ccor


def measure_uvot_flux(galaxy, reg, foreground_reg=[],
                      coi_correct=False, **kwargs):
    """Measure fluxes in UVOT imaging from Hoversten.
    """
    uvotbands = ['w1', 'w2', 'm2']
    fluxes, uncertainties = [], []
    for band in uvotbands:
        hdr, images,  masks = load_uvot_images(galaxy, band, **kwargs)
        cps, exp, resp = images
        # Get a counts image for uncertainties
        counts = cps * (exp * resp)
        unc = np.sqrt(counts) / (exp * resp)
        if coi_correct:
            # Brown CoI correction.  Assuimes 1" pixels
            corcps = brown_coi_correction(cps)
            unc = corcps / cps * unc
            counts = corcps * (exp * resp)
            cps = corcps

        # Get a background value (in c/s) using the bkg mask
        bkg, bkg_unc = np.median(cps * masks[0]), 0.
        # Photometer the BG subtracted image
        flux, ps, units = photometer(cps - bkg, hdr, [reg] + foreground_reg, mef=False)
        # Photometer the variance image (with bg_unc if it exists)
        var, _, _ = photometer(unc**2 + bkg_unc**2, hdr, [reg] + foreground_reg, mef=False)
        # All units are counts/sec
        unc = np.sqrt(var)
        fluxes.append(flux)
        uncertainties.append(unc)

    return fluxes, uncertainties, uvotbands


if __name__ == "__main__":
    # Run parameters
    apname = 'dale'
    imdir = '/Users/bjohnson/Codes/SED/pho/data/'

    # Get a bunch of galaxy names that you want to photometer
    galaxy_names = glob.glob(imdir + '*/')
    galaxy_names = np.sort([g.split(os.path.sep)[-2] for g in galaxy_names])
    galaxy_names = ['ngc7793']

    # Read the aperture data
    if apname == 'dale':
        cat = read_sings_apertures()
        dimensions = ['a', 'b']
    elif apname == 'brown':
        cat = read_brown_coordinates(coordfile=imdir+'brown_coordinates.txt')
        cat = read_brown_apertures(cat=cat, aperturefile=imdir+'brown_apertures.txt')
        dimensions = ['width', 'height']

    # Build output catalog
    fields = (['ra','dec'] + dimensions +
              ['PA', 'flag', 'uvot_w1', 'uvot_w1_unc',
               'uvot_w2', 'uvot_w2_unc', 'uvot_m2', 'uvot_m2_unc',
              ])
    dt = [('name', 'S20')] + [(f, np.float) for f in fields]
    fcat = np.zeros(len(galaxy_names), dtype=np.dtype(dt))

    # Loop over galaxies filling catalog with photometry
    for i, name in enumerate(galaxy_names):
        fcat[i]['name'] = name
        # skip galaxies without aperture info
        if name not in cat:
            fcat[i]['flag'] = -99
            continue
        # Propogate aperture data to the output catalog
        for k in cat[name]:
            try:
                fcat[i][k] = cat[name][k]
            except:
                pass
        # make a ds9 region object from the aperture data for this galaxy
        reg = make_ds9_region(cat[name])
        # Measure the uvot fluxes in the aperture
        flux, unc, uband = measure_uvot_flux(name, reg, imdir=imdir)
        if (len(flux) == 0):
            fcat[i]['flag'] = 1
        for f, u, b in zip(flux, unc, uband):
            fcat[i][b] = f
            fcat[i][b + '_unc'] = u


    h = pyfits.hdu.table.BinTableHDU(fcat)
    h.header['FlUXUNIT'] = 'counts/sec'
    h.header['PAUNIT'] = 'Degrees E of N'
    h.header['APUNIT'] = 'Arcseconds'
    pyfits.writeto('uvot.{}_apertures.flux.fits'.format(apname), fcat, h.header, clobber=True)
