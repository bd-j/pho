import numpy as np
import glob, os

from pho.apertures import *
from pho.photo import *

import astropy.io.fits as pyfits


# -- UVOT AB magnitude zeropoints from Hoversten (priv. comm.) --
uvot_zp_ab = {'w2': 17.30 + 1.755,
              'm2': 16.86 + 1.690,
              'w1': 17.46 + 1.557,
              'uu': 18.35 + 1.031,
              }


def load_uvot_images(galaxy, band, imdir=''):
    """ugh:
    :returns hdr:

    :returns images:
        List of images  [cps, exposure, sensitivity]

    :returns masks:
        Dictionary of mask images, keyed by strings:
        ['bkg', 'reg', 'reg1', 'reg5', 'reg20']
    """
    imroot = os.path.join(imdir, galaxy, galaxy)
    info = {'root': imroot, 'band': band}
    
    names = ['{root}_t{band}.fits', '{root}_t{band}_ex.fits', '{root}_t{band}_lss.fits']
    #print(names[0].format(**info))
    try:
        hdr = pyfits.getheader(names[0].format(**info), 1)
    except:
        hdr = pyfits.getheader(names[0].format(**info), 0)
    images = [pyfits.getdata(n.format(**info)) for n in names]

    masknames = ['bkg', 'reg', 'reg1', 'reg5', 'reg20']
    masks = [pyfits.getdata('{root}_{mask}.fits'.format(root=info['root'], mask=n))
            for n in masknames]

    return hdr, images, dict(zip(masknames, masks))


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
                      coi_correct=False, coi_mask=None, **kwargs):
    """Measure fluxes in UVOT imaging from Hoversten.

    :param coi_mask: (string, default: 'none')
        A switch to turn on masking of the regions affected by coincidence
        loss. One of
        * 'none' -- no masking
        * 'reg1' -- mask all regions affected by coincidence loss at the 1%
                    level or more
        * 'reg5' -- mask all regions affected by coincidence loss at the 5%
                    level or more
        * 'reg20' -- mask all regions affected by coincidence loss at the 20%
                     level or more
        * 'brown' -- correct for co-i loss using the Brown et al. 2012 formula.
    """
    uvotbands = ['w1', 'w2', 'm2']
    fluxes, uncertainties = [], []
    for band in uvotbands:
        hdr, images,  masks = load_uvot_images(galaxy, band, **kwargs)
        cps, exp, resp = images
        # Get a counts image for uncertainties
        counts = cps * (exp * resp)
        unc = np.sqrt(counts) / (exp * resp)
        if coi_mask =='brown':
            # Brown CoI correction.  Assuimes 1" pixels
            corcps = brown_coi_correction(cps)
            unc = corcps / cps * unc
            counts = corcps * (exp * resp)
            cps = corcps
        # Get a background value (in c/s) using the bkg mask
        bkg, bkg_unc = np.nanmedian(cps * masks['bkg']), 0.
        # Choose a mask
        m = masks.get(coi_mask, 1.0)
        print(np.min(m), np.sum(m))#, m.shape)
        
        # Photometer the BG subtracted, masked image
        flux, ps, units = photometer((cps - bkg) * m, hdr, [reg] + foreground_reg, mef=False)
        # Photometer the variance image (with bg_unc if it exists)
        var, _, _ = photometer(m * unc**2 + bkg_unc**2, hdr, [reg] + foreground_reg, mef=False)
        # All units are counts/sec, convert to Jy
        zp = 10**(-0.4 * uvot_zp_ab[band])
        unc = 3631. * zp * np.sqrt(var)
        fluxes.append(3631. * zp * flux)
        uncertainties.append(unc)

    return fluxes, uncertainties, uvotbands


if __name__ == "__main__":
    # Run parameters
    apname = 'dale'
    imdir = '/Users/bjohnson/Codes/SED/pho/data/'
    coi_mask = 'brown' # 'none' | 'reg1' | 'reg5' | 'reg20' | 'brown'

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
        # Skip galaxies without aperture info
        if name not in cat:
            fcat[i]['flag'] = -99
            continue
        # Propogate aperture data to the output catalog
        for k in cat[name]:
            try:
                fcat[i][k] = cat[name][k]
            except:
                pass
        # Make a ds9 region object from the aperture data for this galaxy
        reg = make_ds9_region(cat[name])
        # Measure the uvot fluxes in the aperture
        flux, unc, uband = measure_uvot_flux(name, reg, imdir=imdir, coi_mask=coi_mask)
        # Put the flux in the catalog
        if (len(flux) == 0):
            fcat[i]['flag'] = 1
        for f, u, b in zip(flux, unc, uband):
            bb = 'uvot_{}'.format(b)
            fcat[i][bb] = f
            fcat[i][bb + '_unc'] = u

    # Write out the catalog
    h = pyfits.hdu.table.BinTableHDU(fcat)
    h.header['FLUXUNIT'] = 'Jy'
    h.header['PAUNIT'] = 'Degrees E of N'
    h.header['APUNIT'] = 'Arcseconds'
    h.header['COI_MASK'] = coi_mask
    pyfits.writeto('uvot.{}_apertures.flux.fits'.format(apname), fcat, h.header, clobber=True)
