import numpy as np

from apertures import *
from photo import *


def find_images(imnames, galaxyname):
    possible = []
    for p in imnames:
        if galaxyname.upper() in p.upper():
            possible.append(p)
    return possible
    
def load_uvot_images(imroot):
    """ugh
    """
    names = []

    return [hdr] + [pyfits.getdata(n) for n in names]

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

def measure_uvot_flux(imagenames, reg, foreground_reg=[]):
    """Measure fluxes in background subtracted Herschel PACS imaging.
    """
    uvotbands = ['w1', 'w2', 'm2']
    fluxes, uncertainties, bands = [], [], []
    for imname in imagenames:
        hdr, cps, bkg, exp, resp = load_uvot_images(imname)
        # Get a counts image for uncertainties
        counts = cps * (exp * resp)
        unc = np.sqrt(counts) / (exp * resp)
        if coi_correct:
            # Brown CoI correction
            corcps = brown_coi_correction(cps)
            unc = corcps / cps * unc
            counts = corcps * (exp * resp)

        # Get a background value (in c/s) using the bkg mask
        bg, bg_unc = np.median(cps * bkg), 0.
        # Photometer the BG subtracted image
        flux, ps, units = photometer(cps - bg, hdr, [reg] + foreground_reg, mef=True)
        # Photometer the variance image (with bg_unc if it exists)
        var, _, _ = photometer(unc**2 + bg_unc**2, hdr, [reg] + foreground_reg, mef=True)
        # All units are counts/sec
        unc = np.sqrt(var)
        fluxes.append(flux)
        uncertainties.append(unc)
        bands.append([b for b in pacsbands if b in imname][0])
    return fluxes, uncertainties, bands


if __name__ == "__main__":
    # Read the aperture data
    #cat = read_sings_regions()
    cat = read_brown_coordinates()
    cat = read_brown_apertures(cat=cat)
    apname = 'brown'

    # Get a bunch of image names
    #uvot = glob.glob()

    # Build output catalog
    fields = ['ra','dec','a', 'b', 'PA', 'flag',
              'uvot_w1', 'uvot_w1_unc', 'uvot_w2', 'uvot_w2_unc',
              'uvot_m2', 'uvot_m2_unc',
              ]
    dt = [('name', 'S20')] + [(f, np.float) for f in fields]
        
    galaxy_names = cat.keys()
    # galaxy_names = ['ngc3190']
    fcat = np.zeros(len(galaxy_names), dtype=np.dtype(dt))
    for i, name in enumerate(galaxy_names):
        fcat[i]['name'] = name
        for k in cat[name]:
            try:
                fcat[i][k] = cat[name][k]
            except:
                pass
        # make a ds9 region object from the aperture data for this galaxy
        reg = make_ds9_region(cat[name])
        # Find (background subtracted) images that might be of this galaxy
        thispacs = find_images(uvot, name)
        # Measure the fluxes of these images in the Brown apertures
        flux, unc, uband = measure_uvot_flux(thisuvot, reg)
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
