import numpy as np

from apertures import *
from photo import *


def find_images(imnames, galaxyname):
    possible = []
    for p in imnames:
        if galaxyname.upper() in p.upper():
            possible.append(p)
    return possible
    

def measure_pacs_flux(imagenames, gal_info):
    """Measure fluxes in background subtracted Herschel PACS imaging.
    """
    pacsbands = ['pacs70', 'pacs100', 'pacs160']
    reg = make_ds9_region(gal_info)
    fluxes, uncertainties, bands = [], [], []
    for imname in imagenames:
        image = pyfits.getdata(imname)
        hdr = pyfits.getheader(imname)
        im, unc = image[0:2,:,:]
        flux, ps, units =  photometer(im, hdr, [reg], mef=True)
        var, _, _ = photometer(unc**2, hdr, [reg], mef=True)
        # get flux in Jy, assuming image data in MJy/sr
        # flux = (flux * 1e6) * (ps.prod() * 2.35e-11)
        unc = np.sqrt(var) # * 1e6 * (ps.prod() * 2.35e-11)
        fluxes.append(flux)
        uncertainties.append(unc)
        bands.append([b for b in pacsbands if b in imname][0])
    return fluxes, uncertainties, bands


def measure_spire_flux(imagenames, gal_info):
    """Measure fluxes in background subtracted Herschel SPIRE imaging, given as
    a list of filenames.
    """
    spirebands = ['spire250', 'spire350', 'spire500']
    # make a ds9 region object from the aperture data for this galaxy
    reg = make_ds9_region(gal_info)
    # set up output
    fluxes, uncertainties, bands = [], [], []
    # loop over images
    for imname in imagenames:
        
        # Read the image data and header
        im = pyfits.getdata(imname)
        hdr = pyfits.getheader(imname)
        # Read the uncertainty image.  We assume this is matched pixel by pixel
        # with the flux image.
        uncname = imname.replace('scan.fits', 'scan.unc.fits')
        unc = pyfits.getdata(uncname)
        # Get image flux in the ds9 region
        flux, ps, units =  photometer(im, hdr, [reg])
        # Get image variance in the ds9 region
        var, _, _ = photometer(unc**2, hdr, [reg])
        # Get flux and unc in Jy, assuming image data in MJy/sr
        flux = (flux * 1e6) * (ps.prod() * 2.35e-11)
        unc = (np.sqrt(var) * 1e6) * (ps.prod() * 2.35e-11)
        fluxes.append(flux)
        uncertainties.append(unc)
        # set the band flag for this image
        bands.append([b for b in spirebands if b in imname][0])

    return fluxes, uncertainties, bands


if __name__ == "__main__":
    # Read the brown aperture data
    cat = read_brown_coordinates()
    cat = read_brown_apertures(cat=cat)

    # Get a bunch of image names
    pacs = glob.glob('../imaging/kingfish_pacs_scanam_v17/*fits')
    spire = glob.glob('../imaging/KINGFISH_SPIRE_v3.0_updated/*scan.fits')

    # Build output catalog
    fields = ['ra','dec','height', 'width', 'PA', 'flag',
              'pacs70', 'pacs70_unc', 'pacs100', 'pacs100_unc',
              'pacs160', 'pacs160_unc', 'spire250', 'spire250_unc',
              'spire350', 'spire350_unc', 'spire500', 'spire500_unc']
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
        # Find (background subtracted) images that might be of this galaxy
        thispacs = find_images(pacs, name)
        thisspire = find_images(spire, name)
        # Measure the fluxes of these images in the Brown apertures
        pflux, punc, pband = measure_pacs_flux(thispacs, cat[name])
        sflux, sunc, sband = measure_spire_flux(thisspire, cat[name])
        if (len(pflux) == 0) and (len(sflux) ==0):
            fcat[i]['flag'] = 1
        for f, u, b in zip(pflux+sflux, punc+sunc, pband+sband):
            fcat[i][b] = f
            fcat[i][b + '_unc'] = u


    h = pyfits.hdu.table.BinTableHDU(fcat)
    h.header['FlUXUNIT'] = 'Jy'
    h.header['PAUNIT'] = 'Degrees E of N'
    h.header['APUNIT'] = 'Arcseconds'
    pyfits.writeto('kingfish.brownapertures.flux.fits', fcat, h.header, clobber=True)


