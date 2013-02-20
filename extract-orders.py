import numpy as np
import astropy.io.fits as pyfits
import pyregion
import argparse
import scipy.ndimage as ndi


def extract_single_order(imhdu, wavhdu, orderfilter, ordername, wavbox, wavmask):
    """
    Isolate a single order
    """
    # All pixels that we think are in this order
    ordermask = orderfilter.mask(wavhdu.data.shape)

    # First find wavelengths that ought to fall in the order
    orderwavs = wavhdu.data[ordermask & wavmask]
    if len(orderwavs):
        print "{} : {:.2f}-{:.2f}".format(ordername, orderwavs.min(), orderwavs.max())
    else:
        print "{} : No valid wavelengths found".format(ordername)

    # Second, look at wavelengths in the contiguous wavelength box that we found
    orderwavs = wavhdu.data[wavbox]
    orderwavs = orderwavs[orderwavs > 1000.0]
    if len(orderwavs):
        print "{}*: {:.2f}-{:.2f}".format(ordername, orderwavs.min(), orderwavs.max())
    else:
        print "{}*: No valid wavelengths found".format(ordername)

    # enclosing rectangle around this order
    bbox, = ndi.find_objects(ordermask.astype(int))
    print bbox, wavbox

    # 
    return imhdu.data[bbox]




def extract_orders(specfile, wavfile, regionfile,
                   wavmin=1000.0, wavmax=10000.0):
    """
    Go through all the orders, extracting each one
    """

    imhdu, = pyfits.open(specfile + ".fits")
    wavhdu, = pyfits.open(wavfile + ".fits")
    regions = pyregion.open(regionfile + ".reg")
    filters = regions.get_filter()
    ordernames = [box.attr[1]["text"] for box in regions]

    #
    # Auto-identify contiguous regions  in the wavelength map
    #

    # All pixels that have a valid wavelength
    wavmask = (wavhdu.data >= wavmin) & (wavhdu.data <= wavmax)
    labels, nlabels = ndi.label(wavmask, structure=np.ones((3,3)))
    wavboxes = ndi.find_objects(labels)
    wavboxes.reverse()
    
    for orderfilter, ordername, wavbox in zip(filters, ordernames, wavboxes):
        extract_single_order(imhdu, wavhdu, orderfilter, ordername, wavbox, wavmask)





if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Extract individual spectral orders from 
                       a Keck HIRES image"""
        )
    parser.add_argument("specfile", type=str, 
                        help="""Name of spectral image FITS file (sans extension)"""
                        )
    parser.add_argument("wavfile", type=str, 
                        help="""Name of wavelength FITS file (sans extension)"""
                        )
    parser.add_argument("regionfile", type=str, 
                        help="""Name of DS9 region file containing orders (sans extension)"""
                        )
    
    cmd_args = parser.parse_args()

    extract_orders(**vars(cmd_args))
