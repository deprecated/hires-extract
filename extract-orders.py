import numpy as np
import astropy.io.fits as pyfits
import pyregion
import argparse
import scipy.ndimage as ndi

# Mapping between the boxes that we added by hand and the labels of the contiguous patches
box2label = {
    "51" : 24, "52" : 23, "53" : 22, "54" : 21, "55" : 20, "56" : 19, "57" : 18, 
    "58" : 17, "59" : 16, "60" : 15, "61" : 14, "62" : 14, "63" : 14, "64" : 13, 
    "65" : 12, "66" : 11, "67" : 10, "68" : 9,  "69" : 8,  "70" : 7,  "71" : 6, 
    "72" : 5,  "73" : 4,  "74" : 3,  "75" : 2,  "76" : 1
    }


# Orders up to 72 are well-separated and work with both methods
ordermax = 72

def shrink_box(r, h=8.0):
    """
    Change the height of a pyregion Box region
    """
    r.coord_list[3] = h
    return r


def extract_orders(specfile, wavfile, regionfile,
                   wavmin=1000.0, wavmax=10000.0):
    """
    Go through all the orders, extracting each one
    """

    imhdu, = pyfits.open(specfile + ".fits")
    wavhdu, = pyfits.open(wavfile + ".fits")
    regions = pyregion.open(regionfile + ".reg")
    wide_filters = regions.get_filter()

    # Restrict attention to the center of each slit
    regions = pyregion.ShapeList([shrink_box(region) for region in regions])
    filters = regions.get_filter()
    
    ordernames = [box.attr[1]["text"] for box in regions]

    #
    # Auto-identify contiguous regions  in the wavelength map
    #

    # All pixels that have a valid wavelength
    wavmask = (wavhdu.data >= wavmin) & (wavhdu.data <= wavmax)
    labels, nlabels = ndi.label(wavmask, structure=np.ones((3,3)))

    # save a copy of the labels for debugging
    pyfits.PrimaryHDU(labels).writeto("orders-labels.fits", clobber=True)

    print "Number of order boxes found: ", len(ordernames)
    print "Number of objects found: ", nlabels

    for widefilter, orderfilter, ordername in zip(wide_filters, filters, ordernames):
        # All pixels that we think are in the central part of the slit in this order
        ordermask = orderfilter.mask(wavhdu.data.shape)
        # and the same for the entire order
        widemask = widefilter.mask(wavhdu.data.shape)
        

        # First find wavelengths that ought to fall in the order
        orderwavs = wavhdu.data[ordermask & wavmask]
        if len(orderwavs):
            print "{} : {:.2f}-{:.2f}".format(ordername, orderwavs.min(), orderwavs.max())
        else:
            print "{} : No valid wavelengths found".format(ordername)

        # Second, look at wavelengths in the contiguous wavelength box that we found
        label = box2label[ordername.split()[-1]]
        iorder = int(ordername.split()[-1])
        orderwavs = wavhdu.data[labels == label]
        if len(orderwavs):
            print "Label {}: {:.2f}-{:.2f}".format(label, orderwavs.min(), orderwavs.max())
        else:
            print "{}*: No valid wavelengths found".format(ordername)

        print

        # enclosing rectangle around this entire order
        bbox, = ndi.find_objects(widemask.astype(int))
        
        imorder = imhdu.data[bbox]
        wavorder = wavhdu.data[bbox]
        # Construct a mask of all pixels both methods say should be in this order
        if (iorder <= ordermax):
            # These orders are the easiest to deal with
            m = widemask & (labels == label)
        else:
            m = labels == label
        m = m[bbox]
        print "Number of wavelength pixels found: ", np.sum(m)
        mm = widemask[bbox] # less stringent mask
        ny, nx = wavorder.shape
        # Use a single average wavelength for each column
        meanwav = np.sum(wavorder*m, axis=0) / m.sum(axis=0)
        meanwav = np.vstack( [meanwav]*ny )
        pyfits.PrimaryHDU(imorder*mm).writeto("test-order{}-im.fits".format(iorder), 
                                           clobber=True)
        pyfits.PrimaryHDU(meanwav).writeto("test-order{}-wav.fits".format(iorder), 
                                           clobber=True)


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
