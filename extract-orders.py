import numpy as np
import astropy.io.fits as pyfits
import pyregion
import argparse
import scipy.ndimage as ndi
import skimage.morphology
import os

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

def dilate_mask(mask, margin=1):
    """
    Add an extra margin (default 1 pixel) around all true areas of a logical mask
    """
    return skimage.morphology.dilation(
               mask.astype(np.uint8), 
               skimage.morphology.square(2*margin + 1)
               ).astype(bool)


def parabolic_distorsion(iorder):
    """
    This is most notable in the red, and seems to be zero for orders > 60
    """
    if iorder > 60:
       return 0.0
    else:
       return 2.0*(60.0 - iorder)/(60.0 - 51.0)
    


def extract_orders(specfile, wavfile, regionfile, outdir, 
                   wavmin=1000.0, wavmax=10000.0):
    """
    Go through all the orders, extracting each one
    """

    imhdu, = pyfits.open(specfile + ".fits")
    wavhdu, = pyfits.open(wavfile + ".fits")
    regions = pyregion.open(regionfile + ".reg")

    # Tilt angles from the horizontal
    tilts = [r.coord_list[4] for r in regions]

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

    for widefilter, orderfilter, ordername, tilt in zip(wide_filters, filters, ordernames, tilts):
        # All pixels that we think are in the central part of the slit in this order
        ordermask = orderfilter.mask(wavhdu.data.shape)
        # and the same for the entire order, but adding some padding
        widemask = dilate_mask(widefilter.mask(wavhdu.data.shape), 3)
        

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
        # Add a few more pixels at the top to give equal top/bottom margins
        # We have to do it like this since slice.stop is read-only and tuples are immutable
        start, stop, step = bbox[0].indices(wavhdu.data.shape[0])
        bbox = (slice(start, stop+6), bbox[1])

        imorder = imhdu.data.copy()[bbox]
        wavorder = wavhdu.data.copy()[bbox]
        # Construct a mask of all pixels both methods say should be in this order
        if (iorder <= ordermax):
            # These orders are the easiest to deal with
            m = widemask & (labels == label)
        else:
            m = labels == label
        m = m[bbox]
        print "Number of good wavelength pixels found in order box: ", np.sum(m)
        mm = widemask[bbox] # less stringent mask

        #
        # Remove the horizontal tilt of the orders
        #
        ny, nx = wavorder.shape
        # First the linear tilt
        yshifts = np.arange(nx)*np.tan(np.radians(tilt))
        # Then the parabolic residual distortion 
        yshifts += parabolic_distorsion(iorder)*(2*np.arange(nx).astype(float)/nx - 1.0)**2 
        jshifts = yshifts.astype(int) # required shift of each column
        jshiftset = set(jshifts)
        jtrim = jshifts.max() # Amount to trim off the top of the strip at the end
        jshifts = np.vstack( [jshifts]*ny ) # Expand back to 2D
        for jshift in jshiftset:            # Consider each unique value of jshift
            # Split up into one or more contiguous chunks that have this value of jshift
            chunklabels, nlabels = ndi.label(jshifts == jshift) 
            for chunk in ndi.find_objects(chunklabels):
                # apply the shift to all the arrays
                imorder[chunk] = np.roll(imorder[chunk], -jshift, axis=0)
                wavorder[chunk] = np.roll(wavorder[chunk], -jshift, axis=0)
                m[chunk] = np.roll(m[chunk], -jshift, axis=0)
                mm[chunk] = np.roll(mm[chunk], -jshift, axis=0)
            
        # Use a single average wavelength for each column
        meanwav = np.sum(wavorder*m, axis=0) / m.sum(axis=0)
        meanwav = np.vstack( [meanwav]*ny )
        # Trim the useless space off the top
        imorder = imorder[:-jtrim,:]
        meanwav = meanwav[:-jtrim,:]
        # And save each order to FITS files
        pri = pyfits.PrimaryHDU()
        sci = pyfits.ImageHDU(imorder, name='SCI')
        wav = pyfits.ImageHDU(meanwav, name='WAV')
        outfile = "{}-order{}.fits".format(os.path.split(specfile)[-1], iorder)
        outfile = os.path.join(outdir, outfile)
        pyfits.HDUList([pri, sci, wav]).writeto(outfile, clobber=True)


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
    parser.add_argument("--outdir", "-o", type=str, default="Extract", 
                        help="""Directory for placing the results"""
                        )

    cmd_args = parser.parse_args()

    extract_orders(**vars(cmd_args))
