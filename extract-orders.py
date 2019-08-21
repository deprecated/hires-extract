import numpy as np
import astropy.io.fits as pyfits
import pyregion
import argparse
import scipy.ndimage as ndi
import scipy.interpolate
import skimage.morphology
import os

INTERPOLATE = "linear"

# Mapping between the boxes that we added by hand and
# the labels of the contiguous patches
box2label = {
    51: 24, 52: 23, 53: 22, 54: 21, 55: 20, 56: 19, 57: 18,
    58: 17, 59: 16, 60: 15, 61: 14, 62: 14, 63: 14, 64: 13,
    65: 12, 66: 11, 67: 10, 68: 9,  69: 8,  70: 7,  71: 6,
    72: 5,  73: 4,  74: 3,  75: 2,  76: 1
}


# Orders up to 71 are well-separated and work with both methods
ordermax = 20 # 71


def shrink_box(r, h=8.0):
    """
    Change the height of a pyregion Box region
    """
    r.coord_list[3] = h
    return r


def dilate_mask(mask, margin=1):
    """
    Add an extra margin (default 1 pixel) around all true areas
    of a logical mask
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
                   onlyorders=None, 
                   wavmin=1000.0, wavmax=10000.0):
    """
    Go through all the orders, extracting each one
    """

    imhdu = pyfits.open(specfile + ".fits")[0]
    badpixhdu = pyfits.open(specfile + ".fits")["BADPIX"]
    wavhdu, = pyfits.open(wavfile + ".fits")
    regions = pyregion.open(regionfile + ".reg")

    # Tilt angles from the horizontal
    tilts = [r.coord_list[4] for r in regions]

    # Fractional y shift of box center from nearest grid point
    ypix_fractions = [r.coord_list[1] - int(r.coord_list[1]) for r in regions]

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
    labels, nlabels = ndi.label(wavmask, structure=np.ones((3, 3)))

    # save a copy of the labels for debugging
    pyfits.PrimaryHDU(labels).writeto("orders-labels.fits", clobber=True)

    print("Number of order boxes found:", len(ordernames))
    print("Number of objects found:", nlabels)

    for widefilter, orderfilter, ordername, tilt, ypix_frac in zip(
            wide_filters, filters, ordernames, tilts, ypix_fractions):
        # All pixels that we think are in the central part
        # of the slit in this order
        ordermask = orderfilter.mask(wavhdu.data.shape)
        # and the same for the entire order (adding some padding)
        widemask = dilate_mask(widefilter.mask(wavhdu.data.shape), 3)

        iorder = int(ordername.split()[-1])
        if not (onlyorders is None or iorder in onlyorders):
            # Option for skipping all but some orders
            continue

        # First find wavelengths that ought to fall in the order
        orderwavs = wavhdu.data[ordermask & wavmask]
        if len(orderwavs):
            print("{}: {:.2f}-{:.2f}".format(
                ordername, orderwavs.min(), orderwavs.max()))
        else:
            print("{}: No valid wavelengths found".format(ordername))

        # Second, look at wavelengths in the contiguous wavelength box
        # that we found
        label = box2label[iorder]
        orderwavs = wavhdu.data[labels == label]
        if len(orderwavs):
            print("Label {}: {:.2f}-{:.2f}".format(
                label, orderwavs.min(), orderwavs.max()))
        else:
            print("{}*: No valid wavelengths found".format(ordername))

        print()

        # enclosing rectangle around this entire order
        bbox, = ndi.find_objects(widemask.astype(int))
        # Add a few more pixels at the top to give equal top/bottom margins
        # We have to do it like this since slice.stop is read-only
        # and tuples are immutable
        start, stop, step = bbox[0].indices(wavhdu.data.shape[0])
        bbox = (slice(start, stop+6), bbox[1])

        imorder = imhdu.data.copy()[bbox]
        badpixorder = badpixhdu.data.copy()[bbox]
        wavorder = wavhdu.data.copy()[bbox]
        # Construct a mask of all pixels both methods say
        # should be in this order
        if (iorder <= ordermax):
            # These orders are the easiest to deal with
            m = widemask & (labels == label)
            # This mask is just for the central strip,
            # which is what we need for the wavs
            cm = ordermask & (labels == label)
        else:
            m = labels == label
            cm = labels == label
        m = m[bbox]
        cm = cm[bbox]

        cm = cm & (wavorder > 4000.0)
        print(
            "Number of good wavelength pixels found in order box:",
            np.sum(m), np.sum(cm))
        mm = widemask[bbox]     # less stringent mask

        # Use a single average wavelength for each column
        ny, nx = wavorder.shape
        # print(wavorder[::10, ::200].astype(int))
        # print(cm[::10, ::200])
        # meanwav = np.average(wavorder, axis=0, weights=cm.astype(int))
        meanwav = np.nansum(wavorder*cm.astype(int), axis=0) / np.nansum(cm.astype(int), axis=0)
        # print(meanwav[::200])
        meanwav = np.vstack([meanwav]*ny)

        #
        # Remove the horizontal tilt of the orders
        #
        # First the linear tilt
        yshifts = np.arange(nx)*np.tan(np.radians(tilt))
        # Then the parabolic residual distortion
        yshifts += parabolic_distorsion(iorder)*(
            2*np.arange(nx).astype(float)/nx - 1.0)**2
        # Finally, align the box center to the pixel grid
        yshifts += ypix_frac
        jshifts = yshifts.astype(int)  # required shift of each column
        jshiftset = set(jshifts)
        # Amount to trim off the top of the strip at the end
        jtrim = jshifts.max()
        if INTERPOLATE is not None:
            # These are the grid points we want
            grid_x, grid_y = np.meshgrid(
                np.arange(nx, dtype=np.float), np.arange(ny, dtype=np.float)
            )
            # And these are the coordinates we currently have
            # Note that only the y's change, not the x's
            x, y = grid_x, grid_y - yshifts[np.newaxis, :]
            # Interpolate image onto new grid
            imorder = scipy.interpolate.griddata(
                (x.ravel(), y.ravel()), imorder.ravel(),
                (grid_x, grid_y), method=INTERPOLATE
            )
            # Use nearest-neighbor for the masks,
            # so they don't get converted to reals
            badpixorder = scipy.interpolate.griddata(
                (x.ravel(), y.ravel()), badpixorder.ravel(),
                (grid_x, grid_y), method="nearest"
            )
            m = scipy.interpolate.griddata(
                (x.ravel(), y.ravel()), m.ravel(),
                (grid_x, grid_y), method="nearest"
            )
            mm = scipy.interpolate.griddata(
                (x.ravel(), y.ravel()), mm.ravel(),
                (grid_x, grid_y), method="nearest"
            )
        else:
            jshifts = np.vstack([jshifts]*ny)  # Expand back to 2D
            for jshift in jshiftset:  # Consider each unique value of jshift
                # Split up into one or more contiguous chunks
                # that have this value of jshift
                chunklabels, nlabels = ndi.label(jshifts == jshift)
                for chunk in ndi.find_objects(chunklabels):
                    # apply the shift to all the arrays
                    # (except meanwav, which is constant in y)
                    imorder[chunk] = np.roll(imorder[chunk], -jshift, axis=0)
                    badpixorder[chunk] = np.roll(badpixorder[chunk], -jshift, axis=0)
                    m[chunk] = np.roll(m[chunk], -jshift, axis=0)
                    mm[chunk] = np.roll(mm[chunk], -jshift, axis=0)

        # Trim the useless space off the top
        imorder = imorder[:-jtrim, :]
        badpixorder = badpixorder[:-jtrim, :]
        meanwav = meanwav[:-jtrim, :]
        # And save each order to FITS files
        pri = pyfits.PrimaryHDU()
        sci = pyfits.ImageHDU(imorder, name='SCI')
        bad = pyfits.ImageHDU(badpixorder, name='BAD')
        wav = pyfits.ImageHDU(meanwav, name='WAV')
        outfile = "{}-order{}.fits".format(os.path.split(specfile)[-1], iorder)
        outfile = os.path.join(outdir, outfile)
        pyfits.HDUList([pri, sci, wav, bad]).writeto(outfile, clobber=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Extract individual spectral orders from
                       a Keck HIRES image"""
    )
    parser.add_argument(
        "specfile", type=str,
        help="""Name of spectral image FITS file (sans extension)"""
    )
    parser.add_argument(
        "wavfile", type=str,
        help="""Name of wavelength FITS file (sans extension)"""
    )
    parser.add_argument(
        "regionfile", type=str,
        help="""Name of DS9 region file containing orders (sans extension)"""
    )
    parser.add_argument(
        "--outdir", "-o", type=str, default="Extract",
        help="""Directory for placing the results"""
    )
    parser.add_argument(
        "--onlyorders", type=int, metavar="N", nargs="*", default=None,
        help="""If set, only process these particular orders.  Default is to
        process all the orders."""
    )

    cmd_args = parser.parse_args()

    extract_orders(**vars(cmd_args))
