import argparse
import os
import astropy.io.fits as pyfits


def divide_by_flat(specid, quartzid, extractdir):

    # Utility functions
    def order_file(iorder, specid, suffix1="o"):
        """Return the rectified image of the isolated order
        """
        return os.path.join(
            extractdir,
            "{}{}-order{}.fits".format(specid, suffix1, iorder)
        )

    orders = range(51, 77)
    for iorder in orders:
        print "Processing ", iorder
        hdulist = pyfits.open(order_file(iorder, specid))
        image = hdulist["SCI"].data
        flat_image = pyfits.open(order_file(iorder, quartzid))["SCI"].data
        flat_image = flat_image[13:49,:].mean(axis=0)
        image /= flat_image[None,:]
        hdulist.writeto(order_file(iorder, specid, suffix1="f"), clobber=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Divide by the quartz lamp image to correct for flat-field
        variations"""
    )
    parser.add_argument(
        "specid", type=str,
        help="""Prefix of object FITS spectrum file (e.g., p84)"""
    )
    parser.add_argument(
        "--quartzid", type=str, default="q69",
        help="""Prefix of quartz lamp FITS spectrum file"""
    )
    parser.add_argument(
        "--extractdir", type=str, default="Extract",
        help="""Directory containing the extracted orders"""
    )

    cmd_args = parser.parse_args()

    divide_by_flat(**vars(cmd_args))
