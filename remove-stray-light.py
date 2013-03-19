import argparse
import os
import astropy.io.fits as pyfits
import numpy as np


def remove_stray_light(specid, extractdir):

    # Utility functions
    def order_file(iorder, suffix1="b-cr", suffix2=""):
        """Return the rectified image of the  isolated order

        Optional argument suffix2 can be used to select the
        CR-rejected version, for instance
        """
        return os.path.join(
            extractdir,
            "{}{}-order{}{}.fits".format(specid, suffix1, iorder, suffix2)
        )

    orders = range(51, 77)
    for iorder in orders:
        print "Processing ", iorder
        hdulist = pyfits.open(order_file(iorder))
        image = hdulist["SCI"].data
        ny, nx = image.shape
        x = np.arange(nx)
        # Ad hoc form for the stray light:
        # ... Constant term
        stray0 = 4.0 + 2.0*(iorder - 51)/(76-51)
        # ... plus bump on left
        bump0 = 20.0
        xx = (15.0+x)/30.0
        stray = bump0*xx*np.exp(-xx) + stray0
        image -= stray[None, :]
        hdulist.writeto(order_file(iorder, suffix1="s"), clobber=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Remove stray light in extracted orders from
                       a Keck HIRES image"""
    )
    parser.add_argument(
        "specid", type=str,
        help="""Prefix of original FITS spectrum file (e.g., p84)"""
    )
    parser.add_argument(
        "--extractdir", type=str, default="Extract",
        help="""Directory containing the extracted orders"""
    )

    cmd_args = parser.parse_args()

    remove_stray_light(**vars(cmd_args))
