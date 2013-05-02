import yaml
import argparse
import os
import astropy.io.fits as pyfits
from scipy.special import erf
import numpy as np

DB_FILE_NAME = "quartz-database.json"


def remove_overlap(specid, extractdir, suffix):
    """
    Loop through the orders, dealing with each one
    """
    # Utility functions
    def order_file(iorder, suffix1="s", suffix2=""):
        """Return the rectified image of the  isolated order

        Optional argument suffix2 can be used to select the
        CR-rejected version, for instance
        """
        return os.path.join(
            extractdir,
            "{}{}-order{}{}.fits".format(specid, suffix1, iorder, suffix2)
        )

    def window(y, y1, y2, sigma):
        """Window function with Gaussian broadening"""
        return (0.5 * (erf((y - y1)/sigma) + 1.0) *
                0.5 * (erf((y2 - y)/sigma) + 1.0))

    def other_order(y1_key, w_key, sig_key, j1, j2):
        # Theoretical profile of overlapping order
        win = window(y, pars[y1_key], pars[y1_key] + pars[w_key],
                     pars[sig_key])
        # Only consider pixels where the order is not too weak
        m = win[j1:j2] > 0.01
        if m.any():
            # Find average strength of order in the
            # restricted range close to the edge
            av = np.mean(imcol[j1:j2][m]/win[j1:j2][m])
            # And apply that strength, weighted by the window function,
            # over the entire column
            return av*win
        else:
            return 0.0

    quartz_solutions = yaml.load(open(DB_FILE_NAME))
    orders = sorted(quartz_solutions.keys())
    for iorder in orders:
        print "Processing ", iorder
        pars = quartz_solutions[iorder]
        hdulist = pyfits.open(order_file(iorder, suffix1=suffix))
        image = hdulist["SCI"].data
        ny, nx = image.shape
        y = np.arange(ny) + 1.0  # 1-based scale for FITS y-pixels
        # Fiddle the positions
        if "y1_02" in pars:
            pars["y1_02"] += 0.25
        if "y1_03" in pars:
            pars["y1_03"] += 0.25
        # Loop over each image column
        for i in range(nx):
            imcol = image[:, i]
            # Subtract the linear term but only when it is the quartz lamp
            if specid.startswith("q69"):
                imcol -= pars["const"] + pars["linear"] * y
            # Subtract the previous order at the top
            if "y1_02" in pars:
                imcol -= other_order("y1_02", "w_02", "sigma_02", 52, 55)
            # Subtract the next order at the bottom
            if "y1_03" in pars:
                imcol -= other_order("y1_03", "w_03", "sigma_03", 8, 11)
        hdulist.writeto(order_file(iorder, suffix1="o"), clobber=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Correct the order overlap in extracted orders from
                       a Keck HIRES image"""
    )
    parser.add_argument(
        "specid", type=str,
        help="""Prefix of original FITS spectrum file (e.g., p84)"""
    )
    parser.add_argument(
        "--suffix", type=str, default="s",
        help="""Prefix of original FITS spectrum file (e.g., p84)"""
    )
    parser.add_argument(
        "--extractdir", type=str, default="Extract",
        help="""Directory containing the extracted orders"""
    )

    cmd_args = parser.parse_args()

    remove_overlap(**vars(cmd_args))
