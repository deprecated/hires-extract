import os
import argparse
import astropy.io.fits as pyfits
import numpy as np
from pv_utils import make_grids


def remove_continuum(stampname, vrange, stampdir):
    # Read from the FITS file
    stamp_prefix = os.path.join(stampdir, stampname + "-stamp")
    hdulist = pyfits.open(stamp_prefix + ".fits")
    # Assume that the first HDU in the file is the image
    image = hdulist[0].data
    U, Y = make_grids(hdulist[0].header)
    assert U.shape == image.shape, \
        "Incompatible shapes {} versus {}".format(U.shape, image.shape)
    # Mask that includes only continuum pixels
    mask = (U < vrange[0]) | (U > vrange[1])
    continuum = np.zeros_like(image)
    # Fit a 2nd order polynomial independently to each row of the image
    # and subtract it from the original image
    for u, row, m, cont_row in zip(U, image, mask, continuum):
        coeffs = np.polyfit(u[m], row[m], deg=2)
        cont_row += np.poly1d(coeffs)(u)
        row -= cont_row
    # Save the fitted continuum as an extra HDU in the file
    hdulist.append(pyfits.ImageHDU(continuum, name="CONT"))
    hdulist.writeto(stamp_prefix + "-nc.fits", clobber=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Fit and remove continuum from single line spectrum"""
    )
    parser.add_argument(
        "stampname", type=str,
        help="""Prefix of stamp file (e.g., p84-N_I_5200)"""
    )
    parser.add_argument(
        "--vrange", type=float, nargs=2, default=[-20.0, 60.0],
        help="""Range of velocities to exclude from continuum fit"""
    )
    parser.add_argument(
        "--stampdir", type=str, default="Stamps",
        help="""Directory for placing the results"""
    )

    cmd_args = parser.parse_args()

    remove_continuum(**vars(cmd_args))