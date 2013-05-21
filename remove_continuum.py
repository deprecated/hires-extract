import os
import argparse
import astropy.io.fits as pyfits
import numpy as np
from pv_utils import make_grids


def remove_continuum(stampname, vrange, stampdir):
    # Read from the FITS file
    stamp_prefix = os.path.join(stampdir, stampname + "-stamp")
    hdulist = pyfits.open(stamp_prefix + ".fits")
    # 21 May 2013 - change extracted stamps to have named HDUs for the images
    image = hdulist["SCI"].data
    U, Y = make_grids(hdulist["SCI"].header)
    assert U.shape == image.shape, \
        "Incompatible shapes {} versus {}".format(U.shape, image.shape)
    # Mask that includes only continuum pixels
    mask = (U < vrange[0]) | (U > vrange[1])
    continuum = np.zeros_like(image)
    # This is the std of the entire row due to the uncertainty in the continuum fit
    sigma_r = np.zeros_like(image)
    # Fit a 2nd order polynomial independently to each row of the image
    # and subtract it from the original image
    for u, row, m, cont_row, sig_row in zip(U, image, mask, continuum, sigma_r):
        coeffs = np.polyfit(u[m], row[m], deg=2)
        cont_row += np.poly1d(coeffs)(u)
        row -= cont_row
        # calculate standard error of the mean of the residuals from the fit to this row
        sem = np.std(row[m] - cont_row[m])/np.sqrt(m.sum())
        # increase per-pixel uncertainty accordingly
        sig_row += sem
        
    # Save the fitted continuum as an extra HDU in the file
    hdulist.append(pyfits.ImageHDU(continuum, name="CONT"))
    # Save the per-row sigma as another extra HDU in the file.  This
    # is stored separately from sigma because the individual pixels in
    # the row are not independent, so it must be treated differently
    # in subsequent steps
    hdulist.append(pyfits.ImageHDU(sigma_r, name="SIGR"))
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
