import yaml
import os
import argparse
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt

DB_FILE_NAME = "line-database.json"
# Physical units
LIGHT_SPEED_CGS = 2.99792458e10
KM_CGS = 1.0e5
# Limits of velocity window
VMIN, VMAX = -100.0, 100.0
# Size of pixels in arcsec
PLATE_SCALE = 0.382
# j-limits of the extracted order
J1, J2 = 13, 49


def extract_stamps(specid, extractdir, stampdir):
    """
    For each line in the database,
    produce a postage stamp calibrated 2D spectrum
    """
    # Utility functions
    def order_file(iorder, suffix1="o", suffix2=""):
        """Return the CR-rejected, rectified, de-overlapped
        image of the  isolated order"""
        return os.path.join(
            extractdir,
            "{}{}-order{}{}.fits".format(specid, suffix1, iorder, suffix2)
        )

    def sanitize(s):
        """
        Make a string safe for use as a filename

        Eliminate [ and ], and change SPACE to _
        """
        return s.translate(None, "[]").replace(" ", "_")

    def linsolve(x, y, plotname=None):
        """
        Find the best linear fit: "y = m x + b " for two vectors, x, y.

        Returns coefficients m, b

        If the the optional argument plotname is set,
        then also save a plot showing the fit
        """
        A = np.vstack([x, np.ones(len(x))]).T
        m, b = np.linalg.lstsq(A, y)[0]
        if plotname:
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(x, y - (m*x + b), 'b.')
            ax.plot(x, np.zeros_like(x) + b, 'r-')
            ax.set_xlabel("pixel")
            ax.set_ylabel("velocity - vfit")
            ax.set_title("Linear fit: " + plotname)
            ax.set_ylim(-2.0, 2.0)
            figname = os.path.join(stampdir, "linsolve-" +
                                   specid + "-" + plotname + ".pdf")
            fig.savefig(figname)
        return m, b

    # Read in the line database
    linedb = yaml.load(open(os.path.join(stampdir, DB_FILE_NAME)))
    for lineid, emline in linedb.items():
        inhdulist = pyfits.open(order_file(emline["iorder"]))
        print lineid
        print inhdulist.info()
        image = inhdulist["SCI"].data
        wavs = inhdulist["WAV"].data
        wavs[~np.isfinite(wavs)] = 0.0  # eliminate NaNs
        # Velocities in km/s with respect to rest wavelength of line
        vels = (LIGHT_SPEED_CGS/KM_CGS) * (wavs - emline["wav"])/emline["wav"]
        print vels[wavs != 0.0].max(), vels[wavs != 0.0].min()
        # Extract the window of VMIN -> VMAX
        ny, nx = vels.shape
        v1d = vels[ny/2, :]
        print v1d.max()
        i1 = np.argmin(np.abs(v1d - VMIN))
        i2 = np.argmin(np.abs(v1d - VMAX))
        # print i1, i2
        # print v1d[i1] - VMIN, v1d[i1] - VMAX
        # print v1d[i2] - VMIN, v1d[i2] - VMAX
        vpix = v1d[i1:i2]
        imstamp = image[J1:J2, i1:i2]
        ny, nx = imstamp.shape
        # Determine the WCS parameters
        # Fit a linear function to 1D velocity versus x-index
        xpix = np.arange(nx, dtype=np.float) + 1.0
        try:
            m, b = linsolve(xpix, vpix, plotname=sanitize(lineid))
        except ValueError:
            print "**** Error....  couldn't solve linear  equation"
            print "**** Skipping " + lineid
            continue
        hdustamp = pyfits.PrimaryHDU(imstamp)
        hdustamp.header.update(
            CRPIX1=0.0, CRVAL1=b, CD1_1=m, CD1_2=0.0,
            CTYPE1="VELO", CUNIT1="km/s",
            CRPIX2=0.5*(1 + J2 - J1), CRVAL2=0.0, CD2_1=0.0, CD2_2=PLATE_SCALE,
            CTYPE2="PARAM", CUNIT2="arcsec",
        )
        outfile = os.path.join(
            stampdir, "-".join([specid, sanitize(lineid), "stamp.fits"])
        )
        hdustamp.writeto(outfile, clobber=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Extract stamps of individual spectral lines from
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
    parser.add_argument(
        "--stampdir", type=str, default="Stamps",
        help="""Directory for placing the results"""
    )

    cmd_args = parser.parse_args()

    extract_stamps(**vars(cmd_args))
