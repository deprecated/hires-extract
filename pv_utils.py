"""
General utility routines for dealing with position-velocity (PV) images
"""

import numpy as np


def make_grids(hdr):
    """Create 2D velocity and position grids from WCS info in FITS header

    """
    nx, ny = [hdr["NAXIS{}".format(i)] for i in 1, 2]
    i0, u0, dudi, dudj = [hdr[k] for k in
                          "CRPIX1", "CRVAL1", "CD1_1", "CD1_2"]
    j0, y0, dydi, dydj = [hdr[k] for k in
                          "CRPIX2", "CRVAL2", "CD2_1", "CD2_2"]
    ii, jj = np.meshgrid(np.arange(nx) + 1, np.arange(ny) + 1)
    u = u0 + dudi*(ii - i0) + dudj*(jj - j0)
    y = y0 + dydi*(ii - i0) + dydj*(jj - j0)
    return u, y


if __name__ == "__main__":
    print "TODO: write tests and demo code for this module"
