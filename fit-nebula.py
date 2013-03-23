"""
Fit and remove the nebular components to the line profiles

"""
import argparse
import os
import json
import numpy as np
import astropy.io.fits as pyfits
# import fit_utils
from doublet_utils import multiplet_moments, equivalent_doublet
from pv_utils import make_grids


def find_doublet(moments):
    """
    Convert moments into two delta functions
    """
    components = {"I1": [], "I2": [], "v1": [], "v2": []}
    for intensity, mean, sigma, skewness in zip(
            *[moments[m] for m in "intensity", "mean", "sigma", "skewness"]):
        doublet = equivalent_doublet(intensity, mean, sigma, skewness)
        for k in ["I1", "I2", "v1", "v2"]:
            components[k].append(doublet[k])
    return components


def _find_moments_array(S, U, minfrac=0.05):
    # Maximum brightness in each row
    Smax = S.max(axis=1)[:, None]
    # Mask that deselects pixels that are too faint
    m = S > minfrac*Smax
    M0 = np.sum(m*S, axis=1)
    ubar = np.sum(m*S*U, axis=1)/M0
    # Center the velocities by subtracting off the mean
    Uc = U - ubar[:, None]
    # Variance from centered 2nd moment
    sigma = np.sqrt(np.sum(m*S*Uc**2, axis=1)/M0)
    # Skewness from centered 3rd moment
    skew = np.sum(m*S*Uc**3, axis=1)/(M0*sigma**3)
    return {
        k: v.astype(float).tolist() for k, v in {
            "intensity": M0,
            "mean": ubar,
            "sigma": sigma,
            "skewness": skew,
        }.items()
    }


def _find_moments_multiplet(image, velarray):
    # Re-purpose the routine from doublet_utils to do the dirty
    # work.  It may be necessary to find a more efficient way in
    # the future, but hopefully not.
    moments = {"intensity": [], "mean": [], "sigma": [], "skewness": []}
    for intensities, vels in zip(image, velarray):
        row_moments = multiplet_moments(vels, intensities)
        for key in moments:
            # Cast to python float to facilitate dumping to JSON
            moments[key].append(float(row_moments[key]))
    return moments


def find_moments(image, velarray, method="array", **kwds):
    """Find the first four velocity moments for each row in a 2D
    position-velocity image.  Velocities are provided in a second
    array velarray of the same shape as image.  The velocity axis is
    assumed to be the last axis of the array.  Returns a dict of
    lists: {"intensity": [...], "mean": [...], "sigma": [...],
    "skewness": [...]}.

    """
    return {"multiplet": _find_moments_multiplet,
            "array": _find_moments_array}[method](image, velarray, **kwds)


def save_data(data, prefix, method="json"):
    if method == "json":
        with open(prefix + "-moments.json", "w") as f:
            json.dump(data, f, indent=2)
    elif method == "table":
        raise NotImplementedError
    else:
        raise ValueError("Unrecognised method: " + method)


def fit_nebula(stampname, vrange,
               stampdir="Stamps", extra_suffix="", min_fraction=0.05):
    # Read from the FITS file
    stamp_prefix = os.path.join(
        stampdir, stampname + "-stamp-nc" + extra_suffix)
    hdulist = pyfits.open(stamp_prefix + ".fits")
    # Assume that the first HDU in the file is the image
    image = hdulist[0].data
    U, Y = make_grids(hdulist[0].header)
    # When calculating the rectangular box to use for the moments (but
    # only then), we ignore variations in the velocity scale with
    # y-position
    u = U.mean(axis=0)          # Collapse velocities
    i1 = len(u[u < vrange[0]])
    i2 = len(u[u < vrange[1]])
    data = find_moments(image[:, i1:i2], U[:, i1:i2], minfrac=min_fraction)
    y = Y[:, i1:i2].mean(axis=1)
    data.update(position=y.astype(float).tolist())
    data.update(find_doublet(data))
    save_data(data, stamp_prefix)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Fit and remove nebular component"""
    )
    parser.add_argument(
        "stampname", type=str,
        help="""Prefix of stamp file (e.g., p84-N_I_5200)"""
    )
    parser.add_argument(
        "--vrange", type=float, nargs=2, default=[-10.0, 40.0],
        help="""Range of velocities to use for calculating moments"""
    )
    parser.add_argument(
        "--stampdir", type=str, default="Stamps",
        help="""Directory for placing the results"""
    )
    parser.add_argument(
        "--extra-suffix", type=str, default="",
        help="""Extra suffix for images (e.g., -dd)"""
    )
    parser.add_argument(
        "--min-fraction", type=float, default=0.05,
        help="""Minimum fraction of peak brightness in order that a
        pixel should contribute to the velocity moments"""
    )
    cmd_args = parser.parse_args()
    fit_nebula(**vars(cmd_args))
