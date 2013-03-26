"""
Fit and remove the nebular components to the line profiles

"""
import argparse
import os
import json
import numpy as np
from numpy.polynomial import Chebyshev as T
import astropy.io.fits as pyfits
import lmfit
from doublet_utils import multiplet_moments, equivalent_doublet
from pv_utils import make_grids
from fit_utils import init_single_component, model_minus_data
from fit_utils import model, save_params, YDOMAIN


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


def main(stampname, vrange, ylo, yhi, stampdir="Stamps",
         extra_suffix="", min_fraction=0.05):

    # Step 1: Read data from the FITS file
    stamp_prefix = os.path.join(stampdir, stampname + "-stamp-nc")
    if extra_suffix:
        stamp_prefix = '-'.join([stamp_prefix, extra_suffix])
    # Open FITS file in "update" mode since we will be saving new
    # images to it later
    hdulist = pyfits.open(stamp_prefix + ".fits", mode="update")
    # Assume that the first HDU in the file is the image
    hdr = hdulist[0].header
    image = hdulist[0].data
    U, Y = make_grids(hdr)

    # Step 2: Calculate moments and 2-delta decomposition
    #
    # When calculating the rectangular box to use for the moments (but
    # only then), we ignore variations in the velocity scale with
    # y-position
    u = U.mean(axis=0)          # Collapse velocities
    i1 = len(u[u < vrange[0]])
    i2 = len(u[u < vrange[1]])
    mdata = find_moments(image[:, i1:i2], U[:, i1:i2], minfrac=min_fraction)
    y = Y[:, i1:i2].mean(axis=1)
    mdata.update(position=y.astype(float).tolist())
    mdata.update(find_doublet(mdata))
    save_data(mdata, stamp_prefix)

    # Step 3: Fit polynomials to the BG portion of 2-delta decomposition
    #
    # Mask in y that selects only BG positions
    mask = ((y > ylo[0]) & (y < ylo[1])) | ((y > yhi[0]) & (y < yhi[1]))
    print mask
    print y[mask]
    neb_coeffs = {"I1": None, "I2": None, "v1": None, "v2": None}
    for linepar in neb_coeffs.keys():
        print linepar
        print np.array(mdata[linepar])[mask]
        neb_coeffs[linepar] = T.fit(
            y[mask], np.array(mdata[linepar])[mask], 2,
            domain=YDOMAIN
        )

    # Step 4: initialize 2D model with 2 Gaussian components
    #
    params = lmfit.Parameters()
    # Assume all components have a constant width to start with
    init_single_component(params, "A", i_coeffs=neb_coeffs["I1"],
                          u_coeffs=neb_coeffs["v1"],
                          w_coeffs=[0.0, 0.0, 5.0])
    init_single_component(params, "B", i_coeffs=neb_coeffs["I2"],
                          u_coeffs=neb_coeffs["v2"],
                          w_coeffs=[0.0, 0.0, 5.0])

    # Step 5: fit the model to the BG portion of the 2d data
    #
    # Construct 2d mask of only BG positions and only velocities in
    # the defined window
    Ymask = ((Y > ylo[0]) & (Y < ylo[1])) | ((Y > yhi[0]) & (Y < yhi[1]))
    Umask = (U >= vrange[0]) & (U < vrange[1])
    bgmask = Ymask & Umask
    result = lmfit.minimize(model_minus_data, params,
                            args=(U[bgmask], Y[bgmask], image[bgmask]))
    print result.message

    # Print a verbose report of the error estimates and correlatons
    lmfit.report_errors(params)

    # Step 6: Save everything
    #
    imbg = model(U, Y, params)
    improp = image - imbg

    # Save the resulting BG-subtracted image in-place to the same FITS
    # file as a separate image HDU and do the same with the model
    # nebular BG image and the original guess too.  In each case, try
    # to save to an existing HDU with the right name first
    try:
        hdulist["PROP"].data = improp
    except KeyError:
        hdulist.append(pyfits.ImageHDU(improp, hdr, name="PROP"))

    try:
        hdulist["NEB"].data = imbg
    except KeyError:
        hdulist.append(pyfits.ImageHDU(imbg, hdr, name="NEB"))

    try:
        hdulist["NEB0"].data = improp
    except KeyError:
        hdulist.append(pyfits.ImageHDU(imbg0, hdr, name="NEB0"))

    hdulist.flush()
    # Save the fit parameters as well
    save_params(params, stamp_prefix + "-bgfit.json")


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
        "--ylo", type=float, nargs=2, default=[-6.5, -5.0],
        help="""Range of positions for lower BG sample"""
    )
    parser.add_argument(
        "--yhi", type=float, nargs=2, default=[5.0, 7.0],
        help="""Range of positions for upper BG sample"""
    )
    parser.add_argument(
        "--stampdir", type=str, default="Stamps",
        help="""Directory for placing the results"""
    )
    parser.add_argument(
        "--extra-suffix", type=str, default="",
        help="""Extra suffix for images (e.g., dd)"""
    )
    parser.add_argument(
        "--min-fraction", type=float, default=0.05,
        help="""Minimum fraction of peak brightness in order that a
        pixel should contribute to the velocity moments"""
    )
    cmd_args = parser.parse_args()
    main(**vars(cmd_args))
