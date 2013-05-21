import argparse
import json
import os
import numpy as np
import astropy.io.fits as pyfits
from doublet_utils import multiplet_moments, equivalent_doublet, \
    undoubletify, partition_doublet
DB_FILE_NAME = "multiplet-database.json"
LIGHT_SPEED_KMS = 2.99792458e5


def get_multiplet_data(stampdir):
    """Return the data for all the multiplets"""
    return json.load(open(os.path.join(stampdir, DB_FILE_NAME)))


def get_prefix(specid, lineid, stampdir="Stamps"):
    return os.path.join(stampdir, '-'.join([specid, lineid, "stamp", "nc"]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Deconvolve all the multiplet permitted lines"""
    )
    parser.add_argument(
        "specid", type=str,
        help="""Prefix of original FITS spectrum file (e.g., p84)"""
    )
    parser.add_argument(
        "--stampdir", type=str, default="Stamps",
        help="""Directory for placing the results"""
    )

    parser.add_argument(
        "--method", type=str, default="partition",
        choices=["partition", "moments"],
        help="""Which method to use for reducing the multiplet to a doublet"""
    )
    cmd_args = parser.parse_args()

    all_data = get_multiplet_data(cmd_args.stampdir)

    for lineid, data in all_data.items():
        prefix = get_prefix(cmd_args.specid, lineid, cmd_args.stampdir)
        hdulist = pyfits.open(prefix + ".fits")
        image = hdulist["SCI"].data
        variance = hdulist["SIGMA"].data**2
        du = hdulist["SCI"].header["CD1_1"]

        weights = [A*(2*J + 1) for A, J in zip(data["A_ki"], data["J_k"])]
        wavs = data["Wavelength"]
        moments = multiplet_moments(wavs, weights)
        # convert from wav (Angstrom) to velocity (km/s)
        sigma = moments["sigma"] * LIGHT_SPEED_KMS / moments["mean"]
        vels = [(wav - moments["mean"]) * LIGHT_SPEED_KMS /
                  moments["mean"] for wav in wavs]
        if cmd_args.method == "moments":
            doublet = equivalent_doublet(1.0, 0.0, sigma,  moments["skewness"])
        elif cmd_args.method == "partition":
            doublet = partition_doublet(vels, weights)
        
        newim, newvar = undoubletify(image, doublet["a"], doublet["delta"]/du, variance)
        hdulist["SCI"].data = newim
        hdulist["SIGMA"].data = np.sqrt(newvar)
        # Shift the velocity frame to be centered on the principal component
        for hdu in hdulist:
            if "CRVAL1" in hdu.header:
                hdu.header["CRVAL1"] -= doublet["x1"]
        hdulist.writeto(prefix + "-dd.fits", clobber=True)
