import argparse
import os
import astropy.io.fits as pyfits
import numpy as np

def pixel_noise(specid, extractdir, gain, readout):
    
    # Utility functions
    def order_file(iorder, specid, suffix1="b-cr"):
        """Return the rectified image of the isolated order
        """
        return os.path.join(
            extractdir,
            "{}{}-order{}.fits".format(specid, suffix1, iorder)
        )

    orders = range(51, 77)
    for iorder in orders:
        print "Processing ", iorder
        hdulist = pyfits.open(order_file(iorder, specid), mode='update')
        image = hdulist["SCI"].data
        sig2 = image/gain + (readout/gain)**2
        hdulist.append(pyfits.ImageHDU(np.sqrt(sig2), name="SIGMA"))
        hdulist.flush(output_verify="fix")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Calculate standard deviation of each pixel
        
        Takes into account Poisson noise plus readout noise.
        Should be applied to the bias-subtracted extracted orders. 
        Writes an extra HDU ("SIGMA") to the FITS file, 
        which contains an image of the stdev. 
        """
    )
    parser.add_argument(
        "specid", type=str,
        help="""Prefix of object FITS spectrum file (e.g., p84)"""
    )
    parser.add_argument(
        "--extractdir", type=str, default="Extract",
        help="""Directory containing the extracted orders"""
    )
    parser.add_argument(
        "--gain", type=float, default=4.0,
        help="""CCD inverse gain in elec/ADU"""
    )
    parser.add_argument(
        "--readout", type=float, default=5.0,
        help="""CCD readout noise in electrons"""
    )

    cmd_args = parser.parse_args()
    pixel_noise(**vars(cmd_args))
