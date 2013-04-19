import argparse
import glob
import astropy.io.fits as pyfits


def remove_cosmetic_defects(fits_file_path):
    hdulist = pyfits.open(fits_file_path, mode='update')
    data = hdulist[0].data
    for i in 983, 1149, 1152:
        data[:, i-1] = 0.5*(data[:, i] + data[:, i-2])
    data[:, 106] = 0.5*(data[:, 105] + data[:, 108])
    data[:, 107] = 0.5*(data[:, 105] + data[:, 108])
    data[:, 2027] = 0.5*(data[:, 2026] + data[:, 2029])
    data[:, 2028] = 0.5*(data[:, 2026] + data[:, 2029])
    hdulist.flush()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Remove cosmetic defects from the chips"""
    )
    parser.add_argument(
        "pattern", type=str,
        help="""Glob pattern for files to process (e.g., Keck?/p??b)"""
    )
    cmd_args = parser.parse_args()
    for fits_file in glob.glob(cmd_args.pattern):
        remove_cosmetic_defects(fits_file)
