import argparse
import os
import glob
from remove_continuum import remove_continuum
DB_FILE_NAME = "line-database.json"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Fit and remove continuum from each line"""
    )
    parser.add_argument(
        "specid", type=str,
        help="""Prefix of original FITS spectrum file (e.g., p84)"""
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

    # Read in the line database
    stampfiles = glob.glob(
        os.path.join(cmd_args.stampdir, cmd_args.specid) + "-*-stamp.fits")
    for stampfile in stampfiles:
        stampname = '-'.join(
            os.path.splitext(
                os.path.basename(stampfile)
            )[0].split('-')[:-1]
        )
        print "Removing continuum from ", stampname
        remove_continuum(stampname, cmd_args.vrange, cmd_args.stampdir)
