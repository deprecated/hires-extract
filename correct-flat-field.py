import argparse
import os
import astropy.io.fits as pyfits


def divide_by_flat(specid, quartzid, extractdir, interorder_divide):

    # Utility functions
    def order_file(iorder, specid, suffix1="o"):
        """Return the rectified image of the isolated order
        """
        return os.path.join(
            extractdir,
            "{}{}-order{}.fits".format(specid, suffix1, iorder)
        )

    orders = range(51, 77)
    for iorder in orders:
        print "Processing ", iorder
        hdulist = pyfits.open(order_file(iorder, specid))
        image = hdulist["SCI"].data
        sigma = hdulist["SIGMA"].data
        flat_image = pyfits.open(order_file(iorder, quartzid))["SCI"].data
        flat_image = flat_image[13:49,:].mean(axis=0)
        #  Make sure that the magnitude of the image is not changed too much
        if iorder == orders[0]:
            reddest_order_mean = flat_image.mean() 
        if interorder_divide:
            # Optionally attempt to scale all to the reddest order
            flat_image /= reddest_order_mean
        else:
            # But by default, just scale to mean within this order
            flat_image /= flat_image[-50:-10].mean() 
        image /= flat_image[None,:]
        sigma /= flat_image[None,:]
        hdulist.writeto(order_file(iorder, specid, suffix1="f"), clobber=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="""Divide by the quartz lamp image to correct for flat-field
        variations"""
    )
    parser.add_argument(
        "specid", type=str,
        help="""Prefix of object FITS spectrum file (e.g., p84)"""
    )
    parser.add_argument(
        "--quartzid", type=str, default="q69",
        help="""Prefix of quartz lamp FITS spectrum file"""
    )
    parser.add_argument(
        "--extractdir", type=str, default="Extract",
        help="""Directory containing the extracted orders"""
    )
    parser.add_argument(
        "--interorder-divide", action="store_true", 
        help="""Try to scale all the orders to brightness of the reddest order.
        Note that this is a bit of a stupid thing to do unless we
        think that the quartz lamp has a perfectly flat spectrum,
        which it surely hasn't"""
    )

    cmd_args = parser.parse_args()

    divide_by_flat(**vars(cmd_args))
