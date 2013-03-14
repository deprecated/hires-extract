"""
Construct a parametrized model of the quartz lamp
illumination image.  Fit it to the real image of each order.
"""

import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf
import lmfit
import json


def model_minus_data(params, y, data):
    """Function to minimize"""

    def window(y, y1, y2, sigma):
        """Window function with Gaussian broadening"""
        return (0.5 * (erf((y - y1)/sigma) + 1.0) *
                0.5 * (erf((y2 - y)/sigma) + 1.0))

    def model(y, params):
        """Summation of various window functions"""
        nwin = len(params)/4
        total = np.zeros_like(y)
        total += params["const"].value + params["linear"].value*y
        for iwin in range(1, nwin+1):
            suff = "_{:02d}".format(iwin)
            height = params["height" + suff].value
            y1 = params["y1" + suff].value
            w = params["w" + suff].value
            y2 = y1 + w
            sigma = params["sigma" + suff].value
            total += height*window(y, y1, y2, sigma)
        return total

    return model(y, params) - data


def init_window_params(params, iwin, height, y1,
                       w=36.0, sigma=1.0, vary_height=True, sigmax=2.0):
    """Initial parameters for a single window"""
    suff = "_{:02d}".format(iwin)
    params.add("height" + suff, value=height, min=0.0, vary=vary_height)
    # we always allow y1 to vary
    params.add("y1" + suff, value=y1, max=51.0, min=-1.1*w)
    # Window width is always kept fixed when the window
    # is not wholly within the strip
    params.add("w" + suff, value=w, min=30.0, max=40.0,
               vary=(vary_height and (y1 < 20.0) and (y1 > -5.0)))
    params.add("sigma" + suff, value=sigma,
               min=0.5, max=sigmax)


def prefix(iorder):
    return "Extract/q69b-order{}".format(iorder)

INITIAL_GUESSES = {
    51: ((5000.0, 4.), None, None),
    52: ((5500.0, 4.), None, None),
    53: ((6000.0, 4.), None, None),
    54: ((5800.0, 4.), None, None),
    55: ((5800.0, 4.), None, None),
    56: ((5800.0, 4.), None, None),
    57: ((5800.0, 4.), (5800.0, 48.), None),
    58: ((5800.0, 4.), (5800.0, 48.), None),
    59: ((4300.0, 4.), (4000.0, 47.), (4000.0, -36.)),
    60: ((4000.0, 4.), (4000.0, 46.), (4000.0, -30.)),
    61: ((4000.0, 4.), (4000.0, 42.), (4000.0, -30.)),
    62: ((3200.0, 3.), (4000.0, 41.), (3000.0, -33.)),
    63: ((3000.0, 4.), (3000.0, 41.), (3000.0, -28.)),
    64: ((3000.0, 4.), (3000.0, 40.), (3000.0, -27.)),
    65: ((2500.0, 4.), (2500.0, 40.), (2500.0, -27.)),
    66: ((2000.0, 4.), (2000.0, 39.), (2000.0, -26.)),
    67: ((2000.0, 4.), (2000.0, 39.), (2000.0, -26.)),
    68: ((2000.0, 4.), (2000.0, 38.), (2000.0, -25.)),
    69: ((1500.0, 4.), (1500.0, 38.), (1500.0, -25.)),
    70: ((1500.0, 4.), (1500.0, 37.), (1500.0, -24.)),
    71: ((1500.0, 4.), (1500.0, 36.), (1500.0, -24.)),
    72: ((1000.0, 4.), (1000.0, 35.), (1000.0, -23.)),
    73: ((1000.0, 4.), (1000.0, 34.), (1000.0, -23.)),
    74: ((1000.0, 4.), (1000.0, 33.), (1000.0, -22.)),
    75: ((400.0, 4.),  (400.0, 32.),  (400.0, -22.)),
    76: ((100.0, 4.),  (100.0, 31.),  None),
}

INITIAL_CONST = {74: 40.0, 75: 30.0, 76: 20.0}

if __name__ == "__main__":
    win2_save = (0.0, 0.0)

    orders = {}
    for iorder in range(51, 77):
        # Read in image data and construct 1D profile
        im = pyfits.open(prefix(iorder) + ".fits")["SCI"].data
        ny, nx = im.shape
        y = 1.0 + np.arange(ny)
        if iorder in [60, 61]:
            # Cut out the bad spot in the center
            im = np.hstack((im[:, :930], im[:, 1100:]))
        yprofile = im.mean(axis=1)
        # Mask out the NaNs
        m = np.isfinite(yprofile)
        y, yprofile = y[m], yprofile[m]

        # Initialize model with up to 3 windows
        params = lmfit.Parameters()
        guess_const = INITIAL_CONST.get(iorder, 180.0)
        params.add("const", value=guess_const,
                   min=guess_const/3.0, max=200.0)
        params.add("linear", value=0.0,
                   min=-1.0, max=1.0)
        win1, win2, win3 = INITIAL_GUESSES[iorder]
        height, y1 = win1
        init_window_params(params, 1, height, y1)
        if win2 is not None:
            height, y1 = win2
            # Now always fix win2 height and sigma from the previous order
            height, sigma = win2_save
            if iorder in [63, 65]:
                # More stupid special-casing
                sigmax = 1.3
            else:
                sigmax = 2.0
            init_window_params(params, 2, height, y1,
                               sigma=sigma, vary_height=False, sigmax=sigmax)
        if win3 is not None:
            height, y1 = win3
            if iorder == 59:
                # Special treatment for this one
                # - we cheat and look ahead to order 60
                init_window_params(params, 3, 3821.0, y1,
                                   sigma=1.06, vary_height=False)
            else:
                init_window_params(params, 3, height, y1)

        # Fit model to profile
        result = lmfit.minimize(model_minus_data, params,
                                args=(y, yprofile))

        # Save height of win1 to use as win2 in next order
        win2_save = (params["height_01"].value,
                     params["sigma_01"].value)

        # Write error report
        print "\n\n**** Order ", iorder
        lmfit.report_errors(params)

        # Plot results
        plt.plot(y, yprofile, "o")
        plt.plot(y, yprofile + result.residual)
        plt.xlim(0.0, 53.0)
        plt.ylim(0.0, 1.1*yprofile.max())
        plt.grid()
        plt.xlabel("J pixel (1-based)")
        plt.ylabel("Lamp intensity")
        plt.title("Order {}".format(iorder))
        plt.savefig(prefix(iorder) + "-fit.pdf")
        plt.clf()

        # Save result
        orders[iorder] = dict([(k, v.value) for k, v in params.items()])

    with open("quartz-database.json", "w") as f:
        json.dump(orders, f, indent=4)
