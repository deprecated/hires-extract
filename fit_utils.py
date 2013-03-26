"""
Utility functions for fitting the line profiles

"""
import json
import numpy as np
from scipy.stats import norm
import lmfit
import astropy.io.fits as pyfits
from numpy.polynomial import Chebyshev as T


def gauss(u, area=1.0, u0=0.0, sigma=1.0, du=None):
    """Gaussian profile with total area under the profile of area,
    centered on u0 and with RMS width sigma

    Note that u, area, u0, sigma can all be 2-d arrays

    If optional argument du is present it is the velocity cell size
    and the profile is returned averaged over the cell

    """
    if du is None:
        # straightforward evaluation of the function
        return area*norm.pdf(u, loc=u0, scale=sigma)
    else:
        # average over the velocity cell
        return area*(norm.cdf(u+0.5*du, loc=u0, scale=sigma) -
                     norm.cdf(u-0.5*du, loc=u0, scale=sigma))/du


MAXORDER = 2
NPAR_PER_GAUSSIAN = (1 + MAXORDER)*3
# Maximum extent of y-axis - mapped onto [-1, 1] for Chebyshev polynomials
YDOMAIN = [-7.0, 7.0]


def model(U, Y, params, du=None):
    """Summation of various Gaussian components with parameters that are
    all Chebyshev polynomials of position

    """
    ngauss = len(params)/NPAR_PER_GAUSSIAN
    total = np.zeros_like(U)
    for ABC in "ABCDEFGH"[:ngauss]:
        # Unpack the Chebyshev polynomial coefficients
        # Lowest order polynomial is first in the coefficient list
        i_coeffs = [params["{}_i{}".format(ABC, k)].value
                    for k in range(MAXORDER+1)]
        u_coeffs = [params["{}_u{}".format(ABC, k)].value
                    for k in range(MAXORDER+1)]
        w_coeffs = [params["{}_w{}".format(ABC, k)].value
                    for k in range(MAXORDER+1)]
        i_p = T(i_coeffs, domain=YDOMAIN)
        u_p = T(u_coeffs, domain=YDOMAIN)
        w_p = T(w_coeffs, domain=YDOMAIN)
        total += gauss(U, i_p(Y), u_p(Y), w_p(Y), du)
    return total


def model_minus_data(params, u, y, data):
    """Function to minimize"""
    return model(u, y, params) - data


def init_single_component(params, ABC, i_coeffs, u_coeffs, w_coeffs):
    # k=0 corresponds to the constant term
    # k=n corresponds to the y**n yerm
    for k, coeff in enumerate(i_coeffs):
        params.add("{}_i{}".format(ABC, k), value=coeff)
    for k, coeff in enumerate(u_coeffs):
        params.add("{}_u{}".format(ABC, k), value=coeff)
    # Widths should only vary between about 2 and about 10
    params.add(ABC+"_w0", value=w_coeffs[0], min=2.0, max=10.0)
    # The variation should be even more restricted
    for k, coeff in enumerate(w_coeffs[1:], start=1):
        params.add("{}_w{}".format(ABC, k), value=coeff, min=-1.0, max=1.0)


def save_params(params, fn):
    """Save fit parameters params to filename fn

    """
    with open(fn, "w") as f:
        json.dump({k: v.value for k, v in params.items()}, f, indent=4)


if __name__ == "__main__":
    # Test set-up and save of simple model parameters
    params = lmfit.Parameters()
    init_single_component(params, "A",
                          i_coeffs=[2.0, 0.0, -1.5],
                          u_coeffs=[25.0, 2.0, 1.0],
                          w_coeffs=[5.0, 0.0, 2.0])
    init_single_component(params, "B",
                          i_coeffs=[1.0, 0.0, 0.0],
                          u_coeffs=[0.0, -0.5, 0.5],
                          w_coeffs=[2.0, 0.0, 0.3])
    save_params(params, "gauss-test.json")

    # Test with some sample data
    U, Y = np.meshgrid(np.linspace(-10.0, 50.0, 31),
                       np.linspace(-7.0, 7.0, 29))
    im = model(U, Y, params, du=2.0)
    pyfits.PrimaryHDU(im).writeto("gauss-test.fits", clobber=True)
