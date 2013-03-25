"""
Utility functions for fitting the line profiles

"""
import json
import numpy as np
from scipy.stats import norm
import lmfit


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


def model(U, Y, params, du=None):
    """Summation of various Gaussian components"""
    ngauss = len(params)/NPAR_PER_GAUSSIAN
    total = np.zeros_like(U)
    for ABC in "ABCDEFGH"[:ngauss]:
        i_coeffs = []
        u_coeffs = []
        w_coeffs = []
        # we need the highest power first
        for k in range(MAXORDER, -1, -1):
            # Unpack each component from the model parameters
            suff = "{}".format(k)
            i_coeffs.append(params[ABC + "_i" + suff].value)
            u_coeffs.append(params[ABC + "_u" + suff].value)
            w_coeffs.append(params[ABC + "_w" + suff].value)
        i_p = np.poly1d(i_coeffs)
        u_p = np.poly1d(u_coeffs)
        w_p = np.poly1d(w_coeffs)
        total += gauss(U, i_p(Y), u_p(Y), w_p(Y), du)
    return total


def model_minus_data(params, u, y, data):
    """Function to minimize"""
    return model(u, y, params) - data


def init_single_component(params, ABC, i_coeffs, u_coeffs, w_coeffs):
    # k=0 corresponds to the constant term
    # k=n corresponds to the y**n yerm
    for k, coeff in enumerate(reversed(i_coeffs)):
        params.add("{}_i{}".format(ABC, k), value=coeff)
    for k, coeff in enumerate(reversed(u_coeffs)):
        params.add("{}_u{}".format(ABC, k), value=coeff)
    for k, coeff in enumerate(reversed(w_coeffs)):
        params.add("{}_w{}".format(ABC, k), value=coeff)

if __name__ == "__main__":
    params = lmfit.Parameters()
    init_single_component(params, "A", i_coeffs=[0.0, 0.0, 1.0],
                          u_coeffs=[0.0, 0.0, 25.0], w_coeffs=[0.0, 0.0, 5.0])
    with open("gauss-test.json", "w") as f:
        json.dump({k: v.value for k, v in params.items()}, f, indent=4)














