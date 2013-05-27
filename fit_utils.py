"""
Utility functions for fitting the line profiles

"""
import json
import numpy as np
from scipy.stats import norm
import lmfit
import astropy.io.fits as pyfits
from astropy.table import Table
from numpy.polynomial import Chebyshev as T
import bottleneck as bn


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


def model_minus_data(params, u, y, data, du=None):
    """Function to minimize"""
    return model(u, y, params, du) - data


def model_minus_data_over_sigma(params, u, y, data, sigma, du=None):
    """Another function to minimize"""
    return (model(u, y, params, du) - data)/sigma


def init_single_component(params, ABC, i_coeffs, u_coeffs, w_coeffs, min_width=3.0):
    # k=0 corresponds to the constant term
    # k=n corresponds to the y**n yerm
    eps = 1.e-8
    for k, coeff in enumerate(i_coeffs):
        coeff_id = "{}_i{}".format(ABC, k)
        if coeff is None:
            params.add(coeff_id, value=eps, vary=False)
        else:
            params.add(coeff_id, value=coeff)
    for k, coeff in enumerate(u_coeffs):
        coeff_id = "{}_u{}".format(ABC, k)
        if coeff is None:
            params.add(coeff_id, value=eps, vary=False)
        else:
            params.add(coeff_id, value=coeff)
    # Widths should only vary between about 2 and about 15
    params.add(ABC+"_w0", value=w_coeffs[0], min=min_width, max=15.0)
    # The variation should be even more restricted
    for k, coeff in enumerate(w_coeffs[1:], start=1):
        coeff_id = "{}_w{}".format(ABC, k)
        if coeff is None:
            params.add(coeff_id, value=eps, vary=False)
        else:
            params.add(coeff_id, value=coeff, min=-1.0, max=1.0)


def save_params(params, fn):
    """Save fit parameters params to filename fn

    """
    with open(fn, "w") as f:
        json.dump({k: v.value for k, v in params.items()}, f, indent=4)




def std_from_model_fuzzing(U, Y, params, du=None, nsamp=100,
                           debug_prefix=False, full_covar=False):
    """Estimate std of nebular pv image due to std of fit parameters

    Uses a Monte Carlo simulation of nsamp realizations of the model,
    with parameters drawn from Gaussian distributions around the
    best-fit values, with widths equal to the reported stderror.

    Currently does not attempt to make use of the correlations between
    model parameters, which means we may overestimate the
    uncertainties....

    """
    ny, nu = U.shape
    model_stack = np.empty((nsamp, ny, nu))
    scaled_means = [p.value/find_param_scale(params, n) for n, p in params.items()]
    scaled_covar = calculate_covar_array(params)
    fuzzy_params_stack = []
    # Fill in a stack of nebular models, all fuzzed around the best fit
    for i in range(nsamp):
        fuzzy_params = lmfit.Parameters()
        fuzzy_scaled_values = np.random.multivariate_normal(scaled_means, scaled_covar)
        for (name, param), fuzzy_scaled_value in zip(params.items(), fuzzy_scaled_values):
            if param.vary and param.stderr > 0.0:
                if full_covar:
                    fuzzy_value = fuzzy_scaled_value*find_param_scale(params, name)
                else:
                    fuzzy_value = np.random.normal(param.value, param.stderr)
            else:
                # pegged parameter does not vary
                fuzzy_value = param.value
            # Ensure we do not stray outside of the established bounds
            if param.max:
                fuzzy_value = min(fuzzy_value, param.max)
            if param.min: 
                fuzzy_value = max(fuzzy_value, param.min)
            fuzzy_params.add(name, value=fuzzy_value)
        model_stack[i, :, :] = model(U, Y, fuzzy_params, du)
        fuzzy_params_stack.append({k: v.value for k, v in fuzzy_params.items()})
    if debug_prefix:
        pyfits.PrimaryHDU(model_stack).writeto(
            debug_prefix + "_model_stack.fits", clobber=True)
        with open(debug_prefix + "_model_stack.tab", "w") as f:
            f.write("\n".join(
                Table(fuzzy_params_stack).pformat(max_lines=-1, max_width=-1)))
        # Table(fuzzy_params_stack).write("debug_model_stack.tab", format="ascii")
    return bn.nanstd(model_stack, axis=0)
    

def find_param_scale(params, name):
    """In order to use the covariance array sensibly, we need all the
    parameters to have the same order of magnitude.  Otherwise, the
    smaller of the two is affected a lot more by the covariance.

    So, we use [ABC]_i0 as the scale for the intensities.  And we use
    5 km/s as the scale for velocities and widths.
    """
    if "_i" in name:
        ABC = name[0]
        scale = params[ABC+"_i0"].value
    else:
        scale = 5.0
    return scale


def calculate_covar_array(params):
    """Find the covariance array between each parameter in params

    Calculated from the stderr and correlations that lmfit reports
    """
    npar = len(params)
    parnames = params.keys()
    covar = np.zeros((npar, npar))
    for i, iname in enumerate(parnames):
        iscale = find_param_scale(params, iname)
        covar[i, i] = (params[iname].stderr/iscale)**2
        for j, jname in enumerate(parnames):
            if j != i:
                jscale = find_param_scale(params, jname)
                if params[iname].correl:
                    correl = params[iname].correl.get(jname, 0.0)
                else:
                    correl = 0.0
                covar[i, j] = correl*(
                    params[iname].stderr/iscale)*(
                    params[jname].stderr/jscale)
    return covar


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
