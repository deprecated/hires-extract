"""Utility functions for dealing with multiplets:

    multiplet_moments - Find velocity-like moments 0 to 3 of a
                        multiplet with known component intensities

    equivalent_doublet - Find the unique doublet that has a
                         given set of velocity-like moments 0 to 3

    undoubletify - Deconvolve a doublet by the subtract-and-shift method

"""

import yaml
import numpy as np
from scipy.interpolate import interp1d


def multiplet_moments(wavs, intensities):
    """Find the moments (intensity, mean, sigma, skewness) of a multiplet

    Parameters
    ----------
    wavs : array_like or scalar
         The wavelengths (or frequencies, or velocities, or ...)
         of each component of the multiplet
    intensities : array_like or scalar
         The correponding component intensities

    Returns
    -------
    moments : dict
         {intensity, mean, sigma, skewness}

    """

    ##
    ## This currently uses a single-pass method by calculating all the
    ## non-central moments M_k first.  An alternative would be to
    ## calculate just M_0 and M_1 first, then the mean, and then
    ## subtract the mean from the wavelength list before calculating
    ## the sigma and skewness.  That might be more numerically stable
    ## and less prone to over/underflow.  We can try it if this way
    ## has any problems.
    ##

    # normalize/check input
    wavs, intensities = np.atleast_1d(wavs, intensities)
    assert len(wavs) == len(intensities), (
        "Length mismatch between wavs and intensities")
    if len(wavs) == 1:
        return {"intensity": intensities[0], "mean": wavs[0], "sigma": 0.0,
                "skewness": 0.0}

    # First the "raw" moments, M_0 to M_3:
    M = [np.sum(intensities*wavs**k) for k in range(4)]

    # Now return the central and regularised moments
    intensity = M[0]
    mu = M[1]/M[0]
    sig = np.sqrt((M[2]/M[0]) - mu**2)
    skew = ((M[3]/M[0]) - 3*mu*sig**2 - mu**3)/sig**3

    return {"intensity": intensity, "mean": mu, "sigma": sig, "skewness": skew}


def equivalent_doublet(intensity, mean, sigma, skew):
    """Find the unique doublet of delta functions that has a given
    combination of intensity, mean, sigma, and skewness.

    """
    #
    # This uses the notation of "Moments of a pure doublet" in stis-lv.org
    # In particular, a is always < 1 and delta may be +/-
    #
    # Find the intensity ratio weak/strong
    a = 1.0 - np.abs(skew) * (np.sqrt(1.0 + 0.25*skew**2) - 0.5*np.abs(skew))
    # Find the velocity difference
    delta = skew*sigma*(1.0 + a)/(1.0 - a)
    # Find the velocities of the two components with respect to the mean
    x1 = -a*delta/(1.0 + a)
    x2 = delta/(1.0 + a)
    # Find the intensities of the two components
    I1 = intensity/(1.0 + a)
    I2 = a*I1
    # Just return everything in case it is needed
    return {
        "x1": x1, "x2": x2, "I1": I1, "I2": I2, "a": a, "delta": delta,
        "v1": x1 + mean, "v2": x2 + mean,
    }


def partition_doublet(wavs, intensities):
    """Partition a multiplet into two groups, either side of the largest gap

    Return the doublet that corresponds to treating each group as a
    single component.  Unlike the equivalent_doublet function, each
    component of the doublet has its own internal sigma, but hopefully
    these are much smaller than the width of the entire multiplet.

    """
    assert len(wavs) == len(intensities), (
        "Length mismatch between wavs and intensities")
    # Make sure wavelengths are in ascending order
    wavs, intensities = np.atleast_1d(*zip(*sorted(zip(wavs, intensities))))
    # Find where the split the multiplet in two
    isplit = 1 + np.argmax(np.diff(wavs))
    wavs1, wavs2 = wavs[:isplit], wavs[isplit:]
    intensities1, intensities2 = intensities[:isplit], intensities[isplit:]
    moms1 = multiplet_moments(wavs1, intensities1)
    moms2 = multiplet_moments(wavs2, intensities2)
    
    I1, I2 = moms1["intensity"], moms2["intensity"]
    a = I2/I1
    v1, v2 = moms1["mean"], moms2["mean"]
    delta = v2 - v1
    vmean = (v1 + a*v2)/(1.0 + a)
    x1, x2 = v1 - vmean, v2 - vmean
    if I1 > I2:
        return {
            "x1": x1, "x2": x2, "I1": I1, "I2": I2, "a": a, "delta": delta,
            "v1": v1, "v2": v2, "sig1": moms1["sigma"], "sig2": moms2["sigma"],
            "skew1": moms1["skewness"], "skew2": moms2["skewness"],
        }
    else:
        return {
            "x1": x2, "x2": x1, "I1": I2, "I2": I1,
            "a": 1.0/a, "delta": -delta,
            "v1": v2, "v2": v1, "sig1": moms2["sigma"], "sig2": moms1["sigma"],
            "skew1": moms2["skewness"], "skew2": moms1["skewness"],
        }


def undoubletify(spectrum, a, x0):
    """Deconvolve a doublet line profile using repeated subtraction
    method

    Input arguments:

                    spectrum - 1D or 2D array of line profile
                               Dispersion is along the x axis (last dimension)

                           a - Intensity ratio of the two doublet components
                               (right over left)

                          x0 - Separation in pixels (or fraction) between the
                               two components.

    Returns:
                 newspectrum - Deconvolved spectrum

    """

    def shift_right(array, xshift):
        """Shift an array by xshift pixels to the right along the last axis
        with linear interpolation and infinite zero padding from the
        left

        """
        f = interp1d(x, array, bounds_error=False, fill_value=0.0)
        # The array is shifted right, so the coords must be shifted left
        return f(x - xshift)

    # Make sure that a is less than unity ...
    if a > 1.0:
        # ... by swapping over the two components if necessary
        a = 1./a
        x0 = -x0

    print "Doublet parameters: a = {}, x0 = {}".format(a, x0)

    # Construct an array of nx pixel values along the dispersion direction
    nx = spectrum.shape[-1]
    x = np.arange(nx)
    # We can truncate the series once the entire profile has
    # shifted off the grid
    nterms = abs(int(nx/x0))
    print "Using ", nterms, " terms"
    # Calculate the series
    newspectrum = spectrum
    for n in 1 + np.arange(nterms):
        newspectrum += (-a)**n * shift_right(spectrum, n*x0)

    return newspectrum


def pprint_dict(d):
    """Pretty print a dict that may have numpy scalars in it"""
    return yaml.dump({k: float(v) for k, v in d.items()},
                     default_flow_style=False)


if __name__ == "__main__":
    # Test mechanism on an [O I] triplet with common upper level
    #
    # | J_i | J_k |      Wav |     A_ki |         I |    dV | I dV**2 |
    # |-----+-----+----------+----------+-----------+-------+---------|
    # |   1 |   1 | 6046.233 | 1.04E+06 | 0.3324808 | -7.09 |   16.71 |
    # |   2 |   1 | 6046.438 | 1.74E+06 | 0.5562660 |  3.07 |    5.24 |
    # |   0 |   1 | 6046.495 | 3.48E+05 | 0.1112532 |  5.90 |    3.87 |
    # |-----+-----+----------+----------+-----------+-------+---------|
    # |     |     | 6046.376 | 3128000. | 1.0000000 |  0.00 |    5.08 |

    triplet_wavs = [6046.233, 6046.438, 6046.495]
    triplet_weights = [1.04E+06, 1.74E+06, 3.48E+05]
    light_speed_kms = 2.99792458e5

    moments = multiplet_moments(triplet_wavs, triplet_weights)
    print "O I 6046 multiplet:"
    print "-------------------"
    print pprint_dict(moments)

    intensity = 1.0
    mean = 0.0
    # convert to km/s
    sigma = moments["sigma"] * light_speed_kms / moments["mean"]
    assert (sigma - 5.08) < 0.01, "Sigma = {} is not right".format(sigma)

    doublet = equivalent_doublet(intensity, mean, sigma,  moments["skewness"])
    print "Equivalent doublet"
    print pprint_dict(doublet)

    triplet_vels = [(wav - moments["mean"]) * light_speed_kms /
                    moments["mean"] for wav in triplet_wavs]
    triplet_weights = [w/moments["intensity"] for w in triplet_weights]
    doublet = partition_doublet(triplet_vels, triplet_weights)
    print "Partitioned doublet"
    print pprint_dict(doublet)
    
    # Test mechanism on an [O I] sextuplet, assuming the upper levels are
    # distributed according to statistical weight 2 J + 1

    # | J_i | J_k |         Wav |     A_ki |         I |    dV | I dV**2 |
    # |-----+-----+-------------+----------+-----------+-------+---------|
    # |   1 |   1 |    7001.899 | 1.39E+06 | 0.0834534 | -9.65 |    7.77 |
    # |   1 |   2 |    7001.922 | 2.50E+06 | 0.2501601 | -8.67 |   18.80 |
    # |   2 |   1 |    7002.173 | 9.26E+04 | 0.0055636 |  2.08 |    0.02 |
    # |   2 |   2 |    7002.196 | 8.33E+05 | 0.0834534 |  3.07 |    0.79 |
    # |   2 |   3 |    7002.230 | 3.33E+06 | 0.4662984 |  4.52 |    9.53 |
    # |   0 |   1 |    7002.250 | 1.85E+06 | 0.1110711 |  5.38 |    3.21 |
    # |-----+-----+-------------+----------+-----------+-------+---------|
    # |     |     | 7002.124395 | 9.9956e6 | 1.0000000 |  0.00 |    6.33 |

    sextuplet_wavs = [7001.899, 7001.922, 7002.173,
                      7002.196, 7002.230, 7002.250]
    sextuplet_A21s = [1.39E+06, 2.50E+06, 9.26E+04,
                      8.33E+05, 3.33E+06, 1.85E+06]
    sextuplet_J2s = [1,  2,  1,  2,  3,  1]
    sextuplet_weights = [A*(2*J + 1) for A, J in
                         zip(sextuplet_A21s, sextuplet_J2s)]

    moments = multiplet_moments(sextuplet_wavs, sextuplet_weights)
    print "O I 7002 multiplet:"
    print "-------------------"
    print pprint_dict(moments)

    intensity = 1.0
    mean = 0.0
    # convert to km/s
    sigma = moments["sigma"] * light_speed_kms / moments["mean"]
    assert (sigma - 6.33) < 0.01, "Sigma = {} is not right".format(sigma)
    doublet = equivalent_doublet(intensity, mean, sigma,  moments["skewness"])
    print "Equivalent doublet"
    print pprint_dict(doublet)

    sextuplet_vels = [(wav - moments["mean"]) * light_speed_kms /
                      moments["mean"] for wav in sextuplet_wavs]
    sextuplet_weights = [w/moments["intensity"] for w in sextuplet_weights]
    doublet = partition_doublet(sextuplet_vels, sextuplet_weights)
    print "Partitioned doublet"
    print pprint_dict(doublet)
