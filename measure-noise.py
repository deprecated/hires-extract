"""
Measure noise in the flat field images, and plot noise against signal level
"""
import astropy.io.fits as pyfits
import os
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp

# GAIN = 2.38


def flat_image(specid="q69", folder="Keck1", dark=5.5, full=False):
    """Return bias-subtracted original image

    Optionally, also return exposure time is full=True
    """
    hdu, = pyfits.open(
        os.path.join(folder, specid + ".fits")
    )
    bzero = hdu.header.get("BZERO", 0.0)
    bscale = hdu.header.get("BSCALE", 1.0)
    print bzero, bscale
    scaled_data = (hdu.data - bzero)/bscale - dark
    if full:
        exposure_time = hdu.header.get("EXPOSURE", 0.0)
        return scaled_data, exposure_time
    else:
        return scaled_data


def fsigma(N1, M=1.0, delta0=0.0, R0=1.0):
    """Calculate Poisson and readout noise contributions to sigma(R)/R

    Where R = N1/N2 is the ratio of two "matched" images (N1, N2) and
    sigma(R) is the standard deviation of R.

    The parameters are ADU gain, M (e-/DN), RMS readout noise, delta0
    (e-), and the expectation value of R, R0 (calculated as ratio of
    exposure times).  By "matched", I mean that the images should be
    close to identical to within a constant factor, apart from the
    noise, obviously.

    Note that this function uses the notation that I use in the docs
    (see keck-revisited.org), unlike the rest of this file.  In
    particular, N1 (brightness of the longer exposure) actually
    corresponds to file2 and image2 in plot_statistics()

    """
    return np.sqrt(
        (1.0 + R0)/(N1*M) + 
        (1.0 + R0**2)*(delta0**2)/(N1*M)**2
    )


def robust_statistics(x, y, xedges):
    """Calculate robust estimates of location and scale of a distribution

    Returns (loc, scale) of y, binned by x according to xedges

    Returns vectors of length len(xedges)-1 

    loc is the "average" value, estimated as the trimean

    scale is the width of the distribution, estimated from the
    interquartile range (IQR) and rescaled to be equal to the standard
    deviation, sigma, for a Gaussian distribution

    The IQR and trimean are much less affected by outliers and
    power-law wings than are the regular mean and sigma

    """
    # Interquartile range in units of std for Gaussian profile
    iqr_over_sigma = 2.0*np.sqrt(2.0)*sp.erfinv(0.5)
    trimean, iqr = [], []
    # Loop over bins in x
    for x1, x2 in zip(xedges[:-1], xedges[1:]):
        # Construct a mask that selects this bin
        m = (x >= x1) & (x < x2)
        # Find the 1st, 2nd, 3rd quartiles
        q1 = np.percentile(y[m], 25.0)
        q2 = np.percentile(y[m], 50.0)
        q3 = np.percentile(y[m], 75.0)
        # The trimean is the weighted average of the median and the
        # two quartiles (http://en.wikipedia.org/wiki/Trimean)
        trimean.append(0.25*(q1 + 2*q2 + q3))
        # The interquartile range is the difference between quartiles
        iqr.append(q3 - q1) 
    # Convert lists to numpy arrays before returning
    loc = np.array(trimean)
    # Put scale in units of Gaussian standard deviation
    scale = np.array(iqr)/iqr_over_sigma
    return loc, scale 


def plot_statistics(file1="q69", file2="q110",
                    gain=4.4, dy=0.05, xmax=None, bins=200, robust=True, max_deviation=5.0):

    """Plot noise statistics from ratio of two images
    
    The two images should be as close to identical as possible (e.g.,
    two quartz lamp exposures, or two science exposures of the same
    object at same PA and position), but they may have different
    exposure times.

    """
    if "/" in file1:
        folder, file1 = file1.split("/")
    else:
        folder = "Keck1"
    image1, texp1 = flat_image(file1, folder=folder, full=True)

    if "/" in file2:
        folder, file2 = file2.split("/")
    else:
        folder = "Keck1"
    image2, texp2 = flat_image(file2, folder=folder, full=True)

    ratio = image2/image1
    R0 = texp2/texp1            # Ratio of exposure times

    mean, median = ratio.mean(), np.median(ratio)
    print "Average image ratio: mean, median = ", mean, median
    print "Exposure time ratio: ", R0

    # Flattened versions of the data to analyse: 
    #
    # x-axis is the brightness of the longer exposure image
    x = np.log10(image2.ravel())
    # y-axis is the ratio of the two images divided by the ratio of
    # exposure times
    y = ratio.ravel()/R0

    xmin = 1.5
    if xmax is None:
        xmax = x.max()
    else:
        xmax = np.log10(xmax)
    ymin, ymax = 1 - dy, 1 + dy

    H, xedges, yedges = np.histogram2d(
        x, y, bins=bins, range=[[xmin, xmax], [ymin, ymax]],
        normed=False)
    # Normalize to the sum for each x-value
    H /= H.sum(axis=1)[:, None]

    # Fix-up Nans
    H[~np.isfinite(H)] = 0.0

    # Make a 2d grid of pixel center x, y-values for the histogram array
    ygrid = 0.5*(yedges[:-1] + yedges[1:])
    xgrid = 0.5*(xedges[:-1] + xedges[1:])
    Ygrid, Xgrid = np.meshgrid(ygrid, xgrid)

    # Assuming we know the gain, then we can predict the std of the ratio
    # image as a function of brightness

    readout = 5.0
    Hsig_theory = fsigma(10**xgrid, M=gain, delta0=readout, R0=R0)

    if robust:
        # Restrict calculation to the graph viewport
        m = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
        print y.shape, m.shape, m.sum()
        # Reject any points that are more than 5 sig away from the nominal R0
        m = m & (np.abs(y - 1.0) < max_deviation*fsigma(10**x, gain, readout, R0))
        print y.shape, m.shape, m.sum()
        Hmean, Hsig_obs = robust_statistics(x[m], y[m], xedges)
    else:
        # H is already normalized so that integral along y-axis is unity
        Hmean = np.sum(Ygrid*H, axis=1)
        dY = Ygrid - Hmean[:, None]
        Hsig_obs = np.sqrt(np.sum((dY**2)*H, axis=1))


    np.set_printoptions(precision=6)
    print H.T[::20, ::20]
    print xmin, xmax, ymin, ymax

    # Make a graph of the 2d histogram and the fits
    plt.clf()
    plt.imshow(H.T, extent=[xmin, xmax, ymin, ymax],
               cmap=plt.cm.gray_r,
               origin="low", aspect="auto", interpolation="nearest")
    plt.plot(xgrid, Hmean, label="Mean")
    plt.plot(xgrid, Hmean + Hsig_theory, '--b', 
             label="Theoretical sigma, gain = {:.1f}".format(gain))
    plt.plot(xgrid, Hmean - Hsig_theory, '--b')
    plt.plot(xgrid, Hmean + Hsig_obs, '--r',
             label="Observed sigma")
    plt.plot(xgrid, Hmean - Hsig_obs, '--r')
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel("log_10(Brightness of {})".format(file2))
    plt.ylabel("Ratio ({1}/{0}) / {2:.2f}".format(file1, file2, R0))
    plt.title(

        "Empirical determination of gain from ratio {1}/{0}".format(
            file1, file2))
    plt.legend(loc="lower right", fontsize="small")
    plt.savefig("measure-noise-{}-{}.png".format(file1, file2), 
                dpi=150)
    plt.savefig("measure-noise-{}-{}.pdf".format(file1, file2))

    # Make another graph of sigma against sqrt(1/brightness)
    plt.clf()
    plt.plot(10**xgrid, Hsig_obs, "om", label="Observed", ms=10.0, alpha=0.5)
    for gfac, rfac, style in (1.0/1.5, 1.0, "-r"), (1.0, 1.0, "-g"), (1.0, 0.6, ":g"), (1.5, 1.0, "-b"):
        plt.loglog(10**xgrid, fsigma(10**xgrid, gain*gfac, readout*rfac, R0), style,
                   label="Gain = {:.1f} e-/DN, readout = {:.1f} e-".format(
                       gain*gfac, readout*rfac),
                   lw=3, alpha=0.7
        )
    plt.xlabel("N({}), DN".format(file2))
    plt.ylabel("Standard deviation of N({})/N({})".format(file2, file1))
    plt.xlim(20.0, 3.e4)
    plt.ylim(0.004, 0.5)
    plt.legend(loc="lower left", fontsize="small")
    plt.savefig("measure-gain-{}-{}.png".format(file1, file2), 
                dpi=150)
    plt.savefig("measure-gain-{}-{}.pdf".format(file1, file2))
    

if __name__ == "__main__":
    # Flat fields
    plot_statistics("q69", "q110", xmax=21500.0,
                    dy=1., bins=200, gain=6.0)

    # 244-440
    plot_statistics("p83", "p84", xmax=5000.0,
                    dy=1.0, bins=(100, 100), gain=3.5)

    # 182-413
    plot_statistics("p78", "p77", xmax=5000.0,
                    dy=1.0, bins=(100, 100), gain=3.5)
    plot_statistics("p86", "p87", xmax=5000.0,
                    dy=1.0, bins=(100, 100), gain=3.5)

    # 170-337 - doesn't work at all (as in: crashes)
    # plot_statistics("p71", "p72", xmax=3000.0,
    #                 dy=1.0, bins=50, gain=2.4, robust=True)


    #
    # The technique does not work very well with images from different
    # nights - not similar enough
    #
    
    # 177-341 - different nights - makes two strands
    plot_statistics("Keck2/p59", "p75", xmax=3000.0,
                    dy=1.0, bins=(100, 100), gain=3.5, max_deviation=3.0)

    # plot_statistics("Keck2/p64", "p85", xmax=3000.0,
    #                 dy=1.0, bins=50, gain=2.4, robust=True)











