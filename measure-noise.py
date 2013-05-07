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


def robust_statistics(x, y, xedges):
    """Calculate robust estimates of location and scale of a distribution

    Returns a vector of length len(xedges)-1 

    Returns (loc, scale) of y, binned by x according to xedges

    loc is the "average" value, estimated from the trimean

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
                    gain=4.4, dy=0.05, xmax=None, bins=200, robust=False):

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
    x = image2.ravel()
    # y-axis is the ratio of the two images divided by the ratio of
    # exposure times
    y = ratio.ravel()/R0

    xmin = 0.0
    if xmax is None:
        xmax = image2.max()
    ymin, ymax = 1 - dy, 1 + dy

    # Old plot that no longer is needed
    # every = 1
    # plt.clf()
    # plt.plot(x[::every], y[::every], '.', alpha=0.003)
    # plt.plot([xmin, xmax], [1.0, 1.0])
    # plt.xlim(xmin, xmax)
    # plt.ylim(ymin, ymax)
    # plt.savefig("measure-noise.png")

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

    if robust:
        m = (x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)
        Hmean, Hsig_obs = robust_statistics(x[m], y[m], xedges)
    else:
        # H is already normalized so that integral along y-axis is unity
        Hmean = np.sum(Ygrid*H, axis=1)
        dY = Ygrid - Hmean[:, None]
        Hsig_obs = np.sqrt(np.sum((dY**2)*H, axis=1))

    # Assuming we know the gain, then we can predict the std of the ratio
    # image as a function of brightness

    readout_noise_electrons = 5.0
    Hsig_theory = np.sqrt(
        (1.0 + R0)/(xgrid*gain) + 
        (1.0 + R0**2)*(readout_noise_electrons**2)/(xgrid*gain)**2
    )

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
    plt.xlabel("Brightness of {}".format(file2))
    plt.ylabel("Ratio ({1}/{0}) / {2:.2f}".format(file1, file2, R0))
    plt.title(

        "Empirical determination of gain from ratio {1}/{0}".format(
            file1, file2))
    plt.legend(loc="lower right", fontsize="small")
    plt.savefig("measure-noise-{}-{}.png".format(file1, file2), 
                dpi=150)

    # Make another graph of sigma against sqrt(1/brightness)
    plt.clf()
    plt.plot(1/np.sqrt(xgrid), Hsig_obs, "o", label="Observed")
    for factor in 0.5, 1.0, 1.5:
        plt.loglog(1/np.sqrt(xgrid), Hsig_theory/np.sqrt(factor), "-",
                 label="Gain = {:.1f}".format(gain*factor))
    plt.xlabel("1 / sqrt[ N({}) ]".format(file2))
    plt.ylabel("Standard deviation of N({})/N({})".format(file2, file1))
    plt.xlim(0.003, 1.0)
    plt.ylim(0.003, 1.0)
    plt.legend(loc="lower right", fontsize="small")
    plt.savefig("measure-gain-{}-{}.png".format(file1, file2), 
                dpi=150)
    

if __name__ == "__main__":
    # Flat fields
    plot_statistics("q69", "q110", xmax=21500.0,
                    dy=0.2, bins=200, robust=True)

    # 244-440
    plot_statistics("p83", "p84", xmax=3000.0,
                    dy=0.4, bins=50, gain=2.4, robust=True)

    # 182-413
    plot_statistics("p78", "p77", xmax=1500.0,
                    dy=0.4, bins=50, gain=2.4, robust=True)
    plot_statistics("p86", "p87", xmax=1500.0,
                    dy=0.4, bins=50, gain=2.4, robust=True)

    # 170-337 - doesn't work at all (as in: crashes)
    # plot_statistics("p71", "p72", xmax=3000.0,
    #                 dy=1.0, bins=50, gain=2.4, robust=True)


    #
    # The technique does not work very well with images from different
    # nights - not similar enough
    #
    
    # 177-341 - different nights - makes two strands
    plot_statistics("Keck2/p59", "p75", xmax=3000.0,
                    dy=1.0, bins=50, gain=2.4, robust=True)

    # plot_statistics("Keck2/p64", "p85", xmax=3000.0,
    #                 dy=1.0, bins=50, gain=2.4, robust=True)



# Old stuff
# labels = []
# for iorder in orders:
#     image = pyfits.open(flat_order_file(iorder))["SCI"].data
#     for ystart in ystarts:
#         ywindow = slice(ystart, ystart+ychunksize, None)
#         window = xwindow, ywindow
#         bright.append(image[window].mean())
#         sigma.append(image[window].std())
#         labels.append(str(iorder))
# plt.plot(bright, sigma, '.')
# x = np.arange(0.0, 8000)
# plt.plot(x, np.sqrt(x/GAIN))
# # for x, y, label in zip(bright, sigma, labels):
# #     plt.annotate(label, (x, y))
# plt.xlim(0.0, 8000.0)
# plt.ylim(0.0, 100.0)
# plt.savefig("sigma-vs-signal.pdf")








