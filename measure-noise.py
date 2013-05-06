"""
Measure noise in the flat field images, and plot noise against signal level
"""
import astropy.io.fits as pyfits
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

# GAIN = 2.38


def flat_image(specid="q69", folder="Keck1", dark=5.5):
    """Return bias-subtracted original image
    """
    hdu, = pyfits.open(
        os.path.join(folder, specid + ".fits")
    )
    bzero = hdu.header.get("BZERO", 0.0)
    bscale = hdu.header.get("BSCALE", 1.0)
    print bzero, bscale
    return (hdu.data - bzero)/bscale - dark


def plot_statistics(file1="q69", file2="q110",
                    gain=4.4, dy=0.05, xmax=None, bins=200):
    """Plot noise statistics from ratio of two images
    
    The two images should be as close to identical as possible (e.g.,
    two quartz lamp exposures, or two science exposures of the same
    object at same PA and position), but they may have different
    exposure times.

    """
    image1 = flat_image(file1)
    image2 = flat_image(file2)

    ratio = image2/image1

    mean, median = ratio.mean(), np.median(ratio)
    print "Average image ratio: mean, median = ", mean, median

    R0 = median                 # Guess for ratio of exposure times

    every = 1
    x = image2.ravel()
    y = ratio.ravel()/R0

    xmin = 0.0
    if xmax is None:
        xmax = image2.max()
    ymin, ymax = 1 - dy, 1 + dy
    plt.clf()
    plt.plot(x[::every], y[::every], '.', alpha=0.003)
    plt.plot([xmin, xmax], [1.0, 1.0])
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.savefig("measure-noise.png")

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

    # H is already normalized so that integral along y-axis is unity
    Hmean = np.sum(Ygrid*H, axis=1)
    dY = Ygrid - Hmean[:, None]
    Hsig_obs = np.sqrt(np.sum((dY**2)*H, axis=1))

    # Assuming we know the gain, then we can predict the std of the ratio
    # image as a function of brightness

    Hsig_theory = np.sqrt((1.0 + R0)/(xgrid*gain))

    plt.clf()
    np.set_printoptions(precision=4)
    print H.T[::20, ::20]
    print xmin, xmax, ymin, ymax
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



if __name__ == "__main__":
    plot_statistics("q69", "q110", xmax=21500.0, dy=0.5)
    plot_statistics("p83", "p84", xmax=3000.0,
                    dy=0.5, bins=50, gain=2.4)



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








