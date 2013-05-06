"""
Measure noise in the flat field images, and plot noise against signal level
"""
import astropy.io.fits as pyfits
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

def flat_image(specid="q69", folder="Keck1"):
    """Return bias-subtracted original image
    """
    hdu, = pyfits.open(
        os.path.join(folder, specid + ".fits")
    )
    bzero = hdu.header.get("BZERO", 0.0)
    bscale = hdu.header.get("BSCALE", 1.0)
    print bzero, bscale
    return (hdu.data - bzero)/bscale

# GAIN = 2.38
GAIN = 4.5

image1 = flat_image("q69")
image2 = flat_image("q110")

ratio = image2/image1

mean, median = ratio.mean(), np.median(ratio)
print "Average image ratio: mean, median = ", mean, median

every = 1
x = image2.ravel()
y = ratio.ravel()/median
xmin, xmax = 0.0, 21500.0
ymin, ymax = 1 - 0.05, 1 + 0.05
plt.plot(x[::every], y[::every], '.', alpha=0.003)
plt.plot([xmin, xmax], [1.0, 1.0])
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.savefig("measure-noise.png")

H, xedges, yedges = np.histogram2d(
    x, y, bins=100, range=[[xmin, xmax], [ymin, ymax]], normed=True)
# Normalize to the sum for each x-value
H /= H.sum(axis=1)[:, None]

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

Hsig_theory = np.sqrt((1.0 + median)/(xgrid*GAIN))

plt.clf()
plt.imshow(H.T, extent=[xmin, xmax, ymin, ymax],
           cmap=plt.cm.gray_r,
           origin="low", aspect="auto", interpolation="nearest")
plt.plot(xgrid, Hmean, label="Mean")
plt.plot(xgrid, Hmean + Hsig_theory, '--b', label="Theoretical sigma")
plt.plot(xgrid, Hmean - Hsig_theory, '--b')
plt.plot(xgrid, Hmean + Hsig_obs, '--r', label="Observed sigma")
plt.plot(xgrid, Hmean - Hsig_obs, '--r')
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.xlabel("Brightness of Flat 1")
plt.ylabel("Ratio (Flat 1/Flat 2) / <(Flat 1/Flat 2)>")
plt.legend()
plt.savefig("measure-noise2.png")

sys.exit()


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








