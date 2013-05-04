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

GAIN = 2.09

image1 = flat_image("q69")
image2 = flat_image("q110")

ratio = image2/image1

mean, median = ratio.mean(), np.median(ratio)
print "Average image ratio: mean, median = ", mean, median

every = 1
x = image2.ravel()
y = ratio.ravel()
plt.plot(x[::every], y[::every], '.', alpha=0.003)
plt.plot([0.0, 25000.0], [median, median])
xmin, xmax = 0.0, 25000.0
ymin, ymax = median*(1-0.03), median*(1+0.03)
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.savefig("measure-noise.png")

H, xedges, yedges = np.histogram2d(
    x, y, bins=50, range=[[xmin, xmax], [ymin, ymax]], normed=True)
plt.clf()
plt.imshow(H.T, extent=[xmin, xmax, ymin, ymax],
           cmap=plt.cm.gray_r,
           origin="low", aspect="auto", interpolation="nearest")
plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
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








