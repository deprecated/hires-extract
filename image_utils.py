import numpy as np
import os.path
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import pywcsgrid2
import matplotlib.pyplot as plt

KMS = 1000.0  # m/s


def kms2xpix(hdr, u, y=0.0):
    u = np.atleast_1d(u)
    coords = np.vstack((u*KMS, [y]*len(u))).T
    return pywcs.WCS(hdr).wcs_sky2pix(coords, 0)[:, 0]


def arcsec2ypix(hdr, y, u=0.0):
    y = np.atleast_1d(y)
    coords = np.vstack(([u]*len(y), y)).T
    return pywcs.WCS(hdr).wcs_sky2pix(coords, 0)[:, 1]


def doublet_stages(sname, dname, stampdir="Stamps", drange=None):
    """
    Plot all the stages in the reduction of a doublet
    """
    stampname = os.path.join(stampdir,
                             "-".join([sname, dname, "stamp"]))
    fig = plt.figure(figsize=[16, 7])
    h = None
    for suff, index, win, title in [
            ["", 0, 241, "Original"],
            ["-nc", 0, 242, "Line"],
            ["-nc", "CONT", 246, "Continuum"],
            ["-nc-dd", 0, 243, "Deconvolved"],
            ["-nc-dd", "PROP", 244, "Proplyd"],
            ["-nc-dd", "NEB", 248, "Nebula"]]:
        hdu = pyfits.open(stampname + suff + ".fits")[index]
        if index != "CONT":
            h = hdu.header
        d = hdu.data
        ax = pywcsgrid2.subplot(win, header=h)
        ax.imshow(d, vmin=drange[0], vmax=drange[1],
                  origin="low", cmap=plt.cm.gray_r,
                  interpolation="nearest")
        ax.set_ticklabel2_type("arcsec", scale=1./3600, nbins=8)
        ax.set_ticklabel1_type("absval", scale=1./1000, nbins=7)
        ax.set_xlim(*kms2xpix(h, [-30.0, 60.0]))
        ax.set_title(title)
        ax.grid()
        ax.axis["bottom"].label.set_text(
            r"$V_{\mathrm{top}}\ \mathrm{(km\,s^{-1})}$")

    plt.tight_layout()
    fig.savefig(stampname + "-doublet-stages.pdf")
    plt.close(fig)


def line_stages(sname, linename, fancyname=None,
                stampdir="Stamps", drange=None, sky=0):
    """
    Plot all the stages in the reduction of a normal line
    """
    stampname = os.path.join(stampdir,
                             "-".join([sname, linename, "stamp"]))
    fig = plt.figure(figsize=[12, 7])
    h = None
    for suff, index, win, title in [
            ["", 0, 231, "Original"],
            ["-nc", 0, 232, "Line"],
            ["-nc", "CONT", 235, "Continuum"],
            ["-nc", "PROP", 233, "Proplyd"],
            ["-nc", "NEB", 236, ["", "Sky + "][sky] + "Nebula"]]:
        hdu = pyfits.open(stampname + suff + ".fits")[index]
        if index != "CONT":
            h = hdu.header
        d = hdu.data
        ax = pywcsgrid2.subplot(win, header=h)
        im = ax.imshow(d, vmin=drange[0], vmax=drange[1],
                       origin="low", cmap=plt.cm.gray_r,
                       aspect="auto",
                       interpolation="nearest")
        ax.set_ticklabel2_type("arcsec", scale=1./3600, nbins=8)
        ax.set_ticklabel1_type("absval", scale=1./1000, nbins=7)
        ax.set_xlim(*kms2xpix(h, [-30.0, 60.0]))
        ax.set_title(title)
        ax.grid()
        ax.axis["bottom"].label.set_text(
            r"$V_{\mathrm{top}}\ \mathrm{(km\,s^{-1})}$")

    # Make a separate axis on the second row just for the line label
    # Use 5 times as many horizontal divisions
    ax = plt.subplot(2, 30, 31, frameon=False)
    ax.axis("off")
    ax.text(0.0, 0.5, fancyname or linename,
            fontsize="x-large",
            horizontalalignment='left',
            verticalalignment='center',
            transform=ax.transAxes)

    # And another one for the color bar
    cb_ax = plt.subplot(2, 30, 36)
    # colorbar_axis = plt.axes([0.2, 0.05, 0.02, 0.35])
    cb = plt.colorbar(im, cax=cb_ax, orientation="vertical")
    cb.set_label("Intensity")

    fig.subplots_adjust(bottom=0.08, top=0.95,
                        left=0.05, right=0.99,
                        wspace=0.25, hspace=0.4)
    fig.savefig(stampname + "-line-stages.pdf")
    plt.close(fig)


if __name__ == "__main__":
    # Example plots
    # doublet_stages("p84", "O_I_6046", drange=[-3.0, 180.0])
    # doublet_stages("p84", "O_I_7002", drange=[-3.0, 180.0])
    line_stages("p84", "N_I_5198", "[N I] 5198", drange=[-2.0, 140.0], sky=1)
    line_stages("p84", "O_I_6300", "[O I] 6300", drange=[-10.0, 1600.0], sky=1)
    line_stages("p84", "C_II_6578", "C II 6578", drange=[-2.0, 150.0])
    line_stages("p84", "He_I_S_6678", "He I 6678", drange=[-5.0, 700.0])
    line_stages("p84", "O_I_5577", "[O I] 5577", drange=[-3.0, 250.0], sky=1)
    line_stages("p84", "Si_II_6371", "Si II 6371", drange=[-3.0, 120.0])
    line_stages("p84", "O_III_5007", "[O III] 5007", drange=[-10.0, 10000.0])
    line_stages("p84", "S_II_6731", "[S II] 6731", drange=[-10.0, 5000.0])
    line_stages("p84", "N_II_6548", "[N II] 6548", drange=[-10.0, 10000.0])
    line_stages("p84", "Cl_III_5538", "[Cl III] 5538", drange=[-2.0, 100.0])
    line_stages("p84", "Fe_III_4881", "[Fe III] 4881", drange=[-2.0, 60.0])
    line_stages("p84", "Fe_II_5262", "[Fe II] 5262", drange=[-2.0, 25.0])
    line_stages("p84", "S_III_6312", "[S III] 6312", drange=[-4.0, 330.0])
