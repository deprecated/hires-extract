import numpy as np
import os.path
import scipy.ndimage as ni
import astropy.io.fits as pyfits
import astropy.wcs as pywcs
import pywcsgrid2
import matplotlib.pyplot as plt

KMS = 1000.0  # m/s
NPCONT = 8  # Number of proplyd contours
FWHM = 2.0*np.sqrt(2.0*np.log(2.0))  # Conversion RMS sigma -> FWHM


def kms2xpix(hdr, u, y=0.0):
    u = np.atleast_1d(u)
    coords = np.vstack((u*KMS, [y]*len(u))).T
    return pywcs.WCS(hdr).wcs_sky2pix(coords, 0)[:, 0]


def arcsec2ypix(hdr, y, u=0.0):
    y = np.atleast_1d(y)
    coords = np.vstack(([u]*len(y), y)).T
    return pywcs.WCS(hdr).wcs_sky2pix(coords, 0)[:, 1]


def line_stages(sname, linename, fancyname=None,
                stampdir="Stamps", drange=None, sky=0,
                pmax=None, psmooth=None, doublet=False):
    """
    Plot all the stages in the reduction of a normal or doublet line
    """
    SINGLET_TABLE = [
        ["", 0, 231, "Original"],
        ["-nc", 0, 232, "Line"],
        ["-nc", "CONT", 235, "Continuum"],
        ["-nc", "PROP", 233, "Proplyd"],
        ["-nc", "NEB", 236, ["", "Sky + "][sky] + "Nebula"]]
    DOUBLET_TABLE = [
        ["", 0, 241, "Original"],
        ["-nc", 0, 242, "Line"],
        ["-nc", "CONT", 246, "Continuum"],
        ["-nc-dd", 0, 243, "Deconvolved"],
        ["-nc-dd", "PROP", 244, "Proplyd"],
        ["-nc-dd", "NEB", 248, "Nebula"]]

    if doublet:
        table = DOUBLET_TABLE
        ncolumns = 4
        outid = "-doublet"
    else:
        table = SINGLET_TABLE
        ncolumns = 3
        outid = "-line"

    stampname = os.path.join(stampdir,
                             "-".join([sname, linename, "stamp"]))
    fig = plt.figure(figsize=[4*ncolumns, 7])
    h = None

    for suff, index, win, title in table:
        hdu = pyfits.open(stampname + suff + ".fits")[index]
        if index != "CONT":
            h = hdu.header
        d = hdu.data
        ax = pywcsgrid2.subplot(win, header=h)
        im = ax.imshow(d, vmin=drange[0], vmax=drange[1],
                       origin="low", cmap=plt.cm.gray_r,
                       aspect="auto",
                       interpolation="nearest")
        if index == "PROP" and pmax is not None:
            if psmooth is None:
                pd = d
            else:
                pd = ni.gaussian_filter(d.astype("d"), psmooth/FWHM)
            cs = ax.contour(pd, pmax*2**(-0.5*np.arange(NPCONT)),
                            colors="r", alpha=0.5)
            plt.clabel(cs, cs.levels[::2],
                       inline=True, fmt="%.0d", inline_spacing=1, fontsize=6)
        ax.set_ticklabel2_type("arcsec", scale=1./3600, nbins=8)
        ax.set_ticklabel1_type("absval", scale=1./1000, nbins=7)
        ax.set_xlim(*kms2xpix(h, [-30.0, 60.0]))
        ax.set_title(title)
        ax.grid()
        ax.axis["bottom"].label.set_text(
            r"$V_{\mathrm{top}}\ \mathrm{(km\,s^{-1})}$")

    # Make a separate axis on the second row just for the line label
    # Use 10 times as many horizontal divisions
    ax = plt.subplot(2, 10*ncolumns, 10*ncolumns+1, frameon=False)
    ax.axis("off")
    ax.text(0.0, 0.5, fancyname or linename,
            fontsize="x-large",
            horizontalalignment='left',
            verticalalignment='center',
            transform=ax.transAxes)

    # And another one for the color bar
    cb_ax = plt.subplot(2, 10*ncolumns, 10*ncolumns+6)
    cb = plt.colorbar(im, cax=cb_ax, orientation="vertical")
    cb.set_label("Intensity")

    fig.subplots_adjust(bottom=0.08, top=0.95,
                        left=0.05, right=0.99,
                        wspace=0.25, hspace=0.4)
    fig.savefig(stampname + outid + "-stages.pdf")
    plt.close(fig)


if __name__ == "__main__":
    # Example plots
    if True:
        line_stages("p84", "Fe_II_5159", "[Fe II] 5159",
                    drange=[-2.0, 60.0], pmax=40.0, psmooth=2.0)
        line_stages("p84", "He_I_T_5876", "He I 5876",
                    drange=[-5.0, 1000.0], pmax=600.0)
        line_stages("p84", "N_II_5755", "[N II] 5755",
                    drange=[-10.0, 800.0], pmax=600.0)
        line_stages("p84", "N_II_6583", "[N II] 6583",
                    drange=[-10.0, 60000.0], pmax=32000.0)
        line_stages("p84", "N_II_6548", "[N II] 6548",
                    drange=[-10.0, 15000.0], pmax=8000.0)
        line_stages("p84", "N_I_5200", "[N I] 5200",
                    drange=[-2.0, 100.0], pmax=100.0, sky=1)

    if False:
        # The reserve list
        line_stages("p84", "O_I_6046", "O I 6046", doublet=True,
                    drange=[-3.0, 180.0], pmax=400.0)
        line_stages("p84", "O_I_7002", "O I 7002", doublet=True,
                    drange=[-3.0, 180.0], pmax=400.0)
        line_stages("p84", "Si_II_6347", "Si II 6347",
                    drange=[-3.0, 160.0], pmax=120.0, psmooth=1.0)
        line_stages("p84", "Si_II_6371", "Si II 6371",
                    drange=[-3.0, 120.0], pmax=80.0, psmooth=1.0)
        line_stages("p84", "N_I_5198", "[N I] 5198",
                    drange=[-2.0, 140.0], pmax=160.0, sky=1)
        line_stages("p84", "O_I_6300", "[O I] 6300",
                    drange=[-10.0, 1600.0], sky=1, pmax=1600.0)
        line_stages("p84", "C_II_6578", "C II 6578",
                    drange=[-2.0, 150.0], pmax=80.0, psmooth=2.0)
        line_stages("p84", "He_I_S_6678", "He I 6678",
                    drange=[-5.0, 700.0], pmax=400.0)
        line_stages("p84", "O_I_5577", "[O I] 5577",
                    drange=[-3.0, 250.0], sky=1, pmax=160.0)
        line_stages("p84", "O_III_5007", "[O III] 5007",
                    drange=[-10.0, 10000.0], pmax=4000.0)
        line_stages("p84", "S_II_6716", "[S II] 6716",
                    drange=[-10.0, 3000.0], pmax=1600.0)
        line_stages("p84", "S_II_6731", "[S II] 6731",
                    drange=[-10.0, 5000.0], pmax=3200.0)
        line_stages("p84", "S_III_6312", "[S III] 6312",
                    drange=[-4.0, 330.0], pmax=320.0)
        line_stages("p84", "Cl_III_5518", "[Cl III] 5518",
                    drange=[-2.0, 50.0], pmax=40.0, psmooth=2.0)
        line_stages("p84", "Cl_III_5538", "[Cl III] 5538",
                    drange=[-2.0, 100.0], pmax=40.0, psmooth=2.0)
        line_stages("p84", "Fe_III_5270", "[Fe III] 5270",
                    drange=[-2.0, 40.0], pmax=20.0, psmooth=2.0)
        line_stages("p84", "Fe_III_4881", "[Fe III] 4881",
                    drange=[-2.0, 60.0], pmax=40.0, psmooth=2.0)
        line_stages("p84", "Fe_II_5262", "[Fe II] 5262",
                    drange=[-2.0, 25.0], pmax=40.0, psmooth=2.0)
