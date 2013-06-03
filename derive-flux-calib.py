import numpy as np
from numpy.polynomial import Chebyshev as T
from scipy.stats import norm
import astropy.io.fits as pyfits
import lmfit
import yaml
import argparse
import matplotlib.pyplot as plt

ORDERMIN, ORDERMAX = 51, 76
# HST 1 slit 
OBJ = "p75f"
OBJ2 = "p75c"
# Power law index for continuum derived from STIS spectrum
M = 0.07            
# Normalization at 7000 A
F7000 = 28.0
# Full extent of the slit
J1 = 15
J2 = 47
# How wide a range to mask out around each line
DWAV = 3.0
DWAV_STRONG = 15.0

STRONG_LINES = [
    "H a 6563", "H b 4861", "[N II] 6583", "[N II] 6548",
    "[O III] 4959", "[O III] 5007", # "[S III] 6312", "[O I] 6300", 
    "[S II] 6731", "[S II] 6716", "He I S 6678", "He I T 5876",
]


def model(x, params, npoly=2):
    """Gaussian (centered on x=0) plus Chebychev polynomial

    """
    gsigma = params["gsigma"].value
    gheight = params["gheight"].value
    Tcoeffs = [params["T{}".format(k)].value for k in range(npoly + 1)]
    return gheight*norm.pdf(x, loc=0.0, scale=gsigma) + T(Tcoeffs, domain=[0.0, 1.0])(x)


def model_minus_data(params, x, data, npoly=2):
    """Function to minimize"""
    return model(x, params, npoly) - data

    
def main(npoly):
    orders = range(ORDERMIN, ORDERMAX)
    linedb = yaml.load(open("Stamps/line-database.json"))
    plt.plot([4700.0, 7020.0], [1.0, 1.0], "k", lw=5, alpha=0.3)
    for i in orders:
        f = pyfits.open("Extract/{}-order{}.fits".format(OBJ, i))
        spec = f["SCI"].data[J1:J2,:].mean(axis=0) 
        wav = f["WAV"].data[J1:J2,:].mean(axis=0)
        stis = F7000*(wav/7000.0)**M
        ratio = spec/stis
        # Mask out NaNs
        mask = np.isfinite(wav) & np.isfinite(ratio)
        wavmin, wavmax = np.min(wav[mask]), np.max(wav[mask])
        # Mask out lines
        linewavs = [v['wav'] for k, v in linedb.items()
                    if v['iorder'] == i]
        strongwavs = [v['wav'] for k, v in linedb.items()
                      if v['iorder'] == i and k in STRONG_LINES]
        for wav0 in linewavs:
            inthisline = np.abs(wav-wav0) < 0.5*DWAV
            mask[inthisline] = False
        for wav0 in strongwavs:
            inthisline = np.abs(wav-wav0) < 0.5*DWAV_STRONG
            mask[inthisline] = False
        # And just to be on the safe side, mask out outliers too
        outliers = np.abs(ratio - 1.0) > 0.75 
        mask[outliers] = False

        print i, len(mask), sum(mask)

        # Normalize the wavelength scale
        x = (wav - wavmin)/(wavmax - wavmin)

        # Fit model with lmfit
        params = lmfit.Parameters()
        params.add("gheight", value=0.5, min=0.0, max=1.0)
        params.add("gsigma", value=0.20, vary=False, min=0.01, max=0.5)
        params.add("T0", value=1.0, min=0.3, max=1.5)
        for k in range(npoly):
            params.add("T"+str(k+1), value=0.0, min=-0.3, max=0.3)
        try:
            result = lmfit.minimize(model_minus_data, params,
                                    args=(x[mask], ratio[mask], npoly))
            print "Reduced chi-squared of fit: ", result.redchi
            lmfit.report_errors(params)
        except:
            print "Fit failed"
            continue
        # # Fit third order polynomial
        # try: 
        #     coeffs = np.polyfit(x[mask], ratio[mask], npoly)
        #     p = np.poly1d(coeffs)
        #     print coeffs
        # except:
        #     print "polyfit failed"
        #     continue
        
        plt.plot(wav, ratio, "-", alpha=0.7, lw=0.3)
        plt.plot(wav, model(x, params, npoly), "-k", alpha=0.4, lw=2)
        meanwav = wav[mask].mean()
        plt.text(meanwav, 0.2, "{}".format(i), horizontalalignment='center')

        # save corrected spectra
        correct = model(x, params, npoly)[None, :]
        f["SCI"].data /= correct
        f["SIGMA"].data /= correct
        f.writeto("Extract/{}-order{}.fits".format(OBJ2, i), clobber=True)

    F = plt.gcf()
    F.set_size_inches( (30, 3) )
    plt.grid()
    plt.xlim(4700., 7020)
    plt.ylim(0.0, 2.1)
    plt.xlabel("Wavelength, Angstrom")
    plt.ylabel("F(lam) / lam^{:.2f}".format(M))
    pltname = "derive-flux-calib.png"
    plt.savefig(pltname, bbox_inches='tight', dpi=300)



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--npoly", type=int, default=3, help="Order of polynomial")
    cmd_args = parser.parse_args()

    main(**vars(cmd_args))
