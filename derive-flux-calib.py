import numpy as np
from numpy.polynomial import Chebyshev as T
import astropy.io.fits as pyfits
import lmfit
import yaml

ORDERMIN, ORDERMAX = 51, 76
# HST 1 slit 
OBJ = "p75f"
# Power law index for continuum derived from STIS spectrum
M = 0.07            
# Normalization at 7000 A
F7000 = 28.0
# Full extent of the slit
J1 = 15
J2 = 47
# How wide a range to mask out around each line
DWAV = 5.0

           

def main():
    orders = range(ORDERMIN, ORDERMAX)
    linedb = yaml.load(open("Stamps/line-database.json"))
    for i in orders:
        f = pyfits.open("Extract/{}-order{}.fits".format(OBJ, i))
        spec = f["SCI"].data[J1:J2,:].mean(axis=0) 
        wav = f["WAV"].data[J1:J2,:].mean(axis=0)
        stis = F7000*(wav/7000.0)**M
        ratio = spec/stis
        # Mask out NaNs
        mask = np.isfinite(wav) & np.isfinite(ratio)
        # Mask out lines
        linewavs = [v['wav'] for v in linedb.values()
                    if v['iorder'] == i]
        for wav0 in linewavs:
            inthisline = np.abs(wav-wav0) < 0.5*DWAV
            mask[inthisline] = False
        # And just to be on the safe side, mask out outliers too
        outliers = np.abs(ratio - 1.0) > 0.75 
        mask[outliers] = False

        print i, len(mask), sum(mask)

        # Fit third order polynomial
        # coeffs = np.polyfit(wav, ratio, 3)
        # p = np.poly1d(coeffs)

        # print i, coeffs

        


if __name__ == "__main__":
    main()
