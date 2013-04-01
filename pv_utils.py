"""
General utility routines for dealing with position-velocity (PV) images
"""

import numpy as np
from astropy import coordinates as coord
from astropy import units as u
from pyslalib.slalib import sla_dcs2c, sla_evp, sla_rverot, sla_obs


def make_grids(hdr):
    """Create 2D velocity and position grids from WCS info in FITS header

    """
    nx, ny = [hdr["NAXIS{}".format(i)] for i in 1, 2]
    i0, u0, dudi, dudj = [hdr[k] for k in
                          "CRPIX1", "CRVAL1", "CD1_1", "CD1_2"]
    j0, y0, dydi, dydj = [hdr[k] for k in
                          "CRPIX2", "CRVAL2", "CD2_1", "CD2_2"]
    ii, jj = np.meshgrid(np.arange(nx) + 1, np.arange(ny) + 1)
    u = u0 + dudi*(ii - i0) + dudj*(jj - j0)
    y = y0 + dydi*(ii - i0) + dydj*(jj - j0)
    return u, y


def rvcorr(hdr):
    """Calculate helio- and geo-centric corrections to radial velocities

    The result should be subtracted from the observed (topocentric)
    velocities in order to give velocities in the heliocentric frame.

    """
    # The following page was very helpful in pointing to the relevant
    # SLALIB routines and how to use them
    # http://star-www.rl.ac.uk/docs/sun67.htx/node230.html

    # Object coordinates
    ra = coord.RA(hdr["RA"], u.hour)
    dec = coord.Dec(hdr["DEC"], u.degree)
    # Modified Julian Date
    jdate = float(hdr["MJD"])
    # Sidereal time
    st = coord.Angle((hdr["ST"]), u.hour)
    
    # line-of-sight unit vector to astronomical object
    k_los = sla_dcs2c(ra.radians, dec.radians)
    
    # Velocity and position of earth in barycentric and heliocentric frames
    # Units are AU and AU/s
    vel_bary, pos_bary, vel_hel, pos_hel = sla_evp(jdate, 2000.0)
    
    # Radial velocity correction (km/s) due to helio-geocentric transformation  
    # Positive when earth is moving away from object
    vcorr_hel = u.AU.to(u.km, -np.dot(vel_hel, k_los))
    
    # Long/lat/altitude of observatory (radians, radians, meters)
    obs_id, obs_name, obs_long, obs_lat, obs_height = sla_obs(0, "KECK1")
    
    # Radial velocity correction (km/s) due to geo-topocentric transformation
    # Positive when observatory is moving away from object
    vcorr_geo = sla_rverot(obs_lat, ra.radians, dec.radians, st.radians)

    return vcorr_hel + vcorr_geo


if __name__ == "__main__":
    print "TODO: write tests and demo code for this module"
