"""
Tools for working in XYZ, UVW, and Δμ_δ, Δμ_α' [km/s] coordinates.

Contents:
    calculate_XYZ_given_RADECPLX
"""

###########
# imports #
###########
import numpy as np
from numpy import array as nparr
import astropy.units as u

from astropy.coordinates import Galactocentric
import astropy.coordinates as coord
_ = coord.galactocentric_frame_defaults.set('v4.0')

##########
# config #
##########

VERBOSE = 1
if VERBOSE:
    print(Galactocentric())

#############
# functions #
#############

def calculate_XYZ_given_RADECPLX(ra, dec, plx):
    """
    Given numpy arrays of right ascension, declination, and parallax, calculate
    the corresponding galactic XYZ coordinates.  NOTE: this code converts from
    parallax to distance assuming parallax [arcsec] = 1/distance[pc].  This
    assumption is wrong in the regime of low signal-to-noise parallaxes (very
    faint stars).

    Args:
        ra/dec/plx: np.ndarray of right ascension, declination, and parallax.
        RA/DEC in units of degrees.  Parallax in units of milliarcseconds.

    Returns:
        X, Y, Z: tuple of corresponding physical coordinates.  In this
        coordinate system the Sun is at {X,Y,Z}={-8122,0,+20.8} parsecs.
    """

    # convert from parallax to distance assuming high S/N parallaxes.
    distance = nparr(1/(plx*1e-3))*u.pc

    _c = coord.SkyCoord(
        ra=nparr(ra)*u.deg, dec=nparr(dec)*u.deg, distance=distance,
        pm_ra_cosdec=None, pm_dec=None, radial_velocity=None
    )

    # get XYZ in galactocentric
    c = _c.transform_to(coord.Galactocentric())

    x = c.x.to(u.pc).value
    y = c.y.to(u.pc).value
    z = c.z.to(u.pc).value

    return x,y,z
