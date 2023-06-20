"""
1-d comet model

CITE:
* Zieba+2019 (actual portions of code directly copied!)
* Brogi+2012 (original 1d formalism)

Contents:
    | vcirc
    | P_to_a
    | rc_hat
    | disk_intensity
    | rho
    | impact_param
    | rchordit
"""

import numpy as np, matplotlib.pyplot as plt
from astropy import units as u
from astropy import constants as c
import scipy.optimize as op
import os

@u.quantity_input
def vcirc(m1:u.Msun,m2:u.Mjup,a:u.au)->u.km/u.s:
    """
    Circular orbital velocity of m2 about m1 at distance a

    Args:
        m1, m2: Primary and secondary masses
        a: semimajor axis

    Returns:
        velocity: circular orbital velocity

    >>> import astropy.units as u
    >>> vcirc(1.0 *u.M_sun, 1.0 * u.M_jup, 5.2 * u.au)
    <Quantity 13.06768412 km / s>
    """
    mu = c.G * (m1 + m2)
    vcirc = np.power(mu/a, 0.5)
    return vcirc

@u.quantity_input
def P_to_a(P:u.year, m1:u.M_sun, m2:u.M_jup)->u.au:
    """calculate orbital radius from period

    Args:
        P: orbital period
        m1, m2: Primary and secondary masses

    Returns:
        a: semi-major axis

    >>> import astropy.units as u
    >>> P_to_a(11.86*u.year, 1.0*u.M_sun, 1.0*u.M_jup)
    <Quantity 5.20222482 AU>
    """

    # a^3/P^2 = (G/4pipi) (m1 + m2)
    const = c.G / (4.*np.pi*np.pi)
    mu = m1 + m2
    a3 = P*P*const*mu
    aa = np.power(a3, 1./3.)

    return aa

def rc_hat(rc, a):
    'angular size of chord on star as seen from distance a'
    return np.arcsin(rc/(2*a))

def disk_intensity(thetad, _rc_hat, a, Rstar, ulimb=0.79):
    """
    stellar disk intensity collapsed to 1D
    (Zieba+19 Eq2, or more accurate citation would be Brogi+2012)
    """
    s1 = np.power(np.sin(_rc_hat/2),2) - np.power(np.sin(thetad),2)
    s1_negative = (s1<0)
    s1[s1_negative] = 0. # to prevent error with sqrt below
    s2 = 1 - (a/Rstar)*np.sqrt(s1)
    res = 1-(ulimb*s2)
    res[s1_negative] = 0
    return res

# make an exponential decay absorption for a comet tail
def rho(dtheta, ce, lam):
    'extinction cross section for comet'
    rho = ce * np.exp(-lam*dtheta)
    rho[(dtheta<0)] = 0 # no dust in front of the comet
    return rho

def impact_param(rc, Rstar):
    'impact parameter for a star and orbit around the star'
    d1 = rc/(2*Rstar)
    d2 = 1 - np.power(d1,2.)
    return np.sqrt(d2)

def rchordit(b, Rstar):
    'chord length for b'
    # b2 + (rc/2R)2 = 1
    # rc/2R = sqrt(1-b2)
    d1 = np.sqrt(1.-b*b)
    return d1*2*Rstar

def mu_comet_model(tin, tmid, b_impact, cmax, lam, P, Mstar, Rstar,
                   ulimb, extra=False):

    tin_min = tin.min()
    tin_max = tin.max()
    t_delt = tin_max - tin_min

    # NOTE: lower this number to speed up
    # TODO TODO TODO
    # TODO TODO TODO
    # make a wider range of regularly spaced times so that convolution doesn't have round off error
    t_model = np.linspace(tin_min-t_delt, tin_max+t_delt, 100001)

    # convert time to orbital phase angle
    t_phase = (2*np.pi* (t_model - tmid) / P) * u.radian

    # convert period to semimajor axis
    a = P_to_a(P, Mstar, 0.0 * c.M_sun)

    # calculate size of a chord across the star
    rc = rchordit(b_impact, Rstar)

    # calculate angular size of chord across the star
    _rc_hat = rc_hat(rc, a)

    # calculate intensity profile for a chord across the star
    I = disk_intensity(t_phase, _rc_hat, a, Rstar, ulimb)
    Inorm = I / I.sum()

    # calculate comet tail
    rh = rho(t_phase, cmax, lam)

    # convolve two curves together to get total intensity
    Itot = np.convolve(Inorm, rh, 'same')

    # now use linear interpolator to find values of Itot at the tin times
    Itot_tin = np.interp(tin, t_model, Itot)

    if extra:
        return t_phase, Inorm, rh, 1 - Itot_tin

    return 1 - Itot_tin
