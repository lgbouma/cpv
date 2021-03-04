"""
Contents:
    plot_river
"""
import os, corner, pickle
from glob import glob
from datetime import datetime
import numpy as np, matplotlib.pyplot as plt, pandas as pd, pymc3 as pm
from numpy import array as nparr

import matplotlib as mpl

from astropy import units as u, constants as const
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.table import Table

import matplotlib.patheffects as pe
from matplotlib.ticker import MaxNLocator

from aesthetic.plot import savefig, format_ax, set_style

def plot_river(time, flux, period, t0=None):
    pass
