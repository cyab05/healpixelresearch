import sys
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import hstack
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy import optimize
from scipy import interpolate
from scipy.interpolate import griddata
from scipy.stats import gaussian_kde
from dl import authClient as ac, queryClient as qc
from dl.helpers.utils import convert
from getpass import getpass
import pandas as pd
import healpy as hp

class list_healpixel:
    def __init__(self, nside=32, min_ra=0, max_ra=360, min_dec=-90, max_dec=90):
        self.nside = nside
        self.min_ra = min_ra
        self.max_ra = max_ra
        self.min_dec = min_dec
        self.max_dec = max_dec
        
        self.nside_range = np.arange(0, hp.nside2npix(self.nside))
        self.ra, self.dec = self.pix2radec(self.nside, self.nside_range)
        
        self.pixelsinrange = self.check(self.nside_range, self.ra, self.dec, 
                                        self.min_ra, self.max_ra, self.min_dec, self.max_dec, self.nside)
    
    def pix2radec(self, nside, pixel):
        theta, phi = hp.pix2ang(nside, pixel)
        ra = np.degrees(phi)
        dec = 90 - np.degrees(theta)
        return ra, dec

    def radec2healpix(self, nside, ra, dec):
        theta = np.radians(90 - dec)
        phi = np.radians(ra)
        pixel = hp.ang2pix(nside, theta, phi)
        return pixel

    def check(self, indices, ra, dec, min_ra, max_ra, min_dec, max_dec, nside):
        healpix_pixels = []
        for index in indices:
            if min_ra < ra[index] < max_ra and min_dec < dec[index] < max_dec:
                pixel = self.radec2healpix(nside, ra[index], dec[index])
                healpix_pixels.append(pixel)
        return healpix_pixels

# Example usage
#list_healpixels = ListHealpixels(nside=32, min_ra=0, max_ra=360, min_dec=-90, max_dec=90)
#print(list_healpixels.pixelsinrange)
