import numpy as np
import ephem
from lsst.sims.utils import calcLmstLast

def wrapRA(ra):
    """
    Wrap only RA values into 0-2pi (using mod).
    """
    ra = ra % (2.0*np.pi)
    return ra

def raDecToAltAz(raRad, decRad, longRad, latRad, mjd):
    """
    Numpy-ified the version from  lsst.sims.catalogs.generation.db
    """
    lst,last = calcLmstLast(mjd, longRad)
    haRad = np.radians(last*15.) - raRad
    altRad = np.arcsin(np.sin(decRad)*np.sin(latRad)+np.cos(decRad)*np.cos(latRad)*np.cos(haRad))
    argument = (np.sin(decRad) - np.sin(altRad)*np.sin(latRad))/(np.cos(altRad)*np.cos(latRad))
    # Make sure we don't screw up on machine precision and throw nans.
    argument[np.where(argument > 1.)] = 1.
    argument[np.where(argument < -1.)] = -1.
    azRad = np.arccos(argument)
    azRad = wrapRA(azRad)
    rollOver = np.where(np.sin(haRad) >= 0)
    azRad[rollOver] = 2.*np.pi - azRad[rollOver]

    return altRad, azRad


def altAzToRaDec(altRad, azRad, longRad, latRad, mjd):
    """
    Numpy-ified the version from  lsst.sims.utils.coordinateTransformations
    """
    lst,last = calcLmstLast(mjd, longRad)
    decRad = np.arcsin(np.sin(latRad)*np.sin(altRad)+ np.cos(latRad)*np.cos(altRad)*np.cos(azRad))
    haRad = np.arccos((np.sin(altRad) - np.sin(decRad)*np.sin(latRad))/(np.cos(decRad)*np.cos(latRad)))
    haRad[np.isnan(haRad)] = 0. #catch nans from machine precision roundoff issues going between deg and rad
    raRad = np.radians(last*15.) - haRad
    return raRad, decRad


def mjd2djd(inDate):
    """
    Convert Modified Julian Date to Dublin Julian Date (what pyephem uses).
    """
    if not hasattr(mjd2djd, 'doff'):
        mjd2djd.doff = ephem.Date(0)-ephem.Date('1858/11/17')
    djd = inDate-mjd2djd.doff
    return djd

def robustRMS(array, missing=0.):
    """
    Use the interquartile range to compute a robust approximation of the RMS.
    if passed an array smaller than 2 elements, return missing value
    """
    if np.size(array) < 2:
        rms = missing
    else:
        iqr = np.percentile(array,75)-np.percentile(array,25)
        rms = iqr/1.349 #approximation
    return rms

def ut2Mjd(dateString):
    obs = ephem.Observer()
    obs.date = dateString
    doff = ephem.Date(0)-ephem.Date('1858/11/17')
    mjd = obs.date+doff
    return mjd

def mjd2ut(mjd):
    obs = ephem.Observer()
    doff = ephem.Date(0)-ephem.Date('1858/11/17')
    djd = mjd-doff
    obs.date = djd
    print(obs.date)
