import numpy as np
import lsst.sims.skybrightness as sb
from scipy.stats import binned_statistic

from lsst.sims.selfcal.analysis import healplots
from lsst.sims.utils import _altAzPaFromRaDec
from lsst.sims.maf.utils.telescopeInfo import TelescopeInfo
import healpy as hp
import os
import lsst.sims.photUtils.Bandpass as Bandpass
from lsst.sims.utils import haversine

def healpixels2dist(nside, hp1,hp2):
    """
    Compute the angular distance between 2 healpixels
    """
    lat1, ra1 = hp.pix2ang(nside,hp1)
    lat2, ra2 = hp.pix2ang(nside,hp2)
    dec1 = np.pi/2.0 - lat1
    dec2 = np.pi/2.0 - lat2

    angDist = haversine(ra1,dec1,ra2,dec2)
    return angDist




# Make some validation plots of the model compared to the observations we have

# Select all the near-zenith photometry from the Canon camera
#sqlQ = 'select stars.ra, stars.dec, obs.sky, obs.filter, dates.mjd  from obs,dates,stars where obs.starID = stars.ID and obs.dateID=dates.ID and obs.alt > 85;'
#skydata, mjd = sb.allSkyDB(2744 , sqlQ=sqlQ, filt=None,
#                           dtypes=zip(['ra','dec','sky','filt','mjd'],
#                                      [float,float,float,'|S1',float]))


# Set up the Canon fiters
canonDict = {}
canonFiles = {'R':'red_canon.csv','G':'green_canon.csv','B':'blue_canon.csv'}

path = os.path.join(os.environ.get('SIMS_SKYBRIGHTNESS_DATA_DIR'), 'Canon')
for key in canonFiles.keys():
    data = np.loadtxt(os.path.join(path,canonFiles[key]), delimiter=',',
                      dtype=zip(['wave','throughput'],[float,float]))
    band = Bandpass()
    band.setBandpass(data['wave'], data['throughput'])
    canonDict[key]=band




sqlQ = 'select id,mjd,sunAlt,moonAlt from dates'
dateData,mjd = sb.allSkyDB(2744 , sqlQ=sqlQ, filt=None,dtypes=zip(['dateID', 'mjd', 'moonAlt','sunAlt'],
                                                                  [int,float,float,float]))


skipsize = 10000
indices = np.arange(0,dateData.size, skipsize)

# Maybe bin things on 15-min timescales to cut down the number of calls I need to make to the model
#binsize = 15./60./24.
#edges = np.arange(skydata['mjd'].min(),skydata['mjd'].max()+binsize*2, binsize)

filters = ['R','G','B']
filters = ['R']

sm = sb.SkyModel(mags=False)
telescope = TelescopeInfo('LSST')

nside = 16
# Demand this many stars before trying to fit. This should reject the very cloudy frames
starLimit = 200

# XXX---OK, I'm getting NaNs this way.  Maybe a better way would be to just pull every Nth mjd
# that fits the sun alt limits, and then I can check the number of stars, demand there be over XXX,
# and then make a map, then I can compare the zenith model, and find the model and observation faintest directions.


# Create array to hold the results of
names = ['moonAlt', 'sunAlt', 'obsZenith', 'modelZenith', 'obsDarkestHP',
         'modelDarkestHP', 'obsPointing', 'modelPointing']
types = [float, float,float,float,int,int,float,float]
validationArr = np.zeros(indices.size, dtype=zip(names,types) )


for filterName in filters:
    for i,indx in enumerate(indices):
        skydata,mjd = sb.allSkyDB(dateData[indx]['dateID'], filt=filterName)
        if skydata.size > starLimit:
            airmass = 1./np.cos(np.radians(90.-skydata['alt']))
            skydata = skydata[np.where(airmass < 2.1)]
            alt,az,pa =  _altAzPaFromRaDec(np.radians(skydata['ra']), np.radians(skydata['dec']),
                                           telescope.lon, telescope.lat, mjd)
            skyhp = healplots.healbin(az,alt, skydata['sky'], nside=nside)
            skyhp[np.isnan(skyhp)] = hp.UNSEEN

            sm.setRaDecMjd(np.radians(skydata['ra']), np.radians(skydata['dec']), mjd, degrees=False)
            sm.computeSpec()
            mags = sm.computeMags(canonDict[filterName])
            good = np.where(mags > 0)
            modelhp = healplots.healbin(az[good], alt[good], mags[good], nside=nside)

            modelhp[np.isnan(modelhp)] = hp.UNSEEN

            # ok, now to find the darkest spots and compute angular distance.
            # save data and model zenith sky brightnesses
            # save at some other alt/az combo?
            #validationArr['modelZenith'][i] = modelhp[0]
            #validationArr['obsZenith'][i] = skyhp[0]

            notnan = np.where(skyhp != hp.UNSEEN)
            validationArr['obsDarkestHP'][i] = np.where(skyhp == skyhp[notnan].max() )[0].min()
            notnan =np.where(modelhp != hp.UNSEEN)
            validationArr['modelDarkestHP'][i]  = np.where(modelhp == modelhp[notnan].max() )[0]

            good = np.where(skyhp[0:4] != hp.UNSEEN)[0]
            if good.size >= 1:
                validationArr['obsZenith'][i] = skyhp[0:4][good].mean()
            good = np.where(modelhp[0:4] != hp.UNSEEN)[0]
            if good.size >= 1:
                validationArr['modelZenith'][i] = modelhp[0:4][good].mean()
            validationArr['moonAlt'][i] = np.degrees(sm.moonAlt)
            validationArr['sunAlt'][i] = np.degrees(sm.sunAlt)

angDistances = healpixels2dist(nside,validationArr['obsDarkestHP'], validationArr['modelDarkestHP'])
