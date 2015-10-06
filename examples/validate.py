import numpy as np
import lsst.sims.skybrightness as sb
from scipy.stats import binned_statistic_2d

from lsst.sims.selfcal.analysis import healplots
from lsst.sims.utils import _altAzPaFromRaDec
from lsst.sims.maf.utils.telescopeInfo import TelescopeInfo
import healpy as hp
import os
import lsst.sims.photUtils.Bandpass as Bandpass
from lsst.sims.utils import haversine
import matplotlib.pylab as plt

def robustRMS(data):
    iqr = np.percentile(data,75)-np.percentile(data,25)
    rms = iqr/1.349 #approximation
    return rms


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
dateData,mjd = sb.allSkyDB(2744 , sqlQ=sqlQ, filt=None,dtypes=zip(['dateID', 'mjd','sunAlt', 'moonAlt'],
                                                                  [int,float,float,float]))


skipsize = 10
indices = np.arange(0,dateData.size, skipsize)

# Maybe bin things on 15-min timescales to cut down the number of calls I need to make to the model
#binsize = 15./60./24.
#edges = np.arange(skydata['mjd'].min(),skydata['mjd'].max()+binsize*2, binsize)

filters = ['R','G','B']
sm = sb.SkyModel(mags=False)
telescope = TelescopeInfo('LSST')

nside = 16
# Demand this many stars before trying to fit. This should reject the very cloudy frames
starLimit = 200

# Create array to hold the results of
names = ['moonAlt', 'sunAlt', 'obsZenith', 'modelZenith', 'obsDarkestHP',
         'modelDarkestHP', 'obsBrightestHP',
         'modelBrightestHP', 'obsPointing', 'modelPointing', 'angDistancesFaint', 'angDistancesBright']

# I guess I might as well add the brightest points in the sky as well...

types = [float, float,float,float,int,int,int,int, float,float, float, float, float]
validationArr = np.zeros(indices.size, dtype=zip(names,types) )

# Don't look at Canon above this limit.
amMax = 2.1


npix = hp.nside2npix(nside)
hpid = np.arange(npix)
lat,lon = hp.pix2ang(nside,hpid)
hpra = lon
hpdec = np.pi/2.-lat



for filterName in filters:
    for i,indx in enumerate(indices):
        skydata,mjd = sb.allSkyDB(dateData[indx]['dateID'], filt=filterName)
        if skydata.size > starLimit:
            airmass = 1./np.cos(np.radians(90.-skydata['alt']))
            skydata = skydata[np.where((airmass < amMax) & (airmass >= 1.))]
            alt,az,pa =  _altAzPaFromRaDec(np.radians(skydata['ra']), np.radians(skydata['dec']),
                                           telescope.lon, telescope.lat, mjd)
            skyhp = healplots.healbin(az,alt, skydata['sky'], nside=nside)
            skyhp[np.isnan(skyhp)] = hp.UNSEEN

            good = np.where(skyhp != hp.UNSEEN)
            sm.setRaDecMjd(hpra[good], hpdec[good], mjd, degrees=False, azAlt=True)
            sm.computeSpec()
            modelhp = np.zeros(npix,dtype=float)+hp.UNSEEN
            modelhp[good] = sm.computeMags(canonDict[filterName])

            notnan = np.where(skyhp != hp.UNSEEN)
            validationArr['obsDarkestHP'][i] = np.where(skyhp == skyhp[notnan].max() )[0].min()
            validationArr['obsBrightestHP'][i] = np.where(skyhp == skyhp[notnan].min() )[0].min()
            notnan =np.where(modelhp != hp.UNSEEN)
            validationArr['modelDarkestHP'][i]  = np.where(modelhp == modelhp[notnan].max() )[0]
            validationArr['modelBrightestHP'][i]  = np.where(modelhp == modelhp[notnan].min() )[0]

            good = np.where(skyhp[0:4] != hp.UNSEEN)[0]
            if good.size >= 1:
                validationArr['obsZenith'][i] = skyhp[0:4][good].mean()
            good = np.where(modelhp[0:4] != hp.UNSEEN)[0]
            if good.size >= 1:
                validationArr['modelZenith'][i] = modelhp[0:4][good].mean()
            validationArr['moonAlt'][i] = np.degrees(sm.moonAlt)
            validationArr['sunAlt'][i] = np.degrees(sm.sunAlt)

        else:
            validationArr['moonAlt'][i] = -666

    validationArr['angDistancesFaint'] = np.degrees(healpixels2dist(nside,validationArr['obsDarkestHP'],
                                                                    validationArr['modelDarkestHP']))
    validationArr['angDistancesBright'] = np.degrees(healpixels2dist(nside,validationArr['obsBrightestHP'],
                                                                    validationArr['modelBrightestHP']))

    fig,ax = plt.subplots()

    good = np.where(validationArr['moonAlt'] != -666)
    moonAltBins = np.arange(-90,90+5, 5)
    sunAltBins = np.arange(-90,90+5, 5)

    #medianResid, x_edge,y_edge, bn = binned_statistic_2d(validationArr['moonAlt'][good],
    #                                                     validationArr['sunAlt'][good],
    #                                                     validationArr['angDistancesFaint'][good],
    #                                                     statistic=np.median, bins=20)


    #xx,yy = np.meshgrid(xx,yy)
    #good = ~np.isnan(medianResid)
    #im = ax.imshow(medianResid, extent=[x_edge.min(),x_edge.max(),y_edge.min(),y_edge.max()],
    #               origin='lower', interpolation='nearest')
    im = ax.hexbin(validationArr['moonAlt'][good],validationArr['sunAlt'][good],
                   C=validationArr['angDistancesFaint'][good], reduce_C_function=np.median,
                   gridsize=20)

    ax.set_ylabel('Sun Altitude (degrees)')
    ax.set_xlabel('Moon Altitude (degrees)')
    ax.set_title('filter %s' % filterName)
    ax.axhline(-22, linestyle='--')
    ax.axvline(0, linestyle='--')
    cb = plt.colorbar(im)
    cb.set_label('Median Offset (degrees)')
    fig.savefig('Plots/medianAngDiff_%s_.pdf' % filterName)
    plt.close(fig)


    fig,ax = plt.subplots()

    im = ax.hexbin(validationArr['moonAlt'][good],validationArr['sunAlt'][good],
                   C=validationArr['angDistancesFaint'][good], reduce_C_function=robustRMS,
                   gridsize=20)

    ax.set_ylabel('Sun Altitude (degrees)')
    ax.set_xlabel('Moon Altitude (degrees)')
    ax.set_title('filter %s' % filterName)
    ax.axhline(-22, linestyle='--')
    ax.axvline(0, linestyle='--')

    cb = plt.colorbar(im)
    cb.set_label('Robust-RMS Offset (degrees)')
    fig.savefig('Plots/rmsAngDiff_%s_.pdf' % filterName)
    plt.close(fig)

    fig,ax = plt.subplots()


    im = ax.hexbin(validationArr['moonAlt'][good],validationArr['sunAlt'][good],
                   gridsize=20)
    ax.set_ylabel('Sun Altitude (degrees)')
    ax.set_xlabel('Moon Altitude (degrees)')
    ax.set_title('filter %s, %i total frames ' % (filterName, good[0].size))
    ax.axhline(-22, linestyle='--')
    ax.axvline(0, linestyle='--')

    cb = plt.colorbar(im)
    cb.set_label('Count')
    fig.savefig('Plots/countObs_%s_.pdf' % filterName)
    plt.close(fig)

    np.savez('Plots/valAr_%s.npz' % filterName,validationArr=validationArr)



# Do I need to use the origin='lower' ? YES
#x = np.arange(10)
#y=np.arange(10)
#xx,yy = np.meshgrid(x,y)
#z = xx+yy
#plt.imshow(z,extent=[xx.min(),xx.max(),yy.min(),yy.max()])
