import numpy as np
import lsst.sims.skybrightness as sb
from scipy.stats import binned_statistic


# Make some validation plots of the model compared to the observations we have

# Select all the near-zenith photometry from the Canon camera
sqlQ = 'select stars.ra, stars.dec, obs.sky, obs.filter, dates.mjd  from obs,dates,stars where obs.starID = stars.ID and obs.dateID=dates.ID and obs.alt > 85;'
skydata, mjd = sb.allSkyDB(2744 , sqlQ=sqlQ, filt=None,
                           dtypes=zip(['ra','dec','sky','filt','mjd'],
                                      [float,float,float,'|S1',float]))


# Maybe bin things on 15-min timescales to cut down the number of calls I need to make to the model
binsize = 15./60./24.
edges = np.arange(skydata['mjd'].min(),skydata['mjd'].max()+binsize*2, binsize)

filters = ['R','G','B']
filters = ['R']

sm = sb.SkyModel(mags=True)

for filterName in filters:
    good = np.where(skydata['filt'] == filterName)
    mjds,junk1,junk2 = binned_statistic(skydata['mjd'][good], skydata['mjd'][good], bins=edges)
    skymags,junk1,junk2 = binned_statistic(skydata['mjd'][good], skydata['sky'][good], bins=edges)
    ras, junk1,junk2 = binned_statistic(skydata['mjd'][good], skydata['ra'][good],
                                        bins=edges, statistic='median')
    decs, junk1,junk2 = binned_statistic(skydata['mjd'][good], skydata['dec'][good],
                                        bins=edges, statistic='median')

    empty = np.isnan(mjds+skymags)
    mjds = mjds[~empty]
    skymags = skymags[~empty]
    ras = ras[~empty]
    decs = decs[~empty]

    results = np.zeros(mjds.size, dtype=zip(['mjd','sunAlt', 'moonAlt', 'modelMag', 'obsMag'],[float]*5) )
    results['obsMag'] = skymags
    i=0
    for mjd,ra,dec in zip(mjds,ras,decs):
        sm.setRaDecMjd(np.array([ra]),np.array([dec]), mjd, degrees=True)
        sm.computeSpec()
        mags = sm.computeMags()
        results['modelMag'][i] = mags[0][2] # for r band
        results['sunAlt'][i] = np.degrees(sm.sunAlt)
        results['moonAlt'][i] = np.degrees(sm.moonAlt)
        i += 1
