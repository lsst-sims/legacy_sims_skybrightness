from __future__ import print_function
from builtins import zip
import numpy as np
import ephem
import lsst.sims.skybrightness as sb
import healpy as hp
import os
import sqlalchemy as sqla
from lsst.sims.maf.utils.telescopeInfo import TelescopeInfo
from lsst.sims.skybrightness.utils import wrapRA, mjd2djd, robustRMS
from lsst.sims.utils import _altAzPaFromRaDec, Site, _healbin, ObservationMetaData
from lsst.utils import getPackageDir


# Set up the telescope:
telescope = Site('LSST')
Observatory = ephem.Observer()
Observatory.lat = telescope.latitude_rad
Observatory.lon = telescope.longitude_rad
Observatory.elevation = telescope.height
sun = ephem.Sun()
moon = ephem.Moon()

# Let's bin things up in terms of sun altitude.
sunAlts = np.radians(np.arange(-5, -23, -0.25))
altBin = np.radians(0.1)

nside = 8
npix = hp.nside2npix(nside)
magMap = np.zeros((npix, sunAlts.size), dtype=float)
rmsMap = np.zeros((npix, sunAlts.size), dtype=float)

filterNames = ['R', 'G', 'B']
#filterNames = ['R']
#sunAlts = [sunAlts[5]]

for filterName in filterNames:

    dataPath = getPackageDir('SIMS_SKYBRIGHTNESS_DATA')
    dbAddress = 'sqlite:///'+os.path.join(dataPath, 'photometry', 'skydata.sqlite')

    names = ['mjd', 'ra', 'dec', 'alt', 'starMag', 'sky', 'filter']
    types = [float, float, float, float, float, float, '|S1']
    dtypes = list(zip(names, types))

    engine = sqla.create_engine(dbAddress)
    connection = engine.raw_connection()
    cursor = connection.cursor()

    for i, ack in enumerate(sunAlts):
        q = 'select dates.mjd, stars.ra, stars.dec, obs.alt, obs.starMag, obs.sky, obs.filter from obs,stars,dates where obs.starID = stars.ID and obs.dateID = dates.ID and obs.filter = "%s" and obs.dateID in (select ID from dates where sunAlt >= %f and sunAlt <= %f)' % (filterName, sunAlts[
            i]-altBin, sunAlts[i]+altBin)

        print('Executing:')
        print(q)
        print('%i of %i' % (i, np.size(sunAlts)))

        cursor.execute(q)
        data = cursor.fetchall()
        data = np.asarray(data, dtype=dtypes)

        print('got %i results' % data.size)

        data['ra'] = np.radians(data['ra'])
        data['dec'] = np.radians(data['dec'])

        data.sort(order='mjd')

        umjd = np.unique(data['mjd'])

        left = np.searchsorted(data['mjd'], umjd)
        right = np.searchsorted(data['mjd'], umjd, side='right')

        altaz = np.zeros(data.size, dtype=list(zip(['alt', 'az'], [float]*2)))
        moonAlt = np.zeros(data.size, dtype=float)

        print('computing alts and azs')

        for j, (le, ri, mjd) in enumerate(zip(left, right, umjd)):
            Observatory.date = mjd2djd(mjd)
            sun.compute(Observatory)
            obs_metadata = ObservationMetaData(pointingRA=np.degrees(sun.ra),
                                               pointingDec=np.degrees(sun.dec),
                                               rotSkyPos=np.degrees(0),
                                               mjd=mjd)
            alt, az, pa = _altAzPaFromRaDec(data['ra'][le:ri], data['dec'][le:ri], obs_metadata)
            az = wrapRA(az - sun.az)
            altaz['alt'][le:ri] += alt
            altaz['az'][le:ri] += az
            moon.compute(Observatory)
            moonAlt[le:ri] += moon.alt

        print('making maps')
        good = np.where(moonAlt < 0)
        magMap[:, i] = _healbin(altaz['az'][good], altaz['alt'][good], data['sky'][good],
                                nside=nside, reduceFunc=np.median)
        rmsMap[:, i] = _healbin(altaz['az'][good], altaz['alt'][good], data['sky'][good],
                                nside=nside, reduceFunc=robustRMS)

    print('saving maps')
    np.savez('TwilightMaps/twiMaps_%s.npz' % filterName, magMap=magMap, rmsMap=rmsMap, sunAlts=sunAlts)
