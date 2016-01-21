import numpy as np
import lsst.sims.skybrightness as sb
import healpy as hp
from lsst.sims.utils import Site
import ephem
from lsst.sims.utils import _altAzPaFromRaDec, haversine, Site, ObservationMetaData
from scipy.spatial import cKDTree as kdtree
import sys
import warnings
warnings.filterwarnings("ignore")

# Let's generate some arrays that could be used for lookup/lookahead


def treexyz(ra, dec):
    """Calculate x/y/z values for ra/dec points, ra/dec in radians."""
    # Note ra/dec can be arrays.
    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)
    return x, y, z


# Define avoidance spots
zenith_avoid = np.radians(5.)  # rad
moon_avoid =  np.radians(35.) # rad
venus_avoid =  np.radians(5.) # rad

# Sun twilight limit
sun_alt_limit = np.radians(-12.) # rad

# Telescope pointing altitude limit.
alt_limit = np.radians(20.)

nside = 16
map_size = hp.nside2npix(nside)
dec,ra = hp.pix2ang(nside, np.arange(map_size))
dec = np.pi/2. - dec

leafsize = 100
x,y,z = treexyz(ra, dec)
healTree = kdtree(zip(x,y,z),leafsize=leafsize)

# Need to build a kdtree to quick search for healpixes within a radius

mjd_start = 59580.033829
survey_length = .01 #years
time_step = 10. # minitues

mjds = np.arange(mjd_start,
                 mjd_start + survey_length*365.25+time_step/60./24.,
                 time_step/60./24.)

names=['u','g','r','i','z','y','airmass', 'mask']
types = [float]*7
types.append(int)
results = np.zeros((mjds.size,map_size), dtype=zip(names,types) )
results['mask'] = 1
site = Site('LSST')


moon = ephem.Moon()
venus = ephem.Venus()
sun = ephem.Sun()
doff = ephem.Date(0)-ephem.Date('1858/11/17')

lsstObs = ephem.Observer()
lsstObs.lat = np.radians(site.latitude)
lsstObs.lon = np.radians(site.longitude)
lsstObs.elevation = site.height

sm = sb.SkyModel(mags=True)

filterDict = {'u':0,'g':1,'r':2,'i':3,'z':4,'y':5}
indices = np.arange(mjds.size)
loopSize = mjds.size

for i,mjd,djd in zip(indices,mjds,mjds-doff):
    lsstObs.date = djd
    sun.compute(lsstObs)
    if sun.alt > sun_alt_limit:
        results['mask'][i] = 0
    else:
        moon.compute(djd)
        moon_mask_indices = healTree.query_ball_point(treexyz(moon.ra,moon.dec),moon_avoid)
        venus.compute(djd)
        venus_mask_indices = healTree.query_ball_point(treexyz(venus.ra,venus.dec),venus_avoid)
        results['mask'][i][moon_mask_indices] = 0
        results['mask'][i][venus_mask_indices] = 0

        alt,az,pa = _altAzPaFromRaDec(ra, dec, ObservationMetaData(mjd=mjd,site=site))
        zenith_mask = np.where(alt > np.radians(90.-zenith_avoid))
        results['mask'][i][zenith_mask] = 0

        # Let's only compute airmass for good ones
        high_enough = np.where( (alt > np.radians(alt_limit)))
        too_low = np.where( (alt <= np.radians(alt_limit)))

        results['airmass'][i][high_enough] = 1./np.cos(np.pi/2.-alt[high_enough])
        results['mask'][i][too_low] = 0

        sm.setRaDecMjd(ra[high_enough],dec[high_enough],mjd)
        skyMags = sm.returnMags()
        for key in filterDict.keys():
            results[key][i][high_enough] = skyMags[:,filterDict[key]]
        percentComplete = float(i)/loopSize*100
        if percentComplete == i/loopSize*100:
            sys.stdout.write('%i\% complete' % percentComplete)
            sys.stdout.flush()

sys.stdout.write('\n')



# Can cut off any pixels that are always masked.
collapseMask = np.sum(results['mask'], axis=0)
unmaskedIDs = np.where(collapseMask != 0)[0]
results = results[:,unmaskedIDs]

# Chop off mjds where everything is masked
collapseMask = np.sum(results['mask'], axis=1)
unmaskedIDs = np.where(collapseMask != 0)[0]
results = results[unmaskedIDs,:]
mjds = mjds[unmaskedIDs]

hpids = np.arange(map_size)
hpids = hpids[unmaskedIDs]

np.savez('testsave.npz', precomputed=results, mjds=mjds, hpids=hpids)
