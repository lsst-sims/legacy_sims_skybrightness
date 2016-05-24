import lsst.sims.skybrightness as sb
import lsst.sims.photUtils.Bandpass as Bandpass
import os
import numpy as np
import matplotlib.pylab as plt
import healpy as hp
# from lsst.sims.selfcal.analysis import healplots
from lsst.sims.skybrightness.utils import robustRMS, ut2Mjd, mjd2ut
from lsst.sims.utils import _altAzPaFromRaDec, haversine, calcLmstLast, _raDecFromAltAz, ObservationMetaData, Site
from scipy.interpolate import griddata

plt.rcParams.update({'axes.labelsize': 'x-large'})
plt.rcParams.update({'axes.titlesize': 'x-large'})
plt.rcParams.update({'xtick.labelsize': 'large', 'ytick.labelsize': 'large'})

# Let's recreate the delta m_5 plot from figure 3 in:
# http://xxx.lanl.gov/pdf/1510.07574.pdf
telescope = Site('LSST')
nside = 32
lat, ra = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
dec = np.pi/2-lat

kwargs = dict(twilight=False, zodiacal=False, moon=True, scatteredStar=False, mergedSpec=False)

sm = sb.SkyModel(observatory='LSST', mags=True)#, **kwargs)
mjd = 49353.177645
sm.setRaDecMjd(ra,dec,mjd)
mag = sm.returnMags()
lmst, last = calcLmstLast(mjd,telescope.longitude_rad)

moonRA, moonDec =  _raDecFromAltAz(sm.moonAlt, sm.moonAz, ObservationMetaData(mjd=mjd,site=telescope))

alt, az, pa = _altAzPaFromRaDec(ra,dec, ObservationMetaData(mjd=mjd,site=telescope))
angDist2Moon = np.degrees(haversine(az,alt, sm.moonAz,sm.moonAlt))
ang2 = np.degrees(haversine(ra,dec, moonRA,moonDec))
alt = np.degrees(alt)

mags = -0.5*(np.nanmin(mag['u'])-mag['u'])


#extent = (0,130, 0,90)
extent = (20,120, 20,90)

xs,ys = np.mgrid[extent[0]:extent[1], extent[2]:extent[3]]
resampled = griddata((angDist2Moon, alt), mags, (xs, ys))
notInf = np.where((resampled != np.inf) & (~np.isnan(resampled)))

resampled[notInf] = resampled[notInf]-resampled[notInf].max()

blah=plt.imshow(resampled.T, extent=extent, origin='lower')
#blah = plt.hexbin(angDist2Moon, alt, mags)
cb = plt.colorbar(blah)
plt.xlabel('Angular Distance to Moon (deg)')
plt.ylabel('Altitude (deg)')
plt.title('$u$-band')
plt.ylim([20,90])
plt.xlim([20,120])
cb.set_label('$\Delta m_5$')
cb.solids.set_edgecolor("face")
plt.savefig('deltam5.pdf')
plt.close('all')
