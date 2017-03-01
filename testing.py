import lsst.sims.skybrightness as sb
import numpy as np
import healpy as hp

sm = sb.SkyModel(mags=True)
nside = 16
hpmap = np.zeros(hp.nside2npix(nside))
lat, ra = hp.pix2ang(nside, np.arange(hpmap.size))
dec = np.pi/2-lat

#mjd = 60291.35423611111 -.2
mjd = 49353.177645
sm.setRaDecMjd(ra, dec, mjd)
mags = sm.returnMags()
mags

sm = sb.SkyModel(twilight=True, zodiacal=False, moon=True,
                 airglow=True, lowerAtm=False, upperAtm=False, scatteredStar=False,
                 mergedSpec=True, mags=True)

sm.setRaDecMjd(ra, dec, mjd)
mags = sm.returnMags()
mags


# OK, moon does not appear to be constributing when it should be huge on the linux machine.
# airglow and mergedSpec seem to be fine.
# 


# twilight check
mjd = 49353.177645-.14
sm = sb.SkyModel(twilight=True, zodiacal=False, moon=False,
                 airglow=True, lowerAtm=False, upperAtm=False, scatteredStar=False,
                 mergedSpec=True, mags=True)

sm.setRaDecMjd(ra, dec, mjd)
mags = sm.returnMags()
mags

# Twilight is fine.


