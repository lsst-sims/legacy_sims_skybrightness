import lsst.sims.skybrightness as sb
import numpy as np
import healpy as hp

sm = sb.SkyModel(mags=True)
nside = 16
hpmap = np.zeros(hp.nside2npix(nside))
lat, ra = hp.pix2ang(nside, np.arange(hpmap.size))
dec = np.pi/2-lat

mjd = 60291.35423611111
sm.setRaDecMjd(ra, dec, mjd)
mags = sm.returnMags()

