import numpy as np
import lsst.sims.skybrightness as sb
import healpy as hp
import matplotlib.pylab as plt

# Load up the spectra and plot some examples of each component

comps = ['twilight', 'zodiacal', 'moon','airglow','lowerAtm',
         'upperAtm','scatteredStar', 'mergedSpec']

sm = sb.SkyModel(lowerAtm=True, upperAtm=True, scatteredStar=True)

nside = 8
hpmap = np.zeros(hp.nside2npix(nside))
lat, ra = hp.pix2ang(nside, np.arange(hpmap.size))
dec = np.pi/2-lat
sm.setRaDecMjd(ra,dec,49353.177645, degrees=False, azAlt=True)

# Zodical
for comp in comps:
    setattr(sm,comp,False)
sm.zodiacal = True
sm.computeSpec()

fig, ax = plt.subplots()
ax.semilogy(sm.wave, sm.spec[0], label='airmass = %.1f, $l=$%.1f, $b$=%.1f' % (sm.airmass[0],np.degrees(sm.points['azEclipRelSun'][0]), np.degrees(sm.points['altEclip'][0]) ), alpha=.5)
ax.semilogy(sm.wave, sm.spec[100], label='airmass = %.1f, $l=$%.1f, $b$=%.1f' % (sm.airmass[100], np.degrees(sm.points['azEclipRelSun'][100]), np.degrees(sm.points['altEclip'][100])), alpha=.5)
ax.set_ylim([1e-19,1e-16])
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Flux (erg/s/cm$^2$/nm/arcsec$^2$)')
ax.set_title('Zodiacal')
ax.legend()
fig.savefig('zodiacal.pdf')
plt.close(fig)

# Airglow
for comp in comps:
    setattr(sm,comp,False)
sm.airglow = True
sm.computeSpec()

fig, ax = plt.subplots()
ax.semilogy(sm.wave, sm.spec[0], label='airmass = %.1f, sfu=%i' % (sm.airmass[0], sm.solarFlux),
            alpha=.5)

sm.setRaDecMjd(ra,dec,49353.177645, degrees=False, azAlt=True, solarFlux=200)
sm.computeSpec()

ax.semilogy(sm.wave, sm.spec[0], label='airmass = %.1f, sfu=%i' % (sm.airmass[0], sm.solarFlux),
            alpha=.5)

sm.setRaDecMjd(ra,dec,49353.177645, degrees=False, azAlt=True)

ax.set_ylim([1e-19,1e-15])
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Flux (erg/s/cm$^2$/nm/arcsec$^2$)')
ax.set_title('Airglow')
ax.legend(loc='upper left')
fig.savefig('airglow.pdf')
plt.close(fig)



alpha = 0.5

fig, ax =  plt.subplots(4,1,sharex=True)
# Merged
for comp in comps:
    setattr(sm,comp,False)

sm.lowerAtm = True
sm.computeSpec()
ax[0].semilogy(sm.wave, sm.spec[0], alpha=alpha, label='airmass = %.1f' % sm.airmass[0])
ax[0].semilogy(sm.wave, sm.spec[200], alpha=alpha, label='airmass = %.1f' % sm.airmass[200])
ax[0].legend(loc='upper left')

ax[0].set_title('Lower Atmosphere')

sm.lowerAtm = False
sm.upperAtm = True
sm.computeSpec()
ax[1].plot(sm.wave, sm.spec[0] , alpha=alpha)
ax[1].plot(sm.wave, sm.spec[200] , alpha=alpha)
ax[1].set_title('Upper Atmosphere')

sm.upperAtm = False
sm.scatteredStar = True
sm.computeSpec()
ax[2].semilogy(sm.wave, sm.spec[0] , alpha=alpha)
ax[2].semilogy(sm.wave, sm.spec[200] , alpha=alpha)
ax[2].set_title('Scattered Star Light')
sm.scatteredStar = False
sm.mergedSpec = True
sm.computeSpec()
ax[3].semilogy(sm.wave, sm.spec[0], alpha=alpha)
ax[3].semilogy(sm.wave, sm.spec[200], alpha=alpha)
ax[3].set_title('Merged')


for oneAx in ax:
    oneAx.set_ylim([1e-19,1e-12])

ax[2].set_ylabel('Flux (erg/s/cm$^2$/nm/arcsec$^2$)')

ax[-1].set_xlabel('Wavelength (nm)')
fig.savefig('merged.pdf')
plt.close(fig)

fig, ax =  plt.subplots()
for comp in comps:
    setattr(sm,comp,False)
sm.moon = True

sm.setRaDecMjd(ra,dec,49353.177645, degrees=False, azAlt=True)
sm.computeSpec()

am = sm.airmass[0]
moonAlt = sm.points['moonAltitude'][0]
moonPhase = sm.points['moonSunSep'][0]
moonAz = np.degrees(sm.points['azRelMoon'][0])

ax.semilogy(sm.wave,sm.spec[0], alpha=0.5, label='$X=$%.1f, Moon Alt = %.1f, \n Moon Phase = %.1f, Azimuth to Moon=%.1f' % (am, moonAlt, moonPhase, moonAz))

sm.setRaDecMjd(ra,dec,49373.177645, degrees=False, azAlt=True)
sm.computeSpec()
am = sm.airmass[200]
moonAlt = sm.points['moonAltitude'][200]
moonPhase = sm.points['moonSunSep'][200]
moonAz = np.degrees(sm.points['azRelMoon'][200])

ax.semilogy(sm.wave,sm.spec[200], alpha=0.5, label='$X=$%.1f, Moon Alt = %.1f, \n Moon Phase = %.1f, Azimuth to Moon=%.1f' % (am, moonAlt, moonPhase, moonAz))
ax.legend()
ax.set_title('Moon')
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Flux (erg/s/cm$^2$/nm/arcsec$^2$)')
ax.set_ylim([1e-20,1e-12])

fig.savefig('moon.pdf')
plt.close(fig)
