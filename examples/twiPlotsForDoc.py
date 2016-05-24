import numpy as np
import healpy as hp
import matplotlib.pylab as plt
import matplotlib

plt.rcParams.update({'axes.labelsize': 'x-large'})
plt.rcParams.update({'axes.titlesize': 'x-large'})
plt.rcParams.update({'xtick.labelsize': 'large', 'ytick.labelsize': 'large'})

# Just make a few example plots for the document
filterName = 'R'
twi = np.load('TwilightMaps/twiMaps_'+filterName+'.npz')

sunAlts = twi['sunAlts'].copy()
magMaps = twi['magMap'].copy()
rmsMaps = twi['rmsMap'].copy()

nside = hp.npix2nside(magMaps[:,0].size)
npix = magMaps[:,0].size
hpid = np.arange(magMaps[:,0].size)

lat, az = hp.pix2ang(nside, np.arange(npix))
alt = np.pi/2.-lat
airmass = 1./np.cos(np.pi/2.-alt)


useAlts = [-12.5, -14, -16.5, -18.]


# Make the masked data white
cmap0 = matplotlib.cm.jet
newcm = matplotlib.colors.LinearSegmentedColormap('newcm',cmap0._segmentdata,cmap0.N)
newcm.set_over(newcm(1.0))
newcm.set_under('w')
newcm.set_bad('w')


for i,alt in enumerate(useAlts):
    fig = plt.figure(1)
    good = np.where(sunAlts == np.radians(alt))[0]
    toShow = magMaps[:,good[0]]
    toShow[np.where((airmass < 1) | (airmass > 10))] = hp.UNSEEN
    hp.mollview(toShow, rot=(0,90),
                title=' Sun Alt = %.1f deg' % alt, fig=1, min=4.5,max=9.5,
                unit='%s mags' % filterName, sub=(2,2,i+1), cmap=newcm)

#fig.tight_layout()
fig.savefig('Plots/twiExamples.pdf')
plt.close(fig)
