from builtins import zip
import numpy as np
import glob
from astropy.io import fits
import os
from lsst.sims.photUtils import Sed, Bandpass

dataDir = os.getenv('SIMS_SKYBRIGHTNESS_DATA_DIR')
outDir = os.path.join(dataDir, 'ESO_Spectra/Airglow')

# Read in all the spectra from ESO call and package into a single npz file

files = glob.glob('skytable*.fits')

temp = fits.open(files[0])
wave = temp[1].data['lam'].copy()*1e3

airmasses = []
solarFlux = []
specs = []

for i, filename in enumerate(files):
    fitsfile = fits.open(filename)
    if np.max(fitsfile[1].data['flux']) > 0:
        specs.append(fitsfile[1].data['flux'].copy())
        header = fitsfile[0].header['comment']
        for card in header:
            if 'SKYMODEL.TARGET.AIRMASS' in card:
                airmasses.append(float(card.split('=')[-1]))
            elif 'SKYMODEL.MSOLFLUX' in card:
                solarFlux.append(float(card.split('=')[-1]))


airmasses = np.array(airmasses)
solarFlux = np.array(solarFlux)

nrec = airmasses.size
nwave = wave.size

dtype = [('airmass', 'float'),
         ('solarFlux', 'float'),
         ('spectra', 'float', (nwave)), ('mags', 'float', (6))]
Spectra = np.zeros(nrec, dtype=dtype)
Spectra['airmass'] = airmasses
Spectra['solarFlux'] = solarFlux
Spectra['spectra'] = specs


hPlank = 6.626068e-27  # erg s
cLight = 2.99792458e10  # cm/s

# Convert spectra from ph/s/m2/micron/arcsec2 to erg/s/cm2/nm/arcsec2
Spectra['spectra'] = Spectra['spectra']/(100.**2)*hPlank*cLight/(wave*1e-7)/1e3

# Sort things since this might be helpful later
Spectra.sort(order=['airmass', 'solarFlux'])

# Load LSST filters
throughPath = os.getenv('LSST_THROUGHPUTS_BASELINE')
keys = ['u', 'g', 'r', 'i', 'z', 'y']
nfilt = len(keys)
filters = {}
for filtername in keys:
    bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                    dtype=list(zip(['wave', 'trans'], [float]*2)))
    tempB = Bandpass()
    tempB.setBandpass(bp['wave'], bp['trans'])
    filters[filtername] = tempB

filterWave = np.array([filters[f].calcEffWavelen()[0] for f in keys])

for i, spectrum in enumerate(Spectra['spectra']):
    tempSed = Sed()
    tempSed.setSED(wave, flambda=spectrum)
    for j, filtName in enumerate(keys):
        try:
            Spectra['mags'][i][j] = tempSed.calcMag(filters[filtName])
        except:
            pass

np.savez(os.path.join(outDir, 'airglowSpectra.npz'), wave=wave, spec=Spectra, filterWave=filterWave)
