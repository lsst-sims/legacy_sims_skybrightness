import numpy as np
import glob
import pyfits
import os

dataDir = os.getenv('SIMS_SKYBRIGHTNESS_DATA_DIR')
outDir = os.path.join(dataDir, 'ESO_Spectra/LowerAtm')

# Read in all the spectra from ESO call and package into a single npz file

files = glob.glob('skytable*.fits')

temp = pyfits.open(files[0])
wave = temp[1].data['lam'].copy()*1e3

airmasses = []
nightTimes = []
specs = []

for i,filename in enumerate(files):
    fits = pyfits.open(filename)
    if np.max(fits[1].data['flux']) > 0:
        specs.append(fits[1].data['flux'].copy())
        header = fits[0].header['comment']
        for card in header:
            if 'SKYMODEL.TARGET.AIRMASS' in card:
                airmasses.append(float(card.split('=')[-1]))
            elif 'SKYMODEL.TIME' in card:
                nightTimes.append(float(card.split('=')[-1]))


airmasses = np.array(airmasses)
nigtTimes = np.array(nightTimes)

nrec = airmasses.size
nwave = wave.size

dtype = [('airmass', 'float'),
         ('nightTimes', 'float'),
         ('spectra', 'float', (nwave))]
Spectra = np.zeros(nrec, dtype=dtype)
Spectra['airmass'] = airmasses
Spectra['nightTimes'] = nightTimes
Spectra['spectra'] = specs


hPlank = 6.626068e-27 # erg s
cLight = 2.99792458e10 # cm/s

# Convert spectra from ph/s/m2/micron/arcsec2 to erg/s/cm2/nm/arcsec2
Spectra['spectra'] = Spectra['spectra']/(100.**2)*hPlank*cLight/(wave*1e-7)/1e3

# Sort things since this might be helpful later
Spectra.sort(order=['airmass','nightTimes'])

np.savez(os.path.join(outDir,'Spectra.npz'), wave = wave, spec=Spectra)
