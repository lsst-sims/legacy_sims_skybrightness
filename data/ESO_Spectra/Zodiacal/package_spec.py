from builtins import str
from builtins import zip
import numpy as np
import glob
from astropy.io import fits
import healpy as hp
import os
from lsst.sims.photUtils import Sed, Bandpass

dataDir = os.getenv('SIMS_SKYBRIGHTNESS_DATA_DIR')
outDir = os.path.join(dataDir, 'ESO_Spectra/Zodiacal')

nside = 4

# Read in all the spectra from ESO call and package into a single npz file

files = glob.glob('skytable*.fits')

temp = fits.open(files[0])
wave = temp[1].data['lam'].copy()*1e3

airmasses = []
eclLon = []
eclLat = []
specs = []

for i, filename in enumerate(files):
    fitsname = fits.open(filename)
    if np.max(fitsname[1].data['flux']) > 0:
        specs.append(fitsname[1].data['flux'].copy())
        header = fitsname[0].header['comment']
        for card in header:
            if 'SKYMODEL.TARGET.AIRMASS' in card:
                airmasses.append(float(card.split('=')[-1]))
            elif 'SKYMODEL.ECL.LON' in card:
                eclLon.append(float(card.split('=')[-1]))
            elif 'SKYMODEL.ECL.LAT' in card:
                eclLat.append(float(card.split('=')[-1]))


airmasses = np.array(airmasses)
eclLon = np.array(eclLon)
eclLat = np.array(eclLat)

wrapA = np.where(eclLon < 0.)
eclLon[wrapA] = eclLon[wrapA]+360.

uAM = np.unique(airmasses)
nAM = uAM.size
nwave = wave.size

dtype = [('airmass', 'float'),
         ('hpid', 'int'),
         ('spectra', 'float', (nwave)), ('mags', 'float', (6))]
npix = hp.nside2npix(nside)
Spectra = np.zeros(nAM*npix, dtype=dtype)
for i, am in enumerate(uAM):
    Spectra['airmass'][i*npix:i*npix+npix] = am
    Spectra['hpid'][i*npix:i*npix+npix] = np.arange(npix)


for am, lat, lon, spec in zip(airmasses, eclLat, eclLon, specs):
    hpid = hp.ang2pix(nside, np.radians(lat+90.), np.radians(lon))
    good = np.where((Spectra['airmass'] == am) & (Spectra['hpid'] == hpid))
    Spectra['spectra'][good] = spec.copy()


#Spectra['airmass'] = airmasses
#Spectra['eclLon'] = eclLon
#Spectra['eclLat'] = eclLat
#Spectra['spectra'] = specs


hPlank = 6.626068e-27  # erg s
cLight = 2.99792458e10  # cm/s

# Convert spectra from ph/s/m2/micron/arcsec2 to erg/s/cm2/nm/arcsec2
Spectra['spectra'] = Spectra['spectra']/(100.**2)*hPlank*cLight/(wave*1e-7)/1e3

# Sort things since this might be helpful later
Spectra.sort(order=['airmass', 'hpid'])

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


# span this over multiple files to store in github
nbreak = 3
nrec = np.size(Spectra)

for i in np.arange(nbreak):
    np.savez(os.path.join(outDir, 'zodiacalSpectra_'+str(i)+'.npz'), wave=wave,
             spec=Spectra[i*nrec/nbreak:(i+1)*nrec/nbreak], filterWave=filterWave)
