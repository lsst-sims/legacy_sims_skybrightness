from __future__ import print_function
from builtins import str
from builtins import zip
import numpy as np
import glob
from astropy.io import fits
import os
from lsst.sims.photUtils import Sed, Bandpass

dataDir = os.getenv('SIMS_SKYBRIGHTNESS_DATA_DIR')
outDir = os.path.join(dataDir, 'ESO_Spectra/Moon')

# Read in all the spectra from ESO call and package into a single npz file

files = glob.glob('skytable*.fits')

temp = fits.open(files[0])
moonWave = temp[1].data['lam'].copy()*1e3

moonSpec = []
moonAM = []  # Actually target airmass
moonAlt = []  # altitude of the moon
moonSunSep = []  # moon Phase
moonTargetSep = []
hpid = []
for i, filename in enumerate(files):
    fitsname = fits.open(filename)
    try:
        if np.max(fitsname[1].data['flux']) > 0:
            moonSpec.append(fitsname[1].data['flux'].copy())
            header = fitsname[0].header['comment']
            for card in header:
                if 'SKYMODEL.MOON.SUN.SEP' in card:
                    moonSunSep.append(float(card.split('=')[-1]))
                elif 'SKYMODEL.TARGET.AIRMASS' in card:
                    #moonAM.append( 1./np.cos(np.radians(90.-float(card.split('=')[-1]))) )
                    moonAM.append(float(card.split('=')[-1]))
                elif 'SKYMODEL.MOON.TARGET.SEP' in card:
                    moonTargetSep.append(float(card.split('=')[-1]))
                elif 'SKYMODEL.MOON.ALT' in card:
                    moonAlt.append(float(card.split('=')[-1]))
    except:
        print(filename, ' Failed')


import healpy as hp
from lsst.sims.utils import haversine

nside = 4
lat, az = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
alt = np.pi/2.-lat
airmass = 1./np.cos(np.pi/2.-alt)


# Only need low airmass and then 1/2 to sky
good = np.where((az >= 0) & (az <= np.pi) & (airmass <= 2.6) & (airmass >= 1.))
airmass = airmass[good]
alt = alt[good]
az = az[good]


moonAM = np.array(moonAM)
moonAlt = np.array(moonAlt)
moonSunSep = np.array(moonSunSep)
moonTargetSep = np.array(moonTargetSep)
moonAzDiff = moonTargetSep*0
targetAlt = np.pi/2.-np.arccos(1./moonAM)
# Compute the azimuth difference given the moon-target-seperation
# Let's just do a stupid loop:
for i in np.arange(targetAlt.size):
    possibleDistances = haversine(0., np.radians(moonAlt[i]), az, az*0+targetAlt[i])
    diff = np.abs(possibleDistances - np.radians(moonTargetSep[i]))
    good = np.where(diff == diff.min())
    try:
        moonAzDiff[i] = az[good][0]
        # ok, now I have an alt and az, I can convert that back to a healpix id.

        hpid.append(hp.ang2pix(nside, np.pi/2.-targetAlt[i], moonAzDiff[i]))
    except:
        import pdb
        pdb.set_trace()
    if diff.min() > 1e-5:
        import pdb
        pdb.set_trace()


nrec = moonAM.size
nwave = moonWave.size

dtype = [('hpid', 'int'),
         ('moonAltitude', 'float'),
         ('moonSunSep', 'float'),
         ('spectra', 'float', (nwave)), ('mags', 'float', (6))]
moonSpectra = np.zeros(nrec, dtype=dtype)
moonSpectra['hpid'] = hpid
moonSpectra['moonAltitude'] = moonAlt
moonSpectra['moonSunSep'] = moonSunSep
moonSpectra['spectra'] = moonSpec

hPlank = 6.626068e-27  # erg s
cLight = 2.99792458e10  # cm/s


# Convert spectra from ph/s/m2/micron/arcsec2 to erg/s/cm2/nm/arcsec2
moonSpectra['spectra'] = moonSpectra['spectra']/(100.**2)*hPlank*cLight/(moonWave*1e-7)/1e3

# Sort things since this might be helpful later
moonSpectra.sort(order=['moonSunSep', 'moonAltitude', 'hpid'])

# Crop off the incomplete ones

good = np.where((moonSpectra['moonAltitude'] >= 0) & (moonSpectra['moonAltitude'] < 89.))
moonSpectra = moonSpectra[good]


# Load LSST filters
throughPath = os.getenv('LSST_THROUGHPUTS_BASELINE')
keys = ['u', 'g', 'r', 'i', 'z', 'y']
nfilt = len(keys)
filters = {}
for filtername in keys:
    bp = np.loadtxt(os.path.join(throughPath, 'total_'+filtername+'.dat'),
                    dtype=list(zip(['wave', 'trans'], [float]*2)))
    tempB = Bandpass()
    tempB.setBandpass(bp['wave'], bp['trans'])
    filters[filtername] = tempB

filterWave = np.array([filters[f].calcEffWavelen()[0] for f in keys])

for i, spectrum in enumerate(moonSpectra['spectra']):
    tempSed = Sed()
    tempSed.setSED(moonWave, flambda=spectrum)
    for j, filtName in enumerate(keys):
        try:
            moonSpectra['mags'][i][j] = tempSed.calcMag(filters[filtName])
        except:
            pass


nbreak = 5
indices = np.linspace(0, moonSpectra.size, nbreak+1, dtype=int)

for i in np.arange(nbreak):
    break_index = np.int(np.floor(nrec/nbreak))
    np.savez(os.path.join(outDir, 'moonSpectra_'+str(i)+'.npz'), wave=moonWave,
             spec=moonSpectra[indices[i]:indices[i+1]], filterWave=filterWave)


