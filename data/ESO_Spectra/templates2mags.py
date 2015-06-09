import numpy as np
import os
import glob
from lsst.sims.photUtils import Sed,Bandpass

# Read in the template spectra and compute the ugrizy magnitudes for them and save the results


dirs = ['Airglow','MergedSpec','ScatteredStarLight','Zodiacal',
        'LowerAtm','Moon','UpperAtm']



dataDir =  os.path.join(os.environ.get('SIMS_SKYBRIGHTNESS_DATA_DIR'), 'ESO_Spectra/')


# Load LSST filters
throughPath = os.getenv('LSST_THROUGHPUTS_BASELINE')
keys = ['u','g','r','i','z','y']
nfilt = len(keys)
filters = {}
for filtername in keys:
    bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                    dtype=zip(['wave','trans'],[float]*2 ))
    tempB = Bandpass()
    tempB.setBandpass(bp['wave'],bp['trans'])
    filters[filtername] = tempB

newWave = np.array([filters[f].calcEffWavelen()[0] for f in keys ])


for directory in dirs:
    dataDir =  os.path.join(os.environ.get('SIMS_SKYBRIGHTNESS_DATA_DIR'), 'ESO_Spectra/'+directory)
    filenames = glob.glob(dataDir+'/*.npz')
    if len(filenames) == 1:
        temp = np.load(filenames[0])
        wave = temp['wave'].copy()
        spec = temp['spec'].copy()
    else:
        temp = np.load(filenames[0])
        wave = temp['wave'].copy()
        spec = temp['spec'].copy()
        for filename in filenames[1:]:
            temp = np.load(filename)
            spec = np.append(spec, temp['spec'])
    # Construct the new dtype for the mag array
    dt = []
    for ack in spec.dtype.names:
        if ack != 'spectra':
            dt.append((ack,'<f8'))
    dt.append(('mags', '<f8', (nfilt,)))

    magArr = np.zeros(spec.size, dtype=dt)
    # Copy over all the other info
    for dtname in spec.dtype.names:
        if dtname != 'spectra':
            magArr[dtname] = spec[dtname]

    for i,spectrum in enumerate(spec['spectra']):
        tempSed = Sed()
        tempSed.setSED(wave,flambda=spectrum)
        for j,filtName in enumerate(keys):
            try:
                magArr['mags'][i][j] = tempSed.calcMag(filters[filtName])
            except:
                pass


    np.savez(os.path.join(dataDir, 'mags','mags.npz'), wave=newWave, magArr=magArr)
