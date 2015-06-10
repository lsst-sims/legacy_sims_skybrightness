import numpy as np
import os

dataDir = os.getenv('SIMS_SKYBRIGHTNESS_DATA_DIR')
outDir = os.path.join(dataDir, 'ESO_Spectra/MergedSpec')

# A large number of the background components only depend on Airmass, so we can merge those together

npzs = ['Airglow/airglowSpectra.npz', 'LowerAtm/Spectra.npz',
        'ScatteredStarLight/scatteredStarLight.npz',
        'UpperAtm/Spectra.npz']
files = [os.path.join(dataDir, 'ESO_Spectra', npz) for npz in npzs]


temp = np.load(files[0])

wave = temp['wave'].copy()
spec = temp['spec'].copy()
spec['spectra'] = spec['spectra']*0.
spec['mags'] = spec['mags']*0.

for filename in files:
    restored = np.load(filename)
    spec['spectra'] += restored['spec']['spectra']
    try:
        flux = 10.**(-0.4*(restored['spec']['mags']-np.log10(3631.)))
    except:
        import pdb ; pdb.set_trace()
    flux[np.where(restored['spec']['mags'] == 0.)] = 0.
    spec['mags'] += flux

spec['mags'] = -2.5*np.log10(spec['mags'])+np.log10(3631.)

np.savez(os.path.join(outDir,'mergedSpec.npz'), spec=spec, wave=wave, filterWave=temp['filterWave'])
