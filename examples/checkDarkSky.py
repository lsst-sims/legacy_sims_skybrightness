import lsst.sims.skybrightness as sb
import lsst.sims.photUtils.Bandpass as Bandpass
import numpy as np
import os

# Check and see what the model returns for a dark sky at zenith


# Load LSST filters
throughPath = os.getenv('LSST_THROUGHPUTS_BASELINE')
keys = ['u', 'g', 'r', 'i', 'z', 'y']
overview_vals = [22.9, 22.3, 21.2, 20.5, 19.6, 18.6]
filters = {}
for filtername in keys:
    bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                    dtype=zip(['wave', 'trans'], [float]*2))
    tempB = Bandpass()
    tempB.setBandpass(bp['wave'], bp['trans'])
    filters[filtername] = tempB

# XXX--hmm, the zodiacal seems to be fairly important here!
sm = sb.SkyModel(moon=False, twilight=False)  # , zodiacal=False)
sm2 = sb.SkyModel(moon=False, twilight=False, zodiacal=False)
mjd = 56948.05174
sm.setRaDecMjd(np.array([0.]), np.array([90.]), mjd, azAlt=True, degrees=True)
sm2.setRaDecMjd(np.array([0.]), np.array([90.]), mjd, azAlt=True, degrees=True)
sm.computeSpec()
sm2.computeSpec()

print 'filter  ESO model ESO(no Z)  Overview Paper'
for i, key in enumerate(keys):
    print key+'     %.2f &  %.2f  &  %.2f \\\\' % (sm.computeMags(filters[key])[0], sm2.computeMags(filters[key])[0], overview_vals[i])

# Let's also output the cannon filters while we're at it:
canonFilters = {}
fnames = ['blue_canon.csv', 'green_canon.csv', 'red_canon.csv']
dataDir = os.getenv('SIMS_SKYBRIGHTNESS_DATA_DIR')
# Filter names, from bluest to reddest.
cannonKeys = ['B', 'G', 'R']
for key, fname in zip(cannonKeys, fnames):
    bpdata = np.genfromtxt(os.path.join(dataDir, 'Canon/', fname), delimiter=',',
                           dtype=zip(['wave', 'through'], [float]*2))
    bpTemp = Bandpass()
    bpTemp.setBandpass(bpdata['wave'], bpdata['through'])
    canonFilters[key] = bpTemp

print '----------'
print 'Cannon Filters'
for i, key in enumerate(cannonKeys):
    print key+'     %.2f    ' % (sm.computeMags(canonFilters[key])[0])


# lets loop over a year and try things!
times = np.arange(0, 365.25, .25)

results = np.zeros((times.size, 6), dtype=float)
sm = sb.SkyModel(moon=False, twilight=False, mags=True)
for i, time in enumerate(times):
    sm.setRaDecMjd(np.array([0.]), np.array([90.]), mjd+time, azAlt=True, degrees=True)
    if sm.sunAlt < np.radians(-12.):
        for j, f in enumerate(keys):
            sm.computeSpec()
            results[i, :] = sm.computeMags()[0]

good = np.where(results[:, 0] != 0)
results = results[good[0], :]
print ' filter   median model pm std model,  overview vals'
for i, key in enumerate(keys):
    print '$'+key+'$'+' &    %.2f $\pm$  %.2f  &  %.2f \\\\' % (np.median(results[:, i]), results[:, i].std(), overview_vals[i])
