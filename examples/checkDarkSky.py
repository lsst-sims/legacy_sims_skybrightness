import lsst.sims.skybrightness as sb
import lsst.sims.photUtils.Bandpass as Bandpass
import numpy as np
import os

# Check and see what the model returns for a dark sky at zenith


# Load LSST filters
throughPath = os.getenv('LSST_THROUGHPUTS_BASELINE')
keys = ['u','g','r','i','z','y']
overview_vals = [22.9,22.3,21.2,20.5,19.6,18.6]
filters = {}
for filtername in keys:
    bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                    dtype=zip(['wave','trans'],[float]*2 ))
    tempB = Bandpass()
    tempB.setBandpass(bp['wave'],bp['trans'])
    filters[filtername] = tempB

# XXX--hmm, the zodiacal seems to be fairly important here!
sm = sb.SkyModel( moon=False, twilight=False)#, zodiacal=False)
mjd = 56948.05174
sm.setRaDecMjd(np.array([0.]), np.array([90.]), mjd,azAlt=True, degrees=True)
sm.computeSpec()
print 'filter  ESO model   Overview Paper'
for i,key in enumerate(keys):
    print key+'     %.2f      %.2f ' %(sm.computeMags(filters[key])[0], overview_vals[i])

# Let's also output the cannon filters while we're at it:
canonFilters = {}
fnames = ['blue_canon.csv', 'green_canon.csv','red_canon.csv']
dataDir = os.getenv('SIMS_SKYBRIGHTNESS_DATA_DIR')
# Filter names, from bluest to reddest.
cannonKeys = ['B','G','R']
for key,fname in zip(cannonKeys,fnames):
    bpdata = np.genfromtxt(os.path.join(dataDir, 'Canon/', fname), delimiter=',',
                           dtype=zip(['wave','through'],[float]*2))
    bpTemp = Bandpass()
    bpTemp.setBandpass(bpdata['wave'], bpdata['through'])
    canonFilters[key] = bpTemp

print '----------'
print 'Cannon Filters'
for i,key in enumerate(cannonKeys):
    print key+'     %.2f    ' %(sm.computeMags(canonFilters[key])[0])
