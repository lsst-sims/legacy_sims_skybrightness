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


sm = sb.SkyModel( )
mjd = 56948.05174
sm.setRaDecMjd(np.array([0.]), np.array([90.]), mjd,azAlt=True, degrees=True)
sm.computeSpec()
print 'filter  ESO model   Overview Paper'
for i,key in enumerate(keys):
    print key+'     %.2f      %.2f ' %(sm.computeMags(filters[key])[0], overview_vals[i])

