from __future__ import print_function
import numpy as np
from scipy.optimize import curve_fit
import lsst.sims.skybrightness as sb
import matplotlib.pylab as plt
from lsst.sims.skybrightness.utils import robustRMS
import os
import lsst.sims.photUtils.Bandpass as Bandpass
from lsst.sims.skybrightness import zenithTwilight


# Let's just fit the photodiode data and use the ESO model to get a zeropoint.


def medFilt(x, y, bins=30):
    """
    make bins of x, and then median values of y and x in those bins
    """

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    ack, binEdges = np.histogram(x, bins=bins)
    left = np.searchsorted(x, binEdges, side='left')
    right = np.searchsorted(x, binEdges, side='right')

    xbinned = []
    ybinned = []
    yrms = []
    for i in np.arange(left.size-1):
        xbinned.append(np.median(x[left[i]:left[i+1]]))
        ybinned.append(np.median(y[left[i]:left[i+1]]))
        yrms.append(robustRMS(y[left[i]:left[i+1]]))

    return np.array(xbinned), np.array(ybinned), np.array(yrms)

# Note this is an old npz I have lying around
if 'diode' not in globals():
    print('Loading diode data')
    #data = np.load('/Users/yoachim/gitRepos/stash_skybrigtness/data/photodiode/photodiode.npz')
    data = np.load('/Users/yoachim/gitRepos/sims_skybrightness/data/photodiode/photodiode.npz')
    diode = data['photodiode'].copy()
    data.close()

    good = np.where((diode['sunAlt'] > np.radians(-23.)) & (diode['moonAlt'] < 0.)
                    & (diode['r'] > 0) & (diode['z'] > 0) & (diode['y'] > 0))
    diode = diode[good]

    # Arg, I think the z and y filters might have been mis-labeled and swapped at some point.
    tempy = diode['z'].copy()
    tempz = diode['y'].copy()

    diode['y'] = tempy
    diode['z'] = tempz

# Load up LSST filters
throughPath = os.getenv('LSST_THROUGHPUTS_BASELINE')
keys = ['r', 'z', 'y']
nfilt = len(keys)
filters = {}
for filtername in keys:
    bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                    dtype=zip(['wave', 'trans'], [float]*2))
    tempB = Bandpass()
    tempB.setBandpass(bp['wave'], bp['trans'])
    filters[filtername] = tempB


sm = sb.SkyModel(twilight=False, moon=False, zodiacal=False, mags=True)
sm.setRaDecMjd([0.], [np.radians(90.)], 45000., azAlt=True)
sm.computeSpec()
mags = sm.computeMags()
modelFluxLevels = {}
modelFluxLevels['r'] = sm.spec[0][2]
modelFluxLevels['z'] = sm.spec[0][4]
modelFluxLevels['y'] = sm.spec[0][5]

modelMagLevels = {}
modelMagLevels['r'] = mags[0][2]
modelMagLevels['z'] = mags[0][4]
modelMagLevels['y'] = mags[0][5]


# Make some plots
filterNames = ['r', 'z', 'y']

# Manually set some limits sun Alt limits
altLimits = {'r': -9.5, 'z': -8.4, 'y': -7.5}

fittedParams = {}

fig = plt.figure()
for i, filterName in enumerate(filterNames):
    good = np.where(diode[filterName] < 480.)
    ax = fig.add_subplot(3, 1, i+1)
    ax.plot(np.degrees(diode['sunAlt'][good][::100]),
            -2.5*np.log10(diode[filterName][good][::100]), 'ko', alpha=.01)
    ax.set_xlabel('Sun Altitude (degrees)')
    ax.set_ylabel('Mags')
    ax.set_title(filterName)

    xbinned, ybinned, yrms = medFilt(np.degrees(diode['sunAlt'][good]),
                                     diode[filterName][good])
    magErr = 1.0857*yrms/ybinned
    ax.errorbar(xbinned, -2.5*np.log10(ybinned), yerr=magErr, fmt='yo')
    ax.axvline(x=altLimits[filterName], color='b', linestyle='--')
    ax.set_xlim([-5, -24])
    ax.axvline(x=-12, color='g', linestyle='-')

    goodBinned = np.where(xbinned < altLimits[filterName])
    p0 = np.array([1., 20., -1, -1, 100.])
    fitParams, fitCovariances = curve_fit(zenithTwilight, np.radians(xbinned[goodBinned]),
                                          ybinned[goodBinned],
                                          sigma=yrms[goodBinned], p0=p0)

    ax.plot(xbinned[goodBinned], -2.5 *
            np.log10(zenithTwilight(np.radians(xbinned[goodBinned]), *fitParams)), 'b')

    ax.set_ylim([-1, -8])

    fittedParams[filterName] = fitParams


fig.tight_layout()
fig.savefig('Plots/diode.png')
fig.savefig('Plots/diode.pdf')

plt.close(fig)

print('Best fit exponential = constant, normalized to ESO dark zenith')
print(fittedParams)
