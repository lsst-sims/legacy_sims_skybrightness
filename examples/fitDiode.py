import numpy as np
from scipy.optimize import curve_fit
import lsst.sims.skybrightness as sb
import matplotlib.pylab as plt
from lsst.sims.skybrightness.utils import robustRMS


def expPlusC(xdata, x0,x1, x3):
    """
    Let the sky flux be exponentially declining and hit some constant value
    """
    flux = x0*np.exp( (xdata)*x1)+x3
    return flux



def medFilt(x,y,bins=30):
    """
    make bins of x, and then median values of y and x in those bins
    """

    order = np.argsort(x)
    x = x[order]
    y = y[order]

    ack,binEdges = np.histogram(x,bins=bins)
    left = np.searchsorted(x,binEdges, side='left')
    right = np.searchsorted(x,binEdges,side='right')


    xbinned=[]
    ybinned=[]
    yrms = []
    for i in np.arange(left.size-1):
        xbinned.append(np.median(x[left[i]:left[i+1]]))
        ybinned.append(np.median(y[left[i]:left[i+1]]))
        yrms.append(robustRMS(y[left[i]:left[i+1]]))

    return np.array(xbinned), np.array(ybinned), np.array(yrms)


# Let's just fit the photodiode data and use the ESO model to get a zeropoint.

# Note this is an old npz I have lying around
data = np.load('/Users/yoachim/gitRepos/stash_skybrigtness/data/photodiode/photodiode.npz')
diode = data['photodiode'].copy()
data.close()

good = np.where((diode['sunAlt'] > np.radians(-23.)) & (diode['moonAlt'] < 0.)
                & (diode['r'] > 0) & (diode['z'] > 0) &(diode['y'] > 0) )
diode = diode[good]

# Arg, I think the z and y filters might have been mis-labeled and swapped at some point.
tempy = diode['z'].copy()
tempz = diode['y'].copy()

diode['y'] = tempy
diode['z'] = tempz


# Make some plots
filterNames = ['r','z','y']

# Manually set some limits sun Alt limits
altLimits = {'r':-9.5, 'z':-8.4,'y':-7.5}

fittedParams = {}

fig = plt.figure()
for i,filterName in enumerate(filterNames):
    good = np.where( diode[filterName] < 480. )
    ax = fig.add_subplot(3,1,i+1)
    ax.semilogy(np.degrees(diode['sunAlt'][good][::10]),
                diode[filterName][good][::10], 'ko', alpha=.01 )
    ax.set_xlabel('Sun Altitude (degrees)')
    ax.set_ylabel('Flux (arbitrary) ')
    ax.set_title(filterName)

    xbinned,ybinned,yrms = medFilt(np.degrees(diode['sunAlt'][good]),
                                   diode[filterName][good])
    ax.errorbar(xbinned,ybinned, yerr=yrms, fmt='yo')
    ax.axvline(x=altLimits[filterName], color='b', linestyle='--')

    goodBinned = np.where(xbinned < altLimits[filterName])
    fitParams, fitCovariances = curve_fit(expPlusC, np.radians(xbinned[goodBinned]),
                                          ybinned[goodBinned],
                                          sigma=yrms[goodBinned])
    fittedParams[filterName] = fitParams
    ax.plot(xbinned[goodBinned], expPlusC(np.radians(xbinned[goodBinned]), *fitParams), 'b' )

    # now need to compute what the zeropoint is expected to be


fig.tight_layout()
fig.savefig('Plots/diode.png')
plt.close(fig)
