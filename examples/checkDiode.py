import numpy as np
from scipy.optimize import curve_fit
import lsst.sims.skybrightness as sb
import matplotlib.pylab as plt
from lsst.sims.skybrightness.utils import robustRMS
import os
import lsst.sims.photUtils.Bandpass as Bandpass
from lsst.sims.skybrightness import zenithTwilight


# I don't know why the diode is so totally faint at 12-degree twilight, let's look at other situations


if 'diode' not in globals():
    print 'Loading diode data'
    #data = np.load('/Users/yoachim/gitRepos/stash_skybrigtness/data/photodiode/photodiode.npz')
    data = np.load('/Users/yoachim/gitRepos/sims_skybrightness/data/photodiode/photodiode.npz')
    diode = data['photodiode'].copy()
    data.close()

    # Arg, I think the z and y filters might have been mis-labeled and swapped at some point.
    tempy = diode['z'].copy()
    tempz = diode['y'].copy()

    diode['y'] = tempy
    diode['z'] = tempz


good = np.where( (diode['r'] > 0) & (diode['z'] > 0) &
                 (diode['y'] > 0) &
                 (diode['sunAlt'] < np.radians(-11.)) &
                 (diode['r'] > 20))
diode = diode[good]

subDiode=diode[::1000]


filters = ['r','z','y']
for j,filterName in enumerate(filters):
    fig,ax = plt.subplots(3)

    hb = ax[0].hexbin(np.degrees(subDiode['moonAlt']),np.degrees(subDiode['sunAlt']),
                      -2.5*np.log10(subDiode[filterName]) )
    cb = fig.colorbar(hb, ax=ax[0])
    ax[0].set_xlabel('Moon Alt (degrees)')
    ax[0].set_ylabel('Sun Alt (degrees)')
    ax[0].set_title('Photodiode, %s' % filterName)

    modelVals = np.zeros((subDiode.size,3), dtype=float)

    sm = sb.SkyModel(mags=True)
    for i,mjd in enumerate(subDiode['mjd']):
        sm.setRaDecMjd(0., 89.9, mjd, degrees=True, azAlt=True)
        sm.computeSpec()
        mags = sm.computeMags()
        modelVals[i] = mags[0][np.array([2,4,5])]

    hb = ax[1].hexbin(np.degrees(subDiode['moonAlt']),np.degrees(subDiode['sunAlt']), modelVals[:,j] )
    cb = fig.colorbar(hb, ax=ax[1])
    ax[1].set_xlabel('Moon Alt (degrees)')
    ax[1].set_ylabel('Sun Alt (degrees)')
    ax[1].set_title('Sky Model, %s' % filterName)


    resid = -2.5*np.log10(subDiode[filterName]) - modelVals[:,j]
    resid -= np.median(resid)

    hb = ax[2].hexbin(np.degrees(subDiode['moonAlt']),np.degrees(subDiode['sunAlt']), resid )
    cb = fig.colorbar(hb, ax=ax[2])
    ax[2].set_xlabel('Moon Alt (degrees)')
    ax[2].set_ylabel('Sun Alt (degrees)')
    ax[2].set_title('Diode-Model, %s' % filterName)

    fig.tight_layout()
    fig.savefig('Plots/diodeCheck_%s.pdf' % filterName)
    plt.close(fig)
