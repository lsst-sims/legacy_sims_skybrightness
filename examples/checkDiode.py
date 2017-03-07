from __future__ import print_function
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
    print('Loading diode data')
    #data = np.load('/Users/yoachim/gitRepos/stash_skybrigtness/data/photodiode/photodiode.npz')
    data = np.load('/Users/yoachim/gitRepos/sims_skybrightness/data/photodiode/photodiode.npz')
    diode = data['photodiode'].copy()
    data.close()

    # Arg, I think the z and y filters might have been mis-labeled and swapped at some point.
    tempy = diode['z'].copy()
    tempz = diode['y'].copy()

    diode['y'] = tempy
    diode['z'] = tempz


moonMin = -100.
good = np.where((diode['r'] > 0) & (diode['z'] > 0) &
                (diode['y'] > 0) &
                (diode['sunAlt'] < np.radians(-11.)) &
                (diode['r'] > 20) &
                (diode['moonPhase'] > moonMin))
diode = diode[good]

subDiode = diode[::1000]

gridsize = 20

filters = ['r', 'z', 'y']
pads = {'r': 2.5, 'z': 2, 'y': .5}
for j, filterName in enumerate(filters):
    pad = pads[filterName]
    fig, ax = plt.subplots(2)

    obsMags = -2.5*np.log10(subDiode[filterName])
    dark = np.where((subDiode['moonAlt'] < 0) & (subDiode['sunAlt'] < np.radians(-22.)))
    darkLevel = np.median(obsMags[dark])

    hb = ax[0].hexbin(np.degrees(subDiode['moonAlt']), np.degrees(subDiode['sunAlt']),
                      obsMags, vmin=darkLevel-pad, vmax=darkLevel,
                      gridsize=gridsize)
    cb = fig.colorbar(hb, ax=ax[0])
    ax[0].set_xlabel('Moon Alt (degrees)')
    ax[0].set_ylabel('Sun Alt (degrees)')
    ax[0].set_title('Photodiode, %s' % filterName)

    modelVals = np.zeros((subDiode.size, 3), dtype=float)

    sm = sb.SkyModel(mags=True)
    if not os.path.isfile('diode_%s.npz' % filterName):
        for i, mjd in enumerate(subDiode['mjd']):
            sm.setRaDecMjd(0., 89.9, mjd, degrees=True, azAlt=True)
            sm.computeSpec()
            mags = sm.computeMags()
            modelVals[i] = mags[0][np.array([2, 4, 5])]
        np.savez('diode_%s.npz' % filterName, modelVals=modelVals)
    else:
        ack = np.load('diode_%s.npz' % filterName)
        modelVals = ack['modelVals'].copy()

    darkLevel = np.median(modelVals[:, j][dark])
    hb = ax[1].hexbin(np.degrees(subDiode['moonAlt']), np.degrees(subDiode['sunAlt']), modelVals[:, j],
                      gridsize=gridsize, vmin=darkLevel-pad, vmax=darkLevel)
    cb = fig.colorbar(hb, ax=ax[1])
    ax[1].set_xlabel('Moon Alt (degrees)')
    ax[1].set_ylabel('Sun Alt (degrees)')
    ax[1].set_title('Sky Model, %s' % filterName)

    resid = -2.5*np.log10(subDiode[filterName]) - modelVals[:, j]
    resid -= np.median(resid)

#    hb = ax[2].hexbin(np.degrees(subDiode['moonAlt']),np.degrees(subDiode['sunAlt']),
#                      resid, vmin=-3, vmax=3,
#                      gridsize=gridsize )
#    cb = fig.colorbar(hb, ax=ax[2])
#    ax[2].set_xlabel('Moon Alt (degrees)')
#    ax[2].set_ylabel('Sun Alt (degrees)')
#    ax[2].set_title('Diode-Model, %s' % filterName)

    fig.tight_layout()
    fig.savefig('Plots/diodeCheck_%s.pdf' % filterName)
    plt.close(fig)
