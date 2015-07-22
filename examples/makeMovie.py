import numpy as np
import lsst.sims.skybrightness as sb
import lsst.sims.photUtils.Bandpass as Bandpass
import os
import matplotlib.pylab as plt
import healpy as hp
from lsst.sims.selfcal.analysis import healplots
from lsst.sims.maf.utils.telescopeInfo import TelescopeInfo
from lsst.sims.utils import altAzPaFromRaDec
import subprocess
import sys


def makeMovie(outfileroot, outDir='', ips=10.0, fps=10.0, figformat='png'):
    """
    call ffmpeg to stitch a movie together
    """

    callList = ['ffmpeg', '-r', str(ips), '-i',
                    os.path.join(outDir,'%s_%%04d.png'%(outfileroot)),
                    '-r', str(fps), '-pix_fmt', 'yuv420p', '-crf', '18', '-preset', 'slower',
                    os.path.join(outDir,'%s_%s_%s.mp4' %(outfileroot, str(ips), str(fps)))]
    print 'Attempting to call ffmpeg with:'
    print ' '.join(callList)
    p = subprocess.check_call(callList)

    # OK, let's do a gif too:
    callList = ['ffmpeg','-i',os.path.join(outDir,'%s_%%04d.png'%(outfileroot)),
                '-vf', 'scale=%s:%s' %(str(320),str(-1)), '-t', str(10), '-r', str(10),
                os.path.join(outDir,'%s_%s_%s.gif' %(outfileroot, str(ips), str(fps)))]
    print 'converting to animated gif with:'
    print ' '.join(callList)
    p2 = subprocess.check_call(callList)

telescope = TelescopeInfo('LSST')
sm = sb.SkyModel()

# Let's load up a pile of data and make a movie

canonDict = {}
canonFiles = {'R':'red_canon.csv','G':'green_canon.csv','B':'blue_canon.csv'}

path = os.path.join(os.environ.get('SIMS_SKYBRIGHTNESS_DATA_DIR'), 'Canon')
for key in canonFiles.keys():
    data = np.loadtxt(os.path.join(path,canonFiles[key]), delimiter=',',
                      dtype=zip(['wave','throughput'],[float,float]))
    band = Bandpass()
    band.setBandpass(data['wave'], data['throughput'])
    canonDict[key]=band


band = 'R'
nside = 8

# I can try 1665 to 2506
#select * from dates where sunAlt < -12/180*3.14 and mjd > 56948.258032-1 and mjd < 56948.258032+1 ;
#select * from dates where sunAlt < -12/180*3.14 and mjd > 56948  and mjd < 56961.106458+1 and moonPhase > 50 limit  100
# 8191. 8649
#dateIDs = np.arange(8191, 8694, 2)  #np.arange(1665, 2506, 5) #[2133,2134]
dateIDs = np.arange(28793, 29651, 2)
outDir = 'Movie'

cmin = 18.
cmax = 21.5

zp = 12.11

#zps = []
counter = 0
maxI = np.size(dateIDs)
for i,dateID in enumerate(dateIDs):

    skydata, mjd = sb.allSkyDB(dateID , filt=band)#2744
    airmass = 1./np.cos(np.radians(90.-skydata['alt']))
    skydata = skydata[np.where((airmass < 2.1) & (airmass >= 1.))]

    alt,az,pa =  altAzPaFromRaDec(np.radians(skydata['ra']), np.radians(skydata['dec']),
                                  telescope.lon, telescope.lat, mjd)

    skyhp = healplots.healbin(az, alt, skydata['sky'], nside=nside)

    sm.setRaDecMjd(np.radians(skydata['ra']), np.radians(skydata['dec']), mjd, degrees=False)
    if sm.sunAlt < np.radians(-12.):
        sm.computeSpec()
        mags = sm.computeMags(canonDict[band])

        good = np.where(mags > 0)
        if np.size(good[0]) > 10:
            modelhp = healplots.healbin(az[good], alt[good], mags[good], nside=nside)
            zp = np.median(mags[good] - skydata['sky'][good])
            #zps.append(zp)
            fig = plt.figure(num=1)
            hp.mollview(skyhp+zp, rot=(0,90), sub=(2,1,1), fig=1, title='Cannon '+band+' mjd=%0.2f'%sm.mjd,
                        min=cmin, max=cmax)
            hp.mollview(modelhp, rot=(0,90), sub=(2,1,2), fig=1,
                        title='Model. Sun Alt = %0.1f, moon alt = %0.1f' % (np.degrees(sm.sunAlt), np.degrees(sm.moonAlt)),
                        min=cmin, max=cmax)

            fig.savefig(os.path.join(outDir,'skymovie_%04d' % counter + '.png'))
            plt.close(1)
            counter += 1

    progress = i/float(maxI)*100
    text = "\rprogress = %.1f%%"%progress
    sys.stdout.write(text)
    sys.stdout.flush()


makeMovie('skymovie', outDir='Movie')
