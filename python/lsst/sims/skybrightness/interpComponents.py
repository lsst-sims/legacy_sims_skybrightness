import numpy as np
import os
import glob
import healpy as hp
from lsst.sims.photUtils import Sed,Bandpass
from lsst.sims.skybrightness.twilightFunc import twilightFunc
from scipy.interpolate import interp1d

class BaseSingleInterp(object):
    """
    Base class for sky components that only need to be interpolated on airmass
    """
    def __init__(self, compName=None, sortedOrder=['airmass','nightTimes']):

        dataDir =  os.path.join(os.environ.get('SIMS_SKYBRIGHTNESS_DATA'), 'ESO_Spectra/'+compName)
        filenames = glob.glob(dataDir+'/*.npz')
        if len(filenames) == 1:
            temp = np.load(filenames[0])
            self.wave = temp['wave'].copy()
            self.spec = temp['spec'].copy()
        else:
            temp = np.load(filenames[0])
            self.wave = temp['wave'].copy()
            self.spec = temp['spec'].copy()
            for filename in filenames[1:]:
                temp = np.load(filename)
                self.spec = np.append(self.spec, temp['spec'])

        # Take the log of the spectra in case we want to interp in log space.
        self.logSpec = np.log10(self.spec['spectra'])
        # What order are the dimesions sorted by (from how the .npz was packaged)
        self.sortedOrder = sortedOrder
        self.dimDict = {}
        self.dimSizes = {}
        for dt in self.sortedOrder:
            self.dimDict[dt] = np.unique(self.spec[dt])
            self.dimSizes[dt] = np.size(np.unique(self.spec[dt]))

    def __call__(self, interpPoints):
        """
        given a list/array of airmass values, return a dict with the interpolated
        spectrum at each airmass and the wavelength array.

        Input interpPoints should be sorted
        """

        order = np.argsort(interpPoints, order=self.sortedOrder)

        results = np.zeros( (interpPoints.size, self.spec['spectra'][0].size) ,dtype=float)

        # The model values for the left and right side.
        right = np.searchsorted(self.dimDict['airmass'], interpPoints['airmass'][order])
        left = right-1

        # catch it somewhere if the interp point is outside the model range?
        #inRange = np.where((left >= 0) & (right <= self.dimDict['airmass'].size)  & (left < right) )
        inRange = np.where( (interpPoints['airmass'][order] <= np.max(self.dimDict['airmass'])) &
                            (interpPoints['airmass'][order] >= np.min(self.dimDict['airmass'])))

        left[np.where(left < 0)] = 0
        right[np.where(right >= self.dimDict['airmass'].size)] = self.dimDict['airmass']-1

        # Calc the weights associated with each of those corners
        fullRange = self.dimDict['airmass'][right]-self.dimDict['airmass'][left]
        w1 = (self.dimDict['airmass'][right] - interpPoints['airmass'][order])/fullRange
        w2 = (interpPoints['airmass'][order] - self.dimDict['airmass'][left])/fullRange

        # Catch points that land on a model point
        onPoint = np.where(fullRange == 0)
        w1[onPoint] = 1.
        w2[onPoint] = 0.

        # Little kludge to make up for the fact that each airmass
        # has 3 "time of night" values that we're ignoring.
        nextra = 3

        # XXX--should I use the log spectra?  Make a check and switch back and forth?
        results[order[inRange]] = w1[inRange,np.newaxis]*self.logSpec[left[inRange]*nextra] + \
                                  w2[inRange,np.newaxis]*self.logSpec[right[inRange]*nextra]
        results[order[inRange]] = 10.**results[order[inRange]]
        #results[order] = results
        return {'spec':results, 'wave':self.wave}


class ScatteredStar(BaseSingleInterp):
    """
    Interpolate the spectra caused by scattered starlight.
    """
    def __init__(self, compName='ScatteredStarLight'):
        super(ScatteredStar,self).__init__(compName=compName)

class Airglow(BaseSingleInterp):
    """
    Interpolate the spectra caused by airglow.
    """
    def __init__(self, compName='Airglow'):
        super(Airglow,self).__init__(compName=compName)

class LowerAtm(BaseSingleInterp):
    """
    Interpolate the spectra caused by the lower atmosphere.
    """
    def __init__(self, compName='LowerAtm'):
        super(LowerAtm,self).__init__(compName=compName)

class UpperAtm(BaseSingleInterp):
    """
    Interpolate the spectra caused by the upper atmosphere.
    """
    def __init__(self, compName='UpperAtm'):
        super(UpperAtm,self).__init__(compName=compName)

class MergedSpec(BaseSingleInterp):
    """
    Interpolate the spectra caused by the sum of the scattered starlight, airglow, upper and lower atmosphere.
    """
    def __init__(self, compName='MergedSpec'):
        super(MergedSpec,self).__init__(compName=compName)



class TwilightInterp(object):
    def __init__(self):
        """
        Read the Solar spectrum into a handy object and compute mags in different filters
        """
        dataDir = os.getenv('SIMS_SKYBRIGHTNESS_DATA')

        solarSaved = np.load(os.path.join(dataDir,'solarSpec/solarSpec.npz'))
        self.solarSpec = Sed(wavelen=solarSaved['wave'], flambda=solarSaved['spec'])
        solarSaved.close()

        canonFilters = []
        fnames = ['blue_canon.csv', 'green_canon.csv','red_canon.csv']

        # Filter names, from bluest to reddest.
        self.filterNames = ['B', 'G', 'R']

        for fname in fnames:
            bpdata = np.genfromtxt(os.path.join(dataDir, 'Canon/', fname), delimiter=',',
                                   dtype=zip(['wave','through'],[float]*2))
            bpTemp = Bandpass()
            bpTemp.setBandpass(bpdata['wave'], bpdata['through'])
            canonFilters.append(bpTemp)
        self.effWave = []
        self.solarMag = []
        for cfilter in canonFilters:
            self.effWave.append(cfilter.calcEffWavelen()[0])
            self.solarMag.append(self.solarSpec.calcMag(cfilter))

        self.solarMag = np.array(self.solarMag)

        # MAGIC NUMBERS from fitting the all-sky camera:
        # Code to make in fitTwiSlopesSimul.py
        self.fitResults = {'B': [ 53.20504282,4.85950702e-04,-0.65325829,-1.,-0.69345613],
                           'G': [52.38200428,4.18033020e-04,-0.69706203,-1.,-0.72186434],
                           'R': [50.08252239,2.01774730e-04,-0.6953,-1.,-0.76808688]}


    def __call__(self, interpPoints, maxAM=2.5,
                     limits=[np.radians(-11.), np.radians(-20.)]):
        """
        interpPoints should have airmass, azRelSun, and sunAlt.
        """

        npts = np.size(self.solarSpec.wavelen)
        result = np.zeros((np.size(interpPoints), npts), dtype=float )

        good = np.where( (interpPoints['sunAlt'] >= np.min(limits)) &
                         (interpPoints['sunAlt'] <= np.max(limits)) &
                         (interpPoints['airmass'] <= maxAM) &
                         (interpPoints['airmass'] >= 1.) )[0]

        fluxes = []
        for filterName in self.filterNames:
            fluxes.append( twilightFunc(interpPoints[good],*self.fitResults[filterName]))
        fluxes = np.array(fluxes)

        yvals = 10.**(np.log10(fluxes.T)+self.solarMag/2.5)

        for i,yval in enumerate(yvals):
            interpF = interp1d(self.effWave, yval, bounds_error=False, fill_value=yval[-1])
            ratio = interpF(self.solarSpec.wavelen)
            ratio[np.where(self.solarSpec.wavelen < np.min(self.effWave))] = yval[0]
            result[good[i]] = self.solarSpec.flambda*ratio

        return {'spec':result, 'wave':self.solarSpec.wavelen}



class MoonInterp(BaseSingleInterp):
    """
    Read in the saved Lunar spectra and interpolate.
    """
    def __init__(self, compName='Moon', sortedOrder=['moonSunSep','moonAltitude', 'hpid']):
        super(MoonInterp,self).__init__(compName=compName, sortedOrder=sortedOrder)
        # Magic number from when the templates were generated
        self.nside = 4



    def __call__(self, interpPoints):
        """
        A temporary method that does a stupid loop until I can figure out how to do the proper
        all numpy array slicing
        """

        result = np.zeros( (interpPoints.size, self.spec['spectra'][0].size) ,dtype=float)

        for i,point in enumerate(interpPoints):
            hpids, hweights = hp.get_neighbours(self.nside, np.pi/2.-point['alt'],
                                            point['azRelMoon'] )
            badhp = np.in1d(hpids, self.dimDict['hpid'], invert=True)
            hweights[badhp] = 0.
            norm = np.sum(hweights,axis=0)
            if norm != 0:
                hweights = hweights/norm

                # Find the phase points
                upperPhase = self.spec['moonSunSep'][np.where(self.spec['moonSunSep'] >=
                                                              point['moonSunSep'])].min()
                lowerPhase = self.spec['moonSunSep'][np.where(self.spec['moonSunSep'] <=
                                                              point['moonSunSep'])].max()
                if upperPhase == lowerPhase:
                    phases = [upperPhase]
                    phaseWeights = [1.]
                else:
                    phases = [upperPhase, lowerPhase]
                    phaseWeights = np.abs(point['moonSunSep']-np.array(phases))/(upperPhase-lowerPhase)


                upperMoonAlt = self.spec['moonAltitude'][np.where(self.spec['moonAltitude'] >=
                                                                  point['moonAltitude'])]
                lowerMoonAlt = self.spec['moonAltitude'][np.where(self.spec['moonAltitude'] <=
                                                                  point['moonAltitude'])]
                if (np.size(upperMoonAlt) == 0) | (np.size(lowerMoonAlt) == 0):
                    pass
                else:
                    upperMoonAlt = upperMoonAlt.min()
                    lowerMoonAlt = lowerMoonAlt.max()

                    if upperMoonAlt == lowerMoonAlt:
                        moonAlts = [upperMoonAlt]
                        moonAltWeights = [1.]
                    else:
                        moonAlts = [upperMoonAlt,lowerMoonAlt]
                        moonAltWeights = np.abs(point['moonAltitude']-np.array(moonAlts))/(upperMoonAlt-lowerMoonAlt)

                    for hpid,hweight in zip(hpids, hweights):
                        for phase,phasew in zip(phases,phaseWeights):
                            for moonAlt,maw in zip(moonAlts,moonAltWeights):
                                good = np.where( (self.spec['moonSunSep'] == phase)  &
                                                 (self.spec['moonAltitude'] == moonAlt) &
                                                 (self.spec['hpid'] == hpid))[0]
                                if np.size(good) > 0:
                                    result[i] += hweight*phasew*maw*self.logSpec[good[0]]
                    result[i] = 10.**result[i]


        return {'spec':result, 'wave':self.wave}


    def New__call__(self, interpPoints):
        """
        interpPoints:  numpy array with dtypes of 'moonSunSep','moonAltitude', 'alt', and 'azRelMoon' where the
        azimuth is the azimuth relative to the moon (with a range 0-pi).
        """



        # maybe I should just assume the moonSunSep and moonAlt are all the same?  That could make things a
        # bit easier!

        order = np.argsort(interpPoints, order=self.sortedOrder[:-1])


        hpids, hweights = hp.get_neighbours(self.nside, np.pi/2.-interpPoints['alt'],
                                            interpPoints['azRelMoon'] )

        # Mask any neighbours that are not in the templates
        badhp = np.in1d(hpids, self.dimDict['hpid'], invert=True)
        badhp = badhp.reshape(hpids.shape)

        hweights[badhp] = 0.
        # Renormalize
        norm = np.sum(hweights,axis=0)
        good = np.where(norm > 0.)

        if norm.size == 1:
            norm = np.array([norm])

        hweights[:,good] = hweights[:,good]/norm[good]

        # Need to convert hpid to an index
        origShape = hpids.shape
        hpids = hpids.ravel()
        hporder = np.argsort(hpids)

        hpind = np.searchsorted(self.dimDict['hpid'],hpids[hporder] )
        hpind[hporder] = hpind
        hpind = hpind.reshape(origShape)

        result = np.zeros( (interpPoints.size, self.spec['spectra'][0].size) ,dtype=float)

        inRange = np.where( (norm > 0 ) &
                            (interpPoints['moonSunSep'][order] <= np.max(self.dimDict['moonSunSep'])) &
                            (interpPoints['moonSunSep'][order] >= np.min(self.dimDict['moonSunSep'])) &
                            (interpPoints['moonAltitude'][order] <= np.max(self.dimDict['moonAltitude'])) &
                            (interpPoints['moonAltitude'][order] >= np.min(self.dimDict['moonAltitude'])) )


        # Compute moonSunSep weights
        sepleft = np.searchsorted(self.dimDict['moonSunSep'], interpPoints['moonSunSep'][order])-1
        sepright = np.searchsorted(self.dimDict['moonSunSep'], interpPoints['moonSunSep'][order])
        fullRange = self.dimDict['moonSunSep'][sepright]-self.dimDict['moonSunSep'][sepleft]
        sepWleft = (self.dimDict['moonSunSep'][sepright] - interpPoints['moonSunSep'][order])/fullRange
        sepWright = (interpPoints['moonSunSep'][order] - self.dimDict['moonSunSep'][sepleft])/fullRange


        # Compute moonAltitude weights
        altleft = np.searchsorted(self.dimDict['moonAltitude'], interpPoints['moonAltitude'][order])-1
        altright = np.searchsorted(self.dimDict['moonAltitude'], interpPoints['moonAltitude'][order])
        fullRange = self.dimDict['moonAltitude'][altright]-self.dimDict['moonAltitude'][altleft]
        altWleft = (self.dimDict['moonAltitude'][altright] - interpPoints['moonAltitude'][order])/fullRange
        altWright = (interpPoints['moonAltitude'][order] - self.dimDict['moonAltitude'][altleft])/fullRange



        w1 = hweights[:,inRange[0]]*sepWleft[inRange] * altWleft[inRange]
        w2 = hweights[:,inRange[0]]*sepWright[inRange] * altWleft[inRange]
        w3 = hweights[:,inRange[0]]*sepWleft[inRange] * altWright[inRange]
        w4 = hweights[:,inRange[0]]*sepWright[inRange] * altWright[inRange]

        if interpPoints.size == 1:
            w1 = w1[:,np.newaxis]
            w2 = w2[:,np.newaxis]
            w3 = w3[:,np.newaxis]
            w4 = w4[:,np.newaxis]
        else:
            w1 = w1[:,:,np.newaxis]
            w2 = w2[:,:,np.newaxis]
            w3 = w3[:,:,np.newaxis]
            w4 = w4[:,:,np.newaxis]


        import pdb ; pdb.set_trace()

        result[order[inRange]] += np.sum(w1*self.spec['spectra'][hpind[:,inRange[0]]+altleft[inRange]*self.dimSizes['hpid']+
                                                 sepleft[inRange]*self.dimSizes['hpid']*self.dimSizes['moonAltitude']],
                         axis=0 )
        result[order[inRange]] += np.sum(w2*self.spec['spectra'][hpind[:,inRange[0]]+altleft[inRange]*self.dimSizes['hpid']+
                                                 sepright[inRange]*self.dimSizes['hpid']*self.dimSizes['moonAltitude']],
                         axis=0 )
        result[order[inRange]] += np.sum(w3*self.spec['spectra'][hpind[:,inRange[0]]+altright[inRange]*self.dimSizes['hpid']+
                                                 sepleft[inRange]*self.dimSizes['hpid']*self.dimSizes['moonAltitude']],
                         axis=0 )
        result[order[inRange]] += np.sum(w4*self.spec['spectra'][hpind[:,inRange[0]]+altright[inRange]*self.dimSizes['hpid']+
                                                 sepright[inRange]*self.dimSizes['hpid']*self.dimSizes['moonAltitude']],
                         axis=0 )

        return {'spec':result, 'wave':self.wave}



class ZodiacalInterp(BaseSingleInterp):
    """
    Interpolate the zodiacal light based on the airmass and the healpix ID where
    the healpixels are in ecliptic coordinates, with the sun at ecliptic longitude zero
    """

    def __init__(self, compName='Zodiacal', sortedOrder=['airmass', 'hpid']):
        super(ZodiacalInterp,self).__init__(compName=compName, sortedOrder=sortedOrder)
        self.nside = hp.npix2nside(np.size(np.where(self.spec['airmass'] ==
                                                    np.unique(self.spec['airmass'])[0])[0]))

    def __call__(self, interpPoints):
        """
        Use some np.where mojo to find the templates that surround each interpolation
        point in parameter-space. Then, calculate a biliniear interpolation weight for each model spectrum.
        """

        result = np.zeros( (interpPoints.size, self.spec['spectra'][0].size) ,dtype=float)

        for i,point in enumerate(interpPoints):

            hpids, hweights = hp.get_neighbours(self.nside, np.pi/2.-point['altEclip'],
                                                point['azEclipRelSun'] )

            badhp = np.in1d(hpids, self.dimDict['hpid'], invert=True)
            hweights[badhp] = 0.
            norm = np.sum(hweights,axis=0)
            if (norm != 0) & (point['airmass'] <= self.spec['airmass'].max()) & (point['airmass'] >= self.spec['airmass'].min()) :
                for hpid,hweight in zip(hpids,hweights):
                    hweights = hweights/norm

                # Find the airmass points

                upperAM = self.spec['airmass'][np.where(self.spec['airmass'] >= point['airmass'])].min()
                lowerAM = self.spec['airmass'][np.where(self.spec['airmass'] <= point['airmass'])].max()
                if upperAM == lowerAM:
                    airmasses = [upperAM]
                    amWeights = [1.]
                else:
                    airmasses = [upperAM, lowerAM]
                    amWeights = np.abs(point['airmass']-np.array(airmasses))/(upperAM-lowerAM)

                for hpid,hweight in zip(hpids, hweights):
                    for airmass,amw in zip(airmasses,amWeights):
                        good = np.where( (self.spec['airmass'] == airmass)  &
                                         (self.spec['hpid'] == hpid))[0]
                        if np.size(good) > 0:
                            result[i] += hweight*amw*self.logSpec[good[0]]
                result[i] = 10.**result[i]
        return {'spec':result, 'wave':self.wave}

    def new__call__(self, interpPoints):
        """
        interpPoints should be a numpy array with dtypes of alt, az, and airmass.
        Note the alt should be ecliptic declination, and az should be heliocentric ecliptic azimuth.
        """

        if interpPoints.size > 1:
            order = np.argsort(interpPoints, order=self.sortedOrder[:-1])
            interpPoints = interpPoints[order]

        #should use np.in1d I think for the healpixels
        hpids, hweights = hp.get_neighbours(self.nside, np.pi/2.-interpPoints['altEclip'],
                                            interpPoints['azEclipRelSun'] )
        # Mask any neighbours that are not in the templates
        badhp = np.in1d(hpids, self.dimDict['hpid'], invert=True)
        hpids[badhp] = 0
        hweights[badhp] = 0.
        # Renormalize
        norm = np.sum(hweights,axis=0)
        good = np.where(norm > 0.)

        if norm.size == 1:
            norm = np.array([norm])

        hweights[:,good] = hweights[:,good]/norm[good]

        # Need to convert hpid to an index
        origShape = hpids.shape
        hpids = hpids.ravel()
        hporder = np.argsort(hpids)

        hpind = np.searchsorted(self.dimDict['hpid'],hpids[hporder] )
        hpind[hporder] = hpind
        hpind = hpind.reshape(origShape)

        result = np.zeros( (interpPoints.size, self.spec['spectra'][0].size) ,dtype=float)

        inRange = np.where( (interpPoints['airmass'] >= self.dimDict['airmass'].min()) &
                            (interpPoints['airmass'] <= self.dimDict['airmass'].max() ))

        right = np.searchsorted(self.dimDict['airmass'], interpPoints['airmass'][inRange])
        left = right-1

        # Compute the weights based on how close the airmass is
        fullRange = self.dimDict['airmass'][right]-self.dimDict['airmass'][left]
        wleft = (self.dimDict['airmass'][right] - interpPoints['airmass'][inRange])/fullRange
        wright = (interpPoints['airmass'][inRange] - self.dimDict['airmass'][left])/fullRange

        if np.size(wleft) == 1:
            wleft = np.array([wleft])
            wright = np.array([wright])

        w1 = hweights[:,inRange[0]]*wleft
        w2 = hweights[:,inRange[0]]*wright

        if (np.size(wleft.shape) == 1) & (wleft.shape[0] == 1):
            w1 = w1[:,np.newaxis]
            w2 = w2[:,np.newaxis]
        else:
            w1 = w1[:,:,np.newaxis]
            w2 = w2[:,:,np.newaxis]

        result[inRange] += np.sum(w1*self.spec['spectra'][hpind[:,inRange[0]]+left*self.dimSizes['hpid']],
                         axis=0 )
        result[inRange] += np.sum(w2*self.spec['spectra'][hpind[:,inRange[0]]+right*self.dimSizes['hpid']],
                         axis=0 )

        if interpPoints.size > 1:
            result[order] = result
        return {'spec':result, 'wave':self.wave}
