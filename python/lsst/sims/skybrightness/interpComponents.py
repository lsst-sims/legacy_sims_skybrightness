import numpy as np
import os
import glob
import healpy as hp
from lsst.sims.photUtils import Sed,Bandpass
from lsst.sims.skybrightness.twilightFunc import twilightFunc
from scipy.interpolate import InterpolatedUnivariateSpline, interp1d, RegularGridInterpolator
import os
from lsst.utils import getPackageDir


def id2intid(ids):
    """
    take an array of ids, and convert them to an integer id.
    Handy if you want to put things into a sparse array.
    """
    uids = np.unique(ids)
    order = np.argsort(ids)
    oids = ids[order]
    uintids = np.arange(np.size(uids), dtype=int)
    left = np.searchsorted(oids , uids)
    right = np.searchsorted(oids,uids, side='right')
    intids = np.empty(ids.size, dtype=int)
    for i in range(np.size(left)): intids[left[i]:right[i]]=uintids[i]
    result = intids*0
    result[order] = intids
    return result, uids, uintids


def intid2id(intids, uintids, uids, dtype=int):
    """
    convert an int back to an id
    """
    ids = np.zeros(np.size(intids))

    order = np.argsort(intids)
    ointids = intids[order]
    left = np.searchsorted(ointids,uintids,side='left' )
    right = np.searchsorted(ointids,uintids,side='right' )
    for i,(le,ri) in enumerate(zip(left,right)): ids[le:ri]=uids[i]
    result = np.zeros(np.size(intids), dtype=dtype)
    result[order] = ids

    return result




class BaseSingleInterp(object):
    """
    Base class for sky components that only need to be interpolated on airmass
    """
    def __init__(self, compName=None, sortedOrder=['airmass','nightTimes'], mags=False):
        """
        mags: Rather than the full spectrum, return the LSST ugrizy magnitudes.
        """

        self.mags = mags

        dataDir = os.path.join(getPackageDir('sims_skybrightness_data'), 'ESO_Spectra/'+compName)

        filenames = glob.glob(dataDir+'/*.npz')
        if len(filenames) == 1:
            temp = np.load(filenames[0])
            self.wave = temp['wave'].copy()
            self.filterWave = temp['filterWave'].copy()
            self.spec = temp['spec'].copy()
        else:
            temp = np.load(filenames[0])
            self.wave = temp['wave'].copy()
            self.filterWave = temp['filterWave'].copy()
            self.spec = temp['spec'].copy()
            for filename in filenames[1:]:
                temp = np.load(filename)
                self.spec = np.append(self.spec, temp['spec'])
        # Take the log of the spectra in case we want to interp in log space.
        self.logSpec = np.log10(self.spec['spectra'])
        self.specSize = self.spec['spectra'][0].size

        # What order are the dimesions sorted by (from how the .npz was packaged)
        self.sortedOrder = sortedOrder
        self.dimDict = {}
        self.dimSizes = {}
        for dt in self.sortedOrder:
            self.dimDict[dt] = np.unique(self.spec[dt])
            self.dimSizes[dt] = np.size(np.unique(self.spec[dt]))

    def __call__(self, intepPoints):
        if self.mags:
            return self.interpMag(intepPoints)
        else:
            return self.interpSpec(intepPoints)


    def indxAndWeights(self, points, grid):
        """
        for given 1-D points, find the grid points on either side and return the weights
        assume grid is sorted
        """

        order = np.argsort(points)

        indxL = np.empty(points.size, dtype=int)
        indxR = np.empty(points.size, dtype=int)

        indxR[order] = np.searchsorted(grid, points[order])
        indxL = indxR-1

        fullRange = grid[indxR]-grid[indxL]
        wL = np.zeros(fullRange.size, dtype=float)
        wR = np.ones(fullRange.size, dtype=float)

        good = np.where(fullRange != 0)
        wL[good] = (grid[indxR][good] - points[good])/fullRange[good]
        wR[good] = (points[good] - grid[indxL[good]])/fullRange[good]

        return indxR,indxL,wR,wL


    def _weighting(self, interpPoints, values):
        """
        given a list/array of airmass values, return a dict with the interpolated
        spectrum at each airmass and the wavelength array.

        Input interpPoints should be sorted
        """

        results = np.zeros( (interpPoints.size, np.size(values[0])) ,dtype=float)

        inRange = np.where( (interpPoints['airmass'] <= np.max(self.dimDict['airmass'])) &
                            (interpPoints['airmass'] >= np.min(self.dimDict['airmass'])))
        indxR,indxL,wR,wL = self.indxAndWeights(interpPoints['airmass'][inRange],
                                                self.dimDict['airmass'])

        nextra = 3

        # XXX--should I use the log spectra?  Make a check and switch back and forth?
        results[inRange] = wR[:,np.newaxis]*values[indxR*nextra] + \
                           wL[:,np.newaxis]*values[indxL*nextra]

        return results


    def interpSpec(self, interpPoints):
        result = self._weighting(interpPoints, self.logSpec)
        mask = np.where(result == 0.)
        result = 10.**result
        result[mask]  = 0.
        return {'spec':result, 'wave':self.wave}

    def interpMag(self, interpPoints):
        result = self._weighting(interpPoints, self.spec['mags'])
        mask = np.where(result == 0.)
        result =  10.**(-0.4*(result-np.log10(3631.)))
        result[mask]  = 0.
        return {'spec':result, 'wave':self.filterWave}

class ScatteredStar(BaseSingleInterp):
    """
    Interpolate the spectra caused by scattered starlight.
    """
    def __init__(self, compName='ScatteredStarLight', mags=False):
        super(ScatteredStar,self).__init__(compName=compName, mags=mags)


class LowerAtm(BaseSingleInterp):
    """
    Interpolate the spectra caused by the lower atmosphere.
    """
    def __init__(self, compName='LowerAtm', mags=False):
        super(LowerAtm,self).__init__(compName=compName,mags=mags )

class UpperAtm(BaseSingleInterp):
    """
    Interpolate the spectra caused by the upper atmosphere.
    """
    def __init__(self, compName='UpperAtm', mags=False):
        super(UpperAtm,self).__init__(compName=compName, mags=mags)

class MergedSpec(BaseSingleInterp):
    """
    Interpolate the spectra caused by the sum of the scattered starlight, airglow, upper and lower atmosphere.
    """
    def __init__(self, compName='MergedSpec', mags=False):
        super(MergedSpec,self).__init__(compName=compName, mags=mags)


class Airglow(BaseSingleInterp):
    """
    Interpolate the spectra caused by airglow.
    """
    def __init__(self, compName='Airglow', sortedOrder=['airmass','solarFlux'], mags=False):
        super(Airglow,self).__init__(compName=compName, mags=mags, sortedOrder=sortedOrder)
        self.nSolarFlux = np.size(self.dimDict['solarFlux'])


    def _weighting(self, interpPoints, values):

        results = np.zeros( (interpPoints.size, np.size(values[0])) ,dtype=float)
        # Only interpolate point that lie in the model grid
        inRange = np.where( (interpPoints['airmass'] <= np.max(self.dimDict['airmass'])) &
                            (interpPoints['airmass'] >= np.min(self.dimDict['airmass'])) &
                            (interpPoints['solarFlux'] >= np.min(self.dimDict['solarFlux'])) &
                            (interpPoints['solarFlux'] <= np.max(self.dimDict['solarFlux'])) )
        usePoints = interpPoints[inRange]
        amRightIndex, amLeftIndex, amRightW, amLeftW = self.indxAndWeights(usePoints['airmass'],
                                                                           self.dimDict['airmass'])

        sfRightIndex, sfLeftIndex, sfRightW, sfLeftW = self.indxAndWeights(usePoints['solarFlux'],
                                                                           self.dimDict['solarFlux'])

        for amIndex, amW in zip([amRightIndex,amLeftIndex], [amRightW,amLeftW] ):
            for sfIndex, sfW in zip([sfRightIndex, sfLeftIndex],[sfRightW, sfLeftW] ):
                results[inRange] += amW[:,np.newaxis]*sfW[:,np.newaxis]*values[amIndex*self.nSolarFlux+sfIndex ]

        return results


class TwilightInterp(object):
    def __init__(self, mags=False,
                 darkSkyMags = {'u':22.8, 'g':22.3, 'r':21.2,
                                'i':20.3, 'z':19.3, 'y':18.0,
                                'B':22.35, 'G':21.71, 'R':21.3}):
        """
        Read the Solar spectrum into a handy object and compute mags in different filters
        mags:  If true, only return the LSST filter magnitudes, otherwise return the full spectrum

        darkSkyMags = dict of the zenith dark sky values to be assumed. The twilight fits are
        done relative to the dark sky level.
        """
        # XXX Note the darkSkyMags still should to be averaged over lots of zodiacal values.

        self.mags = mags

        dataDir = getPackageDir('sims_skybrightness_data')

        solarSaved = np.load(os.path.join(dataDir,'solarSpec/solarSpec.npz'))
        self.solarSpec = Sed(wavelen=solarSaved['wave'], flambda=solarSaved['spec'])
        solarSaved.close()

        canonFilters = {}
        fnames = ['blue_canon.csv', 'green_canon.csv','red_canon.csv']

        # Filter names, from bluest to reddest.
        self.filterNames = ['B', 'G', 'R']

        for fname,filterName in zip(fnames,self.filterNames) :
            bpdata = np.genfromtxt(os.path.join(dataDir, 'Canon/', fname), delimiter=',',
                                   dtype=zip(['wave','through'],[float]*2))
            bpTemp = Bandpass()
            bpTemp.setBandpass(bpdata['wave'], bpdata['through'])
            canonFilters[filterName] = bpTemp

        # Tack on the LSST filters
        throughPath =  os.path.join(getPackageDir('throughputs'), 'baseline')
        lsstKeys = ['u', 'g', 'r', 'i', 'z','y']
        for key in lsstKeys:
            bp = np.loadtxt(os.path.join(throughPath, 'filter_'+key+'.dat'),
                            dtype=zip(['wave','trans'],[float]*2 ))
            tempB = Bandpass()
            tempB.setBandpass(bp['wave'],bp['trans'])
            canonFilters[key] = tempB
            self.filterNames.append(key)


        # MAGIC NUMBERS from fitting the all-sky camera:
        # Code to generate values in fitTwiSlopesSimul.py
        # values are of the form:
        # 0: ratio of f^z_12 to f_dark^z
        # 1: slope of curve wrt sun alt
        # 2: airmass term (10^(arg[2]*(X-1)))
        # 3: azimuth term.
        # 4: zenith dark sky flux (erg/s/cm^2)

        # r, z, and y are based on fitting the zenith decay in:
        # fitDiode.py
        # Just assuming the shape parameter fits are similar to the other bands.
        # XXX-- I don't understand why R and r are so different. Or why z is so bright.
        self.fitResults = {'B': [8.42181815e+00,   2.29622121e+01,   2.85862729e-01,
                                 2.99902574e-01,   2.32325117e-04 ],
                           'G': [4.13747072e+00,   2.29416735e+01,   2.97229615e-01,
                                 3.15696066e-01,   4.20961686e-04],
                           'R': [2.73450774e+00,   2.22015053e+01,   2.97825187e-01,
                                 3.28865935e-01,   2.08470485e-04],
                           #'r': [ 0.52247301,  22.51393345, 0.3, 0.3,  54.8812249],
                           'z': [0.74072461,  23.37634241, 0.3, 0.3,  12.88718065],
                           'y': [0.13894689,  23.41098193, 0.3, 0.3,  29.46852266]}


        # XXX-completely arbitrary fudge factor to make things brighter in the blue
        # Just copy the blue and say it's brighter.
        self.fitResults['u'] = [16.,   2.29622121e+01,   2.85862729e-01,
                                2.99902574e-01,   2.32325117e-04 ]

        # Take out any filters that don't have fit results
        self.filterNames = [ key for key in self.filterNames if key in self.fitResults.keys() ]

        self.effWave = []
        self.solarMag = []
        for filterName in self.filterNames :
            self.effWave.append(canonFilters[filterName].calcEffWavelen()[0])
            self.solarMag.append(self.solarSpec.calcMag(canonFilters[filterName]))

        ord = np.argsort(self.effWave)
        self.filterNames = np.array(self.filterNames)[ord]
        self.effWave = np.array(self.effWave)[ord]
        self.solarMag = np.array(self.solarMag)[ord]

        # update the fit results to be zeropointed properly
        for key in self.fitResults:
            f0 = 10.**(-0.4*(darkSkyMags[key]-np.log10(3631.)))
            self.fitResults[key][-1] = f0

        self.solarWave = self.solarSpec.wavelen
        self.solarFlux = self.solarSpec.flambda
        # This one isn't as bad as the model grids, maybe we could get away with computing the magnitudes
        # in the __call__ each time.
        if mags:
            # Load up the LSST filters and convert the solarSpec.flabda and solarSpec.wavelen to fluxes
            throughPath = os.path.join(getPackageDir('throughputs'),'baseline')
            self.lsstFilterNames = ['u','g','r','i','z','y']
            self.lsstEquations = np.zeros((np.size(self.lsstFilterNames),
                                           np.size(self.fitResults['B'])), dtype=float)
            self.lsstEffWave = []

            fits = np.empty((np.size(self.effWave), np.size(self.fitResults['B'])), dtype=float)
            for i,fn in enumerate(self.filterNames):
                fits[i,:] = self.fitResults[fn]

            for filtername in self.lsstFilterNames:
                bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                                dtype=zip(['wave','trans'],[float]*2 ))
                tempB = Bandpass()
                tempB.setBandpass(bp['wave'],bp['trans'])
                self.lsstEffWave.append(tempB.calcEffWavelen()[0] )
            # loop through the parameters and interpolate to new eff wavelengths
            for i in np.arange(self.lsstEquations[0,:].size):
                interp = InterpolatedUnivariateSpline(self.effWave,fits[:,i]) #interp1d(self.effWave,fits[:,i])
                self.lsstEquations[:,i] = interp(self.lsstEffWave)
            # Set the dark sky flux
            for i,filterName in enumerate(self.lsstFilterNames):
                self.lsstEquations[i,-1] = 10.**(-0.4*(darkSkyMags[filterName]-np.log10(3631.)))


    def printFitsUsed(self):
        """
        Print out the fit parameters being used
        """
        print '\\tablehead{\colhead{Filter} & \colhead{$r_{12/z}$} & \colhead{$a$ (1/radians)} & \colhead{$b$ (1/airmass)} & \colhead{$c$ (az term/airmass)} & \colhead{$f_z_dark$ (erg/s/cm$^2$)$\\times 10^8$} & \colhead{m$_z_dark$}}'
        for key in self.fitResults.keys():
            numbers = ''
            for num in  self.fitResults[key]:
                if num > .001:
                    numbers += ' & %.2f' % num
                else:
                    numbers += ' & %.2f' % (num*1e8)
            print key, numbers, ' & ', '%.2f' % (-2.5*np.log10(self.fitResults[key][-1])+np.log10(3631.))



    def __call__(self, intepPoints):
        if self.mags:
            return self.interpMag(intepPoints)
        else:
            return self.interpSpec(intepPoints)



    def interpMag(self, interpPoints, maxAM=2.5,
                     limits=[np.radians(-11.), np.radians(-20.)]):
        npts = np.size(self.lsstEffWave)
        result = np.zeros((np.size(interpPoints), npts), dtype=float )

        good = np.where( (interpPoints['sunAlt'] >= np.min(limits)) &
                         (interpPoints['sunAlt'] <= np.max(limits)) &
                         (interpPoints['airmass'] <= maxAM) &
                         (interpPoints['airmass'] >= 1.) )[0]

        for i,filterName in enumerate(self.lsstFilterNames):
            result[good,i] = twilightFunc(interpPoints[good], *self.lsstEquations[i,:].tolist() )
        #mask = np.where(result == 0.)
        #result =  10.**(-0.4*(result-np.log10(3631.)))
        #result[mask]  = 0.
        return {'spec':result, 'wave':self.lsstEffWave}

    def interpSpec(self, interpPoints, maxAM=2.5,
                     limits=[np.radians(-11.), np.radians(-20.)]):
        """
        interpPoints should have airmass, azRelSun, and sunAlt.
        """

        npts = np.size(self.solarWave)
        result = np.zeros((np.size(interpPoints), npts), dtype=float )

        good = np.where( (interpPoints['sunAlt'] >= np.min(limits)) &
                         (interpPoints['sunAlt'] <= np.max(limits)) &
                         (interpPoints['airmass'] <= maxAM) &
                         (interpPoints['airmass'] >= 1.) )[0]

        # Compute the expected flux in each of the filters that we have fits for
        fluxes = []
        for filterName in self.filterNames:
            fluxes.append( twilightFunc(interpPoints[good],*self.fitResults[filterName]))
        fluxes = np.array(fluxes)

        # ratio of model flux to raw solar flux:
        yvals = fluxes.T/(10.**(-0.4*(self.solarMag-np.log10(3631.)) ))

        # Find wavelengths bluer than cutoff
        blueRegion = np.where(self.solarWave < np.min(self.effWave))


        for i,yval in enumerate(yvals):
            interpF = interp1d(self.effWave, yval, bounds_error=False, fill_value=yval[-1])
            ratio = interpF(self.solarWave)
            interpBlue = InterpolatedUnivariateSpline(self.effWave, yval, k=1)
            ratio[blueRegion] = interpBlue(self.solarWave[blueRegion])
            result[good[i]] = self.solarFlux*ratio

        return {'spec':result, 'wave':self.solarWave}



class MoonInterp(BaseSingleInterp):
    """
    Read in the saved Lunar spectra and interpolate.
    """
    def __init__(self, compName='Moon', sortedOrder=['moonSunSep','moonAltitude', 'hpid'], mags=False):
        super(MoonInterp,self).__init__(compName=compName, sortedOrder=sortedOrder, mags=mags)
        # Magic number from when the templates were generated
        self.nside = 4


    def _weighting(self, interpPoints, values):
        """

        """

        result = np.zeros( (interpPoints.size, np.size(values[0])) ,dtype=float)

        # Check that moonAltitude is in range, otherwise return zero array
        if np.max(interpPoints['moonAltitude']) < np.min(self.dimDict['moonAltitude']):
            return result

        # Find the neighboring healpixels
        hpids, hweights =  hp.get_neighbours(self.nside, np.pi/2.-interpPoints['alt'],
                                                interpPoints['azRelMoon'] )

        badhp = np.in1d(hpids.ravel(), self.dimDict['hpid'], invert=True).reshape(hpids.shape)
        hweights[badhp] = 0.

        norm = np.sum(hweights,axis=0)
        good= np.where(norm != 0.)[0]
        hweights[:,good] = hweights[:,good]/norm[good]

        # Find the neighboring moonAltitude points in the grid
        rightMAs, leftMAs, maRightW, maLeftW = self.indxAndWeights(interpPoints['moonAltitude'],
                                                                   self.dimDict['moonAltitude'])

        # Find the neighboring moonSunSep points in the grid
        rightMss, leftMss, mssRightW, mssLeftW = self.indxAndWeights(interpPoints['moonSunSep'],
                                                                   self.dimDict['moonSunSep'])

        nhpid = self.dimDict['hpid'].size
        nMA = self.dimDict['moonAltitude'].size
        # Convert the hpid to an index.
        tmp = intid2id(hpids.ravel(),  self.dimDict['hpid'],
                          np.arange( self.dimDict['hpid'].size))
        hpindx = tmp.reshape(hpids.shape)
        # loop though the hweights and the moonAltitude weights

        for hpid,hweight in zip(hpindx,hweights):
            for maid,maW in zip([rightMAs,leftMAs],[maRightW,maLeftW] ):
                for mssid,mssW in zip([rightMss,leftMss], [mssRightW,mssLeftW]):
                    weight = hweight*maW*mssW
                    result += weight[:,np.newaxis]*values[mssid*nhpid*nMA+maid*nhpid+hpid]

        return result


class ZodiacalInterp(BaseSingleInterp):
    """
    Interpolate the zodiacal light based on the airmass and the healpix ID where
    the healpixels are in ecliptic coordinates, with the sun at ecliptic longitude zero
    """

    def __init__(self, compName='Zodiacal', sortedOrder=['airmass', 'hpid'], mags=False):
        super(ZodiacalInterp,self).__init__(compName=compName, sortedOrder=sortedOrder, mags=mags)
        self.nside = hp.npix2nside(np.size(np.where(self.spec['airmass'] ==
                                                    np.unique(self.spec['airmass'])[0])[0]))


    def _weighting(self,interpPoints, values):
        """
        interpPoints is a numpy array where interpolation is desired
        values are the model values.
        """
        result = np.zeros( (interpPoints.size, np.size(values[0])) ,dtype=float)

        # Find the neighboring healpixels
        hpids, hweights =  hp.get_neighbours(self.nside, np.pi/2.-interpPoints['altEclip'],
                                                interpPoints['azEclipRelSun'] )

        badhp = np.in1d(hpids.ravel(), self.dimDict['hpid'], invert=True).reshape(hpids.shape)
        hweights[badhp] = 0.

        norm = np.sum(hweights,axis=0)
        good= np.where(norm != 0.)[0]
        hweights[:,good] = hweights[:,good]/norm[good]


        #norm = np.sum(hweights,axis=0)
        #hweights = hweights/norm

        # Find the neighboring airmass points in the grid
        order = np.argsort(interpPoints['airmass'])
        good = np.where( (interpPoints['airmass'][order] >= np.min( self.dimDict['airmass'])) &
                         (interpPoints['airmass'][order] <= np.max( self.dimDict['airmass']))  )
        rightAMs = np.searchsorted(self.dimDict['airmass'], interpPoints[order]['airmass'] )
        leftAMs = rightAMs-1

        # Set the indices that are out of the grid to 0.
        leftAMs[np.where(leftAMs) < 0] = 0
        rightAMs[np.where(rightAMs > self.dimDict['airmass'].size-1)] = 0
        amids = np.array([rightAMs,leftAMs] )

        amWs= np.zeros((2,interpPoints.size), dtype=float)
        amWs[0,order[good]] = (self.dimDict['airmass'][rightAMs[good]]-interpPoints['airmass'][order[good]])/(self.dimDict['airmass'][rightAMs[good]]- self.dimDict['airmass'][leftAMs[good]])
        amWs[1,order[good]] =(interpPoints['airmass'][order[good]]-self.dimDict['airmass'][leftAMs[good]])/(self.dimDict['airmass'][rightAMs[good]]- self.dimDict['airmass'][leftAMs[good]])

        nhpid = self.dimDict['hpid'].size
        # loop though the hweights and the airmass weights
        for hpid,hweight in zip(hpids,hweights):
            for amid,amW in zip(amids, amWs):
                weight = hweight*amW
                result += weight[:,np.newaxis]*values[amid*nhpid+hpid]

        return result
