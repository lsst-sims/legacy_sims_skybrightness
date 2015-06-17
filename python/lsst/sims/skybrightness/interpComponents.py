import numpy as np
import os
import glob
import healpy as hp
from lsst.sims.photUtils import Sed,Bandpass
from lsst.sims.skybrightness.twilightFunc import twilightFunc
from scipy.interpolate import interp1d
import os

class BaseSingleInterp(object):
    """
    Base class for sky components that only need to be interpolated on airmass
    """
    def __init__(self, compName=None, sortedOrder=['airmass','nightTimes'], mags=False):
        """
        mags: Rather than the full spectrum, return the LSST ugrizy magnitudes.
        """

        self.mags = mags

        dataDir =  os.path.join(os.environ.get('SIMS_SKYBRIGHTNESS_DATA_DIR'), 'ESO_Spectra/'+compName)

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

    def _weighting(self, interpPoints, values):
        """
        given a list/array of airmass values, return a dict with the interpolated
        spectrum at each airmass and the wavelength array.

        Input interpPoints should be sorted
        """

        order = np.argsort(interpPoints, order=self.sortedOrder)

        results = np.zeros( (interpPoints.size, np.size(values[0])) ,dtype=float)

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
        results[order[inRange]] = w1[inRange,np.newaxis]*values[left[inRange]*nextra] + \
                                  w2[inRange,np.newaxis]*values[right[inRange]*nextra]

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

class Airglow(BaseSingleInterp):
    """
    Interpolate the spectra caused by airglow.
    """
    def __init__(self, compName='Airglow', mags=False):
        super(Airglow,self).__init__(compName=compName, mags=mags)

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



class TwilightInterp(object):
    def __init__(self, mags=False,
                 darkSkyMags = {'u':22.8, 'g':22.3, 'r':21.2,
                                'i':20.3, 'z':19.3, 'y':18.0,
                                'B':22.35, 'G':21.71, 'R':21.3}):
        """
        Read the Solar spectrum into a handy object and compute mags in different filters

        darkSkyMags = dict of the zenith dark sky values to be assumed. The twilight fits are
        done relative to the dark sky level.
        """
        # XXX Note the darkSkyMags still need to be averaged over lots of zodiacal values.

        dataDir = os.getenv('SIMS_SKYBRIGHTNESS_DATA_DIR')

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

        # Tack on the LSST r z and y filter
        throughPath = os.getenv('LSST_THROUGHPUTS_BASELINE')
        lsstKeys = ['r', 'z','y']
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
        self.fitResults = {'B': [ 6.65202455e+00,   2.31066560e+01,   2.83706875e-01,
                                  3.01164450e-01,   3.02304747e-04 ],
                           'G': [4.05196324e+00,   2.27492146e+01,   3.02730529e-01,
                                 3.13501976e-01,   4.20257609e-04],
                           'R': [1.77783956e+00,   2.17505591e+01,   3.01964448e-01,
                                 3.33575403e-01,   2.98583706e-04],
                           #'r': [ 0.52247301,  22.51393345, 0.3, 0.3,  54.8812249],
                           'z': [0.74072461,  23.37634241, 0.3, 0.3,  12.88718065],
                           'y': [0.13894689,  23.41098193, 0.3, 0.3,  29.46852266]}


        # Take out any filters that don't have fit results
        self.filterNames = [ key for key in self.filterNames if key in self.fitResults.keys() ]

        self.effWave = []
        self.solarMag = []
        for filterName in self.filterNames :
            self.effWave.append(canonFilters[filterName].calcEffWavelen()[0])
            self.solarMag.append(self.solarSpec.calcMag(canonFilters[filterName]))

        self.solarMag = np.array(self.solarMag)

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
            throughPath = os.getenv('LSST_THROUGHPUTS_BASELINE')
            keys = ['u','g','r','i','z','y']
            newSolarWave = []
            newSolarFlux = []
            for filtername in keys:
                bp = np.loadtxt(os.path.join(throughPath, 'filter_'+filtername+'.dat'),
                                dtype=zip(['wave','trans'],[float]*2 ))
                tempB = Bandpass()
                tempB.setBandpass(bp['wave'],bp['trans'])
                newSolarWave.append(tempB.calcEffWavelen()[0])
                mag = self.solarSpec.calcMag(tempB)
                flux = 10.**(-0.4*(mag-np.log10(3631.)))
                newSolarFlux.append(flux)
            self.solarWave = np.array(newSolarWave)
            self.solarFlux = np.array(newSolarFlux)

    def __call__(self, interpPoints, maxAM=2.5,
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

        for i,yval in enumerate(yvals):
            interpF = interp1d(self.effWave, yval, bounds_error=False, fill_value=yval[-1])
            ratio = interpF(self.solarWave)
            ratio[np.where(self.solarWave < np.min(self.effWave))] = yval[0]
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
        A temporary method that does a stupid loop until I can figure out how to do the proper
        all numpy array slicing
        """

        result = np.zeros( (interpPoints.size, np.size(values[0])) ,dtype=float)

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
                                    result[i] += hweight*phasew*maw*values[good[0]]


        return result

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

    def __init__(self, compName='Zodiacal', sortedOrder=['airmass', 'hpid'], mags=False):
        super(ZodiacalInterp,self).__init__(compName=compName, sortedOrder=sortedOrder, mags=mags)
        self.nside = hp.npix2nside(np.size(np.where(self.spec['airmass'] ==
                                                    np.unique(self.spec['airmass'])[0])[0]))

    def _weighting(self, interpPoints, values):
        """
        Use some np.where mojo to find the templates that surround each interpolation
        point in parameter-space. Then, calculate a biliniear interpolation weight for each model spectrum.
        """

        result = np.zeros( (interpPoints.size, np.size(values[0])) ,dtype=float)

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
                            result[i] += hweight*amw*values[good[0]]

        return result

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
