import numpy as np
import lsst.sims.skybrightness as sb
import unittest
import lsst.sims.photUtils.Bandpass as Bandpass
from lsst.utils import getPackageDir
import os
import warnings

class TestSkyModel(unittest.TestCase):

    def testmergedComp(self):
        """
        Test that the 3 components that have been merged return the
        same result if they are computed independently
        """

        sky1 = sb.SkyModel(twilight=False, zodiacal=False,  moon=False,
                         lowerAtm=False, upperAtm=False,
                         airglow=False, scatteredStar=False,
                         mergedSpec=True)
        sky1.setRaDecMjd([36.],[-68.],49353.18, degrees=True)

        sky2 = sb.SkyModel(twilight=False, zodiacal=False,  moon=False,
                         lowerAtm=True, upperAtm=True,
                         airglow=False, scatteredStar=True,
                         mergedSpec=False)
        sky2.setRaDecMjd([36.],[-68.],49353.18, degrees=True)

        dummy, spec1 = sky1.returnWaveSpec()
        dummy, spec2 = sky2.returnWaveSpec()

        np.testing.assert_almost_equal(spec1, spec2)


    def testSetups(self):
        """
        Check that things are the same if the model is set up with
        radecmjd or all the parameters independently
        """

        sm1 = sb.SkyModel()
        sm1.setRaDecMjd([36.],[-68.],49353.18, degrees=True)

        sm2 = sb.SkyModel()
        sm2.setParams( azs=sm1.azs, alts=sm1.alts,
                       moonPhase=sm1.moonPhase,
                       moonAlt=sm1.moonAlt, moonAz=sm1.moonAz,
                       sunAlt=sm1.sunAlt, sunAz=sm1.sunAz,
                       sunEclipLon=sm1.sunEclipLon, eclipLon=sm1.eclipLon,
                       eclipLat=sm1.eclipLat, solarFlux=sm1.solarFlux,
                       degrees=False)

        dummy, spec1 = sm1.returnWaveSpec()
        dummy, spec2 = sm2.returnWaveSpec()

        np.testing.assert_array_equal(spec1, spec2)

        # Check that the degrees kwarg works
        sm2.setParams( azs=np.degrees(sm1.azs), alts=np.degrees(sm1.alts),
                       moonPhase=sm1.moonPhase,
                       moonAlt=np.degrees(sm1.moonAlt), moonAz=np.degrees(sm1.moonAz),
                       sunAlt=np.degrees(sm1.sunAlt), sunAz=np.degrees(sm1.sunAz),
                       sunEclipLon=np.degrees(sm1.sunEclipLon), eclipLon=np.degrees(sm1.eclipLon),
                       eclipLat=np.degrees(sm1.eclipLat), solarFlux=sm1.solarFlux,
                       degrees=True)

        atList = ['azs', 'alts', 'moonPhase', 'moonAlt', 'moonAz', 'sunAlt', 'sunAz',
                  'sunEclipLon', 'eclipLon', 'eclipLat', 'solarFlux']

        # Check each attribute that should match
        for attr in atList:
            np.testing.assert_allclose(getattr(sm1,attr), getattr(sm2,attr))

        # Check the interpolation points
        for name in sm1.points.dtype.names:
            np.testing.assert_allclose(sm1.points[name], sm2.points[name])

        # Check the final output spectra
        np.testing.assert_allclose(sm1.spec, sm2.spec)

    def testMags(self):
        """
        Test that the interpolated mags are similar to mags computed from interpolated spectra
        """

        throughPath = os.path.join(getPackageDir('throughputs'),'baseline')
        filters = ['u','g','r','i','z','y']

        bps = []
        # Load all the LSST bandpasses
        for filterName in filters:
            bp = np.loadtxt(os.path.join(throughPath, 'filter_%s.dat' % filterName),
                            dtype=zip(['wave','trans'],[float]*2 ))
            lsst_bp = Bandpass()
            lsst_bp.setBandpass(bp['wave'], bp['trans'])
            bps.append(lsst_bp)

        # Set up a sky model to interpolate spectra
        sm1 = sb.SkyModel()
        sm1.setRaDecMjd([36.],[-68.],49353.18, degrees=True)
        mags1 = []
        for bp in bps:
            mags1.append(sm1.returnMags(bandpass=bp))
        mags1 = np.array(mags1)

        # Set a sky model to interpolate pre-computed magnitudes directly
        sm2 = sb.SkyModel(mags=True)
        sm2.setRaDecMjd([36.],[-68.],49353.18, degrees=True)
        mag2 = sm2.returnMags()

        # Check that mags computed from spectra match those from interpolation
        rtol = 1e-4
        if np.max(np.abs(mags1-mag2.T)) > rtol:
            warnings.warn('Mags from spectra do not match pre-computed mags. \n Use the repo git@github.com:lsst-sims/sims_skybrightness_fits.git to update templates if throughputs or sims_photUtils have changed.')
        np.testing.assert_allclose(mags1,mag2.T, rtol=rtol)

    def test90Deg(self):
        """
        Make sure we can look all the way to 90 degree altitude.
        """
        mjd =  56973.268218 #56995.22103
        sm = sb.SkyModel(mags=True)
        sm.setRaDecMjd(0.,90.,mjd, degrees=True, azAlt=True)
        mags = sm.returnMags()
        assert(True not in np.isnan(mags))
        assert(True not in np.isnan(sm.spec))

    def testAirglow(self):
        """
        test that the airglow goes up with increasing SFU
        """

        mjd =  56973.268218 #56995.22103
        sm = sb.SkyModel(mags=True)
        sm.setRaDecMjd(0.,90.,mjd, degrees=True, azAlt=True, solarFlux=130.)
        magNormal = sm.returnMags()
        sm.setRaDecMjd(0.,90.,mjd, degrees=True, azAlt=True, solarFlux=50.)
        magFaint = sm.returnMags()
        sm.setRaDecMjd(0.,90.,mjd, degrees=True, azAlt=True, solarFlux=200.)
        magBright = sm.returnMags()

        assert(magNormal[0][2] < magFaint[0][2])
        assert(magNormal[0][2] > magBright[0][2])

    def test_setRaDecAltAzMjd(self):
        """
        Make sure sending in self-computed alt,az works
        """
        sm1 = sb.SkyModel(mags=True)
        sm2 = sb.SkyModel(mags=True)
        ra = np.array([0.,0.,0.])
        dec = np.array([-.1, -.2,-.3])
        mjd = 5900
        sm1.setRaDecMjd(ra,dec,mjd)
        m1 = sm1.returnMags()
        sm2.setRaDecAltAzMjd(ra,dec,sm1.alts,sm1.azs,mjd)
        m2 = sm1.returnMags()

        attrList = ['ra', 'dec', 'alts', 'azs']
        for attr in attrList:
            np.testing.assert_equal(getattr(sm1,attr), getattr(sm2,attr))

        np.testing.assert_allclose(m1,m2, rtol=1e-6)


if __name__=="__main__":
    unittest.main()
