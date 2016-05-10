import numpy as np
import lsst.sims.skybrightness as sb
import unittest
import lsst.sims.photUtils.Bandpass as Bandpass
from lsst.utils import getPackageDir
import os


class TestSkyModel(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        # Load up the spectra just once to speed things up a bit
        self.sm_mags = sb.SkyModel(mags=True)
        self.sm_mags2 = sb.SkyModel(mags=True)
        self.sm_spec = sb.SkyModel(mags=False)
        self.sm_spec2 = sb.SkyModel(mags=False)

    def testmergedComp(self):
        """
        Test that the 3 components that have been merged return the
        same result if they are computed independently
        """

        sky1 = sb.SkyModel(twilight=False, zodiacal=False, moon=False,
                           lowerAtm=False, upperAtm=False,
                           airglow=False, scatteredStar=False,
                           mergedSpec=True)
        sky1.setRaDecMjd([36.], [-68.], 49353.18, degrees=True)

        sky2 = sb.SkyModel(twilight=False, zodiacal=False, moon=False,
                           lowerAtm=True, upperAtm=True,
                           airglow=False, scatteredStar=True,
                           mergedSpec=False)
        sky2.setRaDecMjd([36.], [-68.], 49353.18, degrees=True)

        dummy, spec1 = sky1.returnWaveSpec()
        dummy, spec2 = sky2.returnWaveSpec()

        np.testing.assert_almost_equal(spec1, spec2)

    def testSetups(self):
        """
        Check that things are the same if the model is set up with
        radecmjd or all the parameters independently
        """

        sm1 = self.sm_spec
        sm1.setRaDecMjd([36.], [-68.], 49353.18, degrees=True)

        sm2 = self.sm_spec2
        sm2.setParams(azs=sm1.azs, alts=sm1.alts,
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
        sm2.setParams(azs=np.degrees(sm1.azs), alts=np.degrees(sm1.alts),
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
            np.testing.assert_allclose(getattr(sm1, attr), getattr(sm2, attr))

        # Check the interpolation points
        for name in sm1.points.dtype.names:
            np.testing.assert_allclose(sm1.points[name], sm2.points[name])

        # Check the final output spectra
        np.testing.assert_allclose(sm1.spec, sm2.spec)

    def testMags(self):
        """
        Test that the interpolated mags are similar to mags computed from interpolated spectra
        """

        throughPath = os.path.join(getPackageDir('throughputs'), 'baseline')
        filters = ['u', 'g', 'r', 'i', 'z', 'y']

        bps = []
        for filterName in filters:
            bp = np.loadtxt(os.path.join(throughPath, 'filter_%s.dat' % filterName),
                            dtype=zip(['wave', 'trans'], [float]*2))
            lsst_bp = Bandpass()
            lsst_bp.setBandpass(bp['wave'], bp['trans'])
            bps.append(lsst_bp)

        sm1 = self.sm_spec
        sm1.setRaDecMjd([36.], [-68.], 49353.18, degrees=True)
        mags1 = []
        for bp in bps:
            mags1.append(sm1.returnMags(bandpass=bp))
        mags1 = np.array(mags1)

        sm2 = self.sm_mags
        sm2.setRaDecMjd([36.], [-68.], 49353.18, degrees=True)
        mag2 = sm2.returnMags()
        for i, filtername in enumerate(filters):
            np.testing.assert_allclose(mags1[i, :], mag2[filtername], rtol=1e-4)

    def testGetComputed(self):
        """
        Make sure we can recover computed values.
        """

        sm = self.sm_mags
        sm.setRaDecMjd([36., 36.], [-68., -70.], 49353.18, degrees=True)
        valDict = sm.getComputedVals()

        attrToCheck = ['ra', 'dec', 'alts', 'azs', 'airmass', 'solarFlux', 'moonPhase',
                       'moonAz', 'moonAlt', 'sunAlt', 'sunAz', 'azRelSun', 'moonSunSep',
                       'azRelMoon', 'eclipLon', 'eclipLat', 'moonRA', 'moonDec', 'sunRA',
                       'sunDec', 'sunEclipLon']

        for attr in attrToCheck:
            assert(attr in valDict.keys())
            if np.size(valDict[attr]) > 1:
                np.testing.assert_array_equal(getattr(sm, attr), valDict[attr])
            else:
                assert(getattr(sm, attr) == valDict[attr])

        # Check that things that should be radians are in radian range
        radList = ['ra', 'azs', 'moonAz', 'sunAz', 'azRelSun',
                   'azRelMoon', 'eclipLon', 'moonRA', 'sunRA', 'sunEclipLon']

        for attr in radList:
            if valDict[attr] is not None:
                assert(np.min(valDict[attr]) >= 0)
                assert(np.max(valDict[attr])) <= 2.*np.pi

        # Radians in negative to positive pi range
        radList = ['moonAlt', 'sunAlt', 'alts', 'dec', 'moonDec',
                   'sunDec', 'eclipLat']
        for attr in radList:
            if valDict[attr] is not None:
                assert(np.min(valDict[attr]) >= -np.pi)
                assert(np.max(valDict[attr])) <= np.pi

    def test90Deg(self):
        """
        Make sure we can look all the way to 90 degree altitude.
        """
        mjd = 56973.268218
        sm = self.sm_mags
        sm.setRaDecMjd(0., 90., mjd, degrees=True, azAlt=True)
        mags = sm.returnMags()
        for key in mags.dtype.names:
            assert(True not in np.isnan(mags[key]))
        assert(True not in np.isnan(sm.spec))

    def testAirglow(self):
        """
        test that the airglow goes up with increasing SFU
        """

        mjd = 56973.268218
        sm = self.sm_mags
        sm.setRaDecMjd(0., 90., mjd, degrees=True, azAlt=True, solarFlux=130.)
        magNormal = sm.returnMags()
        sm.setRaDecMjd(0., 90., mjd, degrees=True, azAlt=True, solarFlux=50.)
        magFaint = sm.returnMags()
        sm.setRaDecMjd(0., 90., mjd, degrees=True, azAlt=True, solarFlux=200.)
        magBright = sm.returnMags()

        assert(magNormal[0][2] < magFaint[0][2])
        assert(magNormal[0][2] > magBright[0][2])

    def test_setRaDecAltAzMjd(self):
        """
        Make sure sending in self-computed alt, az works
        """
        sm1 = self.sm_mags
        sm2 = self.sm_mags2
        ra = np.array([0., 0., 0.])
        dec = np.array([-.1, -.2, -.3])
        mjd = 5900
        sm1.setRaDecMjd(ra, dec, mjd)
        m1 = sm1.returnMags()
        sm2.setRaDecAltAzMjd(ra, dec, sm1.alts, sm1.azs, mjd)
        m2 = sm1.returnMags()

        attrList = ['ra', 'dec', 'alts', 'azs']
        for attr in attrList:
            np.testing.assert_equal(getattr(sm1, attr), getattr(sm2, attr))

        for key in m1.dtype.names:
            np.testing.assert_allclose(m1[key], m2[key], rtol=1e-6)


if __name__ == "__main__":
    unittest.main()
