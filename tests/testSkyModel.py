import numpy as np
import lsst.sims.skybrightness as sb
import unittest
import lsst.sims.photUtils.Bandpass as Bandpass
import os

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

        sky1.computeSpec()
        sky2.computeSpec()

        np.testing.assert_almost_equal(sky1.spec, sky2.spec)





    def testSetups(self):
        """
        Check that things are the same if the model is set up with
        redecmjd or all the parameters
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

        sm1.computeSpec()
        sm2.computeSpec()

        np.testing.assert_array_equal(sm1.spec, sm2.spec)

        # Check that the degrees kwarg works
        sm2.setParams( azs=np.degrees(sm1.azs), alts=np.degrees(sm1.alts),
                       moonPhase=sm1.moonPhase,
                       moonAlt=np.degrees(sm1.moonAlt), moonAz=np.degrees(sm1.moonAz),
                       sunAlt=np.degrees(sm1.sunAlt), sunAz=np.degrees(sm1.sunAz),
                       sunEclipLon=np.degrees(sm1.sunEclipLon), eclipLon=np.degrees(sm1.eclipLon),
                       eclipLat=np.degrees(sm1.eclipLat), solarFlux=sm1.solarFlux,
                       degrees=True)
        sm2.computeSpec()

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

        throughPath = os.getenv('LSST_THROUGHPUTS_BASELINE')
        filters = ['u','g','r','i','z','y']

        bps = []
        for filterName in filters:
            bp = np.loadtxt(os.path.join(throughPath, 'filter_%s.dat' % filterName),
                            dtype=zip(['wave','trans'],[float]*2 ))
            lsst_bp = Bandpass()
            lsst_bp.setBandpass(bp['wave'], bp['trans'])
            bps.append(lsst_bp)


        sm1 = sb.SkyModel()
        sm1.setRaDecMjd([36.],[-68.],49353.18, degrees=True)
        sm1.computeSpec()
        mags1 = []
        for bp in bps:
            mags1.append(sm1.computeMags(bandpass=bp))
        mags1 = np.array(mags1)

        sm2 = sb.SkyModel(mags=True)
        sm2.setRaDecMjd([36.],[-68.],49353.18, degrees=True)
        sm2.computeSpec()
        mag2 = sm2.computeMags()


        np.testing.assert_allclose(mags1,mag2.T, rtol=1e-4)


if __name__=="__main__":
    unittest.main()
