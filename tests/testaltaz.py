import numpy as np
import lsst.sims.skybrightness as sb
import unittest
import lsst.utils.tests
from lsst.sims.utils import Site, _raDecFromAltAz, _altAzPaFromRaDec, ObservationMetaData, haversine


class TestAltAz(unittest.TestCase):

    def testradec2altaz(self):
        np.random.seed(42)
        ra = np.random.rand(100)*np.pi*2
        dec = np.random.rand(100)*np.pi-np.pi/2
        site = Site('LSST')
        mjd = 55000
        omd = ObservationMetaData(mjd=mjd, site=site)

        trueAlt, trueAz, pa = _altAzPaFromRaDec(ra, dec, omd)
        fastAlt, fastAz = sb.stupidFast_RaDec2AltAz(ra, dec,
                                                    site.latitude_rad,
                                                    site.longitude_rad, mjd)
        distanceDiff = haversine(trueAz, trueAlt, fastAz, fastAlt)

        degreeTol = 2.  # 2-degree tolerance on the fast transform
        assert(np.degrees(distanceDiff.max()) < degreeTol)

    def testaltaz2radec(self):
        np.random.seed(42)
        az = np.random.rand(100)*np.pi*2
        alt = np.random.rand(100)*np.pi-np.pi/2
        site = Site('LSST')
        mjd = 55000
        omd = ObservationMetaData(mjd=mjd, site=site)

        trueRA, trueDec = _raDecFromAltAz(alt, az, omd)
        fastRA, fastDec = sb.stupidFast_altAz2RaDec(alt, az, site.latitude_rad,
                                                    site.longitude_rad, mjd)
        distanceDiff = haversine(trueRA, trueDec, fastRA, fastDec)
        degreeTol = 2.  # 2-degree tolerance on the fast transform
        assert(np.degrees(distanceDiff.max()) < degreeTol)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
