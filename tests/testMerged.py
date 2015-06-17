import numpy as np
import lsst.sims.skybrightness as sb
import unittest


class TestMerged(unittest.TestCase):

    def testmergedComp(self):
        """
        Test that the 4 components that have been merged return the
        same result if they are computed independently
        """

        sky1 = sb.SkyModel(twilight=False, zodiacal=False,  moon=False,
                         lowerAtm=False, upperAtm=False,
                         airglow=False, scatteredStar=False,
                         mergedSpec=True)
        sky1.setRaDecMjd([36.],[-68.],49353.18, degrees=True)

        sky2 = sb.SkyModel(twilight=False, zodiacal=False,  moon=False,
                         lowerAtm=True, upperAtm=True,
                         airglow=True, scatteredStar=True,
                         mergedSpec=False)
        sky2.setRaDecMjd([36.],[-68.],49353.18, degrees=True)

        sky1.computeSpec()
        sky2.computeSpec()

        np.testing.assert_almost_equal(sky1.spec, sky2.spec)




if __name__=="__main__":
    unittest.main()
