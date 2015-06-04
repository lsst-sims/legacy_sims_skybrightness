import numpy as np
import glob
import ephem
from lsst.sims.skybrightness.skyModel import mjd2djd


# Set up LSST observatory
lsstObs = ephem.Observer()
lsstObs.lat = -0.527868529 #radians of '-30:14:40.7'
lsstObs.lon = -1.2348102646986 #radians of '-70:44:57.9'
lsstObs.elevation = 2662.75 #meters


# Package up all the photodiode data into a single rec array

files = glob.glob('ut*.dat')

mjd = []
rband = []
yband = []
zband =[]

for filename in files:
    with open(filename) as f:
        lines = f.readlines()
    for line in lines:
        ack = line.strip().split(' ')
        if len(ack) > 1:
            mjd.append(float(ack[0]))
            rband.append(float(ack[1]))
            yband.append(ack[2])
            zband.append(ack[3])

names = ['mjd', 'r','y','z', 'sunAlt', 'moonAlt', 'moonPhase']
types = ['float']*len(names)
photodiode = np.zeros(len(mjd), dtype=zip(names,types))
photodiode['mjd'] = mjd
photodiode['r'] = rband
photodiode['y'] = yband
photodiode['z'] = zband

# Tack on the moon altitude and phase and sun altitude
for i, mjd in enumerate(photodiode['mjd']):
    lsstObs.date = mjd2djd(mjd)
    sun = ephem.Sun(lsstObs)
    moon = ephem.Moon(lsstObs)
    photodiode['sunAlt'][i] = sun.alt
    photodiode['moonAlt'][i] = moon.alt
    photodiode['moonPhase'][i] = moon.phase


np.savez('photodiode.npz', photodiode=photodiode)



