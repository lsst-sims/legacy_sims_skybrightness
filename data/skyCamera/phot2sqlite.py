import numpy as np
import glob, sys
import ephem
from lsst.sims.maf.utils.telescopeInfo import TelescopeInfo
from lsst.sims.skybrightness.utils import mjd2djd


# Set up LSST telescope
telescope = TelescopeInfo('LSST')
Observatory = ephem.Observer()
Observatory.lat = telescope.lat
Observatory.lon = telescope.lon
Observatory.elevation = telescope.elev

sun = ephem.Sun()
moon = ephem.Moon()

# Need to read through all the files to make the starID's
bfiles = glob.glob('B/*.txt')
rfiles = glob.glob('R/*.txt')
gfiles = glob.glob('G/*.txt')

allFiles = bfiles + rfiles + gfiles

ra=np.zeros(0, dtype=float)
dec=np.zeros(0, dtype=float)
radec=np.zeros(0, dtype=float)
names = ['mjd', 'ra','dec']
types = [float]*3
mjds = np.zeros(0, dtype=float)


for filename in allFiles:
    data = np.loadtxt(filename, dtype=zip(names,types))
    tempRaDec = np.round(data['dec']*1000)*10+data['ra']
    tempRaDec, uind = np.unique(tempRaDec, return_index=True)

    radec = np.append(radec, tempRaDec)
    ra = np.append(ra, data['ra'][uind])
    dec = np.append(dec, data['dec'][uind])

    radec, uind = np.unique(radec, return_index=True)
    ra = ra[uind]
    dec = dec[uind]
    mjds = np.append(mjds, data['mjd'])
    mjds = np.unique(mjds)


# Generate starID table
starids = np.arange(ra.size)
f = open('starTable.dat', 'w')
for i,rai in enumerate(ra):
    print >>f, '%i,%f,%f,0' % (i,ra[i],dec[i])
f.close()

# Generate mjd table
f = open('mjdTable.dat', 'w')
mjdID = np.arange(mjds.size)
for mjdid,mjd in zip(mjdID,mjds):
    Observatory.date = mjd2djd(mjd)
    sun.compute(Observatory)
    moon.compute(Observatory)
    print >>f, '%i,%f,%f,%f,%f' % (mjdid, mjd, sun.alt, moon.alt, moon.phase)
f.close()


# now to loop through and write the obs table and the mjd table
names = ['mjd', 'ra','dec', 'm','alt', 'dm', 'sky', 'band']
types = [float]*7
types.append('|S1')

obsidMax = 0

f = open('obsTable.dat', 'w')
maxJ = float(len(allFiles))

for j,filename in enumerate(allFiles):
    # Maybe read in a dummy column, set it to the filter and then stack all of these so they can quickly be sorted?
    data = np.genfromtxt(filename, dtype = zip(names,types),
                         usecols=(0,1,2,11,16,19,20,21))
    if data.size > 0:
        data['band'] = filename[0]

        data.sort(order=['mjd'])
        # Look up starIDs and mjdIDs
        starIDs = np.zeros(data.size, dtype=int)
        mjdIDs = np.zeros(data.size, dtype=int)
        obsIDs = np.arange(data.size) + obsidMax
        obsidMax = obsIDs.max()+1

        left = np.searchsorted(mjds, data['mjd']  )
        right = np.searchsorted(mjds, data['mjd'], side='right')

        for i,ack in enumerate(left):
            mjdIDs[i] = mjdID[left[i]:right[i]]

        tempRaDec =  np.round(data['dec']*1000)*10+data['ra']
        ord = np.argsort(tempRaDec)
        tempRaDec = tempRaDec[ord]
        data = data[ord]
        mjdIDs = mjdIDs[ord]


        left = np.searchsorted(radec, tempRaDec)
        right = np.searchsorted(radec, tempRaDec, side='right')
        for i,ack in enumerate(left):
            starIDs[i] = starids[left[i]:right[i]]

        for i,starid in enumerate(starIDs):
            print >>f, '%i,%i,%i,%f,%f,%f,%f,%s' % (obsIDs[i], starIDs[i], mjdIDs[i], data['alt'][i],
                                                    data['m'][i], data['dm'][i], data['sky'][i],
                                                    data['band'][i])
        progress = j/maxJ*100
        text = "\rprogress = %.1f%%"%progress
        sys.stdout.write(text)
        sys.stdout.flush()
f.close()








