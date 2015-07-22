import numpy as np
import lsst.sims.skybrightness as sb
import healpy as hp
from lsst.sims.utils import altAzPaFromRaDec
import healpy as hp
from lsst.sims.maf.utils.telescopeInfo import TelescopeInfo
import ephem
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import lsqr
import matplotlib.pylab as plt

#######
# This uses the all sky camera database, available at:
# http://lsst-web.ncsa.illinois.edu/~coughlin/git-repo/catalognight/

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
    ids = np.empty(np.size(intids))

    order = np.argsort(intids)
    ointids = intids[order]
    left = np.searchsorted(ointids,uintids,side='left' )
    right = np.searchsorted(ointids,uintids,side='right' )
    for i in range(np.size(left)): ids[left[i]:right[i]]=uids[i]
    result = np.zeros(np.size(intids), dtype=dtype)
    result[order] = ids

    return result


# Try and run ubercal on the cannon data.  Hopefully get out a

# Load up the telescope properties, has .lat and .lon
telescope = TelescopeInfo('LSST')

nside = 4

filt = 'R'

starIDs=[]
dateIDs = []
hpIDs = []
starMags= []
starMags_err = []
mjds = []
airmasses = []

altLimit = 10. # Degrees
sunAltLimit = np.radians(-20.)

starBrightLimit = -9.

# get the max dateID
maxID,mjd = sb.allSkyDB(0,'select max(ID) from dates;', dtypes='int')
maxID = np.max(maxID)

minMJD = 56900
minID,mjd = sb.allSkyDB(0,'select min(ID) from dates where mjd > %i;' % minMJD, dtypes='int')

names = ['ra','dec','starID','starMag', 'starMag_err', 'sky', 'filter']
types = [float,float,float, float,float,float,'|S1']
dtypes = zip(names,types)

# Temp to speed things up
maxID = 10000
#maxID= 300

for dateID in np.arange(minID.max(),minID.max()+maxID+1):
    sqlQ = 'select stars.ra, stars.dec, stars.ID, obs.starMag_inst, obs.starMag_err,obs.sky, obs.filter from obs, stars, dates where obs.starID = stars.ID and obs.dateID = dates.ID and obs.filter = "%s" and obs.dateID = %i and obs.starMag_err != 0 and dates.sunAlt < %f and obs.starMag_inst > %f;' % (filt,dateID,sunAltLimit, starBrightLimit)

    # Note that RA,Dec are in degrees
    data,mjd = sb.allSkyDB(dateID, sqlQ=sqlQ, dtypes=dtypes)
    if data.size > 0:
        alt,az,pa = altAzPaFromRaDec(np.radians(data['ra']), np.radians(data['dec']),
                                     telescope.lon, telescope.lat, mjd)

        # Let's trim off any overly high airmass values
        good = np.where(alt > np.radians(altLimit))
        data = data[good]
        hpids = hp.ang2pix(nside, np.pi/2.-alt[good], az[good])
        # Extend the lists
        starIDs.extend(data['starID'].tolist())
        dateIDs.extend([dateID]*np.size(data))
        hpIDs.extend(hpids.tolist())
        starMags.extend(data['starMag'].tolist() )
        starMags_err.extend( data['starMag_err'].tolist())
        mjds.extend([mjd]* np.size(data))
        airmasses.extend((1./np.cos(np.pi/2-alt[good])).tolist() )


# switch to arrays
starIDs = np.array(starIDs)
dateIDs = np.array(dateIDs)
hpIDs = np.array(hpIDs)
starMags = np.array(starMags)
starMags_err = np.array(starMags_err)
mjds = np.array(mjds)
airmasses = np.array(airmasses)

# Clobber the dateID and replace it with mjd's rounded to some timestep
timestep = 15. # minutes
dateIDs = np.round(mjds*24.*60./timestep)
mjds = np.round(mjds*24.*60./timestep)*timestep/24./60.

# Need to construct the patch IDs.  Unique id for each mjd+hp combination
multFactor = 10.**np.ceil(np.log10(np.max(hpIDs)))
patchIDs = dateIDs*multFactor+hpIDs
intPatchIDs, upatchIDs, uintPatchids = id2intid(patchIDs)

#
intStarIDs, ustarIDs,uintStarIDs = id2intid(starIDs)

# Construct and solve the sparse matrix
# Using the simple equation:
# m_ij = m_i + ZP_j
# where m_ij is an observed magnitude
# m_i is the true stellar magnitude
# and ZP_j is the patch the id (a combination of time and alt-az)

nObs = len(starMags)
row1 = np.arange(nObs)

col = np.append(intStarIDs, np.max(intStarIDs)+1 + intPatchIDs)
#col = np.append(col, np.zeros(intStarIDs.size)+ np.max(np.max(intStarIDs)+1 + intPatchIDs)+1 )
data = np.ones(intStarIDs.size, dtype=float)
data = np.append(data/starMags_err, data/starMags_err)
#data = np.append(data, airmasses/starMags_err)
row = np.append(row1,row1)
#row = np.append(row,row1)


b = starMags/starMags_err
A = coo_matrix( (data,(row,col)), shape = (nObs,np.max(col)+1))
A = A.tocsr()
solution = lsqr(A,b,show=True)

patchZP = solution[0][np.max(intStarIDs)+1:]
#airmassK = solution[0][-1]
# Need to back out the resulting patchID and dateID, hpid...

resultPatchIDs = intid2id(np.arange(patchZP.size), uintPatchids, upatchIDs)
resultDateIDs = np.floor(resultPatchIDs/multFactor)
resultHpIDs = resultPatchIDs - resultDateIDs*multFactor

# Here's what the best fit came up with:
resultObsMags = A.dot(solution[0])

residuals = resultObsMags*starMags_err - starMags

# so all the say, patches at helpix #44 should be
hpwanted = 8
zps = patchZP[np.where(resultHpIDs == hpwanted) ]
good = np.where(zps != 0)
#plt.plot(zps[good])

frame = 300
good = np.where( resultDateIDs == frame)
skymap = np.zeros(hp.nside2npix(nside))
skymap[resultHpIDs[good].astype(int)] = patchZP[good]

skymap[np.where(skymap == 0)] = hp.UNSEEN
#hp.mollview(skymap, rot=(0,90))


## XXX--Add a snippet of healpy to convert hpid to alt/az
lat,resultAz = hp.pix2ang(nside, resultHpIDs.astype(int))
resultAlt = np.pi/2.-lat

# Convert resultDateIDs to mjd. I think this should work--presumably the dateIDs and mjds are both increasing?
resultMjds = intid2id(resultDateIDs, np.unique(dateIDs), np.unique(mjds), dtype=float)




# Let's figure out the number of stars per patch:
bins = np.zeros(resultPatchIDs.size*2, dtype=float)
bins[::2] = resultPatchIDs-0.5
bins[1::2] = resultPatchIDs+0.5
starsPerPatch,bins = np.histogram(patchIDs, bins=bins)
starsPerPatch = starsPerPatch[::2]
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(starsPerPatch, patchZP, 'ko', alpha=.1)
ax.set_xlabel('Number of stars per patch')
ax.set_ylabel('Patch zeropoint (mags)')
fig.savefig('Uber/zpDist.png')
plt.close(fig)


fig = plt.figure()
ax = fig.add_subplot(111)
good = np.where(resultHpIDs == 0)
floatzp = np.median(patchZP[good])
sc = ax.scatter(resultMjds[good]-resultMjds[good].min(), patchZP[good]-floatzp,
                c=starsPerPatch[good], edgecolor='none')
cb = fig.colorbar(sc, ax=ax)
cb.set_label('Number of stars')
ax.set_ylabel('Patch Zeropoint (mags)')
ax.set_xlabel('MJD-min(MJD) (days)')
#ax.set_xlim([0,4])
#ax.set_ylim([-1,1])
ax.set_title('nside = %i' % nside)
fig.savefig('Uber/zpEvo.png')
plt.close(fig)

nearZ = np.where(resultAlt > np.radians(75.))

fig = plt.figure()
ax = fig.add_subplot(111)
sc = ax.scatter(resultMjds[nearZ]-resultMjds[nearZ].min(), patchZP[nearZ]-floatzp,
           c=starsPerPatch[nearZ], edgecolor='none')
cb = fig.colorbar(sc, ax=ax)
cb.set_label('Stars Per Patch')
ax.set_title('Near Zenith Points')
ax.set_ylabel('Patch Zeropoint (mags)')
ax.set_xlabel('MJD-min(MJD) (days)')
fig.savefig('Uber/zpEvo_Zenith.png')
plt.close(fig)


# Check the residual distribution
fig = plt.figure()
ax = fig.add_subplot(111)
ax.hist(residuals, bins=80, log=True)
ax.set_xlabel('Model-Observed (mags)')
ax.set_ylabel('#')
fig.savefig('Uber/residualHist.png')
plt.close(fig)


# Now to trim off the nights where the zeropoints go very negative
clipLimit = -0.2
tooLow = np.where( (patchZP-floatzp) < clipLimit)[0]
# Reject nieghbor points
tooLow = np.unique(np.array([tooLow-1, tooLow, tooLow+1]).ravel())
allpts = np.arange(patchZP.size)
good = np.setdiff1d(allpts,tooLow)

clippedPoints = np.in1d(patchIDs,resultPatchIDs[good])




# Hmm, I'm still getting wacky negative outliers (and a negative airmass correction)
# Some possible solutions:
# 1) Run a larger stretch of time and hope that connects things
# 2) combine patches that are near in time!  Looks like they are 71s gaps now.
# Could combine every 20-30 min together.
