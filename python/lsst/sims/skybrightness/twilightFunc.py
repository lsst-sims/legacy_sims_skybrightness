import numpy as np

__all__=['twilightFunc', 'zenithTwilight', 'simpleTwi']




def simpleTwi(xdata, *args):
    """
    Fit a simple slope and constant to many healpixels

    xdata should have keys:
    sunAlt
    hpid

    args:
    0: slope
    1:hpid: magnitudes
    hpid+1:2*hpid: constant offsets
    """

    args = np.array(args)
    hpmax = np.max(xdata['hpid'])
    result = args[xdata['hpid']+1]*np.exp( xdata['sunAlt'] * args[0]) + args[xdata['hpid']+2+hpmax]
    return result



def twilightFunc(xdata, *args):
    """
    xdata: numpy array with columns 'alt', 'az', 'sunAlt' all in radians.
    az should be relative to the sun (i.e., sun is at az zero.

    based on what I've seen, here's my guess for how to fit the twilight:
    args[0] = ratio of (zenith twilight flux at sunAlt = -12) and dark sky zenith flux
    args[1] = zenith dark sky flux
    args[2] = decay slope for all pixels (mags/radian)
    args[3] = airmass term for hemisphere away from the sun.
    args[4] = az,airmass term for hemisphere towards sun
    args[5] = az,airmass^2 term for hemisphere towards sun
    args[6:] = zenith dark sky times constant (optionall)

    """

    ## XXX--I think I might want to promote this to a free parameter to fit.
    amCut = 1.1

    args = np.array(args)
    az = xdata['azRelSun']
    airmass = xdata['airmass']
    sunAlt = xdata['sunAlt']
    flux = np.zeros(az.size, dtype=float)
    away = np.where( (airmass <= amCut) | ((az >= np.pi/2) & (az <= 3.*np.pi/2)))
    towards = np.where( (airmass > amCut) & ((az < np.pi/2) | (az > 3.*np.pi/2)))


    flux = args[0]*args[1]*10.**(args[2]*(sunAlt+np.radians(12.))+args[3]*(airmass-1.))
    flux[towards] *= 10.**(args[4]*np.cos(az[towards])*(airmass[towards]-1.))

    # Adding in an X^2 term seems to help, but now I might just be fitting the zodiacal light?
    # But maybe the zodical is getting mostly absobed in the constant term?
    flux[towards] *= 10.**(args[5]*np.cos(az[towards])*(airmass[towards]-1.)**2 )
    # Adding cos^2 term didn't do much
    #flux[towards] *= 10.**(args[5]*np.cos(az[towards])**2.*(airmass[towards]-1.) )

    # This let's one fit the dark sky background simultaneously.
    # It assumes the dark sky is a function of airmass only. Forced to be args[4] at zenith.
    if np.size(args) >=7:
        #flux += args[4]*np.exp( args[5]*(airmass-1.))
        # Formulate it this way so that it's like adding a constant magnitude
        flux[away] += args[1]*np.exp( args[6:][xdata['hpid'][away]]*(airmass[away]-1.))
        flux[towards] += args[1]*np.exp(args[6:][xdata['hpid'][towards]]*(airmass[towards]-1.))

    return flux



def zenithTwilight(alpha, *args):
    """
    The flux at zenith as a linear combination of a twilight component and a constant:
    alpha = sun altitude (radians)
    args[0] = ratio of (zenith twilight flux at sunAlt = -12) and dark sky zenith flux
    args[1] = zenith dark sky flux
    args[2] = decay slope for all pixels (mags/radian)
    """

    names = ['azRelSun', 'sunAlt', 'airmass', 'hpid']
    types = [float]*(len(names)-1)
    types.append(int)
    xdata = np.zeros(alpha.size, dtype=zip(names,types))
    xdata['sunAlt'] = alpha

    if len(args) < 4:
        flux = twilightFunc(xdata, args[0], args[1], args[2], 0,0,0)
    else:
        flux = twilightFunc(xdata, args[0], args[1], args[2], 0,0,0, args[3])
    #flux = args[0]*args[1]*10.**(args[2]*(alpha+np.radians(12.))) + args[1]
    return flux
