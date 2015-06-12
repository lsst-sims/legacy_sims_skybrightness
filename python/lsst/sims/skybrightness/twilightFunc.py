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
    args[1] = decay slope for all pixels (mags/radian)
    args[2] = airmass term for hemisphere away from the sun. (factor to multiply max brightness at zenith by)
    args[3] = az term for hemisphere towards sun
    args[4] = zenith dark sky flux
    args[5:] = zenith dark sky times constant (optionall)

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


    # Do I need that "1+" in there?  That seems odd. Wait, should the flux have the airmass in an exponent?
    # That probably makes more sense! Maybe the same with the cos term.

    flux = args[0]*args[4]*10.**(args[1]*(sunAlt+np.radians(12.))+args[2]*(1.-airmass))
    flux[towards] *= 10.**(args[3]*np.cos(az[towards])*(1.-airmass[towards]))

    #flux[away] = args[1]*(1.+args[2]*airmass[away])*np.exp(sunAlt[away]*args[0])

    #flux[towards] = args[1]*(1.+args[3]*airmass[towards])*np.exp(sunAlt[towards]*args[0]) \
    #                * (args[4]*np.cos(az[towards])+1)

    # This let's one fit the dark sky background simultaneously.
    # It assumes the dark sky is a function of airmass only. Forced to be args[4] at zenith.
    if np.size(args) >=6:
        flux[away] += args[4]*np.exp( args[5:][xdata['hpid'][away]]*(1.-airmass[away]))
        flux[towards] += args[4]*np.exp(args[5:][xdata['hpid'][towards]]*(1.-airmass[towards]))

    return flux



def zenithTwilight(alpha, *args):
    """
    The flux at zenith as a linear combination of a twilight component and a constant:
    alpha = sun altitude (radians)
    args[0] = ratio of (zenith twilight flux at sunAlt = -12) and dark sky zenith flux
    args[1] = decay slope for all pixels (mags/radian)
    args[2] = airmass term for hemisphere away from the sun. (factor to multiply max brightness at zenith by)
    args[3] = az term for hemisphere towards sun
    args[4] = zenith dark sky flux
    """

    flux = args[0]*args[4]*10.**(args[1]*(alpha+np.radians(12.))) + args[4]
    return flux
