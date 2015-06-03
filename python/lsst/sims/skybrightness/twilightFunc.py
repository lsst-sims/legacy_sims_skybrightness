import numpy as np

def twilightFunc(xdata, *args):
    """
    xdata: numpy array with columns 'alt', 'az', 'sunAlt' all in radians.
    az should be relative to the sun (i.e., sun is at az zero.

    based on what I've seen, here's my guess for how to fit the twilight:
    args[0] = decay slope for all pixels
    args[1] = zenith flux at sunAlt = 0
    args[2] = airmass term for hemisphere away from the sun. (factor to multiply max brightness at zenith by)
    args[3] = airmass term for hemisphere towards sun XXX--depreciated, need to remove.
    args[4] = az term for hemisphere towards sun
    args[5:] = constant offsets per pixel (optionall)

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

    #XXX -- could try to make things blend better by just computing everything
    # with the flux[away], then multiplying the flux[towards] region by the cos term.
    # That would eliminate args[3].

    # Do I need that "1+" in there?  That seems odd. Wait, should the flux have the airmass in an exponent?
    # That probably makes more sense! Maybe the same with the cos term.

    flux = args[1]*np.exp(args[2]*(1.-airmass))*np.exp(sunAlt*args[0])
    flux[towards] *= np.exp(args[4]*np.cos(az[towards])*(1.-airmass[towards]))

    #flux[away] = args[1]*(1.+args[2]*airmass[away])*np.exp(sunAlt[away]*args[0])

    #flux[towards] = args[1]*(1.+args[3]*airmass[towards])*np.exp(sunAlt[towards]*args[0]) \
    #                * (args[4]*np.cos(az[towards])+1)

    if np.size(args) >=6:
        flux[away] += args[5:][xdata['hpid'][away]]
        flux[towards] += args[5:][xdata['hpid'][towards]]
                     #+ args[4]*np.arccos(np.cos(az[towards])) \


    return flux
