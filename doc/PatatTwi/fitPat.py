import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import curve_fit

## Fit a simple exponential + constant to the twilight data scrapped from Patat et al 2006



def expPc(x,a,const, ratio):

    result = ratio*const*10.**((x+np.radians(12.))*a) + const
    return result


files = {'B':'Bband.dat', 'I':'Iband.dat', 'R':'Rband.dat', 'U':'Uband.dat', 'V':'vband.dat'}

names = ['zenthAng', 'sb']
types = [float]*2
dt = zip(names,types)

fig,ax = plt.subplots()

counter = 0

results = np.zeros(len(files), dtype=zip(['filter','a','r12'],['|S1',float,float]))

order = ['U','B','V','R','I']
for i,key in enumerate(order):
    data = np.genfromtxt(files[key], dtype=dt, comments='#', delimiter=',')

    flux = 10.**(-0.4*data['sb'])
    x = np.radians(90.-data['zenthAng'])
    p0 = [22., np.min(flux), 4.]
    popt, pcov = curve_fit(expPc, x, flux, p0=p0, sigma=flux*.05)

    results['filter'][i] = key
    results['a'][i] = popt[0]
    results['r12'][i] = popt[2]

    ax.semilogy(np.degrees(x),flux, label=key)
    ax.semilogy(np.degrees(x),expPc(x,*popt), linestyle='--',
                color='k', alpha=.5)


ax.set_xlim([-6,-22])
ax.set_xlabel('Sun Altitude (degrees)')
ax.set_ylabel('Flux (arbitrary)')
ax.legend()
ax.axvline(-12, linestyle='--', color='red', alpha=.5)
fig.savefig('patatFits.pdf')


for res in results:
    print  '%s %.2f %.2f' % (res['filter'], res['r12'], res['a'])
