import urllib2
import re
import os.path

# Grab the page that has all the dates:
page = urllib2.urlopen('http://lsst-web.ncsa.illinois.edu/~coughlin/allsky/data').read()
dates = []
for i in range(len(page)):
    if (page[i:i+3] == '"ut') & (page[i:i+4] != '"utc'):
        dates.append(page[i+1:i+9])

# Grab the photodiode data for each date
for date in dates:
    if not os.path.isfile(date+'.dat'):
        try:
            print 'downloading %s' % date
            page = urllib2.urlopen('http://lsst-web.ncsa.illinois.edu/~coughlin/allsky/data/' +
                                   date+'/photodiodeplots/photodiode.txt').read()
            f = open(date+'.dat', 'w')
            print >>f, page
            f.close()
        except:
            pass
