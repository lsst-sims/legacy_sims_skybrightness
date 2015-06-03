#!/usr/bin/python

import os
import io
import sys
import string
import urllib
import urllib2
import time
import itertools

# Script to execute the ESO SkyCalc web application http://www.eso.org/observing/etc/skycalc/skycalc.htm in batch mode.
# Prepared on October 7, 2014 by Jakob Vinther (jvinther@eso.org) at European Southern Observatory ESO.
# It is prepared specifically for Peter Yoachim at LSST. Please do not re-ditribute.
#
#
# Tested with python 2.7. Maybe you need to install some packages (see the imports in the beginning of the script).
#
# Here is how to use it:
#(1) Open the SkyCalc web application in a browser with this special URL (this is very important to initialize your script):
#http://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC+SKYCALC.POSTFILE.FLAG.SET=1
#
# (2) Configure the parameters in the web application as you please and submit it.
#
# (3) In the beginning of the output page you will see a link to a file with the POST data. Download this file (default name is post.txt but you can change it after you have donloaded it).
#
# (4) On your local linux machine go to a directory where you want to put the FITS files produced with the sky model.
#
# (5) Make the python script executable with the command:
# chmod +x skycalc.py
#
# (6) Open skycalc.py in an editor and modify the part which is indicated to be edited by users. There you can specify the parameters you want to execute the skycalc with. The script will loop over all combinations of values for the keywords. Note that this can quickly be a huge number, so please try to limit it to a reasonable number  not to stress the web server too much.
#
# The keyword names correspond to the input parameters in the form on the SkyCalc input page. (look at the web page source in the browser, to see what the parameters names are).
#
# (7) Execute the python script with the downloaded post.txt file (or whatever you called it):
# skycalc.py post.txt
#
# (8) You should see some printouts in the shell as the script is progressing through all the parameter value combinations.
#
# (9) Note that normally when the sky model is executed through the browser, some javascript will validate the consistency of parameters before it can be submitted. This validation mechanism is not provided with the script, so you must verify the params yourself. If the parameters lead to an error so that the skymodel cannot succeed, the script will print 'error' next to the name of failed FITS file.
#
# (10) If you plan to execute the script with many and/or big files, could you please let me know (jvinther@eso.org) just before you start it so that I can keep an eye on the web server.
#
# (11) Please do not distribute the python script at this stage.



###################################################################################################################
#
# The file post.txt can be obtained from the skycalc web application with this special URL:
# http://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC+SKYCALC.POSTFILE.FLAG.SET=1
# Fill the parameters and submit - a link to the corresponding post.txt will be provided at the top of the output page.
# Download the file post.txt, modify the combinations below and call this script with:
#
# skycalc.py post.txt
#
# In the following configuartion, the keywords in capitals correspond to input parameters in the SkyCalc input web page.
# The SkyCalc web page has pre-submission validation (with JavaScript). The script does not have much validation, so the user must do it.
# Note that some parameters are limited to certain discrete values and interdependent ranges, as indicated in the SkyCalc input web page.

# The following parameter set is an example executing the SkyCalc web application 3*3*2=18 times.
# It can be keept as comments in this script for reference.
#params = [
#('SKYMODEL.TARGET.AIRMASS',['1.00','1.50','2.0']),
#('SKYMODEL.TIME',          ['1','2','3']),
#('SKYMODEL.SEASON',        ['0','1']),
#]

###################################################################################################################
# Normally the user should only configure the following keywords and lists of parameter values:
###################################################################################################################

params = [
    ('SKYMODEL.TARGET.AIRMASS',['1.0','1.1', '1.2', '1.3', '1.4', '1.50','1.6', '1.7', '1.8', '1.9', '2.0', '2.5']),
    ('SKYMODEL.TIME', ['1','2','3'])
    ]

####################################################################################################################
# Normally do not change anything beyond this
####################################################################################################################


# server and scripts config
server='http://www.eso.org'
locurl='/observing/etc/bin/simu/skycalc'
url=server+locurl
deleter_script_url=server+'/observing/etc/bin/script/rmtmp.py'

# check that only one command-line argument is given
if(len(sys.argv) != 2):
    print 'usage: skycalc.py <post-data-file>'
    print ''
    print 'The file post.txt can be obtained from the skycalc web application with this special URL:'
    print 'http://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC+SKYCALC.POSTFILE.FLAG.SET=1'
    print 'Fill the parameters and submit - a link to the corresponding post.txt is provided at the top of the output page.'
    exit(1)

fname=sys.argv[1] # input filename is first command-line argument
if not os.path.exists(fname):
    print fname+' not found'
    sys.exit(1)

with open(fname) as f:
    lines=f.readlines()

if(len(lines) > 1):
    print 'error: the post data file should have exactly one line.'
    print 'but it has '+len(raw_postdata)
    exit(1)

raw_postdata=lines[0].strip('\n')

# create a dictionary of entries encoded in the string
# e.g.: section ('A=5+B=abc+C=8','+','=') returns  {'A':'5','B':'abc','C':'8'}
def section(line,delim,sub_delim): # create a dictionary of entries encoded in the string
    s=line.split(delim)
    dic={}
    for chunk in s:
        p=chunk.partition(sub_delim)
        dic[p[0]]=p[2]
    return dic

def callEtc(url,d):
    data=urllib.urlencode(d)
    req = urllib2.Request(url, data)
    response = urllib2.urlopen(req)
    the_page = response.readlines()
    return the_page

def cleanUp(tmp_dir,deleter_script_url):
    if(deleter_script_url == ''):
        print 'error, deleter_script_url should be non-empty string'
        print 'deleter_script_url='+deleter_script_url
        return 'empty deleter_script_url'

    if(tmp_dir !=''):
        try:
            deleter_response=urllib2.urlopen(deleter_script_url+'?d='+tmp_dir).read().strip('\n') # remove the temp dir and its contents on the server
            if(deleter_response != 'ok'): # it failed somehow
                return deleter_response
        except URLError,e:
            print 'Could not reach script to delete tmp dir on server: '+tmp_dir
            print str(e)
            return 'Could not reach script to delete tmp dir on server: '+tmp_dir

    return 'ok'

def median(x):
    if len(x)%2 != 0:
        return sorted(x)[len(x)/2]
    else:
        midavg = (sorted(x)[len(x)/2] + sorted(x)[len(x)/2-1])/2.0
        return midavg

# decode encoded characters and clean the POST string
postdata=urllib2.unquote(raw_postdata)
postdata=postdata.replace('++','+')
postdata=postdata.replace('=+','=')

# create a dictionary of entries encoded in the string
d=section(postdata,'&','=')

# string to look for in the returned HTML, this line contains the link to the FITS file
match_str='Download the resulting model spectra as FITS table'

# print keywords and values
print ''
print 'Parameter keywords and list of values (the FITS filenames will be named accordingly):'
print ''
for x in params:
    k, v = x
    print k,
    for vi in v:
        print vi,
    print ''
print ''

for i, x in enumerate(params):
    k, v = x
    params[i] = [(k, val) for val in v]

prod=itertools.product(*params)

# count number of executions
n=1
for p in params: n*=len(p)

# initialize progress counter
cntr=1

# timing
times=[]

t0=time.time()

# loop over all combinations of keywords and values,
# call the skycalc web application,
# from the returned HTML extract the URL to the produced FITS file
# download it and delete it from the server temp dir.
for p in prod:
    d.update(dict(p))

    # construct fits file name
    fits_file_name='skytable'
    for t in p:
        fits_file_name += '_'+t[1]
    fits_file_name += '.fits'
    it=str(cntr) + '/' +str(n) + '  '
    info=it+fits_file_name

    # timing
    callEtc_start= time.time()

    #call the sky model on the server
    lines=callEtc(url,d)

    # timing
    callEtc_end = time.time()
    callEtc_secs=str(round(callEtc_end-callEtc_start, 2))


    deleter_response=''
    fits_status=''
    tmp_dir=''

    # timing
    download_fits_start = time.time()

    for line in lines:
        if(match_str in line):
            fitsUrl=line.partition('<a href=')[2].partition('>')[0] # e.g. '/observing/etc/tmp/AAA4sztU/skytable.fits'

            try:
                urllib.urlretrieve(server+fitsUrl,fits_file_name) # retrieve the fits file and write it to a suitable filename indicating the parameters
                #info+= ' download: '+download_fits_secs

                fits_status='saved'
            except URLError,e:
                print 'FITS file could not be retrieved: '+fits_file_name
                print str(e)
                break

            tmp_dir=line.partition('<a href=')[2].partition('>')[0].partition('/observing/etc/tmp/')[2].partition('/skytable.fits')[0] # e.g. 'AAA4sztU'

            if(fits_status=='saved'):
                break # the line was found, the fits file was received and saved

    # out of inner for loop
    download_fits_end = time.time()
    download_fits_secs=str(round(download_fits_end-download_fits_start, 2))

    cleaned_status=cleanUp(tmp_dir,deleter_script_url)
    if(cleaned_status != 'ok'):
        print info+'An internal error occurred on the web server, please report to jvinther@eso.org. status code:"'+cleaned_status+'"'

    if(fits_status != 'saved'):
        fits_status='error'
        sizestr='------'
    else:
        sizestr= str( round(os.path.getsize(fits_file_name)/1024.0/1024.0,4) )

    #timing
    ts=download_fits_end-callEtc_start
    times.append(ts)
    total_secs=str(round(ts, 2))
    eta=str(round((n-cntr)*median(times),2))

    #print header
    if(cntr ==1):
        for i in ['i/n','FITS file\t\t','status','size/Mb\t','exec/s','download/s','total/s','remaining/s']:
            print str(i)+'\t',
        print ''
        print ''
    # print rows
    for i in [it,fits_file_name+'\t',fits_status,sizestr+'\t',callEtc_secs,download_fits_secs+'\t',total_secs,eta]:
        print str(i)+'\t',
    print ''

    cntr=cntr+1

t1=time.time()
script_time=str(round(t1-t0, 2))
print ''
print 'time for all excecutions: '+script_time+' s'

