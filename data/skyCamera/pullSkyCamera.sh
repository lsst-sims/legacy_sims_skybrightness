#!/bin/sh

wget -cr -nd -np -P B -l 1 -A txt -erobots=off http://lsst-web.ncsa.illinois.edu/~coughlin/git-repo/catalognight/B/ 

# -c continue
# -r recursuve
# -nd don't replicate directory structure
# -np don't go to the parent directory
# -P B put it in directory "B"
# -l 1 only go one level deep
# -a txt only grab txt files
# -erobots=off ignore the robots.txt file on the site

wget -cr -nd -np -P G -l 1 -A txt -erobots=off http://lsst-web.ncsa.illinois.edu/~coughlin/git-repo/catalognight/G/ 

wget -cr -nd -np -P R -l 1 -A txt -erobots=off http://lsst-web.ncsa.illinois.edu/~coughlin/git-repo/catalognight/R/ 

