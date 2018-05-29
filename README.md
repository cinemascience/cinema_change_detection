# CINEMA_CHANGE_DETECTION
## Version 1.0
## Author: Divya Banesh

A 'R' language library that takes as input a Cinema database and produces a new Cinema database with change detection results using Myers et al. This library also produces a change detection png image and json file to be used with the Cinema:Newsfeed.

**Requirements: R v3.1 or higher (https://cran.r-project.org/)**

** R Studio not required**

# Files
### changeDetection.R
The wrapper library that take in a Cinema data.csv file and supporting arguments, applies the change detection algorithm and write the results to a new column of the data.csv file, as a png graph and as a json. 

### run_lm.R
The Myers et al. change detection algorithm using piecewise linear regression to find change points.

### part_sweep.R
Supporting file for change detection algorithm.

### rss_sweep.R
Supporting file for change detection algorithm. 

# Usage

