# CINEMA_CHANGE_DETECTION
## Version 1.0
## Author: Divya Banesh

A R command line tool that takes as input a Cinema database and produces a new Cinema database with change detection results using Myers et al. This tool also produces a change detection png image and json file to be used with the Cinema:Newsfeed.

**Requirements: R v3.1 or higher (https://cran.r-project.org/)**

**R Studio not required**

**Ensure that any proxy settings to download packages are enabled and firewall is properly set**

# Files
### changeDetection.R
The wrapper file that take in a Cinema data.csv file and supporting arguments, applies the change detection algorithm and write the results to a new column of the data.csv file, as a png graph and as a json. 

### run_lm.R
The Myers et al. change detection algorithm using piecewise linear regression to find change points.

### part_sweep.R
Supporting file for change detection algorithm.

### rss_sweep.R
Supporting file for change detection algorithm. 

# Usage

The command line tool is executed as follows: 

> Rscript changeDetection.R -f path/to/data.csv -p path/to/results/folder -r columnName columnValue columnNames columnValue... -x columnNameOfXParameter -y columnNameOfYParameter -cd B alpha delta

- **changeDetection.R** The wrapper file for change detection.
- **-f path/to/data.csv** This file contains the data to apply change detection on.
- **-p path/to/results/folder** This path is where the CCD folder will be located.
- **-r columnName columnValue** If the user wants to limit the data used for change detection, they can supply the algorithm with a columnName and corresponding columnValue to restrict the rows anlayzed to the values given. This can be any number of columns, one value per column.
- **-x columnNameOfXParameter** Identify a column name to iterate over, after restrictions have been taken into account
- **-y columnNameOfYParameter** Identify a column name to apply change detection to, after restrictions have been taken into account.
- **-cd B alpha delta** Supply the B, alpha and delta parameters for the Myers et al. change detection algorithm (see paper for more details).
- **-o outputName** Identify the name to use for outputfiles for this particular run of change detection.

Example:

> Rscript changeDetection.R -f /Users/genericUser/Desktop/metrics/changeDetection/commandLine/project/nyx.cdb/data.csv -p /Users/genericUser/Desktop/metrics/changeDetection/commandLine/project/nyx.cdb -r phi -180 theta -45  iso 2 -x time -y canny -cd 7 0.7 1 -o iso1_Time_vs_Canny
