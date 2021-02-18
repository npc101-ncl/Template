#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 14:00:58 2021

@author: peter
"""

# some tools I've writen to make handing files etc easyer
from python.utilityTools import *
# import the pycotools extention code (you need pycotools installed to
# use this)
from python.pycotoolsHelpers import *
import pandas as pd

# if run using the sigOptTC1.sh file will extract information
# passed from the shell script, cmdDict is a dict containing
# keyed paramiters and cmdFlags the flags
cmdDict, cmdFlags = getCmdLineArgs()

# we extract paramiters from the shell script that called the program
# if it's there, otherwise we set paramiters to defaults.
# the name of this calculation, principally used to indicate which
# directory data is kept in.
if "name" in cmdDict.keys():
    name = cmdDict["name"]
else:
    name = "Hill_coop"

# the mySuperComputer paramiter is just checking the shell script
# to see if the program should asume it's running on the clustor
isSuperComputer = "slurm" in cmdFlags

# load the paramiter set from the previous calculation as well as
# useful paramiters used in the paramiter estimation.
params = loadPick(["data",name,"eqParams.p"], relative=True)
RS = loadPick(["data",name,"runSwitches.p"], relative=True)

# unwrap the paramiters
params = GFID(params)
# filter out the RSS column
params = params[[col for col in params.columns if col!="RSS"]]

# we need to make a version of the paramiter set for each signal
# make an empty list to hold the variations
plist = []
# we loop over the signals
for sig in RS["signalSet"]:
    # we make a copy of the data frame because modifications would
    # be incorperated into the originol.
    temp = params.copy()
    # we loop over the paramiters defining the signal and a them to
    # our copy of the paramiter estimations
    for col, val in sig.items():
        temp[col]=val
    # we append our variation to the list
    plist.append(temp)
# we join the list of data frames ino one big data frame
params = pd.concat(plist, ignore_index=True)

# we creat an opject to haden the model
myModel = modelRunner(antString=RS["antimony_string"],
                      run_dir=resolvePath(["copasiRuns","sigOptTC"],
                                          relative=True))
# we clear the directory associated with the model to stop it filling up
# with huge numbers of files
myModel.clearRunDirectory()
# the new extended set of paramiters has gaps in it, we use this function
# to fill them in using the default values from the model.
params = myModel.preProcessParamEnsam(params)
# we run a colection of simulations for all the difrent paramiter sets
# in params. We get a list of data frames out.
myTC = myModel.runTimeCourse(RS["duration"], intervals=RS["intervals"],
                             adjustParams=params, rocket=isSuperComputer)

# we save the simulation results to the file.
savePick(["data",name,"timeCourses.p"], myTC, relative=True)