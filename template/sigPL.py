#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 16 20:07:45 2021

@author: peter
"""

import site, os, re
import pandas as pd
from python.pycotoolsHelpers import *
from python.utilityTools import *
import pickle
import time, sys

# if run using the sigPL.sh file will extract information
# passed from the shell script, cmdDict is a dict containing
# keyed paramiters and cmdFlags the flags
cmdDict, cmdFlags = getCmdLineArgs()

# set bounds on allowed paramiter values
myUpperBound=10000
myLowerBound=0.00001

# we extract paramiters from the shell script that called the program
# if it's there, otherwise we set paramiters to defaults.
# the mySuperComputer paramiter is just checking the shell script
# to see if the program should asume it's running on the clustor
mySuperComputer = "slurm" in cmdFlags

# the name of this calculation, principally used to indicate which
# directory to save results in.
if "name" in cmdDict.keys():
    name = cmdDict["name"]
else:
    name = "Hill_coop"

# defines the number of points to include in each profile
if "points" in cmdDict.keys():
    points = int(cmdDict["points"])
    if points%2 == 0:
        points = points+1
else:
    points = 3

# defining the number of points per order of magnitude
if "perLog10" in cmdDict.keys():
    perLog10 = float(cmdDict["perLog10"])
else:
    perLog10 = 1

# allows for overide of the methiod used in the profile likelyhood
if "meth" in cmdDict.keys():
    methDict = cmdDict["meth"]
    methDict = [l.split(":") for l in methDict.split(",")]
    methDict = {l[0]:l[1] for l in methDict if len(l)==2}
    for mKey in methDict.keys():
        try:
            methDict[mKey] = float(methDict[mKey])
            if methDict[mKey].is_integer():
                methDict[mKey] = int(methDict[mKey])
        except:
            pass
else:
    methDict = None

# the number of repeats to do of each paramiter estimation
if "depth" in cmdDict.keys():
    myDepth = cmdDict["depth"]
    try:
        myDepth = int(myDepth)
    except:
        myDepth = 3
else:
    myDepth = 3

# open up the paramiter estimation set and the useful paramiter set
# note RS should have had aditional data added to it by sigOptRead.py
newParams = loadPick(['data', name,"eqParams.p"],relative=True)
RS = loadPick(['data', name,'runSwitches.p'],relative=True)

# define the run directory for the profile likelyhood code to keep
# files in and creat it if it doesn't exist
run_dir = resolvePath(['copasiRuns', name+'-sigOpt'],relative=True)
if not os.path.isdir(run_dir):
    os.makedirs(run_dir)

# we are constructing fake experamental data so we generate a data frame
# with just a Time column containing the right number of time points
# for the given duration and 
calData = pd.DataFrame(data={'Time': [i*RS["duration"]/RS["intervals"] 
                                      for i in range(RS["intervals"]+1)]})

# We add a column of 0s for each calibration chanel in our experamental
# data
for chan in RS["calChans"]:
    calData[chan]=0

# create an object to help you work with the model
myModel = modelRunner(RS["antimony_string"], run_dir)

# make a path to save the experamental file used for calibration in
calPath = resolvePath(['copasiRuns', name+'-sigOpt',"myCal.csv"],
                      relative=True)

# make a list of paths all to the same experamental data file
# we can reuse it for each signal because we want to match the
# calibration chanels as close as posable to constant 0
calPaths = [calPath for _ in RS['signalSet']]

# RS["PLSets"] should contain a dictionary of lists of paramiter set
# indexs we want to generate profile likelyhoods for we loop over them
for setType, rows in RS["PLSets"].items():
    for row in rows:
        # we clear the directory to stop file build up
        myModel.clearRunDirectory()
        # we save the all 0 calibration data to file
        calData.to_csv(calPath,index=False)
        # we extract the relivent row of the paramiter estimation
        # as a dictionary
        myOverride = GFID(newParams).iloc[row].to_dict()
        # we filter out the RSS entry from the dictionary
        myOverride = {i:v for i,v in myOverride.items() if i!="RSS"}
        
        # we define the range of difrent variation values we want in the
        # profile likelyhood
        myRange = [10**((i-(points-1)//2)/perLog10)
                   for i in range(points)]
        
        # we define the varables to be estimated
        estVars = RS["oParams"]
        
        # we run a profile likelyhood. This is a colection of paramiter
        # estimations where one estimated value at a time is heald
        # constant at a range of difrent values. The estimated values
        # defined by oParams. Estimations calibrated using the
        # experamental data files pointed to by the paths in
        # calPaths which are associated with the initial conditions in
        # signalSet. The estimations are done relative to the 
        # paramiter set in myOverride which is used as the starting point
        # for tthe estimations and with myDepth copies of eac estimation.
        myPL = myModel.runProfileLikelyhood(calPaths, myRange, estVars,
                                            rocket=mySuperComputer, 
                                            overrideParam=myOverride,
                                            indepToAdd=RS['signalSet'],
                                            upperParamBound=myUpperBound,
                                            lowerParamBound=myLowerBound,
                                            depth=myDepth,
                                            method = methDict)
        
        # the results of the profile liklehood is saved
        savePick(["data",name,"proLick-"+str(row)+".p"], myPL,
                 relative=True)
        
        # we reload the RS varaiable incase it changed durin run time
        RS = loadPick(['data', name,'runSwitches.p'],relative=True)
        
        # we add the data about this particular profile likelyhoods set up
        RS["PL-"+str(row)]={"range":myRange, "points":points,
                            "type":setType}
        # we re save RS
        savePick(['data', name,'runSwitches.p'],RS,relative=True)
