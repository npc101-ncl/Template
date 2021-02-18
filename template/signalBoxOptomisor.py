#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  2 16:27:19 2021

@author: peter
"""

import site, os, re
import pandas as pd
# import the pycotools extention code (you need pycotools installed to
# use this)
from python.pycotoolsHelpers import *
# some tools I've writen to make handing files etc easyer
from python.utilityTools import *
# library for saving python variables to disk
import pickle
import time, sys
import math

# if run using the sigOpt1.sh file will extract information
# passed from the shell script, cmdDict is a dict containing
# keyed paramiters and cmdFlags the flags
cmdDict, cmdFlags = getCmdLineArgs()

# define the paramiters that give the dirent signals. Each
# signal in list corisponds to 1 set of independent conditions
# associated with a difrent 'experament' used in the parameterisation.
signalSet = [{"preSignal":0,   "S_p":0.5*math.pi},
             {"preSignal":1,   "S_p":1.5*math.pi},
             {"preSignal":0.5, "S_f":1/(2*math.pi)},
             {"preSignal":0.5, "S_f":3/(2*math.pi)},
             {"preSignal":0.5, "S_f":9/(2*math.pi)}]

# we set the duration of the experament, the number of intervals in
# the experamental data and the name of the calibration chanels
duration = 20
intervals = 200
calChans = ["Sq2C1","Sq3C1"]

# we are constructing fake experamental data so we generate a data frame
# with just a Time column containing the right number of time points
# for the given duration and 
calData = pd.DataFrame(data={'Time': [i*duration/intervals for i
                                      in range(intervals+1)]})

# We add a column of 0s for each calibration chanel in our experamental
# data
for chan in calChans:
    calData[chan]=0

# we want to try and find an equivalent paramiterisation for each of
# these reaction paramiters. We will look at all posable combinations
# of these posable options
iParamRanges = {"kD1_S":[0.1, 0.5, 1],
                "kD1_V":[0.1, 1, 10],
                "kD1_h":[1, 10, 100],
                "kD1_R":[0.1, 1, 10]}

# we create an empty list to fill with each posable combination of
# the iParamRanges
iParamsDF = []
# we create a dictionary showing the number of posabilities for each
# member of iParamRanges
iParamLens = {k:len(v) for k,v in iParamRanges.items()}
# we initialise k to be 1. k will represent to the total number of
# combinations of iParamRanges
k=1
# we set k by multiplying it by the number of posibilitys for 
# each member of iParamRanges one after the other to get the
# total number of posable combinations.
for i in [v for k,v in iParamLens.items()]:
    k*=i
# we loop through every combination using a number i we will
# link to each posble combination of iParamRanges
for i in range(k):
    # create a new empty dictionary to hold the sate that will
    # be associated with a particular combination
    tDict = {}
    # j will represent the information about the combination we
    # haven't utalised yet, initially just set to i which specifies
    # the combination unequly.
    j = i
    # we loop over the member counts for iParamRanges
    for k, v in iParamLens.items():
        # we use the remainder upon division by the total number
        # of options for a iParamRanges member to select which
        # option should be set for this value of i
        tDict[k] = iParamRanges[k][j%v]
        # we discard the information extracted by the remainder
        # division by setting j to the result of the division with
        # remainder
        j=j//v
    # we add the dictionary represention the combination associated
    # with i to the list
    iParamsDF.append(tDict)
# we convert the list of dictionarys into a data frame represention
# all the relivent combinations of iParamRanges members
iParamsDF = pd.DataFrame(iParamsDF)

# we define the varaibles in the other 2 pathways we want to estimate
oParams = ["kD2_S", "kD2_V", "kD2_h", "kD2_R", "delay2I",
           "kD3_a", "kD3_h", "kD3_R", "delay3I"]

# we set reasionable guesses as to what the uper and lower limit on
# the paramiters should be
myUpperBound=1000
myLowerBound=0.0001

# running on the cluster generally our jobs self terminate after 2 days
# so we give the program a 47 hour deadline to start returning the results
# it has. endTime gives it that cut off time 
secondsToRun = 60*60*47
endTime = time.time()+secondsToRun

# we extract paramiters from the shell script that called the program
# if it's there, otherwise we set paramiters to defaults.
# the mySuperComputer paramiter is just checking the shell script
# to see if the program should asume it's running on the clustor
mySuperComputer = "slurm" in cmdFlags

# allows us to overide the methiod copasi uses for the 1st stage
# of paramiter estimation in the background
if "meth" in cmdDict.keys():
    myMeth = cmdDict["meth"]
else:
    myMeth = "particle_swarm_default"

# the name of this calculation, principally used to indicate which
# directory to save results in.
if "name" in cmdDict.keys():
    name = cmdDict["name"]
else:
    name = "Hill_coop"

# paramiter estimation is an inexact process so we do multipul repeats
# and pick the best (and look to see how close the next few best
# answers were). the default of 3 is rather low but the process is
# intensive
if "copys" in cmdDict.keys():
    myCopyNum = cmdDict["copys"]
else:
    myCopyNum = 3

# we load the antimony string from a text file, the relative arguments
# tells the function to work our the file path from where this script
# is located
antimony_string = loadTxt(["delayAntFile.txt"], relative=True)
# generate a path for the paramiter estimation to take place in
run_dir = resolvePath(['copasiRuns', name+'-sigOpt'],relative=True)
# make a path to save the experamental file used for calibration in
calPath = resolvePath(['copasiRuns', name+'-sigOpt',"myCal.csv"],
                      relative=True)

# create an object to help you work with the model
myModel = modelRunner(antimony_string, run_dir)
# clear the directory to keep from files building up
myModel.clearRunDirectory()
# save the experamental data to file
calData.to_csv(calPath,index=False)

# make a list of paths all to the same experamental data file
# we can reuse it for each signal because we want to match the
# calibration chanels as close as posable to constant 0
calPaths = [calPath for _ in signalSet]

# we run a paramiter estimation search for the paramiters in 
# oParams. Calibrated using the experamental data files pointed
# to by the paths in calPaths which are associated with the
# initial conditions in signalSet. And this estimation is actually
# an ensamble of estimations for each set of initial conditions in
# the data frame iParamsDF and with myCopyNum copies of the estimation
# for each row of iParamsDF.
params = myModel.runParamiterEstimation(calPaths,
                                        copyNum=myCopyNum,
                                        rocket=mySuperComputer,
                                        estimatedVar=oParams,
                                        upperParamBound=myUpperBound,
                                        lowerParamBound=myLowerBound,
                                        method=myMeth,
                                        overrideParam=iParamsDF,
                                        indepToAdd=signalSet,
                                        endTime=endTime) 

# we extract the results of the paramiter estimation and filter out
# the RSS values that estimate the 'goodness' of the fit.
params = GFID(params)
params = params[[col for col in params.columns if col!="RSS"]]

# we repeat the paramiter estimation but this time instead of iParamsDF
# we set the initial condition to be the result of the last estimation.
# this includes the estimated value. Seting randStartVal to False
# ensures the estimation will start with the estimated values the same
# as the output of the last paramiter estimation. This time we set
# copyNum to 1 because we already have the aditional copies in params.
# the methiod use here is a local one that homes in on an answer in a
# much smoother localised proocess. This combination of a paramiter
# estimation methiod that jumps around a lot followed by a local one
# is called global chaser.
params = myModel.runParamiterEstimation(calPaths,
                                        copyNum=1,
                                        rocket=mySuperComputer,
                                        estimatedVar=oParams,
                                        upperParamBound=myUpperBound,
                                        lowerParamBound=myLowerBound,
                                        method={'method':'hooke_jeeves'},
                                        overrideParam=params,
                                        indepToAdd=signalSet,
                                        randStartVal = False,
                                        endTime=endTime)

# we save the output of the paramiter estimation.
savePick(["data",name,"eqParams.p"], params, relative=True)

# we put the useful paramiters into a dict called RS then save it for
# later use
RS={}
RS["signalSet"] = signalSet
RS["duration"] = duration
RS["intervals"] = intervals
RS["calChans"] = calChans
RS["iParamRanges"] = iParamRanges
RS["iParamsDF"] = iParamsDF
RS["oParams"] = oParams
RS["antimony_string"] = antimony_string

savePick(["data",name,"runSwitches.p"], RS, relative=True)