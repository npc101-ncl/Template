#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 13:30:47 2021

@author: peter
"""

from python.utilityTools import *
from python.visualisationTools import *
import os
import pandas as pd

# creats a absolute path relative to this scripts directory
# the scripts figgures are saved here
saveDir = resolvePath(["figures","example"], relative=True)
# ensures the directory exists
if not os.path.isdir(saveDir):
    os.makedirs(saveDir)

# loads results of the parameter estimation
params = loadPick(["data", "example", "param.p"], relative=True)

# initialises object for visualisation with data from paramiter estimation
myPEVis = parameterEstimationVisualiser(params)
# graphs a waterfall plot for the paramiter estimation
myPEVis.waterFall(save=resolvePath(["figures","example","waterfall.png"],
                       relative=True))

# loads result of multipul simulations
timeCourse = loadPick(["data", "example", "timeCourses.p"], relative=True)

# initialise object for visualising multipul simulations with data
myTCVis = timeCourseVisualiser(timeCourse)

# plot simulations number 1, 2 and 3 corisponding to the same paramiter sets
myTCVis.multiPlot(indexSelect=[1,2,3])

# creat path to copy of experamental data for comparision to results
refTC = resolvePath(["refTC.csv"], relative=True)

# read in copy of experamental data
refTC = pd.read_csv(refTC)

# plot 1st 3 simulations only for variable A and B and compare to
# experamental data
myTCVis.multiPlot(indexSelect=[1,2,3], varSelect=["A","B"], wrapNumber=5,
                  compLines=refTC, save = resolvePath(["figures","example",
                                                       "TimeCourse.png"],
                                                      relative=True))

# make a bar chart representing simulation at time 5
myTCVis.barChart(5, indexSelect=[1,2,3], varSelect=["A","B"],
                 compLines=refTC,
                 save = resolvePath(["figures","example",
                                     "TimeCourse.png"], relative=True))

# load profile likelyhood data
myPL = loadPick(["data", "example", "PL.p"], relative=True)

# initialise profile likelyhood visualisor with data
myPLVis = profileLikelyhoodVisualisor(myPL)

# set paramiters to show in visualisation
showVars = myPL[1].keys()

# plot the RSS variation in the profile likelyhood.
myPLVis.plotRSS(showVars, save = resolvePath(["figures","example",
                                              "profileLikelyhood.png"],
                                             relative=True))