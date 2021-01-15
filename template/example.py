#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 17:48:28 2020

@author: peter
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 14:19:16 2020

@author: peter
"""

from python.pycotoolsHelpers import *
from python.utilityTools import *
import os
import pandas as pd

# interprits comandline variables as flags (single words) or keyword
# enteries (patern key:value)
cmdDict, cmdFlags = getCmdLineArgs()

# looks for a comandline flag indicating running on cluster
isOnSlurm = "slurm" in cmdFlags

# looks for keyword entry giving usename of acount running program on
# cluster
if "user" in CLDict.keys():
    userID = CLDict["user"]
else:
    userID = "npc101"

if isOnSlurm: 
    # defines the number of repertitions the profile likelyhood should do to
    # make sure it has the right answer, value greater than one advised for
    # non deterministic methiods
    myDepth=10
    # defines a limit as to how many jobs should be submittered to slurm at a
    # time by the profile likelyhood function.
    myThrottle={"user":userID, "jobLimit":1000}
else:
    myDepth=1
    myThrottle=None

# set number of points to calculate for profile likelyhood
if "adjments" in CLDict.keys():
    adjCount = CLDict["adjments"]
    adjCount=int(adjCount)
else:
    adjCount = 19

# set number of points per factor of 10 in profile likelyhood   
if "perLog10" in CLDict.keys():
    adjPerLog10 = CLDict["perLog10"]
    adjPerLog10=float(adjPerLog10)
else:
    adjPerLog10 = 3

# if not running on a cluster tells python where to find copasise
if not isOnSlurm:
    addCopasiPath("/Applications/copasi")

# the antimony string used in the model
antStr="""
model test_model()
    var A
    var B
    
    R1: ->A; k1
    R2: A->; k2*A*B
    R3: ->B; k3*A
    R4: B->; k4*B
    
    A = 0;
    B = 0;
    
    k1 = 1;
    k2 = 2;
    k3 = 3;
    k4 = 4;
end
"""

# creats a absolute path relative to this scripts directory
# the model's working files are kept here
runDir = resolvePath(["copasiRuns","example1"], relative=True)

# initialises a model based on an antimony string
myModel = modelRunner(antString=antStr, run_dir=runDir)

# runs a simulation optomised for the rocket cluster if on cluster
timeCourse = myModel.runTimeCourse(10, stepSize=0.1, rocket=isOnSlurm)

# selects only the Time A and B outputs of the simulation
timeCourse = timeCourse[["Time","A","B"]]

# creats a absolute path relative to this scripts directory
expFile1 = resolvePath(["copasiRuns","example1","calData1.csv"],
                       relative=True)

# saves the simulation output to csv to be used as experamental data
timeCourse.to_csv(expFile1)

# creat a path for importing the experamental data into the visualisation
# code
refTC = resolvePath(["refTC.csv"], relative=True)

# save time course as csv for use in visualisation code
timeCourse.to_csv(refTC)

# run simulation under the asumption A starts at 1
timeCourse = myModel.runTimeCourse(10, stepSize=0.1, rocket=isOnSlurm,
                                   adjustParams=pd.DataFrame([{"A":1}]))

# selects only the Time A and B outputs of the simulation
timeCourse = timeCourse[["Time","A","B"]]


expFile2 = resolvePath(["copasiRuns","example1","calData2.csv"],
                       relative=True)

# save the simulation as second experamental data set.
timeCourse.to_csv(expFile2)

# define the starting asumptions / paramiter overides for each experamental
# data set
indepCond = pd.DataFrame([{}, {"A":1}])
# fill in any blanks based on values in model
indepCond = myModel.preProcessParamEnsam(indepCond)

# run paramiter extimation on 2 data sets
params = myModel.runParamiterEstimation([expFile1, expFile2],
                                        estimatedVar = ["k1", "k2", "k3",
                                                        "k4", "A", "B"],
                                        copyNum = 100,
                                        rocket = isOnSlurm,
                                        indepToAdd = indepCond,
                                        method="particle_swarm_default")

# unwrap results to give data frame
params = GFID(params)
# purge RSS column from data frame
params = params[[col for col in params.columns if "RSS"!=col]]

# perform 2nd paramiter estimation using output of first paramiter estimation
# as a starting point
params = myModel.runParamiterEstimation([expFile1, expFile2],
                                        estimatedVar = ["k1", "k2", "k3",
                                                        "k4", "A", "B"],
                                        copyNum = 1,
                                        rocket = isOnSlurm,
                                        indepToAdd = indepCond,
                                        overrideParam = params,
                                        randStartVal = False,
                                        method={"method":"hooke_jeeves"})

# save results as a pickel file
savePick(["data", "example", "param.p"], params, relative=True)

# unwrap results to give data frame
params = GFID(params)
# purge RSS column from data frame
params = params[[col for col in params.columns if "RSS"!=col]]

# run colection of simulations using the output of the paramiter estimation
# to overide the models paramiters / initial conditions
timeCourse = myModel.runTimeCourse(10, stepSize=0.1, rocket=isOnSlurm,
                                   adjustParams=params, genReactions=True)

# save results as pickel file
savePick(["data", "example", "timeCourses.p"], timeCourse, relative=True)

# define range of adjustments to use in profile likelyhoods.
adjRange = [10**((x-(adjCount-1)//2)/adjPerLog10)
            for x in range(adjCount)]

# run profile likelyhood on best result from paramiter estimation
myPL = myModel.runProfileLikelyhood([expFile1, expFile2], adjRange,
                                    ["k1", "k2", "k3", "k4", "A", "B"],
                                    rocket = isOnSlurm,
                                    overrideParam = params.iloc[0].to_dict(),
                                    indepToAdd = indepCond, depth = myDepth,
                                    method = {"method":"hooke_jeeves"},
                                    throttle = myThrottle)

# save results to pickel file
savePick(["data", "example", "PL.p"], myPL, relative=True)