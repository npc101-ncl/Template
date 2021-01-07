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

cmdDict, cmdFlags = getCmdLineArgs()

isOnSlurm = "slurm" in cmdFlags

if "user" in CLDict.keys():
    userID = CLDict["user"]
else:
    userID = "npc101"

if isOnSlurm: 
    myDepth=10
    myThrottle={"user":userID, "jobLimit":1000}
else:
    myDepth=1
    myThrottle=None

if "adjments" in CLDict.keys():
    adjCount = CLDict["adjments"]
    adjCount=int(adjCount)
else:
    adjCount = 19
    
if "perLog10" in CLDict.keys():
    adjPerLog10 = CLDict["perLog10"]
    adjPerLog10=float(adjPerLog10)
else:
    adjPerLog10 = 3

if not isOnSlurm:
    addCopasiPath("/Applications/copasi")

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

runDir = resolvePath(["copasiRuns","example1"], relative=True)

myModel = modelRunner(antString=antStr, run_dir=runDir)

timeCourse = myModel.runTimeCourse(10, stepSize=0.1, rocket=isOnSlurm)

timeCourse = timeCourse[["Time","A","B"]]

expFile1 = resolvePath(["copasiRuns","example1","calData1.csv"],
                       relative=True)

timeCourse.to_csv(expFile1)

timeCourse = myModel.runTimeCourse(10, stepSize=0.1, rocket=isOnSlurm,
                                   adjustParams=pd.DataFrame([{"A":1}]))

timeCourse = timeCourse[["Time","A","B"]]

expFile2 = resolvePath(["copasiRuns","example1","calData2.csv"],
                       relative=True)

timeCourse.to_csv(expFile2)

indepCond = pd.DataFrame([{}, {"A":1}])
indepCond = myModel.preProcessParamEnsam(indepCond)

params = myModel.runParamiterEstimation([expFile1, expFile2],
                                        estimatedVar = ["k1", "k2", "k3",
                                                        "k4", "A", "B"],
                                        copyNum = 100,
                                        rocket = isOnSlurm,
                                        indepToAdd = indepCond,
                                        method="particle_swarm_default")

params = GFID(params)
params = params[[col for col in params.columns if "RSS"!=col]]

params = myModel.runParamiterEstimation([expFile1, expFile2],
                                        estimatedVar = ["k1", "k2", "k3",
                                                        "k4", "A", "B"],
                                        copyNum = 1,
                                        rocket = isOnSlurm,
                                        indepToAdd = indepCond,
                                        overrideParam = params,
                                        randStartVal = False,
                                        method={"method":"hooke_jeeves"})

savePick(["data", "example", "param.p"], params, relative=True)

params = GFID(params)
params = params[[col for col in params.columns if "RSS"!=col]]

timeCourse = myModel.runTimeCourse(10, stepSize=0.1, rocket=isOnSlurm,
                                   adjustParams=params, genReactions=True)

savePick(["data", "example", "timeCourses.p"], timeCourse, relative=True)

adjRange = [10**((x-(adjCount-1)//2)/adjPerLog10)
            for x in range(adjCount)]

myPL = myModel.runProfileLikelyhood([expFile1, expFile2], adjRange,
                                    ["k1", "k2", "k3", "k4", "A", "B"],
                                    rocket = isOnSlurm,
                                    overrideParam = params.iloc[0].to_dict(),
                                    indepToAdd = indepCond, depth = myDepth,
                                    method = {"method":"hooke_jeeves"},
                                    throttle = myThrottle)

savePick(["data", "example", "PL.p"], myPL, relative=True)