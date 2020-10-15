#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 14:19:16 2020

@author: peter
"""

from python.pycotoolsHelpers import *
import os
import pandas as pd
from python.utilityTools import *

working_dir = os.path.abspath('')
data_dir = os.path.join(working_dir, 'data')
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)
dataPath = os.path.join(data_dir, 'test8.csv')

addCopasiPath("/Applications/copasi")

test_string="""
model test_model()
    var A
    var B
    var C
    var D
    
    R1: A->B; k1*A
    R2: B->C; k2*B
    R3: C->D; k3*C
    
    A = 1;
    B = 0;
    C = 0;
    D = 0;
    
    k1 = 0;
    k2 = 0;
    k3 = 0;
end
"""

test_data = pd.DataFrame({'Time':[1],'D':[1]}, columns= ['Time','D'])

test_data.to_csv(dataPath)

findVar = ["k1","k2","k3"]

myModel = modelRunner(test_string,os.path.join(working_dir, 'copasiRuns',
                                               'testRun8'))

myModel.clearRunDirectory()
params = myModel.runParamiterEstimation([dataPath],copyNum=3,
                                        estimatedVar=findVar)

print(params)
params = GFID(params)
params = params[findVar]

myModel.clearRunDirectory()
params = myModel.runParamiterEstimation([dataPath],copyNum=1,
                                        estimatedVar=findVar,
                                        method=
                                        {'method':'hooke_jeeves'},
                                        overrideParam=params,
                                        randStartVal=False)

print(params)