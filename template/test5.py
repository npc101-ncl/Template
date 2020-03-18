#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 22:37:46 2020

@author: peter
"""

from python.pycotoolsHelpers import *
import os, re
import pandas as pd
import pickle

working_dir = os.path.abspath('')
data_dir = os.path.join(working_dir, 'data')
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)
dataPath = os.path.join(data_dir, 'test5.csv')

addCopasiPath("/Applications/copasi")

test_string="""
model test_model()
    var A
    
    R1: ->A; k
    R2: A->; B*A
    
    k=1;
    B=1;
end
"""

test_data = pd.DataFrame({'Time':[0,1],'A':[1,1]}, columns= ['Time','A'])
test_data.to_csv(dataPath)

myModel = modelRunner(test_string,os.path.join(working_dir, 'copasiRuns',
                                               'testRun5'))

paramOverride = pd.DataFrame({'k':[1,2]}, columns= ['k'])

params1 = myModel.runParamiterEstimation(dataPath,estimatedVar=["A","B"],
                                         copyNum=2,
                                         overrideParam=paramOverride)

print(params1)