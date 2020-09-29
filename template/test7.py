#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 16:24:36 2020

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
dataPath = os.path.join(data_dir, 'test7.csv')

addCopasiPath("/Applications/copasi")

test_string="""
model test_model()
    var A
    
    R1: ->A; k
    R2: A->; l*A
    
    k=1;
    l=1;
end
"""

test_data = pd.DataFrame({'Time':[0,1],'A':[1,1]}, columns= ['Time','A'])
test_data.to_csv(dataPath)

myModel = modelRunner(test_string,os.path.join(working_dir, 'copasiRuns',
                                               'testRun7'))

paramOverride = pd.DataFrame({'k':[1,2]}, columns= ['k'])

myModel.clearRunDirectory()

profile, refParam = myModel.runProfileLikelyhood([dataPath], [0.9,1.1],
                                                 ["k","l"],
                                                 overrideParam=paramOverride)

print(profile)
print(refParam)