#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 00:05:09 2019

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
dataPath = os.path.join(data_dir, 'test2.csv')

#addCopasiPath("/Applications/copasi")

test_string="""
model test_model()
    var A
    
    R1: ->A; 1
    R2: A->; 1
end
"""

test_data = pd.DataFrame({'Time':[0],'A':[1]}, columns= ['Time','A'])
test_data.to_csv(dataPath)

myModel = modelRunner(test_string,os.path.join(working_dir, 'copasiRuns',
                                               'testRun2'))

params0 = myModel.runParamiterEstimation(dataPath,copyNum=2,rocket=True)

reactions = extractAntReactions(test_string)

params1 = myModel.runParamiterEstimation(dataPath,estimatedVar=["A"],
                                           copyNum=2,rocket=True)

data = myModel.runTimeCourse(1,adjustParams=params1[list(params1)[0]])

print(data)
print(params0)
print(params1)

file = open(os.path.join(data_dir,'test1data.p'), 'wb')
pickle.dump(data, file)
file.close()

file = open(os.path.join(data_dir,'test1params0.p'), 'wb')
pickle.dump(params0, file)
file.close()

file = open(os.path.join(data_dir,'test1params1.p'), 'wb')
pickle.dump(params1, file)
file.close()