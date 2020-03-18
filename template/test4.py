#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 13:53:33 2020

@author: peter
"""

from python.pycotoolsHelpers import *
import pandas as pd
import os

working_dir = os.path.abspath('')
data_dir = os.path.join(working_dir, 'data')
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)
dataPath = os.path.join(data_dir, 'test4.csv')

addCopasiPath("/Applications/copasi")
rocket=False

test_string="""
model test_model()
    var A
    
    R1: ->A; myPar
    R2: A->; 1
    
    myPar = 1
end
"""

test_data = pd.DataFrame({'Time':[0],'A':[1]}, columns= ['Time','A'])
test_data.to_csv(dataPath)

myModel = modelRunner(test_string,os.path.join(working_dir, 'copasiRuns',
                                               'testRun4'))

myParam = pd.DataFrame([{"myPar":2},{"myPar":3}])

params0 = myModel.runParamiterEstimation(dataPath,copyNum=2,rocket=rocket,
                                         overrideParam=None)