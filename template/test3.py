#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 10 15:20:35 2020

@author: peter
"""

from python.pycotoolsHelpers import *
import pandas as pd
import os

working_dir = os.path.abspath('')
data_dir = os.path.join(working_dir, 'data')
if not os.path.isdir(data_dir):
    os.makedirs(data_dir)
dataPath = os.path.join(data_dir, 'test1.csv')

addCopasiPath("/Applications/copasi")

test_string="""
model test_model()
    var A
    
    R1: ->A; 1
    R2: A->; 1
end
"""

myModel = modelRunner(test_string,os.path.join(working_dir, 'copasiRuns',
                                               'testRun3'))

