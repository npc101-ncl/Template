#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 14:53:16 2019

@author: peter
"""

from python.visualisationTools import *
import os, re
import pandas as pd
import pickle

working_dir = os.path.abspath('')
data_dir = os.path.join(working_dir, 'data')

f = open(os.path.join(data_dir,'test1data.p'), "rb" )
data = pickle.load(f)
file.close()

f = open(os.path.join(data_dir,'test1params0.p'), "rb" )
params0 = pickle.load(f)
file.close()

f = open(os.path.join(data_dir,'test1params1.p'), "rb" )
params1 = pickle.load(f)
file.close()

print(data)
print(params0)
print(params1)

myTCV = timeCourseVisualiser(data)

myTCV.multiPlot()

PEVis = parameterEstimationVisualiser([params0,params1])

PEVis.violinPlot(indexSelect=[0,1])

PEVis.waterFall()