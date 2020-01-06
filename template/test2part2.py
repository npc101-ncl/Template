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
f.close()

f = open(os.path.join(data_dir,'test1params0.p'), "rb" )
params0 = pickle.load(f)
f.close()

f = open(os.path.join(data_dir,'test1params1.p'), "rb" )
params1 = pickle.load(f)
f.close()

f = open(os.path.join(data_dir,'test1data0.p'), "rb" )
data0 = pickle.load(f)
f.close()

f = open(os.path.join(data_dir,'test1data1.p'), "rb" )
data1 = pickle.load(f)
f.close()

print(data0)
print(data1)
print(params0)
print(params1)

myTCV0 = timeCourseVisualiser(data0)

myTCV0.multiPlot()

myTCV1 = timeCourseVisualiser(data1)

myTCV1.multiPlot()

PEVis = parameterEstimationVisualiser([params0,params1])

PEVis.violinPlot(indexSelect=[0,1])

PEVis.waterFall()