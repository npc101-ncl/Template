#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 00:12:03 2019

@author: peter
"""

import os, re
from pycotools3 import model, tasks, viz
import pandas as pd

rocket = True

if not re.search("/Applications/copasi", os.environ["PATH"]):
    os.environ["PATH"] += os.pathsep + "/Applications/copasi"

working_dir = os.path.abspath('')

copasi_filename = os.path.join(working_dir, "testNoRun.cps")

test_string="""
model test_model()
    var A
    
    R1: ->A; 1
    R2: A->; 1
end
"""
    
theModel = model.loada(test_string, copasi_filename)

theTimeCourse = tasks.TimeCourse(theModel, run=False)

shellString = """#!/bin/bash
# Example SLURM job script for serial (non-parallel) jobs
#
#
# Tell SLURM if you want to be emailed when your job starts, ends, etc.
# Currently mail can only be sent to addresses @ncl.ac.uk
#
#SBATCH --mail-type=ALL
#SBATCH --mail-user=peter.clark@ncl.ac.uk
#


"""

if rocket:
    shellString = shellString+"CopasiSE "+copasi_filename
    shellPath = os.path.join(working_dir, "scratch1part2.sh")
    f = open(shellPath, 'w')
    f.write(shellString)
    f.close()
    os.system("sbatch "+shellPath)
else:
    os.system("CopasiSE "+copasi_filename)

while True:
    try:
        myTable = viz.Parse(theTimeCourse)
        break
    except:
        continue

print(myTable)