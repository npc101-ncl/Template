#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 13:27:25 2020

@author: peter
"""

import os, glob, re
import pandas, numpy
import matplotlib.pyplot as plt
import seaborn
from pycotools3 import model, tasks, viz

if not re.search("/Applications/copasi", os.environ["PATH"]):
        os.environ["PATH"] += os.pathsep + "/Applications/copasi"

working_directory = os.path.abspath('')

antimony_string = """
model simple_parameter_estimation()
    compartment Cell = 1;

    A in Cell;
    B in Cell;
    C in Cell;

    // reactions
    R1: A => B ; Cell * k1 * A;
    R2: B => A ; Cell * k2 * B * C;
    R3: B => C ; Cell * k3 * B;
    R4: C => B ; Cell * k4 * C;

    // initial concentrations
    A = 100;
    B = 1;
    C = 1;

    // reaction parameters
    k1 = 0.1;
    k2 = 0.1;
    k3 = 0.1;
    k4 = 0.1;
end
"""

copasi_file = os.path.join(working_directory, 'example_model.cps')

## build model
mod = model.loada(antimony_string, copasi_file)

assert isinstance(mod, model.Model)

## simulate some data, returns a pandas.DataFrame
data = mod.simulate(0, 20, 1)

## write data to file
experiment_filename = os.path.join(working_directory, 'experiment_data.txt')
data.to_csv(experiment_filename)

with tasks.ParameterEstimation.Context(
    mod, experiment_filename,
    context='pl', parameters='gm'
) as context:
    context.set('method', 'hooke_jeeves')
    context.set('pl_lower_bound', 1000)
    context.set('pl_upper_bound', 1000)
    context.set('pe_number', 25) # number of steps in each profile likelihood
    context.set('run_mode', True)
    config = context.get_config()
    
myPE = tasks.ParameterEstimation(config)

data = viz.Parse(myPE).data