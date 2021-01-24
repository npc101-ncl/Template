#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 13:59:15 2021

@author: peter
"""

import tellurium as te

# the antimony string used in the model
antStr="""
model test_model()
    var A
    var B
    
    R1: ->A; k1
    R2: A->; k2*A*B
    R3: ->B; k3*A
    R4: B->; k4*B
    
    A = 0;
    B = 0;
    
    k1 = 1;
    k2 = 2;
    k3 = 3;
    k4 = 4;
end
"""

# load model
r = te.loada(antStr)

# simulate from 0 to 50 with 100 steps
r.simulate(0, 10, 100)
# plot the simulation
r.plot()