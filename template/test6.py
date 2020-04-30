#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 10:29:03 2020

@author: peter
"""

from python.pycotoolsHelpers import *
import os, re
import pandas as pd
import pickle
from libsbgnpy import render, utils
from libsbgnpy import libsbgn
from libsbgnpy import Extension, Language

working_dir = os.path.abspath('')

addCopasiPath("/Applications/copasi")

test_string="""
model test_model()
    var A
    
    R1: ->A; myPar
    R2: A->; A
    
    myPar = 1
end
"""
class networkVisualiser:
    def __init__(self,sbgnPath,nodeBindings,edgeBindings,timeCourse=None):
        self.sbgn = utils.read_from_file(sbgnPath)
        self.map = sbgn.get_map()
        self.map.set_language(Language.PD)
        glyphs = self.map.get_glyph()
        validIDs = [g.get_id() for g in glyphs]
        processes = {g.get_id():[p.get_id() for p in g.get_port()] for
                     g in glyphs if g.get_class()=="process"}
        arcs = self.map.get_arc()
        self.edges = {g:[a.get_id() for a in arcs 
                         if (a.get_source() in p) or (a.get_target() in p)]
                      for g, p in processes.items()}
        self.edges = {reaction:self.edges[sbgnID] for reaction, sbgnID in
                      edgeBindings if sbgnID in self.edges.keys()}
        self.nodes = {metabolite:sbgnID for metabolite, sbgnID in
                      nodeBindings if sbgnID in validIDs} 
        self.timeCourse = timeCourse
        
    def sbgnLineStyle(self, ids, intensity, tab="    "):
        if intensity>1 or intensity<-1:
            lineColour = "ffff"
        elif intensity>=0:
            lineColour = hex(int(intensity*(16**4-1)))[2:]
        else:
            lineColour = hex(int(-intensity*(16**4-1)))[2:]
        while len(lineColour)<4:
            lineColour = "0"+lineColour
        lineColour = "#00" + lineColour
        if intensity<0:
            lineColour = lineColour+"80"
            width="1"
        elif intensity>1:
            width="4"
        else:
            width="2"
        gTag = ('<g stroke="' + lineColour + '" stroke-width="' + width +
                '" />')
        if isinstance(ids,list):
            if len(ids)==1:
                sTag = ids
            else:
                sTag = " ".join(ids)+" "
        else:
            sTag = ids
        sTag = '<style idList="' + sTag + '">'
        return sTag+"\n"+tab+gTag+"\n</style>"
    
    def sbgnFillStyle(self, ids, intensity, tab="    "):
        if intensity>1 or intensity<-1:
            fillColour = "0000"
        elif intensity>=0:
            fillColour = hex(int((1-intensity)*(16**4-1)))[2:]
        else:
            fillColour = hex(int((1+intensity)*(16**4-1)))[2:]
        while len(fillColour)<4:
            fillColour = "0"+fillColour
        fillColour = "#" + fillColour + "ff"
        if intensity<0:
            fillColour = fillColour+"80"
            width="1"
        elif intensity>1:
            width="4"
        else:
            width="2"
        gTag = ('<g fill="' + fillColour + '" stroke-width="' + width +
                '" />')
        if isinstance(ids,list):
            if len(ids)==1:
                sTag = ids
            else:
                sTag = " ".join(ids)+" "
        else:
            sTag = ids
        sTag = '<style idList="' + sTag + '">'
        return sTag+"\n"+tab+gTag+"\n</style>"
    
    def genImage(self, frame, fPath):
        if self.timeCourse is None:
            return False
        if frame<0 or frame>=len(self.timeCourse):
            return False
        styleStr = """<renderInformation id="example" programName="SBML Layout" programVersion="3.0"
         xmlns="http://projects.eml.org/bcb/sbml/render/level2">
            <listOfStyles>
        """
        styleStr = styleStr + "\n"
        for index, value in self.timeCourse.iloc[frame].items():
            if index in self.nodes.keys():
                styleStr = styleStr + self.sbgnFillStyle(self.nodes[index],
                                                         value) + "\n"
            if index in self.edges.keys():
                for target in self.edges[index]:
                    styleStr = styleStr + self.sbgnLineStyle(target,
                                                             value) + "\n"
        styleStr = styleStr + """</listOfStyles>
        </renderInformation>"""
        extension = Extension(extString)
        self.map.set_extension(extension)
        render.render_sbgn(self.sbgn, image_file=fPath, file_format="png")

myModel = modelRunner(test_string,os.path.join(working_dir, 'copasiRuns',
                                               'testRun6'))            

myModel.genReactionAntString()

data = myModel.runTimeCourse(1)
data2 = myModel.runTimeCourse(1,genReactions=True)

sbgn = utils.read_from_file(os.path.join(working_dir, "test6.sbgn"))

map = sbgn.get_map()
map.set_language(Language.PD)

extString = """<renderInformation id="example" programName="SBML Layout" programVersion="3.0"
 xmlns="http://projects.eml.org/bcb/sbml/render/level2">
    <listOfStyles>
"""
extString = (extString +"\n"+sbgnLineStyle("R2-react0",0.5) +"\n"+
             sbgnFillStyle("sa0",0.5) +"\n"+ """</listOfStyles>
</renderInformation>""")

extension = Extension(extString)
map.set_extension(extension)

#print(utils.write_to_string(sbgn))

render.render_sbgn(sbgn, image_file=os.path.join(working_dir, "test6.png"),
                   file_format="png")

"https://libsbgn-python.readthedocs.io/en/latest/examples.html#Create-SBGN"