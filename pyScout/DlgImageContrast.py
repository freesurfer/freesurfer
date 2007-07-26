
import math
import copy
import vtk

from Tkinter import *
import tkSimpleDialog

import Globals

import ImagePipeline

def UpdateLUT():
    Globals.imagePipeline.lut = Globals.imagePipeline.lutBuilder.Build()
    Globals.imagePipeline.UpdateLUT()
    Globals.renWin.Render()

class DlgImageContrast(tkSimpleDialog.Dialog):

    def body(self, master):
        """define the controls of the widget """

        builder = Globals.imagePipeline.lutBuilder
        # backup original LUT
        self.lutBackup = copy.deepcopy(builder)
        
        scalarRange = Globals.imagePipeline.volume.GetScalarRange()
        
        minVal = DoubleVar()
        minVal.set( builder.rangeMin )
        Scale(master,
              from_=scalarRange[0],to=scalarRange[1],
              variable = minVal, orient="horizontal",
              command=self.SetMinValue,
              label="Minimum").pack(side="top",
                                    fill="x",
                                    expand="false")

        maxVal = DoubleVar()
        maxVal.set( builder.rangeMax )
        Scale(master,
              from_=scalarRange[0],to=scalarRange[1],
              variable = maxVal, orient="horizontal",
              command=self.SetMaxValue,
              label="Maximum").pack(side="top",
                                    fill="x",
                                    expand="false")

        # use logarithmic scale for the sigma
        sigmaVal = DoubleVar()
        sigmaVal.set(
            math.log(builder.sigma) )
        Scale(master,
              from_=-5, to=5,
              resolution=.1,
              variable=sigmaVal, orient="horizontal",
              command=self.SetSigma,
              label="Sigma").pack(side="top",
                                  fill="x",
                                  expand="false")

    def SetMinValue(self,val):
        Globals.imagePipeline.lutBuilder.rangeMin = float(val)
        UpdateLUT()
    def SetMaxValue(self,val):
        Globals.imagePipeline.lutBuilder.rangeMax = float(val)
        UpdateLUT()
    def SetSigma(self,val):
        Globals.imagePipeline.lutBuilder.sigma = \
                                               math.exp( float(val))
        UpdateLUT()
        
    def ok(self,event=None):
        tkSimpleDialog.Dialog.cancel(self,event)

    def cancel(self,event=None):
        Globals.imagePipeline.lutBuilder = self.lutBackup
        UpdateLUT()
        tkSimpleDialog.Dialog.cancel(self,event)
    
