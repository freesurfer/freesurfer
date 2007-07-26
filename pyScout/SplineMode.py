
import time

import Globals
import ImagePipeline

from SplineSet import *
from Laplace import Locator, WrapperPolyData

import Tkinter
import vtk
from vtk.tk.vtkTkRenderWindowInteractor import \
     vtkTkRenderWindowInteractor

def GetEventWorldPoint():
    """ uses the image pipeline to infer the position in global coordinates """
    
    iren = Globals.renWin.GetInteractor()
    x,y = iren.GetEventPosition()

    cam = Globals.ren.GetActiveCamera()
    viewFocus = cam.GetFocalPoint()
    Globals.ren.SetWorldPoint(viewFocus[0],viewFocus[1],viewFocus[2], 1.0)
    Globals.ren.WorldToDisplay()
    z = Globals.ren.GetDisplayPoint()[2]
    
    Globals.ren.SetDisplayPoint(x,y,z)
    Globals.ren.DisplayToWorld()
    world = Globals.ren.GetWorldPoint()
    
    return world[:3]

def ResetStyle():
    istyle = vtk.vtkInteractorStyle()
    iren = Globals.renWin.GetInteractor()
    iren.SetInteractorStyle( istyle )


#############################

class ModeInterface:
    def __init__(self,caller=None):
        pass
    def Start(self):
        pass
    def Finish(self):
        pass

class Navigation(ModeInterface):

    def Start(self):
        istyle = vtk.vtkInteractorStyleImage()
        istyle.AddObserver("RightButtonPressEvent",self.PlacePointer)
        iren = Globals.renWin.GetInteractor()
        iren.SetInteractorStyle(istyle)


    def PlacePointer(self,a,b):
        w = GetEventWorldPoint()
        Globals.crossHair.SetPosition(w[0],w[1])
        Globals.renWin.Render()

class AddSpline(ModeInterface):

    def __init__(self,caller):
        self.caller = caller

    def Finish(self):
        ResetStyle()
        self.caller.imageObservers.remove(self.UpdateObserver)

    def UpdateObserver(self):
        self.spline.Hide()
        self.ResetSpline()

    def Start(self):
        self.caller.imageObservers.append(self.UpdateObserver)
        self.ResetSpline()

        istyle = vtk.vtkInteractorStyleImage()
        istyle.AddObserver("LeftButtonReleaseEvent", self.AddPoint)
        istyle.AddObserver("RightButtonPressEvent", self.AddSpline)
        iren = Globals.renWin.GetInteractor()
        iren.SetInteractorStyle(istyle)

    def AddPoint(self,a,b):
        world = GetEventWorldPoint()
        self.spline.AddPoint(world)
        Globals.renWin.Render()
        
    def AddSpline(self,a,b):
        self.AddPoint(a,b)
        self.spline.rendering.profile.GetProperty().SetColor(1,1,1)
        Globals.objectSet.AddSpline(self.spline)
        self.ResetSpline()
        Globals.renWin.Render()

    def ResetSpline(self):
        self.spline = IndexedSpline(1)
        self.spline.sliceIndex = Globals.imagePipeline.GetSliceIndex()
        self.spline.orientation = Globals.imagePipeline.GetOrientation()
        self.spline.rendering.profile.GetProperty().SetColor(0,1,1)


class CmdDeleteSpline(ModeInterface):

    def __init__(self,caller):
        self.caller = caller
        self.caller.mode = None
        
    def Start(self):
        self.caller.imageObservers.append(self.UpdateObserver)
        self.locator = None
        self.ResetLocator()
        
        istyle = vtk.vtkInteractorStyle()
        istyle.AddObserver("LeftButtonPressEvent", self.Delete)
        iren = Globals.renWin.GetInteractor()
        iren.SetInteractorStyle(istyle)

    def Delete(self,a,b):
        if not len(Globals.objectSet.visibleSplines): return

        wp = GetEventWorldPoint()
        index = self.locator.FindClosestPoint( wp[0],wp[1],wp[2] )

        Globals.objectSet.DeleteVisibleSpline(index)
        Globals.renWin.Render()

        self.Finish()

    def Finish(self):
        ResetStyle()
        self.caller.imageObservers.remove(self.UpdateObserver)
        self.caller.ChangeMode()

    def UpdateObserver(self):
        self.ResetLocator()

    def ResetLocator(self):

        self.locator = Locator()
        
        splines = Globals.objectSet.visibleSplines
        if not len(splines): return
        
        wrapper = WrapperPolyData()
        counter = 0
        for s in splines:
            startCounter = counter
            for pt in s.points:
                wrapper.InsertPoint( counter,
                                     pt[0],pt[1],pt[2]
                                     )
                counter += 1
            wrapper.InsertNextCell( len(s.points) )
            for i in range(startCounter,counter):
                wrapper.InsertCellPoint(i)

        self.locator.SetInputData( wrapper )
            
