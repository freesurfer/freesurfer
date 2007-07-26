

import Globals
import ImagePipeline
import Spline
import CrossHair

from Region import Region
from Laplace import Locator, WrapperPolyData

import Tkinter
import vtk
from vtk.tk.vtkTkRenderWindowInteractor import \
     vtkTkRenderWindowInteractor

from DlgRegionAddInfo import DlgRegionAddInfo
from DlgProfileInspector import DlgProfileInspector
from DlgRegionHighlights import DlgRegionHighlights

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


##################################

class ModeInterface:
    def Start(self):
        pass
    def Finish(self):
        pass


class RegionAwareMode(ModeInterface):
    """ generic mode that keeps track of the displayed regions """

    def __init__(self,caller):
        self.caller = caller
        caller.imageObservers.append(self.UpdateObservers)
        self.locator = None

    def SelectRegion(self,a,b):
        pass

    def Start(self):
        istyle = vtk.vtkInteractorStyleImage()
        istyle.AddObserver("LeftButtonPressEvent",self.SelectRegion)
        iren = Globals.renWin.GetInteractor()
        iren.SetInteractorStyle(istyle)

        self.UpdateObservers()

    def Finish(self):
        self.caller.imageObservers.remove(self.UpdateObservers)

    def UpdateObservers(self):
        " callback for when the slice is changed "
        
        # use the visible regions
        regions = Globals.objectSet.visibleRegions

        if len(regions)==0:
            self.locator = None
            return

        self.locator = Locator()
        wrapper = WrapperPolyData()

        counter = 0
        for r in regions:
            l = r.GetPoints()
            startCounter = counter
            for pt in l:
                wrapper.InsertPoint( counter, pt[0], pt[1], pt[2])
                counter += 1
            wrapper.InsertNextCell( len(l) )
            for i in range(startCounter,counter): wrapper.InsertCellPoint(i)

        self.locator.SetInputData( wrapper )

    def GetSelectedRegion(self,a,b):
        if self.locator == None: return None

        wp = GetEventWorldPoint()
        regionId = self.locator.FindClosestPoint( \
            wp[0], wp[1], wp[2] )

        return regionId

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

class TestNavigation:
    """ default mode
    also base class for all other modes"""

    def Start(self):
        istyle = vtk.vtkInteractorStyleImage()
        istyle.AddObserver("LeftButtonPressEvent",self.Left)
        istyle.AddObserver("RightButtonPressEvent", self.Right)
        iren = Globals.renWin.GetInteractor()
        iren.SetInteractorStyle(istyle)

        self.spline = Spline.PlaneSpline()
        self.renSpline = Spline.RenderSpline(1)
        self.splineExists = 0

    def Left(self,a,b):

        if not self.splineExists: self.AddPoint()
        else: self.SelectPoint()

    def AddPoint(self):
        world = GetEventWorldPoint()
        self.spline.AddPoint(world)
        print " number of points = " , len(self.spline.points)
        self.renSpline.SetInput( self.spline )
        Globals.renWin.Render()

    def SelectPoint(self):
        print " selecting point "
        a = self.locator.FindClosestPoint( GetEventWorldPoint() )
        print " selected point ", a
        self.spline.SetPointPosition( a, GetEventWorldPoint())
        self.renSpline.SetPoints( self.spline )
        Globals.renWin.Render()

    def Right(self,a,b):
        self.splineExists = 1

        # build locator
        self.locator = vtk.vtkPointLocator()
        self.locator.SetDataSet( self.spline.GetVtkPolyData() )    


class AddRegion(ModeInterface):

    def __init__(self,caller):
        caller.imageObservers.append(self.UpdateObserver)
        self.region = Region()
        self.renSpline = Spline.RenderSpline()
        self.InitCurves()

        self.caller = caller
    def Finish(self):
        self.caller.imageObservers.remove(self.UpdateObserver)
        self.dlg.cancel()

    def UpdateObserver(self):
        self.Reset()
        self.dlg.reset()
        self.dlg.SetRegionInformation(Globals.imagePipeline.GetSliceIndex(),
                                      Globals.imagePipeline.GetOrientation()
                                      )
        Globals.renWin.Render()
    def Start(self):
        istyle = vtk.vtkInteractorStyleImage()
        istyle.AddObserver("LeftButtonPressEvent",self.AddPoint)
        iren = Globals.renWin.GetInteractor()
        iren.SetInteractorStyle(istyle)

        self.dlg = DlgRegionAddInfo(Globals.root,
                                    self.Next,
                                    self.Reset,
                                    self.Del)
        self.dlg.SetRegionInformation(Globals.imagePipeline.GetSliceIndex(),
                                      Globals.imagePipeline.GetOrientation()
                                      )

        # define a state list
        self.stateChain = [ self.AddZeroCurve,
                            self.AddOneCurve
                            ]
        self.stateIndex = 0

    def Next(self):
        if self.stateChain[self.stateIndex]():
            # pass onto next state
            self.stateIndex = ( 1 + self.stateIndex ) % len(self.stateChain)
            return 1
        return 0

    def Reset(self):
        self.stateIndex = 0
        self.region.Hide()
        self.region = Region()
        self.InitCurves()
        Globals.renWin.Render()

    def Del(self):
        pts = self.spline.points
        del pts[-1]
        self.spline = Spline.PlaneSpline(pts)
        self.renSpline.SetInput(self.spline)
        Globals.renWin.Render()

    def AddPoint(self,a,b):
        """ adds a point to the current spline """
        world = GetEventWorldPoint()
        self.spline.AddPoint(world)
        self.renSpline.SetInput( self.spline )
        Globals.renWin.Render()

    def AddZeroCurve(self):
        if len(self.spline.points) > 0 :
            self.region.SetZeroPoints(self.spline)
            self.InitCurves()
            Globals.renWin.Render()
            return 1
        return 0

    def AddOneCurve(self):
        if len(self.spline.points) > 0:
            self.region.SetOnePoints(self.spline)
            self.InitCurves()
            Globals.renWin.Render()

            # commit the region and reinit it
            Globals.objectSet.AddRegion( self.region )
            self.region = Region()
            return 1
        return 0

    def InitCurves(self):
        self.renSpline.Hide()
        self.spline = Spline.PlaneSpline()
        self.renSpline = Spline.RenderSpline(1)


class Laplace(RegionAwareMode):
    
    def SelectRegion(self,a,b):

        regionId = self.GetSelectedRegion(a,b)
        if regionId == None: return
        region = Globals.objectSet.visibleRegions[regionId]

        region.ApplyLaplace()
        # update profiles

        region.ComputeProfiles()
        region.ShowProfiles()
        Globals.renWin.Render()
        

class ProfileInspector(ModeInterface):
    
    def __init__(self,caller):
        self.caller = caller
        self.dlg = DlgProfileInspector(Globals.root,
                                       self.SelectProfile,
                                       self.SelectPoint )

        self.pointId = None
        self.glyph = None
        self.profile = None
        self.region = None
        
    def Start(self):
        self.caller.imageObservers.append(self.UpdateObservers)
        istyle = vtk.vtkInteractorStyleImage()
        istyle.AddObserver("LeftButtonPressEvent",self.ClickRegion)
        iren = Globals.renWin.GetInteractor()
        iren.SetInteractorStyle(istyle)
        
        self.dic = {}
        self.BuildCellLocator()
        
        
    def Finish(self):
        self.caller.imageObservers.remove(self.UpdateObservers)
        self.UpdateObservers()
        self.dlg.cancel()
        Globals.renWin.Render()
        
    def ClickRegion(self,a,b):

        if not self.locator: return
        
        wp = GetEventWorldPoint()
        cellId = self.locator.FindClosestPoint(
            wp[0],wp[1],wp[2] )
        elts = self.dic[cellId]
        region = elts[0]
        profileId = elts[1]
        ptId = self.GetClosestProfilePoint(
            wp,
            region.profiles[profileId] )
        self.dlg.SetSelectedRegion(region,
                                   profileId,
                                   ptId)
        self.SelectRegion(region)
        self.SelectProfile(profileId)
        self.SelectPoint(ptId)
        
        # update dlg info
        self.dlg.SetSelectedRegion(self.region,
                                   profileId, ptId)
        Globals.renWin.Render()

    def UpdateObservers(self):
        
        # rebuild locator
        self.dic = {}
        self.BuildCellLocator()
        
        # clear selected data and restore default states
        self.ClearSelectedPoint()
        self.ClearSelectedProfile()
        self.ClearSelectedRegion()
            
        # clear GUI data
        self.dlg.SetSelectedRegion()

    def ClearSelectedPoint(self):
        if self.glyph:
            self.glyph.Hide()
            self.glyph = None
            self.pointId = None
    def ClearSelectedProfile(self):
        if self.profile:
            self.profile.profile.GetProperty().SetColor(1,1,1)
            self.profile = None
    def ClearSelectedRegion(self):
        if self.region:
            self.region.rendering.SetDefaultColors()
            self.region = None

    def SelectPoint(self,ptId):
        self.ClearSelectedPoint()
        if self.region:
            if not self.glyph:
                self.glyph = CrossHair.CrossHairWidget(1,
                                                       Globals.renderProps["sphereSize"])
            self.pointId = ptId
            pt = self.profile.handles.points[ptId]
            self.glyph.SetPosition(pt[0],
                                   pt[1])
        Globals.renWin.Render()
    def SelectProfile(self,profile):
        self.ClearSelectedProfile()
        if self.region:
            self.profile = self.region.rendering.profiles[
                profile]
            self.profile.profile.GetProperty().SetColor(0,1,0)
            if self.pointId:
                self.SelectPoint(self.pointId)
        Globals.renWin.Render()
    def SelectRegion(self,region):
        self.ClearSelectedPoint()
        self.ClearSelectedProfile()
        self.ClearSelectedRegion()
        self.region = region
        self.region.rendering.SetSelected()

    def GetClosestProfilePoint(self,worldPoint,
                               spline):
        pointLocator = vtk.vtkPointLocator()
        pointLocator.SetDataSet( spline.GetVtkPolyData() )
        
        return pointLocator.FindClosestPoint( worldPoint )

    def BuildCellLocator(self):

        regions = Globals.objectSet.visibleRegions
        if not len(regions):
            self.locator = None
            return

        self.locator = Locator()
        wrapper = WrapperPolyData()

        pointCounter = 0
        cellCounter = 0
        for r in regions:
            localCounter = 0
            for p in r.profiles:
                startCounter = cellCounter
                for pt in p.points:
                    wrapper.InsertPoint( pointCounter,
                                         pt[0],pt[1],pt[2] )
                    pointCounter += 1
                wrapper.InsertNextCell( len(p.points) )
                for i in range(startCounter, pointCounter):
                    wrapper.InsertCellPoint(i)
                item = [r,localCounter]
                self.dic[cellCounter] = item
                cellCounter += 1
                localCounter += 1

        self.locator.SetInputData( wrapper )

class Highlight(RegionAwareMode):

    def __init__(self,caller):
        RegionAwareMode.__init__(self,caller)

        self.dlg = DlgRegionHighlights(Globals.root,
                                       self.Update,
                                       self.GoToRegion)

        self.region = None
        self.profileIndices = None
        self.renProfiles = None
        self.caller = caller

    def Finish(self):
        RegionAwareMode.Finish(self)
        #  class specific operations
        #     - cancel the dialog
        #     - destroy special glyphs, if any
        self.dlg.cancel()
        self.ClearProfiles()
        Globals.renWin.Render()

    def UpdateObservers(self):
        """ adds functionality to that of the basic RegionAware mode
        current region needs to be deselected when changing """

        RegionAwareMode.UpdateObservers(self)

        self.ClearProfiles()
        self.dlg.ClearData()
        #print " nullifying region"
        self.region = None

            
    def ClearProfiles(self):

        if self.renProfiles:
            for rp in self.renProfiles:
                rp.profile.GetProperty().SetColor(1,1,1)
            self.renProfiles = None

    def ShowProfiles(self):
        """ highlights profiles associated with the current region """

        self.ClearProfiles()
        if self.region and self.profileIndices:
            self.renProfiles = [ self.region.rendering.profiles[p] \
                                 for p in self.profileIndices ]
            for rp in self.renProfiles:
                rp.profile.GetProperty().SetColor(
                    Globals.preferences["highlightRed"],
                    Globals.preferences["highlightGreen"],
                    Globals.preferences["highlightBlue"] )
        Globals.renWin.Render()
        
    def GoToRegion(self,regionId):
        """ will change the slice for a certain region"""

        # search through the list of regions and see if the value exists
        allRegions = Globals.objectSet.GetAllRegions()
        region = None
        for r in allRegions:
            if r.index == regionId:
                region = r
                break
        if not region: return

        # set the view accordingly
        # to prevent dialog data from being cleared, set no update flag
        self.caller.SetOrientation( region.orientation )
        self.caller.SetSlice( region.sliceIndex )
        self.caller.scaleSlice.set(region.sliceIndex)
        print " setting region"
        #self.region = region
        #self.dlg.SetRegion(regionId)


    def SelectRegion(self,a,b):

        regionId = self.GetSelectedRegion(a,b)
        if regionId == None: return

        self.region = Globals.objectSet.visibleRegions[regionId]
        self.dlg.SetRegion(self.region.index)

    def Update(self,indices):
        """ will do the actual highlighting of the profiles"""

        # todo - should this be here?
        if self.region==None: return
        
        # check to see if the indices are valid
        tmpIndices = [int(i) for i in indices]
        self.profileIndices = [ pidx for pidx in tmpIndices
                                if pidx>=0 and pidx<len(self.region.profiles) ]
        self.ShowProfiles()

class DeleteRegion(RegionAwareMode):

    def SelectRegion(self,a,b):
        print " got here"
        regionId = self.GetSelectedRegion(a,b)
        if regionId == None: return
        
        Globals.objectSet.DeleteVisibleRegion(regionId)
        Globals.renWin.Render()
        self.UpdateObservers()
        
