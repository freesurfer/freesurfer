
import vtk

import Laplace

from Spline import *
import Globals
from ImagePipeline import ImagePipeline

class RegionRendering:

    def __init__(self):

        self.zero = RenderSpline(1)
        self.one = RenderSpline(1)
        self.cBegin = RenderSpline(0)
        self.cEnd = RenderSpline(0)

        self.SetDefaultColors()
        self.profiles = []

    def UpdateZero(self,pzero):
        self.zero.SetInput(pzero)
        self.UpdateEnds()

    def UpdateOne(self,pone):
        self.one.SetInput(pone)
        self.UpdateEnds()

    def HideBounds(self):
        for rs in (self.zero,self.one,self.cBegin,self.cEnd):
            rs.Hide()
    def HideProfiles(self):
        for rs in self.profiles:
            rs.Hide()
    def HideAll(self):
        self.HideBounds()
        self.HideProfiles()
    def ShowBounds(self):
        for rs in (self.zero,self.one,self.cBegin,self.cEnd):
            rs.Show()
    def ShowProfiles(self):
        for rs in self.profiles:
            rs.Show()
    def ShowAll(self):
        self.ShowBounds()
        self.ShowProfiles()

    def UpdateProfiles(self, lstProfiles):
        for p in lstProfiles:
            self.profiles.append(RenderSpline(0,p))

    def SetSelected(self):
        self.zero.profile.GetProperty().SetColor(.5,0,0)
        self.one.profile.GetProperty().SetColor(0,0,.5)
        for rs in (self.cBegin,self.cEnd):
            rs.profile.GetProperty().SetColor(0,.5,0)

    def SetDefaultColors(self):
        self.zero.profile.GetProperty().SetColor(1,0,0)
        self.one.profile.GetProperty().SetColor(0,0,1)
        for rs in (self.cBegin,self.cEnd):
            rs.profile.GetProperty().SetColor(0,1,0)

    def ClearProfiles(self):
        self.profiles = []

    def UpdateEnds(self):
        if not ( len(self.zero.handles.points) and
                 len(self.one.handles.points) ): return
        
        s = PlaneSpline()
        s.AddPoint( self.zero.handles.points[0] )
        s.AddPoint( self.one.handles.points[0] )

        self.cBegin.SetInput( s )

        s = PlaneSpline()
        s.AddPoint( self.zero.handles.points[-1] )
        s.AddPoint( self.one.handles.points[-1] )

        self.cEnd.SetInput( s )

    def GetRenderedSplines(self):
        splines = self.profiles[:]
        splines.extend([self.zero,self.one,self.cBegin,self.cEnd])
        return splines

    def UpdateProperties(self):
        for s in self.GetRenderedSplines():
            s.UpdateProperties()
    
class Region:

    def __init__(self,zero=PlaneSpline(),
                 one=PlaneSpline()):

        self.index = -1
        self.sliceIndex = Globals.imagePipeline.GetSliceIndex()
        self.orientation = Globals.imagePipeline.GetOrientation()

        self.zero = zero
        self.one = one
        self.profiles = []

        # Laplace structures
        self.tracer = None
        self.profileSpacing = None
        self.numberOfPoints = None

        self.rendering = RegionRendering()
        if len(self.zero.points) and len(self.one.points):
            self.FixOrientation()
            self.rendering.UpdateZero(self.zero)
            self.rendering.UpdateOne(self.one)

    def FixOrientation(self):
        # fix orientation if necessary
        if PointDistanceSqr(self.zero.points[0],
                            self.one.points[0]) > \
                            PointDistanceSqr(self.zero.points[0],
                                             self.one.points[-1]):
            self.one.points.reverse()

    def GetPoints(self):
        l = self.zero.points[:]
        l.extend( self.one.points )
        return l

    def SetZeroPoints(self,zero):
        self.zero = zero
        # give the member as argument - orientation
        self.rendering.UpdateZero(self.zero)
    def SetOnePoints(self,one):
        self.one = one
        self.FixOrientation()
        # give the member as argument - orientation
        self.rendering.UpdateOne(self.one)
    def Render(self):
        self.rendering.ShowBounds()
        Globals.renWin.Render()

    def Hide(self):
        self.rendering.HideAll()
    def Show(self):
        self.rendering.ShowAll()

    def ApplyLaplace(self):

        resolution = Globals.preferences["laplaceResolution"] \
                     * Globals.referenceSize

        dist = resolution / 3.0
        # resample the data to fit the resolution
        resampleZero = PlaneSpline(self.zero.ComputeRepresentation(dist))
        resampleOne = PlaneSpline(self.one.ComputeRepresentation(dist))

        bufBegin = PlaneSpline()
        bufBegin.AddPoint( resampleZero.points[0] )
        bufBegin.AddPoint( resampleOne.points[0] )
        resampleBegin = PlaneSpline(bufBegin.ComputeRepresentation(dist))

        bufEnd = PlaneSpline()
        bufEnd.AddPoint( resampleZero.points[-1] )
        bufEnd.AddPoint( resampleOne.points[-1] )
        resampleEnd = PlaneSpline(bufEnd.ComputeRepresentation(dist))

        splines = [ resampleZero, resampleOne, resampleBegin, resampleEnd ]
        # build the polydata
        wrapper = Laplace.WrapperPolyData()
        counter = 0
        for spline in splines:
            startCounter = counter
            for pt in spline.points:
                wrapper.InsertPoint(counter,
                                    pt[0], pt[1], pt[2] )
                counter += 1
            wrapper.InsertNextCell( len(spline.points) )
            for i in range(startCounter, counter): wrapper.InsertCellPoint(i)

        # update the tracer
        self.tracer = Laplace.ComputeSolution( wrapper,
                                               10, 10,
                                               resolution,
                                               int(Globals.preferences["laplaceConvergence"])
                                               )

    def ComputeProfiles(self):

        if not self.tracer: return
        self.profiles = []
        dv = Laplace.DoubleVector()
        dv.append(.5)

        
        vml = Laplace.ComputeIsolines( dv, self.tracer,
                                       self.zero.points[0][0],
                                       self.zero.points[0][1] )

        profiles = Laplace.ComputeProfiles( 10, \
                                            Globals.preferences["profileSpacing"] \
                                            * Globals.referenceSize,
                                            vml[0], self.tracer )

        for line in profiles:
            bufSpline = PlaneSpline()
            for pt in line:
                pt.append(0)
                bufSpline.AddPoint(pt)
            newSpline = PlaneSpline( \
                bufSpline.Resample( Globals.preferences["profilePoints"] ) )
            self.profiles.append(newSpline)

    def ShowProfiles(self):
        self.rendering.UpdateProfiles(self.profiles)
        Globals.renWin.Render()


def CommitSplineRenderingProperties():
    regions = Globals.objectSet.GetAllRegions()
    for r in regions:
        r.rendering.UpdateProperties()

    Globals.renWin.Render()


class RegionSet:

    def __init__(self):
        self.dic = {}
        self.visibleRegions = []
        self.nextIndex = 0

    def AddRegion(self,region,index=None):
        if index:
            region.index = index
            self.nextIndex = max(index+1, self.nextIndex)
        else:
            region.index = self.nextIndex
            self.nextIndex += 1

        if self.dic.has_key(region.sliceIndex):
            self.dic[region.sliceIndex].append(region)
        else:
            self.dic[region.sliceIndex] = [region]

        self.visibleRegions.append(region)

    def UpdateFromGUI(self):
        """ only display regions matching the current position """

        for r in self.visibleRegions:
            r.Hide()

        si = Globals.imagePipeline.GetSliceIndex()
        ori = Globals.imagePipeline.GetOrientation()

        regions = self.GetRegions(si,ori)
        self.visibleRegions = regions
        
        for r in regions:
            r.Show()

    def GetRegions(self,sliceIndex,ori):

        regions = []

        if self.dic.has_key(sliceIndex):
            for r in self.dic[sliceIndex]:
                if r.orientation == ori:
                    regions.append(r)

        return regions

    def GetAllRegions(self):
        regions = []
        for k in self.dic:
            rlst = self.dic[k]
            regions.extend(rlst)
        return regions

    def GetIndices(self,si,ori):
        kSmaller = None
        kLarger = None
        # create a list of available indices
        for k in self.dic.keys():
            found = 0
            index = 0
            regions = self.dic[k]
            while (not found) and (index<len(regions)):
                if regions[index].orientation == ori:
                    found = 1
                index += 1
            if found:
                if k>si:
                    if not kLarger: kLarger=k
                    elif k<kLarger: kLarger=k
                elif k<si:
                    if not kSmaller: kSmaller=k
                    elif k>kSmaller: kSmaller=k
                    
        return kSmaller,kLarger

    def UpdateAllProfiles(self):
        regions = self.GetAllRegions()
        
        for r in regions:
            r.rendering.ClearProfiles()
            r.ComputeProfiles()
            r.rendering.UpdateProfiles(r.profiles)
            r.Hide()
        for r in self.visibleRegions:
            r.Show()
    
    def DeleteVisibleRegion(self,regionId):
        region = self.visibleRegions[regionId]
        if region == None: return

        sliceIndex = region.sliceIndex
        region.Hide()
        self.dic[sliceIndex].remove(region)
        if len(self.dic[sliceIndex])==0:
            del self.dic[sliceIndex]
        self.UpdateFromGUI()

