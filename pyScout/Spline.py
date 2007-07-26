
import math
import vtk

import Globals

def PointDistanceSqr(p1,p2):
    dsum = 0.0
    for i in range(len(p1)):
        dsum += (p1[i]-p2[i])*(p1[i]-p2[i])
    return dsum

class PlaneSpline:
    """ 2D spline = ordered list of points """

    def __init__(self,points=None):
        self.points = []
        if points:
            self.points = points[:]

    def AddPoint(self,point):
        self.points.append(point)
    def SetPointPosition(self,id, newPosition):
        self.points[id] = newPosition

    def GetVtkPoints(self):
        vtkPoints = vtk.vtkPoints()
        vtkPoints.SetNumberOfPoints( len(self.points))
        count = 0
        for pt in self.points:
            vtkPoints.SetPoint(count,pt)
            count += 1
        return vtkPoints
    def GetVtkPolyData(self):
        vtkPoints = self.GetVtkPoints()
        line = vtk.vtkCellArray()
        line.InsertNextCell( len(self.points) )
        for i in range( len(self.points) ):
            line.InsertCellPoint(i)

        data = vtk.vtkPolyData()
        data.SetPoints(vtkPoints)
        data.SetLines(line)

        return data
    def ComputeRepresentation(self, maxDist):
        """ uses cardinal splines to compute compute an interpolated representation of the line set"""

        alen = self.GetApproximateLength()
        factor = math.ceil(alen/maxDist/float(len(self.points)))
        spoints = self.Oversample(factor )
        return spoints

    def __MeasureLen(self,spoints):

        def sqr(x):
            return x*x
        dsum = 0
        for i in range(1,len(spoints)):
            dsum += math.sqrt( sqr( spoints[i-1][0] - spoints[i][0] ) +
                               sqr( spoints[i-1][1] - spoints[i][1] )
                               )
        return dsum

    def GetApproximateLength(self):
        """ computes the approximate length of the SPLINE representation obtained from the set of points"""

        factor = 2
        spoints = self.Oversample( )
        crtLen = self.__MeasureLen(spoints)

        oldLen = crtLen *2.0
        while  math.fabs(oldLen-crtLen)/oldLen > .01:
            factor *= 2
            spoints = self.Oversample( factor )
            oldLen = crtLen
            crtLen = self.__MeasureLen(spoints)

        return crtLen

    def Resample(self, numberOfOutputPoints):

        sx = vtk.vtkCardinalSpline()
        sy = vtk.vtkCardinalSpline()
        counter = 0
        for pt in self.points:
            sx.AddPoint(counter, pt[0])
            sy.AddPoint(counter, pt[1])
            counter += 1
        numberOfInputPoints = len(self.points)

        ratio = (numberOfInputPoints-1.0) / (numberOfOutputPoints-1.0)
        spoints = []
        for i in range(numberOfOutputPoints):
            t = ratio * i
            pt = [ sx.Evaluate(t), sy.Evaluate(t), 0]
            spoints.append(pt)

        return spoints

    def Oversample(self, factor=2):
        return self.Resample( int(factor * len(self.points)) )

    def GetIntensities(self,selector):
        probe = vtk.vtkProbeFilter()
        volume = Globals.imagePipeline.volume
        
        m = vtk.vtkMatrix4x4()
        # populate the matrix
        m.DeepCopy( selector.GetDirectionCosines() )
        axesOrigin = selector.GetAxesOrigin()
        m.SetElement(0,3, axesOrigin[0])
        m.SetElement(1,3, axesOrigin[1])
        m.SetElement(2,3, axesOrigin[2])
        
        # use the selector to project points in the right spatial position
        volSpline = PlaneSpline()
        for pt in self.points:
            wpt = m.MultiplyPoint( [pt[0],pt[1],0,1] )
            volSpline.AddPoint(wpt[0:3])

        polyData = volSpline.GetVtkPolyData()
        probe.SetInput(polyData)
        probe.SetSource(volume)
        probe.Update()
        
        lstValues = []
        vtkValues = probe.GetOutput().GetPointData().GetScalars()
        for i in range(vtkValues.GetNumberOfTuples()):
            lstValues.append( vtkValues.GetComponent(i,0) )
            
        return str(lstValues)[1:-1]


class RenderSpline:
    """ handles the rendering of a spline in the VTK window """

    def __del__(self):
        self.Hide()

    def __init__(self, showHandles=0, planeSpline=None):
        self.showHandles = showHandles

        self.glyph = None
        if self.showHandles: self.SetupGlyph()
        self.SetupProfile()

        self.handles = PlaneSpline()
        if planeSpline: self.SetInput(planeSpline)

        self.isVisible = 1

    def UpdateProperties(self):
        
        self.tubes.SetRadius( Globals.renderProps["tubeSize"]
                              * Globals.referenceSize )
        self.profile.GetProperty().SetOpacity( \
            Globals.renderProps["opacity"] )
        if self.showHandles:
            self.sphere.SetRadius( \
                Globals.renderProps["sphereSize"]
                * Globals.referenceSize )
            self.sphere.Update()
            self.glyphPoints.Modified()
            self.glyph.GetProperty().SetOpacity( \
                Globals.renderProps["opacity"] )

    def SetupGlyph(self):
        self.handleData = vtk.vtkPolyData()
        self.sphere = vtk.vtkSphereSource()
        self.sphere.SetRadius( \
            Globals.renderProps["sphereSize"]
            * Globals.referenceSize )
        self.sphere.SetPhiResolution( \
            Globals.renderProps["spherePhiResolution"] )
        self.sphere.SetThetaResolution( Globals.renderProps["sphereThetaResolution"])

        self.glyphPoints = vtk.vtkGlyph3D()
        self.glyphPoints.SetInput( self.handleData )
        self.glyphPoints.SetSource( self.sphere.GetOutput() )

        self.glyphMapper = vtk.vtkPolyDataMapper()
        self.glyphMapper.SetInputConnection( self.glyphPoints.GetOutputPort() )

        self.glyph = vtk.vtkActor()
        self.glyph.SetMapper( self.glyphMapper )

        Globals.ren.AddActor( self.glyph )
        self.glyph.GetProperty().SetOpacity( \
            Globals.renderProps["opacity"] )

    def SetupProfile(self):
        self.profileData = vtk.vtkPolyData()
        
        self.tubes = vtk.vtkTubeFilter()
        self.tubes.SetNumberOfSides( Globals.renderProps["profileNumberOfSides"] )
        self.tubes.SetInput( self.profileData )
        self.tubes.SetRadius( \
            Globals.renderProps["tubeSize"]
            * Globals.referenceSize )

        self.profileMapper = vtk.vtkPolyDataMapper()
        self.profileMapper.SetInputConnection( self.tubes.GetOutputPort() )

        self.profile = vtk.vtkActor()
        self.profile.SetMapper( self.profileMapper )

        Globals.ren.AddActor( self.profile )
        self.profile.GetProperty().SetOpacity( \
            Globals.renderProps["opacity"] )

    def SetInput(self, planeSpline):
        self.handles = planeSpline

        if self.showHandles:
            tmpHandleData = self.handles.GetVtkPolyData()
            self.handleData.SetPoints( tmpHandleData.GetPoints() )
            self.handleData.SetLines(  tmpHandleData.GetLines() )

        if len(self.handles.points) > 1:
            profileSpline = PlaneSpline(self.handles.ComputeRepresentation(.5))
            tmpProfileData = profileSpline.GetVtkPolyData()
            #pdb.set_trace()
            self.profileData.SetPoints( tmpProfileData.GetPoints() )
            self.profileData.SetLines(  tmpProfileData.GetLines() )

    def Show(self):
        if self.isVisible: return

        if self.glyph: Globals.ren.AddActor(self.glyph)
        if self.profile: Globals.ren.AddActor(self.profile)
        self.isVisible = 1

    def Hide(self):
        
        if not self.isVisible: return

        if self.glyph: Globals.ren.RemoveActor(self.glyph)
        if self.profile: Globals.ren.RemoveActor(self.profile)
        self.isVisible = 0


