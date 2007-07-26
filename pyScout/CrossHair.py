
import vtk

import Globals

class CrossHairWidget:
    """ navigation cursor
     the coordinates of this actor will be in 2D, i.e. the last coordinate will always be zero"""

    def __init__(self, profile=0, radius=2):

        if not radius:
            self.radius = Globals.renderProps["sphereSize"] * 1.1
        else: self.radius = radius
        self.z = 0

        self.data = vtk.vtkPolyData()


        self.sphere = vtk.vtkSphereSource()
        self.sphere.SetRadius( Globals.referenceSize * self.radius )
        self.sphere.SetPhiResolution( Globals.renderProps["spherePhiResolution"] )
        self.sphere.SetThetaResolution( Globals.renderProps["sphereThetaResolution"] )

        self.glyphPoints = vtk.vtkGlyph3D()
        self.glyphPoints.SetInput( self.data )
        self.glyphPoints.SetSource( self.sphere.GetOutput() )

        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputConnection( self.glyphPoints.GetOutputPort() )

        self.actor = vtk.vtkActor()
        self.actor.SetMapper(self.mapper)
        if profile:
            self.actor.GetProperty().SetColor(0,0,1)
        else:
            self.actor.GetProperty().SetColor(1,1,0)

        Globals.ren.AddActor(self.actor)

    def SetPosition(self, x,y):

        self.x = x
        self.y = y

        point = vtk.vtkPoints()
        point.SetNumberOfPoints(1)
        point.SetPoint(0, x,y,0)
        self.data.SetPoints(point)

    def Hide(self):
        Globals.ren.RemoveActor(self.actor)
    def Show(self):
        Globals.ren.AddActor(self.actor)

    def UpdateProperties(self):
        self.radius = Globals.renderProps["sphereSize"] * 1.1
        self.sphere.SetRadius( Globals.referenceSize * self.radius )
        self.sphere.Update()
        self.glyphPoints.Modified()
        
