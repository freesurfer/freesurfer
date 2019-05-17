import numpy as np
import vtk

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from vtk.util.numpy_support import numpy_to_vtk

from . import LookupTable


def view_surface(surface, colors=None, labels=None, lut=None, overlay=None, cmap='spectral'):
    # create surface mesh
    mesh = poly(surface.vertices, surface.faces)

    if labels is not None:
        if lut is None:
            lut = LookupTable.from_list(labels)
        colors = [lut.color(l) for l in labels]

    if colors is not None:
        if not isinstance(colors, vtk.vtkArray):
            colors = numpy_to_vtk(colors, deep=1, array_type=3)
        mesh.GetPointData().SetScalars(colors)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(mesh)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    view([actor])


def view(actors, background=(1, 1, 1)):
    # create a renderer, render window, and an interactor
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(background)
    for actor in actors:
        renderer.AddActor(actor)

    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindow.SetWindowName('FreeSurfer Viewer')
    size = np.array(renderWindow.GetScreenSize()) * 0.75
    renderWindow.SetSize(size.astype('int'))

    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()

    renderWindow.Render()
    renderWindowInteractor.Start()


def poly(vertices, faces):
    # init points
    points = vtk.vtkPoints()
    for vertex in vertices:
        points.InsertNextPoint(vertex)

    # init faces
    triangles = vtk.vtkCellArray()
    for face in faces:
        tri = vtk.vtkTriangle()
        tri.GetPointIds().SetId(0, face[0])
        tri.GetPointIds().SetId(1, face[1])
        tri.GetPointIds().SetId(2, face[2])
        triangles.InsertNextCell(tri)

    # create a polydata object
    poly = vtk.vtkPolyData()
    poly.SetPoints(points)
    poly.SetPolys(triangles)
    return poly
