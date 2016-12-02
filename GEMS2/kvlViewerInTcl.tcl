###
### from /home/koen/software/VTK/Graphics/Testing/Tcl/TestRectilinearGridToTetrahedra.tcl
###
### and /home/koen/software/VTK/Examples/VisualizationAlgorithms/Tcl/ClipCow.tcl
###
### ./Medical3 /home/koen/software/VTKData/Data/headsq/quarter

package require vtk
package require vtkinteraction
package require vtktesting


set sagittalSliceNumber 32
set axialSliceNumber 46
set coronalSliceNumber 32

vtkMetaImageReader  imageReader
  imageReader SetFileName "/home/koen/software/Atlas3D-KDevelop/bin/alphaImage.mhd"
  imageReader Update

set dimensions [[imageReader GetOutput] GetDimensions]
set origin [[imageReader GetOutput] GetOrigin]
set spacing [[imageReader GetOutput] GetSpacing]

vtkMetaImageReader  overlayReader
  overlayReader SetFileName "/home/koen/software/Atlas3D-KDevelop/bin/overlay.mhd"
  overlayReader Update

# Read the grid
vtkDataSetReader  gridReader
  gridReader SetFileName "/home/koen/software/Atlas3D-KDevelop/bin/grid.vtk"
  gridReader Update
  [ gridReader GetOutput ] GetClassName


# An outline provides context around the data.
#
vtkOutlineFilter outlineData
  #outlineData SetInputConnection [imageReader GetOutputPort]
  #outlineData SetInputConnection [gridReader GetOutputPort]
  outlineData SetInput [gridReader GetOutput]
vtkPolyDataMapper mapOutline
  #mapOutline SetInputConnection [outlineData GetOutputPort ]
  mapOutline SetInput [outlineData GetOutput ]
vtkActor outline
  outline SetMapper mapOutline
  [outline GetProperty] SetColor 0 0 0

# Start by creating a black/white lookup table.
vtkLookupTable bwLut
  bwLut SetTableRange 0.0 1.0 
  bwLut SetSaturationRange 0 0
  bwLut SetHueRange 0 0
  bwLut SetValueRange 0 1
  bwLut Build 

# Also create a color lookup table for the overlay
vtkLookupTable  overlayLookup
   # input range of image; outside values will be clamped
   overlayLookup SetTableRange 0 255
   # specify color here
   overlayLookup SetHueRange 0.75 0.75
   # very colorful
   overlayLookup SetSaturationRange 1 1
   # very bright
   overlayLookup SetValueRange 1 1 
   overlayLookup SetAlphaRange 0 1
   overlayLookup Build



#
vtkImageMapToColors saggitalColors
  saggitalColors SetInputConnection [imageReader GetOutputPort]
  saggitalColors SetLookupTable bwLut
  saggitalColors PassAlphaToOutputOn
vtkImageMapToColors saggitalOverlayColors
  saggitalOverlayColors SetInputConnection [ overlayReader GetOutputPort ]
  saggitalOverlayColors SetLookupTable overlayLookup
  saggitalOverlayColors PassAlphaToOutputOn
vtkImageBlend  saggitalBlender
  saggitalBlender SetInput 0 [saggitalColors GetOutput]
  saggitalBlender SetInput 1 [saggitalOverlayColors GetOutput]
  saggitalBlender SetOpacity 1 0.8
vtkImageActor saggital
  saggital SetInput [saggitalBlender GetOutput]
  saggital SetDisplayExtent $sagittalSliceNumber $sagittalSliceNumber 0 [expr [lindex $dimensions 1] - 1 ] 0 [expr [lindex $dimensions 2] - 1 ]

# Create the second (axial) plane of the three planes. We use the
# same approach as before except that the extent differs.
vtkImageMapToColors axialColors
  axialColors SetInputConnection [imageReader GetOutputPort]
  axialColors SetLookupTable bwLut 
vtkImageMapToColors axialOverlayColors
  axialOverlayColors SetInputConnection [ overlayReader GetOutputPort ]
  axialOverlayColors SetLookupTable overlayLookup
  axialOverlayColors PassAlphaToOutputOn
vtkImageBlend  axialBlender
  axialBlender SetInput 0 [axialColors GetOutput]
  axialBlender SetInput 1 [axialOverlayColors GetOutput]
  axialBlender SetOpacity 1 0.8
vtkImageActor axial
  axial SetInput [axialBlender GetOutput]
  axial SetDisplayExtent 0 [expr [lindex $dimensions 0] - 1 ] 0 [expr [lindex $dimensions 1] - 1 ] $axialSliceNumber $axialSliceNumber 

# Create the third (coronal) plane of the three planes. We use 
# the same approach as before except that the extent differs.
vtkImageMapToColors  coronalColors
  coronalColors SetInputConnection [imageReader GetOutputPort]
  coronalColors SetLookupTable bwLut
vtkImageActor coronal
  coronal SetInput [coronalColors GetOutput]
  coronal SetDisplayExtent 0 [expr [lindex $dimensions 0] - 1 ] $coronalSliceNumber $coronalSliceNumber 0 [expr [lindex $dimensions 2] - 1 ] 

# It is convenient to create an initial view of the data. The
# FocalPoint and Position form a vector direction. Later on
# (ResetCamera() method) this vector is used to position the camera
# to look at the data in this direction.
vtkCamera aCamera
  aCamera SetViewUp 0 0 1
  #aCamera SetPosition 0 1 0
  aCamera SetPosition 1 1 1
  aCamera SetFocalPoint 0 0 0
  #aCamera ComputeViewPlaneNormal

vtkRenderer aRenderer
  aRenderer SetActiveCamera aCamera
  aRenderer AddActor outline
  aRenderer AddActor saggital
  aRenderer AddActor axial
  aRenderer AddActor coronal

# An initial camera view is created.  The Dolly() method moves 
# the camera towards the FocalPoint, thereby enlarging the image.
  #aRenderer SetActiveCamera aCamera
  #aRenderer Render
  #aRenderer ResetCamera
  #aCamera Dolly 1.5

# Set a background color for the renderer and set the size of the
# render window (expressed in pixels).
  aRenderer SetBackground 1 1 1

# Note that when camera movement occurs (as it does in the Dolly()
# method), the clipping planes often need adjusting. Clipping planes
# consist of two planes: near and far along the view direction. The 
# near plane clips out objects in front of the plane; the far plane
# clips out objects behind the plane. This way only what is drawn
# between the planes is actually rendered.
  #aRenderer ResetCameraClippingRange









# The (implicit) planes is used to do the cutting
vtkPlane  axialPlane
  axialPlane SetOrigin 0 0 [expr $axialSliceNumber * [lindex $spacing 2] ] 
  axialPlane SetNormal 0 0 1
vtkPlane  coronalPlane
  coronalPlane SetOrigin 0 [expr $coronalSliceNumber * [lindex $spacing 1] ] 0
  coronalPlane SetNormal 0 1 0
vtkPlane  sagittalPlane
  sagittalPlane SetOrigin [expr $sagittalSliceNumber * [lindex $spacing 0] ] 0 0
  sagittalPlane SetNormal 1 0 0

# vtkImplicitBoolean combinedPlanes
#     #combinedPlanes SetOperationTypeToIntersection
#     #combinedPlanes SetOperationTypeToUnionOfMagnitudes
#     combinedPlanes SetOperationTypeToUnion
#     #combinedPlanes SetOperationTypeToDifference
#     combinedPlanes AddFunction axialPlane
#     combinedPlanes AddFunction coronalPlane
#     combinedPlanes AddFunction sagittalPlane

vtkBox  box
  box SetBounds [expr $sagittalSliceNumber * [lindex $spacing 0] ] [expr [lindex $dimensions 0] * [lindex $spacing 0] ]      [expr $coronalSliceNumber * [lindex $spacing 1] ] [expr [lindex $dimensions 1] * [lindex $spacing 1] ]      [expr $axialSliceNumber * [lindex $spacing 2] ] [expr [lindex $dimensions 2] * [lindex $spacing 2] ]



# Extract edges, clip them, and make them cylindres
vtkExtractEdges  edgeExtracter
  edgeExtracter SetInput [ gridReader GetOutput ]


vtkClipPolyData  axialEdgeClipper
  #axialEdgeClipper SetInputConnection [ edgeExtracter GetOutputPort ]
  axialEdgeClipper SetInput [ edgeExtracter GetOutput ]
  axialEdgeClipper SetClipFunction axialPlane
  #clipper GenerateClipScalarsOn
  #clipper GenerateClippedOutputOn
  axialEdgeClipper SetValue 0.0
vtkTubeFilter axialEdgeTuber
   #axialEdgeTuber SetInputConnection [ axialEdgeClipper GetOutputPort ]
   axialEdgeTuber SetInput [ axialEdgeClipper GetOutput ]
   axialEdgeTuber SetRadius 0.1
   axialEdgeTuber SetNumberOfSides 6


# Extract triangles
vtkDataSetSurfaceFilter  triangleExtracter
#vtkGeometryFilter  meshSurfaces
  #meshSurfaces SetInputConnection 0 [FormMesh GetOutputPort]
  triangleExtracter SetInput [ gridReader GetOutput ]


vtkClipPolyData  axialTriangleClipper
  axialTriangleClipper SetInputConnection [triangleExtracter GetOutputPort]
  axialTriangleClipper SetClipFunction axialPlane
  #clipper GenerateClipScalarsOn
  #clipper GenerateClippedOutputOn
  axialTriangleClipper SetValue 0.0


# # Let's try to actually clip the unstructured grid
# vtkClipDataSet  clipper
#   clipper SetInput [ gridReader GetOutput ]
#   clipper SetClipFunction combinedPlanes
#   #clipper GenerateClipScalarsOn
#   #clipper GenerateClippedOutputOn
#   clipper SetValue 0.0
# vtkDataSetSurfaceFilter  clippedSurface
#  clippedSurface SetInputConnection 0 [clipper GetOutputPort]



vtkCutter axialCutter
  #cutter SetInput [ meshSurfaces GetOutput ]
  axialCutter SetInput [ gridReader GetOutput ]
  axialCutter SetCutFunction axialPlane
  #cutter GenerateCutScalarsOff
  #cutter SetSortByToSortByCell
#   cutter Update
#   set cutResult [cutter GetOutput]
#   eval $cutResult GetNumberOfPolys



# OpenGL stuff for the tethrahedral mesh
vtkPolyDataMapper axialEdgeMapper
  axialEdgeMapper SetInput [ axialEdgeTuber GetOutput]
vtkActor axialEdgeActor
  axialEdgeActor SetMapper axialEdgeMapper
  [axialEdgeActor GetProperty] SetColor 1 0 0
  [axialEdgeActor GetProperty] SetSpecularColor 1 1 1
  [axialEdgeActor GetProperty] SetSpecular 0.3
  [axialEdgeActor GetProperty] SetSpecularPower 20
  [axialEdgeActor GetProperty] SetAmbient 0.2
  [axialEdgeActor GetProperty] SetDiffuse 0.8

vtkPolyDataMapper axialTriangleMapper
  #mapSurface SetInput [clippedSurface GetOutput]
  axialTriangleMapper SetInput [axialTriangleClipper GetOutput]
vtkProperty backProp
  backProp SetDiffuseColor 0.5 0.5 0.5
vtkActor axialTriangleActor
  axialTriangleActor SetMapper axialTriangleMapper
  axialTriangleActor SetBackfaceProperty backProp
  [axialTriangleActor GetProperty] SetColor 0.7 0.4 0.4
  [axialTriangleActor GetProperty] SetSpecularColor 1 1 1
  [axialTriangleActor GetProperty] SetSpecular 0.3
  [axialTriangleActor GetProperty] SetSpecularPower 20
  [axialTriangleActor GetProperty] SetAmbient 0.2
  [axialTriangleActor GetProperty] SetDiffuse 0.8
  #[edgeActor GetProperty] SetRepresentationToWireframe




#aRenderer AddActor surfaceActor
## OpenGL stuff for the extracted plane
vtkPolyDataMapper axialCutMapper
 axialCutMapper SetInput [ axialCutter GetOutput]
vtkActor axialCutActor
  axialCutActor SetMapper axialCutMapper
 [axialCutActor GetProperty] SetColor 0 1 0
 [axialCutActor GetProperty] SetSpecularColor 1 1 1
 [axialCutActor GetProperty] SetSpecular 0.3
 [axialCutActor GetProperty] SetSpecularPower 20
 [axialCutActor GetProperty] SetAmbient 0.2
 [axialCutActor GetProperty] SetDiffuse 0.8



aRenderer AddActor axialEdgeActor
aRenderer AddActor axialCutActor


vtkRenderer ren
  ren AddActor axialEdgeActor
  ren AddActor axialTriangleActor
  ren SetBackground 0 0 0
  ren AddActor axialCutActor
#  ren ResetCamera
#  [ren GetActiveCamera] Zoom 1
#  [ren GetActiveCamera] SetPosition 1.73906 12.7987 -0.257808
#  [ren GetActiveCamera] SetViewUp 0.992444 0.00890284 -0.122379
#  [ren GetActiveCamera] SetClippingRange 9.36398 15.0496
#  ren ResetCamera
#  [ren GetActiveCamera] Zoom 1
#  [ren GetActiveCamera] SetPosition 1.0 -6.0 4.0
##  [ren GetActiveCamera] SetViewUp 0.992444 0.00890284 -0.122379
#  [ren GetActiveCamera] SetClippingRange 0 15.0496

  ren SetActiveCamera aCamera
  #aRenderer Render
  ren ResetCamera


#
vtkCamera  axialCamera
  axialCamera SetFocalPoint [expr ( [lindex $dimensions 0] -1 ) / 2.0 ] [expr ( [lindex $dimensions 1] -1 ) / 2.0 ] $axialSliceNumber
  axialCamera ParallelProjectionOn
  axialCamera SetParallelScale [expr ( [lindex $dimensions 0] -1 ) / 2.0 ];
  axialCamera SetPosition [expr ( [lindex $dimensions 0] -1 ) / 2.0 ] [expr ( [lindex $dimensions 1] -1 ) / 2.0 ] 0
  axialCamera SetViewUp 0 1 0
vtkRenderer  axialRenderer
  axialRenderer SetActiveCamera axialCamera
  axialRenderer AddActor axial
  #axialRenderer AddActor axialCutActor
  axialRenderer SetLayer 0
vtkRenderer  axialRenderer2
  axialRenderer2 SetActiveCamera axialCamera
  axialRenderer2 AddActor axialCutActor
  axialRenderer2 SetLayer 1


# Render window
vtkRenderWindow renWin
  renWin SetSize 1800 600
  renWin AddRenderer ren
  renWin AddRenderer aRenderer
  renWin SetNumberOfLayers 2
  renWin AddRenderer axialRenderer
  renWin AddRenderer axialRenderer2
  ren SetViewport 0 0 0.33 1
  aRenderer SetViewport 0.33 0 0.66 1
  axialRenderer SetViewport 0.66 0 1 1
  axialRenderer2 SetViewport 0.66 0 1 1
vtkRenderWindowInteractor iren
  set istyle [vtkInteractorStyleSwitch istyleswitch]
  iren SetInteractorStyle $istyle
  $istyle SetCurrentStyleToTrackballCamera
  iren SetRenderWindow renWin


proc Cut {axialSliceNumber coronalSliceNumber sagittalSliceNumber} {
  global spacing
  global dimensions
  axialPlane SetOrigin 0 0 [expr $axialSliceNumber * [lindex $spacing 2] ]
  coronalPlane SetOrigin 0 [expr $coronalSliceNumber * [lindex $spacing 1] ] 0
  sagittalPlane SetOrigin [expr $sagittalSliceNumber * [lindex $spacing 0] ] 0 0
  box SetBounds [expr $sagittalSliceNumber * [lindex $spacing 0] ] [expr [lindex $dimensions 0] * [lindex $spacing 0] ]      [expr $coronalSliceNumber * [lindex $spacing 1] ] [expr [lindex $dimensions 1] * [lindex $spacing 1] ]      [expr $axialSliceNumber * [lindex $spacing 2] ] [expr [lindex $dimensions 2] * [lindex $spacing 2] ]

  saggital SetDisplayExtent $sagittalSliceNumber $sagittalSliceNumber 0 [expr [lindex $dimensions 1] - 1 ] 0 [expr [lindex $dimensions 2] - 1 ]
  axial SetDisplayExtent 0 [expr [lindex $dimensions 0] - 1 ] 0 [expr [lindex $dimensions 1] - 1 ] $axialSliceNumber $axialSliceNumber 
  coronal SetDisplayExtent 0 [expr [lindex $dimensions 0] - 1 ] $coronalSliceNumber $coronalSliceNumber 0 [expr [lindex $dimensions 2] - 1 ] 

  renWin Render
}

# for {set i 1} {$i<1000} {incr i 1} {
#  set v [expr $i / 1000.0 * [lindex $dimensions 0] ]
#  Cut $v $v $v
# }
# for {set i 1000} {$i>0} {incr i -1} {
#  set v [expr $i / 1000.0 * [lindex $dimensions 0] ]
#  Cut $v $v $v
# }

# for {set i 0} {$i<99} {incr i 1} {
#  Cut $i 50 50
#  after 100
# }
# for {set i 99} {$i>0} {incr i -1} {
#  Cut $i 50 50
#  after 100
# }


Cut $axialSliceNumber $coronalSliceNumber $sagittalSliceNumber

iren AddObserver UserEvent {wm deiconify .vtkInteract}
iren Initialize
wm withdraw .

