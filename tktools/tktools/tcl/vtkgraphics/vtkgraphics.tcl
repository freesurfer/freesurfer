package require -exact vtkfiltering 5.6

if {[info commands ::vtk::init::require_package] != ""} {
  if {[::vtk::init::require_package vtkGraphicsTCL 5.6]} {
    package provide vtkgraphics 5.6
  }
} else {
  if {[info commands vtkAxes]  != "" ||
    [::vtk::load_component vtkGraphicsTCL] == ""} {
    package provide vtkgraphics 5.6
  }
}
