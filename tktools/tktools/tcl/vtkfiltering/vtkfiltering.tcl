package require -exact vtkcommon 5.6

if {[info commands ::vtk::init::require_package] != ""} {
  if {[::vtk::init::require_package vtkFilteringTCL 5.6]} {
    package provide vtkfiltering 5.6
  }
} else {
  if {[info commands vtkCardinalSpline] != "" ||
    [::vtk::load_component vtkFilteringTCL] == ""} {
    package provide vtkfiltering 5.6
  }
}
