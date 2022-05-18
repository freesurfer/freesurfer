package require -exact vtkfiltering 5.6

if {[info commands ::vtk::init::require_package] != ""} {
  if {[::vtk::init::require_package vtkIOTCL 5.6]} {
    package provide vtkio 5.6
  }
} else {
  if {[info commands vtkBMPReader] != "" ||
    [::vtk::load_component vtkIOTCL] == ""} {
    package provide vtkio 5.6
  }
}
