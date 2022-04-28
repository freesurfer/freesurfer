package require -exact vtkhybrid 5.6

if {[info commands ::vtk::init::require_package] != ""} {
  if {[::vtk::init::require_package vtkWidgetsTCL 5.6]} {
    package provide vtkwidgets 5.6
  }
} else {
  if {[info commands vtkBorderWidget] != "" ||
    [::vtk::load_component vtkWidgetsTCL] == ""} {
    package provide vtkwidgets 5.6
  }
}
