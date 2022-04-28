package require -exact vtkfiltering 5.6

if {[info commands ::vtk::init::require_package] != ""} {
  if {[::vtk::init::require_package vtkImagingTCL 5.6]} {
    package provide vtkimaging 5.6
  }
} else {
  if {[info commands vtkImageFFT] != "" ||
    [::vtk::load_component vtkImagingTCL] == ""} {
    package provide vtkimaging 5.6
  }
}
