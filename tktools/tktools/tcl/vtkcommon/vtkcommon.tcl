package require -exact vtkbase 5.6

if {[info commands ::vtk::init::require_package] != ""} {
  if {[::vtk::init::require_package vtkCommonTCL 5.6]} {
    package provide vtkcommon 5.6
  }
} else {
  puts stderr "Warning: Your TCLLIBPATH points to the VTK source tree. \
 Support for this is deprecated in VTK 4.2, and will be removed in a\
 future version.  Please point TCLLIBPATH to\
 {VTK_BINARY_DIR}/Wrapping/Tcl<config>, where <config>\
 is the build configuration.  The <config> value is empty for most\
 platforms, and /Debug, /Release, etc. for Visual Studio builds. \
 After this is done, you no longer need to set PATH or LD_LIBRARY_PATH\
 to point to the VTK build tree."
  if {[info commands vtkObject] != "" ||
    [::vtk::load_component vtkCommonTCL] == ""} {
    package provide vtkcommon 5.6
  }
}
