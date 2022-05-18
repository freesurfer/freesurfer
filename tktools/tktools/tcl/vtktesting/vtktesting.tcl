package require -exact vtkcommon 5.6

foreach s {colors mccases backdrop grab} {
  source [file join [file dirname [info script]] "${s}.tcl"]
}

package provide vtktesting 5.6
