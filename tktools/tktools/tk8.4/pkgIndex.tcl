if {[package vcompare [package provide Tcl] 8.4] != 0} { return }
package ifneeded Tk 8.4 [list load [file join $dir .. libtk8.4.so] Tk]
