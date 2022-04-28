# -*-mode: tcl; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: CObjView.tcl,v 1.2.2.1 2001/11/04 05:10:08 idiscovery Exp $
#
# Tix Demostration Program
#
# This sample program is structured in such a way so that it can be
# executed from the Tix demo program "widget": it must have a
# procedure called "RunSample". It should also have the "if" statment
# at the end of this file so that it can be run as a standalone
# program using tixwish.

# This program demonstrates the use of the CObjView (Canvas Object
# View) class.
#


proc RunSample {w} {
    label $w.lab  -justify left -text \
"Click on the buttons to add or delete canvas
objects randomally. Notice the scrollbars automatically
adjust to include all objects in the scroll-region."

    pack $w.lab -anchor c -padx 10 -pady 6
    tixCObjView $w.c
    pack $w.c -expand yes -fill both -padx 4 -pady 2
    button $w.add -command "CVDemo_Add $w.c"    -text Add    -width 6
    button $w.del -command "CVDemo_Delete $w.c" -text Delete -width 6
    pack $w.add $w.del -side left -padx 20 -pady 10 -anchor c -expand yes
}

set cvdemo_counter 0
proc CVDemo_Add {cov} {
    global cvdemo_counter

    # Generate four pseudo random numbers (x,y,w,h) to define the coordinates
    # of a rectangle object in the canvas.
    #
    set colors {red green blue white black gray yellow}

    set px [expr [lindex [time update] 0] + $cvdemo_counter]
    set py [expr [lindex [time update] 0] + $cvdemo_counter]
    set pw [expr [lindex [time update] 0] + $cvdemo_counter]
    set ph [expr [lindex [time update] 0] + $cvdemo_counter]
    set pc [expr [lindex [time update] 0] + $cvdemo_counter]

    set x [expr (20 - ($px % 37))   * 10]
    set y [expr (10 - ($py % 23))  * 10]
    set w [expr ($pw % 17)  * 10]
    set h [expr ($ph % 17)  * 10]

    # Create the canvas object
    #
    $cov subwidget canvas create rectangle $x $y [expr $x+$w] [expr $y+$h] \
	-fill [lindex $colors [expr $pc % [llength $colors]]]

    # Call the adjustscrollregion command to set the scroll bars so that all
    # objects are included in the scroll-region
    #
    $cov adjustscrollregion

    # This number acts as the seed for the next round of randomization.
    #
    set cvdemo_counter [expr ($px % 37)]
}

proc CVDemo_Delete {cov} {
    set px [lindex [time update] 0]
    set w [$cov subwidget canvas]
    set items [$w find withtag all]

    if [string compare $items ""] {
	# There are items in the canvas, randomally delete one of them
	# and re-adjust the scroll-region
	#
	set toDelete [expr $px % [llength $items]]
	$w delete [lindex $items $toDelete]

	$cov adjustscrollregion
    }
}

if {![info exists tix_demo_running]} {
    wm withdraw .
    set w .demo
    toplevel $w; wm transient $w ""
    RunSample $w
    bind $w <Destroy> exit
}
