#!/bin/sh
# the next line restarts using wish \
exec wish "$0" "$@"


load [file dirname [info script]]/test_ScubaFrame[info sharedlibextension] Test_scubaframe

set gCurWindowID 0
proc GetNewWindowID { } {
    global gCurWindowID
    incr gCurWindowID
    return $gCurWindowID;
}

proc CreateWindow { } {

    set windowID [GetNewWindowID]

    set ww         .ww$windowID
    set fwTop      $ww.fwTop
    set twMain     $fwTop.twMain

    set fwControls  $ww.fwControls
    set bwNewWindow $fwControls.bwNewWindow
    set bw11       $fwControls.bw11
    set bw22       $fwControls.bw22
    set bw44       $fwControls.bw44
    set bw13       $fwControls.bw13

    toplevel $ww

    frame $fwTop
    togl $twMain -width 300 -height 300 -rgba true -ident $windowID

    frame $fwControls
    button $bwNewWindow -text "New Window" \
	-command CreateWindow

    button $bw11 -text "11" \
	-command "SetViewConfiguration $windowID c11"
    button $bw22 -text "22" \
	-command "SetViewConfiguration $windowID c22"
    button $bw44 -text "44" \
	-command "SetViewConfiguration $windowID c44"
    button $bw13 -text "13" \
	-command "SetViewConfiguration $windowID c13"


    bind $twMain <Motion> "%W MouseMotionCallback %x %y %b"
    bind $twMain <ButtonPress> "%W MouseDownCallback %x %y %b"
    bind $twMain <ButtonRelease> "%W MouseUpCallback %x %y %b"
    bind $twMain <KeyRelease> "%W KeyUpCallback %x %y %K"
    bind $twMain <KeyPress> "%W KeyDownCallback %x %y %K"
    bind $twMain <Enter> "focus $twMain"

    pack $twMain -fill both -expand true
    pack $bwNewWindow $bw11 $bw22 $bw44 $bw13 -side left

    pack $fwTop -fill both -expand true
    pack $fwControls -fill x -expand true

    puts "tcl: created window $windowID"
}


CreateWindow


wm withdraw .
