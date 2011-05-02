#!/bin/sh
# the next line restarts using wish \
exec wish "$0" "$@"

##
## test_ToglManager.tcl
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2011/05/02 21:02:44 $
##    $Revision: 1.4 $
##
## Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
##
## Terms and conditions for use, reproduction, distribution and contribution
## are found in the 'FreeSurfer Software License Agreement' contained
## in the file 'LICENSE' found in the FreeSurfer distribution, and here:
##
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
##
## Reporting: freesurfer@nmr.mgh.harvard.edu
##
## General inquiries: freesurfer@nmr.mgh.harvard.edu
## Bug reports: analysis-bugs@nmr.mgh.harvard.edu
##

load [file dirname [info script]]/test_ToglManager[info sharedlibextension] Test_toglmanager

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

    toplevel $ww

    frame $fwTop
    togl $twMain -width 300 -height 300 -rgba true -ident $windowID

    frame $fwControls
    button $bwNewWindow -text "New Window" \
	-command CreateWindow

    bind $twMain <Motion> "%W MouseMotionCallback %x %y %b"
    bind $twMain <ButtonPress> "%W MouseDownCallback %x %y %b"
    bind $twMain <ButtonRelease> "%W MouseUpCallback %x %y %b"
    bind $twMain <KeyRelease> "%W KeyUpCallback %x %y %K"
    bind $twMain <KeyPress> "%W KeyDownCallback %x %y %K"
    bind $twMain <Enter> "focus $twMain"
    bind $twMain <FocusOut> "%W ExitCallback"

    pack $twMain -fill both -expand true
    pack $bwNewWindow -side left

    pack $fwTop -fill both -expand true
    pack $fwControls -fill x -expand true

    puts "tcl: created window $windowID"
}


CreateWindow


wm withdraw .
