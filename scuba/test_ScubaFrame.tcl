#!/bin/sh
# the next line restarts using wish \
exec wish "$0" "$@"

##
## test_ScubaFrame.tcl
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2011/05/02 21:02:44 $
##    $Revision: 1.5 $
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


load [file dirname [info script]]/test_ScubaFrame[info sharedlibextension] Test_scubaframe

set gCurWindowID 0
proc GetNewWindowID { } {
    global gCurWindowID
    incr gCurWindowID
    return $gCurWindowID;
}

proc DoGetViewIDFromFrameColRowDlog { iFrameID } {

    set ww       .dlog
    set fwTop    $ww.fwTop
    set lwColRow $fwTop.lwColRow
    set ewCol    $fwTop.ewCol
    set ewRow    $fwTop.ewRow
    set bwOK     $fwTop.bwOK

    toplevel .dlog
    
    frame $fwTop

    label $lwColRow -text "Col/Row"
    entry $ewCol -textvariable gCol -width 4
    entry $ewRow -textvariable gRow -width 4
    button $bwOK -text "OK" \
	-command "DoGetViewIDFromFrameColRow $iFrameID \$gCol \$gRow; destroy $ww"

    grid $lwColRow -column 0 -row 0 
    grid $ewRow    -column 1 -row 0 
    grid $ewCol    -column 2 -row 0 
    grid $bwOK     -column 0 -row 1 -sticky e -columnspan 3

    pack $fwTop
}

proc DoGetViewIDFromFrameColRow { iFrameID iCol iRow } {
    set ID [GetViewIDFromFrameColRow $iFrameID $iCol $iRow]
    puts "ID for $iCol $iRow is $ID"
}

proc DoGetSelectedViewID { iFrameID } {
    set ID [GetSelectedViewID $iFrameID]
    puts "ID for selected view is $ID"
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

    set fwTestButtons $ww.fwTestButtons
    set bwVID      $fwTestButtons.bwVID
    set bwSVID     $fwTestButtons.bwSVID

    toplevel $ww

    frame $fwTop
    togl $twMain -width 300 -height 300 -rgba true -ident $windowID

    frame $fwControls
    button $bwNewWindow -text "New Window" \
	-command CreateWindow
    button $bw11 -text "11" \
	-command "SetFrameViewConfiguration $windowID c1"
    button $bw22 -text "22" \
	-command "SetFrameViewConfiguration $windowID c22"
    button $bw44 -text "44" \
	-command "SetFrameViewConfiguration $windowID c44"
    button $bw13 -text "13" \
	-command "SetFrameViewConfiguration $windowID c13"

    frame $fwTestButtons
    button $bwVID -text "Get ID from Frame Row" \
	-command "DoGetViewIDFromFrameColRowDlog $windowID"
    button $bwSVID -text "Get Selected ID" \
	-command "DoGetSelectedViewID $windowID"


    bind $twMain <Motion> "%W MouseMotionCallback %x %y %b"
    bind $twMain <ButtonPress> "%W MouseDownCallback %x %y %b"
    bind $twMain <ButtonRelease> "%W MouseUpCallback %x %y %b"
    bind $twMain <KeyRelease> "%W KeyUpCallback %x %y %K"
    bind $twMain <KeyPress> "%W KeyDownCallback %x %y %K"
    bind $twMain <Enter> "focus $twMain"

    pack $twMain -fill both -expand true
    pack $bwNewWindow $bw11 $bw22 $bw44 $bw13 $bwVID $bwSVID -side left


    pack $fwTop -fill both -expand true
    pack $fwControls $fwTestButtons -fill x -expand true -side top

    puts "tcl: created window $windowID"

    SetFrameViewConfiguration $windowID c1
}


CreateWindow


wm withdraw .
