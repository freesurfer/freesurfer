#!/bin/sh
# the next line restarts using wish \
exec wish "$0" "$@"

##
## test_ScubaLayers.tcl
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


load [file dirname [info script]]/test_ScubaLayers[info sharedlibextension] Test_scubalayers

set gCurWindowID 0
proc GetNewWindowID { } {
    global gCurWindowID
    incr gCurWindowID
    return $gCurWindowID;
}

proc SetLayersInViews {} {

    global frameID layerID

    set err [catch { set cRows [GetNumberOfRowsInFrame $frameID] } sResult]
    if { 0 != $err } { puts $sResult; exit }
    for { set nRow 0 } { $nRow < $cRows } { incr nRow } {
	set err [catch { 
	    set cCols [GetNumberOfColsAtRowInFrame $frameID $nRow]
	} sResult]
	if { 0 != $err } { puts $sResult; exit }
	for { set nCol 0 } { $nCol < $cCols } { incr nCol } {
	    
	    set err [catch { 
		set viewID [GetViewIDFromFrameColRow $frameID $nCol $nRow] 
	    } sResult]
	    if { 0 != $err } { puts $sResult; exit }
	    puts "Got viewID from frame $frameID at c$nCol, r$nRow = $viewID"
	    
	    set err [catch { AddLayerToView $viewID $layerID 0 } sResult]
	    if { 0 != $err } { puts $sResult; exit }
	    puts "Added layer $layerID to view $viewID at level 0"
	}
    }
}

proc CreateWindow { } {

    global frameID

    set windowID [GetNewWindowID]
    set frameID $windowID

    set ww          .ww$windowID
    set fwTop       $ww.fwTop
    set twMain      $fwTop.twMain

    set bwNewWindow $fwTop.bwNewWindow
    set bw11        $fwTop.bw11
    set bw22        $fwTop.bw22
    set bw13        $fwTop.bw13

    toplevel $ww

    frame $fwTop
    togl $twMain -width 300 -height 300 -rgba true -ident $windowID

    button $bwNewWindow -text "New Window" \
	-command CreateWindow
    button $bw11 -text "11" \
	-command "SetFrameViewConfiguration $windowID c1; SetLayersInViews"
    button $bw22 -text "22" \
	-command "SetFrameViewConfiguration $windowID c22; SetLayersInViews"
    button $bw13 -text "13" \
	-command "SetFrameViewConfiguration $windowID c13; SetLayersInViews"

    bind $twMain <Motion> "%W MouseMotionCallback %x %y %b"
    bind $twMain <ButtonPress> "%W MouseDownCallback %x %y %b"
    bind $twMain <ButtonRelease> "%W MouseUpCallback %x %y %b"
    bind $twMain <KeyRelease> "%W KeyUpCallback %x %y %K"
    bind $twMain <KeyPress> "%W KeyDownCallback %x %y %K"
    bind $twMain <Enter> "focus $twMain"

    grid $twMain      -column 0 -row 0 -columnspan 4 -sticky news
    grid $bwNewWindow -column 0 -row 1
    grid $bw11        -column 1 -row 1
    grid $bw22        -column 2 -row 1
    grid $bw13        -column 3 -row 1

    grid columnconfigure $fwTop 0 -weight 1
    grid rowconfigure $fwTop 0 -weight 1
    grid rowconfigure $fwTop 1 -weight 0

    pack $fwTop -fill both -expand true

    puts "tcl: created window $windowID"
}


CreateWindow

set err [catch { set volumeID [MakeDataCollection Volume] } sResult]
if { 0 != $err } { puts $sResult; exit }
puts "Made volume $volumeID"

set fnMRI /Users/kteich/work/subjects/bert/mri/T1
if { [info exists env(SUBJECTS_DIR)] } {
    set fnMRI $env(SUBJECTS_DIR)/bert/mri/T1
}

set err [catch { SetVolumeCollectionFileName $volumeID $fnMRI } sResult]
if { 0 != $err } { puts $sResult; exit }
puts "Set volume $volumeID filename to $fnMRI"


set err [catch { set layerID [MakeLayer 2DMRI] } sResult]
if { 0 != $err } { puts $sResult; exit }
puts "Made layer $layerID"
 
set err [catch { Set2DMRILayerVolumeCollection $layerID $volumeID } sResult]
if { 0 != $err } { puts $sResult; exit }
puts "Set layer $layerID volume to $volumeID"

set err [catch { SetLayerLabel $layerID "bert T1" } sResult]
if { 0 != $err } { puts $sResult; exit }
puts "Set layer $layerID label to bert T1"


SetFrameViewConfiguration $frameID c1
SetLayersInViews



wm withdraw .
