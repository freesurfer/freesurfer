#!/bin/sh
# the next line restarts using wish \
exec wish "$0" "$@"

package require Tix

load [file dirname [info script]]/scuba[info sharedlibextension] scuba

# Also look for tkUtils.tcl.
foreach sSourceFileName { tkUtils.tcl tkcon.tcl } {
    set lPath [list "." "$env(DEV)/scripts" "$env(MRI_DIR)/lib/tcl"]
    set bFound 0
    foreach sPath $lPath {
       if { $bFound == 0 } {
	    set sFullFileName [ file join $sPath $sSourceFileName ]
	    set nErr [catch { source $sFullFileName } sResult]
	    if { $nErr == 0 } {
		puts "Reading $sFullFileName"
		set bFound 1;
	    }
	}
    }    
    if { $bFound == 0 } {
	puts "Couldn't load $sSourceFileName: Not found in $lPath"
    }
}

# gTool
#   current - current selected tool (nav,)

# gaWidget
#   window - main window
#   tkconWindow - tkcon window
#   menuBar - menu bar
#     col, row - grid position
#   toolBar - tool bar
#     col, row - grid position
#   scubaFrame - scuba frame object
#     col, row - grid position
#   labelArea - label area
#     col, row - grid position

# gView
#   config - 

proc dputs { isMsg } {
    global gbDebugOutput
    if { $gbDebugOutput } {
	puts "scuba: $isMsg"
    }
}


set gNextFrameID 0
proc GetNewFrameID { } {

    global gNextFrameID
    set frameID $gNextFrameID
    incr gNextFrameID
    return $frameID;
}

proc BuildShortcutDirsList {} {
    global glShortcutDirs env
    set glShortcutDirs {}
    if { [info exists env(SUBJECTS_DIR)] } {
	lappend glShortcutDirs $env(SUBJECTS_DIR)
    }
    if { [info exists env(FREESURFER_DATA)] } {
	lappend glShortcutDirs $env(FREESURFER_DATA)
    }
    if { [info exists env(FREESURFER_HOME)] } {
	lappend glShortcutDirs $env(FREESURFER_HOME)
    }
    if { [info exists env(PWD)] } {
	lappend glShortcutDirs $env(PWD)
    }
    if { [info exists env(FSDEV_TEST_DATA)] } {
	lappend glShortcutDirs $env(FSDEV_TEST_DATA)
    }
}


proc AddDirToShortcutDirsList { iDir } {

    global glShortcutDirs
    foreach dir $glShortcutDirs {
	if { $iDir == $dir } { return }
    }
    lappend glShortcutDirs $iDir
}

proc GetDefaultFileLocation { iType } {
    global gsaDefaultLocation 
    global env
    if { [info exists gsaDefaultLocation($iType)] == 0 } {
	switch $iType {
	    LoadVolume {
		if { [info exists env(SUBJECTS_DIR)] } {
		    set gsaDefaultLocation($iType) $env(SUBJECTS_DIR)
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}	       
	    }
	    default { 
		if { [info exists env(SUBJECTS_DIR)] } {
		    set gsaDefaultLocation($iType) $env(SUBJECTS_DIR)
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}	       
	    }
	}
    }
    return $gsaDefaultLocation($iType)
}

proc SetDefaultFileLocation { iType isValue } {
    global gsaDefaultLocation
    if { [string range $isValue 0 0] == "/" } {
	set gsaDefaultLocation($iType) $isValue
    }
}

proc SetSubjectName { isSubject } {
    global gSubject
    global env
    
    # TODO: make sure this subject exists in the subject directory
    
    set gSubject(name) $isSubject
    set gSubject(homeDir) [file join $env(SUBJECTS_DIR) $isSubject]
    set gSubject(subjectsDir) $env(SUBJECTS_DIR)
}

proc FindFile { ifn } {
    global gSubject

    set fn $ifn

    # If this is not a full path...
    if { [file pathtype $fn] != "absolute" } {

	# If it's partial, and if we have a subject name...
	if { [info exists gSubject(homeDir)] } {
	    
	    # Check a couple of the common places to find files.
	    lappend lfn [file join $gSubject(homeDir) mri $ifn]
	    lappend lfn [file join $gSubject(homeDir) label $ifn]
	    lappend lfn [file join $gSubject(homeDir) surf $ifn]
	    foreach fnTest $lfn {
		if { [file exists $fnTest] } {
		    set fn $fnTest
		    break
		}
	    }
	}
    }

    
    # Make sure this file exists
    if { [file exists $fn] == 0 } {
	tkuErrorDlog "Couldn't find file '$fn'"
	return ""
    }

    # Make sure it is readable
    if { [file readable $fn] == 0 } {
	tkuErrorDlog "File $fn isn't readable"
	return ""
    }

    return $fn
}

proc ExtractLabelFromFileName { ifnData } {
    global gbDebugOutput

    set sSeparator [string range [file join " " " "] 1 1]

    set sSubject ""
    set sData ""
    set sLabel ""

    # Look for 'subjects' to see if we have a subjects path.
    set bFoundSubjectName 0
    if { [string first subjects $ifnData] != -1 } {

	# First look for subjects. If found
	set nBegin [string first subjects $ifnData]
	set nEnd 0
	if { $nBegin != -1 } {
	    incr nBegin 9   ; #skip past subjects/
	    #look for the next separator
	    set nEnd [string first $sSeparator $ifnData $nBegin]
	    if { $nEnd == -1 } { ; # if / not found, just use the whole rest
		set nEnd [string length $ifnData]
	    } else {
		incr nEnd -1 ; # end is at sep now so go back one
	    }
	    set bFoundSubjectName 1
	} 
    }
    
    # look for a ***/SUBJECTNAME/mri/*** pattern. first look
    # for mri, this-1 will be the end.
    if { ! $bFoundSubjectName } {
	
	set nEnd [string first mri $ifnData]
	if { $nEnd != -1 } {
	    incr nEnd -2 ; # go back to before the separator
	    set nBegin $nEnd ; # go backwards until we hit another sep.
	    while { [string range $ifnData $nBegin $nBegin] != $sSeparator &&
		    $nBegin > 0 } {
		incr nBegin -1
	    }
	    if { $nBegin != 0 } {
		incr nBegin ; # skip seprator
		set bFoundSubjectName 1
	    }
	}
	
    } 
    
    # That's it for subject name.
    if { $bFoundSubjectName } {
	set sSubject [string range $ifnData $nBegin $nEnd]
    } else {
	# still not found, just use nothing.
	set sSubject ""
    }
    
    # Data name is between mri/ and the next slash.
    set bFoundDataName 0
    set nBegin [string first mri $ifnData]
    set nEnd 0
    if { $nBegin != -1 } {
	incr nBegin 4
	set nEnd [string first / $ifnData $nBegin]
	if { $nEnd == -1 } {
	    set nEnd [string length $ifnData]
	    set bFoundDataName 1
	}
    }

    if { $bFoundDataName } {
	set sData [string range $ifnData $nBegin $nEnd]
	set sData [file rootname $sData] ; # Remove any file suffixes.
    } else {
	# Not found, just use file name without suffix and path.
	set sData [file rootname [file tail $ifnData]]
    }
    
    set sLabel [string trimleft [string trimright "$sSubject $sData"]]

    return $sLabel
}

set ksImageDir   "$env(FREESURFER_HOME)/lib/images/"
proc LoadImages {} {

    global ksImageDir
    
    set sFileErrors ""
    foreach sImageName { icon_edit_label icon_edit_volume 
	icon_navigate icon_edit_ctrlpts icon_edit_parc icon_line_tool 
	icon_view_single icon_view_multiple icon_view_31 
	icon_cursor_goto icon_cursor_save 
	icon_main_volume icon_aux_volume icon_linked_cursors 
	icon_arrow_up icon_arrow_down icon_arrow_left icon_arrow_right 
	icon_arrow_cw icon_arrow_ccw 
	icon_arrow_expand_x icon_arrow_expand_y 
	icon_arrow_shrink_x icon_arrow_shrink_y 
	icon_orientation_coronal icon_orientation_horizontal 
	icon_orientation_sagittal 
	icon_zoom_in icon_zoom_out 
	icon_brush_square icon_brush_circle icon_brush_3d 
	icon_surface_main icon_surface_original icon_surface_pial 
	icon_snapshot_save icon_snapshot_load 
	icon_marker_crosshair icon_marker_diamond 
	icon_stopwatch } {

	set fnImage [file join $ksImageDir $sImageName.gif]
	if { [catch {image create photo $sImageName -file $fnImage} \
	      sResult] != 0 } {
	    set sFileErrors "$sFileErrors $fnImage"
	}
    }

    if { $sFileErrors != "" } {
	tkuFormattedErrorDlog "Error Loading Images" \
	    "Couldn't load some images." \
	    "Couldn't find the following images: $sFileErrors"
    }
}


proc MakeMenuBar { ifwTop } {

    global gaMenu
    set fwMenuBar     $ifwTop.fwMenuBar
    set gaMenu(file)  $fwMenuBar.mbwFile

    frame $fwMenuBar -border 2 -relief raised

    tkuMakeMenu -menu $gaMenu(file) -label "File" -items {
	{command "Load Volume..." { DoLoadVolumeDlog } }
	{separator}
	{command "Quit" { Quit } }
    }

    pack $gaMenu(file) -side left

    return $fwMenuBar
}


proc MakeToolBar { ifwTop } {
    global gTool
    global gView

    set fwToolBar     $ifwTop.fwToolBar

    frame $fwToolBar -border 2 -relief raised
    
    tkuMakeToolbar $fwToolBar.fwTool \
	-allowzero false \
	-radio true \
	-variable gTool(current) \
	-command {ToolBarWrapper} \
	-buttons {
	    { -type image -name nav -image icon_navigate } 
	}

    set gTool(current) nav

    tkuMakeToolbar $fwToolBar.fwView \
	-allowzero false \
	-radio true \
	-variable gView(config) \
	-command {ToolBarWrapper} \
	-buttons {
	    { -type image -name c1 -image icon_view_single }
	    { -type image -name c22 -image icon_view_multiple }
	    { -type image -name c13 -image icon_view_31 }
	}

    set gView(config) c1

    pack $fwToolBar.fwTool $fwToolBar.fwView -side left

    return $fwToolBar
}

proc ToolBarWrapper { isName iValue } {
    if { $iValue == 1 } {
	switch $isName {
	    nav {

	    }
	    c1 - c22 - c13 {
		SetFrameViewConfiguration [GetMainFrameID] $isName
	    }
	}
    }
}

proc GetMainFrameID {} {
    global gFrameWidgetToID
    global gaWidget
    if { ![info exists gaWidget(scubaFrame)] } {
	return 0
    }
    return $gFrameWidgetToID($gaWidget(scubaFrame))
}

proc MakeScubaFrame { ifwTop } {

    global gFrameWidgetToID

    set fwScuba $ifwTop.fwScuba
    
    set frameID [GetNewFrameID]
    togl $fwScuba -width 300 -height 300 -rgba true -ident $frameID

    bind $fwScuba <Motion> \
	"%W MouseMotionCallback %x %y %b; ScubaMouseMotionCallback %x %y %b"
    bind $fwScuba <ButtonPress> "%W MouseDownCallback %x %y %b"
    bind $fwScuba <ButtonRelease> "%W MouseUpCallback %x %y %b"
    bind $fwScuba <KeyRelease> "%W KeyUpCallback %x %y %K"
    bind $fwScuba <KeyPress> "%W KeyDownCallback %x %y %K"
    bind $fwScuba <Enter> "focus $fwScuba"

    set gFrameWidgetToID($fwScuba) $frameID

    return $fwScuba
}

proc ScubaMouseMotionCallback { inX inY iButton } {

    set err [catch { set viewID [GetViewIDAtFrameLocation [GetMainFrameID] $inX $inY] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    set err [catch { set labelValues [GetLabelValuesSet $viewID cursor] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    UpdateLabelArea $labelValues
}

proc MakeLabelArea { ifwTop } {
    global gaWidget

    set fwLabelArea     $ifwTop.fwLabelArea

    frame $fwLabelArea -border 2 -relief raised

    set gaWidget(labelArea,labelValueWidgets) {}
    set gaWidget(labelArea,numberOfLabels) 0

    return $fwLabelArea
}

proc ShowHideLabelArea { ibShow } {
    global gaWidget

    if { $ibShow } {
	grid $gaWidget(labelArea) \
	    -column $gaWidget(labelArea,column) -row $gaWidget(labelArea,row)
    } else {
	grid remove $gaWidget(labelArea)
    }
}

proc UpdateLabelArea { ilLabelValues } {
    global glLabelValues
    set glLabelValues $ilLabelValues

    DrawLabelArea
}

proc DrawLabelArea {} {
    global gaWidget
    global glLabelValues

    foreach fw $gaWidget(labelArea,labelValueWidgets) {
#	grid forget $fw
    }

    set nLabel 0
    foreach lLabelValue $glLabelValues {
	
	set label [lindex $lLabelValue 0]
	set value [lindex $lLabelValue 1]

	set ewLabel $gaWidget(labelArea).ewLabel$nLabel
	set ewValue $gaWidget(labelArea).ewValue$nLabel
	
	if { $nLabel >= $gaWidget(labelArea,numberOfLabels) } {

	    tkuMakeNormalLabel $ewLabel -label $label
	    tkuMakeNormalLabel $ewValue -label $value
	    
	    grid $ewLabel -column 0 -row $nLabel -sticky ew
	    grid $ewValue -column 1 -row $nLabel -sticky ew
	    
	    lappend gaWidget(labelArea,labelValueWidgets) $ewLabel
	    lappend gaWidget(labelArea,labelValueWidgets) $ewValue

	} else {
	    
	    $ewLabel.lw config -text $label
	    $ewValue.lw config -text $value
	}

	incr nLabel
    }

    grid columnconfigure $gaWidget(labelArea) 0 -weight 1
    grid columnconfigure $gaWidget(labelArea) 1 -weight 1

    if { $nLabel > $gaWidget(labelArea,numberOfLabels) } {
	set gaWidget(labelArea,numberOfLabels) $nLabel
    }
}

proc MakeVolumeCollection { ifnVolume } {

    set err [catch { set colID [MakeDataCollection Volume] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    set err [catch { SetVolumeCollectionFileName $colID $ifnVolume } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    # Get a good name for the collection.
    set sLabel [ExtractLabelFromFileName $ifnVolume]
    
    set err [catch { SetCollectionLabel $colID $sLabel } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    return $colID
}

proc Make2DMRILayer { isLabel } {

    set err [catch { set layerID [MakeLayer 2DMRI] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    set err [catch { SetLayerLabel $layerID $isLabel } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

     return $layerID
}

proc AddLayerToAllViewsInFrame { iFrameID iLayerID } {

    set err [catch { set cRows [GetNumberOfRowsInFrame $iFrameID] } sResult]
    if { 0 != $err } { tkuErrorDlog "$sResult"; return }
    
    for { set nRow 0 } { $nRow < $cRows } { incr nRow } {
	
	set err [catch { 
	    set cCols [GetNumberOfColsAtRowInFrame $iFrameID $nRow]
	} sResult]
	if { 0 != $err } { tkuErrorDlog $sResult; return }
	
	for { set nCol 0 } { $nCol < $cCols } { incr nCol } {
	    
	    set err [catch { 
		set viewID [GetViewIDFromFrameColRow $iFrameID $nCol $nRow] 
	    } sResult]
	    if { 0 != $err } { tkuErrorDlog $sResult; return }

	    set err [catch { AddLayerToView $viewID $iLayerID 0 } sResult]
	    if { 0 != $err } { tkuErrorDlog $sResult; return }
       }
   }
}

proc LoadVolume { ifnVolume ibCreateLayer iFrameIDToAdd } {

    set fnVolume [FindFile $ifnVolume]

    set err [catch { set colID [MakeVolumeCollection $fnVolume] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    if { $ibCreateLayer } {

	set err [catch { set layerID [Make2DMRILayer volume] } sResult]
	if { 0 != $err } { tkuErrorDlog $sResult; return }

	set err [catch { SetVolumeCollection $layerID $colID } sResult]
	if { 0 != $err } { tkuErrorDlog $sResult; return }
	
	if { $iFrameIDToAdd != -1 } {
	    AddLayerToAllViewsInFrame $iFrameIDToAdd $layerID
	}
    }

    # Add this directory to the shortcut dirs if it isn't there already.
    AddDirToShortcutDirsList [file dirname $ifnVolume]
}

proc DoLoadVolumeDlog {} {
    global glShortcutDirs

    tkuDoFileDlog -title "Load Volume" \
	-prompt1 "Load Volume: " \
	-defaultdir1 [GetDefaultFileLocation LoadVolume] \
	-shortcuts $glShortcutDirs \
	-okCmd { 
	    LoadVolume %s1 1 $gFrameWidgetToID($gaWidget(scubaFrame))
	}

    # Future options.
#       -type2 checkbox \
#	-prompt2 "Automatically create new layer" \
#	-defaultvalue2 1 \
#	-type3 checkbox \
#	-prompt3 "Automatically add new layer to all views" \
#	-defaultvalue3 1 \
}


# Look at our command line args. For some we will want to process and
# exit without bringing up all our windows. For some, we need to bring
# up our windows first. So cache those in lCommands and we'll execute
# them later.
set lCommands {}
set nArg 0
while { $nArg < $argc } {
    set sArg [lindex $argv $nArg]
    set sOption [string range $sArg [expr [string last "-" $sArg]+1] end]
    switch $sOption {
	m - volume {
	    incr nArg
	    set fnVolume [lindex $argv $nArg]
	    lappend lCommands "LoadVolume $fnVolume 1 [GetMainFrameID]"
	}
	j - subject {
	    incr nArg
	    set sSubject [lindex $argv $nArg]
	    lappend lCommands "SetSubjectName $sSubject"
	}
	help - default {
	    if {$sOption != "help"} {puts "Option $sOption not recognized."}
	    puts ""
	    puts "scuba"
	    exit
	}
    }
    incr nArg
}

BuildShortcutDirsList
LoadImages


# Make the tkcon window. This must be done at this scope because the
# tkcon.tcl script needs access to some global vars.
set av $argv
set argv "" 
toplevel .dummy
::tkcon::Init
tkcon attach main
wm geometry .tkcon -10-10
destroy .dummy
set argv $av

set gaWidget(tkconWindow) .tkcon


# Make the main window.
set gaWidget(window) .main
toplevel $gaWidget(window)

# Make the areas in the window.
set gaWidget(menuBar) [MakeMenuBar $gaWidget(window)]
set gaWidget(toolBar) [MakeToolBar $gaWidget(window)]
set gaWidget(scubaFrame) [MakeScubaFrame $gaWidget(window)]
set gaWidget(labelArea) [MakeLabelArea $gaWidget(window)]

# Set the grid coords of our areas and the grid them in.
set gaWidget(menuBar,column)    0; set gaWidget(menuBar,row)    0
set gaWidget(toolBar,column)    0; set gaWidget(toolBar,row)    1
set gaWidget(scubaFrame,column) 0; set gaWidget(scubaFrame,row) 2
set gaWidget(labelArea,column)  0; set gaWidget(labelArea,row)  3

grid $gaWidget(menuBar) \
    -column $gaWidget(menuBar,column) -row $gaWidget(menuBar,row) \
    -sticky ew
grid $gaWidget(toolBar) \
    -column $gaWidget(toolBar,column) -row $gaWidget(toolBar,row) \
    -sticky ew
grid $gaWidget(scubaFrame) \
    -column $gaWidget(scubaFrame,column) -row $gaWidget(scubaFrame,row) \
    -sticky news
grid $gaWidget(labelArea) \
    -column $gaWidget(labelArea,column) -row $gaWidget(labelArea,row) \
    -sticky ew

grid columnconfigure $gaWidget(window) 0 -weight 1
grid rowconfigure $gaWidget(window) 0 -weight 0
grid rowconfigure $gaWidget(window) 1 -weight 0
grid rowconfigure $gaWidget(window) 2 -weight 1
grid rowconfigure $gaWidget(window) 3 -weight 0

wm withdraw .

# Now execute all the commands we cached before.
foreach command $lCommands {
    eval $command
}

tkuFinish


proc test_ExtractLabelFromFileName {} {

    foreach { fn sLabel } { 
	/path/to/freesurfer/subjects/SUBJECT/mri/DATA
	"SUBJECT DATA"
	/path/to/some/data/SUBJECT/mri/DATA.mgh
	"SUBJECT DATA"
	/path/to/some/data/DATA.mgh
	"DATA"
    } {
	set sFoundLabel [ExtractLabelFromFileName $fn]
	if { $sFoundLabel != $sLabel } {
	    puts "$fn FAILED, was $sFoundLabel"
	}
    }

}

test_ExtractLabelFromFileName
