
package require Tix


# This will probably be phased out, but for now. Check to see if the
# PrintAllCommands command works. If so, we are running from a binary
# and don't need to load the scuba shared lib. If not, load it.
set err [catch { PrintAllCommands } sResult]
if { $err } {
    load [file dirname [info script]]/libscuba[info sharedlibextension] scuba
}

# gTool
#   current - current selected tool (nav,)

# gaWidget
#   window - main window
#   tkcon - tkcon frame
#   scubaFrame,ID - frame widget for frame ID
#   menuBar - menu bar
#     col, row - grid position
#   toolBar - tool bar
#     col, row - grid position
#   scubaFrame - scuba frame object
#     col, row - grid position
#   labelArea - label area
#     col, row - grid position
#   subjectsLoader
#     subjectsMenu - menu of subjects, value is index gaSubject(nameList)
#     volumeMenu - list of currect subj's mri vols, value is full filename
#     surfaceMenu - list of current subj's surfs, value is full filename
#   layerProperties
#     menu - layer menu, value is layer ID
#   viewProperties
#     menu - view menu, value is view ID
#     drawLevelMenuN - menu of layers at draw level N, value is layer ID
#     transformMenu - menu of transforms, value is transform ID
#   transformProperties
#     menu - transform menu, value is transform ID

# gSubject
#   nameList - directory listing of SUBJECTS_DIR

# gaMenu

# gaSubject

# gaFrame
#   n - id of frame
#     viewConfig
#     toolID

# gaView
#   current
#     id
#     menuIndex
#     id
#     col
#     row
#     linked
#     drawN
#     transformID
#     inPlane

# gaLayer
#   current - currently displayed layer in props panel
#     id
#     type
#     label
#     opacity
#     colorMapMethod  - 2DMRI only
#     sampleMethod - 2DMRI only
#     brightness - 2DMRI only
#     contrast - 2DMRI only
#   idList - list of IDs in layer props listbox

# gaTool
#   n - id of tool
#     mode

set gbDebugOutput false
proc dputs { isMsg } {
    global gbDebugOutput
    if { $gbDebugOutput } {
	puts "scuba: $isMsg"
    }
}


set gNextFrameID 0
proc GetNewFrameID { } {
    dputs "GetNewFrameID   "

    global gNextFrameID
    set frameID $gNextFrameID
    incr gNextFrameID
    return $frameID;
}

proc BuildShortcutDirsList {} {
    dputs "BuildShortcutDirsList  "

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
    dputs "AddDirToShortcutDirsList  $iDir  "


    global glShortcutDirs
    foreach dir $glShortcutDirs {
	if { $iDir == $dir } { return }
    }
    lappend glShortcutDirs $iDir
}

proc GetDefaultFileLocation { iType } {
    dputs "GetDefaultFileLocation  $iType  "

    global gsaDefaultLocation 
    global env
    global gSubject
    if { [info exists gsaDefaultLocation($iType)] == 0 } {
	switch $iType {
	    LoadVolume {
		if { [info exists env(SUBJECTS_DIR)] } {
		    if { [info exists gSubject(homeDir)] } {
			set gsaDefaultLocation($iType) $gSubject(homeDir)/mri
		    } else {
			set gsaDefaultLocation($iType) $env(SUBJECTS_DIR)
		    }
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}	       
	    }
	    SaveLabel - LoadLabel {
		if { [info exists env(SUBJECTS_DIR)] } {
		    if { [info exists gSubject(homeDir)] } {
			set gsaDefaultLocation($iType) $gSubject(homeDir)/label
		    } else {
			set gsaDefaultLocation($iType) $env(SUBJECTS_DIR)
		    }
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}	       
	    }
	    LUT {
		if { [info exists env(FREESURFER_HOME)] } {
		    set gsaDefaultLocation($iType) $env(FREESURFER_HOME)
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
    dputs "SetDefaultFileLocation  $iType $isValue  "

    global gsaDefaultLocation
    if { [string range $isValue 0 0] == "/" } {
	set gsaDefaultLocation($iType) $isValue
    }
}

proc SetSubjectName { isSubject } {
    dputs "SetSubjectName  $isSubject  "

    global gSubject
    global gaWidget
    global env
    
    # Make sure this subject exists in the subject directory
    if { ![info exists env(SUBJECTS_DIR)] } { 
	tkuErrorDlog "SUBJECTS_DIR environment variable not set."
	return
    }
    if { ![file isdirectory $env(SUBJECTS_DIR)/$isSubject] } { 
	tkuErrorDlog "Subject $isSubject doesn't exist."
	return
    }

    # Set some info.
    set gSubject(name) $isSubject
    set gSubject(homeDir) [file join $env(SUBJECTS_DIR) $isSubject]
    set gSubject(subjectsDir) $env(SUBJECTS_DIR)

    # Select it in the subjects loader.
    SelectSubjectInSubjectsLoader $isSubject
	
}

proc FindFile { ifn } {
    dputs "FindFile  $ifn  "

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
    dputs "ExtractLabelFromFileName  $ifnData  "

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
    
    # Volume data name is between mri/ and the next slash.
    set bFoundDataName 0
    set nBegin [string first mri $ifnData]
    set nEnd 0
    if { $nBegin != -1 } {
	incr nBegin 4
	set nEnd [string first / $ifnData $nBegin]
	if { $nEnd == -1 } {
	    set nEnd [string length $ifnData]
	    set sData [string range $ifnData $nBegin $nEnd]
	    set sData [file rootname $sData] ; # Remove any file suffixes.
	    set bFoundDataName 1
	}
    }

    # Surface data name is after surf/.
    set nBegin [string first surf$sSeparator $ifnData]
    if { $nBegin != -1 } {
	incr nBegin 5
	set sData [string range $ifnData $nBegin end]
	set bFoundDataName 1
    }

    if { ! $bFoundDataName } {
	# Not found, just use file name without suffix and path.
	set sData [file rootname [file tail $ifnData]]
    }
    
    set sLabel [string trimleft [string trimright "$sSubject $sData"]]

    return $sLabel
}

set ksImageDir   "$env(FREESURFER_HOME)/lib/images/"
proc LoadImages {} {
    dputs "LoadImages  "


    global ksImageDir
    
    set sFileErrors ""
    foreach sImageName { icon_edit_label icon_edit_volume icon_draw_line
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
    dputs "MakeMenuBar  $ifwTop  "

    global gaMenu
    global gaView
    global gaWidget
    
    set fwMenuBar     $ifwTop.fwMenuBar
    set gaMenu(file)  $fwMenuBar.mbwFile
    set gaMenu(edit)  $fwMenuBar.mbwEdit
    set gaMenu(view)  $fwMenuBar.mbwView

    frame $fwMenuBar -border 2 -relief raised

    tkuMakeMenu -menu $gaMenu(file) -label "File" -items {
	{command "Load Volume..." { DoLoadVolumeDlog } }
	{command "Load Label..." { DoLoadLabelDlog } }
	{separator}
	{command "Save Label..." { DoSaveLabelDlog } }
	{command "Export Labels as Segmenation..." { DoExportLabelsDlog } }
	{command "Quit:Alt Q" { Quit } }
    }

    pack $gaMenu(file) -side left

    tkuMakeMenu -menu $gaMenu(edit) -label "Edit" -items {
	{command "Nothing to undo" 
	    { UndoOrRedo; RedrawFrame [GetMainFrameID]; UpdateUndoMenuItem } }
	{separator}
	{command "Preferences..." { DoPrefsDlog } }
    }

    pack $gaMenu(edit) -side left

    tkuMakeMenu -menu $gaMenu(view) -label "View" -items {
	{check "Flip Left/Right" { SetViewFlipLeftRightYZ $gaView(current,id) $gaView(flipLeftRight) } gaView(flipLeftRight) }
	{check "Coordinate Overlay" { SetPreferencesValue DrawCoordinateOverlay $gaView(coordOverlay); RedrawFrame [GetMainFrameID] } gaView(coordOverlay) }
	{check "Center Crosshair Overlay" { SetPreferencesValue DrawCenterCrosshairOverlay $gaView(crosshairOverlay); RedrawFrame [GetMainFrameID] } gaView(crosshairOverlay) }
	{check "Show Console:Alt N" { ShowHideConsole $gaView(tkcon,visible) } gaView(tkcon,visible) }
    }

    set gaView(flipLeftRight) [GetPreferencesValue ViewFlipLeftRight]
    set gaView(tkcon,visible) [GetPreferencesValue ShowConsole]
    set gaView(coordOverlay)  [GetPreferencesValue DrawCoordinateOverlay]
    set gaView(crosshairOverlay)  [GetPreferencesValue DrawCenterCrosshairOverlay]

    pack $gaMenu(view) -side left

    return $fwMenuBar
}

proc UpdateUndoMenuItem {} {
    global gaMenu

    # Try to get an Undo title. This may fail if there has been no
    # undoable action yet.
    set sLabel "Nothing to undo"
    catch { set sLabel [GetUndoTitle] }

    tkuSetMenuItemName $gaMenu(edit) 1 $sLabel
}

proc MakeToolBar { ifwTop } {
    dputs "MakeToolBar  $ifwTop  "

    global gaTool
    global gaFrame

    set fwToolBar     $ifwTop.fwToolBar

    frame $fwToolBar -border 2 -relief raised

    tkuMakeToolbar $fwToolBar.fwTool \
	-allowzero false \
	-radio true \
	-variable gaTool($gaFrame([GetMainFrameID],toolID),mode) \
	-command {ToolBarWrapper} \
	-buttons {
	    { -type image -name navigation -image icon_navigate } 
	    { -type image -name voxelEditing -image icon_edit_volume } 
	    { -type image -name roiEditing -image icon_edit_label } 
	    { -type image -name straightLine -image icon_line_tool } 
	    { -type image -name edgeLine -image icon_draw_line } 
	}

    set gaTool($gaFrame([GetMainFrameID],toolID),mode) navigation

    tkuMakeToolbar $fwToolBar.fwView \
	-allowzero false \
	-radio true \
	-variable gaFrame([GetMainFrameID],viewConfig) \
	-command {ToolBarWrapper} \
	-buttons {
	    { -type image -name c1 -image icon_view_single }
	    { -type image -name c22 -image icon_view_multiple }
	    { -type image -name c13 -image icon_view_31 }
	}

    set gaFrame([GetMainFrameID],viewConfig) c1

    tkuMakeToolbar $fwToolBar.fwInPlane \
	-allowzero false \
	-radio true \
	-variable gaView(current,inPlane) \
	-command {ToolBarWrapper} \
	-buttons {
	    { -type image -name x -image icon_orientation_sagittal }
	    { -type image -name y -image icon_orientation_coronal }
	    { -type image -name z -image icon_orientation_horizontal }
	}

    set gaView(current,inPlane) x

    pack $fwToolBar.fwTool $fwToolBar.fwView $fwToolBar.fwInPlane \
	-side left

    return $fwToolBar
}

proc ToolBarWrapper { isName iValue } {
    dputs "ToolBarWrapper  $isName $iValue  "

    global gaLayer
    global gaFrame
    global gaROI
    global gaTool
    global gaView

    if { $iValue == 1 } {
	switch $isName {
	    navigation - voxelEditing - roiEditing - straightLine - edgeLine {
		SetToolMode $gaFrame([GetMainFrameID],toolID) $isName
		SelectToolInToolProperties $isName
	    }
	    c1 - c22 - c13 {
		SetFrameViewConfiguration [GetMainFrameID] $isName
		UpdateViewList
	    }
	    x - y - z {
		SetViewInPlane [GetSelectedViewID [GetMainFrameID]] $isName
		set gaView(current,inPlaneInc) \
   	          [GetViewInPlaneMovementIncrement $gaView(current,id) $isName]
		RedrawFrame [GetMainFrameID]
	    }
	    grayscale - heatScale - lut {
		Set2DMRILayerColorMapMethod \
		    $gaLayer(current,id) $gaLayer(current,colorMapMethod)
		RedrawFrame [GetMainFrameID]
	    }
	    nearest - trilinear - sinc - magnitude {
		Set2DMRILayerSampleMethod \
		    $gaLayer(current,id) $gaLayer(current,sampleMethod)
		RedrawFrame [GetMainFrameID]
	    }
	    structure - free {
		SetROIType $gaROI(current,id) $gaROI(current,type)
		RedrawFrame [GetMainFrameID]
	    }
	    circle - square {
		SetToolBrushShape $gaFrame([GetMainFrameID],toolID) \
		    $gaTool(current,brushShape)
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
    dputs "MakeScubaFrame  $ifwTop  "

    global gFrameWidgetToID
    global gaFrame
    global gaTool
    global gaWidget

    set fwScuba $ifwTop.fwScuba
    
    set frameID [GetNewFrameID]
    togl $fwScuba -width 512 -height 512 -rgba true -ident $frameID -time 20

    bind $fwScuba <Motion> \
	"%W MouseMotionCallback %x %y %b; ScubaMouseMotionCallback %x %y %b"
    bind $fwScuba <ButtonPress> \
	"%W MouseDownCallback %x %y %b; ScubaMouseDownCallback %x %y %b"
    bind $fwScuba <ButtonRelease> "%W MouseUpCallback %x %y %b"
    bind $fwScuba <KeyRelease> "%W KeyUpCallback %x %y %K"
    bind $fwScuba <KeyPress> "%W KeyDownCallback %x %y %K"
    bind $fwScuba <Enter> "focus $fwScuba"

    set gaWidget(scubaFrame,$frameID) $fwScuba
    set gFrameWidgetToID($fwScuba) $frameID

    set gaFrame($frameID,toolID) [GetToolIDForFrame $frameID]
    set gaTool($frameID,mode) [GetToolMode $gaFrame($frameID,toolID)]

    return $fwScuba
}

proc MakeScubaFrameBindings { iFrameID } {
    dputs "MakeScubaFrameBindings  $iFrameID  "

    global gaWidget

    set fwScuba $gaWidget(scubaFrame,$iFrameID)

    set sKeyInPlaneX [GetPreferencesValue KeyInPlaneX]
    set sKeyInPlaneY [GetPreferencesValue KeyInPlaneY]
    set sKeyInPlaneZ [GetPreferencesValue KeyInPlaneZ]
    set sKeyCycleView [GetPreferencesValue KeyCycleViewsInFrame]

    bind $fwScuba <Key-$sKeyInPlaneX> {
	set gaView(current,inPlane) x
	SetViewInPlane [GetSelectedViewID [GetMainFrameID]] $gaView(current,inPlane)
    }
    bind $fwScuba <Key-$sKeyInPlaneY> {
	set gaView(current,inPlane) y
	SetViewInPlane [GetSelectedViewID [GetMainFrameID]] $gaView(current,inPlane)
    }
    bind $fwScuba <Key-$sKeyInPlaneZ> { 
	set gaView(current,inPlane) z
	SetViewInPlane [GetSelectedViewID [GetMainFrameID]] $gaView(current,inPlane)
    }
    bind $fwScuba <Key-$sKeyCycleView> {
	CycleCurrentViewInFrame [GetMainFrameID]
	set viewID [GetSelectedViewID [GetMainFrameID]]
	SelectViewInViewProperties $viewID
	RedrawFrame [GetMainFrameID]
    }
}

proc ScubaMouseMotionCallback { inX inY iButton } {
    dputs "ScubaMouseMotionCallback  $inX $inY $iButton  "


    set err [catch { 
	set viewID [GetViewIDAtFrameLocation [GetMainFrameID] $inX $inY] 
    } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    set err [catch { 
	set labelValues [GetLabelValuesSet $viewID cursor] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    UpdateLabelArea $labelValues
}

proc ScubaMouseDownCallback { inX inY iButton } {
    dputs "ScubaMouseDownCallback  $inX $inY $iButton  "

    global gaView

    set viewID [GetSelectedViewID [GetMainFrameID]]
    if { $viewID != $gaView(current,id) } {
	SelectViewInViewProperties $viewID
    }
}

proc GetPreferences {} {
    global gaPrefs

    foreach sKey {
	KeyMoveViewLeft
	KeyMoveViewRight
	KeyMoveViewUp 
	KeyMoveViewDown 
	KeyMoveViewIn 
	KeyMoveViewOut
	KeyZoomViewIn 
	KeyZoomViewOut
    } {
	set gaPrefs($sKey) [GetPreferencesValue $sKey]
    }
}

proc SetPreferences {} {
    global gaPrefs

    foreach sKey {
	KeyMoveViewLeft
	KeyMoveViewRight
	KeyMoveViewUp 
	KeyMoveViewDown 
	KeyMoveViewIn 
	KeyMoveViewOut
	KeyZoomViewIn 
	KeyZoomViewOut
    } {
	SetPreferencesValue $sKey $gaPrefs($sKey)
    }
}

proc Quit {} {
    dputs "Quit  "
    global gaView
    global gaWidget

    # Set our prefs values and save our prefs.
    SetPreferencesValue ViewFlipLeftRight $gaView(flipLeftRight)
    SetPreferencesValue ShowConsole $gaView(tkcon,visible)
    SetPreferencesValue DrawCoordinateOverlay $gaView(coordOverlay)
    SetPreferencesValue DrawCenterCrosshairOverlay $gaView(crosshairOverlay)
    SetPreferences
    SaveGlobalPreferences

    exit
}

# INTERFACE CREATION ==================================================

proc MakeLabelArea { ifwTop } {
    dputs "MakeLabelArea  $ifwTop  "

    global gaWidget

    set fwLabelArea     $ifwTop.fwLabelArea

    frame $fwLabelArea -border 2 -relief raised

    set gaWidget(labelArea,labelValueWidgets) {}
    set gaWidget(labelArea,numberOfLabels) 0

    return $fwLabelArea
}

proc MakeNavigationArea { ifwTop } {
    dputs "MakeNavigationArea  $ifwTop  "

    global gaWidget

    set fwNavArea $ifwTop.fwNavArea
    set fwPan     $fwNavArea.fwPan
    set fwPlane   $fwNavArea.fwPlane
    set fwCenter  $fwNavArea.fwCenter
    set fwInPlane $fwCenter.fwInPlane
    set fwZoom    $fwNavArea.fwZoom
    
    frame $fwNavArea

    frame $fwPan
    button $fwPan.bwUp    -image icon_arrow_up -command { MoveUp }
    button $fwPan.bwDown  -image icon_arrow_down -command { MoveDown }
    button $fwPan.bwLeft  -image icon_arrow_left -command { MoveLeft }
    button $fwPan.bwRight -image icon_arrow_right -command { MoveRight }
    grid $fwPan.bwUp    -column 1 -row 0
    grid $fwPan.bwDown  -column 1 -row 2
    grid $fwPan.bwLeft  -column 0 -row 1
    grid $fwPan.bwRight -column 2 -row 1

    frame $fwPlane
    button $fwPlane.bwIn   -image icon_arrow_up -command { MoveIn }
    button $fwPlane.bwOut  -image icon_arrow_down -command { MoveOut }

    pack $fwPan

    return $fwNavArea
}

proc MakeSubjectsLoaderPanel { ifwTop } {
    dputs "MakeSubjectsLoaderPanel  $ifwTop  "

    global gaWidget
    global gaSubject

    set fwTop  $ifwTop.fwSubjects
    set fwData $fwTop.fwData

    frame $fwTop

    frame $fwData -relief ridge -border 2

    tixOptionMenu $fwData.fwMenu \
	-label "Subject:" \
	-variable gaSubject(current,menuIndex) \
	-command { SubjectsLoaderSubjectMenuCallback }
    set gaWidget(subjectsLoader,subjectsMenu) $fwData.fwMenu

    grid $fwData.fwMenu  -column 0 -row 0 -columnspan 2 -sticky ew

    tixOptionMenu $fwData.volumesMenu \
	-label "Volumes:"
    set gaWidget(subjectsLoader,volumeMenu) $fwData.volumesMenu

    button $fwData.volumesButton \
	-text "Load" \
	-command {LoadVolumeFromSubjectsLoader [$gaWidget(subjectsLoader,volumeMenu) cget -value]}

    grid $fwData.volumesMenu   -column 0 -row 1 -sticky ew
    grid $fwData.volumesButton -column 1 -row 1 -sticky e


    tixOptionMenu $fwData.surfacesMenu \
	-label "Surfaces:"
    set gaWidget(subjectsLoader,surfaceMenu) $fwData.surfacesMenu

    button $fwData.surfacesButton \
	-text "Load" \
	-command {LoadSurfaceFromSubjectsLoader [$gaWidget(subjectsLoader,surfaceMenu) cget -value]}
    
    grid $fwData.surfacesMenu   -column 0 -row 2 -sticky ew
    grid $fwData.surfacesButton -column 1 -row 2 -sticky e


    grid columnconfigure $fwData 0 -weight 1
    grid columnconfigure $fwData 1 -weight 0
    
    pack $fwData -side top -expand yes -fill x

    return $fwTop
}

proc MakePropertiesPanel { ifwTop } {
    dputs "MakePropertiesPanel  $ifwTop  "

    global gaPanel
    global gaWidget

    set fwTop  $ifwTop.fwProps

    frame $fwTop
    set fwButtons  $fwTop.fwButtons
    set fwPanels   $fwTop.fwPanels
    
    frame $fwButtons -relief raised -border 2

    tkuMakeToolbar $fwButtons.tbwPanelsTop \
	-allowzero 1 -radio 1 \
	-variable gaPanel(currentTop) \
	-command PanelBarWrapper \
	-buttons {
	    {-type text -name subjectsLoader -label "Subjects"}
	    {-type text -name viewProperties -label "Views"}
	    {-type text -name layerProperties -label "Layers"}
	    {-type text -name toolProperties -label "Tools"}
	}
    tkuMakeToolbar $fwButtons.tbwPanelsBottom \
	-allowzero 1 -radio 1 \
	-variable gaPanel(currentBottom) \
	-command PanelBarWrapper \
	-buttons {
	    {-type text -name collectionProperties -label "Data"}
	    {-type text -name transformProperties -label "Transforms"}
	    {-type text -name lutProperties -label "Color LUTs"}
	}
    
    pack $fwButtons.tbwPanelsTop $fwButtons.tbwPanelsBottom -side top \
	-fill x -expand yes
    
    frame $fwPanels

    pack [frame $fwPanels.fwPanel]

    set gaWidget(collectionProperties) \
	[MakeDataCollectionsPropertiesPanel $fwPanels.fwPanel]
    set gaWidget(layerProperties) \
	[MakeLayerPropertiesPanel $fwPanels.fwPanel]
    set gaWidget(viewProperties) \
	[MakeViewPropertiesPanel $fwPanels.fwPanel]
    set gaWidget(subjectsLoader) \
	[MakeSubjectsLoaderPanel $fwPanels.fwPanel]
    set gaWidget(transformProperties) \
	[MakeTransformsPanel $fwPanels.fwPanel]
    set gaWidget(lutProperties) \
	[MakeLUTsPanel $fwPanels.fwPanel]
    set gaWidget(toolProperties) \
	[MakeToolsPanel $fwPanels.fwPanel]

    set gaPanel(currentTop) subjectsLoader
    PanelBarWrapper subjectsLoader 1

    pack $fwButtons $fwPanels -side top -fill x -expand yes
    
    return $fwTop
}

proc PanelBarWrapper { isName iValue } {
    dputs "PanelBarWrapper  $isName $iValue  "

    global gaPanel
    global gaWidget

    if { $iValue == 0 } {
	pack forget $gaWidget($isName)
    }
    if { $iValue == 1 } {
	pack $gaWidget($isName)
	switch $isName {
	    subjectsLoader - viewProperties - 
	    layerProperties - toolProperties {
		set gaPanel(currentBottom) ""
	    }
	    collectionProperties - transformProperties - lutProperties {
		set gaPanel(currentTop) ""
	    }
	}
    }
}

proc MakeDataCollectionsPropertiesPanel { ifwTop } {
    dputs "MakeDataCollectionsPropertiesPanel  $ifwTop  "

    global gaWidget
    global gaCollection
    global gaROI
    global glShortcutDirs

    set fwTop        $ifwTop.fwCollectionsProps
    set fwProps      $fwTop.fwProps
    set fwROIs       $fwTop.fwROIs
    set fwCommands   $fwTop.fwCommands

    frame $fwTop

    frame $fwProps -relief ridge -border 2
    set fwMenu         $fwProps.fwMenu
    set fwPropsCommon  $fwProps.fwPropsCommon
    set fwPropsVolume  $fwProps.fwPropsVolume
    set fwPropsSurface $fwProps.fwPropsSurface

    tixOptionMenu $fwMenu \
	-label "Data Collection:" \
	-variable gaCollection(current,menuIndex) \
	-command { CollectionPropertiesMenuCallback }
    set gaWidget(collectionProperties,menu) $fwMenu

    frame $fwPropsCommon
    tkuMakeActiveLabel $fwPropsCommon.ewID \
	-variable gaCollection(current,id) -width 2
    tkuMakeActiveLabel $fwPropsCommon.ewType \
	-variable gaCollection(current,type) -width 10
    tkuMakeEntry $fwPropsCommon.ewLabel \
	-variable gaCollection(current,label) \
	-command {SetCollectionLabel $gaCollection(current,id) $gaCollection(current,label); UpdateCollectionList} \
	-notify 1
    set gaWidget(collectionProperties,labelEntry) $fwPropsCommon.ewLabel

    grid $fwPropsCommon.ewID      -column 0 -row 0               -sticky nw
    grid $fwPropsCommon.ewType    -column 1 -row 0               -sticky new
    grid $fwPropsCommon.ewLabel   -column 0 -row 1 -columnspan 2 -sticky we


    frame $fwPropsVolume
    tkuMakeFileSelector $fwPropsVolume.fwVolume \
	-variable gaCollection(current,fileName) \
	-text "Volume file name:" \
	-shortcutdirs [list $glShortcutDirs] \
	-command {SetVolumeCollectionFileName $gaCollection(current,id) $gaCollection(current,fileName); RedrawFrame [GetMainFrameID]}
    
    grid $fwPropsVolume.fwVolume -column 0 -row 0 -sticky ew
    set gaWidget(collectionProperties,volume) $fwPropsVolume


    frame $fwPropsSurface
    tkuMakeFileSelector $fwPropsSurface.fwSurface \
	-variable gaCollection(current,fileName) \
	-text "Surface file name:" \
	-shortcutdirs [list $glShortcutDirs] \
	-command {SetSurfaceCollectionFileName $gaCollection(current,id) $gaCollection(current,fileName); RedrawFrame [GetMainFrameID]}
    
    grid $fwPropsSurface.fwSurface -column 0 -row 0 -sticky ew
    set gaWidget(collectionProperties,surface) $fwPropsSurface


    frame $fwROIs -relief ridge -border 2
    tixOptionMenu $fwROIs.menu \
	-label "Current ROI:" \
	-variable gaROI(current,menuIndex) \
	-command { ROIPropertiesMenuCallback }
    set gaWidget(roiProperties,menu) $fwROIs.menu

    tkuMakeActiveLabel $fwROIs.ewID \
	-variable gaROI(current,id) -width 2
    tkuMakeEntry $fwROIs.ewLabel \
	-variable gaROI(current,label) \
	-command {SetROILabel $gaROI(current,id) $gaROI(current,label); UpdateROIList} \
	-notify 1
    set gaWidget(roiProperties,labelEntry) $fwROIs.ewLabel

    tkuMakeToolbar $fwROIs.tbwType \
	-allowzero 0 -radio 1 \
	-variable gaROI(current,type) \
	-command ToolBarWrapper \
	-buttons {
	    {-type text -name free -label "Free"}
	    {-type text -name structure -label "Structure"}
	}
    
    tixOptionMenu $fwROIs.mwLUT \
	-label "LUT:" \
	-command "ROIPropertiesLUTMenuCallback"
    set gaWidget(roiProperties,lutMenu) $fwROIs.mwLUT

    tixScrolledListBox $fwROIs.lbStructure \
	-scrollbar auto \
	-browsecmd ROIPropertiesStructureListBoxCallback
   $fwROIs.lbStructure subwidget listbox configure -selectmode single
   set gaWidget(roiProperties,structureListBox) $fwROIs.lbStructure

    tkuMakeActiveLabel $fwROIs.ewStructure \
	-variable gaROI(current,structureLabel)
    
   
    # hack, necessary to init color pickers first time
    set gaROI(current,redColor) 0
    set gaROI(current,greenColor) 0
    set gaROI(current,blueColor) 0
    
    tkuMakeColorPickers $fwROIs.cpFree \
	-pickers {
	    {-label "Free Color:" 
		-redVariable   gaROI(current,redColor) 
		-greenVariable gaROI(current,greenColor)
		-blueVariable  gaROI(current,blueColor)
		-command {SetROIColor $gaROI(current,id) $gaROI(current,redColor) $gaROI(current,greenColor) $gaROI(current,blueColor); RedrawFrame [GetMainFrameID]}}
	}
    set gaWidget(roiProperties,freeColor) $fwROIs.cpFree

    grid $fwROIs.menu        -column 0 -columnspan 2 -row 0 -sticky new
    grid $fwROIs.ewID        -column 0 -row 1   -sticky nw
    grid $fwROIs.ewLabel     -column 1 -row 1   -sticky new
    grid $fwROIs.tbwType     -column 0 -columnspan 2 -row 2   -sticky new
    grid $fwROIs.mwLUT       -column 0 -columnspan 2 -row 3   -sticky new
    grid $fwROIs.lbStructure -column 0 -columnspan 2 -row 4   -sticky new
    grid $fwROIs.ewStructure -column 0 -columnspan 2 -row 5   -sticky new
    grid $fwROIs.cpFree      -column 0 -columnspan 2 -row 6   -sticky new

    grid columnconfigure $fwROIs 0 -weight 0
    grid columnconfigure $fwROIs 1 -weight 1

    frame $fwCommands
    button $fwCommands.bwMakeROI -text "Make New ROI" \
	-command { set roiID [NewCollectionROI $gaCollection(current,id)]; SetROILabel $roiID "New ROI"; UpdateROIList; SelectROIInROIProperties $roiID }
    pack $fwCommands.bwMakeROI -expand yes -fill x


    grid $fwMenu        -column 0 -row 0 -sticky news
    grid $fwPropsCommon -column 0 -row 1 -sticky news

    grid $fwProps    -column 0 -row 1 -sticky news
    grid $fwROIs     -column 0 -row 3 -sticky news
    grid $fwCommands -column 0 -row 4 -sticky news

    return $fwTop
}

proc MakeToolsPanel { ifwTop } {
    dputs "MakeToolsPanel  $ifwTop  "

    global gaWidget
    global gaTool

    set fwTop        $ifwTop.fwToolsProps
    set fwProps      $fwTop.fwProps

    frame $fwTop

    frame $fwProps -relief ridge -border 2

    set fwMenu               $fwProps.fwMenu
    set fwPropsCommon        $fwProps.fwPropsCommon
    set fwPropsBrush         $fwProps.fwPropsBrush
    set fwPropsFill          $fwProps.fwPropsFill

    tixOptionMenu $fwMenu \
	-label "Tools:" \
	-variable gaTool(current,menuIndex) \
	-command { ToolPropertiesMenuCallback }
    set gaWidget(toolProperties,menu) $fwMenu

    FillMenuFromList $fwMenu \
	{ navigation voxelEditing roiEditing straightLine edgeLine }  "" \
	{ "Navigation" "Voxel Editing" "ROI Editing" "Straight Line" "Edge Line"} false


    frame $fwPropsCommon


    tixLabelFrame $fwPropsBrush \
	-label "Brush Options" \
	-labelside acrosstop \
	-options { label.padX 5 }

    set fwPropsBrushSub [$fwPropsBrush subwidget frame]
    
    tkuMakeToolbar $fwPropsBrushSub.tbwBrushShape \
	-allowzero false \
	-radio true \
	-variable gaTool(current,brushShape) \
	-command { ToolBarWrapper } \
	-buttons {
	    {-type text -name circle -label "Circle"}
	    {-type text -name square -label "Square"}
	}

    tkuMakeSliders $fwPropsBrushSub.swRadius -sliders {
	{-label "Radius" -variable gaTool(current,radius) 
	    -min 1 -max 20 -entry true
	    -command {SetToolBrushRadius $gaFrame([GetMainFrameID],toolID) $gaTool(current,radius)} }
    }

    tkuMakeCheckboxes $fwPropsBrushSub.cb3D \
	-font [tkuNormalFont] \
	-checkboxes { 
	    {-type text -label "Work in 3D" 
		-variable gaTool(current,brush3D)
		-command {SetToolBrush3D $gaFrame([GetMainFrameID],toolID) $gaTool(current,brush3D)}}
	}

    grid $fwPropsBrushSub.tbwBrushShape -column 0 -row 0 -sticky ew
    grid $fwPropsBrushSub.swRadius      -column 0 -row 1 -sticky ew
    grid $fwPropsBrushSub.cb3D          -column 0 -row 2 -sticky ew

    set gaWidget(toolProperties,brush) $fwPropsBrush


    tixLabelFrame $fwPropsFill \
	-label "Fill Options" \
	-labelside acrosstop \
	-options { label.padX 5 }

    set fwPropsFillSub [$fwPropsFill subwidget frame]

    tkuMakeCheckboxes $fwPropsFillSub.cbFillOptions \
	-font [tkuNormalFont] \
	-checkboxes { 
	    {-type text -label "Stop at other ROIs" 
		-variable gaTool(current,stopROI)
		-command {SetToolFloodStopAtROIs $gaFrame([GetMainFrameID],toolID) $gaTool(current,stopROI)}}
	    {-type text -label "Stop at lines" 
		-variable gaTool(current,stopLines)
		-command {SetToolFloodStopAtLines $gaFrame([GetMainFrameID],toolID) $gaTool(current,stopLines)}}
	}
    
    tkuMakeSliders $fwPropsFillSub.swFuzziness -sliders {
	{-label "Fill Fuzziness" -variable gaTool(current,fuzziness) 
	    -min 0 -max 50 -entry true
	    -command {SetToolFloodFuzziness $gaFrame([GetMainFrameID],toolID) $gaTool(current,fuzziness)}}
	{-label "Fill Max Distance" -variable gaTool(current,maxDistance) 
	    -min 0 -max 100 -entry true
	    -command {SetToolFloodMaxDistance $gaFrame([GetMainFrameID],toolID) $gaTool(current,maxDistance)}}
    }

    tkuMakeCheckboxes $fwPropsFillSub.cb3D \
	-font [tkuNormalFont] \
	-checkboxes { 
	    {-type text -label "Work in 3D" 
		-variable gaTool(current,flood3D)
		-command {SetToolFlood3D $gaFrame([GetMainFrameID],toolID) $gaTool(current,flood3D)}}
	}

    grid $fwPropsFillSub.cbFillOptions -column 0 -row 0 -sticky ew
    grid $fwPropsFillSub.swFuzziness   -column 0 -row 1 -sticky ew
    grid $fwPropsFillSub.cb3D          -column 0 -row 2 -sticky ew

    set gaWidget(toolProperties,fill) $fwPropsFill

    grid $fwMenu        -column 0 -row 0 -sticky news
    grid $fwPropsCommon -column 0 -row 1 -sticky news

    grid $fwProps    -column 0 -row 0 -sticky news

    return $fwTop
}


proc MakeLayerPropertiesPanel { ifwTop } {
    dputs "MakeLayerPropertiesPanel  $ifwTop  "

    global gaWidget
    global gaLayer
    global glShortcutDirs

    set fwTop        $ifwTop.fwLayerProps
    set fwProps      $fwTop.fwProps

    frame $fwTop

    frame $fwProps -relief ridge -border 2
    set fwMenu        $fwProps.fwMenu
    set fwPropsCommon $fwProps.fwPropsCommon
    set fwProps2DMRI  $fwProps.fwProps2DMRI
    set fwProps2DMRIS $fwProps.fwProps2DMRIS

    tixOptionMenu $fwMenu \
	-label "Layer:" \
	-variable gaLayer(current,menuIndex) \
	-command { LayerPropertiesMenuCallback }
    set gaWidget(layerProperties,menu) $fwMenu


    frame $fwPropsCommon
    tkuMakeActiveLabel $fwPropsCommon.ewID \
	-variable gaLayer(current,id) -width 2
    tkuMakeActiveLabel $fwPropsCommon.ewType \
	-variable gaLayer(current,type) -width 5
    tkuMakeEntry $fwPropsCommon.ewLabel \
	-variable gaLayer(current,label) \
	-command {SetLayerLabel $gaLayer(current,id) $gaLayer(current,label); UpdateLayerList} \
	-notify 1
    set gaWidget(layerProperties,labelEntry) $fwPropsCommon.ewLabel
    tkuMakeSliders $fwPropsCommon.swOpacity -sliders {
	{-label "Opacity" -variable gaLayer(current,opacity) 
	    -min 0 -max 1 -resolution 0.1
	    -command {SetLayerOpacity $gaLayer(current,id) $gaLayer(current,opacity); RedrawFrame [GetMainFrameID]}}
    }

    grid $fwPropsCommon.ewID      -column 0 -row 0               -sticky nw
    grid $fwPropsCommon.ewType    -column 1 -row 0               -sticky new
    grid $fwPropsCommon.ewLabel   -column 0 -row 1 -columnspan 2 -sticky we
    grid $fwPropsCommon.swOpacity -column 0 -row 2 -columnspan 2 -sticky we


    frame $fwProps2DMRI
    tkuMakeToolbar $fwProps2DMRI.tbwColorMapMethod \
	-allowzero 0 -radio 1 \
	-variable gaLayer(current,colorMapMethod) \
	-command ToolBarWrapper \
	-buttons {
	    {-type text -name grayscale -label "Grayscale"}
	    {-type text -name heatScale -label "Heat scale"}
	    {-type text -name lut -label "LUT"}
	}

    tixOptionMenu $fwProps2DMRI.mwLUT \
	-label "LUT:" \
	-command "LayerPropertiesLUTMenuCallback"
    set gaWidget(layerProperties,lutMenu) \
	$fwProps2DMRI.mwLUT

    tkuMakeToolbar $fwProps2DMRI.tbwSampleMethod \
	-allowzero 0 -radio 1 \
	-variable gaLayer(current,sampleMethod) \
	-command ToolBarWrapper \
	-buttons {
	    {-type text -name nearest -label "Nearest"}
	    {-type text -name trilinear -label "Trilinear"}
	    {-type text -name sinc -label "Sinc"}
	    {-type text -name magnitude -label "Mag"}
	}
    tkuMakeCheckboxes $fwProps2DMRI.cbwClearZero \
	-font [tkuNormalFont] \
	-checkboxes { 
	    {-type text -label "Draw 0 values clear" 
		-variable gaLayer(current,clearZero) 
		-command {Set2DMRILayerDrawZeroClear $gaLayer(current,id) $gaLayer(current,clearZero); RedrawFrame [GetMainFrameID]} }
	}
    tkuMakeSliders $fwProps2DMRI.swBC -sliders {
	{-label "Brightness" -variable gaLayer(current,brightness) 
	    -min 1 -max 0 -resolution 0.01 
	    -command {Set2DMRILayerBrightness $gaLayer(current,id) $gaLayer(current,brightness); RedrawFrame [GetMainFrameID]}}
	{-label "Contrast" -variable gaLayer(current,contrast) 
	    -min 0 -max 30 -resolution 1
	    -command {Set2DMRILayerContrast $gaLayer(current,id) $gaLayer(current,contrast); RedrawFrame [GetMainFrameID]}}
    }
    tkuMakeSliders $fwProps2DMRI.swMinMax -sliders {
	{-label "Min" -variable gaLayer(current,minVisibleValue) 
	    -min 0 -max 1 -entry 1
	    -command {Set2DMRILayerMinVisibleValue $gaLayer(current,id) $gaLayer(current,minVisibleValue); RedrawFrame [GetMainFrameID]}}
	{-label "Max" -variable gaLayer(current,maxVisibleValue) 
	    -min 0 -max 1 -entry 1
	    -command {Set2DMRILayerMaxVisibleValue $gaLayer(current,id) $gaLayer(current,maxVisibleValue); RedrawFrame [GetMainFrameID]}}
    }
    set gaWidget(layerProperties,minMaxSliders) $fwProps2DMRI.swMinMax

    tkuMakeSliders $fwProps2DMRI.swROIOpacity -sliders {
	{-label "ROI Opacity" -variable gaLayer(current,roiOpacity) 
	    -min 0 -max 1 -resolution 0.1 
	    -command {Set2DMRILayerROIOpacity $gaLayer(current,id) $gaLayer(current,roiOpacity); RedrawFrame [GetMainFrameID]}}
    }


    grid $fwProps2DMRI.tbwColorMapMethod -column 0 -row 0 -sticky ew
    grid $fwProps2DMRI.mwLUT             -column 0 -row 1 -sticky ew
    grid $fwProps2DMRI.cbwClearZero      -column 0 -row 2 -sticky ew
    grid $fwProps2DMRI.tbwSampleMethod   -column 0 -row 3 -sticky ew
    grid $fwProps2DMRI.swBC              -column 0 -row 4 -sticky ew
    grid $fwProps2DMRI.swMinMax          -column 0 -row 5 -sticky ew
    grid $fwProps2DMRI.swROIOpacity      -column 0 -row 6 -sticky ew
    set gaWidget(layerProperties,2DMRI) $fwProps2DMRI

    # hack, necessary to init color pickers first time
    set gaLayer(current,redLineColor) 0
    set gaLayer(current,greenLineColor) 0
    set gaLayer(current,blueLineColor) 0

    frame $fwProps2DMRIS
    tkuMakeColorPickers $fwProps2DMRIS.cpColors \
	-pickers {
	    {-label "Line Color" -redVariable gaLayer(current,redLineColor) 
		-greenVariable gaLayer(current,greenLineColor)
		-blueVariable gaLayer(current,blueLineColor)
		-command {Set2DMRISLayerLineColor $gaLayer(current,id) $gaLayer(current,redLineColor) $gaLayer(current,greenLineColor) $gaLayer(current,blueLineColor); RedrawFrame [GetMainFrameID]}}
	}
    set gaWidget(layerProperties,lineColorPickers) $fwProps2DMRIS.cpColors

    grid $fwProps2DMRIS.cpColors        -column 0 -row 1 -sticky ew
    set gaWidget(layerProperties,2DMRIS) $fwProps2DMRIS


    grid $fwMenu        -column 0 -row 0 -sticky news
    grid $fwPropsCommon -column 0 -row 1 -sticky news

    grid $fwMenu -column 0 -row 0 -sticky new
    grid $fwProps -column 0 -row 1 -sticky news

    return $fwTop
}

proc MakeViewPropertiesPanel { ifwTop } {
    dputs "MakeViewPropertiesPanel  $ifwTop  "

    global gaWidget
    global gaView

    set fwTop        $ifwTop.fwViewProps
    set fwProps      $fwTop.fwProps

    frame $fwTop

    frame $fwProps -relief ridge -border 2

    tixOptionMenu $fwProps.fwMenu \
	-label "View:" \
	-variable gaView(current,menuIndex) \
	-command { ViewPropertiesMenuCallback }
    set gaWidget(viewProperties,menu) $fwProps.fwMenu

    tkuMakeActiveLabel $fwProps.ewID \
	-label "ID: " \
	-variable gaView(current,id) -width 2
    tkuMakeActiveLabel $fwProps.ewCol \
	-label "Column: " \
	-variable gaView(current,col) -width 2
    tkuMakeActiveLabel $fwProps.ewRow \
	-label "Row: " \
	-variable gaView(current,row) -width 2

    tkuMakeCheckboxes $fwProps.cbwLinked \
	-checkboxes {
	    {-type text -label "Linked" -variable gaView(current,linked)
		-command {SetViewLinkedStatus $gaView(current,id) $gaView(current,linked)} }
	}

    for { set nLevel 0 } { $nLevel < 10 } { incr nLevel } {
	tixOptionMenu $fwProps.mwDraw$nLevel \
	    -label "Draw Level $nLevel:" \
	    -variable gaView(current,draw$nLevel) \
	    -command "ViewPropertiesDrawLevelMenuCallback $nLevel"
	set gaWidget(viewProperties,drawLevelMenu$nLevel) \
	    $fwProps.mwDraw$nLevel
    }

    tixOptionMenu $fwProps.mwTransform \
	-label "Transform:" \
	-variable gaView(current,transformID) \
	-command "ViewPropertiesTransformMenuCallback"
    set gaWidget(viewProperties,transformMenu) \
	$fwProps.mwTransform
    
    button $fwProps.bwCopyLayers -text "Copy Layers to Other Views" \
	-command { CopyViewLayersToAllViewsInFrame [GetMainFrameID] $gaView(current,id) }

    tkuMakeEntry $fwProps.ewInPlaneInc \
	-label "In-plane key increment" \
	-font [tkuNormalFont] \
	-variable gaView(current,inPlaneInc) \
	-command {SetViewInPlaneMovementIncrement $gaView(current,id) $gaView(current,inPlane) $gaView(current,inPlaneInc) } \
	-notify 1
    set gaWidget(viewProperties,inPlaneInc) $fwProps.ewInPlaneInc
    
    grid $fwProps.fwMenu    -column 0 -row 0 -sticky nw
    grid $fwProps.ewID      -column 0 -row 1 -sticky nw
    grid $fwProps.ewCol     -column 1 -row 1 -sticky e
    grid $fwProps.ewRow     -column 2 -row 1 -sticky e
    grid $fwProps.cbwLinked -column 0 -row 2 -sticky ew -columnspan 3
    for { set nLevel 0 } { $nLevel < 10 } { incr nLevel } {
	grid $fwProps.mwDraw$nLevel \
	    -column 0 -row [expr $nLevel + 3] -sticky ew -columnspan 3
    }
    grid $fwProps.mwTransform  -column 0 -row 14 -sticky ew -columnspan 3
    grid $fwProps.bwCopyLayers -column 0 -row 15 -sticky ew -columnspan 3
    grid $fwProps.ewInPlaneInc -column 0 -row 16 -sticky ew -columnspan 3
    
    grid $fwProps -column 0 -row 0 -sticky news

    return $fwTop
}

proc MakeTransformsPanel { ifwTop } {
    dputs "MakeTransformsPanel  $ifwTop  "

    global gaWidget
    global gaTransform

    set fwTop      $ifwTop.fwTransforms
    set fwProps    $fwTop.fwProps
    set fwCommands $fwTop.fwCommands
 
    frame $fwTop

    frame $fwProps -relief ridge -border 2

    tixOptionMenu $fwProps.fwMenu \
	-label "Transform:" \
	-variable gaTransform(current,menuIndex) \
	-command { TransformPropertiesMenuCallback }
    set gaWidget(transformProperties,menu) $fwProps.fwMenu

    grid $fwProps.fwMenu -column 0 -row 0 -columnspan 4 -sticky new

    tkuMakeEntry $fwProps.ewLabel \
	-variable gaTransform(current,label) \
	-notify 1 \
	-command {SetTransformLabel $gaTransform(current,id) $gaTransform(current,label); UpdateTransformList} 
    set gaWidget(transformProperties,labelEntry) $fwProps.ewLabel
    
    grid $fwProps.ewLabel -column 0 -row 1 -columnspan 4 -sticky ew

    for { set nRow 0 } { $nRow < 4 } { incr nRow } {
	for { set nCol 0 } { $nCol < 4 } { incr nCol } {

	    tkuMakeEntry $fwProps.ewValue$nCol-$nRow \
		-width 6 \
		-variable gaTransform(current,value$nCol-$nRow) \
		-command { UpdateCurrentTransformValueList } \
		-notify 1
	    set gaWidget(transformProperties,value$nCol-$nRow) \
		$fwProps.ewValue$nCol-$nRow
	    
	    grid $fwProps.ewValue$nCol-$nRow \
		-column $nCol -row [expr $nRow + 2] -sticky ew
	}
    }

    button $fwProps.bwSetTransform -text "Set Values" \
	-command { SetTransformValues $gaTransform(current,id) $gaTransform(current,valueList); ClearSetTransformValuesButton }
    set gaWidget(transformProperties,setValuesButton) $fwProps.bwSetTransform

    grid $fwProps.bwSetTransform -column 0 -row 6 -columnspan 4 -sticky ew

    frame $fwCommands
    button $fwCommands.bwMakeTransform -text "Make New Transform" \
	-command { set transformID [MakeNewTransform]; SetTransformLabel $transformID "New Transform"; UpdateTransformList; SelectTransformInTransformProperties $transformID }

    pack $fwCommands.bwMakeTransform -expand yes -fill x

    pack $fwProps $fwCommands -side top -expand yes -fill x

    return $fwTop
}

proc MakeLUTsPanel { ifwTop } {
    dputs "MakeLUTsPanel  $ifwTop  "

    global gaWidget
    global gaLUT
    global glShortcutDirs

    set fwTop      $ifwTop.fwLUTs
    set fwProps    $fwTop.fwProps
    set fwCommands $fwTop.fwCommands
 
    frame $fwTop

    frame $fwProps -relief ridge -border 2

    tixOptionMenu $fwProps.fwMenu \
	-label "LUT:" \
	-variable gaLUT(current,menuIndex) \
	-command { LUTPropertiesMenuCallback }
    set gaWidget(lutProperties,menu) $fwProps.fwMenu

    tkuMakeEntry $fwProps.ewLabel \
	-variable gaLUT(current,label) \
	-notify 1 \
	-command {SetColorLUTLabel $gaLUT(current,id) $gaLUT(current,label); UpdateLUTList} 
    set gaWidget(lutProperties,labelEntry) $fwProps.ewLabel

    tkuMakeFileSelector $fwProps.fwLUT \
	-variable gaLUT(current,fileName) \
	-text "File name:" \
	-shortcutdirs [list $glShortcutDirs] \
	-command {SetColorLUTFileName $gaLUT(current,id) $gaLUT(current,fileName); UpdateLUTList; RedrawFrame [GetMainFrameID]}
    
    grid $fwProps.fwMenu  -column 0 -row 0 -sticky ew
    grid $fwProps.ewLabel -column 0 -row 1 -sticky ew
    grid $fwProps.fwLUT   -column 0 -row 2 -sticky ew

    frame $fwCommands
    button $fwCommands.bwMakeLUT -text "Make New LUT" \
	-command { set lutID [MakeNewColorLUT]; SetColorLUTLabel $lutID "New LUT"; UpdateLUTList; SelectLUTInLUTProperties $lutID }

    pack $fwCommands.bwMakeLUT -expand yes -fill x

    pack $fwProps $fwCommands -side top -expand yes -fill x

    return $fwTop
}


# DATA COLLECTION PROPERTIES FUNCTIONS =====================================

proc CollectionPropertiesMenuCallback { iColID } {
    dputs "CollectionPropertiesMenuCallback  $iColID  "

    SelectCollectionInCollectionProperties $iColID
}

proc SelectCollectionInCollectionProperties { iColID } {
    dputs "SelectCollectionInCollectionProperties  $iColID  "

    global gaWidget
    global gaCollection

    # Unpack the type-specific panels.
    grid forget $gaWidget(collectionProperties,volume)
    grid forget $gaWidget(collectionProperties,surface)

    # Get the general collection properties from the collection and
    # load them into the 'current' slots.
    set gaCollection(current,id) $iColID
    set gaCollection(current,type) [GetCollectionType $iColID]
    set gaCollection(current,label) [GetCollectionLabel $iColID]
    tkuRefreshEntryNotify $gaWidget(collectionProperties,labelEntry)

    # Make sure that this is the item selected in the menu. Disale the
    # callback and set the value of the menu to collection ID. Then
    # reenable the callback.
    $gaWidget(collectionProperties,menu) config -disablecallback 1
    $gaWidget(collectionProperties,menu) config -value $iColID
    $gaWidget(collectionProperties,menu) config -disablecallback 0
    
    # Do the type specific stuff.
    switch $gaCollection(current,type) {
	Volume { 
	    # Pack the type panel.
	    grid $gaWidget(collectionProperties,volume) \
		-column 0 -row 2 -sticky news

	    # Get the type specific properties.
	    set gaCollection(current,fileName) \
		[GetVolumeCollectionFileName $iColID]
	}
	Surface {
	    # Pack the type panel.
	    grid $gaWidget(collectionProperties,surface) \
		-column 0 -row 2 -sticky news

	    # Get the type specific properties.
	    set gaCollection(current,fileName)\
		[GetSurfaceCollectionFileName $iColID]
	}
    }

    # Rebuild the ROI list.
    UpdateROIList
}

# This builds the data collection ID list and populates the menu that
# selects the current collection in the collection props panel.  It
# should be called whenever a collection is created or deleted.
proc UpdateCollectionList {} {
    dputs "UpdateCollectionList  "

    global gaWidget
    global gaCollection

    # Get the layer ID list and build the menu.
    set gaCollection(idList) [GetDataCollectionIDList]
    FillMenuFromList $gaWidget(collectionProperties,menu) \
	$gaCollection(idList) "GetCollectionLabel %s" {} false

    # Reselect the current collection.
    if { [info exists gaCollection(current,id)] && 
	 $gaCollection(current,id) >= 0 } {
	SelectCollectionInCollectionProperties $gaCollection(current,id)
    }

    UpdateROIList
}

# ROI PROPERTIES FUNCTIONS ===============================================

proc ROIPropertiesMenuCallback { iROIID } {
    dputs "ROIPropertiesMenuCallback  $iROIID  "

    SelectROIInROIProperties $iROIID
}

proc SelectROIInROIProperties { iROIID } {
    dputs "SelectROIInROIProperties  $iROIID  "

    global gaWidget
    global gaCollection
    global gaROI

    SelectCollectionROI $gaCollection(current,id) $iROIID

    # Get the general ROI properties from the ROI and load them into
    # the 'current' slots.
    set gaROI(current,id) $iROIID
    set gaROI(current,label) [GetROILabel $iROIID]
    tkuRefreshEntryNotify $gaWidget(roiProperties,labelEntry)
    set gaROI(current,type) [GetROIType $iROIID]
    set gaROI(current,lutID) [GetROILUTID $iROIID]
    set gaROI(current,structure) [GetROIStructure $iROIID]
    set lColor [GetROIColor $iROIID]
    set gaROI(current,redColor) [lindex $lColor 0]
    set gaROI(current,greenColor) [lindex $lColor 1]
    set gaROI(current,blueColor) [lindex $lColor 2]
    tkuUpdateColorPickerValues $gaWidget(roiProperties,freeColor)

    # Make sure that this is the item selected in the menu. Disable the
    # callback and set the value of the menu to collection ID. Then
    # reenable the callback.
    $gaWidget(roiProperties,menu) config -disablecallback 1
    $gaWidget(roiProperties,menu) config -value $iROIID
    $gaWidget(roiProperties,menu) config -disablecallback 0

    # Show the right LUT in the listbox.
    SelectLUTInROIProperties $gaROI(current,lutID)

    SelectStructureInROIProperties $gaROI(current,structure)
}

proc ClearROIInROIProperties {} {
    dputs "ClearROIInROIProperties  "

    global gaWidget
    global gaROI

    # Clear the stuff in the current slots.
    set gaROI(current,id) -1
    set gaROI(current,label) ""
    tkuRefreshEntryNotify $gaWidget(roiProperties,labelEntry)
    
    # Clear the listbox.
    $gaWidget(roiProperties,structureListBox) subwidget listbox \
	delete 0 end
}

proc ROIPropertiesLUTMenuCallback { iLUTID } {
    dputs "ROIPropertiesLUTMenuCallback  $iLUTID  "

    SelectLUTInROIProperties $iLUTID
}


proc SelectLUTInROIProperties { iLUTID } {
    dputs "SelectLUTInROIProperties  $iLUTID  "

    global gaROI
    global gaWidget
    global gaCollection

    # Set the ROI data if we can.
    catch {
	set gaROI(current,lutID) $iLUTID
	SetROILUTID $gaROI(current,id) $gaROI(current,lutID)
    }

    # Clear the listbox.
    $gaWidget(roiProperties,structureListBox) subwidget listbox \
	delete 0 end

    # Put the entries in the list box.
    set cEntries [GetColorLUTNumberOfEntries $gaROI(current,lutID)]
    for { set nEntry 0 } { $nEntry < $cEntries } { incr nEntry } {
	catch {
	    set sLabel "$nEntry: [GetColorLUTEntryLabel $gaROI(current,lutID) $nEntry]"
	    $gaWidget(roiProperties,structureListBox) subwidget listbox \
		insert end $sLabel
	}
    }

    # Make sure the right menu item is selected.
    $gaWidget(roiProperties,lutMenu) config -disablecallback 1
    $gaWidget(roiProperties,lutMenu) config -value $iLUTID
    $gaWidget(roiProperties,lutMenu) config -disablecallback 0

    SelectStructureInROIProperties $gaROI(current,structure)
}

proc ROIPropertiesStructureListBoxCallback {} {
    dputs "ROIPropertiesStructureListBoxCallback  "

    global gaWidget

    set nStructure [$gaWidget(roiProperties,structureListBox) \
			subwidget listbox curselection]

    SelectStructureInROIProperties $nStructure
}


proc SelectStructureInROIProperties { inStructure } {
    dputs "SelectStructureInROIProperties  $inStructure  "

    global gaROI
    global gaWidget
    
    # Set value in ROI.
    catch {
	set gaROI(current,structure) $inStructure
	SetROIStructure $gaROI(current,id) $gaROI(current,structure)
	RedrawFrame [GetMainFrameID]

	# Set the label.
	set gaROI(current,structureLabel) "$inStructure: [GetColorLUTEntryLabel $gaROI(current,lutID) $inStructure]"
    }
    
    # Make sure the structure is highlighted and visible in the listbox.
    catch {
	$gaWidget(roiProperties,structureListBox) subwidget listbox \
	    selection clear 0 end
	$gaWidget(roiProperties,structureListBox) subwidget listbox \
	    selection set $gaROI(current,structure)
	$gaWidget(roiProperties,structureListBox) subwidget listbox \
	    see $gaROI(current,structure)
    }

}

# This builds the roi ID list based on the current data collection and
# populates the menu that selects the current roi in the
# collection/roi props panel.  It should be called whenever a
# collection or roi is created or deleted or when a new collection is
# selected.
proc UpdateROIList {} {
    dputs "UpdateROIList  "

    global gaWidget
    global gaCollection
    global gaROI

    if { [info exists gaCollection(current,id)] &&
	 $gaCollection(current,id) >= 0} {

	# Get the roi ID list and build the menu.
	set gaROI(idList) [GetROIIDListForCollection $gaCollection(current,id)]
	FillMenuFromList $gaWidget(roiProperties,menu) \
	    $gaROI(idList) "GetROILabel %s" {} false

	# Reselect the current ROI. If it doesn't exist in this
	# collection, get a new roi ID from the list we got from the
	# collection. If there aren't any, clear the properties.
	if { [info exists gaROI(current,id)] } {
	    if { [lsearch $gaROI(idList) $gaROI(current,id)] == -1 } {
		if { [llength $gaROI(idList)] > 0 } {
		    set gaROI(current,id) [lindex $gaROI(idList) 0]
		} else {
		    ClearROIInROIProperties
		    return
		}
	    }
	    SelectROIInROIProperties $gaROI(current,id)
	}
    }
}

# LAYER PROPERTIES FUNCTIONS ===========================================

proc LayerPropertiesMenuCallback { iLayerID } {
    dputs "LayerPropertiesMenuCallback  $iLayerID  "

    SelectLayerInLayerProperties $iLayerID
}

proc LayerPropertiesLUTMenuCallback { iLUTID } {
    dputs "LayerPropertiesLUTMenuCallback  $iLUTID  "

    global gaLayer
    
    # Set the LUT in this layer and redraw.
    Set2DMRILayerColorLUT $gaLayer(current,id) $iLUTID
    RedrawFrame [GetMainFrameID]
}

proc SelectLayerInLayerProperties { iLayerID } {
    dputs "SelectLayerInLayerProperties  $iLayerID  "

    global gaWidget
    global gaLayer

    # Unpack the type-specific panels.
    grid forget $gaWidget(layerProperties,2DMRI)
    grid forget $gaWidget(layerProperties,2DMRIS)

    # Get the general layer properties from the specific layer and
    # load them into the 'current' slots.
    set gaLayer(current,id) $iLayerID
    set gaLayer(current,type) [GetLayerType $iLayerID]
    set gaLayer(current,label) [GetLayerLabel $iLayerID]
    set gaLayer(current,opacity) [GetLayerOpacity $iLayerID]
    tkuRefreshEntryNotify $gaWidget(layerProperties,labelEntry)

    # Make sure that this is the item selected in the menu. Disale the
    # callback and set the value of the menu to the layer ID. Then
    # reenable the callback.
    $gaWidget(layerProperties,menu) config -disablecallback 1
    $gaWidget(layerProperties,menu) config -value $iLayerID
    $gaWidget(layerProperties,menu) config -disablecallback 0
    
    # Do the type specific stuff.
    switch $gaLayer(current,type) {
	2DMRI { 
	    # Pack the type panel.
	    grid $gaWidget(layerProperties,2DMRI) -column 0 -row 2 -sticky news

	    # Configure the length of the value sliders.
	    set gaLayer(current,minValue) [Get2DMRILayerMinValue $iLayerID]
	    set gaLayer(current,maxValue) [Get2DMRILayerMaxValue $iLayerID]
	    tkuUpdateSlidersRange $gaWidget(layerProperties,minMaxSliders) \
		$gaLayer(current,minValue) $gaLayer(current,maxValue)

	    # Get the type specific properties.
	    set gaLayer(current,colorMapMethod) \
		[Get2DMRILayerColorMapMethod $iLayerID]
	    set gaLayer(current,clearZero) \
		[Get2DMRILayerDrawZeroClear $iLayerID]
	    set gaLayer(current,sampleMethod) \
		[Get2DMRILayerSampleMethod $iLayerID]
	    set gaLayer(current,brightness) [Get2DMRILayerBrightness $iLayerID]
	    set gaLayer(current,contrast) [Get2DMRILayerContrast $iLayerID]
	    set gaLayer(current,lutID) [Get2DMRILayerColorLUT $iLayerID]
	    set gaLayer(current,minVisibleValue) \
		[Get2DMRILayerMinVisibleValue $iLayerID]
	    set gaLayer(current,maxVisibleValue) \
		[Get2DMRILayerMaxVisibleValue $iLayerID]
	    set gaLayer(current,roiOpacity) \
		[Get2DMRILayerROIOpacity $iLayerID]

	    # Set the LUT menu.
	    $gaWidget(layerProperties,lutMenu) config -disablecallback 1
	    $gaWidget(layerProperties,lutMenu) config \
		-value $gaLayer(current,lutID)
	    $gaWidget(layerProperties,lutMenu) config -disablecallback 0    
	}
	2DMRIS {
	    # Pack the type panel.
	    grid $gaWidget(layerProperties,2DMRIS) \
		-column 0 -row 2 -sticky news

	    # Get the type specific properties.
	    set lColor [Get2DMRISLayerLineColor $iLayerID]
	    set gaLayer(current,redLineColor) [lindex $lColor 0]
	    set gaLayer(current,greenLineColor) [lindex $lColor 1]
	    set gaLayer(current,blueLineColor) [lindex $lColor 2]

	    # Configure color selector.
	    tkuUpdateColorPickerValues \
		$gaWidget(layerProperties,lineColorPickers)
	}
    }
}

# This builds the layer ID list and populates the menu that selects
# the current layer in the layer props panel, and the menus in the
# view props panel. It should be called whenever a layer is created or
# deleted, or when a lyer is added to or removed from a view.
proc UpdateLayerList {} {
    dputs "UpdateLayerList  "

    global gaLayer
    global gaWidget
    global gaView

    # We have two jobs here. First we need to populate the menu that
    # selects the current layer in the layer props panel. Then we need
    # to populate all the level-layer menus in the view props
    # panel. First do the layer props.

    # Get the layer ID list and build the menu.
    set gaLayer(idList) [GetLayerIDList]
    FillMenuFromList $gaWidget(layerProperties,menu) $gaLayer(idList) \
	"GetLayerLabel %s" {} false

    # Reselect the current layer.
    if { [info exists gaLayer(current,id)] && 
	 $gaLayer(current,id) >= 0 } {
	SelectLayerInLayerProperties $gaLayer(current,id)
    }


    # Populate the menus in the view props draw level menus.
    for { set nLevel 0 } { $nLevel < 10 } { incr nLevel } {

	FillMenuFromList $gaWidget(viewProperties,drawLevelMenu$nLevel) \
	    $gaLayer(idList) "GetLayerLabel %s" {} true
    }

    # Make sure the right layers are selected in the view draw level
    # menus.
    UpdateCurrentViewProperties
}

# VIEW PROPERTIES FUNCTIONS =============================================

proc ViewPropertiesMenuCallback { iViewID } {
    dputs "ViewPropertiesMenuCallback  $iViewID  "

    SelectViewInViewProperties $iViewID
}

proc ViewPropertiesDrawLevelMenuCallback { iLevel iLayerID } {
    dputs "ViewPropertiesDrawLevelMenuCallback  $iLevel $iLayerID  "

    global gaView
    global gaLayer
    
    # Set the layer in this view and redraw.
    SetLayerInViewAtLevel $gaView(current,id) $iLayerID $iLevel
    RedrawFrame [GetMainFrameID]

    # Get the new inplane inc value if necessary.
    set gaView(current,inPlaneInc) \
	[GetViewInPlaneMovementIncrement $gaView(current,id) $gaView(current,inPlane)]
}

proc ViewPropertiesTransformMenuCallback { iTransformID } {
    dputs "ViewPropertiesTransformMenuCallback  $iTransformID  "

    global gaView
    global gaTransform
    
    # Set the transform in this view and redraw.
    SetViewTransform $gaView(current,id) $iTransformID
    RedrawFrame [GetMainFrameID]
}

proc SelectViewInViewProperties { iViewID } {
    dputs "SelectViewInViewProperties  $iViewID  "

    global gaWidget
    global gaView

    SetSelectedViewID [GetMainFrameID] $iViewID

    # Get the general view properties from the specific view and
    # load them into the 'current' slots.
    set gaView(current,id) $iViewID
    set gaView(current,col) [GetColumnOfViewInFrame [GetMainFrameID] $iViewID]
    set gaView(current,row) [GetRowOfViewInFrame [GetMainFrameID] $iViewID]
    set gaView(current,linked) [GetViewLinkedStatus $iViewID]
    set gaView(current,transformID) [GetViewTransform $iViewID]
    set gaView(current,inPlane) [GetViewInPlane $iViewID]
    set gaView(current,inPlaneInc) \
	[GetViewInPlaneMovementIncrement $iViewID $gaView(current,inPlane)]
    tkuRefreshEntryNotify $gaWidget(viewProperties,inPlaneInc)

    for { set nLevel 0 } { $nLevel < 10 } { incr nLevel } {
	$gaWidget(viewProperties,drawLevelMenu$nLevel) \
	    config -disablecallback 1
	set gaView(current,draw$nLevel) \
	    [GetLayerInViewAtLevel $iViewID $nLevel]
	$gaWidget(viewProperties,drawLevelMenu$nLevel) \
	    config -disablecallback 0
    }
    
    # Make sure that this is the item selected in the menu. Disale the
    # callback and set the value of the menu to the view ID. Then
    # reenable the callback.
    $gaWidget(viewProperties,menu) config -disablecallback 1
    $gaWidget(viewProperties,menu) config -value $iViewID
    $gaWidget(viewProperties,menu) config -disablecallback 0

    # Select the right transform in the transform menu
    $gaWidget(viewProperties,transformMenu) config -disablecallback 1
    $gaWidget(viewProperties,transformMenu) config \
	-value $gaView(current,transformID)
    $gaWidget(viewProperties,transformMenu) config -disablecallback 0    

    UpdateCurrentViewProperties
}

# This gets the layers at each level of the currently selected view
# and makes sure the draw level menus are set properly. Call it
# whenever a layer has been set in the current view.
proc UpdateCurrentViewProperties {} {
    dputs "UpdateCurrentViewProperties  "

    global gaWidget
    global gaView
    global gaLayer

    for { set nLevel 0 } { $nLevel < 10 } { incr nLevel } {

	# Get the current value of this layer.
	set layerID [GetLayerInViewAtLevel $gaView(current,id) $nLevel]

	# Disable callback.
	$gaWidget(viewProperties,drawLevelMenu$nLevel) \
	    config -disablecallback 1

	# Set the layer id in the draw level.
	$gaWidget(viewProperties,drawLevelMenu$nLevel) config -value $layerID

	# Renable the callback.
	$gaWidget(viewProperties,drawLevelMenu$nLevel) \
	    config -disablecallback 0
    }
}


# This builds the view ID list from the current view configuration and
# populates the menu that selects the view in the view props panel. It
# should be called every time the view configuration changes.
proc UpdateViewList {} {
    dputs "UpdateViewList  "

    global gaView
    global gaWidget

    set gaView(idList) {}
    set gaView(labelList) {}
    
    # Build the ID list.
    set err [catch { set cRows [GetNumberOfRowsInFrame [GetMainFrameID]] } sResult]
    if { 0 != $err } { tkuErrorDlog "$sResult"; return }
    for { set nRow 0 } { $nRow < $cRows } { incr nRow } {
	
	set err [catch { 
	    set cCols [GetNumberOfColsAtRowInFrame [GetMainFrameID] $nRow]
	} sResult]
	if { 0 != $err } { tkuErrorDlog $sResult; return }
	
	for { set nCol 0 } { $nCol < $cCols } { incr nCol } {

	    set err [catch { 
		set viewID [GetViewIDFromFrameColRow [GetMainFrameID] $nCol $nRow] 
	    } sResult]
	    if { 0 != $err } { tkuErrorDlog $sResult; return }

	    lappend gaView(idList) $viewID

	    lappend gaView(labelList) "[GetColumnOfViewInFrame [GetMainFrameID] $viewID], [GetRowOfViewInFrame [GetMainFrameID] $viewID]"
	}
    }

    FillMenuFromList $gaWidget(viewProperties,menu) $gaView(idList) \
	"" $gaView(labelList) false
}

# TOOL PROPERTIES FUNCTIONS =============================================

proc ToolPropertiesMenuCallback { iTool } {
    SelectToolInToolProperties $iTool
}

proc SelectToolInToolProperties { iTool } {
    dputs "SelectToolInToolProperties  $iTool"

    global gaWidget
    global gaTool
    global gaFrame

    # Unpack the type-specific panels.
    grid forget $gaWidget(toolProperties,brush)
    grid forget $gaWidget(toolProperties,fill)

    # Get the general layer properties from the specific layer and
    # load them into the 'current' slots.
    set gaTool(current,id) $gaFrame([GetMainFrameID],toolID)
    set gaTool(current,type) $iTool

    # Make sure that this is the item selected in the menu. Disale the
    # callback and set the value of the menu to the layer ID. Then
    # reenable the callback.
    $gaWidget(toolProperties,menu) config -disablecallback 1
    $gaWidget(toolProperties,menu) config -value $iTool
    $gaWidget(toolProperties,menu) config -disablecallback 0
    
    # Set this tool in the toolbar too.
    set gaTool($gaFrame([GetMainFrameID],toolID),mode) $gaTool(current,type)

    # Do the type specific stuff. Pack the relevant tool panels.
    switch $gaTool(current,type) {
	roiEditing - voxelEditing { 
	    # ROI and voxel editing get brush and fill stuff.
	    grid $gaWidget(toolProperties,brush) -column 0 -row 2 -sticky ew

	    set gaTool(current,brushShape) \
		[GetToolBrushShape $gaTool(current,id)]
	    set gaTool(current,radius) \
		[GetToolBrushRadius $gaTool(current,id)]
	    set gaTool(current,brush3D) \
		[GetToolBrush3D $gaTool(current,id)]

	    grid $gaWidget(toolProperties,fill)  -column 0 -row 3 -sticky ew

	    set gaTool(current,stopROI) \
		[GetToolFloodStopAtROIs $gaTool(current,id)]
	    set gaTool(current,stopLines) \
		[GetToolFloodStopAtLines $gaTool(current,id)]
	    set gaTool(current,maxDistance) \
		[GetToolFloodMaxDistance $gaTool(current,id)]
	    set gaTool(current,flood3D) \
		[GetToolFlood3D $gaTool(current,id)]
	}
    }
}

# SUBJECTS LOADER FUNCTIONS =============================================

proc SubjectsLoaderSubjectMenuCallback { inSubject } {
    dputs "SubjectsLoaderSubjectMenuCallback  $inSubject  "

    global gaSubject

    # Get the name at this index in the nameList, then select that
    # subject.
    set gaSubject(current) [lindex $gaSubject(nameList) $inSubject]
    SelectSubjectInSubjectsLoader $gaSubject(current)
}

proc SelectSubjectInSubjectsLoader { isSubject } {
    dputs "SelectSubjectInSubjectsLoader  $isSubject  "

    global gaWidget
    global gaSubject
    global env

    # Make sure we know this subject.
    set nSubject [lsearch $gaSubject(nameList) $isSubject]
    if { $nSubject == -1 } {
	tkuErrorDlog "Subject $isSubject doesn't exist."
	return
    }

    # Make sure that this is the item selected in the menu. Disale the
    # callback and set the value of the menu to the index of this
    # subject name in the subject name list. Then reenable the callback.
    $gaWidget(subjectsLoader,subjectsMenu) config -disablecallback 1
    $gaWidget(subjectsLoader,subjectsMenu) config -value $nSubject
    $gaWidget(subjectsLoader,subjectsMenu) config -disablecallback 0

    # We need to populate the data menus for this subject.  Empty them
    # first.

    set lEntries [$gaWidget(subjectsLoader,volumeMenu) entries]
    foreach entry $lEntries { 
	$gaWidget(subjectsLoader,volumeMenu) delete $entry
    }

    # For volumes, look for all the $sSubject/mri/ subdirs except
    # transforms.and tmp. Make sure they have COR-.info files in them.
    set lContents [dir -full $env(SUBJECTS_DIR)/$isSubject/mri]
    foreach sItem $lContents {
	if { [file isdirectory $env(SUBJECTS_DIR)/$isSubject/mri/$sItem] &&
	  [file exists $env(SUBJECTS_DIR)/$isSubject/mri/$sItem/COR-.info]} {
	    set sVolume [string trim $sItem /]
	    if { "$sVolume" != "transforms" &&
		 "$sVolume" != "tmp" } {
		$gaWidget(subjectsLoader,volumeMenu) add \
		    command "$env(SUBJECTS_DIR)/$isSubject/mri/$sItem" \
		    -label $sVolume
	    }
	}
    }
    
    set lEntries [$gaWidget(subjectsLoader,surfaceMenu) entries]
    foreach entry $lEntries { 
	$gaWidget(subjectsLoader,surfaceMenu) delete $entry
    }
    # For surfaces, look for all the $sSubject/surf/{l,r}h files.
    set lContents [dir -full $env(SUBJECTS_DIR)/$isSubject/surf]
    foreach sItem $lContents {
	if { [string range $sItem 0 1] == "lh" ||
	     [string range $sItem 0 1] == "rh" } {
	    $gaWidget(subjectsLoader,surfaceMenu) add \
		command "$env(SUBJECTS_DIR)/$isSubject/surf/$sItem" \
		-label $sItem
	}
    }
    
}

proc LoadVolumeFromSubjectsLoader { isVolume } {
    dputs "LoadVolumeFromSubjectsLoader  $isVolume  "

    global gaSubject

    LoadVolume "$isVolume" 1 [GetMainFrameID]
}

proc LoadSurfaceFromSubjectsLoader { isSurface } {
    dputs "LoadSurfaceFromSubjectsLoader  $isSurface  "

    global gaSubject

    LoadSurface "$isSurface" 1 [GetMainFrameID]
}

# Builds the subject nameList by looking in SUBJECTS_DIR.
proc UpdateSubjectList {} {
    dputs "UpdateSubjectList  "

    global gaSubject
    global gaWidget
    global env

    set gaSubject(nameList) {}

    # Disable the menu.
    $gaWidget(subjectsLoader,subjectsMenu) config -disablecallback 1

    # Build the ID list. Go through and make sure each is a
    # directory. Trim slashes.
    set lContents [dir -full $env(SUBJECTS_DIR)]
    foreach sItem $lContents {
	if { [file isdirectory $env(SUBJECTS_DIR)/$sItem] } {
	    lappend gaSubject(nameList) [string trim $sItem /]
	}
    }

    # Empty the current subject menu.
    set lEntries [$gaWidget(subjectsLoader,subjectsMenu) entries]
    foreach entry $lEntries { 
	$gaWidget(subjectsLoader,subjectsMenu) delete $entry
    }
    
    # Add the entries from the subject name list to the menu.
    set nSubject 0
    foreach sSubject $gaSubject(nameList) {
	$gaWidget(subjectsLoader,subjectsMenu) add command $nSubject -label $sSubject
	incr nSubject
    }

    # Reenable the menu.
    $gaWidget(subjectsLoader,subjectsMenu) config -disablecallback 0

    # If we don't have a subject select, select the first one.
    if { ![info exists gaSubject(current,id)] } {
	SelectSubjectInSubjectsLoader [lindex $gaSubject(nameList) 0]
    }
}

# TRANSFORM PROPERTIES FUNCTIONS =========================================

proc TransformPropertiesMenuCallback { iTransformID } {
    dputs "TransformPropertiesMenuCallback  $iTransformID  "

    SelectTransformInTransformProperties $iTransformID
}

proc SelectTransformInTransformProperties { iTransformID } {
    dputs "SelectTransformInTransformProperties  $iTransformID  "

    global gaWidget
    global gaTransform

    # Get the tranforms properties from the transofmr layer and
    # load them into the 'current' slots.
    set gaTransform(current,id) $iTransformID
    set gaTransform(current,label) [GetTransformLabel $iTransformID]
    set gaTransform(current,valueList) [GetTransformValues $iTransformID]
    tkuRefreshEntryNotify $gaWidget(transformProperties,labelEntry)

    # Set the invidual values from the value list.
    for { set nRow 0 } { $nRow < 4 } { incr nRow } {
	for { set nCol 0 } { $nCol < 4 } { incr nCol } {
	    set gaTransform(current,value$nCol-$nRow) \
		[lindex $gaTransform(current,valueList) \
		     [expr ($nRow * 4) + $nCol]]
	    tkuRefreshEntryNotify \
		$gaWidget(transformProperties,value$nCol-$nRow)
	}
    }

     
    # Make sure that this is the item selected in the menu. Disale the
    # callback and set the value of the menu to the transform ID. Then
    # reenable the callback.
    $gaWidget(transformProperties,menu) config -disablecallback 1
    $gaWidget(transformProperties,menu) config -value $iTransformID
    $gaWidget(transformProperties,menu) config -disablecallback 0
}


# This builds the transform ID list and populates the menu that selects
# the current transform in the transform props panel, and the menu in the
# view props panel. It should be called whenever a transform is created or
# deleted, or when a lyer is added to or removed from a view.
proc UpdateTransformList {} {
    dputs "UpdateTransformList  "

    global gaTransform
    global gaWidget
    global gaView

    # Get the transform ID list.
    set gaTransform(idList) [GetTransformIDList]

    # First rebuild the transform list in the transform props panel.
    FillMenuFromList $gaWidget(transformProperties,menu) $gaTransform(idList) \
	"GetTransformLabel %s" {} false

    # Reselect the current transformProperties.
    if { [info exists gaTransform(current,id)] && 
	 $gaTransform(current,id) >= 0 } {
	SelectTransformInTransformProperties $gaTransform(current,id)
    }

    # Now rebuild the transform list in the view props panel.
    FillMenuFromList $gaWidget(viewProperties,transformMenu) \
	$gaTransform(idList) "GetTransformLabel %s" {} false
}

proc UpdateCurrentTransformValueList {} {
    dputs "UpdateCurrentTransformValueList  "

    global gaTransform
    global gaWidget

    set gaTransform(current,valueList) {}

    for { set nRow 0 } { $nRow < 4 } { incr nRow } {
	for { set nCol 0 } { $nCol < 4 } { incr nCol } {
	    lappend gaTransform(current,valueList) \
		$gaTransform(current,value$nCol-$nRow)
	}
    }

    # Change the set button to red to remind the user to click that button.
    $gaWidget(transformProperties,setValuesButton) config -fg red
}

proc ClearSetTransformValuesButton {} {
    dputs "ClearSetTransformValuesButton  "

    global gaWidget

    # Change the set button to normal.
    $gaWidget(transformProperties,setValuesButton) config -fg black
}

# COLOR LUT PROPERTIES FUNCTIONS =========================================

proc LUTPropertiesMenuCallback { iLUTID } {
    dputs "LUTPropertiesMenuCallback  $iLUTID  "

    SelectLUTInLUTProperties $iLUTID
}

proc SelectLUTInLUTProperties { iLUTID } {
    dputs "SelectLUTInLUTProperties  $iLUTID  "

    global gaWidget
    global gaLUT

    # Get the lut properties and load them into the 'current' slots.
    set gaLUT(current,id) $iLUTID
    set gaLUT(current,label) [GetColorLUTLabel $iLUTID]
    tkuRefreshEntryNotify $gaWidget(lutProperties,labelEntry)
    set gaLUT(current,fileName) [GetColorLUTFileName $iLUTID]

     
    # Make sure that this is the item selected in the menu. Disale the
    # callback and set the value of the menu to the lut ID. Then
    # reenable the callback.
    $gaWidget(lutProperties,menu) config -disablecallback 1
    $gaWidget(lutProperties,menu) config -value $iLUTID
    $gaWidget(lutProperties,menu) config -disablecallback 0
}


# This builds the lut ID list and populates the menu that selects
# the current lut in the lut props panel, and the menu in the
# layer props panel. It should be called whenever a transform is created or
# deleted.
proc UpdateLUTList {} {
    dputs "UpdateLUTList  "

    global gaLUT
    global gaWidget

    # Get the lut ID list.
    set gaLUT(idList) [GetColorLUTIDList] 

    # First rebuild the lut list in the lut props panel.
    FillMenuFromList $gaWidget(lutProperties,menu) $gaLUT(idList) \
	"GetColorLUTLabel %s" {} false

    # Reselect the current lutProperties.
    if { [info exists gaLUT(current,id)] && $gaLUT(current,id) >= 0 } {
	SelectLUTInLUTProperties $gaLUT(current,id)
    }

    # Now rebuild the lut list in the layer props panel.
    FillMenuFromList $gaWidget(layerProperties,lutMenu) $gaLUT(idList) \
	"GetColorLUTLabel %s" {} false

    # Rebuild the list in the ROI props.
    FillMenuFromList $gaWidget(roiProperties,lutMenu) $gaLUT(idList) \
	"GetColorLUTLabel %s" {} false
}

# LABEL AREA FUNCTIONS ==================================================

proc ShowHideLabelArea { ibShow } {
    dputs "ShowHideLabelArea  $ibShow  "

    global gaWidget

    if { $ibShow } {
	grid $gaWidget(labelArea) \
	    -column $gaWidget(labelArea,column) -row $gaWidget(labelArea,row)
    } else {
	grid remove $gaWidget(labelArea)
    }
}

proc UpdateLabelArea { ilLabelValues } {
    dputs "UpdateLabelArea  $ilLabelValues  "

    global glLabelValues
    set glLabelValues $ilLabelValues

    DrawLabelArea
}

proc DrawLabelArea {} {
    dputs "DrawLabelArea  "

    global gaWidget
    global glLabelValues

    set nLabel 0
    foreach lLabelValue $glLabelValues {
	
	set label [lindex $lLabelValue 0]
	set value [lindex $lLabelValue 1]

	set fw $gaWidget(labelArea).fw$nLabel
	set ewLabel $fw.ewLabel
	set ewValue $fw.ewValue
	
	if { $nLabel >= $gaWidget(labelArea,numberOfLabels) } {

	    frame $fw
	    
	    tkuMakeNormalLabel $ewLabel -label $label -width 20
	    tkuMakeNormalLabel $ewValue -label $value -width 20
	    
	    pack $ewLabel $ewValue -side left -anchor w
	    pack $fw

	} else {
	    
	    $ewLabel.lw config -text $label
	    $ewValue.lw config -text $value
	}

	incr nLabel
    }

    grid columnconfigure $gaWidget(labelArea) 0 -weight 1
    grid columnconfigure $gaWidget(labelArea) 1 -weight 1

    # Save the new number of labels. The way nLabel is incremented
    # here, it's equal to the number of labels we're displaying. If
    # it's less than the last total..
    if { $nLabel > $gaWidget(labelArea,numberOfLabels) } {
	set gaWidget(labelArea,numberOfLabels) $nLabel

    } elseif { $nLabel < $gaWidget(labelArea,numberOfLabels) } {
	
	# Start with the nLabel before this one and go to our total,
	# setting the labels blank. Don't save the new total because
	# that's really the number of labels we have _created_, and we
	# don't want to recreate them.
	for { set nDelLabel $nLabel }  \
	    { $nDelLabel < $gaWidget(labelArea,numberOfLabels) } \
	    { incr nDelLabel } {
		
		set fw $gaWidget(labelArea).fw$nDelLabel
		$fw.ewLabel.lw config -text ""
		$fw.ewValue.lw config -text ""
	    }

    }
}

# PREFS DIALOG =========================================================

proc DoPrefsDlog {} {
    global gaDialog
    global gaPrefs
    
    set wwDialog .prefs
    if { [tkuCreateDialog $wwDialog "Preferences" {-borderwidth 10}] } {

	set fwKeys     $wwDialog.fwKeys
	set fwButtons  $wwDialog.fwButtons

	frame $fwKeys

	foreach {sKey sLabel} {
	    KeyMoveViewLeft   "Move View Left"
	    KeyMoveViewRight  "Move View Right"
	    KeyMoveViewUp     "Move View Up"
	    KeyMoveViewDown   "Move View Down"
	    KeyMoveViewIn     "Move View In"
	    KeyMoveViewOut    "Move View Out"
	    KeyZoomViewIn     "Zoom View In"
	    KeyZoomViewOut    "Zoom View Out"
	} {
   
	    tkuMakeEntry $fwKeys.fw$sKey \
		-label $sLabel -width 8 -font [tkuNormalFont] \
		-labelwidth 20 \
		-variable gaPrefs($sKey) \
		-command "SetPreferencesValue $sKey \$gaPrefs($sKey)" \
		-notify 1

	    pack $fwKeys.fw$sKey \
		-side top -fill x -expand yes
	   
	}
	

	tkuMakeCloseButton $fwButtons $wwDialog
	
	pack $fwKeys $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5
    }
}

# DIALOGS ==============================================================

proc NewDlog { args } {
    global gaDialog

    # set default arguments for all fields
    set aArgs(-title) ""
    set aArgs(-text) ""
    set aArgs(-meter) false
    set aArgs(-buttons) {}

    # Set arg items.
    array set aArgs $args

    set wwDialog .dialog
    if { [tkuCreateDialog $wwDialog $aArgs(-title) {-borderwidth 10}] } {

	set fwText     $wwDialog.fwText
	set fwMeter    $wwDialog.fwMeter
	set fwButtons  $wwDialog.fwButtons

	set gaDialog(current,text) $aArgs(-text)
	tkuMakeActiveLabel $fwText -variable gaDialog(current,text)

	set gaDialog(current,meterText) "Progress: 0%"
	tkuMakeActiveLabel $fwMeter -variable gaDialog(current,meterText)

	frame $fwButtons
	set nButton 0
	foreach btn $aArgs(-buttons) {
	    button $fwButtons.fw$nButton -text $btn \
		-command "DlogCallback \"$btn\""
	    pack $fwButtons.fw$nButton \
		-side left -anchor e 
	    incr nButton
	}

	pack $fwText $fwMeter $fwButtons -side top -fill x -expand yes

	if { 0 } {
	    after idle [format {
		update idletasks
		wm minsize %s %d [winfo reqheight %s]
		wm geometry %s =%dx[winfo reqheight %s]
	    } $wwDialog 400 $wwDialog $wwDialog 400 $wwDialog]
	}

	set gaDialog(current,callBacks) {}
	set gaDialog(current,window) $wwDialog
	set gaDialog(current,percent) 0
    }
}

proc DlogCallback { isButton } {
    global gaDialog

    lappend gaDialog(current,callBacks) $isButton
}

proc UpdateDlog { args } {
    global gaDialog

    # set default arguments for all fields
    set aArgs(-text) $gaDialog(current,text)
    set aArgs(-percent) $gaDialog(current,percent)

    # Set arg items.
    array set aArgs $args

    set gaDialog(current,text) $aArgs(-text)
    set gaDialog(current,percent) $aArgs(-percent)

    set gaDialog(current,meterText) "Progress: $gaDialog(current,percent)"
}

proc CheckDlogForButtons {} {
    global gaDialog

    set lCallbacks $gaDialog(current,callBacks)
    set gaDialog(current,callBacks) {}
    return $lCallbacks
}

proc CloseDlog {} {
    global gaDialog

    tkuCloseDialog $gaDialog(current,window)
}

# VIEW CONFIGURATION ==================================================

proc SetLayerInAllViewsInFrame { iFrameID iLayerID } {
    dputs "SetLayerInAllViewsInFrame  $iFrameID $iLayerID  "


    # For each view...
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

	    # Get the first unused draw level and add the layer to the
	    # view at that level.
	    set err [catch { 
		set level [GetFirstUnusedDrawLevelInView $viewID] } sResult]
	    if { 0 != $err } { tkuErrorDlog $sResult; return }

	    set err [catch {
		SetLayerInViewAtLevel $viewID $iLayerID $level } sResult]
	    if { 0 != $err } { tkuErrorDlog $sResult; return }
       }
    }
    
    UpdateLayerList
}

proc FillMenuFromList { imw ilEntries iLabelFunction ilLabels ibNone  } {

    # Disable callback.
    $imw config -disablecallback 1
    
    # Delete all the entries and add ones for all the IDs in the
    # ID list. Also add a command for 'none' with index of -1.
    set lEntries [$imw entries]
    foreach entry $lEntries { 
	$imw delete $entry
    }

    if { $ibNone } {
	$imw add command -1 -label "None"
    }

    set nEntry 0
    foreach entry $ilEntries {

	if { $iLabelFunction != "" } {
	    regsub -all %s $iLabelFunction $entry sCommand
	    set sLabel [eval $sCommand]
	}

	if { $ilLabels != {} } {
	    set sLabel [lindex $ilLabels $nEntry]
	}

	$imw add command $entry -label $sLabel
	incr nEntry
    }
    
    # Renable the callback.
    $imw config -disablecallback 0
}

proc ShowHideConsole { ibShow } {
    global gaWidget
    global gaView

    if { $ibShow } {

	grid $gaWidget(tkcon) -sticky ews \
	    -column $gaWidget(tkcon,column) -row $gaWidget(tkcon,row) \
	    -columnspan 2

    } else {
	
	grid forget $gaWidget(tkcon)
    }

    # Make sure our visible var is set correctly.
    set gaView(tkcon,visible) $ibShow
}

# DATA LOADING =====================================================

proc MakeVolumeCollection { ifnVolume } {
    dputs "MakeVolumeCollection  $ifnVolume  "


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

proc MakeSurfaceCollection { ifnSurface } {
    dputs "MakeSurfaceCollection  $ifnSurface  "


    set err [catch { set colID [MakeDataCollection Surface] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    set err [catch { SetSurfaceCollectionFileName $colID $ifnSurface } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    # Get a good name for the collection.
    set sLabel [ExtractLabelFromFileName $ifnSurface]
    
    set err [catch { SetCollectionLabel $colID $sLabel } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    return $colID
}

proc Make2DMRILayer { isLabel } {
    dputs "Make2DMRILayer  $isLabel  "


    set err [catch { set layerID [MakeLayer 2DMRI] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    set err [catch { SetLayerLabel $layerID $isLabel } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    UpdateLayerList

    return $layerID
}

proc Make2DMRISLayer { isLabel } {
    dputs "Make2DMRISLayer  $isLabel  "


    set err [catch { set layerID [MakeLayer 2DMRIS] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    set err [catch { SetLayerLabel $layerID $isLabel } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    UpdateLayerList

    return $layerID
}

proc LoadVolume { ifnVolume ibCreateLayer iFrameIDToAdd } {
    dputs "LoadVolume  $ifnVolume $ibCreateLayer $iFrameIDToAdd  "


    set fnVolume [FindFile $ifnVolume]

    set err [catch { set colID [MakeVolumeCollection $fnVolume] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    if { $ibCreateLayer } {

	set sLabel [ExtractLabelFromFileName $fnVolume]

	set layerID [Make2DMRILayer "$sLabel"]

	set err [catch {
	    Set2DMRILayerVolumeCollection $layerID $colID } sResult]
	if { 0 != $err } { tkuErrorDlog $sResult; return }
	
	if { $iFrameIDToAdd != -1 } {
	    SetLayerInAllViewsInFrame $iFrameIDToAdd $layerID
	}

	UpdateCollectionList
	SelectCollectionInCollectionProperties $colID
	SelectLayerInLayerProperties $layerID
    }

    # Create a new ROI for this collection.
    set roiID [NewCollectionROI $colID]
    SetROILabel $roiID "New ROI"
    SetROIType $roiID free
    SetROIColor $roiID 0 0 255
    UpdateROIList
    SelectROIInROIProperties $roiID
    

    # Add this directory to the shortcut dirs if it isn't there already.
    AddDirToShortcutDirsList [file dirname $ifnVolume]
}

proc LoadSurface { ifnSurface ibCreateLayer iFrameIDToAdd } {
    dputs "LoadSurface  $ifnSurface $ibCreateLayer $iFrameIDToAdd  "


    set fnSurface [FindFile $ifnSurface]

    set err [catch { set colID [MakeSurfaceCollection $fnSurface] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    if { $ibCreateLayer } {

	set sLabel [ExtractLabelFromFileName $fnSurface]

	set layerID [Make2DMRISLayer "$sLabel"]

	set err [catch {
	    Set2DMRISLayerSurfaceCollection $layerID $colID } sResult]
	if { 0 != $err } { tkuErrorDlog $sResult; return }
	
	if { $iFrameIDToAdd != -1 } {
	    SetLayerInAllViewsInFrame $iFrameIDToAdd $layerID
	}

	UpdateCollectionList
	SelectCollectionInCollectionProperties $colID
	SelectLayerInLayerProperties $layerID
    }

    # Add this directory to the shortcut dirs if it isn't there already.
    AddDirToShortcutDirsList [file dirname $ifnSurface]
}

proc DoLoadVolumeDlog {} {
    dputs "DoLoadVolumeDlog  "

    global glShortcutDirs

    tkuDoFileDlog -title "Load Volume" \
	-prompt1 "Load Volume: " \
	-defaultdir1 [GetDefaultFileLocation LoadVolume] \
	-type2 checkbox \
	-prompt2 "Automatically add new layer to all views" \
	-defaultvalue2 1 \
	-shortcuts $glShortcutDirs \
	-okCmd { 
	    set frameID -1
	    if { %s2 } {
		set frameID $gFrameWidgetToID($gaWidget(scubaFrame))
	    }
	    LoadVolume %s1 1 $frameID
	}
}

proc DoSaveLabelDlog {} {
    dputs "DoSaveLabelDlog  "

    global glShortcutDirs
    global gaROI
    global gaCollection

    if { [info exists gaCollection(current,id)] &&
	 $gaCollection(current,id) >= -1 } {

	if { [info exists gaROI(current,id)] &&
	     $gaROI(current,id) >= -1 } {
	    
	    tkuDoFileDlog -title "Save Label" \
		-prompt1 "Will save ROI \"$gaROI(current,label)\" from data collection \"$gaCollection(current,label)\"" \
		-type1 note \
		-prompt2 "Save Label: " \
		-defaultdir2 [GetDefaultFileLocation SaveLabel] \
		-shortcuts $glShortcutDirs \
		-okCmd { 
		    set err [catch {
			WriteVolumeROIToLabel $gaCollection(current,id) \
			    $gaROI(current,id) %s1
		    } sResult]
		    if { 0 != $err } { tkuErrorDlog $sResult }
		}
	    
	} else {
	    tkuErrorDlog "No ROI is selected. Make sure you have selected one in the Data Collections panel."
	}
	
    } else {
	
	 tkuErrorDlog "There are no data collections. Load some data and try again."

     }
 }


proc DoLoadLabelDlog {} {
    dputs "DoLoadLabelDlog  "

    global glShortcutDirs
    global gaCollection

    if { [info exists gaCollection(current,id)] &&
	 $gaCollection(current,id) >= -1 } {

	tkuDoFileDlog -title "Load Label" \
	    -prompt1 "Will create a new ROI for collection \"$gaCollection(current,label)\"" \
	    -type1 note \
	    -prompt2 "Load Label: " \
	    -defaultdir2 [GetDefaultFileLocation LoadLabel] \
	    -shortcuts $glShortcutDirs \
	    -okCmd { 
		set err [catch {
		    set roiID [NewVolumeROIFromLabel \
				   $gaCollection(current,id) %s2]
		} sResult]
		if { 0 != $err } { 
		    tkuErrorDlog $sResult 
		} else {
		    SetROILabel $roiID [file tail %s2]
		    UpdateROIList
		    SelectROIInROIProperties $roiID
		    RedrawFrame [GetMainFrameID]
		}
	    }
	
    } else {
	
	tkuErrorDlog "There are no data collections. Load some data and try again."
	
     }
 }


 proc DoExportLabelsDlog {} {
     dputs "DoExportLabelsDlog"

     global glShortcutDirs
     global gaCollection

     if { [info exists gaCollection(current,id)] &&
	  $gaCollection(current,id) >= -1 } {

	 tkuDoFileDlog -title "Export Segmentation" \
	     -prompt1 "Will export all structure ROIs from data collection \"$gaCollection(current,label)\"" \
	     -type1 note \
	     -prompt2 "Save Volume: " \
	     -defaultdir2 [GetDefaultFileLocation ExportSegmentation] \
	     -shortcuts $glShortcutDirs \
	     -okCmd { 
		 set err [catch {
		     WriteVolumeROIsToSegmentation \
			 $gaCollection(current,id) %s1
		 } sResult]
		 if { 0 != $err } { tkuErrorDlog $sResult }
	     }

     } else {
	 tkuErrorDlog "There are no data collections. Load some data and try again."
     }

}


# MAIN =============================================================

# Look at our command line args. For some we will want to process and
# exit without bringing up all our windows. For some, we need to bring
# up our windows first. So cache those in lCommands and we'll execute
# them later.
set lCommands {}

set argc [GetArgc]
set argv [GetArgv]

set nArg 0
while { $nArg < $argc } {
    set sArg [lindex $argv $nArg]
    set sOption [string range $sArg [expr [string last "-" $sArg]+1] end]
    switch $sOption {
	v - volume {
	    incr nArg
	    set fnVolume [lindex $argv $nArg]
	    lappend lCommands "LoadVolume $fnVolume 1 [GetMainFrameID]"
	}
	f - surface {
	    incr nArg
	    set fnSurface [lindex $argv $nArg]
	    lappend lCommands "LoadSurface $fnSurface 1 [GetMainFrameID]"
	}
	s - subject {
	    incr nArg
	    set sSubject [lindex $argv $nArg]
	    lappend lCommands "SetSubjectName $sSubject"
	}
	
	help - default {
	    if {$sOption != "help"} {puts "Option $sOption not recognized."}
	    puts ""
	    puts "Usage: scuba \[OPTION\]..."
	    puts "Data viewer for the FreeSurfer package."
	    puts ""
	    puts "Options:"
	    puts "-s, --subject SUBJECT Set the subject for this session. Environment variable "
	    puts "                      SUBJECTS_DIR should be set."
	    puts "-v, --volume FILE     Load a volume file. Can be a file name or a subdir in"
	    puts "                      the subject's directory specified with -s."
	    puts "-f, --surface FILE    Load a surface file Can be a file name or a subdir in"
	    puts "                      the subject's directory specified with -s."
	    exit
	}
    }
    incr nArg
}


# Source our support files.
foreach sSourceFileName { tkUtils.tcl tkcon.tcl } {
    set lPath [list "$env(PWD)"]
    if { [info exists env(DEV)] } {
	lappend lPath "$env(DEV)/scripts"
    }
    if { [info exists env(FREESURFER_HOME)] } {
	lappend lPath "$env(MRI_DIR)/lib/tcl"
    }
    set bFound 0
    foreach sPath $lPath {
	if { $bFound == 0 } {
	    set sFullFileName [ file join $sPath $sSourceFileName ]
	    if { [file exists $sFullFileName] } { 
		set nErr [catch { source $sFullFileName } sResult]
		if { $nErr != 0 } {
		    puts "Error sourcing $sFullFileName: $sResult"
		} else {
		    set bFound 1
		}
	    }
	}
    }    
    if { $bFound == 0 } {
	puts "Couldn't load $sSourceFileName: Not found in $lPath"
    }
}


# Do some startup stuff.
BuildShortcutDirsList
LoadImages


# Make the main window.
set gaWidget(window) .main
toplevel $gaWidget(window)

# Make the tkcon panel. This must be done at this scope because the
# tkcon.tcl script needs access to some global vars.
set gaWidget(tkcon) [frame $gaWidget(window).tkcon -height 40]
::tkcon::Init -root $gaWidget(window).tkcon -showmenu 0 -embed 1 -exec ""
tkcon attach main


# Make the areas in the window. Make the scuba frame first because it
# inits stuff that is needed by other areas.
set gaWidget(scubaFrame) [MakeScubaFrame $gaWidget(window)]
set gaWidget(menuBar) [MakeMenuBar $gaWidget(window)]
set gaWidget(toolBar) [MakeToolBar $gaWidget(window)]
set gaWidget(labelArea) [MakeLabelArea $gaWidget(window)]
set gaWidget(properties) [MakePropertiesPanel $gaWidget(window)]

# Set the grid coords of our areas and the grid them in.
set gaWidget(menuBar,column)    0; set gaWidget(menuBar,row)    0
set gaWidget(toolBar,column)    0; set gaWidget(toolBar,row)    1
set gaWidget(scubaFrame,column) 0; set gaWidget(scubaFrame,row) 2
set gaWidget(labelArea,column)  0; set gaWidget(labelArea,row)  3
set gaWidget(properties,column) 1; set gaWidget(properties,row) 2
set gaWidget(tkcon,column)      0; set gaWidget(tkcon,row)      4

grid $gaWidget(menuBar) -sticky ew -columnspan 2 \
    -column $gaWidget(menuBar,column) -row $gaWidget(menuBar,row)
grid $gaWidget(toolBar) -sticky ew -columnspan 2 \
    -column $gaWidget(toolBar,column) -row $gaWidget(toolBar,row)
grid $gaWidget(scubaFrame) \
    -column $gaWidget(scubaFrame,column) -row $gaWidget(scubaFrame,row) \
    -sticky news
grid $gaWidget(labelArea) \
    -column $gaWidget(labelArea,column) -row $gaWidget(labelArea,row) \
    -sticky ew
grid $gaWidget(properties) -sticky n \
    -column $gaWidget(properties,column) -row $gaWidget(properties,row) \
    -rowspan 2

grid columnconfigure $gaWidget(window) 0 -weight 1
grid columnconfigure $gaWidget(window) 1 -weight 0
grid rowconfigure $gaWidget(window) 0 -weight 0
grid rowconfigure $gaWidget(window) 1 -weight 0
grid rowconfigure $gaWidget(window) 2 -weight 1
grid rowconfigure $gaWidget(window) 3 -weight 0
grid rowconfigure $gaWidget(window) 4 -weight 0

wm withdraw .

# Let tkUtils finish up.
tkuFinish

# Make the default color LUTs.
foreach fnLUT {tkmeditColorsCMA tkmeditParcColorsCMA surface_labels.txt Simple_surface_labels2002.txt} {
    if { [file exists $env(FREESURFER_HOME)/$fnLUT] } {
	set lutID [MakeNewColorLUT]
	SetColorLUTLabel $lutID "$fnLUT"
	SetColorLUTFileName $lutID [file join $env(FREESURFER_HOME) $fnLUT]
    }
}


# Make the default transform if necessary.
set transformList [GetTransformIDList]
if { [llength $transformList] == 0 } {
    set transformID [MakeNewTransform]
    SetTransformLabel $transformID "Identity"
}


# Set default view configuration and update/initialize the
# menus. Select the view to set everything up.
SetFrameViewConfiguration [GetMainFrameID] c1
UpdateCollectionList
UpdateLayerList
UpdateViewList
UpdateSubjectList
UpdateTransformList
UpdateLUTList
SelectViewInViewProperties 0
SelectTransformInTransformProperties 0
SelectLUTInLUTProperties 0
SelectToolInToolProperties navigation

ShowHideConsole $gaView(tkcon,visible)

GetPreferences

MakeScubaFrameBindings [GetMainFrameID]

# Now execute all the commands we cached before.
foreach command $lCommands {
    eval $command
}




bind $gaWidget(window) <Alt-Key-q> "Quit"
bind $gaWidget(window) <Alt-Key-n> {
    if { $gaView(tkcon,visible) } {
	set gaView(tkcon,visible) 0
    } else {
	set gaView(tkcon,visible) 1
    }
    ShowHideConsole $gaView(tkcon,visible)
}
bind $gaWidget(window) <Key-1> {set gaTool(current,radius) 1; SetToolBrushRadius $gaFrame([GetMainFrameID],toolID) $gaTool(current,radius) } 
bind $gaWidget(window) <Key-2> {set gaTool(current,radius) 2; SetToolBrushRadius $gaFrame([GetMainFrameID],toolID) $gaTool(current,radius) }
bind $gaWidget(window) <Key-3> {set gaTool(current,radius) 3; SetToolBrushRadius $gaFrame([GetMainFrameID],toolID) $gaTool(current,radius) }
bind $gaWidget(window) <Key-4> {set gaTool(current,radius) 4; SetToolBrushRadius $gaFrame([GetMainFrameID],toolID) $gaTool(current,radius) }
bind $gaWidget(window) <Key-5> {set gaTool(current,radius) 5; SetToolBrushRadius $gaFrame([GetMainFrameID],toolID) $gaTool(current,radius) }
bind $gaWidget(window) <Key-6> {set gaTool(current,radius) 6; SetToolBrushRadius $gaFrame([GetMainFrameID],toolID) $gaTool(current,radius) }
bind $gaWidget(window) <Key-7> {set gaTool(current,radius) 7; SetToolBrushRadius $gaFrame([GetMainFrameID],toolID) $gaTool(current,radius) }
bind $gaWidget(window) <Key-8> {set gaTool(current,radius) 8; SetToolBrushRadius $gaFrame([GetMainFrameID],toolID) $gaTool(current,radius) }
bind $gaWidget(window) <Key-9> {set gaTool(current,radius) 9; SetToolBrushRadius $gaFrame([GetMainFrameID],toolID) $gaTool(current,radius) }


