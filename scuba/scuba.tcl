
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

# gSubject

# gaMenu

# gView
#  inPlane - inPlane of current view

# gaSubject

# gaFrame
#   n - id of frame
#     viewConfi

# gaView


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
    global gsaDefaultLocation
    if { [string range $isValue 0 0] == "/" } {
	set gsaDefaultLocation($iType) $isValue
    }
}

proc SetSubjectName { isSubject } {
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
    set gaMenu(view)  $fwMenuBar.mbwView

    frame $fwMenuBar -border 2 -relief raised

    tkuMakeMenu -menu $gaMenu(file) -label "File" -items {
	{command "Load Volume..." { DoLoadVolumeDlog } }
	{separator}
	{command "Quit" { Quit } }
    }

    pack $gaMenu(file) -side left

    tkuMakeMenu -menu $gaMenu(view) -label "View" -items {
	{check "Flip Left/Right" { SetViewFlipLeftRightYZ $gaView(current,id) $gaView(flipLeftRight) } gaView(flipLeftRight) }
    }

    pack $gaMenu(view) -side left

    return $fwMenuBar
}


proc MakeToolBar { ifwTop } {
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
	-variable gView(inPlane) \
	-command {ToolBarWrapper} \
	-buttons {
	    { -type image -name r -image icon_orientation_sagittal }
	    { -type image -name a -image icon_orientation_coronal }
	    { -type image -name s -image icon_orientation_horizontal }
	}

    set gView(inPlane) r

    pack $fwToolBar.fwTool $fwToolBar.fwView $fwToolBar.fwInPlane \
	-side left

    return $fwToolBar
}

proc ToolBarWrapper { isName iValue } {
    global gaLayer
    global gaFrame
    if { $iValue == 1 } {
	switch $isName {
	    navigation - voxelEditing {
		SetToolMode $gaFrame([GetMainFrameID],toolID) $isName
	    }
	    c1 - c22 - c13 {
		SetFrameViewConfiguration [GetMainFrameID] $isName
		UpdateViewList
	    }
	    r - a - s {
		SetViewInPlane [GetSelectedViewID [GetMainFrameID]] $isName
		RedrawFrame [GetMainFrameID]
	    }
	    grayscale - heatScale - lut {
		Set2DMRILayerColorMapMethod \
		    $gaLayer(current,id) $gaLayer(current,colorMapMethod)
		RedrawFrame [GetMainFrameID]
	    }
	    nearest - trilinear - sinc {
		Set2DMRILayerSampleMethod \
		    $gaLayer(current,id) $gaLayer(current,sampleMethod)
		RedrawFrame [GetMainFrameID]
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
    global gaFrame
    global gaTool

    set fwScuba $ifwTop.fwScuba
    
    set frameID [GetNewFrameID]
    togl $fwScuba -width 512 -height 512 -rgba true -ident $frameID

    bind $fwScuba <Motion> \
	"%W MouseMotionCallback %x %y %b; ScubaMouseMotionCallback %x %y %b"
    bind $fwScuba <ButtonPress> "%W MouseDownCallback %x %y %b"
    bind $fwScuba <ButtonRelease> "%W MouseUpCallback %x %y %b"
    bind $fwScuba <KeyRelease> "%W KeyUpCallback %x %y %K"
    bind $fwScuba <KeyPress> "%W KeyDownCallback %x %y %K"
    bind $fwScuba <Enter> "focus $fwScuba"

    set gFrameWidgetToID($fwScuba) $frameID

    set gaFrame($frameID,toolID) [GetToolIDForFrame $frameID]
    set gaTool($frameID,mode) [GetToolMode $gaFrame($frameID,toolID)]

    return $fwScuba
}

proc ScubaMouseMotionCallback { inX inY iButton } {

    set err [catch { 
	set viewID [GetViewIDAtFrameLocation [GetMainFrameID] $inX $inY] 
    } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    set err [catch { 
	set labelValues [GetLabelValuesSet $viewID cursor] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    UpdateLabelArea $labelValues
}

proc Quit {} {

    SaveGlobalPreferences
    exit
}

# INTERFACE CREATION ==================================================

proc MakeLabelArea { ifwTop } {
    global gaWidget

    set fwLabelArea     $ifwTop.fwLabelArea

    frame $fwLabelArea -border 2 -relief raised

    set gaWidget(labelArea,labelValueWidgets) {}
    set gaWidget(labelArea,numberOfLabels) 0

    return $fwLabelArea
}

proc MakeNavigationArea { ifwTop } {
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
    global gaWidget
    global gaSubject

    set fwTop  $ifwTop.fwSubjects
    set fwMenu $fwTop.fwMenu
    set fwData $fwTop.fwData

    frame $fwTop

    frame $fwMenu
    tixOptionMenu $fwMenu.menu \
	-label "Subject:" \
	-variable gaSubject(current,menuIndex) \
	-command { SubjectsLoaderSubjectMenuCallback }
    set gaWidget(subjectsLoader,subjectsMenu) $fwMenu.menu
    pack $fwMenu.menu

    frame $fwData
    tixOptionMenu $fwData.volumesMenu \
	-label "Volumes:"
    set gaWidget(subjectsLoader,volumeMenu) $fwData.volumesMenu

    button $fwData.volumesButton \
	-text "Load" \
	-command {LoadVolumeFromSubjectsLoader [$gaWidget(subjectsLoader,volumeMenu) cget -value]}

    grid $fwData.volumesMenu   -column 0 -row 0 -sticky ew
    grid $fwData.volumesButton -column 1 -row 0 -sticky e


    tixOptionMenu $fwData.surfacesMenu \
	-label "Surfaces:"
    set gaWidget(subjectsLoader,surfaceMenu) $fwData.surfacesMenu

    button $fwData.surfacesButton \
	-text "Load" \
	-command {LoadSurfaceFromSubjectsLoader [$gaWidget(subjectsLoader,surfaceMenu) cget -value]}
    
    grid $fwData.surfacesMenu   -column 0 -row 1 -sticky ew
    grid $fwData.surfacesButton -column 1 -row 1 -sticky e


    grid columnconfigure $fwData 0 -weight 1
    grid columnconfigure $fwData 1 -weight 0
    
    pack $fwMenu $fwData -side top -expand yes -fill x

    return $fwTop
}

proc MakePropertiesPanel { ifwTop } {
    global gaWidget

    set fwTop  $ifwTop.fwProps
    
    tixNoteBook $fwTop
    
    $fwTop add layerPanel -label Layers
    $fwTop add viewPanel -label Views
    $fwTop add subjectsLoader -label Subjects
    $fwTop add transformPanel -label Transforms

    set gaWidget(layerProperties) \
	[MakeLayerPropertiesPanel [$fwTop subwidget layerPanel]]
    set gaWidget(viewProperties) \
	[MakeViewPropertiesPanel [$fwTop subwidget viewPanel]]
    set gaWidget(subjectsLoader) \
	[MakeSubjectsLoaderPanel [$fwTop subwidget subjectsLoader]]
    set gaWidget(transformProperties) \
	[MakeTransformsPanel [$fwTop subwidget transformPanel]]

    pack $gaWidget(layerProperties)
    pack $gaWidget(viewProperties)
    pack $gaWidget(subjectsLoader)
    pack $gaWidget(transformProperties)

    return $fwTop
}

proc MakeLayerPropertiesPanel { ifwTop } {
    global gaWidget
    global gaLayer
    global glShortcutDirs

    set fwTop        $ifwTop.fwLayerProps
    set fwMenu       $fwTop.fwMenu
    set fwProps      $fwTop.fwProps

    frame $fwTop

    frame $fwMenu
    tixOptionMenu $fwMenu.menu \
	-label "Layer:" \
	-variable gaLayer(current,menuIndex) \
	-command { LayerPropertiesMenuCallback }
    set gaWidget(layerProperties,menu) $fwMenu.menu
    pack $fwMenu.menu

    frame $fwProps
    set fwPropsCommon $fwProps.fwPropsCommon
    set fwProps2DMRI  $fwProps.fwProps2DMRI
    set fwProps2DMRIS $fwProps.fwProps2DMRIS

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
    tkuMakeFileSelector $fwProps2DMRI.fwLUT \
	-text "LUT" \
	-variable gaLayer(current,fileName) \
	-shortcutdirs [list $glShortcutDirs] \
	-defaultdir [GetDefaultFileLocation LUT] \
	-command {Set2DMRILayerFileLUTFileName $gaLayer(current,id) $gaLayer(current,fileName); RedrawFrame [GetMainFrameID]}
    tkuMakeToolbar $fwProps2DMRI.tbwSampleMethod \
	-allowzero 0 -radio 1 \
	-variable gaLayer(current,sampleMethod) \
	-command ToolBarWrapper \
	-buttons {
	    {-type text -name nearest -label "Nearest"}
	    {-type text -name trilinear -label "Trilinear"}
	    {-type text -name sinc -label "Sinc"}
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

    grid $fwProps2DMRI.tbwColorMapMethod -column 0 -row 0 -sticky ew
    grid $fwProps2DMRI.fwLUT             -column 0 -row 1 -sticky ew
    grid $fwProps2DMRI.cbwClearZero      -column 0 -row 2 -sticky ew
    grid $fwProps2DMRI.tbwSampleMethod   -column 0 -row 3 -sticky ew
    grid $fwProps2DMRI.swBC              -column 0 -row 4 -sticky ew
    grid $fwProps2DMRI.swMinMax          -column 0 -row 5 -sticky ew
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


    grid $fwPropsCommon -column 0 -row 0 -sticky news

    grid $fwMenu -column 0 -row 0 -sticky new
    grid $fwProps -column 0 -row 1 -sticky news

    return $fwTop
}

proc MakeViewPropertiesPanel { ifwTop } {
    global gaWidget
    global gaView

    set fwTop        $ifwTop.fwViewProps
    set fwMenu       $fwTop.fwMenu
    set fwProps      $fwTop.fwProps

    frame $fwTop

    frame $fwMenu
    tixOptionMenu $fwMenu.menu \
	-label "View:" \
	-variable gaView(current,menuIndex) \
	-command { ViewPropertiesMenuCallback }
    set gaWidget(viewProperties,menu) $fwMenu.menu
    pack $fwMenu.menu

    frame $fwProps
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
    
    grid $fwProps.ewID      -column 0 -row 0 -sticky nw
    grid $fwProps.ewCol     -column 1 -row 0 -sticky e
    grid $fwProps.ewRow     -column 2 -row 0 -sticky e
    grid $fwProps.cbwLinked -column 0 -row 1 -sticky ew -columnspan 3
    for { set nLevel 0 } { $nLevel < 10 } { incr nLevel } {
	grid $fwProps.mwDraw$nLevel \
	    -column 0 -row [expr $nLevel + 2] -sticky ew -columnspan 3
    }
    grid $fwProps.mwTransform  -column 0 -row 13 -sticky ew -columnspan 3
    grid $fwProps.bwCopyLayers -column 0 -row 14 -sticky ew -columnspan 3
    
    grid $fwMenu -column 0 -row 0 -sticky new
    grid $fwProps -column 0 -row 1 -sticky news

    return $fwTop
}

proc MakeTransformsPanel { ifwTop } {
    global gaWidget
    global gaTransform

    set fwTop      $ifwTop.fwSubjects
    set fwMenu     $fwTop.fwMenu
    set fwProps    $fwTop.fwProps
    set fwCommands $fwTop.fwCommands
 
    frame $fwTop

    frame $fwMenu
    tixOptionMenu $fwMenu.menu \
	-label "Transform:" \
	-variable gaTransform(current,menuIndex) \
	-command { TransformPropertiesMenuCallback }
    set gaWidget(transformProperties,menu) $fwMenu.menu
    pack $fwMenu.menu

    
    frame $fwProps
    tkuMakeEntry $fwProps.ewLabel \
	-variable gaTransform(current,label) \
	-notify 1 \
	-command {SetTransformLabel $gaTransform(current,id) $gaTransform(current,label); UpdateTransformList} 
    set gaWidget(transformProperties,labelEntry) $fwProps.ewLabel
    
    grid $fwProps.ewLabel -column 0 -row 0 -columnspan 4 -sticky ew

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
		-column $nCol -row [expr $nRow + 1] -sticky ew
	}
    }

    button $fwProps.bwSetTransform -text "Set Values" \
	-command { SetTransformValues $gaTransform(current,id) $gaTransform(current,valueList); ClearSetTransformValuesButton }
    set gaWidget(transformProperties,setValuesButton) $fwProps.bwSetTransform

    grid $fwProps.bwSetTransform -column 0 -row 5 -columnspan 4 -sticky ew

    frame $fwCommands
    button $fwCommands.bwMakeTransform -text "Make New Transform" \
	-command { set transformID [MakeNewTransform]; SetTransformLabel $transformID "New Transform"; UpdateTransformList; SelectTransformInTransformProperties $transformID }

    pack $fwCommands.bwMakeTransform -expand yes -fill x

    pack $fwMenu $fwProps $fwCommands -side top -expand yes -fill x

    return $fwTop
}


# LAYER PROPERTIES FUNCTIONS ===========================================

proc LayerPropertiesMenuCallback { inLayer } {
    global gaLayer

    # Get the ID at this index in the idList, then select that layer.
    set layerID [lindex $gaLayer(idList) $inLayer]
    SelectLayerInLayerProperties $layerID
}

proc SelectLayerInLayerProperties { iLayerID } {
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
    # callback and set the value of the menu to the index of this
    # layer ID in the layer ID list. Then reenable the callback.
    $gaWidget(layerProperties,menu) config -disablecallback 1
    $gaWidget(layerProperties,menu) config \
	-value [lsearch $gaLayer(idList) $iLayerID]
    $gaWidget(layerProperties,menu) config -disablecallback 0
    
    # Do the type specific stuff.
    switch $gaLayer(current,type) {
	2DMRI { 
	    # Pack the type panel.
	    grid $gaWidget(layerProperties,2DMRI) -column 0 -row 1 -sticky news

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
	    set gaLayer(current,fileName) \
		[Get2DMRILayerFileLUTFileName $iLayerID]
	    set gaLayer(current,minVisibleValue) \
		[Get2DMRILayerMinVisibleValue $iLayerID]
	    set gaLayer(current,maxVisibleValue) \
		[Get2DMRILayerMaxVisibleValue $iLayerID]
	}
	2DMRIS {
	    # Pack the type panel.
	    grid $gaWidget(layerProperties,2DMRIS) \
		-column 0 -row 1 -sticky news

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
    global gaLayer
    global gaWidget
    global gaView

    # We have two jobs here. First we need to populate the menu that
    # selects the current layer in the layer props panel. Then we need
    # to populate all the level-layer menus in the view props
    # panel. First do the layer props.

    # Get the layer ID list.
    set err [catch { set gaLayer(idList) [GetLayerIDList] } sResult]
    if { $err } { 
	set gaLayer(idList) {} 
    }

    # Disable the menu callback.
    $gaWidget(layerProperties,menu) config -disablecallback 1

    # Get all the entries, delete them, then add commands for all the
    # IDs in the layer ID list.
    set lEntries [$gaWidget(layerProperties,menu) entries]
    foreach entry $lEntries { 
	$gaWidget(layerProperties,menu) delete $entry
    }
    foreach id $gaLayer(idList) {
	$gaWidget(layerProperties,menu) add command $id \
	    -label "$id: [GetLayerLabel $id]"
    }
    # Renable the menu.
    $gaWidget(layerProperties,menu) config -disablecallback 0

    # Reselect the current layer.
    if { [info exists gaLayer(current,id)] && 
	 $gaLayer(current,id) >= 0 } {
	SelectLayerInLayerProperties $gaLayer(current,id)
    }


    # Populate the menus in the view props draw level menus.
    for { set nLevel 0 } { $nLevel < 10 } { incr nLevel } {

	# Disable callback.
	$gaWidget(viewProperties,drawLevelMenu$nLevel) \
	    config -disablecallback 1

	# Delete all the entries and add ones for all the IDs in the
	# ID list. Also add a command for 'none' with index of -1.
	set lEntries [$gaWidget(viewProperties,drawLevelMenu$nLevel) entries]
	foreach entry $lEntries { 
	    $gaWidget(viewProperties,drawLevelMenu$nLevel) delete $entry
	}
	$gaWidget(viewProperties,drawLevelMenu$nLevel) \
	    add command -1 -label "None"
	foreach id $gaLayer(idList) {
	    $gaWidget(viewProperties,drawLevelMenu$nLevel) add command $id \
		-label "$id: [GetLayerLabel $id]"
	}

	# Renable the callback.
	$gaWidget(viewProperties,drawLevelMenu$nLevel) \
	    config -disablecallback 0
    }

    # Make sure the right layers are selected in the view draw level
    # menus.
    UpdateCurrentViewProperties
}

# VIEW PROPERTIES FUNCTIONS =============================================

proc ViewPropertiesMenuCallback { iViewID } {
    global gaView

    SelectViewInViewProperties $iViewID
}

proc ViewPropertiesDrawLevelMenuCallback { iLevel inLayer } {
    global gaView
    global gaLayer
    
    # If we didn't get -1, find the layer ID from the list of
    # indices. Otherwise we'll set it to ID -1.
    set layerID -1
    if { $inLayer > -1 } {
	set layerID [lindex $gaLayer(idList) [expr $inLayer]]
    }

    # Set the layer in this view and redraw.
    SetLayerInViewAtLevel $gaView(current,id) $layerID $iLevel
    RedrawFrame [GetMainFrameID]
}

proc ViewPropertiesTransformMenuCallback { inTransform } {
    global gaView
    global gaTransform
    
    # Find the transform ID from the list of indices.
    set transformID [lindex $gaTransform(idList) [expr $inTransform]]

    # Set the transform in this view and redraw.
    SetViewTransform $gaView(current,id) $transformID
    RedrawFrame [GetMainFrameID]
}

proc SelectViewInViewProperties { iViewID } {
    global gaWidget
    global gaView

    # Get the general view properties from the specific view and
    # load them into the 'current' slots.
    set gaView(current,id) $iViewID
    set gaView(current,col) [GetColumnOfViewInFrame [GetMainFrameID] $iViewID]
    set gaView(current,row) [GetRowOfViewInFrame [GetMainFrameID] $iViewID]
    set gaView(current,linked) [GetViewLinkedStatus $iViewID]
 
    for { set nLevel 0 } { $nLevel < 10 } { incr nLevel } {
	set gaView(current,draw$nLevel) \
	    [GetLayerInViewAtLevel $iViewID $nLevel]
    }
    
    # Make sure that this is the item selected in the menu. Disale the
    # callback and set the value of the menu to the index of this
    # view ID in the view ID list. Then reenable the callback.
    $gaWidget(viewProperties,menu) config -disablecallback 1
    $gaWidget(viewProperties,menu) config \
	-value [lsearch $gaView(idList) $iViewID]
    $gaWidget(viewProperties,menu) config -disablecallback 0

    UpdateCurrentViewProperties
}

# This gets the layers at each level of the currently selected view
# and makes sure the draw level menus are set properly. Call it
# whenever a layer has been set in the current view.
proc UpdateCurrentViewProperties {} {
    global gaWidget
    global gaView
    global gaLayer

    for { set nLevel 0 } { $nLevel < 10 } { incr nLevel } {

	# Get the current value of this layer.
	set layerID [GetLayerInViewAtLevel $gaView(current,id) $nLevel]
	set gaView(current,draw$nLevel) $layerID

	# Disable callback.
	$gaWidget(viewProperties,drawLevelMenu$nLevel) \
	    config -disablecallback 1

	# Find the index of the layer ID at this draw level in the
	# view, and set the menu appropriately.
	$gaWidget(viewProperties,drawLevelMenu$nLevel) config \
	    -value [lsearch $gaLayer(idList) $layerID]

	# Renable the callback.
	$gaWidget(viewProperties,drawLevelMenu$nLevel) \
	    config -disablecallback 0
    }
}


# This builds the view ID list from the current view configuration and
# populates the menu that selects the view in the view props panel. It
# should be called every time the view configuration changes.
proc UpdateViewList {} {
    global gaView
    global gaWidget

    set gaView(idList) {}

    # Disable the menu.
    $gaWidget(viewProperties,menu) config -disablecallback 1

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
	}
    }

    # Empty the current view list.
    set lEntries [$gaWidget(viewProperties,menu) entries]
    foreach entry $lEntries { 
	$gaWidget(viewProperties,menu) delete $entry
    }
    
    # Add the entries from the view ID list to the menu.
    foreach id $gaView(idList) {
	set sLabel "[GetColumnOfViewInFrame [GetMainFrameID] $id], [GetRowOfViewInFrame [GetMainFrameID] $id]"
	$gaWidget(viewProperties,menu) add command $id -label $sLabel
    }

    # Reenable the menu.
    $gaWidget(viewProperties,menu) config -disablecallback 0
}

# SUBJECTS LOADER FUNCTIONS =============================================

proc SubjectsLoaderSubjectMenuCallback { inSubject } {
    global gaSubject

    # Get the name at this index in the nameList, then select that
    # subject.
    set gaSubject(current) [lindex $gaSubject(nameList) $inSubject]
    SelectSubjectInSubjectsLoader $gaSubject(current)
}

proc SelectSubjectInSubjectsLoader { isSubject } {
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
    global gaSubject

    LoadVolume "$isVolume" 1 [GetMainFrameID]
}

proc LoadSurfaceFromSubjectsLoader { isSurface } {
    global gaSubject

    LoadSurface "$isSurface" 1 [GetMainFrameID]
}

# Builds the subject nameList by looking in SUBJECTS_DIR.
proc UpdateSubjectList {} {
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

proc TransformPropertiesMenuCallback { inTransform } {
    global gaTransform

    # Get the ID at this index in the idList, then select that transform.
    set transformID [lindex $gaTransform(idList) $inTransform]
    SelectTransformInTransformProperties $transformID
}

proc SelectTransformInTransformProperties { iTransformID } {
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
    # callback and set the value of the menu to the index of this
    # transform ID in the transform ID list. Then reenable the callback.
    $gaWidget(transformProperties,menu) config -disablecallback 1
    $gaWidget(transformProperties,menu) config \
	-value [lsearch $gaTransform(idList) $iTransformID]
    $gaWidget(transformProperties,menu) config -disablecallback 0
}


# This builds the transform ID list and populates the menu that selects
# the current transform in the transform props panel, and the menu in the
# view props panel. It should be called whenever a transform is created or
# deleted, or when a lyer is added to or removed from a view.
proc UpdateTransformList {} {
    global gaTransform
    global gaWidget
    global gaView

    # Get the transform ID list.
    set err [catch { set gaTransform(idList) [GetTransformIDList] } sResult]
    if { $err } { 
	set gaTransform(idList) {} 
    }

    # First rebuild the transform list in the transform props panel.

    # Disable the menu callback.
    $gaWidget(transformProperties,menu) config -disablecallback 1

    # Get all the entries, delete them, then add commands for all the
    # IDs in the transform ID list.
    set lEntries [$gaWidget(transformProperties,menu) entries]
    foreach entry $lEntries { 
	$gaWidget(transformProperties,menu) delete $entry
    }
    foreach id $gaTransform(idList) {
	$gaWidget(transformProperties,menu) add command $id \
	    -label "$id: [GetTransformLabel $id]"
    }
    # Renable the menu.
    $gaWidget(transformProperties,menu) config -disablecallback 0

    # Reselect the current transformProperties.
    if { [info exists gaTransform(current,id)] && 
	 $gaTransform(current,id) >= 0 } {
	SelectTransformInTransformProperties $gaTransform(current,id)
    }

    # Now rebuild the transform list in the view props panel.

    # Disable callback.
    $gaWidget(viewProperties,transformMenu) \
	config -disablecallback 1
    
    # Delete all the entries and add ones for all the IDs in the
    # ID list.
    set lEntries [$gaWidget(viewProperties,transformMenu) entries]
    foreach entry $lEntries { 
	$gaWidget(viewProperties,transformMenu) delete $entry
    }
    foreach id $gaTransform(idList) {
	$gaWidget(viewProperties,transformMenu) add command $id \
	    -label "$id: [GetTransformLabel $id]"
    }
    
    # Renable the callback.
    $gaWidget(viewProperties,transformMenu) \
	    config -disablecallback 0
}

proc UpdateCurrentTransformValueList {} {
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
    global gaWidget

    # Change the set button to normal.
    $gaWidget(transformProperties,setValuesButton) config -fg black
}

# LABEL AREA FUNCTIONS ==================================================

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

    if { $nLabel > $gaWidget(labelArea,numberOfLabels) } {
	set gaWidget(labelArea,numberOfLabels) $nLabel
    }
}

# VIEW CONFIGURATION ==================================================

proc SetLayerInAllViewsInFrame { iFrameID iLayerID } {

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

# DATA LOADING =====================================================

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

proc MakeSurfaceCollection { ifnSurface } {

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

    set err [catch { set layerID [MakeLayer 2DMRI] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    set err [catch { SetLayerLabel $layerID $isLabel } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    UpdateLayerList

    return $layerID
}

proc Make2DMRISLayer { isLabel } {

    set err [catch { set layerID [MakeLayer 2DMRIS] } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    set err [catch { SetLayerLabel $layerID $isLabel } sResult]
    if { 0 != $err } { tkuErrorDlog $sResult; return }

    UpdateLayerList

    return $layerID
}

proc LoadVolume { ifnVolume ibCreateLayer iFrameIDToAdd } {

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

	SelectLayerInLayerProperties $layerID
    }

    # Add this directory to the shortcut dirs if it isn't there already.
    AddDirToShortcutDirsList [file dirname $ifnVolume]
}

proc LoadSurface { ifnSurface ibCreateLayer iFrameIDToAdd } {

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

	SelectLayerInLayerProperties $layerID
    }

    # Add this directory to the shortcut dirs if it isn't there already.
    AddDirToShortcutDirsList [file dirname $ifnSurface]
}

proc DoLoadVolumeDlog {} {
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

    # Future options.
#       -type2 checkbox \
#	-prompt2 "Automatically create new layer" \
#	-defaultvalue2 1 \
}


# MAIN =============================================================

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
	r - surface {
	    incr nArg
	    set fnSurface [lindex $argv $nArg]
	    lappend lCommands "LoadSurface $fnSurface 1 [GetMainFrameID]"
	}
	j - subject {
	    incr nArg
	    set sSubject [lindex $argv $nArg]
	    lappend lCommands "SetSubjectName $sSubject"
	}
	s - segmentation - seg {
	    
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

# Do some startup stuff.
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
set gaWidget(properties,column)  1; set gaWidget(properties,row) 2

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
    -column $gaWidget(properties,column) -row $gaWidget(properties,row)

grid columnconfigure $gaWidget(window) 0 -weight 1
grid columnconfigure $gaWidget(window) 1 -weight 0
grid rowconfigure $gaWidget(window) 0 -weight 0
grid rowconfigure $gaWidget(window) 1 -weight 0
grid rowconfigure $gaWidget(window) 2 -weight 1
grid rowconfigure $gaWidget(window) 3 -weight 0

wm withdraw .

# Let tkUtils finish up.
tkuFinish

# Set default view configuration and update/initialize the
# menus. Select the view to set everything up.
SetFrameViewConfiguration [GetMainFrameID] c1
UpdateLayerList
UpdateViewList
UpdateSubjectList
UpdateTransformList
SelectViewInViewProperties 0
SelectTransformInTransformProperties 0

# Now execute all the commands we cached before.
foreach command $lCommands {
    eval $command
}
