##
## tkmedit.tcl
##
##
## Copyright (C) 2000-2011, CorTechs Labs, Inc. (La Jolla, CA) and
## The General Hospital Corporation (Boston, MA).
##
## Terms and conditions for use, reproduction, distribution and contribution
## are found in the 'FreeSurfer/CorTechs Software License Agreement' contained
## in the file 'license.cortechs.txt' found in the FreeSurfer distribution,
## and here:
##
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferCorTechsLicense
##
## Reporting: freesurfer@nmr.mgh.harvard.edu
##

source $env(FREESURFER_HOME)/tktools/tkm_common.tcl

# This function finds a file from a list of directories.
proc FindFile { ifnFile ilDirs } {
    foreach sPath $ilDirs {
	set sFullFileName [ file join $sPath $ifnFile ]
	if { [file readable $sFullFileName] } {
	    puts "Reading $sFullFileName"
	    return $sFullFileName
	}
    }
    puts "Couldn't find $ifnFile: Not in $ilDirs"
    return ""
}

# Try to get some default script locations from environment variables.
set sDefaultScriptsDir ""
catch { set sDefaultScriptsDir "$env(FREESURFER_HOME)/lib/tcl" }
set sTkToolsDir ""
catch { set sTkToolsDir "$env(FREESURFER_HOME)/tktools" }
set sTkmeditScriptsDir ""
catch { set sTkmeditScriptsDir "$env(TKMEDIT_SCRIPTS_DIR)" }
set sFsgdfDir ""
catch { set sFsgdfDir "$env(FSGDF_DIR)" }


# Source the tkm_common.tcl and tkm_wrappers.tcl
set fnCommon \
    [FindFile tkm_common.tcl \
	 [list $sTkmeditScriptsDir "." "../scripts" $sDefaultScriptsDir $sTkToolsDir]]
if { [string compare $fnCommon ""] == 0 } { exit }
source $fnCommon

set fnWrappers \
    [FindFile tkm_wrappers.tcl \
	 [list $sTkmeditScriptsDir "." "../scripts" $sDefaultScriptsDir $sTkToolsDir]]
if { [string compare $fnWrappers ""] == 0 } { exit }
source $fnWrappers

set fnFsgdf \
    [FindFile fsgdfPlot.tcl \
	 [list $sFsgdfDir "." "../scripts" $sDefaultScriptsDir $sTkToolsDir]]
if { [string compare $fnFsgdf ""] == 0 } { exit }
source $fnFsgdf


# constants
set ksWindowName "TkMedit Tools"
set ksImageDir   "$env(FREESURFER_HOME)/lib/images/"

# mri_tOrientation
set mri_tOrientation_Coronal    0
set mri_tOrientation_Horizontal 1
set mri_tOrientation_Sagittal   2

# DspA_tDisplayFlag
set glDisplayFlag { \
  flag_AuxVolume \
  flag_Anatomical \
  flag_Cursor \
  flag_MainSurface \
  flag_OriginalSurface \
  flag_PialSurface \
  flag_InterpolateSurfaceVertices \
  flag_DisplaySurfaceVertices \
  flag_ControlPoints \
  flag_Selection \
  flag_FunctionalOverlay \
  flag_FunctionalColorScaleBar \
  flag_MaskToFunctionalOverlay \
  flag_MaskFunctionalOverlayToAux \
  flag_HistogramPercentChange \
  flag_SegmentationVolumeOverlay \
  flag_AuxSegmentationVolumeOverlay \
  flag_SegLabelVolumeCount \
  flag_DTIOverlay \
  flag_VectorField \
  flag_FocusFrame \
  flag_UndoVolume \
  flag_Axes \
  flag_MaxIntProj \
  flag_HeadPoints \
  flag_VerboseGCADump
}
set nFlagIndex 1
foreach flag $glDisplayFlag {
    set gnFlagIndex($flag) $nFlagIndex
    incr nFlagIndex 
}
set glActiveFlags {}

# DspA_tTool
set DspA_tTool_Navigate   0
set DspA_tTool_Select     1
set DspA_tTool_Edit       2
set DspA_tTool_EditSeg   3
set DspA_tTool_CtrlPts    4
set DspA_tTool_Line       5

# DspA_tBrush
set DspA_tBrush_EditOne 0
set DspA_tBrush_EditTwo 1

set ksaBrushString(0) "Button 2"
set ksaBrushString(1) "Button 3"

# DspA_tBrushTarget
set DspA_tBrushTarget(main)    0
set DspA_tBrushTarget(aux)     1
set DspA_tBrushTarget(mainaux) 2

# DspA_tBrushShape
set DspA_tBrushShape_Circle 0
set DspA_tBrushShape_Square 1

# DspA_tBrushMode
set DspA_tBrushMode(set)   0
set DspA_tBrushMode(clone) 1

# DspA_tMarker
set DspA_tMarker_Crosshair 0
set DspA_tMarker_Diamond   1

# view presets
set MWin_tLinkPolicy_None                  0
set MWin_tLinkPolicy_MultipleOrientations  1
set MWin_tLinkPolicy_Mosaic                2

set tViewPreset_Single   0
set tViewPreset_Multiple 1
set tViewPreset_Mosaic   2

set ksaViewPresetString(0) "single"
set ksaViewPresetString(1) "multiple"
set ksaViewPresetString(2) "mosaic"

# tkm_tVolumeType
set tkm_tVolumeType(main) 0
set tkm_tVolumeType(aux)  1

# tkm_tSegType
set tkm_tSegType(main) 0
set tkm_tSegType(aux)  1

# tkm_tVolumeTarget
set tkm_tVolumeTarget_MainAna 0
set tkm_tVolumeTarget_AuxAna  1
set tkm_tVolumeTarget_MainSeg 2
set tkm_tVolumeTarget_AuxSeg  3

set ksaDisplayedVolumeString(0) "main"
set ksaDisplayedVolumeString(1) "aux"

# tkm_tSurfaceType
set tkm_tSurfaceType(main) 0
set tkm_tSurfaceType(aux)  1

set ksaSurfaceTypeString(0) "Main"
set ksaSurfaceTypeString(1) "Aux"

# Surf_tVertexSet
set Surf_tVertexSet(main)      0
set Surf_tVertexSet(original)  1
set Surf_tVertexSet(pial) 2

set ksaSurfaceVertexSetString(0) "Main Surface"
set ksaSurfaceVertexSetString(1) "Original Surface"
set ksaSurfaceVertexSetString(2) "Pial Surface"

# tFunctionalVolume
set tFunctionalVolume_Overlay    0
set tFunctionalVolume_TimeCourse 1

# mri_tCoordSpace
set mri_tCoordSpace_VolumeIdx  0
set mri_tCoordSpace_SurfaceRAS 1
set mri_tCoordSpace_RAS        2
set mri_tCoordSpace_Talairach  3

# Volm_tSampleType
set Volm_tSampleType(nearest)   0
set Volm_tSampleType(trilinear) 1
set Volm_tSampleType(sinc)      2

# Volm_tResampleMethod
set Volm_tResampleMethod(RAS)   0
set Volm_tResampleMethod(slice) 1

# FunD_tRegistrationType
set FunD_tRegistration(file) 0
set FunD_tRegistration(find) 1
set FunD_tRegistration(identity) 2
set FunD_tRegistration(noneNeeded) 3

set ksaLinkedCursorString(0) notlinked
set ksaLinkedCursorString(1) linked

# current location
set gOrientation 0
set gbLinkedCursor 1
set gbLinkedCursorString $ksaLinkedCursorString($gbLinkedCursor)
set gnVolX(cursor) 0
set gnVolY(cursor) 0
set gnVolZ(cursor) 0
set gnVolX(mouseover) 0
set gnVolY(mouseover) 0
set gnVolZ(mouseover) 0
set gnVolSlice 0
set gnZoomLevel 0

# for tool setting and buttons
set gTool $DspA_tTool_Select

set gDisplayIntermediateResults 1

# seg edit brush
set gSegBrush(color) 0
set gSegBrush(3d) 0
set gSegBrush(fuzzy) 0
set gSegBrush(distance) 0
set gSegBrush(src) $tkm_tVolumeType(main)
set glSegEditColors {}

# flood select params
set gFloodSelectParams(3d) 0
set gFloodSelectParams(src) $tkm_tVolumeType(main)
set gFloodSelectParams(fuzzy) 0
set gFloodSelectParams(distance) 0


# for view preset setting and buttons
set gDisplayCols 1
set gDisplayRows 1
set gViewPreset $tViewPreset_Single
set gViewPresetString $ksaViewPresetString($gViewPreset)

# display flags
foreach flag $glDisplayFlag {
    set gbDisplayFlag($flag) 0
}

# tool bar frames
set gfwaToolBar(main)  ""
set gfwaToolBar(nav)   ""
set gfwaToolBar(recon) ""

# labels
set glLabel { \
  kLabel_Coords_Vol \
  kLabel_Coords_Vol_RAS \
  kLabel_Coords_Vol_Scanner \
  kLabel_Coords_Vol_MNI \
  kLabel_Coords_Vol_Tal \
  kLabel_Value_Vol \
  kLabel_Value_Aux \
  kLabel_Coords_Func \
  kLabel_Coords_Func_RAS \
  kLabel_Value_Func \
  kLabel_Label_SegLabel \
  kLabel_Label_AuxSegLabel \
  kLabel_Label_Head \
  kLabel_SurfaceDistance \
  kLabel_LineLength }
foreach label $glLabel {
    set gfwaLabel($label,cursor) ""
    set gfwaLabel($label,mouseover) ""
}

set gsaLabelContents(kLabel_Coords_Vol,name) "Volume index"
set gsaLabelContents(kLabel_Coords_Vol_RAS,name) "Volume RAS"
set gsaLabelContents(kLabel_Coords_Vol_Scanner,name) "Volume Scanner RAS"
set gsaLabelContents(kLabel_Coords_Vol_MNI,name) "MNI Talairach"
set gsaLabelContents(kLabel_Coords_Vol_Tal,name) "Talairach"
set gsaLabelContents(kLabel_Coords_Func,name) "Functional index"
set gsaLabelContents(kLabel_Coords_Func_RAS,name) "Functional RAS"
set gsaLabelContents(kLabel_Value_Func,name) "Functional value"
set gsaLabelContents(kLabel_Label_SegLabel,name) "Sgmtn label"
set gsaLabelContents(kLabel_Label_AuxSegLabel,name) "Aux Sgmtn label"
set gsaLabelContents(kLabel_Label_Head,name) "Head Point"
set gsaLabelContents(kLabel_SurfaceDistance,name) "Surface Distance"
set gsaLabelContents(kLabel_LineLength,name) "Line Length"

foreach label $glLabel {
    set gsaLabelContents($label,value,cursor) "none"
    set gsaLabelContents($label,value,mouseover) "none"
}  


# brush info
set gBrushInfo(target)   $DspA_tBrushTarget(main)
set gBrushInfo(radius)   1
set gBrushInfo(shape)    $DspA_tBrushShape_Square
set gBrushInfo(3d)       true
set gBrushInfo(3dfill)   false
set gBrushInfo(fuzzy)    0
set gBrushInfo(distance) 0

foreach tool "$DspA_tBrush_EditOne $DspA_tBrush_EditTwo" {
    set gEditBrush($tool,mode)  0
    set gEditBrush($tool,new)  0
    set gEditBrush($tool,cloneSource)  0
    set gEditBrush($tool,high) 0
    set gEditBrush($tool,low)  0
}

# cursor
set gCursor(color,red) 0
set gCursor(color,green) 0
set gCursor(color,blue) 0
set gCursor(shape) $DspA_tMarker_Crosshair

# surface
foreach surface "$tkm_tSurfaceType(main) $tkm_tSurfaceType(aux)" {
    foreach set "$Surf_tVertexSet(main) $Surf_tVertexSet(original) \
                 $Surf_tVertexSet(pial)" {
	set gSurface($surface,$set,width) 0
	set gSurface($surface,$set,color,red)   0
	set gSurface($surface,$set,color,green) 0
	set gSurface($surface,$set,color,blue)  0
    }
}

set gbUseRealRAS 0

# general volume info
foreach volume "$tkm_tVolumeType(main) $tkm_tVolumeType(aux)" {
    set gVolume($volume,colorScale,brightness) 0
    set gVolume($volume,colorScale,contrast) 0
    set gVolume($volume,colorScale,min) 0
    set gVolume($volume,colorScale,max) 0

    set gVolume($volume,minValue) 0
    set gVolume($volume,maxValue) 0

    set gVolume($volume,sampleType) $Volm_tSampleType(nearest)

    set gVolume($volume,resampleMethod) $Volm_tResampleMethod(RAS)
}

# initialize global vars
set gsSubjectDirectory "/"
set gsSegmentationColorTable ""
set gbVolumeDirty 0
set gbAuxVolumeDirty 0
set gbTalTransformPresent 0
set gsSurfaceHemi($tkm_tSurfaceType(main)) lh
set gsSurfaceHemi($tkm_tSurfaceType(aux)) rh

# determine the list of shortcut dirs for the file dlog boxes
proc BuildShortcutDirsList {} {
    global glShortcutDirs gsSubjectDirectory env
    set glShortcutDirs {}
    if { [info exists env(SUBJECTS_DIR)] } {
	lappend glShortcutDirs $env(SUBJECTS_DIR)
    }
    if { [info exists gsSubjectDirectory] } {
	lappend glShortcutDirs $gsSubjectDirectory
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
BuildShortcutDirsList

proc BuildVolumeLoadList {} {
    global gaLoadMenu
    global gsSubjectDirectory

    # Clear the menus.
    $gaLoadMenu(main) delete 0 end
    $gaLoadMenu(aux) delete 0 end
	
    # Go through looking for file names.
    foreach {fn label} {
	brain/COR-.info "Load brain COR"
	brain/brain.mgz "Load brain.mgz"
	brain/brain.mgh "Load brain.mgh"
	filled/COR-.info "Load filled COR"
	filled/filled.mgz "Load filled.mgz"
	filled/filled.mgh "Load filled.mgh"
	orig/COR-.info "Load orig COR"
	orig/orig.mgz "Load orig.mgz"
	orig/orig.mgh "Load orig.mgh"
	T1/COR-.info "Load T1 COR"
	T1/T1.mgz "Load T1.mgz"
	T1/T1.mgh "Load T1.mgh"
	wm/COR-.info "Load wm COR"
	wm/wm.mgz "Load wm.mgz"
	wm/wm.mgh "Load wm.mgh"
    } {

	# Create the full path using the subject dir and mri/. Then
	# look fir the file. If it exists, make a menu command for it.
	set fnFull [file join $gsSubjectDirectory mri $fn]
	if { [file exists $fnFull] } {
	    $gaLoadMenu(main) add command \
		-command "LoadVolume $fnFull" \
		-label $label
	    $gaLoadMenu(aux) add command \
		-command "LoadAuxVolume $fnFull" \
		-label $label
	}
    }
}

# ========================================================= UPDATES FROM MEDIT

proc UpdateLinkedCursorFlag { ibLinked } {
    global gbLinkedCursor 
    set gbLinkedCursor $ibLinked
}

proc UpdateVolumeCursor { iSet inX inY inZ } {
    global gnVolX gnVolY gnVolZ gsaLabelContents
    set gsaLabelContents(kLabel_Coords_Vol,value,$iSet) \
      "$inX $inY $inZ"
    # set the volume coords
    set gnVolX($iSet) $inX
    set gnVolY($iSet) $inY
    set gnVolZ($iSet) $inZ
}

proc UpdateVolumeSlice { inSlice } {
    global gnVolSlice
    set gnVolSlice $inSlice
}

proc UpdateRASCursor { iSet ifX ifY ifZ } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Coords_Vol_RAS,value,$iSet) \
      "$ifX $ifY $ifZ"
}

proc UpdateTalCursor { iSet ifX ifY ifZ } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Coords_Vol_Tal,value,$iSet) \
      "$ifX $ifY $ifZ"
}

proc UpdateScannerCursor { iSet ifX ifY ifZ } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Coords_Vol_Scanner,value,$iSet) \
      "$ifX $ifY $ifZ"
}

proc UpdateMNICursor { iSet ifX ifY ifZ } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Coords_Vol_MNI,value,$iSet) \
      "$ifX $ifY $ifZ"
}

proc UpdateVolumeName { isName } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Value_Vol,name) $isName
}

proc UpdateVolumeValue { iSet inValue } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Value_Vol,value,$iSet) $inValue
}

proc UpdateAuxVolumeName { isName } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Value_Aux,name) $isName
}

proc UpdateAuxVolumeValue { iSet inValue } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Value_Aux,value,$iSet) $inValue
}

proc UpdateSegLabel { iSet isLabel } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Label_SegLabel,value,$iSet) $isLabel
}

proc UpdateAuxSegLabel { iSet isLabel } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Label_AuxSegLabel,value,$iSet) $isLabel
}

proc UpdateHeadPointLabel { iSet isLabel } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Label_Head,value,$iSet) $isLabel
}

proc UpdateFunctionalValue { iSet ifValue } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Value_Func,value,$iSet) $ifValue
}

proc UpdateFunctionalCoords { iSet inX inY inZ } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Coords_Func,value,$iSet) \
      "$inX $inY $inZ"
}

proc UpdateFunctionalRASCoords { iSet inX inY inZ } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_Coords_Func_RAS,value,$iSet) \
      "$inX $inY $inZ"
}

proc UpdateSurfaceDistance { iSet ifDistance } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_SurfaceDistance,value,$iSet) $ifDistance
}

proc UpdateLineLength { iSet ifLength } {
    global gsaLabelContents
    set gsaLabelContents(kLabel_LineLength,value,$iSet) $ifLength
}

proc UpdateZoomLevel { inLevel } { 
    global gnZoomLevel
    set gnZoomLevel $inLevel
}

proc UpdateOrientation { iOrientation } {
    global gOrientation
    set gOrientation $iOrientation
}

proc UpdateDisplayFlag { iFlagIndex ibValue } {
    global gbDisplayFlag glDisplayFlag gnFlagIndex
    global glActiveFlags
    foreach flag $glDisplayFlag {
  if { $gnFlagIndex($flag) == $iFlagIndex } {
      set gbDisplayFlag($flag) $ibValue

      # put or remove the flag from a list of active flags.
      # this will only work if anyone is listening for our flags,
      # i.e. a toolbar
      set nIndex [lsearch $glActiveFlags $flag]
      if { $ibValue == 0 } {
    if { $nIndex >= 0 } {
        catch {set glActiveFlags [lreplace $glActiveFlags $nIndex $nIndex]} sResult
    }
      } else {
    if { $nIndex == -1 } {
        catch {lappend glActiveFlags $flag} sResult
    }
      }
  }
    }
}

proc UpdateTool { iTool } {
    global gTool gToolString
    set gTool $iTool
}

proc UpdateBrushTarget { iTarget } {
    global gBrushInfo
    set gBrushInfo(target) $iTarget
}

proc UpdateBrushShape { inRadius iShape ib3D } {
    global gBrushInfo
    set gBrushInfo(radius) $inRadius
    set gBrushInfo(shape)  $iShape
    set gBrushInfo(3d)     $ib3D
}

proc UpdateBrushInfo { inBrush inLow inHigh inNewValue inMode inCloneSource } {
    global gEditBrush
    set gEditBrush($inBrush,low)          $inLow
    set gEditBrush($inBrush,high)         $inHigh
    set gEditBrush($inBrush,new)          $inNewValue
    set gEditBrush($inBrush,mode)         $inMode
    set gEditBrush($inBrush,cloneSource)  $inCloneSource
}

proc UpdateAnatomicalFillInfo { ib3D iFuzzy iDistance } {
    global gBrushInfo
    set gBrushInfo(3dfill)   $ib3D
    set gBrushInfo(fuzzy)    $iFuzzy
    set gBrushInfo(distance) $iDistance
}

proc UpdateCursorColor { ifRed ifGreen ifBlue } {
    global gCursor

    # Need to go from 0-1 based colors to 256 based colors.
    set gCursor(color,red) [expr round($ifRed * 255)]
    set gCursor(color,blue) [expr round($ifBlue * 255)]
    set gCursor(color,green) [expr round($ifGreen * 255)]
}

proc UpdateCursorShape { iShape } {
    global gCursor
    set gCursor(shape) $iShape
}

proc UpdateSurfaceLineWidth { iSurface iSet inWidth } {
    global gSurface
    set gSurface($iSurface,$iSet,width) $inWidth
}

proc UpdateSurfaceLineColor { iSurface iSet ifRed ifGreen ifBlue } {
    global gSurface
    
    # Need to go from 0-1 based colors to 256 based colors.
    set gSurface($iSurface,$iSet,color,red) [expr round($ifRed * 255)]
    set gSurface($iSurface,$iSet,color,green) [expr round($ifGreen * 255)]
    set gSurface($iSurface,$iSet,color,blue) [expr round($ifBlue * 255)]
}

proc UpdateUseRealRAS { ibUseRealRAS } {
    global gbUseRealRAS

    set gbUseRealRAS $ibUseRealRAS
}

proc UpdateSegBrushInfo { inColor ib3D iSrc iFuzzy iDistance } {
    global gSegBrush
    global glSegEditColors

    set oldSelection  $gSegBrush(color)

    set gSegBrush(color)   $inColor
    set gSegBrush(3d)      $ib3D
    set gSegBrush(src)     $iSrc
    set gSegBrush(fuzzy)   $iFuzzy
    set gSegBrush(sitance) $iDistance

    # if the seg brush info dialog box is open, we want to select the
    # item with the index of the seg brush color. do all this in a catch
    # because if the dialog is not open, this will fail. 
    catch {
	
	# We got a structure index, but we need to find the
	# corresponding list index in our glSegEditColors. To do this,
	# search the glSegEditColors for the color, and take the
	# index/2.
	set nListIndex [expr [lsearch -exact $glSegEditColors $inColor] / 2]
	
	set fwColor [.wwEditSegBrushInfoDlog.lfwColor subwidget frame].fwColor
	$fwColor subwidget listbox selection clear $oldSelection
	$fwColor subwidget listbox selection set $nListIndex
	$fwColor subwidget listbox see $nListIndex
    } sResult
}

proc UpdateSegmentationVolumeAlpha { ifAlpha } {
    global gfSegmentationVolumeAlpha
    set gfSegmentationVolumeAlpha $ifAlpha
}

proc UpdateFloodSelectParams { ib3D iSrc iFuzzy iDistance } {
    global gFloodSelectParams

    set gFloodSelectParams(3d) $ib3D
    set gFloodSelectParams(src) $iSrc
    set gFloodSelectParams(fuzzy) $iFuzzy
    set gFloodSelectParams(distance) $iDistance
}

proc UpdateVolumeColorScaleInfo { inVolume inBrightness inContrast
				  inMin inMax } {
    global gVolume

    set gVolume($inVolume,colorScale,brightness) $inBrightness
    set gVolume($inVolume,colorScale,contrast) $inContrast
    set gVolume($inVolume,colorScale,min) $inMin
    set gVolume($inVolume,colorScale,max) $inMax
}

proc UpdateVolumeSampleType { inVolume iType } {
    global gVolume
    set gVolume($inVolume,sampleType) $iType
}

proc UpdateVolumeResampleMethod { inVolume iMethod } {
    global gVolume
    set gVolume($inVolume,resampleMethod) $iMethod
}

proc UpdateDTIVolumeAlpha { ifAlpha } {
    global gfDTIVolumeAlpha
    set gfDTIVolumeAlpha $ifAlpha
}

proc UpdateTimerStatus { ibOn } {
    global gbTimerOn
    set gbTimerOn $ibOn
}

proc UpdateVolumeDirty { ibDirty } {
    global gbVolumeDirty
    set gbVolumeDirty $ibDirty
}

proc UpdateAuxVolumeDirty { ibDirty } {
    global gbAuxVolumeDirty
    set gbAuxVolumeDirty $ibDirty
}

proc UpdateSubjectDirectory { isSubjectDir } {
    global gsSubjectDirectory
    set gsSubjectDirectory $isSubjectDir
    BuildShortcutDirsList
    BuildVolumeLoadList
}

proc UpdateSegmentationColorTable { isColorTable } {
    global gsSegmentationColorTable
    set gsSegmentationColorTable $isColorTable
}

proc UpdateVolumeIsConformed { ibIsConformed } {
    global gfwaToolBar
    global gControlPointsMenuCommand

    # Enable or diable the controlpoints tool button and menu item.
    if { $ibIsConformed } {
	[$gfwaToolBar(main).fwTools.tbw subwidget 4] configure -state normal
	$gControlPointsMenuCommand(menu) entryconfigure \
	    $gControlPointsMenuCommand(entry) -state normal
	tkm_SetEnableGroupStatus tMenuGroup_ControlPoints 1
    } else {
	[$gfwaToolBar(main).fwTools.tbw subwidget 4] configure -state disabled
	$gControlPointsMenuCommand(menu) entryconfigure \
	    $gControlPointsMenuCommand(entry) -state disabled
	tkm_SetEnableGroupStatus tMenuGroup_ControlPoints 0
    }
}

proc SendDisplayFlagValue { iFlag } {
    global gnFlagIndex gbDisplayFlag
    SetDisplayFlag $gnFlagIndex($iFlag) $gbDisplayFlag($iFlag)
}

proc SendLinkedCursorValue { } {
    global gbLinkedCursor
    SetLinkedCursorFlag $gbLinkedCursor
}

proc SendSurfaceInformation { } {
    global tkm_tSurfaceType Surf_tVertexSet
    global gSurface

    foreach surface "$tkm_tSurfaceType(main) $tkm_tSurfaceType(aux)" {
	foreach set "$Surf_tVertexSet(main) $Surf_tVertexSet(original) \
                     $Surf_tVertexSet(pial)" {

	    SetSurfaceLineWidth $surface $set $gSurface($surface,$set,width)
	    
	    # Need to go from 256-based colors to 0-1 colors
	    SetSurfaceLineColor $surface $set \
		[expr $gSurface($surface,$set,color,red).0 / 255.0]\
		[expr $gSurface($surface,$set,color,green).0 / 255.0] \
		[expr $gSurface($surface,$set,color,blue).0 / 255.0]
	}
    }
}

proc SendCursorConfiguration {} {
    global gCursor
	    
    # Need to go from 256-based colors to 0-1 colors
    SetCursorColor \
	[expr $gCursor(color,red).0 / 255.0] \
	[expr $gCursor(color,green).0 / 255.0] \
	[expr $gCursor(color,blue).0 / 255.0]
    SetCursorShape $gCursor(shape)
}

proc SendVolumeSampleType { iVolume } {
    global gVolume
    SetVolumeSampleType $iVolume $gVolume($iVolume,sampleType)
}

proc SendVolumeResampleMethod { iVolume } {
    global gVolume
    SetVolumeResampleMethod $iVolume $gVolume($iVolume,resampleMethod)
}

proc SendUseRealRAS { } {
    global gbUseRealRAS
    SetUseRealRAS $gbUseRealRAS
}

proc UpdateVolumeValueMinMax { iVolume iMin iMax } {
    global gVolume
    global gInterface

    set gVolume($iVolume,minValue) $iMin
    set gVolume($iVolume,maxValue) $iMax

    # Change the value that the button sets the value to. (This will
    # fail silently if the dlog is not open.)
    catch { $gInterface(colorScaleDlog,minValueButton,$iVolume) 
	config -cmd "set gVolume($iVolume,colorScale,min) $gVolume($iVolume,minValue); SendVolumeMinMax" }
    catch { $gInterface(colorScaleDlog,maxValueButton,$iVolume) 
	config -cmd "set gVolume($iVolume,colorScale,max) $gVolume($iVolume,maxValue); SendVolumeMaxMax" }
}

proc UpdateSurfaceHemi { iSurfaceType isHemi } {
    global gsSurfaceHemi
    set gsSurfaceHemi($iSurfaceType) $isHemi
}

# =============================================================== DIALOG BOXES

proc GetDefaultLocation { iType } {
    global gsaDnefaultLocation 
    global gsSubjectDirectory gsSegmentationColorTable env
    global gsSurfaceHemi tkm_tSurfaceType
    if { [info exists gsaDefaultLocation($iType)] == 0 } {
	switch $iType {
	    LoadVolume - LoadAuxVolume - SaveVolumeAs - SaveAuxVolumeAs -
	    LoadSegmentation - LoadAuxSegmentation - SaveSegmentationAs -
	    SaveAuxSegmentationAs - ExportChangedSegmentationVolume -
	    ExportAuxChangedSegmentationVolume {
		if { $gsSubjectDirectory != "/" } {
		    set gsaDefaultLocation($iType) $gsSubjectDirectory/mri
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}	       
	    }
	    LoadVolumeDisplayTransform - LoadAuxVolumeDisplayTransform  { 
		if { $gsSubjectDirectory != "/" } {
		    set gsaDefaultLocation($iType) \
			$gsSubjectDirectory/mri/transforms
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}
	    }
	    SaveLabelAs - LoadLabel - ImportSurfaceAnnotation - 
	    WriteLineLabel { 
		if { $gsSubjectDirectory != "/" } {
		    set gsaDefaultLocation($iType) $gsSubjectDirectory/label
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}
	    }
	    LoadPialSurface {
		if { $gsSubjectDirectory != "/" } {
		    set gsaDefaultLocation($iType) $gsSubjectDirectory/surf
		    set fn [file join $gsaDefaultLocation($iType) \
				$gsSurfaceHemi($tkm_tSurfaceType(main)).pial]
		    if { [file exists $fn] } {
			set gsaDefaultLocation($iType) $fn
		    }
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}
	    }
	    LoadOriginalSurface {
		if { $gsSubjectDirectory != "/" } {
		    set gsaDefaultLocation($iType) $gsSubjectDirectory/surf
		    set fn [file join $gsaDefaultLocation($iType) \
				$gsSurfaceHemi($tkm_tSurfaceType(main)).orig]
		    if { [file exists $fn] } {
			set gsaDefaultLocation($iType) $fn
		    }
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}
	    }
	    LoadMainSurface - LoadMainAuxSurface - 
	    LoadOriginalAuxSurface - LoadPialAuxSurface -
	    WriteSurfaceValues { 
		if { $gsSubjectDirectory != "/" } {
		    set gsaDefaultLocation($iType) $gsSubjectDirectory/surf
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}
	    }
	    LoadHeadPts_Points { 
		set gsaDefaultLocation($iType) [exec pwd] 
	    }
	    LoadHeadPts_Transform { 
		set gsaDefaultLocation($iType) [exec pwd] 
	    }
	    Segmentation_ColorTable { 
		if { $gsSegmentationColorTable != "" } {
		    set gsaDefaultLocation($iType) $gsSegmentationColorTable
		} elseif { [info exists env(FREESURFER_HOME)] } {
		    set gsaDefaultLocation($iType) $env(FREESURFER_HOME)/
		} else {
		    set gsaDefaultLocation($iType) $gsSubjectDirectory 
		}
	    }
	    LoadFunctional-overlay - LoadFunctional-timecourse -
	    SpecifyRegistration-overlay - SpecifyRegistration-timecourse {
		set gsaDefaultLocation($iType) [exec pwd]
	    }
	    LoadGCA_Volume - SaveGCA {
		if { [info exists env(FREESURFER_HOME)] } {
		    set gsaDefaultLocation($iType) $env(FREESURFER_HOME)/average
		} else {
		    set gsaDefaultLocation($iType) $gsSubjectDirectory
		}
	    }
	    LoadGCA_Transform {
		if { $gsSubjectDirectory != "/" } {
		    set gsaDefaultLocation($iType) \
			$gsSubjectDirectory/mri/transforms
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}
		 }
	    SaveTIFF {
		if { $gsSubjectDirectory != "/" &&
		     [file isdirectory $gsSubjectDirectory/tiff] } {
		    set gsaDefaultLocation($iType) $gsSubjectDirectory/tiff
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}
	    }
	    default { 
		if { $gsSubjectDirectory != "/" } {
		    set gsaDefaultLocation($iType) $gsSubjectDirectory 
		} else {
		    set gsaDefaultLocation($iType) [exec pwd]
		}
	    }
	}
    }

    # If the file or path doesn't exist, just give pwd instead.
    if { ![file exists $gsaDefaultLocation($iType)] } {
	set gsaDefaultLocation($iType) [exec pwd]
    }

    # If location is a directory, make sure that the last char is a
    # slash.
    if { [file isdirectory $gsaDefaultLocation($iType)] } {
	if { [string range $gsaDefaultLocation($iType) end end] != "/" } {
	    set gsaDefaultLocation($iType) $gsaDefaultLocation($iType)/
	}
    }

    return $gsaDefaultLocation($iType)
}
proc SetDefaultLocation { iType isValue } {
    global gsaDefaultLocation
    if { [string range $isValue 0 0] == "/" } {
	set gsaDefaultLocation($iType) $isValue
    }
}
set tDlogSpecs(LoadVolume) [list \
  -title "Load Volume" \
  -prompt1 "Load Volume:" \
  -note1 "The volume file (or COR-.info for COR volumes)" \
  -entry1 [list GetDefaultLocation LoadVolume] \
  -default1 [list GetDefaultLocation LoadVolume] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadVolume %s1; SetDefaultLocation LoadVolume %s1} ]
set tDlogSpecs(LoadAuxVolume) [list \
  -title "Load Aux Volume" \
  -prompt1 "Load Volume:" \
  -note1 "The volume file (or COR-.info for COR volumes)" \
  -entry1 [list GetDefaultLocation LoadAuxVolume] \
  -default1 [list GetDefaultLocation LoadAuxVolume] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadAuxVolume %s1; SetDefaultLocation LoadAuxVolume %s1} ]
set tDlogSpecs(LoadGCA) [list \
  -title "Load GCA" \
  -prompt1 "Load Classifier Array:" \
  -note1 "The GCA file (*.gca)" \
  -entry1 [list GetDefaultLocation LoadGCA_Volume] \
  -default1 [list GetDefaultLocation LoadGCA_Volume] \
  -presets1 $glShortcutDirs \
  -prompt2 "Load Transform:" \
  -note2 "The file containing the transform to the atlas space" \
  -entry2 [list GetDefaultLocation LoadGCA_Transform] \
  -default2 [list GetDefaultLocation LoadGCA_Transform] \
  -presets2 $glShortcutDirs \
  -okCmd {LoadGCA %s1 %s2; \
  SetDefaultLocation LoadGCA_Volume %s1; \
  SetDefaultLocation LoadGCA_Transform %s2} ]
set tDlogSpecs(SaveGCA) [list \
  -title "Save GCA" \
  -prompt1 "Save Classifier Array:" \
  -note1 "The GCA file (*.gca)" \
  -entry1 [list GetDefaultLocation SaveGCA] \
  -default1 [list GetDefaultLocation SaveGCA] \
  -presets1 $glShortcutDirs \
  -okCmd {SaveGCA %s1; SetDefaultLocation SaveGCA %s1} ]
set tDlogSpecs(SaveVolumeAs) [list \
  -title "Save Main Volume As" \
  -prompt1 "Save COR Volume:" \
  -type1 dir \
  -note1 "The file name to save, or the directory in which to write the COR volume files" \
  -entry1 [list GetDefaultLocation SaveVolumeAs] \
  -default1 [list GetDefaultLocation SaveVolumeAs] \
  -presets1 $glShortcutDirs \
  -okCmd {SaveVolumeAs 0 %s1; SetDefaultLocation SaveVolumeAs %s1} ]
set tDlogSpecs(SaveAuxVolumeAs) [list \
  -title "Save Aux Volume As" \
  -prompt1 "Save COR Volume:" \
  -type1 dir \
  -note1 "The file name to save, or the directory in which to write the COR volume files" \
  -entry1 [list GetDefaultLocation SaveVolumeAs] \
  -default1 [list GetDefaultLocation SaveVolumeAs] \
  -presets1 $glShortcutDirs \
  -okCmd {SaveVolumeAs 1 %s1; SetDefaultLocation SaveVolumeAs %s1} ]
set tDlogSpecs(LoadVolumeDisplayTransform) [list \
  -title "Load Transform" \
  -prompt1 "Load Transform File:" \
  -note1 "The .lta or .xfm file containing the transform to load" \
  -entry1 [list GetDefaultLocation LoadVolumeDisplayTransform] \
  -default1 [list GetDefaultLocation LoadVolumeDisplayTransform] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadVolumeDisplayTransform 0 %s1; \
  SetDefaultLocation  LoadVolumeDisplayTransform %s1} ]
set tDlogSpecs(LoadAuxVolumeDisplayTransform) [list \
  -title "Load Aux Transform" \
  -prompt1 "Load Transform File:" \
  -note1 "The .lta or .xfm file containing the transform to load" \
  -entry1 [list GetDefaultLocation LoadAuxVolumeDisplayTransform] \
  -default1 [list GetDefaultLocation LoadAuxVolumeDisplayTransform] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadVolumeDisplayTransform 1 %s1; \
  SetDefaultLocation LoadAuxVolumeDisplayTransform %s1} ]
set tDlogSpecs(SaveLabelAs) [list \
  -title "Save Label As" \
  -prompt1 "Save Label:" \
  -note1 "The file name of the label to save" \
  -entry1 [list GetDefaultLocation SaveLabelAs] \
  -default1 [list GetDefaultLocation SaveLabelAs] \
  -presets1 $glShortcutDirs \
  -okCmd {SaveLabel %s1; SetDefaultLocation SaveLabelAs %s1} ]
set tDlogSpecs(LoadLabel) [list \
  -title "Load Label" \
  -prompt1 "Load Label:" \
  -note1 "The file name of the label to load" \
  -entry1 [list GetDefaultLocation LoadLabel] \
  -default1 [list GetDefaultLocation LoadLabel] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadLabel %s1; SetDefaultLocation LoadLabel %s1; RedrawScreen} ]
set tDlogSpecs(LoadMainSurface) [list \
  -title "Load Main Surface" \
  -prompt1 "Load Surface:" \
  -note1 "The file name of the surface to load" \
  -entry1 [list GetDefaultLocation LoadMainSurface] \
  -default1 [list GetDefaultLocation LoadMainSurface] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadMainSurface %s1; SetDefaultLocation LoadMainSurface %s1} ]
set tDlogSpecs(LoadOriginalSurface) [list \
  -title "Load Original Surface" \
  -prompt1 "Load Surface:" \
  -note1 "The file name of the surface to load" \
  -entry1 [list GetDefaultLocation LoadOriginalSurface] \
  -default1 [list GetDefaultLocation LoadOriginalSurface] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadOriginalSurface %s1; \
  SetDefaultLocation LoadOriginalSurface %s1} ]
set tDlogSpecs(LoadPialSurface) [list \
  -title "Load Pial Surface" \
  -prompt1 "Load Surface:" \
  -note1 "The file name of the surface to load" \
  -entry1 [list GetDefaultLocation LoadPialSurface] \
  -default1 [list GetDefaultLocation LoadPialSurface] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadPialSurface %s1; \
  SetDefaultLocation LoadPialSurface %s1} ]
set tDlogSpecs(LoadSurfaceAnnotation) [list \
  -title "Load Annotation" \
  -prompt1 "Load Annoation:" \
  -note1 "The file name of the annotation to load (*.annot)" \
  -entry1 [list GetDefaultLocation LoadSurfaceAnnotation] \
  -default1 [list GetDefaultLocation LoadSurfaceAnnotation] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadSurfaceAnnotation 0 %s1; \
  SetDefaultLocation LoadSurfaceAnnotation %s1} ]
set tDlogSpecs(LoadMainAuxSurface) [list \
  -title "Load Aux Main Surface" \
  -prompt1 "Load Surface:" \
  -note1 "The file name of the surface to load" \
  -entry1 [list GetDefaultLocation LoadMainAuxSurface] \
  -default1 [list GetDefaultLocation LoadMainAuxSurface] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadMainSurface 1 %s1; SetDefaultLocation LoadMainAuxSurface %s1} ]
set tDlogSpecs(LoadOriginalAuxSurface) [list \
  -title "Load Aux Original Surface" \
  -prompt1 "Load Surface:" \
  -note1 "The file name of the surface to load" \
  -entry1 [list GetDefaultLocation LoadOriginalAuxSurface] \
  -default1 [list GetDefaultLocation LoadOriginalAuxSurface] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadOriginalSurface 1 %s1; \
  SetDefaultLocation LoadOriginalAuxSurface %s1} ]
set tDlogSpecs(LoadPialAuxSurface) [list \
  -title "Load Aux Pial Surface" \
  -prompt1 "Load Surface:" \
  -note1 "The file name of the surface to load" \
  -entry1 [list GetDefaultLocation LoadPialAuxSurface] \
  -default1 [list GetDefaultLocation LoadPialAuxSurface] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadPialSurface 1 %s1; \
  SetDefaultLocation LoadPialAuxSurface %s1} ]
set tDlogSpecs(LoadAuxSurfaceAnnotation) [list \
  -title "Load Annotation" \
  -prompt1 "Load Annoation:" \
  -note1 "The file name of the annotation to load (*.annot)" \
  -entry1 [list GetDefaultLocation LoadAuxSurfaceAnnotation] \
  -default1 [list GetDefaultLocation LoadAuxSurfaceAnnotation] \
  -presets1 $glShortcutDirs \
  -okCmd {LoadSurfaceAnnotation 1 %s1; \
  SetDefaultLocation LoadAuxSurfaceAnnotation %s1} ]
set tDlogSpecs(WriteSurfaceValues) [list \
  -title "Write Surface Values" \
  -prompt1 "Save Values As:" \
  -note1 "The file name of the values file to write" \
  -entry1 [list GetDefaultLocation WriteSurfaceValues] \
  -default1 [list GetDefaultLocation WriteSurfaceValues] \
  -presets1 $glShortcutDirs \
  -okCmd {WriteSurfaceValues %s1; \
  SetDefaultLocation WriteSurfaceValues %s1} ]
set tDlogSpecs(ImportSurfaceAnnotation) [list \
  -title "Import Surface Annotation" \
  -prompt1 "Load Annotation File:" \
  -note1 "A .annot file containing the annotation data" \
  -entry1 [list GetDefaultLocation ImportSurfaceAnnotation] \
  -default1 [list GetDefaultLocation ImportSurfaceAnnotation] \
  -presets1 $glShortcutDirs \
  -prompt2 "Load Color Table:" \
  -note2 "The file containing the colors and ROI definitions" \
  -entry2 [list GetDefaultLocation Segmentation_ColorTable] \
  -default2 [list GetDefaultLocation Segmentation_ColorTable] \
  -presets2 $glShortcutDirs \
  -okCmd {ImportSurfaceAnnotationToSegmentation 0 %s1 %s2; \
  SetDefaultLocation ImportSegmentation_Volume %s1; \
  SetDefaultLocation Segmentation_ColorTable %s2} ]
set tDlogSpecs(PrintTimeCourse) [list \
  -title "Print Time Course" \
  -prompt1 "Save Summary As:" \
  -note1 "The file name of the text summary to create" \
  -entry1 [list GetDefaultLocation PrintTimeCourse] \
  -default1 [list GetDefaultLocation PrintTimeCourse] \
  -presets1 $glShortcutDirs \
  -okCmd {TimeCourse_PrintSelectionRangeToFile %s1; \
  SetDefaultLocation PrintTimeCourse %s1} ]
set tDlogSpecs(SaveTimeCourseToPS) [list \
  -title "Save Time Course" \
  -prompt1 "Save Time Course As:" \
  -note1 "The file name of the PostScript file to create" \
  -entry1 [list GetDefaultLocation SaveTimeCourseToPS] \
  -default1 [list GetDefaultLocation SaveTimeCourseToPS] \
  -presets1 $glShortcutDirs \
  -okCmd {TimeCourse_SaveGraphToPS %s1; \
  SetDefaultLocation SaveTimeCourseToPS %s1} ]
set tDlogSpecs(NewSegmentation) [list \
  -title "New Segmentation" \
  -prompt1 "Load Color Table:" \
  -note1 "The file containing the colors and ROI definitions" \
  -entry1 [list GetDefaultLocation Segmentation_ColorTable] \
  -default1 [list GetDefaultLocation Segmentation_ColorTable] \
  -presets1 $glShortcutDirs \
  -okCmd {NewSegmentationVolume 0 0 %s1; \
  SetDefaultLocation Segmentation_ColorTable %s1} ]
set tDlogSpecs(NewAuxSegmentation) [list \
  -title "New Aux Segmentation" \
  -prompt1 "Load Color Table:" \
  -note1 "The file containing the colors and ROI definitions" \
  -entry1 [list GetDefaultLocation Segmentation_ColorTable] \
  -default1 [list GetDefaultLocation Segmentation_ColorTable] \
  -presets1 $glShortcutDirs \
  -okCmd {NewSegmentationVolume 1 1 %s1; \
  SetDefaultLocation Segmentation_ColorTable %s1} ]
set tDlogSpecs(LoadSegmentation) [list \
  -title "Load Segmentation" \
  -prompt1 "Load Volume:" \
  -note1 "The volume file (or COR-.info for COR volumes)" \
  -entry1 [list GetDefaultLocation LoadSegmentation] \
  -default1 [list GetDefaultLocation LoadSegmentation] \
  -presets1 $glShortcutDirs \
  -prompt2 "Load Color Table:" \
  -note2 "The file containing the colors and ROI definitions" \
  -entry2 [list GetDefaultLocation Segmentation_ColorTable] \
  -default2 [list GetDefaultLocation Segmentation_ColorTable] \
  -presets2 $glShortcutDirs \
  -okCmd {LoadSegmentationVolume 0 %s1 %s2; \
  SetDefaultLocation ImportSegmentation_Volume %s1; \
  SetDefaultLocation Segmentation_ColorTable %s2} ]
set tDlogSpecs(LoadAuxSegmentation) [list \
  -title "Load Aux Segmentation" \
  -prompt1 "Load Volume:" \
  -note1 "The volume file (or COR-.info for COR volumes)" \
  -entry1 [list GetDefaultLocation LoadAuxSegmentation] \
  -default1 [list GetDefaultLocation LoadAuxSegmentation] \
  -presets1 $glShortcutDirs \
  -prompt2 "Load Color Table:" \
  -note2 "The file containing the colors and ROI definitions" \
  -entry2 [list GetDefaultLocation Segmentation_ColorTable] \
  -default2 [list GetDefaultLocation Segmentation_ColorTable] \
  -presets2 $glShortcutDirs \
  -okCmd {LoadSegmentationVolume 1 %s1 %s2; \
  SetDefaultLocation ImportSegmentation_Volume %s1; \
  SetDefaultLocation Segmentation_ColorTable %s2} ]
set tDlogSpecs(SaveSegmentationAs) [list \
  -title "Save Segmentation As" \
  -prompt1 "Save COR Volume:" \
  -type1 dir \
  -note1 "The directory in which to write the COR volume files" \
  -entry1 [list GetDefaultLocation SaveSegmentationAs] \
  -default1 [list GetDefaultLocation SaveSegmentationAs] \
  -presets1 $glShortcutDirs \
  -okCmd {SaveSegmentationVolume 0 %s1; \
  SetDefaultLocation SaveSegmentationAs %s1} ]
set tDlogSpecs(SaveAuxSegmentationAs) [list \
  -title "Save Aux Segmentation As" \
  -prompt1 "Save COR Volume:" \
  -type1 dir \
  -note1 "The directory in which to write the COR volume files" \
  -entry1 [list GetDefaultLocation SaveSegmentationAs] \
  -default1 [list GetDefaultLocation SaveSegmentationAs] \
  -presets1 $glShortcutDirs \
  -okCmd {SaveSegmentationVolume 1 %s1; \
  SetDefaultLocation SaveSegmentationAs %s1} ]
set tDlogSpecs(ExportChangedSegmentationVolume) [list \
  -title "Save Changed Segmentation Values As" \
  -prompt1 "Save COR Volume:" \
  -type1 dir \
  -note1 "The directory in which to write the COR volume files" \
  -entry1 [list GetDefaultLocation SaveSegmentationAs] \
  -default1 [list GetDefaultLocation SaveSegmentationAs] \
  -presets1 $glShortcutDirs \
  -okCmd {ExportChangedSegmentationVolume 0 %s1; \
  SetDefaultLocation ExportChangedSegmentationVolume %s1} ]
set tDlogSpecs(ExportAuxChangedSegmentationVolume) [list \
  -title "Save Aux Changed Segmentation Values As" \
  -prompt1 "Save COR Volume:" \
  -type1 dir \
  -note1 "The directory in which to write the COR volume files" \
  -entry1 [list GetDefaultLocation SaveSegmentationAs] \
  -default1 [list GetDefaultLocation SaveSegmentationAs] \
  -presets1 $glShortcutDirs \
  -okCmd {ExportChangedSegmentationVolume 1 %s1; \
  SetDefaultLocation ExportChangedSegmentationVolume %s1} ]
set tDlogSpecs(LoadHeadPts) [list \
  -title "Load Head Points" \
  -prompt1 "Load Head Points:" \
  -note1 "The file name of the .hpts head points file" \
  -entry1 [list GetDefaultLocation LoadHeadPts_Points] \
  -default1 [list GetDefaultLocation LoadHeadPts_Points] \
  -presets1 $glShortcutDirs \
  -prompt2 "Load Transform File:" \
  -note2 "The file name of the .trans transform file" \
  -entry2 [list GetDefaultLocation LoadHeadPts_Transform] \
  -default2 [list GetDefaultLocation LoadHeadPts_Transform] \
  -presets2 $glShortcutDirs \
  -okCmd {LoadHeadPts %s1 %s2; \
  SetDefaultLocation LoadHeadPts_Points %s1; \
  SetDefaultLocation LoadHeadPts_Transform %s2} ]
set tDlogSpecs(SaveRGB) [list \
  -title "Save RGB" \
  -prompt1 "Save RGB File:" \
  -note1 "The file name of the RGB picture to create" \
  -entry1 [list GetDefaultLocation SaveRGB] \
  -default1 [list GetDefaultLocation SaveRGB] \
  -presets1 $glShortcutDirs \
  -okCmd {SaveRGB %s1; SetDefaultLocation SaveRGB %s1} ]
set tDlogSpecs(SaveTIFF) [list \
  -title "Save TIFF" \
  -prompt1 "Save TIFF File:" \
  -note1 "The file name of the TIFF picture to create" \
  -entry1 [list GetDefaultLocation SaveTIFF] \
  -default1 [list GetDefaultLocation SaveTIFF] \
  -presets1 $glShortcutDirs \
  -okCmd {SaveTIFF %s1; SetDefaultLocation SaveTIFF %s1} ]
set tDlogSpecs(WriteLineLabel) [list \
  -title "Write Line Label" \
  -prompt1 "Save Label As:" \
  -note1 "The file name of the label file to write" \
  -entry1 [list GetDefaultLocation WriteLineLabel] \
  -default1 [list GetDefaultLocation WriteLineLabel] \
  -presets1 $glShortcutDirs \
  -okCmd {WriteLineToLabel %s1; \
  SetDefaultLocation WriteLineLabel %s1} ]
set tDlogSpecs(SaveGDFPlotToPS) [list \
  -title "Save Group Plot" \
  -prompt1 "Save Plot As:" \
  -note1 "The file name of the PostScript file to create" \
  -default1 [list GetDefaultLocation SaveGDFPlotToPS] \
  -entry1 [list GetDefaultLocation SaveGDFPlotToPS] \
  -presets1 $glShortcutDirs \
  -okCmd {
      SetDefaultLocation SaveGDFPlotToPS %s1;
      FsgdfPlot_SaveToPostscript [GetGDFID] %s1  
  }]
set tDlogSpecs(SaveGDFPlotToTable) [list \
   -title "Save Group Data" \
   -prompt1 "Save Plot As:" \
   -note1 "The file name of the table to create" \
   -default1 [list GetDefaultLocation SaveGDFPlotToTable] \
   -entry1 [list GetDefaultLocation SaveGDFPlotToTable] \
   -presets1 $glShortcutDirs \
   -okCmd {
       SetDefaultLocation SaveGDFPlotToTable %s1;
       FsgdfPlot_SaveToTable [GetGDFID] %s1  
   }]

proc DoFileDlog { which } {
    global tDlogSpecs
    tkm_DoFileDlog $tDlogSpecs($which)
}

proc FindVertex { inVertex } {

    global Surf_tVertexSet
    global gFindingSurface
    
    if { $Surf_tVertexSet(main)   == $gFindingSurface } {
	GotoMainVertex $inVertex
    }
    if { $Surf_tVertexSet(original)  == $gFindingSurface } {
	GotoOriginalVertex $inVertex
    }
    if { $Surf_tVertexSet(pial) == $gFindingSurface } {
	GotoPialVertex $inVertex
    }
}

proc DoLoadFunctionalDlog { isType } {

    global gDialog gaLinkedVars
    global gaScalarValueID gsaLabelContents
    global glShortcutDirs
    global FunD_tRegistration
    global gfnFunctional gsFuncLoadType 
    global gRegistrationType gfnFunctionalRegistration

    set wwDialog .wwLoadFunctionalDlog

    set knWidth 400

    set gsFuncLoadType $isType

    set sTitle ""
    set sPrompt ""
    if { $gsFuncLoadType == "overlay" } {
	set sTitle "Load Overlay"
	set sPrompt "Load Overlay:"
    } elseif { $gsFuncLoadType == "timecourse" } {
	set sTitle "Load Time Course"
	set sPrompt "Load Time Course:"
    }
    
    # try to create the dlog...
    if { [Dialog_Create $wwDialog $sTitle {-borderwidth 10}] } {
	
	set fwFile             $wwDialog.fwFile
	set fwFileNote         $wwDialog.fwFileNote
	set fwRegistration     $wwDialog.fwRegistration
	set fwButtons          $wwDialog.fwButtons
	
	set gfnFunctional [GetDefaultLocation LoadFunctional-$gsFuncLoadType]
	tkm_MakeFileSelector $fwFile $sPrompt gfnFunctional \
	    [list GetDefaultLocation LoadFunctional-$gsFuncLoadType] \
	    $glShortcutDirs

	[$fwFile.ew subwidget entry] icursor end
	
	tkm_MakeSmallLabel $fwFileNote \
	    "The volume file (or COR-.info for COR volumes)" 400

	tixLabelFrame $fwRegistration \
	    -label "Registration" \
	    -labelside acrosstop \
	    -options { label.padX 5 }
	
	set fwRegSub           [$fwRegistration subwidget frame]

	set fwRegNone          $fwRegSub.fwRegNone
	set fwRegFile          $fwRegSub.fwRegFile
	set fwRegFileName      $fwRegSub.fwRegFileName
	set fwRegFind          $fwRegSub.fwRegFind
	set fwRegIdentity      $fwRegSub.fwRegIdentity

	# The bit of code in the radio buttons disables the file entry
	# field when the file radio button is not clicked.
	tkm_MakeRadioButton $fwRegNone "No registration needed" \
	    gRegistrationType $FunD_tRegistration(noneNeeded) "set state disabled; if { \[set gRegistrationType\] == $FunD_tRegistration(file)} { set state normal }; $fwRegFileName.ew config -state \$state; $fwRegFileName.bw config -state \$state"
	tkm_MakeRadioButton $fwRegFile "Specify registration file" \
	    gRegistrationType $FunD_tRegistration(file) "set state disabled; if { \[set gRegistrationType\] == $FunD_tRegistration(file) } { set state normal }; $fwRegFileName.ew config -state \$state; $fwRegFileName.bw config -state \$state"
	tkm_MakeFileSelector $fwRegFileName "register.dat file:" \
	    gfnFunctionalRegistration \
	    [list GetDefaultLocation LoadFunctional-$gsFuncLoadType] \
	    $glShortcutDirs
	tkm_MakeRadioButton $fwRegFind "Find registration in data directory" \
	    gRegistrationType $FunD_tRegistration(find) "set state disabled; if { \[set gRegistrationType\] == $FunD_tRegistration(file)} { set state normal }; $fwRegFileName.ew config -state \$state; $fwRegFileName.bw config -state \$state"
	tkm_MakeRadioButton $fwRegIdentity "Calculate identity matrix" \
	    gRegistrationType $FunD_tRegistration(identity) "set state disabled; if { \[set gRegistrationType\] == $FunD_tRegistration(file)} { set state normal }; $fwRegFileName.ew config -state \$state; $fwRegFileName.bw config -state \$state"
	set gRegistrationType $FunD_tRegistration(file)

	pack $fwRegNone $fwRegFile $fwRegFileName $fwRegFind $fwRegIdentity \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5

	# buttons.
        tkm_MakeCancelOKButtons $fwButtons $wwDialog \
	    {set fnFunctional $gfnFunctional; 
	      SetDefaultLocation LoadFunctional-$gsFuncLoadType $gfnFunctional;
		DoLoadFunctional $gsFuncLoadType $gfnFunctional $gRegistrationType $gfnFunctionalRegistration }
	
	pack $fwFile $fwFileNote $fwRegistration $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5
	
	# after the next idle, the window will be mapped. set the min
	# width to our width and the min height to the mapped height.
	after idle [format {
	    update idletasks
	    wm minsize %s %d [winfo reqheight %s]
	    wm geometry %s =%dx[winfo reqheight %s]
	} $wwDialog $knWidth $wwDialog $wwDialog $knWidth $wwDialog] 
    }
}

proc DoLoadFunctional { isType ifnVolume iRegistrationType ifnRegister } {

    if { [string match $isType overlay] } {
	LoadFunctionalOverlay \
	    $ifnVolume $iRegistrationType $ifnRegister
    } elseif { [string match $isType timecourse] } {
	LoadFunctionalTimeCourse \
	    $ifnVolume $iRegistrationType $ifnRegister
    }
}

proc DoLoadDTIDlog {} {
    global gDialog glShortcutDirs
    global nColorX nColorY nColorZ

    set wwDialog .wwLoadDTIDlog

    set knWidth 400

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load DTI Volume" {-borderwidth 10}] } {

	set fwEVFile           $wwDialog.fwEVFile
	set fwEVFileNote       $wwDialog.fwEVFileNote
	set fwFAFile           $wwDialog.fwFAFile
	set fwFAFileNote       $wwDialog.fwFAFileNote
	set lwColor            $wwDialog.lwColor
	set fwColorTable       $wwDialog.fwColorTable
	set fwButtons          $wwDialog.fwButtons
	
	set sEVFileName ""
	tkm_MakeFileSelector $fwEVFile "Load DTI Vector Volume:" sEVFileName \
	    [list GetDefaultLocation LoadDTIVolume] \
	    $glShortcutDirs

	[$fwEVFile.ew subwidget entry] icursor end
	
	tkm_MakeSmallLabel $fwEVFileNote "The DTI vector volume to load" 400
	
	set sFAFileName ""
	tkm_MakeFileSelector $fwFAFile "Load DTI FA Volume:" sFAFileName \
	    [list GetDefaultLocation LoadDTIVolume] \
	    $glShortcutDirs

	[$fwFAFile.ew subwidget entry] icursor end
	
	tkm_MakeSmallLabel $fwFAFileNote "The DTI FA volume to load" 400
	
	tkm_MakeNormalLabel $lwColor "Color orientation:"

	frame $fwColorTable

	set nRow 1
	foreach color {Red Green Blue} {
	    tkm_MakeSmallLabel $fwColorTable.lw$color "$color"
	    grid $fwColorTable.lw$color -column 0 -row $nRow
	    incr nRow
	}
	set nColumn 1
	foreach axis {X Y Z} {
	    tkm_MakeSmallLabel $fwColorTable.lw$axis "$axis"
	    grid $fwColorTable.lw$axis -column $nColumn -row 0

	    set nRow 1
	    foreach nColor {0 1 2} {
		radiobutton $fwColorTable.rb$axis-$nColor \
		    -variable nColor$axis -value $nColor
		grid $fwColorTable.rb$axis-$nColor \
		    -column $nColumn -row $nRow
		incr nRow
	    }
	    incr nColumn
	}
        set nColorX 0
        set nColorY 1
        set nColorZ 2
	

	# ok and cancel buttons.
	tkm_MakeCancelOKButtons $fwButtons $wwDialog \
	    { LoadDTIVolume $sEVFileName $sFAFileName $nColorX $nColorY $nColorZ;
	    SetDefaultLocation LoadDTIVolume $sEVFileName }
	
	pack $fwEVFile $fwEVFileNote $fwFAFile $fwFAFileNote \
	    $lwColor $fwColorTable $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5

	# after the next idle, the window will be mapped. set the min
	# width to our width and the min height to the mapped height.
	after idle [format {
	    update idletasks
	    wm minsize %s %d [winfo reqheight %s]
	    wm geometry %s =%dx[winfo reqheight %s]
	} $wwDialog $knWidth $wwDialog $wwDialog $knWidth $wwDialog] 
    }
}

proc DoLoadGDF {} {

    global gDialog gaLinkedVars
    global glShortcutDirs
    global gfnGDF
    global gGDFRegType gfnGDFReg
    global FunD_tRegistration

    set wwDialog .wwLoadGDFDlog

    set knWidth 400

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load GDF" {-borderwidth 10}] } {
	
	set fwFile             $wwDialog.fwFile
	set fwFileNote         $wwDialog.fwFileNote
	set fwRegistration     $wwDialog.fwRegistration
	set fwButtons          $wwDialog.fwButtons
	
	set gfnGDF [GetDefaultLocation LoadGDF]
	tkm_MakeFileSelector $fwFile "Load GDF: " gfnGDF \
	    [list GetDefaultLocation LoadGDF] \
	    $glShortcutDirs

	[$fwFile.ew subwidget entry] icursor end
	
	tkm_MakeSmallLabel $fwFileNote \
	    "The .fsgd file" 400

	tixLabelFrame $fwRegistration \
	    -label "Registration" \
	    -labelside acrosstop \
	    -options { label.padX 5 }
	
	set fwRegSub           [$fwRegistration subwidget frame]

	set fwRegFile          $fwRegSub.fwRegFile
	set fwRegFileName      $fwRegSub.fwRegFileName
	set fwRegFind          $fwRegSub.fwRegFind
	set fwRegIdentity      $fwRegSub.fwRegIdentity

	# The bit of code in the radio buttons disables the file entry
	# field when the file radio button is not clicked.
	tkm_MakeRadioButton $fwRegFile "Specify registration file" \
	    gGDFRegType $FunD_tRegistration(file) "set state disabled; if { \[set gGDFRegType\] == $FunD_tRegistration(file) } { set state normal }; $fwRegFileName.ew config -state \$state; $fwRegFileName.bw config -state \$state"
	tkm_MakeFileSelector $fwRegFileName "register.dat file:" \
	    gfnGDFReg \
	    [list GetDefaultLocation LoadGDF] \
	    $glShortcutDirs
	tkm_MakeRadioButton $fwRegFind "Find registration in data directory" \
	    gGDFRegType $FunD_tRegistration(find) "set state disabled; if { \[set gGDFRegType\] == $FunD_tRegistration(file)} { set state normal }; $fwRegFileName.ew config -state \$state; $fwRegFileName.bw config -state \$state"
	tkm_MakeRadioButton $fwRegIdentity "Calculate identity matrix" \
	    gGDFRegType $FunD_tRegistration(identity) "set state disabled; if { \[set gGDFRegType\] == $FunD_tRegistration(file)} { set state normal }; $fwRegFileName.ew config -state \$state; $fwRegFileName.bw config -state \$state"
	set gGDFRegType $FunD_tRegistration(file)

	pack $fwRegFile $fwRegFileName $fwRegFind $fwRegIdentity \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5

	# buttons.
        tkm_MakeCancelOKButtons $fwButtons $wwDialog \
	    {set fnGDF $gfnGDF; 
	      SetDefaultLocation LoadGDF $gfnGDF;
		LoadGDF $gfnGDF $gGDFRegType $gfnGDFReg }
	
	pack $fwFile $fwFileNote $fwRegistration $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5
	
	# after the next idle, the window will be mapped. set the min
	# width to our width and the min height to the mapped height.
	after idle [format {
	    update idletasks
	    wm minsize %s %d [winfo reqheight %s]
	    wm geometry %s =%dx[winfo reqheight %s]
	} $wwDialog $knWidth $wwDialog $wwDialog $knWidth $wwDialog] 
    }
}

proc DoSaveDlog {} {

    global gDialog

    set wwDialog .wwSaveDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save Main Volume" {-borderwidth 10}] } {

	set fwMain    $wwDialog.fwMain
	set fwButtons $wwDialog.fwButtons
	
	# prompt
	tkm_MakeNormalLabel $fwMain\
	    "Are you sure you wish to save changes to the main volume?"
	
	# ok and cancel buttons.
	tkm_MakeCancelOKButtons $fwButtons $wwDialog \
	    { SaveVolume 0 }
	
	pack $fwMain $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5
    }
}

proc DoAuxSaveDlog {} {

    global gDialog

    set wwDialog .wwAuxSaveDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save Aux Volume" {-borderwidth 10}] } {

	set fwMain    $wwDialog.fwMain
	set fwButtons $wwDialog.fwButtons
	
	# prompt
	tkm_MakeNormalLabel $fwMain\
	    "Are you sure you wish to save changes to the aux volume?"
	
	# ok and cancel buttons.
	tkm_MakeCancelOKButtons $fwButtons $wwDialog \
	    { SaveVolume 1 }
	
	pack $fwMain $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5
    }
}

proc DoAskSaveChangesDlog {} {

  global gDialog

	set wwDialog .wwSaveDlog

	# try to create the dlog...
	if { [Dialog_Create $wwDialog "Save Changes" {-borderwidth 10}] } {

		set fwMain    $wwDialog.fwMain
		set fwButtons $wwDialog.fwButtons
	
		# prompt
		tkm_MakeNormalLabel $fwMain \
		"Do you wish to save changes to the main volume?"
	
		# ok and cancel buttons.
		tkm_MakeButtons $fwButtons {
			{text "Save and Quit" {SaveVolume 0; QuitMedit} }
			{text "Just Quit" { QuitMedit } }
			{text "Don't Quit" { Dialog_Close .wwSaveDlog} } 
		}

		pack $fwMain $fwButtons \
		-side top       \
		-expand yes     \
		-fill x         \
		-padx 5         \
		-pady 5
	}
}	

proc DoAskQuitDlog {} {

  global gDialog

	set wwDialog .wwQuitDlog

	# try to create the dlog...
	if { [Dialog_Create $wwDialog "Quit" {-borderwidth 10}] } {

		set fwMain    $wwDialog.fwMain
		set fwButtons $wwDialog.fwButtons
	
		# prompt
		tkm_MakeNormalLabel $fwMain "Do you wish to quit?"
	
		# ok and cancel buttons.
		tkm_MakeButtons $fwButtons {
			{text "Quit" {QuitMedit} }
			{text "Cancel" { Dialog_Close .wwQuitDlog} } 
		}

		pack $fwMain $fwButtons \
		-side top       \
		-expand yes     \
		-fill x         \
		-padx 5         \
		-pady 5
	}
}	

proc DoBrushInfoDlog {} {

    global gDialog
    global DspA_tBrushShape_Square DspA_tBrushShape_Circle
    global DspA_tBrushTarget
    global gBrushInfo

    set wwDialog .wwBrushInfoDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Brush Info" {-borderwidth 10}] } {
  
	set fwTop                $wwDialog.fwTop
	set fwShape              $fwTop.fwShape
	set fwButtons            $wwDialog.fwButtons
	
	frame $fwTop
	
	set fwRadiusScale        $fwTop.fwRadiusScale
	set fwTargetLabel        $fwTop.fwTargetLabel
	set fwMain               $fwTop.fwMain
	set fwAux                $fwTop.fwAux
	set fwMainAux            $fwTop.fwMainAux
	set fwShapeLabel         $fwTop.fwShapeLabel
	set fwCircle             $fwTop.fwCircle
	set fwSquare             $fwTop.fwSquare
	set fw3DCheckbox         $fwTop.fw3DCheckbox
	
 	# target radio buttons
	tkm_MakeNormalLabel $fwTargetLabel "Target"
	tkm_MakeRadioButton $fwMain "Main volume" \
	    gBrushInfo(target) $DspA_tBrushTarget(main) \
	    "SetBrushConfiguration"
	tkm_MakeRadioButton $fwAux "Aux volume" \
	    gBrushInfo(target) $DspA_tBrushTarget(aux) \
	    "SetBrushConfiguration"
	tkm_MakeRadioButton $fwMainAux "Main and aux volume" \
	    gBrushInfo(target) $DspA_tBrushTarget(mainaux) \
	    "SetBrushConfiguration"

	# radius
	tkm_MakeSliders $fwRadiusScale { \
	       { {"Radius"} gBrushInfo(radius) 1 20 100 "SetBrushConfiguration" 1 } }
	
 	# shape radio buttons
	tkm_MakeNormalLabel $fwShapeLabel "Shape"
	tkm_MakeRadioButton $fwCircle "Circle" \
	    gBrushInfo(shape) $DspA_tBrushShape_Circle \
	    "SetBrushConfiguration"
	tkm_MakeRadioButton $fwSquare "Square" \
	    gBrushInfo(shape) $DspA_tBrushShape_Square \
	    "SetBrushConfiguration"
	
	# 3d checkbox
	tkm_MakeCheckboxes $fw3DCheckbox y { \
	      { text "3D" gBrushInfo(3d) "SetBrushConfiguration" } }
	
	# pack them in a column
	pack $fwTargetLabel $fwMain $fwAux $fwMainAux \
	    $fwRadiusScale $fwShapeLabel $fwCircle \
	    $fwSquare $fw3DCheckbox             \
	    -side top                           \
	    -anchor w                           \
	    -expand yes                         \
	    -fill x
	
	# buttons. 
	tkm_MakeCloseButton $fwButtons $wwDialog
	
	pack $fwTop $fwButtons \
	    -side top       \
	    -expand yes     \
    -fill x
    }
}

proc DoEditBrushInfoDlog {} {

    global gDialog
    global ksaBrushString
    global DspA_tBrush_EditOne DspA_tBrush_EditTwo
    global DspA_tBrushMode
    global tkm_tVolumeType
    global gEditBrush
    global gVolume
    global gBrushInfo

    set wwDialog .wwEditBrushInfoDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Volume Brush Info" {-borderwidth 10}] } {
	
	set fwTop                $wwDialog.fwTop
	set fwInfo               $fwTop.fwInfo
	set fwFillInfo           $fwTop.fwFillInfo
	set fwButtons            $wwDialog.fwButtons
	
	frame $fwTop
	
	tixNoteBook $fwInfo
	foreach tool "$DspA_tBrush_EditOne $DspA_tBrush_EditTwo" {
	    
	    $fwInfo add pane$tool -label $ksaBrushString($tool)
	    
	    set fw [$fwInfo subwidget pane$tool]
	    set fwThresh         $fw.fwThresh
	    set fwMode           $fw.fwMode
	    set fwNewValue       $fw.fwNewValue
	    set fwCloneSrc       $fw.fwCloneSrc
	    set fwDefaults       $fw.fwDefaults
	    
	    # This is a hack. Even though the volume min/max of a COR volume
	    # might be 0-253, 255 is still a valid value, so we'll change
	    # the max here to allow that.
	    set min $gVolume(0,minValue)
	    set max $gVolume(0,maxValue)
	    if { $max < 255 } { set max 255 }
	    
	    # low, high, and new value sliders
	    tkm_MakeSliders $fwThresh \
		[list \
		     [list {"Low"} gEditBrush($tool,low) $min $max \
			  200 "SetEditBrushConfiguration" 1] \
		     [list {"High"} gEditBrush($tool,high) $min $max \
			  200 "SetEditBrushConfiguration" 1 ]]
	    
	    tkm_MakeRadioButtons $fwMode y "Mode" gEditBrush($tool,mode) \
		[list \
		     [list text "New Value" $DspA_tBrushMode(set) \
			  "SetBrushInfoModeParamsState set $fwCloneSrc $fwNewValue; SetEditBrushConfiguration"] \
		     [list text "Clone" $DspA_tBrushMode(clone) \
			  "SetBrushInfoModeParamsState clone $fwCloneSrc $fwNewValue; SetEditBrushConfiguration"]]
	    
	    tkm_MakeEntry $fwNewValue "New Value" gEditBrush($tool,new) \
		6 "SetEditBrushConfiguration"
	    
	    tkm_MakeRadioButtons $fwCloneSrc y "Clone Source" \
		gEditBrush($tool,cloneSource) \
		[list \
		     [list text "Main Volume" $tkm_tVolumeType(main) \
			  "SetEditBrushConfiguration"] \
		     [list text "Aux Volume"  $tkm_tVolumeType(aux) \
			  "SetEditBrushConfiguration"]]
	    
	    
	    
	    # defaults button
	    tkm_MakeButtons $fwDefaults \
		[list \
		     [list text "Restore Defaults" "SetBrushInfoToDefaults $tool"]]
	    
	    # pack them in a column
	    pack $fwThresh $fwMode $fwNewValue $fwCloneSrc $fwDefaults \
		-side top                           \
		-anchor w                           \
		-expand yes                         \
		-fill x
	}
	
	
	# fill characteristics
	tixLabelFrame $fwFillInfo \
	    -label "Fill Parameters" \
	    -labelside acrosstop \
	    -options { label.padX 5 }
	
	set fwFillInfoSub [$fwFillInfo subwidget frame]
	set fwFill      $fwFillInfoSub.fwFill
	set fw3D        $fwFill.fw3D
	set fwSliders   $fwFill.fwSliders
	set fwNote      $fwFill.fwNotes
	
	frame $fwFill

	# 3d
	tkm_MakeCheckboxes $fw3D y { 
	    { text "3D" gBrushInfo(3dfill) "SendFillAnatomicalInfo" } }

	# fuzziness and max distance
	tkm_MakeSliders $fwSliders { 
	    { "Fuzziness" gBrushInfo(fuzzy) 
		0 255 50 "SendFillAnatomicalInfo" 1 0.1 } 
	    {  "\"Max Distance\"" gBrushInfo(distance) 
		0 255 50 "SendFillAnatomicalInfo" 1 0.1 } }
	tkm_MakeSmallLabel $fwNote "enter 0 for no limit"
	
	pack $fw3D $fwSliders $fwNote $fwFill \
	    -side top \
	    -expand yes \
	    -fill x


	# buttons. 
	tkm_MakeCloseButton $fwButtons $wwDialog
	
	pack $fwTop $fwInfo $fwFillInfo $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x
    }
}

proc SetBrushInfoModeParamsState { iWhich cloneBase setBase } {

    set cloneState disabled
    set cloneTextColor gray
    set setState disabled
    set setTextColor gray
    if { $iWhich == "clone" } { 
	set cloneState normal 
	set cloneTextColor black
    }
    if { $iWhich == "set" } { 
	set setState normal 
	set setTextColor black
    }

    [$cloneBase subwidget label] config -fg $cloneTextColor
    [$cloneBase subwidget frame].rb0 config -state $cloneState
    [$cloneBase subwidget frame].rb1 config -state $cloneState
    [$cloneBase subwidget frame].lw0 config -state $cloneState
    [$cloneBase subwidget frame].lw1 config -state $cloneState

    $setBase.lwLabel config -fg $setTextColor
    $setBase.ewEntry config -state $setState
}

proc DoSegmentationVolumeDisplayInfoDlog { } {

    global gDialog
    global gfSegmentationVolumeAlpha

    set wwDialog .wwSegmentationVolumeDisplay

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Segmentation Display" {-borderwidth 10}] } {
  
  set fwSlider             $wwDialog.fwSlider
  set fwButtons            $wwDialog.fwButtons
  
  # alpha
  tkm_MakeSliders $fwSlider { \
    { {"Overlay Alpha"} gfSegmentationVolumeAlpha \
    0 1 80 "SetSegmentationVolumeConfiguration" 1 0.1 } }

  # buttons. 
  tkm_MakeCloseButton $fwButtons $wwDialog 

  pack $fwSlider $fwButtons \
    -side top       \
    -expand yes     \
    -fill x
   }

}

proc DoDTIVolumeDisplayInfoDlog { } {

    global gDialog
    global gfDTIVolumeAlpha

    set wwDialog .wwDTIVolumeDisplay

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "DTI Display" {-borderwidth 10}] } {
  
  set fwSlider             $wwDialog.fwSlider
  set fwButtons            $wwDialog.fwButtons
  
  # alpha
  tkm_MakeSliders $fwSlider { \
    { {"Overlay Alpha"} gfDTIVolumeAlpha \
    0 1 80 "SetDTIVolumeConfiguration" 1 0.1 } }

  # buttons. 
  tkm_MakeCloseButton $fwButtons $wwDialog 

  pack $fwSlider $fwButtons \
    -side top       \
    -expand yes     \
    -fill x
   }

}

proc DoRecomputeSegmentation {} {

    global gDialog

    set wwDialog .wwRecomputeSegmentation

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Recompute Segmentation" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set fwButtons $wwDialog.fwButtons
  set fwCheckbox $wwDialog.fwCheckbox


  # check button
  tkm_MakeCheckboxes $fwCheckbox h { \
          { text "Display Intermediate Results" \
           gDisplayIntermediateResults \
           { "SetSegmentationDisplayStatus" } \
           "Will show each iteration of Gibbs ICM algorithm" } \
       }

  # prompt
  tkm_MakeNormalLabel $fwMain\
    "Do you wish to recompute the segmentation?"

  # ok and cancel buttons.
  tkm_MakeButtons $fwButtons { \
    {text "Try" {RecomputeSegmentation 0} } \
    {text "Revert" {RestorePreviousSegmentation 0} } \
    {text "Cancel" { Dialog_Close .wwRecomputeSegmentation }}
    {text "Update Means" { Dialog_Close .wwRecomputeSegmentation } }}

  pack $fwMain $fwButtons $fwCheckbox \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}


proc DoCursorInfoDlog { } {

    global gDialog
    global gCursor

    set wwDialog .wwCursorInfoDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Cursor Info" {-borderwidth 10}] } {
  
	set fwTop                $wwDialog.fwTop 
	set fwColor              $fwTop.fwColor
	set lfwShape             $fwTop.lfwShape
	set fwButtons            $wwDialog.fwButtons
	
	frame $fwTop
	
	# color
	tkm_MakeColorPickers $fwColor \
	    [list [list "Color" gCursor(color,red) gCursor(color,green) \
		       gCursor(color,blue) "SendCursorConfiguration"]]

	# shape
	tixLabelFrame $lfwShape \
	    -label "Shape" \
	    -labelside acrosstop \
	    -options { label.padX 5 }
	
	set fwShapeSub           [$lfwShape subwidget frame]
	set fwShape              $fwShapeSub.fwShape
	
	tkm_MakeToolbar $fwShape \
	    1 \
	    gCursor(shape) \
	    {puts "cursor toolbar"} {
		{ image 0 icon_marker_crosshair "Crosshair" }
		{ image 1 icon_marker_diamond "Diamond" } }
	
	pack $fwShape \
	    -side left \
	    -anchor w
	
	# buttons. 
	tkm_MakeApplyCloseButtons $fwButtons $wwDialog SendCursorConfiguration
	
	pack $fwTop $fwColor $lfwShape $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x
    }
}

proc DoSurfaceInfoDlog { } {

    global gDialog
    global gSurface
    global Surf_tVertexSet
    global tkm_tSurfaceType
    global ksaSurfaceVertexSetString
    global ksaSurfaceTypeString
    global gbUseRealRAS

    set wwDialog .wwSurfaceInfoDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Surface Info" {-borderwidth 10}] } {
	
	set fwTop     $wwDialog.fwTop
	set fwButtons $wwDialog.fwButtons
	
	frame $fwTop

	tkm_MakeCheckboxes $fwTop.cbwUseRealRAS x {
	    { text "Use real RAS" gbUseRealRAS "SendUseRealRAS" } }
	
	tkm_MakeBigLabel $fwTop.lwColors "Colors"

	set lPickers {}
	foreach surface "$tkm_tSurfaceType(main) $tkm_tSurfaceType(aux)" {
	    foreach set "$Surf_tVertexSet(main) $Surf_tVertexSet(original) \
                         $Surf_tVertexSet(pial)" {
		
		lappend lPickers \
    [list "$ksaSurfaceTypeString($surface) $ksaSurfaceVertexSetString($set)" \
	     gSurface($surface,$set,color,red) \
	     gSurface($surface,$set,color,green) \
	     gSurface($surface,$set,color,blue) \
	     "SendSurfaceInformation"]
	    }
	}		
	tkm_MakeColorPickers $fwTop.fwColors $lPickers

	tkm_MakeBigLabel $fwTop.lwWidths "Widths"

	set lSliders {}
	foreach surface "$tkm_tSurfaceType(main) $tkm_tSurfaceType(aux)" {
	    foreach set "$Surf_tVertexSet(main) $Surf_tVertexSet(original) \
                         $Surf_tVertexSet(pial)" {

		lappend lSliders \
 [list "\"$ksaSurfaceTypeString($surface) $ksaSurfaceVertexSetString($set)\"" \
                gSurface($surface,$set,width) 1 20 80 \
                "SendSurfaceInformation" 1]
	    }
	}
	tkm_MakeSliders $fwTop.fwWidths $lSliders
	
	pack $fwTop.cbwUseRealRAS $fwTop.lwColors $fwTop.fwColors \
	    $fwTop.lwWidths $fwTop.fwWidths \
	    -side top \
	    -anchor w
	
	# buttons. 
	tkm_MakeCloseButton $fwButtons $wwDialog
	
	pack $fwTop $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x
    }
    
}

proc DoResolveUseRealRASDlog { } {

    global gDialog
    global gbUseRealRAS

    set wwDialog .wwResolveUseRealRASDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Use Real RAS" {-borderwidth 10}] } {
	
	set fwTop     $wwDialog.fwTop
	set fwButtons $wwDialog.fwButtons
	
	frame $fwTop

	tkm_MakeNormalLabel $fwTop.lwExplanation \
	    "The surface you are loading uses a different\nRAS coordinate space interpretation than you have\nspecified (or is the default setting)."

	if { $gbUseRealRAS } {
	    set sUseRealRAS_0 "Use surface RAS (surface setting)"
	    set sUseRealRAS_1 "Use real RAS (current setting)"
	} else {
	    set sUseRealRAS_0 "Use surface RAS (current setting)"
	    set sUseRealRAS_1 "Use real RAS (surface setting)"
	}

	tkm_MakeRadioButtons $fwTop.rbwUseRealRAS y "" gbUseRealRAS \
	    [list \
		 [list text $sUseRealRAS_0 0 ""] \
		 [list text $sUseRealRAS_1 1 ""] ] 

	pack $fwTop.lwExplanation $fwTop.rbwUseRealRAS \
	    -side top \
	    -anchor w
	
	# buttons. 
	tkm_MakeCancelOKButtons $fwButtons $wwDialog "SendUseRealRAS"
	
	pack $fwTop $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x
    }
    
}

proc DoEditSegBrushInfoDlog { } {

    global gDialog 
    global gSegBrush glSegEditColors
    global tkm_tVolumeTarget_MainAna
    global tkm_tVolumeTarget_AuxAna
    global tkm_tVolumeTarget_MainSeg
    global tkm_tVolumeTarget_AuxSeg

    set wwDialog .wwEditSegBrushInfoDlog
    
    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Segmentation Brush Info" {-borderwidth 10}] } {
	
	set lfwColor     $wwDialog.lfwColor
	set lfwFill      $wwDialog.lfwFill
	set fwButtons    $wwDialog.fwButtons
	
	# color
	tixLabelFrame $lfwColor \
	    -label "Color" \
	    -labelside acrosstop \
	    -options { label.padX 5}
	
	set fwColorSub           [$lfwColor subwidget frame]
	set fwColor              $fwColorSub.fwColor
	
	tixScrolledListBox $fwColor -scrollbar auto\
	    -browsecmd SendSegBrushInfo
	
	# go thru the list of entry names and insert each into the listbox
	$fwColor subwidget listbox configure -selectmode single
	set nLength [llength $glSegEditColors]
	foreach {structure name} $glSegEditColors {
	    $fwColor subwidget listbox insert end "$structure: $name"
	}
	
	# select the one with the index of the seg brush color
	$fwColor subwidget listbox selection set $gSegBrush(color)
	$fwColor subwidget listbox see $gSegBrush(color)
	
	pack $fwColor \
	    -side top \
	    -expand yes \
	    -fill both
	
	# fill characteristics
	tixLabelFrame $lfwFill \
	    -label "Fill Parameters" \
	    -labelside acrosstop \
	    -options { label.padX 5 }
	
	set fwFillSub       [$lfwFill subwidget frame]
	set fwFill          $fwFillSub.fwFill
	set fw3D            $fwFill.fw3D
	set fwLabel         $fwFill.fwLabel
	set fwMainSrc       $fwFill.fwMainSrc
	set fwAuxSrc        $fwFill.fwAuxSrc
	set fwSegSrc        $fwFill.fwSegSrc
	set fwSliders       $fwFill.fwSliders
	set fwDistanceNote  $fwFill.fwDistanceNote
	
	frame $fwFill
	
	# 3d
	tkm_MakeCheckboxes $fw3D y {
	    { text "3D" gSegBrush(3d) "SendSegBrushInfo" } }
	
	# source radios
	tkm_MakeNormalLabel $fwLabel "Use as source:"
	tkm_MakeRadioButton $fwMainSrc "Main Anatomical" \
	    gSegBrush(src) $tkm_tVolumeTarget_MainAna "SendSegBrushInfo"
	tkm_MakeRadioButton $fwAuxSrc "Aux Anatomical" \
	    gSegBrush(src) $tkm_tVolumeTarget_AuxAna "SendSegBrushInfo"
	tkm_MakeRadioButton $fwSegSrc "Segmentation" \
	    gSegBrush(src) $tkm_tVolumeTarget_MainSeg "SendSegBrushInfo"
	
	# fuzziness and max distance
	tkm_MakeSliders $fwSliders {
	    { "Fuzziness" gSegBrush(fuzzy) 
		0 255 50 "SendSegBrushInfo" 1 0.1 } 
	    {  "\"Max Distance\"" gSegBrush(distance) 
		0 255 50 "SendSegBrushInfo" 1 0.1 } }
	tkm_MakeSmallLabel $fwDistanceNote "enter 0 for no limit"
	
	
	pack $fw3D $fwLabel $fwMainSrc $fwAuxSrc \
	    $fwSegSrc $fwSliders $fwDistanceNote $fwFill \
	    -side top \
	    -expand yes \
	    -fill x
	
	# close button
	tkm_MakeCloseButton $fwButtons $wwDialog 
	
	grid $lfwColor  -column 0 -row 0 -sticky news
	grid $lfwFill   -column 0 -row 1 -sticky ews
	grid $fwButtons -column 0 -row 2 -sticky ews

	grid columnconfigure $wwDialog 0 -weight 1
	grid rowconfigure $wwDialog 0 -weight 1
	grid rowconfigure $wwDialog 1 -weight 0
	grid rowconfigure $wwDialog 2 -weight 0

    }
}

proc DoEditFloodSelectParamsDlog { } {

    global gDialog 
    global gFloodSelectParams
    global tkm_tVolumeTarget_MainAna
    global tkm_tVolumeTarget_AuxAna
    global tkm_tVolumeTarget_MainSeg
    global tkm_tVolumeTarget_AuxSeg

    set wwDialog .wwEditFloodSelectParams

   # try to create the dlog...
    if { [Dialog_Create $wwDialog "Flood Select" {-borderwidth 10}] } {

	
	set fwMain          $wwDialog.fwMain
	set fwButtons       $wwDialog.fwButtons

	set fw3D            $fwMain.fw3D
	set fwLabel         $fwMain.fwLabel
	set fwMainSrc       $fwMain.fwMainSrc
	set fwAuxSrc        $fwMain.fwAuxSrc
	set fwSegSrc        $fwMain.fwSegSrc
	set fwAuxSegSrc     $fwMain.fwAuxSegSrc
	set fwSliders       $fwMain.fwSliders
	set fwDistanceNote  $fwMain.fwDistanceNote
	
	frame $fwMain

	# 3d
	tkm_MakeCheckboxes $fw3D y {
	    { text "3D" gFloodSelectParams(3d) "SendFloodSelectParams" } }
	
	# source radios
	tkm_MakeNormalLabel $fwLabel "Use as source:"
	tkm_MakeRadioButton $fwMainSrc "Main Anatomical" \
	    gFloodSelectParams(src) $tkm_tVolumeTarget_MainAna \
	    "SendFloodSelectParams"
	tkm_MakeRadioButton $fwAuxSrc "Aux Anatomical" \
	    gFloodSelectParams(src) $tkm_tVolumeTarget_AuxAna \
	    "SendFloodSelectParams"
	tkm_MakeRadioButton $fwSegSrc "Main Segmentation" \
	    gFloodSelectParams(src) $tkm_tVolumeTarget_MainSeg \
	    "SendFloodSelectParams"
	tkm_MakeRadioButton $fwAuxSegSrc "Aux Segmentation" \
	    gFloodSelectParams(src) $tkm_tVolumeTarget_AuxSeg \
	    "SendFloodSelectParams"
	
	# fuzziness and max distance
	tkm_MakeSliders $fwSliders { 
	    { "Fuzziness" gFloodSelectParams(fuzzy) 
		0 255 50 "SendFloodSelectParams" 1 0.1 } 
	    {  "\"Max Distance\"" gFloodSelectParams(distance) 
		0 255 50 "SendFloodSelectParams" 1 0.1 } }
	tkm_MakeSmallLabel $fwDistanceNote "enter 0 for no limit"
	
	
	pack $fw3D $fwLabel $fwMainSrc $fwAuxSrc \
	    $fwSegSrc $fwAuxSegSrc $fwSliders $fwDistanceNote \
	    -side top \
	    -expand yes \
	    -fill x

	# close button
	tkm_MakeCloseButton $fwButtons $wwDialog 
	
	pack $fwMain $fwButtons \
	    -side top \
	    -expand yes \
	    -fill x
    }
}


proc DoVolumeColorScaleInfoDlog { } {
    global gInterface
    global gDialog
    global gVolume

    set wwDialog .wwVolumeColorScaleInfoDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog \
	      "Brightness / Contrast" {-borderwidth 10}] } {

	set fwMainSliders $wwDialog.fwMainSliders
	set fwMainMin     $wwDialog.fwMainMin
	set fwMainMax     $wwDialog.fwMainMax
	set fwAuxSliders  $wwDialog.fwAuxSliders
	set fwAuxMin      $wwDialog.fwAuxMin
	set fwAuxMax      $wwDialog.fwAuxMax
	set fwButtons     $wwDialog.fwButtons
	
	# brightness and contrast for main volume.
	tkm_MakeSliders $fwMainSliders [list \
		[list {"Brightness"} gVolume(0,colorScale,brightness) \
		     1 0 100 "SendVolumeBrightnessContrast" 1 0.01] \
		[list {"Contrast"} gVolume(0,colorScale,contrast) \
		     0 30 100 "SendVolumeBrightnessContrast" 1] ]

	# Min and max buttons for main volume. An entry for setting
	# the value, and a button to set the value to the volume's
	# value min or max.
	frame $fwMainMin
	tkm_MakeEntry $fwMainMin.ew "Main Min value" \
	    gVolume(0,colorScale,min) 10 "SendVolumeMinMax"
	button $fwMainMin.bw -text "Set to $gVolume(0,minValue)" \
	    -command "set gVolume(0,colorScale,min) $gVolume(0,minValue); SendVolumeMinMax"
	pack $fwMainMin.ew $fwMainMin.bw \
	    -side left

	frame $fwMainMax
	tkm_MakeEntry $fwMainMax.ew "Main Max value" \
	    gVolume(0,colorScale,max) 10 "SendVolumeMinMax"
	button $fwMainMax.bw -text "Set to $gVolume(0,maxValue)" \
	    -command "set gVolume(0,colorScale,max) $gVolume(0,maxValue); SendVolumeMinMax"
	pack $fwMainMax.ew $fwMainMax.bw \
	    -side left

	# And same for the aux volume.
	tkm_MakeSliders $fwAuxSliders [list \
	        [list {"Aux Brightness"} gVolume(1,colorScale,brightness) \
		     1 0 100 "SendVolumeBrightnessContrast" 1 0.01] \
		[list {"Aux Contrast"} gVolume(1,colorScale,contrast)  \
		     0 30 100 "SendVolumeBrightnessContrast" 1] ]

	frame $fwAuxMin
	tkm_MakeEntry $fwAuxMin.ew "Aux Min value" \
	    gVolume(1,colorScale,min) 10 "SendVolumeMinMax"
	button $fwAuxMin.bw -text "Set to $gVolume(1,minValue)" \
	    -command "set gVolume(1,colorScale,min) $gVolume(1,minValue); SendVolumeMinMax"
	pack $fwAuxMin.ew $fwAuxMin.bw \
	    -side left

	frame $fwAuxMax
	tkm_MakeEntry $fwAuxMax.ew "Aux Max value" \
	    gVolume(1,colorScale,max) 10 "SendVolumeMinMax"
	button $fwAuxMax.bw -text "Set to $gVolume(1,maxValue)" \
	    -command "set gVolume(1,colorScale,max) $gVolume(1,maxValue); SendVolumeMinMax"
	pack $fwAuxMax.ew $fwAuxMax.bw \
	    -side left

	# Save these slider var names so we can update them if the
	# values change.
	set gInterface(colorScaleDlog,brightnessSlider,0) $fwMainSliders.sw0
	set gInterface(colorScaleDlog,contrastSlider,0)   $fwMainSliders.sw1
	set gInterface(colorScaleDlog,minValueButton,0)   $fwMainMin.bw
	set gInterface(colorScaleDlog,maxValueButton,0)   $fwMainMax.bw
	set gInterface(colorScaleDlog,brightnessSlider,1) $fwAuxSliders.sw4
	set gInterface(colorScaleDlog,contrastSlider,1)   $fwAuxSliders.sw5
	set gInterface(colorScaleDlog,minValueButton,1)   $fwAuxMin.bw
	set gInterface(colorScaleDlog,maxValueButton,1)   $fwAuxMax.bw

	# buttons
	tkm_MakeCloseButton $fwButtons $wwDialog
	
	pack $fwMainSliders $fwMainMin $fwMainMax \
	    $fwAuxSliders $fwAuxMin $fwAuxMax \
	    $fwButtons  \
	    -side top    \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5
    }
}


proc DoThresholdDlog {} {

    global gDialog
    
    set wwDialog .wwThresholdDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Threshold" {-borderwidth 10}] } {

  set fwLabel       $wwDialog.fwLabel
  set fwAbove       $wwDialog.fwAbove
  set fwBelow       $wwDialog.fwBelow
  set fwSliders     $wwDialog.fwSliders
  set fwButtons     $wwDialog.fwButtons

  # label 
  tkm_MakeNormalLabel $fwLabel "Change all values"
  
  # direction radios
  tkm_MakeRadioButton $fwAbove "above" bAbove 1
  tkm_MakeRadioButton $fwBelow "below" bAbove 0

  # threshold value
  tkm_MakeSliders $fwSliders { \
    { {"this value"} nThreshold 0 255 200 "" 1 } \
    { {"to this value"} nNewValue 0 255 200 "" 1 } }

  # pack them in a column
  pack $fwLabel $fwAbove $fwBelow $fwSliders \
    -side top                \
    -anchor w                \
    -expand yes              \
    -fill x

  # buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { ThresholdVolume $nThreshold $bAbove $nNewValue }

  pack  $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoRotateVolumeDlog {} {

    global gDialog
    
    set wwDialog .wwRotateVolumeDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Rotate Volume" {-borderwidth 10}] } {

  set fwDegrees     $wwDialog.fwRotateDegrees
  set fwDirection   $wwDialog.fwRotateDirection
  set fwX           $wwDialog.fwRotateX
  set fwY           $wwDialog.fwRotateY
  set fwZ           $wwDialog.fwRotateZ
  set fwButtons     $wwDialog.fwButtons

  set fRotateDegrees 0
  set sRotateDirection x

  # degrees
  tkm_MakeEntry \
    $fwDegrees "Degrees" \
    fRotateDegrees 5


  # direction radios
  tkm_MakeNormalLabel $fwDirection "Around anatomical axis:"
  # these are switched to match the internal representation of
  # the cor structure. x is ears, y is nose, z is thru top of head.
  tkm_MakeRadioButton $fwX "X (Ear to ear, perpendicular to Sagittal plane)" sRotateDirection x
  tkm_MakeRadioButton $fwY "Y (Back of head to nose, perpendicular to Coronal plane)" sRotateDirection z
  tkm_MakeRadioButton $fwZ "Z (Neck to top of head, perpendicular to Horizontal plane)" sRotateDirection y

  # buttons. 
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { RotateVolume $sRotateDirection $fRotateDegrees }

  pack $fwDegrees $fwDirection $fwX $fwY $fwZ $fwButtons \
    -side top       \
    -anchor w       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}
proc DoFlipVolumeDlog {} {

    global gDialog
    
    set bFlipX 0
    set bFlipY 0
    set bFlipZ 0
    set wwDialog .wwFlipDialog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Flip Volume" {-borderwidth 10}] } {

  set fwFlip    $wwDialog.fwFlip
  set fwButtons $wwDialog.fwButtons

  # flip checks
  # these are switched to match the internal representation of
  # the cor structure. x is ears, y is nose, z is thru top of head.
  tkm_MakeCheckboxes $fwFlip y { \
    { text "Flip around middle Sagittal plane" bFlipX {} "" } \
    { text "Flip around middle Horizontal plane" bFlipY {} "" } \
    { text "Flip around middle Coronal plane" bFlipZ {} "" } }

  # buttons. 
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { FlipVolume $bFlipX $bFlipY $bFlipZ }

  pack $fwFlip $fwButtons \
    -side top       \
    -anchor w       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoRegisterHeadPtsDlog {} {

    global gDialog
    
    set wwDialog .wwRegisterHeadPtsDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Register Head Points" {-borderwidth 10}] } {

  set fwTop        $wwDialog.fwTop
  set lfwTranslate $fwTop.lfwTranslate
  set lfwRotate    $fwTop.lfwRotate
  set fwButtons    $fwTop.fwButtons

  frame $fwTop

  # make the label frames and get their subs.
  tixLabelFrame $lfwTranslate \
    -label "Translate" \
    -labelside acrosstop \
    -options { label.padX 5 }

  set fwTranslateSub     [$lfwTranslate subwidget frame]
  set fwTranslateButtons $fwTranslateSub.fwTranslateButtons
  set fwTranslateAmt     $fwTranslateSub.fwTranslateAmt

  # puts buttons in them
  tkm_MakeButtons $fwTranslateButtons { \
    { image icon_arrow_up \
    "TranslateHeadPts $fTranslateDistance y" } \
    { image icon_arrow_down \
    "TranslateHeadPts -$fTranslateDistance y" } \
    { image icon_arrow_left \
    "TranslateHeadPts $fTranslateDistance x" } \
    { image icon_arrow_right \
    "TranslateHeadPts -$fTranslateDistance x" } }

  tkm_MakeEntryWithIncDecButtons \
    $fwTranslateAmt "Distance" \
    fTranslateDistance \
    {} \
    0.5

  pack $fwTranslateButtons $fwTranslateAmt \
    $lfwTranslate \
    -side top                           \
    -anchor w                           \
    -expand yes                         \
    -fill x

  # rotate frame
  tixLabelFrame $lfwRotate \
    -label "Rotate" \
    -labelside acrosstop \
    -options { label.padX 5 }

  set fwRotateSub     [$lfwRotate subwidget frame]
  set fwRotateButtons $fwRotateSub.fwRotateButtons
  set fwRotateAmt     $fwRotateSub.fwRotateAmt

  # puts buttons in them
  tkm_MakeButtons $fwRotateButtons { \
    { image icon_arrow_ccw \
    "RotateHeadPts $fRotateDegrees z" } \
    { image icon_arrow_cw \
    "RotateHeadPts -$fRotateDegrees z" } }

  tkm_MakeEntryWithIncDecButtons \
    $fwRotateAmt "Degrees" \
    fRotateDegrees \
    {} \
    0.5

  pack $fwRotateButtons $fwRotateAmt \
    $lfwRotate \
    -side top                     \
    -anchor w                     \
    -expand yes                   \
    -fill x

  # just a close button here
  tkm_MakeButtons $fwButtons { \
    { text "Close" {Dialog_Close .wwRegisterHeadPtsDlog} } }

  pack $fwButtons \
    -side right \
    -anchor e

  pack $fwTop
    }
}

proc DoFindVertexDlog { iSurface } {

    global gDialog
    global gFindingSurface

    set gFindingSurface $iSurface
    set nVertex 0
    set wwDialog .wwFindVertexDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Find Vertex" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set lwPrompt  $fwMain.lwPrompt
  set ewName    $fwMain.ewName

  set fwButtons $wwDialog.fwButtons
  set bwOK      $fwButtons.bwOK
  set bwCancel  $fwButtons.bwCancel

  # prompt and entry field
  tkm_MakeEntry $fwMain "Find vertex number:" nVertex 6
  
  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { FindVertex $nVertex } {}

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoRegisterOverlayDlog {} {

    global gDialog
    global fScaleFactor

    set wwDialog .wwRegisterOverlayDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Register Functional Overlay" {-borderwidth 10}] } {

  set fwTop        $wwDialog.fwTop
  set lfwTranslate $fwTop.lfwTranslate
  set lfwRotate    $fwTop.lfwRotate
  set lfwScale     $fwTop.lfwScale
  set fwButtons    $fwTop.fwButtons
  
  frame $fwTop

  # make the label frames and get their subs.
  tixLabelFrame $lfwTranslate \
    -label "Translate" \
    -labelside acrosstop \
    -options { label.padX 5 }

  set fwTranslateSub     [$lfwTranslate subwidget frame]
  set fwTranslateButtons $fwTranslateSub.fwTranslateButtons
  set fwTranslateAmt     $fwTranslateSub.fwTranslateAmt

  # puts buttons in them
  tkm_MakeButtons $fwTranslateButtons { \
    { image icon_arrow_up \
    "TranslateOverlayRegistration $fTranslateDistance y" } \
    { image icon_arrow_down \
    "TranslateOverlayRegistration -$fTranslateDistance y" } \
    { image icon_arrow_left \
    "TranslateOverlayRegistration $fTranslateDistance x" } \
    { image icon_arrow_right \
    "TranslateOverlayRegistration -$fTranslateDistance x" } }

  tkm_MakeEntryWithIncDecButtons \
    $fwTranslateAmt "Distance" \
    fTranslateDistance \
    {} \
    0.5

  pack $fwTranslateButtons $fwTranslateAmt \
    $lfwTranslate \
    -side top                           \
    -anchor w                           \
    -expand yes                         \
    -fill x

  # rotate frame
  tixLabelFrame $lfwRotate \
    -label "Rotate" \
    -labelside acrosstop \
    -options { label.padX 5 }

  set fwRotateSub     [$lfwRotate subwidget frame]
  set fwRotateButtons $fwRotateSub.fwRotateButtons
  set fwRotateAmt     $fwRotateSub.fwRotateAmt

  # puts buttons in them
  tkm_MakeButtons $fwRotateButtons { \
    { image icon_arrow_ccw \
    "RotateOverlayRegistration $fRotateDegrees z" } \
    { image icon_arrow_cw \
    "RotateOverlayRegistration -$fRotateDegrees z" } }

  tkm_MakeEntryWithIncDecButtons \
    $fwRotateAmt "Degrees" \
    fRotateDegrees \
    {} \
    0.5

  pack $fwRotateButtons $fwRotateAmt \
    $lfwRotate \
    -side top                     \
    -anchor w                     \
    -expand yes                   \
    -fill x
  
  # scale frame
  tixLabelFrame $lfwScale \
    -label "Scale" \
    -labelside acrosstop \
    -options { label.padX 5 }

  set fwScaleSub     [$lfwScale subwidget frame]
  set fwScaleButtons $fwScaleSub.fwScaleButtons
  set fwScaleAmt     $fwScaleSub.fwScaleAmt

  # puts buttons in them
  tkm_MakeButtons $fwScaleButtons { \
    { image icon_arrow_expand_x \
    "ScaleOverlayRegistration $fScaleFactor x" } \
    { image icon_arrow_shrink_x \
    "ScaleOverlayRegistration [expr 1.0 / $fScaleFactor] x" }
    { image icon_arrow_expand_y \
    "ScaleOverlayRegistration $fScaleFactor y" } \
    { image icon_arrow_shrink_y \
    "ScaleOverlayRegistration [expr 1.0 / $fScaleFactor] y" } }

  set fScaleFactor 1.0
  tkm_MakeEntryWithIncDecButtons \
    $fwScaleAmt "Scale Factor" \
    fScaleFactor \
    {} \
    0.05

  pack $fwScaleButtons $fwScaleAmt \
    $lfwScale \
    -side top                           \
    -anchor w                           \
    -expand yes                         \
    -fill x

  # just a close button here
  tkm_MakeButtons $fwButtons { \
    { text "Close" {Dialog_Close .wwRegisterOverlayDlog} } }

  pack $fwButtons \
    -side right \
    -anchor e

  pack $fwTop
    }
}

proc DoSaveRGBSeriesDlog {} {

    global gDialog
    global gnVolSlice

    set sDir ""
    set sPrefix ""
    set nBegin $gnVolSlice
    set nEnd $gnVolSlice
    set wwDialog .wwSaveRGBSeriesDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save RGB Series" {-borderwidth 10}] } {

	set fwDir     $wwDialog.fwDir
	set fwPrefix  $wwDialog.fwPrefix
	set fwDirection $wwDialog.fwDirection
	set fwBegin   $wwDialog.fwBegin
	set fwEnd     $wwDialog.fwEnd
	set fwButtons $wwDialog.fwButtons
	
	# the directory field
	tkm_MakeDirectorySelector $fwDir \
	    "Directory to save files in:" sDir
	
	# the file prefix
	tkm_MakeEntry $fwPrefix "Prefix:" sPrefix
	
	# begin and end slices
	tkm_MakeEntryWithIncDecButtons $fwBegin \
	    "From slice" nBegin {} 1
	tkm_MakeEntryWithIncDecButtons $fwEnd \
	    "To slice" nEnd {} 1
	
	# ok and cancel buttons.
	tkm_MakeCancelOKButtons $fwButtons $wwDialog \
	    "SaveRGBSeries \$sDir/\$sPrefix \$nBegin \$nEnd"
	
	pack $fwDir $fwPrefix $fwBegin $fwEnd $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5
    }
}

proc DoSaveTIFFSeriesDlog {} {

    global gDialog
    global gnVolSlice

    set sDir ""
    set sPrefix ""
    set nBegin $gnVolSlice
    set nEnd $gnVolSlice
    set wwDialog .wwSaveTIFFSeriesDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save TIFF Series" {-borderwidth 10}] } {

	set fwDir     $wwDialog.fwDir
	set fwPrefix  $wwDialog.fwPrefix
	set fwDirection $wwDialog.fwDirection
	set fwBegin   $wwDialog.fwBegin
	set fwEnd     $wwDialog.fwEnd
	set fwButtons $wwDialog.fwButtons
	
	# the directory field
	tkm_MakeDirectorySelector $fwDir \
	    "Directory to save files in:" sDir
	
	# the file prefix
	tkm_MakeEntry $fwPrefix "Prefix:" sPrefix
	
	# begin and end slices
	tkm_MakeEntryWithIncDecButtons $fwBegin \
	    "From slice" nBegin {} 1
	tkm_MakeEntryWithIncDecButtons $fwEnd \
	    "To slice" nEnd {} 1
	
	# ok and cancel buttons.
	tkm_MakeCancelOKButtons $fwButtons $wwDialog \
	    "SaveTIFFSeries \$sDir/\$sPrefix \$nBegin \$nEnd"
	
	pack $fwDir $fwPrefix $fwBegin $fwEnd $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5
    }
}

proc DoLabelWriterHelperDlog {} {

    global gDialog glShortcutDirs
    global sDir sPrefix sNumber sSuffix
    global sNamingMethod
    global bClearLabel
    global bSelectLine
    global nLabelFileNameLine sNextLabelName

    set sDir [GetDefaultLocation SaveLabelAs]
    set sPrefix "line-"
    set sNumber 0
    set sSuffix ".label"

    set wwDialog .wwLabelWriterHelperDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Label Writer Helper" {-borderwidth 10}] } {

	set fwDir       $wwDialog.fwDir

	set fwOptions   $wwDialog.fwOptions

	set fwUseAutoInc   $wwDialog.fwUseAutoInc
	set fwUseAutoIncCB $fwUseAutoInc.fwUseAutoIncCB
	set fwAutoIncFileName $fwUseAutoInc.fwAutoIncFileName
	set fwPrefix       $fwAutoIncFileName.fwPrefix
	set fwNumber       $fwAutoIncFileName.fwNumber
	set fwSuffix       $fwAutoIncFileName.fwSuffix

	set fwUseFile          $wwDialog.fwUseFile
	set fwUseFileCB        $fwUseFile.fwUseFileCB
	set fwUseFileSelector  $fwUseFile.fwUseFileSelector
	set fwNextFileName     $fwUseFile.fwNextFileName
	set fwButtons   $wwDialog.fwButtons

	# the directory field
	tkm_MakeDirectorySelector $fwDir \
	    "Directory to save files in:" sDir

	# Options.
	tkm_MakeCheckboxes $fwOptions y {
	    { text "Select Line Before Writing" bSelectLine {} }
	    { text "Clear Label After Writing" bClearLabel {} }
	}

	# For using an auto incrementing file name.
	frame $fwUseAutoInc -bd 1 -relief ridge
	tkm_MakeRadioButtons $fwUseAutoIncCB y "" sNamingMethod {
	    { text "Use Incrementing Number for Label Name" autoInc {} }
	}
	frame $fwAutoIncFileName
	tkm_MakeEntry $fwPrefix "File Name:" sPrefix 10
	tkm_MakeEntry $fwNumber "" sNumber 4
	tkm_MakeEntry $fwSuffix "" sSuffix 10
	

	# For getting label names from a file.
	set nLabelFileNameLine 0
	frame $fwUseFile -bd 1 -relief ridge
	tkm_MakeRadioButtons $fwUseFileCB y "" sNamingMethod {
	    { text "Get Label Names from File" useFile {} }
	}
	tkm_MakeFileSelector $fwUseFileSelector "File" fnLabelSource \
	    [list GetDefaultLocation LoadLabel] $glShortcutDirs
	tkm_MakeActiveLabel $fwNextFileName "Next name" sNextLabelName 30

	# ok and cancel buttons.
	tkm_MakeApplyCloseButtons $fwButtons $wwDialog {
	    if { $bSelectLine } { AddLineToSelection }
	    set fnLabel ""
	    if { "$sNamingMethod" == "autoInc" } {
		set fnLabel $sDir/$sPrefix$sNumber$sSuffix
		set sNumber [expr $sNumber + 1]
 	    } elseif { "$sNamingMethod" == "useFile" } {
		if { "$fnLabelSource" == "" } {
		    ErrorDlog "Please specify a source file first."
		    return
		} else {
		    set r [catch { 
			set fnLabel [GetNthLineFromFile \
					 $fnLabelSource $nLabelFileNameLine]
			incr nLabelFileNameLine
			set sNextLabelName [GetNthLineFromFile \
					 $fnLabelSource $nLabelFileNameLine]
		    } sError]
		    if { $r != 0 } { ErrorDlog "$sError" }
		}
	    }
	    SaveLabel $fnLabel
	    if { $bClearLabel } { ClearSelection }
	    RedrawAll
	}

	# Set the initial naming method.
	set sNamingMethod autoInc

	pack $fwPrefix $fwNumber $fwSuffix -side left

	pack $fwUseAutoIncCB $fwAutoIncFileName \
	    -side top       \
	    -expand yes     \
	    -fill x         \

	pack $fwUseFileCB $fwUseFileSelector $fwNextFileName \
	    -side top       \
	    -expand yes     \
	    -fill x         \

	pack $fwDir  $fwOptions $fwUseAutoInc $fwUseFile $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5
    }
}

proc GetNthLineFromFile { ifn inLine } {

    # Try to open the file.
    set rOpen [catch {set fp [open $ifn r]} sError]
    if { $rOpen != 0 } {
	error "Couldn't open file $ifn."
    }

    # Go through the file up to inLine, rewinding to the beginning if
    # we need to.
    set nLine 0
    set sLine ""
    while { $nLine <= $inLine } {
	gets $fp sLine
	if { [eof $fp] } {
	    seek $fp 0 start
	    gets $fp sLine
	}
	incr nLine
    }
	
    close $fp

    return $sLine
}

proc DoGotoPointDlog {} {

    global gDialog
    global mri_tCoordSpace_VolumeIdx mri_tCoordSpace_RAS 
    global mri_tCoordSpace_Talairach mri_tCoordSpace_SurfaceRAS
    global gnVolX gnVolY gnVolZ
    global gbTalTransformPresent
    
    set wwDialog .wwGotoPointDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Goto Point" {-borderwidth 10}] } {

  set fwLabel       $wwDialog.fwLabel
  set fwCoordSpace  $wwDialog.fwCoordSpace
  set fwVolumeIdx   $fwCoordSpace.fwVolumeIdx
  set fwRAS         $fwCoordSpace.fwSurfaceRAS
  set fwSurfaceRAS  $fwCoordSpace.fwRAS
  set fwTalCoords   $fwCoordSpace.fwTalCoords
  set fwWhere       $wwDialog.fwWhere
  set fwX           $fwWhere.fwX
  set fwY           $fwWhere.fwY
  set fwZ           $fwWhere.fwZ
  set fwButtons     $wwDialog.fwButtons

  set fX $gnVolX(cursor)
  set fY $gnVolY(cursor)
  set fZ $gnVolZ(cursor)
  set coordSpace $mri_tCoordSpace_VolumeIdx

  # coord space radios
  tkm_MakeNormalLabel $fwLabel "Coordinate space:"
  frame $fwCoordSpace
  tkm_MakeRadioButton $fwVolumeIdx "Volume Index" \
    coordSpace $mri_tCoordSpace_VolumeIdx
  tkm_MakeRadioButton $fwSurfaceRAS "Surface RAS" \
    coordSpace $mri_tCoordSpace_SurfaceRAS
  tkm_MakeRadioButton $fwRAS "RAS" \
    coordSpace $mri_tCoordSpace_RAS
  pack $fwLabel $fwVolumeIdx $fwSurfaceRAS $fwRAS \
    -side left

  # pack tal coords if we got 'em
  if { $gbTalTransformPresent == 1 } {

      tkm_MakeRadioButton $fwTalCoords "Talairach" \
        coordSpace $mri_tCoordSpace_Talairach
      pack $fwTalCoords \
        -side left
  }

  # x y z fields
  frame $fwWhere
  tkm_MakeEntry $fwX "X" fX 5
  tkm_MakeEntry $fwY "Y" fY 5
  tkm_MakeEntry $fwZ "Z" fZ 5
  pack $fwX $fwY $fwZ \
    -side left

  # buttons. 
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { SetCursor $coordSpace $fX $fY $fZ }

  pack $fwLabel $fwCoordSpace $fwWhere $fwButtons \
    -side top       \
    -anchor w       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoAverageSurfaceVertexPositionsDlog {} {

    global gDialog

    set nNumAverages 10
    set wwDialog .wwFindVertexDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Average Vertex Positions" {-borderwidth 10}] } {

	set fwMain    $wwDialog.fwMain
	set lwPrompt  $fwMain.lwPrompt
	set ewName    $fwMain.ewName
	
	set fwButtons $wwDialog.fwButtons
	set bwOK      $fwButtons.bwOK
	set bwCancel  $fwButtons.bwCancel
	
	# prompt and entry field
	tkm_MakeEntry $fwMain "Number of averages:" nNumAverages 6
  
	# ok and cancel buttons.
	tkm_MakeCancelOKButtons $fwButtons $wwDialog \
	    { AverageSurfaceVertexPositions $nNumAverages } {}
	
	pack $fwMain $fwButtons \
	    -side top       \
	    -expand yes     \
	    -fill x         \
	    -padx 5         \
	    -pady 5
    }
}

proc DoEditHeadPointLabelDlog {} {

    global gDialog
    global gsaLabelContents
    
    set wwDialog .wwEditHeadPointLabelDlog
    set sHeadPointLabel $gsaLabelContents(kLabel_Label_Head,value,cursor)

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Edit Head Point Label" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain

  set fwButtons $wwDialog.fwButtons
  set bwOK      $fwButtons.bwOK
  set bwCancel  $fwButtons.bwCancel

  # prompt and entry field
  tkm_MakeEntry $fwMain "Change label to:" sHeadPointLabel 20
  
  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { SetSelectedHeadPointLabel $sHeadPointLabel } {}

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}



proc SetSegmentationDisplayStatus { } {

    global gDisplayIntermediateResults
    SetGCADisplayStatus $gDisplayIntermediateResults
}
proc SetBrushConfiguration { } {
    global gBrushInfo

    SetBrushTarget $gBrushInfo(target)
    SetBrushShape $gBrushInfo(radius) $gBrushInfo(shape) $gBrushInfo(3d)
}

proc SetEditBrushConfiguration { } {

    global gEditBrush
    global DspA_tBrush_EditOne DspA_tBrush_EditTwo

    foreach tool "$DspA_tBrush_EditOne $DspA_tBrush_EditTwo" {
	SetBrushInfo $tool \
	    $gEditBrush($tool,low) \
	    $gEditBrush($tool,high) \
	    $gEditBrush($tool,new) \
	    $gEditBrush($tool,mode) \
	    $gEditBrush($tool,cloneSource)
    }
}

proc SendFillAnatomicalInfo { } {
    global gBrushInfo
    SetAnatomicalFillInfo \
	$gBrushInfo(3dfill) $gBrushInfo(fuzzy) $gBrushInfo(distance)

}

proc SetSegmentationVolumeConfiguration {} {

    global gfSegmentationVolumeAlpha
    SetSegmentationAlpha $gfSegmentationVolumeAlpha
}

proc SetDTIVolumeConfiguration {} {

    global gfDTIVolumeAlpha
    SetDTIAlpha $gfDTIVolumeAlpha
}

proc SendSegBrushInfo {} {
    global gSegBrush
    global glSegEditColors

    # Get the selected item from the scrolled box in the dlog. 
    set nSelection [[.wwEditSegBrushInfoDlog.lfwColor subwidget frame].fwColor subwidget listbox curselection]

    # Now we need to find the corresponding segmentation index in our
    # list of eit colors. The list is in the format {structure name
    # structure name...} so it will be the (nSelection*2)th element.
    set nStructure [lindex $glSegEditColors [expr $nSelection * 2]]

    SetSegBrushInfo $nStructure $gSegBrush(3d) \
      $gSegBrush(src) $gSegBrush(fuzzy) \
      $gSegBrush(distance)
}

proc SendFloodSelectParams {} {
    global gFloodSelectParams

    SetFloodSelectParams $gFloodSelectParams(3d) \
      $gFloodSelectParams(src) $gFloodSelectParams(fuzzy) \
      $gFloodSelectParams(distance)
}

proc SendVolumeColorScale { } {

    global tkm_tVolumeType
    global gVolume

    foreach volume "$tkm_tVolumeType(main) $tkm_tVolumeType(aux)" {
	SetVolumeColorScale $volume \
	    $gVolume($volume,colorScale,brightness) \
	    $gVolume($volume,colorScale,contrast) \
	    $gVolume($volume,colorScale,min) \
	    $gVolume($volume,colorScale,max) 
    }
}

proc SendVolumeBrightnessContrast { } {

    global tkm_tVolumeType
    global gVolume

    foreach volume "$tkm_tVolumeType(main) $tkm_tVolumeType(aux)" {
	SetVolumeBrightnessContrast $volume \
	    $gVolume($volume,colorScale,brightness) \
	    $gVolume($volume,colorScale,contrast)
    }
}

proc SendVolumeMinMax { } {

    global tkm_tVolumeType
    global gVolume

    foreach volume "$tkm_tVolumeType(main) $tkm_tVolumeType(aux)" {
	SetVolumeMinMax $volume \
	    $gVolume($volume,colorScale,min) \
	    $gVolume($volume,colorScale,max)
    }
}

# ======================================================== INTERFACE MODIFIERS


proc ShowLabel { isLabel ibShow } {

    global gbShowLabel
    PackLabel $isLabel cursor $ibShow 
    PackLabel $isLabel mouseover $ibShow
    set gbShowLabel($isLabel) $ibShow
}

proc PackLabel { isLabel iSet ibShow } {

    global glLabel gfwaLabel

    # find the label index in our list.
    set nLabel [lsearch -exact $glLabel $isLabel]
    if { $nLabel == -1 } {
  puts "Couldn't find $isLabel\n"
  return;
    }

    # are we showing or hiding?
    if { $ibShow == 1 } {

  # go back and try to pack it after the previous labels
  set lTemp [lrange $glLabel 0 [expr $nLabel - 1]]
  set lLabelsBelow ""
  foreach element $lTemp {
      set lLabelsBelow [linsert $lLabelsBelow 0 $element]
  }
  foreach label $lLabelsBelow {
      if {[catch { pack $gfwaLabel($isLabel,$iSet) \
        -after $gfwaLabel($label,$iSet)    \
        -side top                \
        -anchor w } sResult] == 0} {
    return;
      }
  }
  
  # if that fails, go forward and try to pack it before the later labels
  set lLabelsAbove [lrange $glLabel [expr $nLabel + 1] [llength $glLabel]]
  foreach label $lLabelsAbove {
      if {[catch { pack $gfwaLabel($isLabel,$iSet)  \
        -before $gfwaLabel($label,$iSet)    \
        -side top                  \
        -anchor w } sResult] == 0} {
    return;
      }
  }

  # must be the first one. just pack it.
  catch { pack $gfwaLabel($isLabel,$iSet)  \
    -side top                \
    -anchor w } sResult

    } else {
  
  # else just forget it
  pack forget $gfwaLabel($isLabel,$iSet)
    } 
}

proc ShowVolumeCoords { ibShow } {
    ShowLabel kLabel_Coords_Vol $ibShow
}

proc ShowRASCoords { ibShow } {
    ShowLabel kLabel_Coords_Vol_RAS $ibShow
}
 
proc ShowTalCoords { ibShow } {
    global gbTalTransformPresent
    set gbTalTransformPresent $ibShow
    ShowLabel kLabel_Coords_Vol_Tal $ibShow
}

proc ShowAuxValue { ibShow } {
    ShowLabel kLabel_Value_Aux $ibShow
}

proc ShowSegLabel { ibShow } {
    ShowLabel kLabel_Label_SegLabel $ibShow
}

proc ShowAuxSegLabel { ibShow } {
    ShowLabel kLabel_Label_AuxSegLabel $ibShow
}

proc ShowHeadPointLabel { ibShow } {
    ShowLabel kLabel_Label_Head $ibShow
    tkm_SetEnableGroupStatus tMenuGroup_HeadPoints $ibShow
}

proc ShowFuncCoords { ibShow } {
    ShowLabel kLabel_Coords_Func $ibShow
}

proc ShowFuncValue { ibShow } {
    ShowLabel kLabel_Value_Func $ibShow
}

proc ClearSegColorTable { } {

    global glSegEditColors
    set glSegEditColors {}
}

proc AddSegColorTableEntry { inIndex isString } {

    global glSegEditColors
    lappend glSegEditColors $inIndex "$isString"
}

# =============================================================== VIEW PRESETS

proc SetViewPreset { iPreset } {

    global gViewPreset gViewPresetString ksaViewPresetString
    global tViewPreset_Single tViewPreset_Multiple tViewPreset_Mosaic
    global MWin_tLinkPolicy_None MWin_tLinkPolicy_MultipleOrientations
    global MWin_tLinkPolicy_Mosaic

    
    if { [catch {set gViewPreset $iPreset} sResult] == 1 } {
#  puts "caught: $sResult"
    }
    

    if { $iPreset == $tViewPreset_Single } {
  SetDisplayConfig 1 1 $MWin_tLinkPolicy_None
    }
    if { $iPreset == $tViewPreset_Multiple } {
  SetDisplayConfig 2 2 $MWin_tLinkPolicy_MultipleOrientations
    }
    if { $iPreset == $tViewPreset_Mosaic } {
  SetDisplayConfig 4 4 $MWin_tLinkPolicy_Mosaic
    }

}

# ========================================================= BUILDING INTERFACE

proc CreateWindow { iwwTop } {

    global ksWindowName
    frame $iwwTop
    wm title . $ksWindowName
    wm withdraw .
}

proc ToggleAndSendDisplayFlag { iFlag } {
    global gbDisplayFlag

    set gbDisplayFlag($iFlag) [expr 1 - $gbDisplayFlag($iFlag)]
    SendDisplayFlagValue $iFlag
}

proc SetAndSendDisplayFlag { iFlag iValue } {
    global gbDisplayFlag

    set gbDisplayFlag($iFlag) $iValue
    SendDisplayFlagValue $iFlag
}

proc MakeKeyBindings { iwTop } {

    global Surf_tVertexSet
    global gnZoomLevel

    bind $iwTop <Control-Key-1> \
	{SetAndSendDisplayFlag flag_AuxVolume 0}
    bind $iwTop <Control-Key-2> \
	{SetAndSendDisplayFlag flag_AuxVolume 1}
    bind $iwTop <Control-Key-m> \
	{ToggleAndSendDisplayFlag flag_MainSurface}
    bind $iwTop <Alt-Key-m> \
	{FindNearestSurfaceVertex $Surf_tVertexSet(main)}
    bind $iwTop <Control-Key-o> \
	{ToggleAndSendDisplayFlag flag_OriginalSurface}
    bind $iwTop <Alt-Key-o> \
	{FindNearestSurfaceVertex $Surf_tVertexSet(original)}
    bind $iwTop <Control-Key-p> \
	{ToggleAndSendDisplayFlag flag_PialSurface}
    bind $iwTop <Alt-Key-p> \
	{FindNearestSurfaceVertex $Surf_tVertexSet(pial)}
    bind $iwTop <Control-Key-v> \
	{ToggleAndSendDisplayFlag flag_DisplaySurfaceVertices}
    bind $iwTop <Control-Key-i> \
	{ToggleAndSendDisplayFlag flag_InterpolateSurfaceVertices}
    bind $iwTop <Control-Key-f> \
	{ToggleAndSendDisplayFlag flag_FunctionalColorScaleBar}
    bind $iwTop <Control-Key-g> \
	{ToggleAndSendDisplayFlag flag_SegmentationVolumeOverlay}
    bind $iwTop <Alt-Key-g> \
	{ToggleAndSendDisplayFlag flag_AuxSegmentationVolumeOverlay}
    bind $iwTop <Alt-Key-d> \
	{ToggleAndSendDisplayFlag flag_DTIOverlay}
    bind $iwTop <Control-Key-s> \
	{ToggleAndSendDisplayFlag flag_Selection}
    bind $iwTop <Control-Key-t> \
	{ToggleAndSendDisplayFlag flag_ControlPoints}
    bind $iwTop <Control-Key-c> \
	{ToggleAndSendDisplayFlag flag_Cursor}
    bind $iwTop <Control-Key-plus> \
	{SetZoomLevelWrapper [expr $gnZoomLevel / 2] }
    bind $iwTop <Control-Key-minus> \
	{SetZoomLevelWrapper [expr $gnZoomLevel / 2] }
    bind $iwTop <Control-q> \
	{ AllowSaveThenQuit }
    bind $iwTop <Control-b> \
	{ DoVolumeColorScaleInfoDlog }
    bind $iwTop <Control-u> \
	{ DoSurfaceInfoDlog }

    bind $iwTop <Alt-Key-plus> \
	{SetZoomLevelWrapper 32 }
    bind $iwTop <Alt-Key-minus> \
	{SetZoomLevelWrapper 1 }

    bind $iwTop <Key-n> \
	{SetTool $DspA_tTool_Navigate}
    bind $iwTop <Key-s> \
	{SetTool $DspA_tTool_Select}
    bind $iwTop <Key-a> \
	{SetTool $DspA_tTool_Edit}
    bind $iwTop <Key-g> \
	{SetTool $DspA_tTool_EditSeg}
    bind $iwTop <Key-t> \
	{SetTool $DspA_tTool_CtrlPts}
    bind $iwTop <Key-l> \
	{SetTool $DspA_tTool_Line}

}

proc CreateMenuBar { ifwMenuBar } {

    global mri_tOrientation_Sagittal mri_tOrientation_Horizontal 
    global mri_tOrientation_Coronal
    global DspA_tTool_Navigate DspA_tTool_Select
    global DspA_tTool_Edit DspA_tTool_EditSeg  DspA_tTool_CtrlPts 
    global gDisplayCols gDisplayRows gViewPreset
    global tViewPreset_Single tViewPreset_Multiple tViewPreset_Mosaic    
    global gTool
    global gbShowToolBar gbShowLabel
    global glDisplayFlag gbDisplayFlag
    global gControlPointsMenuCommand

    set mbwFile   $ifwMenuBar.mbwFile
    set mbwEdit   $ifwMenuBar.mbwEdit
    set mbwView   $ifwMenuBar.mbwView
    set mbwTools  $ifwMenuBar.mbwTools

    frame $ifwMenuBar -border 2 -relief raised

    # file menu button
    tkm_MakeMenu $mbwFile "File" {
	{ command
	    "Load Main Volume..."
	    {DoFileDlog LoadVolume} }
	{ command
	    "Save Main Volume"
	    DoSaveDlog
	    tMenuGroup_DirtyAnatomicalVolume }
	{ command
	    "Save Main Volume As..."
	    {DoFileDlog SaveVolumeAs}
	    tMenuGroup_DirtyAnatomicalVolume }
	{ cascade "Aux Volume" {
	    { command
		"Load Aux Volume..."
		{DoFileDlog LoadAuxVolume} }
	    { command
		"Save Aux Volume"
		DoAuxSaveDlog
		tMenuGroup_DirtyAuxAnatomicalVolume }
	    { command
		"Save Aux Volume As..."
		{DoFileDlog SaveAuxVolumeAs}
		tMenuGroup_DirtyAuxAnatomicalVolume }
	    { command
		"Unload Aux Volume"
		{UnloadVolume 1} 
	        tMenuGroup_AuxVolumeOptions }
	}}
	{ separator }
	{ command
	    "Load Main Surface..."
	    {DoFileDlog LoadMainSurface} }
	{ cascade "Load Surface Configuration..." {
	    { command
		"Original Vertices"
		{DoFileDlog LoadOriginalSurface}
		tMenuGroup_SurfaceLoading }
	    { command
		"Pial Vertices "
		{DoFileDlog LoadPialSurface}
		tMenuGroup_SurfaceLoading } 
	}}
	{ command
	    "Load Annotation"
	    {DoFileDlog LoadSurfaceAnnotation}
	    tMenuGroup_SurfaceLoading }
	{ command
	    "Unload Surface"
	    {UnloadSurface 0}
	    tMenuGroup_SurfaceLoading }
	{ command
	    "Write Surface Values..."
	    {DoFileDlog WriteSurfaceValues}
	    tMenuGroup_SurfaceLoading }
	{ cascade "Aux Surface" {
	    { command
		"Load Aux Main Surface..."
		{DoFileDlog LoadMainAuxSurface}
		tMenuGroup_SurfaceLoading }
	    { cascade "Load Aux Surface Configuration..." {
		{ command
		    "Original Vertices"
		    {DoFileDlog LoadOriginalAuxSurface}
		    tMenuGroup_SurfaceLoading }
		{ command
		    "Pial Vertices "
		    {DoFileDlog LoadPialAuxSurface}
		    tMenuGroup_SurfaceLoading }
	    }}
	    { command 
		"Load Aux Annotation"
		{DoFileDlog LoadAuxSurfaceAnnotation}
		tMenuGroup_SurfaceLoading }
	    { command
		"Unload Aux Surface"
		{UnloadSurface 1}
		tMenuGroup_SurfaceLoading }
	}}
	{ separator }
	{ command
	    "Load Overlay Data..."
	    {DoLoadFunctionalDlog overlay} }
	{ command
	    "Load Time Course Data..."
	    {DoLoadFunctionalDlog timecourse} }
	{ command
	    "Save Overlay Registration"
	    Overlay_SaveRegistration
	    tMenuGroup_Registration }
	{ separator }
	{ command
	    "New Segmentation..."
	    {DoFileDlog NewSegmentation} }
	{ command
	    "Load Segmentation..."
	    {DoFileDlog LoadSegmentation} }
	{ command
	    "Import Surface Annotation as Segmentation..."
	    {DoFileDlog ImportSurfaceAnnotation}
	    tMenuGroup_SurfaceLoading }
	{ command
	    "Save Segmentation"
	    "SaveSegmentationVolume 0"
	    tMenuGroup_SegmentationOptions }
	{ command
	    "Save Segmentation As..."
	    {DoFileDlog SaveSegmentationAs}
	    tMenuGroup_SegmentationOptions }
	{ command
	    "Save Changed Values"
	    {DoFileDlog ExportChangedSegmentationVolume}
	    tMenuGroup_SegmentationOptions }
	{ cascade "Aux Segmentation" {
	    { command
		"New Aux Segmentation..."
		{DoFileDlog NewAuxSegmentation} }
	    
	    { command
		"Load Aux Segmentation..."
		{DoFileDlog LoadAuxSegmentation} }
	    
	    { command
		"Save Aux Segmentation"
		"SaveSegmentationVolume 1"
		tMenuGroup_AuxSegmentationOptions }
	    
	    { command
		"Save Aux Segmentation As..."
		{DoFileDlog SaveAuxSegmentationAs}
		tMenuGroup_AuxSegmentationOptions }
	    { command
		"Save Aux Changed Values"
		{DoFileDlog ExportAuxChangedSegmentationVolume}
		tMenuGroup_SegmentationOptions }
	}}
	{ separator }
	{ cascade "Transforms" {
	    { command
		"Load Transform for Main Volume..."
		{DoFileDlog LoadVolumeDisplayTransform} }
	    { command
		"Load Transform for Aux Volume..."
		{DoFileDlog LoadAuxVolumeDisplayTransform} }
	    { command
		"Unload Transform for Main Volume"
		{UnloadVolumeDisplayTransform 0}
		tMenuGroup_VolumeMainTransformLoadedOptions }
	    { command
		"Unload Transform for Aux Volume"
		{UnloadVolumeDisplayTransform 1}
		tMenuGroup_VolumeAuxTransformLoadedOptions }
	}}
	{ cascade "Label" {
	    { command
		"Load Label..."
		{DoFileDlog LoadLabel} }
	    { command
		"Save Label As..."
		{DoFileDlog SaveLabelAs} }
	}}
	{ cascade "GCA" {
	    { command
		"Load GCA"
		{DoFileDlog LoadGCA} }
	    { command
		"Save GCA"
		{DoFileDlog SaveGCA}
		tMenuGroup_GCAOptions }
	    { command
		"Unload GCA"
		{UnloadGCA}
		tMenuGroup_GCAOptions }
	}}
	{ cascade "Head Points" {
	    { command
		"Load Head Points..."
		{DoFileDlog LoadHeadPts} }
	    { command
		"Save Head Point Transform"
		WriteHeadPointsTransform
		tMenuGroup_HeadPoints }
	    { command
		"Save Head Points"
		WriteHeadPointsFile
		tMenuGroup_HeadPoints }
	}}
	{ cascade "DTI" {
	    { command
		"Load DTI Volumes..."
		  {DoLoadDTIDlog} }
	}}
	{ cascade "GDF" {
	    { command
		"Load GDF..."
		  {DoLoadGDF} }
	}}
	{ command
	    "Save Control Points"
	    WriteControlPointFile
	    tMenuGroup_ControlPoints }
	{ separator }
	{ command
	    "Quit:Ctrl Q"
	    AllowSaveThenQuit } 
    }
   
    # edit menu 
    tkm_MakeMenu $mbwEdit "Edit" {
	{ command
	    "Undo Last Edit:Ctrl Z"
	    UndoLastEdit }
	{ separator }
	{ command
	    "Take Snapshot of Volume"
	    SnapshotVolume }
	{ command
	    "Restore Volume to Snapshot"
	    RestoreVolumeFromSnapshot }
	{ separator }
	{ command
	    "Clear Selection / Label"
	    "ClearSelection; RedrawAll" }
	{ command
	    "Clear Undo Volume"
	    ClearUndoVolume } 
    }

    # view menu
    tkm_MakeMenu $mbwView "View" {
	{ cascade "View Configurations" {
	    { radio 
		"Single View"
		"SetViewPreset $tViewPreset_None"
		gViewPreset
		0 }
	    { radio 
		"Multiple Orientations"
		"SetViewPreset $tViewPreset_Multiple"
		gViewPreset
		1 }
	    { radio 
		"Mosaic"
		"SetViewPreset $tViewPreset_Mosaic"
		gViewPreset
		2 } 
	}}
	{ separator }
	{ cascade "Tool Bars" {
	    { check
		"Main"
		"ShowToolBar main $gbShowToolBar(main)"
		gbShowToolBar(main) }
	    { check
		"Navigation"
		"ShowToolBar nav $gbShowToolBar(nav)"
		gbShowToolBar(nav) }
	    { check
		"Reconstruction"
		"ShowToolBar recon $gbShowToolBar(recon)"
		gbShowToolBar(recon) }
	}}
	{ cascade "Information" {
	    { check
		"Volume Index Coordinates"
		"ShowLabel kLabel_Coords_Vol $gbShowLabel(kLabel_Coords_Vol)"
		gbShowLabel(kLabel_Coords_Vol) }
	    { check
		"Volume RAS Coordinates"
		"ShowLabel kLabel_Coords_Vol_RAS $gbShowLabel(kLabel_Coords_Vol_RAS)"
		gbShowLabel(kLabel_Coords_Vol_RAS) }
	    { check
		"Volume Scanner Coordinates"
		"ShowLabel kLabel_Coords_Vol_Scanner $gbShowLabel(kLabel_Coords_Vol_Scanner)"
		gbShowLabel(kLabel_Coords_Vol_Scanner) }
	    { check
		"MNI Coordinates"
		"ShowLabel kLabel_Coords_Vol_MNI $gbShowLabel(kLabel_Coords_Vol_MNI)"
		gbShowLabel(kLabel_Coords_Vol_MNI) }
	    { check
		"Talairach Coordinates"
		"ShowLabel kLabel_Coords_Vol_Tal $gbShowLabel(kLabel_Coords_Vol_Tal)"
		gbShowLabel(kLabel_Coords_Vol_Tal) }
	    { check
		"Volume Value"
		"ShowLabel kLabel_Value_Vol $gbShowLabel(kLabel_Value_Vol)"
		gbShowLabel(kLabel_Value_Vol) }
	    { check
		"Aux Volume Value"
		"ShowLabel kLabel_Value_Aux $gbShowLabel(kLabel_Value_Aux)"
		gbShowLabel(kLabel_Value_Aux)
		tMenuGroup_AuxVolumeOptions }
	    { check
		"Functional Overlay Index Coordinates"
		"ShowLabel kLabel_Coords_Func $gbShowLabel(kLabel_Coords_Func)"
		gbShowLabel(kLabel_Coords_Func)
		tMenuGroup_OverlayOptions }
	    { check
		"Functional Overlay RAS Coordinates"
		"ShowLabel kLabel_Coords_Func_RAS $gbShowLabel(kLabel_Coords_Func_RAS)"
		gbShowLabel(kLabel_Coords_Func_RAS)
		tMenuGroup_OverlayOptions }
	    { check
		"Functional Overlay Value"
		"ShowLabel kLabel_Value_Func $gbShowLabel(kLabel_Value_Func)"
		gbShowLabel(kLabel_Value_Func)
		tMenuGroup_OverlayOptions }
	    { check
		"Segmentation Label"
		"ShowLabel kLabel_Label_SegLabel $gbShowLabel(kLabel_Label_SegLabel)"
		gbShowLabel(kLabel_Label_SegLabel)
		tMenuGroup_SegmentationOptions }
	    { check
		"Aux Segmentation Label"
		"ShowLabel kLabel_Label_AuxSegLabel $gbShowLabel(kLabel_Label_AuxSegLabel)"
		gbShowLabel(kLabel_Label_AuxSegLabel)
		tMenuGroup_AuxSegmentationOptions }
	    { check
		"Head Point Label"
		"ShowLabel kLabel_Label_Head $gbShowLabel(kLabel_Label_Head)"
		gbShowLabel(kLabel_Label_Head)
		tMenuGroup_HeadPoints }
	    { check
		"Surface Distance"
		"ShowLabel kLabel_SurfaceDistance $gbShowLabel(kLabel_SurfaceDistance)"
		gbShowLabel(kLabel_SurfaceDistance)
		tMenuGroup_SurfaceLoading} 
	    { check
		"Line Length"
		"ShowLabel kLabel_LineLength $gbShowLabel(kLabel_LineLength)"
		gbShowLabel(kLabel_LineLength)
		} 
	}}
	{ separator }
	{ cascade "Configure..." {
	    { command
		"Brightness / Contrast...:Ctrl B"
		DoVolumeColorScaleInfoDlog }
	    { command
		"Cursor..."
		DoCursorInfoDlog }
	    { command
		"Surface...:Ctrl U"
		DoSurfaceInfoDlog
		tMenuGroup_SurfaceViewing }
	    { command
		"Functional Overlay..."
		Overlay_DoConfigDlog
		tMenuGroup_OverlayOptions }
	    { command
		"Time Course Graph..."
		TimeCourse_DoConfigDlog
		tMenuGroup_TimeCourseOptions }
	    { command
		"Segmentation Display..."
		DoSegmentationVolumeDisplayInfoDlog
		tMenuGroup_SegmentationOptions } 
	    { command
		"DTI Display..."
		DoDTIVolumeDisplayInfoDlog
		tMenuGroup_DTIOptions } 
	}}
	{ cascade "Anatomical Sampling" {
	    { cascade "Main Volume" {
		{ radio 
		    "Nearest Neighbor"
		    "SendVolumeSampleType 0"
		    gVolume(0,sampleType)
		    0 }
		{ radio 
		    "Trilinear"
		    "SendVolumeSampleType 0"
		    gVolume(0,sampleType)
		    1 }
		{ radio 
		    "Sinc"
		    "SendVolumeSampleType 0"
		    gVolume(0,sampleType)
		    2 }
	    }}
	    { cascade "Aux Volume" {
		{ radio 
		    "Nearest Neighbor"
		    "SendVolumeSampleType 1"
		    gVolume(1,sampleType)
		    0 }
		{ radio 
		    "Trilinear"
		    "SendVolumeSampleType 1"
		    gVolume(1,sampleType)
		    1 }
		{ radio 
		    "Sinc"
		    "SendVolumeSampleType 1"
		    gVolume(1,sampleType)
		    2 }
	    }}
	}}
	{ cascade "Anatomical Resampling" {
	    { cascade "Main Volume" {
		{ radio 
		    "RAS"
		    "SendVolumeResampleMethod 0"
		    gVolume(0,resampleMethod)
		    0 }
		{ radio 
		    "Slice"
		    "SendVolumeResampleMethod 0"
		    gVolume(0,resampleMethod)
		    1 }
	    }}
	    { cascade "Aux Volume" {
		{ radio 
		    "RAS"
		    "SendVolumeResampleMethod 1"
		    gVolume(1,resampleMethod)
		    0 }
		{ radio 
		    "Slice"
		    "SendVolumeResampleMethod 1"
		    gVolume(1,resampleMethod)
		    1 }
	    }}
	}}
	{ separator }
	{ check 
	    "Anatomical Volume:Ctrl A"
	    "SendDisplayFlagValue flag_Anatomical"
	    gbDisplayFlag(flag_Anatomical) }
	{ radio 
	    "Main Volume:Ctrl 1"
	    "SendDisplayFlagValue flag_AuxVolume"
	    gbDisplayFlag(flag_AuxVolume) 
	    0 }
	{ radio
	    "Aux Volume:Ctrl 2"
	    "SendDisplayFlagValue flag_AuxVolume"
	    gbDisplayFlag(flag_AuxVolume)
	    1
	    tMenuGroup_AuxVolumeOptions }
	{ separator }
	{ check
	    "Maximum Intensity Projection"
	    "SendDisplayFlagValue flag_MaxIntProj"
	    gbDisplayFlag(flag_MaxIntProj) 
	}
	{ check
	    "Main Surface:Ctrl M"
	    "SendDisplayFlagValue flag_MainSurface"
	    gbDisplayFlag(flag_MainSurface) 
	    tMenuGroup_SurfaceViewing }
	{ check
	    "Original Surface:Ctrl O"
	    "SendDisplayFlagValue flag_OriginalSurface"
	    gbDisplayFlag(flag_OriginalSurface) 
	    tMenuGroup_OriginalSurfaceViewing }
	{ check
	    "Pial Surface:Ctrl P"
	    "SendDisplayFlagValue flag_PialSurface"
	    gbDisplayFlag(flag_PialSurface) 
	    tMenuGroup_PialSurfaceViewing }
	{ check
	    "Surface Vertices:Ctrl V"
		"SendDisplayFlagValue flag_DisplaySurfaceVertices"
	    gbDisplayFlag(flag_DisplaySurfaceVertices) 
	    tMenuGroup_SurfaceViewing }
	{ check
	    "Interpolate Surface Vertices:Ctrl I"
	    "SendDisplayFlagValue flag_InterpolateSurfaceVertices"
	    gbDisplayFlag(flag_InterpolateSurfaceVertices) 
	    tMenuGroup_SurfaceViewing }
	{ check
	    "Functional Overlay:Ctrl F"
	    "SendDisplayFlagValue flag_FunctionalOverlay"
	    gbDisplayFlag(flag_FunctionalOverlay) 
	    tMenuGroup_OverlayOptions }
	{ check
	    "Functional Color Scale Bar"
	    "SendDisplayFlagValue flag_FunctionalColorScaleBar"
	    gbDisplayFlag(flag_FunctionalColorScaleBar) 
	    tMenuGroup_OverlayOptions }
	{ check
	    "Mask to Functional Overlay"
	    "SendDisplayFlagValue flag_MaskToFunctionalOverlay"
	    gbDisplayFlag(flag_MaskToFunctionalOverlay) 
	    tMenuGroup_OverlayOptions }
	{ check
	    "Mask Functional Overlay to Aux Volume"
	    "SendDisplayFlagValue flag_MaskFunctionalOverlayToAux"
	    gbDisplayFlag(flag_MaskFunctionalOverlayToAux) 
	    tMenuGroup_OverlayOptions }
	{ check
	    "Show Histogram Percent Change"
	    "SendDisplayFlagValue flag_HistogramPercentChange"
	    gbDisplayFlag(flag_HistogramPercentChange) 
	    tMenuGroup_VLIOptions }
	{ check
	    "Segmentation Overlay:Ctrl G"
	    "SendDisplayFlagValue flag_SegmentationVolumeOverlay"
	    gbDisplayFlag(flag_SegmentationVolumeOverlay) 
	    tMenuGroup_SegmentationOptions }
	{ check
	    "Aux Segmentation Overlay:Alt G"
	    "SendDisplayFlagValue flag_AuxSegmentationVolumeOverlay"
	    gbDisplayFlag(flag_AuxSegmentationVolumeOverlay) 
	    tMenuGroup_AuxSegmentationOptions }
	{ check
	    "Segmentation Label Volume Count"
	    "SendDisplayFlagValue flag_SegLabelVolumeCount"
	    gbDisplayFlag(flag_SegLabelVolumeCount) 
	    tMenuGroup_SegmentationOptions }
	{ check
	    "DTI Overlay:Alt D"
	    "SendDisplayFlagValue flag_DTIOverlay"
	    gbDisplayFlag(flag_DTIOverlay) 
	    tMenuGroup_DTIOptions }
	{ check
	    "Selection / Label:Ctrl S"
	    "SendDisplayFlagValue flag_Selection"
	    gbDisplayFlag(flag_Selection) }
	{ check
	    "Head Points"
	    "SendDisplayFlagValue flag_HeadPoints"
	    gbDisplayFlag(flag_HeadPoints) 
	    tMenuGroup_HeadPoints }
	{ check
	    "Control Points:Ctrl T"
	    "SendDisplayFlagValue flag_ControlPoints"
	    gbDisplayFlag(flag_ControlPoints) }
	{ check
	    "Cursor:Ctrl C"
	    "SendDisplayFlagValue flag_Cursor"
	    gbDisplayFlag(flag_Cursor) }
	{ check
	    "Axes"
	    "SendDisplayFlagValue flag_Axes"
	    gbDisplayFlag(flag_Axes) }
	{ check
	    "Edited Voxels"
	    "SendDisplayFlagValue flag_UndoVolume"
	    gbDisplayFlag(flag_UndoVolume) }
    }

    # tools menu
    tkm_MakeMenu $mbwTools "Tools" {
	{ radio
	    "Navigate:N"
	    "SetTool $DspA_tTool_Navigate"
	    gTool
	    0 }
	{ radio
	    "Select Voxels:S"
	    "SetTool $DspA_tTool_Select"
	    gTool
	    1 }
	{ radio
	    "Edit Voxels:A"
	    "SetTool $DspA_tTool_Edit"
	    gTool
	    2 }
	{ radio
	    "Edit Segmentation:G"
	    "SetTool $DspA_tTool_EditSeg"
	    gTool
	    3 }
	{ radio
	    "Edit Ctrl Pts:T"
	    "SetTool $DspA_tTool_CtrlPts"
	    gTool
	    4 }
	{ radio
	    "Line:L"
	    "SetTool $DspA_tTool_Line"
	    gTool
	    5 }
	{ separator }
	{ command
	    "Configure Brush Info..."
	    DoBrushInfoDlog }
	{ command
	    "Configure Volume Brush..."
	    DoEditBrushInfoDlog }
	{ command
	    "Configure Segmentation Brush..."
	    DoEditSegBrushInfoDlog
	    tMenuGroup_SegmentationOptions }
	{ command
	    "Configure Flood Select..."
	    DoEditFloodSelectParamsDlog }
	{ separator }
	{ command
	    "Save Point"
	    SendCursor }
	{ command
	    "Goto Saved Point"
	    ReadCursor }
	{ separator }
	{ command
	    "Goto Point..."
	    DoGotoPointDlog }
	{ command
	    "Goto Center of Selection"
	    SetCursorToCenterOfSelection }
	{ separator }
	{ command
	    "Add Line to Selection"
	    AddLineToSelection }
	{ command
	    "Write Line to Label..."
	    {DoFileDlog WriteLineLabel} }
	{ command
	    "Label Writer Helper..."
	    "DoLabelWriterHelperDlog" }
	{ separator }
	{ cascade "Volume" {
	    { command
		"Threshold Volume..."
		DoThresholdDlog }
	    { command
		"Flip Volume..."
		DoFlipVolumeDlog }
	    { command
		"Rotate Volume..."
		DoRotateVolumeDlog }
	    { command
		"Smart Cut"
		SmartCutAtCursor } } }
	{ cascade "Surface" {
	    { command 
		"Show Nearest Main Vertex:Alt M"
		ShowNearestMainVertex
		tMenuGroup_SurfaceViewing }
	    { command
		"Show Nearest Original Vertex:Alt O"
		ShowNearestOriginalVertex
		tMenuGroup_OriginalSurfaceViewing }
	    { command
		"Show Nearest Pial Vertex:Alt P"
		ShowNearestPialVertex
		tMenuGroup_PialSurfaceViewing }
	    { separator }
	    { command 
		"Show Nearest Main Surface Edge"
		ShowNearestInterpolatedMainVertex
		tMenuGroup_SurfaceViewing }
	    { command
		"Show Nearest Original Surface Edge"
		ShowNearestInterpolatedOriginalVertex
		tMenuGroup_OriginalSurfaceViewing }
	    { command
		"Show Nearest Pial Surface Edge"
		ShowNearestInterpolatedPialVertex
		tMenuGroup_PialSurfaceViewing }
	    { separator }
	    { command
		"Find Main Vertex..."
		{ DoFindVertexDlog $Surf_tVertexSet(main) } 
		tMenuGroup_SurfaceViewing }
	    { command
		"Find Original Vertex..."
		{ DoFindVertexDlog $Surf_tVertexSet(original) }
		tMenuGroup_OriginalSurfaceViewing }
	    { command
		"Find Pial Vertex..."
		{ DoFindVertexDlog $Surf_tVertexSet(pial) }
		tMenuGroup_PialSurfaceViewing }
	    { separator }
	    { command
		"Set Vertex Distance at Cursor"
		{ SetSurfaceDistanceAtCursor }
		tMenuGroup_SurfaceLoading }
	    { cascade "Set MRI Value in Surface..." {
		{ command
		    "Set MRI Value in Surface at Main Vertex Set"
		    { SetMRIValueAtCursorInSurface $Surf_tVertexSet(main)}
		    tMenuGroup_SurfaceLoading }
		{ command
		    "Set MRI Value in Surface at Original Vertex Set"
		    { SetMRIValueAtCursorInSurface $Surf_tVertexSet(original)}
		    tMenuGroup_SurfaceLoading }
		{ command
		    "Set MRI Value in Surface at Pial Vertex Set"
		    { SetMRIValueAtCursorInSurface $Surf_tVertexSet(pial)}
		    tMenuGroup_SurfaceLoading }
	    }}
	    { command
		"Average Vertex Positions..."
		{ DoAverageSurfaceVertexPositionsDlog }
		tMenuGroup_SurfaceViewing } 
	}}

	{ cascade "fMRI" {
	    { command
		"Select Contiguous Voxels by Func Value"
		{ SelectVoxelsByFuncValue 0 }
		tMenuGroup_OverlayOptions }
	    { command
		"Select Contiguous Voxels by Threshold"
		{ SelectVoxelsByFuncValue 2 }
		tMenuGroup_OverlayOptions }
	    { command
		"Select Functional Voxel"
		{ SelectVoxelsByFuncValue 1 }
		tMenuGroup_OverlayOptions }
	    { command
		"Register Functional Overlay..."
		{ DoRegisterOverlayDlog }
		tMenuGroup_Registration }
	    { command
		"Restore Overlay Registration"
		{ Overlay_RestoreRegistration }
		tMenuGroup_Registration }
	    { command
		"Set Registration to Identity"
		{ Overlay_SetRegistrationToIdentity }
		tMenuGroup_Registration }
	    { command
		"Graph Current Selection"
		{ GraphSelectedRegion }
		tMenuGroup_TimeCourseOptions }
	    { command
		"Print Time Course Summary to File.."
		{ DoFileDlog PrintTimeCourse }
		tMenuGroup_TimeCourseOptions }
	    { command
		"Save Time Course Graph to Postscript File.."
		{ DoFileDlog SaveTimeCourseToPS }
		tMenuGroup_TimeCourseOptions }
	}}
	{ cascade "Segmentation" {
	    { command
		"Recompute Segmentation"
		DoRecomputeSegmentation
		tMenuGroup_GCAOptions }
	    { command
		"Select Main Segmentation Label At Cursor"
		"SelectSegLabelAtCursor $tkm_tSegType(main)"
		tMenuGroup_SegmentationOptions }
	    { command
		"Select Aux Segmentation Label At Cursor"
		"SelectSegLabelAtCursor $tkm_tSegType(aux)"
		tMenuGroup_SegmentationOptions }
	    { check
		"Verbose GCA Display"
		"SendDisplayFlagValue flag_VerboseGCADump"
		gbDisplayFlag(flag_VerboseGCADump) 
		tMenuGroup_GCAOptions }
	}}
	{ cascade "Head Points" {
	    { command
		"Restore Head Points"
		RestoreHeadPts
		tMenuGroup_HeadPoints }
	    { command
		"Edit Current Head Point Label..."
		DoEditHeadPointLabelDlog
		tMenuGroup_HeadPoints }
	    { command
		"Register Head Points..."
		DoRegisterHeadPtsDlog
		tMenuGroup_HeadPoints }
	}}
	{ cascade "Group" {
	    { command "Save Plotted Data to Table"
		{ DoFileDlog SaveGDFPlotToTable }
		tMenuGroup_GDFLoaded }
	    { command "Save Plot to Postscript File"
		{ DoFileDlog SaveGDFPlotToPS }
		tMenuGroup_GDFLoaded }
	}}
	{ command
	    "Save RGB..."
	    {DoFileDlog SaveRGB} }
	{ command
	    "Save RGB Series..."
	    DoSaveRGBSeriesDlog }
	{ command
	    "Save TIFF..."
	    {DoFileDlog SaveTIFF} }
	{ command
	    "Save TIFF Series..."
	    DoSaveTIFFSeriesDlog }
    }

    # Save the location of the Control Points tool menu item.
    set gControlPointsMenuCommand(menu) $mbwTools.mw
    set gControlPointsMenuCommand(entry) 5
	
    pack $mbwFile $mbwEdit $mbwView $mbwTools \
      -side left
}

proc CreateCursorFrame { ifwTop } {

    global gbLinkedCursor 

    set fwLabel             $ifwTop.fwMainLabel
    set fwLinkCheckbox      $ifwTop.fwLinkCheckbox
    set fwLabels            $ifwTop.fwLabels

    frame $ifwTop

    # the label that goes at the top of the frame
    tkm_MakeBigLabel $fwLabel "Cursor"

    # make the labels
    CreateLabelFrame $fwLabels cursor

    # pack the subframes in a column. 
    pack $fwLabel $fwLabels \
      -side top             \
      -anchor w

}

proc CreateMouseoverFrame { ifwTop } {

    set fwLabel             $ifwTop.fwMainLabel
    set fwLabels            $ifwTop.fwLabels

    frame $ifwTop

    # the label that goes at the top of the frame
    tkm_MakeBigLabel $fwLabel "Mouse"

    # make the labels
    CreateLabelFrame $fwLabels mouseover

    # pack the subframes in a column. 
    pack $fwLabel $fwLabels \
      -side top             \
      -anchor w
}

proc CreateLabelFrame { ifwTop iSet } {

    global glLabel gfwaLabel gsaLabelContents
    global mri_tCoordSpace_VolumeIdx
    global mri_tCoordSpace_SurfaceRAS
    global mri_tCoordSpace_RAS
    global mri_tCoordSpace_Talairach

    frame $ifwTop

    # create the frame names
    foreach label $glLabel {
	set gfwaLabel($label,$iSet) $ifwTop.fw$label
    }
    
    # create two active labels in each label frame
    foreach label $glLabel {
	frame $gfwaLabel($label,$iSet)
	set fwLabel $gfwaLabel($label,$iSet).fwLabel
	set fwValue $gfwaLabel($label,$iSet).fwValue

	tkm_MakeActiveLabel $fwLabel "" gsaLabelContents($label,name) 14

	if { $label == "kLabel_Coords_Vol" && $iSet == "cursor" } {
	    tkm_MakeEntry $fwValue "" gsaLabelContents($label,value,$iSet) 18 \
		"SetCursorFromLabelContents $label $iSet $mri_tCoordSpace_VolumeIdx"
	} elseif { $label == "kLabel_Coords_Vol_RAS" && $iSet == "cursor" } {
	    tkm_MakeEntry $fwValue "" gsaLabelContents($label,value,$iSet) 18 \
		"SetCursorFromLabelContents $label $iSet $mri_tCoordSpace_SurfaceRAS"
	} elseif { $label == "kLabel_Coords_Vol_Scanner" && $iSet == "cursor" } {
	    tkm_MakeEntry $fwValue "" gsaLabelContents($label,value,$iSet) 18 \
		"SetCursorFromLabelContents $label $iSet $mri_tCoordSpace_RAS"
	} elseif { $label == "kLabel_Coords_Vol_Tal" && $iSet == "cursor" } {
	    tkm_MakeEntry $fwValue "" gsaLabelContents($label,value,$iSet) 18 \
		"SetCursorFromLabelContents $label $iSet $mri_tCoordSpace_Talairach"
	} else {
	    tkm_MakeActiveLabel $fwValue "" gsaLabelContents($label,value,$iSet) 18
	}

	pack $fwLabel $fwValue \
	    -side left \
	    -anchor w
    }
    
    ShowLabel kLabel_Coords_Vol_RAS 1
    ShowLabel kLabel_Coords_Vol_Tal 1
    ShowLabel kLabel_Value_Vol 1
}

proc SetCursorFromLabelContents { iLabel iSet iSpace } {
    global gsaLabelContents
    
    # Get the input string.
    set sCoords $gsaLabelContents($iLabel,value,$iSet)

    # [-+]? matches the leading - or +
    # \d+ matches a series of digits like 12
    # \d+\.\d+ matches floating point numbers like 12.34
    set sFiltered [regexp -inline -all -- {[-+]?\d+|[-+]?\d+\.\d+} $sCoords]

    # Make sure we have three elements.
    if { [llength $sFiltered] != 3 } {
	ErrorDlog {Invalid coordinate string. Make sure there are three numbers.}
	return;
    }

    # Set the cursor.
    SetCursor $iSpace \
	[lindex $sFiltered 0] [lindex $sFiltered 1] [lindex $sFiltered 2]
}

proc CreateToolBar { ifwToolBar } {
    
    global gaLoadMenu
    global gfwaToolBar
    global gTool gViewPreset gDisplayedVolume
    global gbLinkedCursor
    global gBrushInfo DspA_tBrushShape_Square DspA_tBrushShape_Circle
    global gOrientation gnZoomLevel gnVolSlice
    global glActiveFlags gnFlagIndex

    frame $ifwToolBar

    # main toolbar
    set gfwaToolBar(main)  $ifwToolBar.fwMainBar
    set fwTools            $gfwaToolBar(main).fwTools
    set fwSurfaces         $gfwaToolBar(main).fwSurfaces
    set fwVolumeToggles    $gfwaToolBar(main).fwVolumeToggles
    set fwSegVolume        $gfwaToolBar(main).fwSegVolume
    set fwOverlayVolume    $gfwaToolBar(main).fwOverlayVolume
    set fwScreenShot       $gfwaToolBar(main).fwScreenShot

    frame $gfwaToolBar(main) -border 2 -relief raised
    
    tkm_MakeToolbar $fwTools \
      1 \
      gTool \
      UpdateToolWrapper {
	  { image 0 icon_navigate "Navigate Tool (n)" }
	  { image 1 icon_edit_label "Select Voxels Tool (s)" }
	  { image 2 icon_edit_volume "Edit Voxels Tool (a)" }
	  { image 3 icon_edit_parc "Edit Segmentation Tool (g)" }
	  { image 4 icon_edit_ctrlpts "Edit Ctrl Pts Tool (t)" }
	  { image 5 icon_line_tool "Line Tool (l)" } 
      }

    # Control-left-click to show settings for the tool.
    bind [$fwTools.tbw subwidget 2] <Control-Button-1> {
	DoBrushInfoDlog
	DoEditBrushInfoDlog
    }
    bind [$fwTools.tbw subwidget 3] <Control-Button-1> {
	DoBrushInfoDlog
	DoEditSegBrushInfoDlog
    }
    
    tkm_MakeCheckboxes $fwSurfaces h {
	{ image icon_surface_main gbDisplayFlag(flag_MainSurface)
	    "SendDisplayFlagValue flag_MainSurface" "Show Main Surface" }
	{ image icon_surface_original gbDisplayFlag(flag_OriginalSurface)
	    "SendDisplayFlagValue flag_OriginalSurface"
	    "Show Original Surface" }
	{ image icon_surface_pial gbDisplayFlag(flag_PialSurface)
	    "SendDisplayFlagValue flag_PialSurface"
	    "Show Pial Surface" } 
    }
    tkm_AddCheckboxToEnableGroup tMenuGroup_SurfaceViewing $fwSurfaces.cb0
    tkm_AddCheckboxToEnableGroup tMenuGroup_OriginalSurfaceViewing \
	$fwSurfaces.cb1
    tkm_AddCheckboxToEnableGroup tMenuGroup_PialSurfaceViewing $fwSurfaces.cb2
    # Control-left-click to show the settings, control-right-click to
    # open a load dialog box. Note that for the control button 1
    # click, we have to toggle flag, since even control-clicking will
    # trigger the main cmd, so we need to reverse it.
    bind $fwSurfaces.cb0 <Control-Button-1> {
	if { [tkm_IsGroupEnabled tMenuGroup_SurfaceViewing] } {
	    set gbDisplayFlag(flag_MainSurface) \
		[expr !$gbDisplayFlag(flag_MainSurface)] 
	    DoSurfaceInfoDlog
	}
    }
    bind $fwSurfaces.cb1 <Control-Button-1> {
	if { [tkm_IsGroupEnabled tMenuGroup_OriginalSurfaceViewing] } {
	    set gbDisplayFlag(flag_OriginalSurface) \
		[expr !$gbDisplayFlag(flag_OriginalSurface)] 
	    DoSurfaceInfoDlog
	}
    }
    bind $fwSurfaces.cb2 <Control-Button-1> {
	if { [tkm_IsGroupEnabled tMenuGroup_PialSurfaceViewing] } {
	    set gbDisplayFlag(flag_PialSurface) \
		[expr !$gbDisplayFlag(flag_PialSurface)] 
	    DoSurfaceInfoDlog
	}
    }
    bind $fwSurfaces.cb0 <Control-Button-3> "DoFileDlog LoadMainSurface"
    bind $fwSurfaces.cb1 <Control-Button-3> "DoFileDlog LoadOriginalSurface"
    bind $fwSurfaces.cb2 <Control-Button-3> "DoFileDlog LoadPialSurface"

    # Volume buttons. Main button, you can always control left click
    # to bring up options and control right click to see load
    # dialog. Aux button, only enables when aux volume is enabled, but
    # has same behavior.
    tkm_MakeToolbar $fwVolumeToggles \
	1 \
	gbDisplayFlag(flag_AuxVolume) \
	UpdateVolumeToggleWrapper {
	    { image 0 icon_main_volume "Show Main Volume" }
	    { image 1 icon_aux_volume "Show Aux Volume" } 
	}
    bind [$fwVolumeToggles.tbw subwidget 0] <Control-Button-1> \
	DoVolumeColorScaleInfoDlog
    bind [$fwVolumeToggles.tbw subwidget 0] <Control-Button-3> \
	"DoFileDlog LoadVolume"
    tkm_AddCheckboxToEnableGroup tMenuGroup_AuxVolumeOptions \
	[$fwVolumeToggles.tbw subwidget 1]
    bind [$fwVolumeToggles.tbw subwidget 1] <Control-Button-3> \
	"DoFileDlog LoadAuxVolume"
    bind [$fwVolumeToggles.tbw subwidget 1] <Control-Button-1> {
	if { [tkm_IsGroupEnabled tMenuGroup_AuxVolumeOptions] } {
	    set gbDisplayFlag(flag_AuxVolume) \
		[expr !$gbDisplayFlag(flag_AuxVolume)] 
	    DoVolumeColorScaleInfoDlog
	}
    }

    # Make and attach the load menus to the buttons. These will be
    # build in BuildVolumeLoadList when the subject dir is set.
    set gaLoadMenu(main) [menu .mmLoadMenu]
    bind [$fwVolumeToggles.tbw subwidget 0] <Button-3> \
	"tk_popup $gaLoadMenu(main) %X %Y"
    set gaLoadMenu(aux) [menu .mmAuxLoadMenu]
    bind [$fwVolumeToggles.tbw subwidget 1] <Button-3> \
	"tk_popup $gaLoadMenu(aux) %X %Y"


    # Segmentation checkbox. Control-left click shows options,
    # control-right click shows load dlog.
    tkm_MakeCheckboxes $fwSegVolume h {
	{ image icon_seg_volume gbDisplayFlag(flag_SegmentationVolumeOverlay)
	    "SendDisplayFlagValue flag_SegmentationVolumeOverlay" 
	    "Show Segmentation" }
    }
    tkm_AddCheckboxToEnableGroup tMenuGroup_SegmentationOptions \
	$fwSegVolume.cb0
    bind $fwSegVolume.cb0 <Control-Button-1> {
	if { [tkm_IsGroupEnabled tMenuGroup_SegmentationOptions] } {
	    set gbDisplayFlag(flag_SegmentationVolumeOverlay) \
		[expr !$gbDisplayFlag(flag_SegmentationVolumeOverlay)] 
	    DoSegmentationVolumeDisplayInfoDlog
	}
    }
    bind $fwSegVolume.cb0 <Control-Button-3> "DoFileDlog LoadSegmentation"
    
    # Overlay checkbox. Control-left click shows options,
    # control-right click shows load dlog.
    tkm_MakeCheckboxes $fwOverlayVolume h {
	{ image icon_overlay gbDisplayFlag(flag_FunctionalOverlay)
	    "SendDisplayFlagValue flag_FunctionalOverlay" 
	    "Show Functional Overlay" }
	{ image icon_color_scalebar gbDisplayFlag(flag_FunctionalColorScaleBar)
	    "SendDisplayFlagValue flag_FunctionalColorScaleBar" 
	    "Show Functional Color Bar" }
    }
    tkm_AddCheckboxToEnableGroup tMenuGroup_OverlayOptions \
	$fwOverlayVolume.cb0
    tkm_AddCheckboxToEnableGroup tMenuGroup_OverlayOptions \
	$fwOverlayVolume.cb1
    bind $fwOverlayVolume.cb0 <Control-Button-1> {
	if { [tkm_IsGroupEnabled tMenuGroup_OverlayOptions] } {
	    set gbDisplayFlag(flag_FunctionalOverlay) \
		[expr !$gbDisplayFlag(flag_FunctionalOverlay)] 
	    Overlay_DoConfigDlog
	}
    }
    bind $fwOverlayVolume.cb0 <Control-Button-3> "DoLoadFunctionalDlog overlay"

    tkm_MakeButtons $fwScreenShot { 
	{ image icon_camera { DoFileDlog SaveTIFF } "Save TIFF" } }
    bind $fwScreenShot.bw0 <Control-Button> { DoFileDlog SaveRGB; break }
    
    pack $fwTools $fwSurfaces $fwVolumeToggles \
	$fwSegVolume $fwOverlayVolume $fwScreenShot \
	-side left \
	-anchor w \
	-padx 5
    
    # navigation toolbar
    set gfwaToolBar(nav)   $ifwToolBar.fwNavBar
    set fwOrientation      $gfwaToolBar(nav).fwOrientation
    set fwViews            $gfwaToolBar(nav).fwViews
    set fwCurSlice         $gfwaToolBar(nav).fwCurSlice
    set fwZoomButtons      $gfwaToolBar(nav).fwZoomButtons
    set fwZoomLevel        $gfwaToolBar(nav).fwZoomLevel
    set fwHome             $gfwaToolBar(nav).fwHome
    set fwPoint            $gfwaToolBar(nav).fwPoint
    set fwLinkedCursor     $gfwaToolBar(nav).fwLinkedCursor

    frame $gfwaToolBar(nav) -border 2 -relief raised

    if { 1 } {
	tkm_MakeToolbar $fwOrientation \
	    1 \
	    gOrientation \
	    UpdateOrientationWrapper {
		{ image 0 icon_orientation_coronal "Coronal View" }
		{ image 1 icon_orientation_horizontal "Horizontal View" }
		{ image 2 icon_orientation_sagittal "Sagittal View" } }
    } else {
	menubutton $fwOrientation \
	    -relief raised -border 2 \
	    -image icon_orientation_coronal_pdm \
	    -menu $fwOrientation.mw
	menu $fwOrientation.mw
	$fwOrientation.mw add command \
	    -image icon_orientation_coronal \
	    -command "SetOrientation 0; $fwOrientation config -image icon_orientation_coronal_pdm"
	$fwOrientation.mw add command \
	    -image icon_orientation_horizontal \
	    -command "SetOrientation 1; $fwOrientation config -image icon_orientation_horizontal_pdm"
	$fwOrientation.mw add command \
	    -image icon_orientation_sagittal \
	    -command "SetOrientation 2; $fwOrientation config -image icon_orientation_sagittal_pdm"
    }
    
    if { 0 } {
	tkm_MakeToolbar $fwViews \
	    1 \
	    gViewPreset \
	    UpdateViewPresetWrapper {
		{ image 0 icon_view_single "Single View" } 
		{ image 1 icon_view_multiple "Multiple Views" } 
		{ image 2 icon_view_mosaic "Mosaic" } 
	    }
    } else {
	menubutton $fwViews \
	    -relief raised -border 2 \
	    -image icon_view_single_pdm \
	    -menu $fwViews.mw
	menu $fwViews.mw
	$fwViews.mw add command \
	    -image icon_view_single \
	    -command "SetViewPreset 0; $fwViews config -image icon_view_single_pdm"
	$fwViews.mw add command \
	    -image icon_view_multiple \
	    -command "SetViewPreset 1; $fwViews config -image icon_view_multiple_pdm"
	$fwViews.mw add command \
	    -image icon_view_mosaic \
	    -command "SetViewPreset 2; $fwViews config -image icon_view_mosaic_pdm"
    }

    tkm_MakeEntryWithIncDecButtons $fwCurSlice "Slice" gnVolSlice \
      { SetSlice $gnVolSlice } 1

    tkm_MakeButtons $fwZoomButtons { \
      { image icon_zoom_out \
      { SetZoomLevelWrapper [expr $gnZoomLevel / 2] } "Zoom Out" } \
      { image icon_zoom_in \
      { SetZoomLevelWrapper [expr $gnZoomLevel * 2] } "Zoom In" } }
    
    tkm_MakeEntryWithIncDecButtons $fwZoomLevel "Zoom" gnZoomLevel \
      { SetZoomLevelWrapper } 1
    
    tkm_MakeButtons $fwHome { \
      { image icon_home {SetCursorToCenterOfVolume 0} "Restore View to Home" } }

    tkm_MakeButtons $fwPoint { \
      { image icon_cursor_save {SendCursor} "Save Point" } \
      { image icon_cursor_goto {ReadCursor} "Goto Saved Point" } }

    tkm_MakeCheckboxes $fwLinkedCursor h {
	{ image icon_linked_cursors gbLinkedCursor \
	      "SendLinkedCursorValue" "Link Cursors" } }
    
    pack $fwOrientation $fwViews $fwCurSlice $fwZoomButtons \
	$fwZoomLevel $fwHome $fwPoint $fwLinkedCursor \
	-side left \
	-anchor w \
	-padx 5
    
    # recon toolbar
    set gfwaToolBar(recon) $ifwToolBar.fwBrushBar
    set fwShape            $gfwaToolBar(recon).fwShape
    set fw3D               $gfwaToolBar(recon).fw3D
    set fwRadius           $gfwaToolBar(recon).fwRadius
    set fwSnapshot         $gfwaToolBar(recon).fwSnapshot

    frame $gfwaToolBar(recon) -border 2 -relief raised

    tkm_MakeToolbar $fwShape \
      1 \
      gBrushInfo(shape) \
      UpdateShapeWrapper { \
      { image 0 icon_brush_circle "Circular Brush" } \
      { image 1 icon_brush_square "Square Brush" } }

    tkm_MakeCheckboxes $fw3D h { \
      { image icon_brush_3d gBrushInfo(3d) \
      "SetBrushConfiguration" "Activate 3D Brush" } }

    tkm_MakeEntryWithIncDecButtons $fwRadius "Radius" gBrushInfo(radius) \
      { UpdateBrushConfigurationWrapper } 1

    tkm_MakeButtons $fwSnapshot { \
      { image icon_snapshot_save \
      { SnapshotVolume } "Take Snapshot of Volume" } \
      { image icon_snapshot_load \
      { RestoreVolumeFromSnapshot } "Restore Volume from Snapshot" } }

    pack $fwShape $fw3D $fwRadius $fwSnapshot \
      -side left \
      -anchor w \
      -padx 5
}

proc ShowToolBar { isWhich ibShow } {

    global gfwaToolBar gbShowToolBar

    if { $ibShow == 1 } {   

  if { [catch { pack $gfwaToolBar($isWhich) \
    -side top \
    -fill x \
    -expand yes \
    -after $gfwaToolBar(main) } sResult] == 1 } {
      
      pack $gfwaToolBar($isWhich) \
        -side top \
        -fill x \
        -expand yes
  }
  
    } else {

  pack forget $gfwaToolBar($isWhich)
    }

    set gbShowToolBar($isWhich) $ibShow
}


proc CreateImages {} {

    global ksImageDir

    foreach image_name { icon_edit_label icon_edit_volume 
	icon_navigate icon_edit_ctrlpts icon_edit_parc icon_line_tool
	icon_view_single icon_view_multiple icon_view_mosaic
	icon_view_single_pdm icon_view_multiple_pdm icon_view_mosaic_pdm
	icon_overlay icon_color_scalebar
	icon_cursor_goto icon_cursor_save
	icon_main_volume icon_aux_volume icon_seg_volume icon_linked_cursors
	icon_arrow_up icon_arrow_down icon_arrow_left icon_arrow_right
	icon_arrow_cw icon_arrow_ccw
	icon_arrow_expand_x icon_arrow_expand_y
	icon_arrow_shrink_x icon_arrow_shrink_y
	icon_orientation_coronal icon_orientation_horizontal
	icon_orientation_sagittal
	icon_orientation_coronal_pdm icon_orientation_horizontal_pdm
	icon_orientation_sagittal_pdm
	icon_zoom_in icon_zoom_out
	icon_brush_square icon_brush_circle icon_brush_3d
	icon_surface_main icon_surface_original icon_surface_pial 
	icon_snapshot_save icon_snapshot_load
	icon_marker_crosshair icon_marker_diamond
	icon_stopwatch icon_camera icon_home } {
	
	if { [catch {image create photo $image_name -file \
			 [file join $ksImageDir $image_name.gif]} \
		  sResult] != 0 } {
	    dputs "Error loading $image_name:"
	    dputs $sResult
	}
    }
}

proc UpdateVolumeToggleWrapper { iValue ibStatus } {
    SendDisplayFlagValue flag_AuxVolume
}

proc UpdateOrientationWrapper { iValue ibStatus } {
    if { $ibStatus == 1 } {
  SetOrientation $iValue
    }
}

proc UpdateShapeWrapper { iValue ibStatus } {
    if { $ibStatus == 1 } {
  SetBrushConfiguration
    }
}

proc UpdateBrushConfigurationWrapper { iValue } {
    SetBrushConfiguration
}

proc UpdateCursorInfoWrapper { iValue ibStatus } {
    global gCursor
    if { $ibStatus == 1 } {
	gCursor(shape) = $iValue
	SendCursorConfiguration
    }   
}

proc UpdateToolWrapper { iValue ibStatus } {
    if { $ibStatus == 1 } {
	SetTool $iValue
    }
}

proc UpdateViewPresetWrapper { iValue ibStatus } {
    if { $ibStatus == 1 } {
	SetViewPreset $iValue
    }
}

proc UpdateLinkedCursorWrapper { iValue ibStatus } {
    SetLinkedCursorFlag $iValue
}

proc SetZoomLevelWrapper { inLevel } {
    global gnZoomLevel
    set gnZoomLevel $inLevel
    SetZoomLevel $gnZoomLevel
}


# ================================================================== BAR CHART


# works just by passing it the following arguments:
#   -title <string> : the title of the window and graph
#   -xAxisTitle <string> : the title of the x axis
#   -yAxisTitle <string> : the title of the x axis
#   -label1 <string> : the label in the legend for the first element
#   -label2 <string> : the label in the legend for the second element
#   -values1 <list> : list of values for the first element
#   -values2 <list> : list of values for the second element
#   -xAxisLabels <list> : the list of labels for the x axis
# note that the number of elements in -values1, -values2, and -xAxisLabels
# should be the same.
proc BarChart_Draw { args } {

    global glsXAxisLabels
    global kNormalFont

    # default values
    set tArgs(-title) ""
    set tArgs(-xAxisTitle) ""
    set tArgs(-yAxisTitle) ""
    set tArgs(-xAxisLabels) ""
    set tArgs(-label1) ""
    set tArgs(-label2) ""
    set tArgs(-values1) ""
    set tArgs(-values2) ""

    # get the params
    array set tArgs $args
    
    # find an unused window name.
    set nSuffix 0
    set wwChartWindow .bcw0
#    while { [winfo exists $wwChartWindow] } {
#  incr nSuffix
#  set wwChartWindow .bcw$nSuffix
#    }


    # if the window doesn't exist already, make it.
    set bcw $wwChartWindow.bcw
    if { [winfo exists $wwChartWindow] == 0 } {

	# create window and set its size.
	toplevel $wwChartWindow
	wm geometry $wwChartWindow 600x800
	
	# create the chart. configure the x axis to call BarChart_GetXLabel
	# to get its labels. create two empty elements.
	blt::barchart $bcw -barmode aligned
	pack $bcw -expand true -fill both
	$bcw axis configure x \
	    -command { BarChart_GetXLabel } \
	    -rotate 90 \
	    -tickfont $kNormalFont
	$bcw element create V1
	$bcw element create V2
    }

    # set the window and chart title.
    wm title $wwChartWindow $tArgs(-title)
    $bcw config -title $tArgs(-title)
    
    # set the x axis labels.
    set glsXAxisLabels($bcw) $tArgs(-xAxisLabels)

    # set the label titles.
    $bcw axis config x -title $tArgs(-xAxisTitle)
    $bcw axis config y -title $tArgs(-yAxisTitle)
    
    # create a vector of indices for the elements. these are used
    # as indices into the x axis labels list.
    blt::vector vX
    vX seq 1 [llength $glsXAxisLabels($bcw)]

    # set the data in the two elements.
    $bcw element config V1 -label $tArgs(-label1) \
      -ydata $tArgs(-values1) -xdata vX -fg blue -bg blue
    $bcw element config V2 -label $tArgs(-label2) \
      -ydata $tArgs(-values2) -xdata vX -fg red -bg red
}


proc BarChart_GetXLabel { iwwTop ifValue } {

    global glsXAxisLabels

    if { [info exists glsXAxisLabels($iwwTop)] == 0 } {
  puts "error: labels for $iwwTop don't exist"
  return $ifValue
    }

    set nIndex [expr round($ifValue)]
    incr nIndex -1
    set sName [lindex $glsXAxisLabels($iwwTop) $nIndex]
    return $sName
}

# ================================================================== FUNCTIONS

proc AllowSaveThenQuit {} {
    global gbVolumeDirty

    if { $gbVolumeDirty } {
	    DoAskSaveChangesDlog;
    } else {
	    DoAskQuitDlog;
    }
}

proc SaveRGBSeries { isPrefix inBegin inEnd } {

    global gnVolX gnVolY gnVolZ gOrientation
    global mri_tOrientation_Sagittal mri_tOrientation_Horizontal 
    global mri_tOrientation_Coronal
    global mri_tCoordSpace_VolumeIdx

    dputs "SaveRGBSeries $isPrefix $inBegin $inEnd"

    # determine which way we're going
    if { $inBegin < $inEnd } {
	set nIncr 1 
    } else {
	set nIncr -1
    }
    
    set nX $gnVolX(cursor)
    set nY $gnVolY(cursor)
    set nZ $gnVolZ(cursor)
    for { set nSlice $inBegin } { $nSlice <= $inEnd } { incr nSlice $nIncr } {
	
	switch $gOrientation {
	    2 { set nX $nSlice }
	    1 { set nY $nSlice }
	    0 { set nZ $nSlice }
	}
	
	SetCursor $mri_tCoordSpace_VolumeIdx $nX $nY $nZ
	RedrawScreen
	SaveRGB $isPrefix[format "%03d" $nSlice].rgb
    }
}

proc SaveTIFFSeries { isPrefix inBegin inEnd } {

    global gnVolX gnVolY gnVolZ gOrientation
    global mri_tOrientation_Sagittal mri_tOrientation_Horizontal 
    global mri_tOrientation_Coronal
    global mri_tCoordSpace_VolumeIdx

    dputs "SaveTIFFSeries $isPrefix $inBegin $inEnd"

    # determine which way we're going
    if { $inBegin < $inEnd } {
	set nIncr 1 
    } else {
	set nIncr -1
    }
    
    set nX $gnVolX(cursor)
    set nY $gnVolY(cursor)
    set nZ $gnVolZ(cursor)
    for { set nSlice $inBegin } { $nSlice <= $inEnd } { incr nSlice $nIncr } {
	
	switch $gOrientation {
	    2 { set nX $nSlice }
	    1 { set nY $nSlice }
	    0 { set nZ $nSlice }
	}
	
	SetCursor $mri_tCoordSpace_VolumeIdx $nX $nY $nZ
	RedrawScreen
	SaveTIFF $isPrefix[format "%03d" $nSlice].tif
    }
}

proc ErrorDlog { isMsg } {

    global gwwTop

    tk_messageBox -type ok \
      -icon error \
      -message $isMsg \
      -title "Error" \
      -parent $gwwTop
}

proc FormattedErrorDlog { isTitle isMsg isDesc } {

    global gDialog
    global kLabelFont kNormalFont kSmallFont

    set wwDialog .wwFormattedErrorDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Error" {-borderwidth 10}] } {

  set fwText       $wwDialog.fwText
  set fwButtons    $wwDialog.fwButtons

  text $fwText -width 40 \
    -height 10 \
    -spacing3 10 \
    -relief flat \
    -wrap word
  $fwText insert end "Error: $isTitle \n" {tTitle}
  $fwText insert end "$isMsg \n" {tMsg}
  $fwText insert end "$isDesc \n" {tDesc}
  $fwText tag configure tTitle -font $kLabelFont
  $fwText tag configure tMsg -font $kNormalFont
  $fwText tag configure tDesc -font $kNormalFont
  $fwText configure -state disabled

  # button.
  tkm_MakeCloseButton $fwButtons $wwDialog

  pack $fwText $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }

}

proc AlertDlog { isMsg } {

    global gwwTop

    tk_messageBox -type ok \
      -icon info \
      -message $isMsg \
      -title "Note" \
      -parent $gwwTop
}

# ======================================================================= MAIN

CreateImages

# build the window
set wwTop        .w
set fwMenuBar    $wwTop.fwMenuBar
set fwToolBar    $wwTop.fwToolBar
set fwLeft       $wwTop.fwLeft
set fwRight      $wwTop.fwRight
set fwCursor     $fwLeft.fwCursor

CreateWindow         $wwTop
MakeKeyBindings      .

frame $fwLeft

CreateMenuBar        $fwMenuBar
CreateToolBar        $fwToolBar
CreateCursorFrame    $fwCursor
CreateMouseoverFrame $fwRight

# pack the window
pack $fwMenuBar $fwToolBar \
  -side top    \
  -expand true \
  -fill x      \
  -anchor w

pack $fwCursor \
  -side top         \
  -expand true      \
  -fill x           \
  -anchor nw

pack $fwLeft $fwRight \
  -side left    \
  -padx 3       \
  -pady 3       \
  -expand true  \
  -fill x       \
  -fill y       \
  -anchor nw

pack $wwTop

# start out with the main bar enabled
ShowToolBar main 1
ShowToolBar nav 1

# look for environment variable settings to automatically show toolbars
foreach toolbar {main nav recon} {
    catch {
	if { $env(TKMEDIT_TOOLBAR_[string toupper $toolbar]) == 1 } {
	    ShowToolBar $toolbar 1
	}
	if { $env(TKMEDIT_TOOLBAR_[string toupper $toolbar]) == 0 } {
	    ShowToolBar $toolbar 0
	}
    }
}

FsgdfPlot_Init

# lets us execute scripts from the command line but only after the
# window is open
after idle { catch { ExecuteQueuedScripts } }

dputs "Successfully parsed tkmedit.tcl"


# now try parsing the prefs files. first look in
# $FREESURFER_HOME/tktools/tkmedit_init.tcl, then
# $SUBJECTS_DIR/scripts/tkmedit_init.tcl, then
# $subject/scripts/tkmedit_init.tcl, then
# ~/scripts/tkmedit_init.tcl.
set lUserScripts {}
if { [info exists env(FREESURFER_HOME)] } {
    lappend lUserScripts $env(FREESURFER_HOME)/tktools/tkmedit_init.tcl
}
if { [info exists env(SUBJECTS_DIR)] } {
    lappend lUserScripts $env(SUBJECTS_DIR)/scripts/tkmedit_init.tcl
}
lappend lUserScripts $gsSubjectDirectory/scripts/tkmedit_init.tcl
lappend lUserScripts ~/tkmedit_init.tcl

foreach fnUserScript $lUserScripts {
    if { [file exists $fnUserScript] } {
	catch { 
	    dputs "Reading $fnUserScript"
	    source $fnUserScript
	    dputs "Successfully parsed $fnUserScript"
	}
    }
}

