#! /usr/bin/tixwish

# our graph code
source $env(MRI_DIR)/lib/tcl/tkm_graph.tcl
source $env(MRI_DIR)/lib/tcl/wrappers.tcl
source $env(MRI_DIR)/lib/tcl/tkm_dialog.tcl
source $env(MRI_DIR)/lib/tcl/tkm_wrappers.tcl

# hungarian notation notes:
# prefix  type                   library (if applicable)
# ww      window widget          tk
# fw      frame widget           tk
# slw     selector widget        tix
# cw      control widget         tix
# lw      label widget           tk
# lfw     label frame widget     tix
# vw      canvas widget          tk
# gw      graph widget           emu_graph
# ew      entry widget           tk
# lew     label entry widget     tix
# rw      radio button widget    tk
# cbw     checkbutton widget     tk
# bw      button widget          tk
# sw      scale (slider) widget  tk
# mbw     menu button widget     tk
# mw      menu widget            tk
# dw      dialog box widget      tix
# nw      notebook widget        tix
# cbxw    combobox widget        tix
# l       list                   tcl
# n       integer / number
# s       string
# f       float
# g       global

# function naming:
# functions prefixed with Set or xxxx_Set are c functions that the script
# calls to change variables. functions prefixed with Update are functions
# that the c code calls when variables are changed. when the user changes
# value from the tk interface, the change is 'suggested' with a Set command.
# the c code will accept or reject the change and then call the Update
# function with the new (or old) value.

# constants
set ksWindowName "TkMedit Tools"
set ksImageDir   "$env(MRI_DIR)/lib/images/"

# mri_tOrientation
set mri_tOrientation_Coronal    0
set mri_tOrientation_Horizontal 1
set mri_tOrientation_Sagittal   2

# DspA_tDisplayFlag
set DspA_tDisplayFlag_AuxVolume                   1
set DspA_tDisplayFlag_Cursor                      2
set DspA_tDisplayFlag_MainSurface                 3
set DspA_tDisplayFlag_OriginalSurface             4
set DspA_tDisplayFlag_CanonicalSurface            5
set DspA_tDisplayFlag_InterpolateSurfaceVertices  6
set DspA_tDisplayFlag_DisplaySurfaceVertices      7
set DspA_tDisplayFlag_ControlPoints               8
set DspA_tDisplayFlag_Selection                   9
set DspA_tDisplayFlag_FunctionalOverlay          10
set DspA_tDisplayFlag_MaskToFunctionalOverlay    11
set DspA_tDisplayFlag_ROIGroupOverlay            12
set DspA_tDisplayFlag_FocusFrame                 13
set DspA_tDisplayFlag_Axes                       14
set DspA_tDisplayFlag_MaxIntProj                 15
set DspA_tDisplayFlag_HeadPoints                 16

# DspA_tTool
set DspA_tTool_Navigate   0
set DspA_tTool_Select     1
set DspA_tTool_Edit       2
set DspA_tTool_CtrlPts    3

set ksaToolString(0) "navigate"
set ksaToolString(1) "select"
set ksaToolString(2) "edit"
set ksaToolString(3) "ctrlpts"

# DspA_tBrush
set DspA_tBrush_EditOne 0
set DspA_tBrush_EditTwo 1

set ksaBrushString(0) "Button 2"
set ksaBrushString(1) "Button 3"

# DspA_tBrushShape
set DspA_tBrushShape_Circle 0
set DspA_tBrushShape_Square 1

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
set tkm_tVolumeType_Main 0
set tkm_tVolumeType_Aux  1

set ksaDisplayedVolumeString(0) "main"
set ksaDisplayedVolumeString(1) "aux"

# Surf_tVertexSet
set Surf_tVertexSet_Main      0
set Surf_tVertexSet_Original  1
set Surf_tVertexSet_Canonical 2

# tFunctionalVolume
set tFunctionalVolume_Overlay    0
set tFunctionalVolume_TimeCourse 1

# mri_tCoordSpace
set mri_tCoordSpace_VolumeIdx 0
set mri_tCoordSpace_RAS       1
set mri_tCoordSpace_Talairach 2

set ksaLinkedCursorString(0) notlinked
set ksaLinkedCursorString(1) linked

# our global vars
set gOrientation 0
set gbLinkedCursor 1
set gbLinkedCursorString $ksaLinkedCursorString($gbLinkedCursor)
set gnVolX 0
set gnVolY 0
set gnVolZ 0
set gsVolCoords ""
set gfwVolCoords ""
set gfRASX 0
set gfRASY 0
set gfRASZ 0
set gsRASCoords ""
set gfwRASCoords ""
set gfTalX 0
set gfTalY 0
set gfTalZ 0
set gsTalCoords ""
set gfwTalCoords ""
set gsVolName ""
set gnVolValue 0
set gfwVolValue ""
set gsAuxVolName ""
set gnAuxVolValue 0
set gfwAuxVolValue ""
set gsROIGroupLabel ""
set gfwROIGroupLabel ""
set gsHeadPointLabel ""
set gfwHeadPointLabel ""
set gsFuncCoords ""
set gfwFuncCoords ""
set gfFuncValue 0
set gfwFuncValue ""
set gnZoomLevel 0
set gTool $DspA_tTool_Select
set gToolString $ksaToolString($gTool)
set gDisplayedVolume 0
set gDisplayedVolumeString $ksaDisplayedVolumeString($gDisplayedVolume)
set gDisplayCols 1
set gDisplayRows 1
set gViewPreset $tViewPreset_Single
set gViewPresetString $ksaViewPresetString($gViewPreset)
set gbCursor 0
set gbSelection 0
set gbControlPoints 0
set gbFunctional 0
set gbMaskFunctional 0
set gbROIGroup 0
set gbAxes 0
set gbMaxIntProj 0
set gbMainSurface 0
set gbOriginalSurface 0
set gbCanonicalSurface 0
set gbDisplaySurfaceVertices 0
set gbInterpolateSurfaceVertices 0
set gbHeadPoints 0
set gfwaToolBar(main)  ""
set gfwaToolBar(brush) ""

# brush info
set gBrush(radius)         1
set gBrush(shape)          $DspA_tBrushShape_Circle
set gBrush(3d)             true
foreach tool "$DspA_tBrush_EditOne $DspA_tBrush_EditTwo" {
    set gBrush(low,$tool)  0
    set gBrush(high,$tool) 0
    set gBrush(new,$tool)  0
}

# volume color scale
set gfVolumeColorScaleThresh($tkm_tVolumeType_Main) 0
set gfVolumeColorScaleSquash($tkm_tVolumeType_Main) 0
set gfVolumeColorScaleThresh($tkm_tVolumeType_Aux) 0
set gfVolumeColorScaleSquash($tkm_tVolumeType_Aux) 0

set gbVolumeDirty 0
set gbTalTransformPresent 0

# ========================================================= UPDATES FROM MEDIT

proc UpdateLinkedCursorFlag { ibLinked } {

    global gbLinkedCursor gbLinkedCursorString ksaLinkedCursorString
    set gbLinkedCursor $ibLinked
#    set gbLinkedCursorString $ksaLinkedCursorString($gbLinkedCursor)
}

proc UpdateVolumeCursor { inX inY inZ } {

    global gnVolX gnVolY gnVolZ gsVolCoords
    set gnVolX $inX
    set gnVolY $inY
    set gnVolZ $inZ

    # update the cursor string
    set gsVolCoords "($gnVolX, $gnVolY, $gnVolZ)"
}

proc UpdateRASCursor { ifX ifY ifZ } {

    global gfRASX gfRASY gfRASZ gsRASCoords
    set gfRASX $ifX
    set gfRASY $ifY
    set gfRASZ $ifZ

    # update the cursor string
    set gsRASCoords "($gfRASX, $gfRASY, $gfRASZ)"
}

proc UpdateTalCursor { ifX ifY ifZ } {

    global gfTalX gfTalY gfTalZ gsTalCoords
    set gfTalX $ifX
    set gfTalY $ifY
    set gfTalZ $ifZ

    # update the cursor string
    set gsTalCoords "($gfTalX, $gfTalY, $gfTalZ)"
}

proc UpdateVolumeName { isName } {

    global gsVolName
    set gsVolName $isName
}

proc UpdateVolumeValue { inValue } {

    global gnVolValue
    set gnVolValue $inValue
}

proc UpdateAuxVolumeName { isName } {

    global gsAuxVolName
    set gsAuxVolName $isName
}

proc UpdateAuxVolumeValue { inValue } {

    global gnAuxVolValue
    set gnAuxVolValue $inValue
}

proc UpdateROILabel { isLabel } {

    global gsROIGroupLabel
    set gsROIGroupLabel $isLabel
}

proc UpdateHeadPointLabel { isLabel } {

    global gsHeadPointLabel
    set gsHeadPointLabel $isLabel
}

proc UpdateFunctionalValue { ifValue } {

    global gfFuncValue
    set gfFuncValue $ifValue
}

proc UpdateFunctionalCoords { inX inY inZ } {

    global gsFuncCoords

    # update the cursor string
    set gsFuncCoords "($inX, $inY, $inZ)"
}

proc UpdateZoomLevel { inLevel } { 

    global gnZoomLevel
    set gnZoomLevel $inLevel
}

proc UpdateOrientation { iOrientation } {

    global gOrientation
    set gOrientation $iOrientation
}

proc UpdateDisplayFlag { iFlag ibValue } {

    # tcl doesn't make it easy to do constants.
    global DspA_tDisplayFlag_AuxVolume DspA_tDisplayFlag_Cursor
    global DspA_tDisplayFlag_MainSurface DspA_tDisplayFlag_OriginalSurface
    global DspA_tDisplayFlag_CanonicalSurface 
    global DspA_tDisplayFlag_InterpolateSurfaceVertices
    global DspA_tDisplayFlag_DisplaySurfaceVertices 
    global DspA_tDisplayFlag_ControlPoints
    global DspA_tDisplayFlag_Selection DspA_tDisplayFlag_FunctionalOverlay
    global DspA_tDisplayFlag_MaskToFunctionalOverlay
    global DspA_tDisplayFlag_ROIGroupOverlay
    global DspA_tDisplayFlag_Axes DspA_tDisplayFlag_MaxIntProj
    global DspA_tDisplayFlag_HeadPoints
    global tkm_tVolumeType_Main tkm_tVolumeType_Aux

    global gDisplayedVolume gDisplayedVolumeString ksaDisplayedVolumeString
    global gbCursor gbSelection gbControlPoints
    global gbFunctional gbMainSurface gbOriginalSurface gbCanonicalSurface
    global gbDisplaySurfaceVertices gbInterpolateSurfaceVertices
    global gbROIGroup gbAxes gbMaxIntProj gbHeadPoints gbMaskFunctional
    

    if { $DspA_tDisplayFlag_AuxVolume == $iFlag } {
  if { $ibValue == 0 } {
      set gDisplayedVolume $tkm_tVolumeType_Main
#      set gDisplayedVolumeString $ksaDisplayedVolumeString($tkm_tVolumeType_Main)
  } else {
      set gDisplayedVolume $tkm_tVolumeType_Aux
#      set gDisplayedVolumeString $ksaDisplayedVolumeString($tkm_tVolumeType_Aux)
  }
    }
    if { $DspA_tDisplayFlag_Cursor == $iFlag } {
  set gbCursor $ibValue
    }
    if { $DspA_tDisplayFlag_MainSurface == $iFlag } {
  set gbMainSurface $ibValue
    }
    if { $DspA_tDisplayFlag_OriginalSurface == $iFlag } {
  set gbOriginalSurface $ibValue
    }
    if { $DspA_tDisplayFlag_CanonicalSurface == $iFlag } {
  set gbCanonicalSurface $ibValue
    }
    if { $DspA_tDisplayFlag_InterpolateSurfaceVertices == $iFlag } {
  set gbInterpolateSurfaceVertices $ibValue
    }
    if { $DspA_tDisplayFlag_DisplaySurfaceVertices == $iFlag } {
  set gbDisplaySurfaceVertices $ibValue
    }
    if { $DspA_tDisplayFlag_ControlPoints == $iFlag } {
  set gbControlPoints $ibValue
    }
    if { $DspA_tDisplayFlag_Selection == $iFlag } {
  set gbSelection $ibValue
    }
    if { $DspA_tDisplayFlag_FunctionalOverlay == $iFlag } {
  set gbFunctional $ibValue
    }
    if { $DspA_tDisplayFlag_MaskToFunctionalOverlay == $iFlag } {
  set gbMaskFunctional $ibValue
    }
    if { $DspA_tDisplayFlag_ROIGroupOverlay == $iFlag } {
  set gbROIGroup $ibValue
    }
    if { $DspA_tDisplayFlag_Axes == $iFlag } {
  set gbAxes $ibValue
    }
    if { $DspA_tDisplayFlag_MaxIntProj == $iFlag } {
  set gbMaxIntProj $ibValue
    }
    if { $DspA_tDisplayFlag_HeadPoints == $iFlag } {
  set gbHeadPoints $ibValue
    }
}

proc UpdateTool { iTool } {

    global gTool ksaToolString gToolString
    set gTool $iTool
#    set gToolString $ksaToolString($gTool)
    
}

proc UpdateBrushShape { inRadius iShape ib3D } {

    global gBrush

    set gBrush(radius) $inRadius
    set gBrush(shape)  $iShape
    set gBrush(3d)     $ib3D
}

proc UpdateBrushInfo { inBrush inLow inHigh inNewValue } {

    global gBrush

    set gBrush(low,$inBrush)  $inLow
    set gBrush(high,$inBrush) $inHigh
    set gBrush(new,$inBrush)  $inNewValue
}

proc UpdateVolumeColorScaleInfo { inVolume inThresh inSquash } {

    global gfVolumeColorScaleThresh gfVolumeColorScaleSquash 

    set gfVolumeColorScaleThresh($inVolume) $inThresh
    set gfVolumeColorScaleSquash($inVolume) $inSquash
}

proc UpdateVolumeDirty { ibDirty } {

    global gbVolumeDirty
    set gbVolumeDirty $ibDirty
}

# =============================================================== DIALOG BOXES

proc DoLoadVolumeDlog {} {

    global gDialog

    set sVolumeName ""
    
    set wwDialog .wwLoadVolumeDialog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Volume" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set fwButtons $wwDialog.fwButtons

  # prompt and entry field
  tkm_MakeDirectorySelector $fwMain "Volume directory:" sVolumeName

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    "tkm_UpdateDirectorySelectorVariable $fwMain; \
    LoadVolume $sVolumeName" {}

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
   }
}

proc DoLoadAuxVolumeDlog {} {

    global gDialog

    set sAuxVolumeName ""
    
    set wwDialog .wwLoadAuxVolumeDialog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Aux Volume" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set fwButtons $wwDialog.fwButtons

  # prompt and entry field
  tkm_MakeDirectorySelector $fwMain "Aux volume directory:" \
    sAuxVolumeName

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    "tkm_UpdateDirectorySelectorVariable $fwMain;
    LoadAuxVolume \$sAuxVolumeName " {}

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
   }
}

proc DoSaveVolumeAsDlog {} {

    global gDialog

    set sVolumeDir ""
    set wwDialog .wwSaveVolumeAsDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save Volume As..." {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set fwButtons $wwDialog.fwButtons

  # prompt and entry field
  tkm_MakeDirectorySelector $fwMain "Directory to save volume in:" \
    sVolumeDir

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    "tkm_UpdateDirectorySelectorVariable $fwMain; \
    SaveVolumeAs \$sVolumeName" {}

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoSaveLabelDlog {} {

    global gDialog

    set sLabelName ""
    set wwDialog .wwSaveLabelDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save Label" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set fwButtons $wwDialog.fwButtons

  # prompt and entry field
  tkm_MakeFileSelector $fwMain "Save this label as:" sLabelName

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    "tkm_UpdateFileSelectorVariable $fwMain; \
    SaveLabel \$sLabelName" {}

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoLoadLabelDlog {} {

    global gDialog
    
    set wwDialog .wwLoadLabelDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Label" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set fwButtons $wwDialog.fwButtons

  # prompt and entry field
  tkm_MakeFileSelector $fwMain "Label file:" sFileName

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    "tkm_UpdateFileSelectorVariable $fwMain; \
    LoadLabel \$sFileName" {}

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc LoadSurface { isSurfaceName } {

    global Surf_tVertexSet_Main Surf_tVertexSet_Original
    global Surf_tVertexSet_Canonical
    global gLoadingSurface

    if { $Surf_tVertexSet_Main == $gLoadingSurface } { 
  LoadMainSurface $isSurfaceName
    }

    if { $Surf_tVertexSet_Original == $gLoadingSurface } {
  LoadOriginalSurface $isSurfaceName
    }
    
    if { $Surf_tVertexSet_Canonical == $gLoadingSurface } {
  LoadCanonicalSurface $isSurfaceName
    }
}

proc DoLoadSurfaceDlog { iSurface } {

    global gDialog
    global gLoadingSurface
 
    set gLoadingSurface $iSurface
    set sSurfaceName ""
    set wwDialog .wwLoadSurfaceDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Surface" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set fwButtons $wwDialog.fwButtons

  # prompt and entry field
  tkm_MakeFileSelector $fwMain "Surface file:" sSurfaceName

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    "tkm_UpdateFileSelectorVariable $fwMain; \
    LoadSurface \$sSurfaceName" {}

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc FindVertex { inVertex } {

    global Surf_tVertexSet_Main Surf_tVertexSet_Original
    global Surf_tVertexSet_Canonical
    global gFindingSurface

    if { $Surf_tVertexSet_Main   == $gFindingSurface } {
  GotoMainVertex $inVertex
    }
    if { $Surf_tVertexSet_Original  == $gFindingSurface } {
  GotoOriginalVertex $inVertex
    }
    if { $Surf_tVertexSet_Canonical == $gFindingSurface } {
  GotoCanonicalVertex $inVertex
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

proc LoadFunctionalVolume { isPath isStem } {

    global tFunctionalVolume_Overlay tFunctionalVolume_TimeCourse
    global gLoadingVolume

    if { $gLoadingVolume == $tFunctionalVolume_Overlay } {
  Overlay_LoadData $isPath $isStem
    }
    if { $gLoadingVolume == $tFunctionalVolume_TimeCourse } {
  TimeCourse_LoadData $isPath $isStem
    }
}

proc DoLoadFunctionalDlog { iVolume } {

    global tFunctionalVolume_Overlay tFunctionalVolume_TimeCourse
    global gDialog
    global gLoadingVolume

    set gLoadingVolume $iVolume
    set sPath ""
    set sStem ""
    set wwDialog .wwLoadFunctionalDlog

    set sWindowName "Load Functional Volume"

    # try to create the dlog...
    if { [Dialog_Create $wwDialog $sWindowName {-borderwidth 10}] } {

  set fwPath    $wwDialog.fwPath
  set fwStem    $wwDialog.fwStem
  set fwButtons $wwDialog.fwButtons

  # stem and path
  tkm_MakeDirectorySelector $fwPath \
    "Enter the directory the data is in:" sPath
  tkm_MakeEntry $fwStem "Enter the stem of the data:" sStem 10

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    "tkm_UpdateDirectorySelectorVariable $fwPath; \
    LoadFunctionalVolume \$sPath \$sStem"

  pack $fwPath $fwStem $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
   }
}

proc DoPrintTimeCourseDlog {} {

    global gDialog

    set sSummaryName ""
    set wwDialog .wwPrintTimeCourseDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Print Time Course" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set fwButtons $wwDialog.fwButtons

  # prompt and entry field
  tkm_MakeFileSelector $fwMain "Time course summary file:" sSummaryName

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    "tkm_UpdateFileSelectorVariable $fwMain; \
    TimeCourse_PrintSelectionRangeToFile \$sSummaryName"

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

proc DoSaveDlog {} {

    global gDialog

    set wwDialog .wwSaveDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save Volume" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set fwButtons $wwDialog.fwButtons

  # prompt
  tkm_MakeNormalLabel $fwMain\
    "Are you sure you wish to save changes to the volume?"

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { SaveVolume }

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
    global ksaBrushString
    global DspA_tBrushShape_Square DspA_tBrushShape_Circle
    global DspA_tBrush_EditOne DspA_tBrush_EditTwo
    global gBrush

    set wwDialog .wwBrushInfoDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Brush Size" {-borderwidth 10}] } {
  
  set fwTop                $wwDialog.fwTop
  set fwShape              $fwTop.fwShape
  set fwInfo               $fwTop.fwInfo
  set fwButtons            $wwDialog.fwButtons
  
  frame $fwTop
  
  tixLabelFrame $fwShape \
    -label "Shape" \
    -labelside acrosstop \
    -options { label.padX 5 }
  
  set fwShapeSubFrame [$fwShape subwidget frame]
  set fwRadiusScale        $fwShapeSubFrame.fwRadiusScale
  set fwShapeLabel         $fwShapeSubFrame.fwShapeLabel
  set fwCircle             $fwShapeSubFrame.fwCircle
  set fwSquare             $fwShapeSubFrame.fwSquare
  set fw3DCheckbox         $fwShapeSubFrame.fw3DCheckbox
  
  # radius
  tkm_MakeSlider $fwRadiusScale "Radius" gBrush(radius) 1 20 50 "" 1
    
  # shape radio buttons
  tkm_MakeNormalLabel $fwShapeLabel "Shape"
  tkm_MakeRadioButton $fwCircle "Circle" \
    gBrush(shape) $DspA_tBrushShape_Circle ""
  tkm_MakeRadioButton $fwSquare "Square" \
    gBrush(shape) $DspA_tBrushShape_Square ""

  # 3d checkbox
  tkm_MakeCheckbox $fw3DCheckbox "3D" gBrush(3d) ""
  
  # pack them in a column
  pack $fwRadiusScale $fwShapeLabel $fwCircle \
    $fwSquare $fw3DCheckbox             \
    -side top                           \
    -anchor w                           \
    -expand yes                         \
    -fill x
  
  tixLabelFrame $fwInfo \
    -label "Edit Tool" \
    -labelside acrosstop \
    -options { label.padX 5 }
  
  set fwInfoSubFrame [$fwInfo subwidget frame]
  set nbwInfo $fwInfoSubFrame.nbwInfo

  tixNoteBook $nbwInfo
  foreach tool "$DspA_tBrush_EditOne $DspA_tBrush_EditTwo" {

      $nbwInfo add pane$tool -label $ksaBrushString($tool)

      set fw [$nbwInfo subwidget pane$tool]
      set fwLowScale        $fw.fwLowScale
      set fwHighScale       $fw.fwHighScale
      set fwNewValueScale   $fw.fwNewValueScale
      set fwDefaults        $fw.fwDefaults

      # low, high, and new value sliders
      tkm_MakeSlider $fwLowScale "Low" gBrush(low,$tool) \
        0 255 100 "" 1
      tkm_MakeSlider $fwHighScale "High" gBrush(high,$tool) \
        0 255 100 "" 1
      tkm_MakeSlider $fwNewValueScale {"New Value"} gBrush(new,$tool) \
        0 255 100 "" 1

      # defaults button
      tkm_MakeButtons $fwDefaults \
        { {text "Restore Defaults" "SetBrushInfoToDefaults $tool"} }

      # pack them in a column
      pack $fwLowScale $fwHighScale $fwNewValueScale $fwDefaults \
        -side top                           \
        -anchor w                           \
        -expand yes                         \
        -fill x
  }

  pack $nbwInfo

  # pack the sides
  pack $fwShape $fwInfo \
    -side left \
    -padx 2 \
    -anchor nw 

  # buttons. 
  tkm_MakeApplyCloseButtons $fwButtons $wwDialog \
    { SetBrushConfiguration } 

  pack $fwTop $fwButtons \
    -side top       \
    -expand yes     \
    -fill x
   }
}

proc DoVolumeColorScaleInfoDlog { } {

    global gDialog
    global gfVolumeColorScaleThresh gfVolumeColorScaleSquash 
    global gfSavedVolumeColorScaleThresh gfSavedVolumeColorScaleSquash 

    set wwDialog .wwVolumeColorScaleInfoDlog

    set gfSavedVolumeColorScaleSquash(0) $gfVolumeColorScaleSquash(0)
    set gfSavedVolumeColorScaleThresh(0) $gfVolumeColorScaleThresh(0)
    set gfSavedVolumeColorScaleSquash(1) $gfVolumeColorScaleSquash(1)
    set gfSavedVolumeColorScaleThresh(1) $gfVolumeColorScaleThresh(1)


    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Brightness / Contrast" {-borderwidth 10}] } {

  set fwBrightness $wwDialog.fwBrightness
  set fwContrast   $wwDialog.fwContrast
  set fwAuxBrightness $wwDialog.fwAuxBrightness
  set fwAuxContrast   $wwDialog.fwAuxContrast
  set fwButtons    $wwDialog.fwButtons

  
  # brightness and contrast sliders
  tkm_MakeSlider $fwBrightness "Brightness" \
    gfVolumeColorScaleThresh(0) \
    1 0 100 "" 1 0.01
  tkm_MakeSlider $fwContrast "Contrast" \
    gfVolumeColorScaleSquash(0) \
    0 20 100 "" 1

  # aux brightness and contrast sliders
  tkm_MakeSlider $fwAuxBrightness "\"Aux Brightness\"" \
    gfVolumeColorScaleThresh(1) \
    1 0 100 "" 1 0.01
  tkm_MakeSlider $fwAuxContrast "\"Aux Contrast\"" \
    gfVolumeColorScaleSquash(1) \
    0 20 100 "" 1

  # pack them in a column
  pack $fwBrightness $fwContrast          \
    $fwAuxBrightness $fwAuxContrast \
    -side top                \
    -anchor w                \
    -expand yes              \
    -fill x
  
  # buttons
  tkm_MakeApplyCloseButtons $fwButtons $wwDialog \
    { SetVolumeColorScale 0 \
    $gfVolumeColorScaleThresh(0) \
    $gfVolumeColorScaleSquash(0);\
    SetVolumeColorScale 1 \
    $gfVolumeColorScaleThresh(1) \
    $gfVolumeColorScaleSquash(1) }

  pack $fwButtons  \
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
  set fwThresh      $wwDialog.fwThresh
  set fwNewValue    $wwDialog.fwNewValue
  set fwButtons     $wwDialog.fwButtons

  # label 
  tkm_MakeNormalLabel $fwLabel "Change all values"
  
  # direction radios
  tkm_MakeRadioButton $fwAbove "above" bAbove 1
  tkm_MakeRadioButton $fwBelow "below" bAbove 0

  # threshold value
  tkm_MakeSlider $fwThresh "\"this value\"" nThreshold 0 255 200 "" 1

  #$fwBelow new value
  tkm_MakeSlider $fwNewValue "\"to this value\"" nNewValue 0 255 200 "" 1

  # pack them in a column
  pack $fwLabel $fwAbove $fwBelow \
    $fwThresh $fwNewValue     \
    -side top                \
    -anchor w                \
    -expand yes              \
    -fill x

  # buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { Threshold $nThreshold $bAbove $nNewValue }

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
  set fwAxesRads    $wwDialog.fwAxesRads
  set fwX           $fwAxesRads.fwRotateX
  set fwY           $fwAxesRads.fwRotateY
  set fwZ           $fwAxesRads.fwRotateZ
  set fwButtons     $wwDialog.fwRotateButtons

  set fRotateDegrees 0
  set sRotateDirection x

  # degrees
  tkm_MakeEntryWithIncDecButtons \
    $fwDegrees "Degrees" \
    fRotateDegrees \
    {} \
    0.5

  # direction radios
  frame $fwAxesRads
  tkm_MakeNormalLabel $fwDirection "Around anatomical axis:"
  tkm_MakeRadioButton $fwX "X" sRotateDirection x
  tkm_MakeRadioButton $fwY "Y" sRotateDirection y
  tkm_MakeRadioButton $fwZ "Z" sRotateDirection z
  pack $fwX $fwY $fwZ  \
    -side left       \
    -anchor w       \
    -expand yes     \
    -fill x 

  # buttons. 
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { RotateVolume $sRotateDirection $fRotateDegrees }

  pack $fwDegrees $fwDirection $fwAxesRads $fwButtons \
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

  set fwFlipX   $wwDialog.fwFlipX
  set fwFlipY   $wwDialog.fwFlipY
  set fwFlipZ   $wwDialog.fwFlipZ
  set fwButtons $wwDialog.fwButtons

  # flip checks
  tkm_MakeCheckbox $fwFlipX "Flip around X axis" bFlipX {}
  tkm_MakeCheckbox $fwFlipY "Flip around Y axis" bFlipY {}
  tkm_MakeCheckbox $fwFlipZ "Flip around Z axis" bFlipZ {}

  # buttons. 
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { FlipVolume $bFlipX $bFlipY $bFlipZ }

  pack $fwFlipX $fwFlipY $fwFlipZ $fwButtons \
    -side top       \
    -anchor w       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoLoadHeadPtsDlog {} {

    global gDialog

    set sHeadPtsName ""
    set sTransformName ""
    set wwDialog .wwLoadHeadPtsDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Head Points" {-borderwidth 10}] } {

  set fwHeadPts    $wwDialog.fwHeadPts
  set fwTransform  $wwDialog.fwTransform
  set fwButtons    $wwDialog.fwButtons

  # prompt and entry fields
  tkm_MakeFileSelector $fwHeadPts "Head points file:" sHeadPtsName
  tkm_MakeFileSelector $fwTransform "Transform file:" sTransformName
  
  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    "tkm_UpdateFileSelectorVariable $fwHeadPts; \
    tkm_UpdateFileSelectorVariable $fwTransform; \
    LoadHeadPts \$sHeadPtsName \$sTransformName"

  pack $fwHeadPts $fwTransform $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoRotateHeadPtsDlog {} {

    global gDialog
    
    set wwDialog .wwRotateHeadPtsDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Rotate" {-borderwidth 10}] } {

  set fwDegrees     $wwDialog.fwRotateDegrees
  set fwDirection   $wwDialog.fwRotateDirection
  set fwAxesRads    $wwDialog.fwAxesRads
  set fwX           $fwAxesRads.fwRotateX
  set fwY           $fwAxesRads.fwRotateY
  set fwZ           $fwAxesRads.fwRotateZ
  set fwButtons     $wwDialog.fwRotateButtons

  set fRotateDegrees 0
  set sRotateDirection x

  # degrees
  tkm_MakeEntryWithIncDecButtons \
    $fwDegrees "Degrees" \
    fRotateDegrees \
    {} \
    0.5

  # direction radios
  frame $fwAxesRads
  tkm_MakeNormalLabel $fwDirection "Around axis:"
  tkm_MakeRadioButton $fwX "X" sRotateDirection x
  tkm_MakeRadioButton $fwY "Y" sRotateDirection y
  tkm_MakeRadioButton $fwZ "Z" sRotateDirection z
  pack $fwX $fwY $fwZ  \
    -side left       \
    -anchor w       \
    -expand yes     \
    -fill x

  # buttons. 
  tkm_MakeApplyCloseButtons $fwButtons $wwDialog \
    { RotateHeadPts $fRotateDegrees $sRotateDirection }

  pack $fwDegrees $fwDirection $fwAxesRads \
    $fwButtons \
    -side top       \
    -anchor w       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoTranslateHeadPtsDlog {} {

    global gDialog
    
    set wwDialog .wwTranslateHeadPtsDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Translate" {-borderwidth 10}] } {

  set fwDistance    $wwDialog.fwTransDistance
  set fwDirection   $wwDialog.fwTransDirection
  set fwX           $wwDialog.fwTransX
  set fwY           $wwDialog.fwTransY
  set fwButtons     $wwDialog.fwTransButtons

  set fTransDirection 0
  set sTransDirection x

  # degrees
  tkm_MakeEntryWithIncDecButtons \
    $fwDistance "Distance" \
    fTransDistance \
    {} \
    1

  # direction radios
  tkm_MakeNormalLabel $fwDirection "Direction:"
  tkm_MakeRadioButton $fwX "X" sTransDirection x
  tkm_MakeRadioButton $fwY "Y" sTransDirection y

  # buttons. 
  tkm_MakeApplyCloseButtons $fwButtons $wwDialog \
    { TranslateHeadPts $fTransDistance $sTransDirection }

  pack $fwDistance $fwDirection $fwX $fwY $fwButtons \
    -side top       \
    -anchor w       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoLoadParcellationDlog {} {

    global gDialog

    set sVolume ""
    set sColorFile ""
    set wwDialog .wwLoadROIGroupDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Parcellation" {-borderwidth 10}] } {

  set fwPath      $wwDialog.fwPath
  set fwColorFile $wwDialog.fwColorFile
  set fwButtons   $wwDialog.fwButtons

  # volume prompt and entry field
  tkm_MakeDirectorySelector $fwPath \
    "Volume directory" sVolume

  # color prompt and entry field
  tkm_MakeFileSelector $fwColorFile \
    "Color look up table file:" sColorFile

  # ok and cancel buttons. 
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    "tkm_UpdateDirectorySelectorVariable $fwPath; \
    tkm_UpdateFileSelectorVariable $fwColorFile; \
    LoadParcellationVolume \$sVolume \$sColorFile"

  pack $fwPath $fwColorFile $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoSaveRGBDlog {} {

    global gDialog

    set sRGBName ""
    set wwDialog .wwSaveRGBDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save RGB" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set fwButtons $wwDialog.fwButtons

  # prompt and entry field
  tkm_MakeFileSelector $fwMain \
    "RGB file:" sRGBName

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    "tkm_UpdateFileSelectorVariable $fwMain; \
    SaveRGB \$sRGBName"

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoGotoPointDlog {} {

    global gDialog
    global mri_tCoordSpace_VolumeIdx mri_tCoordSpace_RAS 
    global mri_tCoordSpace_Talairach
    global gnVolX gnVolY gnVolZ
    global gbTalTransformPresent
    
    set wwDialog .wwGotoPointDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Goto Point" {-borderwidth 10}] } {

  set fwLabel       $wwDialog.fwLabel
  set fwCoordSpace  $wwDialog.fwCoordSpace
  set fwVolumeIdx   $fwCoordSpace.fwVolumeIdx
  set fwRAS         $fwCoordSpace.fwRAS
  set fwTalCoords   $fwCoordSpace.fwTalCoords
  set fwWhere       $wwDialog.fwWhere
  set fwX           $fwWhere.fwX
  set fwY           $fwWhere.fwY
  set fwZ           $fwWhere.fwZ
  set fwButtons     $wwDialog.fwButtons

  set fX $gnVolX
  set fY $gnVolY
  set fZ $gnVolZ
  set coordSpace $mri_tCoordSpace_VolumeIdx

  # coord space radios
  tkm_MakeNormalLabel $fwLabel "Coordinate space:"
  frame $fwCoordSpace
  tkm_MakeRadioButton $fwVolumeIdx "Volume Index" \
    coordSpace $mri_tCoordSpace_VolumeIdx
  tkm_MakeRadioButton $fwRAS "RAS" \
    coordSpace $mri_tCoordSpace_RAS
  pack $fwLabel $fwVolumeIdx $fwRAS \
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

proc DoEditHeadPointLabelDlog {} {

    global gDialog
    global gsHeadPointLabel
    
    set wwDialog .wwEditHeadPointLabelDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Edit Head Point Label" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain

  set fwButtons $wwDialog.fwButtons
  set bwOK      $fwButtons.bwOK
  set bwCancel  $fwButtons.bwCancel

  # prompt and entry field
  tkm_MakeEntry $fwMain "Change label to:" gsHeadPointLabel 20
  
  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { SetSelectedHeadPointLabel $gsHeadPointLabel } {}

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc SetBrushConfiguration { } {

    global gBrush
    global DspA_tBrush_EditOne DspA_tBrush_EditTwo

    SetBrushShape $gBrush(radius) $gBrush(shape) $gBrush(3d)
    foreach tool "$DspA_tBrush_EditOne $DspA_tBrush_EditTwo" {
  SetBrushInfo $tool $gBrush(low,$tool) $gBrush(high,$tool) $gBrush(new,$tool)
    }
    
}

# ======================================================== INTERFACE MODIFIERS

proc ShowVolumeCoords { ibShow } {

    global gfwVolCoords gfwRASCoords gfwTalCoords gfwVolValue

    if { $ibShow == 1 } {
  
  # try to pack before ras coords
  if {[catch { pack $gfwVolCoords  \
    -before $gfwRASCoords    \
    -side top                \
    -anchor w } sResult] } {
      
      # failed, try to pack before tal coords
      if {[catch { pack $gfwVolCoords \
        -before $gfwTalCoords   \
        -side top               \
        -anchor w } sResult] } {
    
    # failed, try vol value
    if {[catch { pack $gfwVolCoords  \
      -before $gfwVolValue     \
      -side top                \
      -anchor w } sResult] } {
        
        # damn it
        return;
    }
      }
  }
      
    } else {
  pack forget $gfwVolCoords
    }
}

proc ShowRASCoords { ibShow } {

    global gfwVolCoords gfwRASCoords gfwTalCoords gfwVolValue

    if { $ibShow == 1 } {

  # try packing before tal coords
  if {[catch { pack $gfwRASCoords            \
    -before $gfwTalCoords \
    -side top             \
    -anchor w } sResult] } {

      # failed, try vol value
      if {[catch { pack $gfwRASCoords  \
        -before $gfwVolValue     \
        -side top                \
        -anchor w } sResult] } {
    
    return;
      }
  }

    } else {
  pack forget $gfwRASCoords
    }
}

proc ShowTalCoords { ibShow } {

    global gfwVolCoords gfwRASCoords gfwTalCoords gfwVolValue
    global gbTalTransformPresent

    # this is used in the gotoPoint dialog
    set gbTalTransformPresent $ibShow

    if { $ibShow == 1 } {
  if {[catch { pack $gfwTalCoords            \
    -before $gfwVolValue  \
    -side top             \
    -anchor w } sResult] } {

      return;
  }
    } else {
  pack forget $gfwTalCoords
    }
}

proc ShowAuxValue { ibShow } {

    global gfwVolValue gfwAuxVolValue

    if { $ibShow == 1 } {
  if {[catch { pack $gfwAuxVolValue  \
    -after $gfwVolValue      \
    -side top                \
    -anchor w } sResult] } {

      puts "ShowAuxValue: pack failed: $sResult"
      return;
  }
    } else {
  pack forget $gfwAuxVolValue
    }
}

proc ShowROILabel { ibShow } {

    global gfwVolValue gfwAuxVolValue gfwROIGroupLabel

    if { $ibShow == 1 } {
  if {[catch { pack $gfwROIGroupLabel  \
    -after $gfwAuxVolValue      \
    -side top                \
    -anchor w } sResult] } {
      if {[catch { pack $gfwROIGroupLabel  \
        -after $gfwVolValue      \
        -side top                \
        -anchor w } sResult] } {
    puts "ShowROILabel pack failed: $sResult"
    return;
      }
  }
    } else {
  pack forget $gfwROIGroupLabel
    }
}

proc ShowHeadPointLabel { ibShow } {

    global gfwVolValue gfwAuxVolValue gfwROIGroupLabel gfwHeadPointLabel

    if { $ibShow == 1 } {
  if {[catch { pack $gfwHeadPointLabel  \
    -after $gfwROIGroupLabel      \
    -side top                \
    -anchor w } sResult] } {
      if {[catch { pack $gfwHeadPointLabel  \
        -after $gfwAuxVolValue      \
        -side top                \
        -anchor w } sResult] } {
    if {[catch { pack $gfwHeadPointLabel  \
      -after $gfwVolValue      \
      -side top                \
      -anchor w } sResult] } {
        puts "ShowHeadPointLabel pack failed: $sResult"
        return;
    }
      }
  }
    } else {
  pack forget $gfwHeadPointLabel
    }

    # set menu items status
    tkm_SetMenuItemGroupStatus tMenuGroup_HeadPoints $ibShow
}

proc ShowFuncCoords { ibShow } {

    global gfwAuxVolValue gfwVolValue gfwFuncCoords gfFuncValue

    if { $ibShow == 1 } {
  if {[catch { pack $gfwFuncCoords            \
    -before $gfwFuncValue  \
    -side top             \
    -anchor w } sResult] } {
      
      if {[catch { pack $gfwFuncCoords  \
        -after $gfwAuxVolValue      \
        -side top                \
        -anchor w } sResult] } {
    
    if {[catch { pack $gfwFuncCoords  \
      -after $gfwVolValue      \
      -side top                \
      -anchor w } sResult] } {

        puts "ShowFuncCoords: pack failed: $sResult"
        return;
    }
      }
  }

    } else {
  pack forget $gfwFuncCoords
    }
}

proc ShowFuncValue { ibShow } {

    global gfwVolValue gfwAuxVolValue gfwFuncValue gfwFuncCoords
    global DspA_tDisplayFlag_FunctionalOverlay

    if { $ibShow == 1 } {
  if {[catch { pack $gfwFuncValue  \
    -after $gfwFuncCoords      \
    -side top                \
    -anchor w } sResult] } {
      
      if {[catch { pack $gfwFuncValue  \
        -after $gfwAuxVolValue      \
        -side top                \
        -anchor w } sResult] } {
    
    if {[catch { pack $gfwFuncValue  \
      -after $gfwVolValue      \
      -side top                \
      -anchor w } sResult] } {
        
        puts "ShowFuncValue: pack failed: $sResult"
        return;
    }
      }
  }

  # show item
  SetDisplayFlag $DspA_tDisplayFlag_FunctionalOverlay 1

    } else {
  pack forget $gfwFuncValue
    }
}

# =============================================================== VIEW PRESETS

proc SetViewPreset { iPreset } {

    global gViewPreset gViewPresetString ksaViewPresetString
    global tViewPreset_Single tViewPreset_Multiple tViewPreset_Mosaic
    global MWin_tLinkPolicy_None MWin_tLinkPolicy_MultipleOrientations
    global MWin_tLinkPolicy_Mosaic

    set gViewPreset $iPreset
#    set gViewPresetString $ksaViewPresetString($gViewPreset)

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
}

proc CreateMenuBar { ifwMenuBar } {

    global mri_tOrientation_Sagittal mri_tOrientation_Horizontal 
    global mri_tOrientation_Coronal
    global DspA_tTool_Navigate DspA_tTool_Select
    global DspA_tTool_Edit DspA_tTool_CtrlPts 
    global gnVolX gnVolY gnVolZ
    global gDisplayCols gDisplayRows gViewPreset
    global tViewPreset_Single tViewPreset_Multiple tViewPreset_Mosaic    
    global gbCursor gbSelection gbControlPoints gbFunctional gbMaskFunctional
    global gbMainSurface gbOriginalSurface gbCanonicalSurface
    global gbDisplaySurfaceVertices gbInterpolateSurfaceVertices
    global gbROIGroup gbAxes gbMaxIntProj gbHeadPoints
    global gTool

    set mbwFile   $ifwMenuBar.mbwFile
    set mbwEdit   $ifwMenuBar.mbwEdit
    set mbwView   $ifwMenuBar.mbwView
    set mbwTools  $ifwMenuBar.mbwTools

    frame $ifwMenuBar -border 2 -relief raised

    # file menu button
    tkm_MakeMenu $mbwFile "File" { \
      { command \
      "Load Main Volume..." \
      DoLoadVolumeDlog } \
      \
      { command \
      "Load Aux Volume..." \
      DoLoadAuxVolumeDlog } \
      \
      { command \
      "Save Volume" \
      DoSaveDlog } \
      \
      { command \
      "Save Volume As..." \
      DoSaveVolumeAsDlog } \
      \
      { separator } \
      \
      { command \
      "Load Main Surface..." \
      "DoLoadSurfaceDlog $Surf_tVertexSet_Main" } \
      \
      { cascade "Load Vertex Set..." { \
      { command \
      "Original Verticies" \
      "DoLoadSurfaceDlog $Surf_tVertexSet_Original" \
      tMenuGroup_SurfaceLoading } \
      \
      { command \
      "Pial Verticies " \
      "DoLoadSurfaceDlog $Surf_tVertexSet_Canonical" \
      tMenuGroup_SurfaceLoading } } } \
      \
      { command \
      "Unload Surface" \
      UnloadSurface \
      tMenuGroup_SurfaceLoading } \
      \
      { separator } \
      \
      { command \
      "Load Overlay Data..." \
      "DoLoadFunctionalDlog $tFunctionalVolume_Overlay" } \
      \
      { command \
      "Load Time Course Data..." \
      "DoLoadFunctionalDlog $tFunctionalVolume_TimeCourse" } \
      \
      { command \
      "Save Overlay Registration" \
      Overlay_SaveRegistration \
      tMenuGroup_Registration } \
      \
      { separator } \
      \
      { command \
      "Load Parcellation..." \
      DoLoadParcellationDlog } \
      \
      { command \
      "Load ROI Group..." \
      DoLoadParcellationDlog } \
      \
      { separator } \
      \
      { command \
      "Load Label..." \
      DoLoadLabelDlog } \
      \
      { command \
      "Save Label As..." \
      DoSaveLabelDlog } \
      \
      { separator } \
      \
      { command \
      "Load Head Points..." \
      DoLoadHeadPtsDlog } \
      \
      { command \
      "Save Head Point Transform" \
      WriteHeadPointsTransform \
      tMenuGroup_HeadPoints } \
      \
      { command \
      "Save Head Points" \
      WriteHeadPointsFile \
      tMenuGroup_HeadPoints } \
      \
      { separator } \
      \
      { command \
      "Save Control Points" \
      WriteControlPointFile } \
      \
      { separator } \
      \
      { command \
      "Quit" \
      QuitMedit } }
   
    # edit menu 
    tkm_MakeMenu $mbwEdit "Edit" { \
      { command \
      "Undo Last Edit" \
      UndoLastEdit } \
      \
      { separator } \
      \
      { command \
      "Take Snapshot of Volume" \
      SnapshotVolume } \
      \
      { command \
      "Restore Volume from Snapshot" \
      RestoreVolumeFromSnapshot } \
      \
      { separator } \
      \
      { command \
      "Clear Label Selection" \
      ClearSelection } \
      \
      { command \
      "Clear Control Points Selection" \
      DeselectAllControlPoints }   }

    # view menu
    tkm_MakeMenu $mbwView "View" { \
      { cascade \
      "View Configurations" { \
      { radio  \
      "Single View" \
      "SetViewPreset $tViewPreset_None" \
      gViewPreset \
      0 } \
      \
      { radio  \
      "Multiple Orientations" \
      "SetViewPreset $tViewPreset_Multiple" \
      gViewPreset \
      1 } \
      \
      { radio  \
      "Mosaic" \
      "SetViewPreset $tViewPreset_Mosaic" \
      gViewPreset \
      2 } } } \
      \
      { separator } \
      \
      { check \
      "Brush Tool Bar" \
      "ShowToolBar brush $gbShowToolBar" \
      gbShowToolBar } \
      \
      { separator } \
      \
      { cascade "Configure..." { \
      { command \
      "Brightness / Contrast..." \
      DoVolumeColorScaleInfoDlog } \
      \
      { command \
      "Functional Overlay..." \
      Overlay_DoConfigDlog \
      tMenuGroup_OverlayOptions } \
      \
      { command \
      "Time Course Graph..." \
      TimeCourse_DoConfigDlog \
      tMenuGroup_TimeCourseOptions } } }\
      \
      { separator } \
      \
      { radio  \
      "Main Volume" \
      "SetDisplayFlag $DspA_tDisplayFlag_AuxVolume \
      $tkm_tVolumeType_Main"  \
      gDisplayedVolume  \
      0 } \
      \
      { radio \
      "Aux Volume" \
      "SetDisplayFlag $DspA_tDisplayFlag_AuxVolume \
      $tkm_tVolumeType_Aux"  \
      gDisplayedVolume  \
      1 } \
      \
      { separator } \
      \
      { check \
      "Maximum Intensity Projection" \
      "SetDisplayFlag $DspA_tDisplayFlag_MaxIntProj \
      \$gbMaxIntProj" \
      gbMaxIntProj } \
      \
      { check \
      "Main Surface" \
      "SetDisplayFlag $DspA_tDisplayFlag_MainSurface \$gbMainSurface" \
      gbMainSurface \
      tMenuGroup_SurfaceViewing } \
      \
      { check \
      "Original Surface" \
      "SetDisplayFlag $DspA_tDisplayFlag_OriginalSurface \
      \$gbOriginalSurface"  \
      gbOriginalSurface \
      tMenuGroup_OriginalSurfaceViewing } \
      \
      { check \
      "Pial Surface" \
      "SetDisplayFlag $DspA_tDisplayFlag_CanonicalSurface \
      \$gbCanonicalSurface"  \
      gbCanonicalSurface \
      tMenuGroup_CanonicalSurfaceViewing } \
      \
      { check \
      "Surface Vertices" \
      "SetDisplayFlag $DspA_tDisplayFlag_DisplaySurfaceVertices \
      \$gbDisplaySurfaceVertices" \
      gbDisplaySurfaceVertices \
      tMenuGroup_SurfaceViewing } \
      \
      { check \
      "Interpolate Surface Vertices" \
      "SetDisplayFlag $DspA_tDisplayFlag_InterpolateSurfaceVertices \
      \$gbInterpolateSurfaceVertices"  \
      gbInterpolateSurfaceVertices \
      tMenuGroup_SurfaceViewing } \
      \
      { check \
      "Functional Overlay" \
      "SetDisplayFlag $DspA_tDisplayFlag_FunctionalOverlay \
      \$gbFunctional" \
      gbFunctional \
      tMenuGroup_OverlayOptions } \
      \
      { check \
      "Mask to Functional Overlay" \
      "SetDisplayFlag $DspA_tDisplayFlag_MaskToFunctionalOverlay \
      \$gbMaskFunctional" \
      gbMaskFunctional \
      tMenuGroup_OverlayOptions } \
      \
      { check \
      "ROI Group Overlay" \
      "SetDisplayFlag $DspA_tDisplayFlag_ROIGroupOverlay \
      \$gbROIGroup" \
      gbROIGroup } \
      \
      { check \
      "Selection / Label" \
      "SetDisplayFlag $DspA_tDisplayFlag_Selection \$gbSelection" \
      gbSelection } \
      \
      { check \
      "Head Points" \
      "SetDisplayFlag $DspA_tDisplayFlag_HeadPoints \
      \$gbHeadPoints" \
      gbHeadPoints \
      tMenuGroup_HeadPoints } \
      \
      { check \
      "Control Points" \
      "SetDisplayFlag $DspA_tDisplayFlag_ControlPoints \
      \$gbControlPoints" \
      gbControlPoints } \
      \
      { check \
      "Cursor" \
      "SetDisplayFlag $DspA_tDisplayFlag_Cursor \$gbCursor" \
      gbCursor } \
      \
      { check \
      "Axes" \
      "SetDisplayFlag $DspA_tDisplayFlag_Axes \
      \$gbAxes" \
      gbAxes } }

#      { radio \
#      "Navigate" \
#      "SetTool $DspA_tTool_Navigate" \
#      gTool \
#      0 } \
#      \
    # tools menu
    tkm_MakeMenu $mbwTools "Tools" { \
      { radio \
      "Select Voxels" \
      "SetTool $DspA_tTool_Select" \
      gTool \
      1 } \
      \
      { radio \
      "Edit Voxels" \
      "SetTool $DspA_tTool_Edit" \
      gTool \
      2 } \
      \
      { radio \
      "Select Ctrl Pts" \
      "SetTool $DspA_tTool_CtrlPts" \
      gTool \
      3 } \
      \
      { separator } \
      \
      { command \
      "Configure Brush..." \
      DoBrushInfoDlog } \
      \
      { separator } \
      \
      { command \
      "Save Point" \
      SendCursor } \
      \
      { command \
      "Goto Saved Point" \
      ReadCursor } \
      \
      { separator } \
      \
      { command \
      "Goto Point..." \
      DoGotoPointDlog } \
      \
      { separator } \
      \
      { cascade "Volume" { \
      { command \
      "Threshold Volume..." \
      DoThresholdDlog } \
      \
      { command \
      "Flip Volume..." \
      DoFlipVolumeDlog } \
      \
      { command \
      "Rotate Volume..." \
      DoRotateVolumeDlog } } } \
      \
      { cascade "Surface" { \
      { command  \
      "Show Nearest Main Vertex" \
      ShowNearestMainVertex \
      tMenuGroup_SurfaceViewing }\
      \
      { command \
      "Show Nearest Original Vertex" \
      ShowNearestOriginalVertex \
      tMenuGroup_OriginalSurfaceViewing } \
      \
      { command \
      "Show Nearest Pial Vertex" \
      ShowNearestCanonicalVertex \
      tMenuGroup_CanonicalSurfaceViewing } \
      \
      { command \
      "Find Main Vertex..." \
      { DoFindVertexDlog $Surf_tVertexSet_Main }  \
      tMenuGroup_SurfaceViewing } \
      \
      { command \
      "Find Original Vertex..." \
      { DoFindVertexDlog $Surf_tVertexSet_Original } \
      tMenuGroup_OriginalSurfaceViewing } \
      \
      { command \
      "Find Pial Vertex..." \
      { DoFindVertexDlog $Surf_tVertexSet_Canonical } \
      tMenuGroup_CanonicalSurfaceViewing } } } \
      \
      { cascade "fMRI" { \
      { command \
      "Register Functional Overlay" \
      { DoRegisterOverlayDlog }\
      tMenuGroup_Registration } \
      \
      { command \
      "Restore Overlay Registration" \
      { Overlay_RestoreRegistration }\
      tMenuGroup_Registration } \
      \
      { command \
      "Set Registration to Identity" \
      { Overlay_SetRegistrationToIdentity }\
      tMenuGroup_Registration } \
      \
      { command \
      "Graph Current ROI Time Course" \
      { GraphSelectedRegion }\
      tMenuGroup_TimeCourseOptions } \
      \
      { command \
      "Print Time Course Summary to File.." \
      { DoPrintTimeCourseDlog } \
      tMenuGroup_TimeCourseOptions } } } \
      \
      { cascade "ROI Group" { \
      { command "Edit ROI Label" {} } \
      { command "Select Current ROI" {} } } } \
      \
      { cascade "Head Points" { \
      { command \
      "Restore Head Points" \
      RestoreHeadPts \
      tMenuGroup_HeadPoints } \
      \
      { command \
      "Edit Current Head Point Label..." \
      DoEditHeadPointLabelDlog \
      tMenuGroup_HeadPoints } \
      \
      { command \
      "Translate Head Points..." \
      DoTranslateHeadPtsDlog \
      tMenuGroup_HeadPoints } \
      \
      { command \
      "Rotate Head Points..." \
      DoRotateHeadPtsDlog \
      tMenuGroup_HeadPoints } } } \
      \
      { cascade "Control Points" { \
      { command \
      "New Control Point" \
      { NewControlPoint } } \
      \
      { command \
      "Delete Selected Control Point" \
      DeleteSelectedControlPoints } } } \
      \
      { command \
      "Save RGB..." \
      DoSaveRGBDlog } }


    pack $mbwFile $mbwEdit $mbwView $mbwTools \
      -side left
}

proc CreateCursorFrame { ifwCursor } {

    global gfwVolCoords gfwRASCoords gfwTalCoords gfwVolValue
    global gfwAuxVolValue gfwFuncValue
    global gbLinkedCursor gsVolCoords gsRASCoords 
    global gsTalCoords gsVolName gnVolValue  gsFuncCoords
    global gsAuxVolName gnAuxVolValue gfwFuncCoords gfFuncValue
    global gsROIGroupLabel gfwROIGroupLabel
    global gsHeadPointLabel gfwHeadPointLabel

    frame $ifwCursor

    set fwLabel             $ifwCursor.fwLabel
    set lwLabel             $fwLabel.lwLabel
    
    set fwLinkCheckbox      $ifwCursor.fwLinkCheckbox
    set cwLink              $fwLinkCheckbox.cwLink
    
    set gfwVolCoords        $ifwCursor.fwVolCoords
    set gfwRASCoords        $ifwCursor.fwRASCoords
    set gfwTalCoords        $ifwCursor.fwTalCoords

    set gfwVolValue         $ifwCursor.fwVolValue
    set fwVolValueLabel     $gfwVolValue.fwVolValueLabel
    set fwVolValue          $gfwVolValue.fwVolValue

    set gfwAuxVolValue      $ifwCursor.fwAuxVolValue
    set fwAuxVolValueLabel  $gfwAuxVolValue.fwAuxVolValueLabel
    set fwAuxVolValue       $gfwAuxVolValue.fwAuxVolValue

    set gfwROIGroupLabel $ifwCursor.fwROIGroupLabel

    set gfwHeadPointLabel   $ifwCursor.fwHeadPointLabel

    set gfwFuncCoords       $ifwCursor.fwFuncCoords

    set gfwFuncValue        $ifwCursor.fwFuncValue

    # the label that goes at the top of the frame
    tkm_MakeBigLabel $fwLabel "Cursor"

    # link checkbox
    tkm_MakeCheckbox $fwLinkCheckbox "Linked Cursors" \
      gbLinkedCursor { SetLinkedCursorFlag $gbLinkedCursor }

    # the volume coords
    tkm_MakeActiveLabel $gfwVolCoords "Volume" gsVolCoords

    # the RAS coords
    tkm_MakeActiveLabel $gfwRASCoords "RAS" gsRASCoords

    # the Talairach coords
    tkm_MakeActiveLabel $gfwTalCoords "Talairach" gsTalCoords

    # the volume value
    frame $gfwVolValue
    tkm_MakeActiveLabel $fwVolValueLabel "" gsVolName
    tkm_MakeActiveLabel $fwVolValue "" gnVolValue
    pack $fwVolValueLabel $fwVolValue -side left \
      -anchor w

    # the aux volume value
    frame $gfwAuxVolValue
    tkm_MakeActiveLabel $fwAuxVolValueLabel "" gsAuxVolName
    tkm_MakeActiveLabel $fwAuxVolValue "" gnAuxVolValue
    pack $fwAuxVolValueLabel $fwAuxVolValue -side left \
      -anchor w

    # the parcellation label
    tkm_MakeActiveLabel $gfwROIGroupLabel "ROI: " gsROIGroupLabel

    # the head point label
    tkm_MakeActiveLabel $gfwHeadPointLabel "Head Point: " gsHeadPointLabel

    # the Func coords
    tkm_MakeActiveLabel $gfwFuncCoords "Overlay" gsFuncCoords

    # the functional value
    tkm_MakeActiveLabel $gfwFuncValue "Overlay value" gfFuncValue

    # pack the subframes in a column. don't pack vol coords, aux value,
    # or func value. these can be explicity shown later.
    pack $fwLabel $fwLinkCheckbox $gfwRASCoords \
      $gfwTalCoords $gfwVolValue          \
      -side top                           \
      -anchor w

}

proc CreateDisplayFrame { ifwDisplay } {

    global mri_tOrientation_Sagittal mri_tOrientation_Horizontal 
    global mri_tOrientation_Coronal
    global gOrientation gnVolX gnVolY gnVolZ gnZoomLevel
    global gb3D gbCursor gbSelection gbControlPoints gbFunctional
    global gbMainSurface gbOriginalSurface gbCanonicalSurface
    global gbDisplaySurfaceVertices gbInterpolateSurfaceVertices

    frame $ifwDisplay

    set fwLabel            $ifwDisplay.fwLabel
    set fwSagittalButton   $ifwDisplay.fwSagittalButton
    set fwSagittalScale    $ifwDisplay.fwSagittalScale
    set fwHorizontalButton $ifwDisplay.fwHorizontalButton
    set fwHorizontalScale  $ifwDisplay.fwHorizontalScale
    set fwCoronalButton    $ifwDisplay.fwCoronalButton
    set fwCoronalScale     $ifwDisplay.fwCoronalScale
    set fwZoomLabel        $ifwDisplay.fwZoomLabel
    set fwZoomScale        $ifwDisplay.fwZoomScale

    # the label that goes at the top of the frame
    tkm_MakeBigLabel $fwLabel "Display"

    # sagittal orientation radio button
    tkm_MakeRadioButton $fwSagittalButton \
      "Sagittal (R -- L)" gOrientation $mri_tOrientation_Sagittal \
      { SetOrientation $gOrientation }

    # sagittal slice slider
    tkm_MakeSlider $fwSagittalScale "" gnVolX 0 255 200 \
      { SetCursor $mri_tCoordSpace_VolumeIdx $gnVolX $gnVolY $gnVolZ } \
      1

    # horizontal orientaiton radio button
    tkm_MakeRadioButton $fwHorizontalButton \
      "Horizontal (Sup -- Inf)" gOrientation \
      $mri_tOrientation_Horizontal \
      { SetOrientation $gOrientation }

    # horiontal slice slider
    tkm_MakeSlider $fwHorizontalScale "" gnVolY 0 255 200 \
      { SetCursor $mri_tCoordSpace_VolumeIdx $gnVolX $gnVolY $gnVolZ } \
      1

    # coronal orientaiton radio button
    tkm_MakeRadioButton $fwCoronalButton \
      "Coronal (Post -- Ant)" gOrientation \
      $mri_tOrientation_Coronal \
      { SetOrientation $gOrientation }

    # coronal slice slider
    tkm_MakeSlider $fwCoronalScale "" gnVolZ 0 255 200 \
      { SetCursor $mri_tCoordSpace_VolumeIdx $gnVolX $gnVolY $gnVolZ } \
      1

    # zoom label
    tkm_MakeNormalLabel $fwZoomLabel "Zoom (Out -- In)"

    # zoom slice slider
    tkm_MakeSlider $fwZoomScale "" gnZoomLevel 1 16 200 \
      { SetZoomLevel $gnZoomLevel } \
      1

    # pack them all together in a column
    pack $fwLabel $fwSagittalButton $fwSagittalScale   \
      $fwHorizontalButton $fwHorizontalScale     \
      $fwCoronalButton $fwCoronalScale           \
      $fwZoomLabel $fwZoomScale                  \
      -side top                                  \
      -anchor w                                  \
      -expand yes                                \
      -fill x
}

proc CreateToolBar { ifwToolBar } {

    global gfwaToolBar
    global gDisplayedVolumeString ksaDisplayedVolume
    global gToolString ksaToolString
    global gViewPresetString ksaViewPresetString
    global gbLinkedCursor
    global gBrush DspA_tBrushShape_Square DspA_tBrushShape_Circle

    frame $ifwToolBar

    set gfwaToolBar(main)  $ifwToolBar.fwMainBar
    set fwTools            $gfwaToolBar(main).fwTools
    set fwViews            $gfwaToolBar(main).fwViews
    set fwPoint            $gfwaToolBar(main).fwPoint
    set fwVolumeToggles    $gfwaToolBar(main).fwVolumeToggles
    set fwLinkedCursor     $gfwaToolBar(main).fwLinkedCursor

    set gfwaToolBar(brush) $ifwToolBar.fwBrushBar
    set fwCircle           $gfwaToolBar(brush).fwCircle
    set fwSquare           $gfwaToolBar(brush).fwSquare
    set fw3D               $gfwaToolBar(brush).fw3D
    set fwRadius           $gfwaToolBar(brush).fwRadius

    frame $gfwaToolBar(main) -border 2 -relief raised
    
    tkm_MakeToolbar $fwTools \
      1 \
      gToolString \
      UpdateToolWrapper { \
      { image select icon_edit_label } \
      { image edit icon_edit_volume } \
      { image ctrlpts icon_control_point } }

    tkm_MakeToolbar $fwViews \
      1 \
      gViewPresetString \
      UpdateViewPresetWrapper { \
      { image single icon_view_single } \
      { image multiple icon_view_multiple } \
      { image mosaic icon_view_mosaic } }
    
    tkm_MakeButtons $fwPoint { \
      { image icon_cursor_save {SendCursor} } \
      { image icon_cursor_goto {ReadCursor} } }

    tkm_MakeToolbar $fwVolumeToggles \
      1 \
      gDisplayedVolumeString \
      UpdateVolumeToggleWrapper { \
      { image main icon_main_volume } \
      { image aux icon_aux_volume } }

#    image create photo icon_linked_cursors -file icon_linked_cursors.gif
    
    #    tkm_MakeToolbar $fwLinkedCursor \
      #      0 \
      #      gbLinkedCursorString \
      #      UpdateLinkedCursorWrapper { \
      #      { image linked icon_linked_cursors } }
    
    frame $gfwaToolBar(brush) -border 2 -relief raised

    tkm_MakeRadioButton $fwCircle "Circle" \
    gBrush(shape) $DspA_tBrushShape_Circle "SetBrushConfiguration"
    tkm_MakeRadioButton $fwSquare "Square" \
      gBrush(shape) $DspA_tBrushShape_Square "SetBrushConfiguration"

    tkm_MakeCheckbox $fw3D "3D" gBrush(3d) "SetBrushConfiguration"

    tkm_MakeSlider $fwRadius "Radius" gBrush(radius) 1 20 80 "SetBrushConfiguration" 1

    pack $fwTools $fwViews $fwPoint $fwVolumeToggles \
      -side left \
      -anchor w \
      -padx 5

    pack $fwCircle $fwSquare $fw3D $fwRadius \
      -side left \
      -anchor w \
      -padx 5

    pack $gfwaToolBar(main) $gfwaToolBar(brush) \
      -side top \
      -fill x \
      -expand yes
}

proc ShowToolBar { isWhich ibShow } {

    global gfwaToolBar

    if { $ibShow == 1 } {   

  pack $gfwaToolBar($isWhich) \
      -side top \
      -fill x \
      -expand yes \
      -after $gfwaToolBar(main)

    } else {

  pack forget $gfwaToolBar($isWhich)
    }
}


proc CreateImages {} {

    global ksImageDir

    foreach image_name { icon_edit_label icon_edit_volume \
      icon_control_point \
      icon_view_single icon_view_multiple icon_view_mosaic \
      icon_cursor_goto icon_cursor_save \
      icon_main_volume icon_aux_volume \
      icon_arrow_up icon_arrow_down icon_arrow_left icon_arrow_right \
      icon_arrow_cw icon_arrow_ccw \
      icon_arrow_expand_x icon_arrow_expand_y \
      icon_arrow_shrink_x icon_arrow_shrink_y } {

  image create photo  $image_name -file \
    [ file join $ksImageDir $image_name.gif ]
    }
}

proc UpdateVolumeToggleWrapper { iValue ibStatus } {

    global DspA_tDisplayFlag_AuxVolume
    global ksaDisplayedVolumeString tkm_tVolumeType_Main tkm_tVolumeType_Aux

    if { $ibStatus == 1 } {
  foreach volume "$tkm_tVolumeType_Aux $tkm_tVolumeType_Main" {
      if { [string compare $iValue $ksaDisplayedVolumeString($volume)] == 0 } {
    SetDisplayFlag $DspA_tDisplayFlag_AuxVolume $volume
      }
  }
    }
}

proc UpdateToolWrapper { iValue ibStatus } {

    global ksaToolString
    global DspA_tTool_Navigate DspA_tTool_Select DspA_tTool_Edit
    global DspA_tTool_CtrlPts 

    if { $ibStatus == 1 } {
  foreach tool "$DspA_tTool_Navigate $DspA_tTool_Select \
    $DspA_tTool_Edit $DspA_tTool_CtrlPts" {
      if { [string compare $iValue $ksaToolString($tool)] == 0 } {
    SetTool $tool
      }
  }
    }
}

proc UpdateViewPresetWrapper { iValue ibStatus } {

    global ksaViewPresetString
    global tViewPreset_Single tViewPreset_Multiple tViewPreset_Mosaic

    if { $ibStatus == 1 } {
  foreach preset "$tViewPreset_Single $tViewPreset_Multiple $tViewPreset_Mosaic" {
      if { [string compare $iValue $ksaViewPresetString($preset)] == 0 } {
    SetViewPreset $preset
      }
  }
    }
}

proc UpdateLinkedCursorWrapper { iValue ibStatus } {

    SetLinkedCursorFlag $iValue
}

# = ====================================================================== MAIN

fixcolors

CreateImages

# build the window
set wwTop        .w
set fwMenuBar    $wwTop.fwMenuBar
set fwToolBar    $wwTop.fwToolBar
set fwLeft       $wwTop.fwLeft
set fwRight      $wwTop.fwRight
set fwCursor     $fwLeft.fwCursor

CreateWindow         $wwTop

frame $fwLeft

CreateMenuBar        $fwMenuBar
CreateToolBar        $fwToolBar
CreateCursorFrame    $fwCursor
CreateDisplayFrame   $fwRight

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

# start out with the brush toolbar disabled
ShowToolBar brush 0

proc ErrorDlog { isMsg } {

    global gwwTop

    tk_messageBox -type ok \
      -icon error \
      -message $isMsg \
      -title "Error" \
      -parent $gwwTop
}

proc AlertDlog { isMsg } {

    global gwwTop

    tk_messageBox -type ok \
      -icon info \
      -message $isMsg \
      -title "Note" \
      -parent $gwwTop
}

