#! /usr/bin/wish

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

# DspA_tOrientation
set DspA_tOrientation_Coronal    0
set DspA_tOrientation_Horizontal 1
set DspA_tOrientation_Sagittal   2

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
set DspA_tDisplayFlag_ParcellationOverlay        11
set DspA_tDisplayFlag_FocusFrame                 12
set DspA_tDisplayFlag_Axes                       13
set DspA_tDisplayFlag_MaxIntProj                 14

# DspA_tTool
set DspA_tTool_Select     0
set DspA_tTool_Edit       1
set DspA_tTool_CtrlPts    2
set DspA_tTool_CustomEdit 3

# DspA_tBrushShape
set DspA_tBrushShape_Circle 0
set DspA_tBrushShape_Square 1

# MWin_tDisplayConfiguration
set MWin_tDisplayConfiguration_1x1 0
set MWin_tDisplayConfiguration_2x2 1

# MWin_tVolumeType
set MWin_tVolumeType_Main 0
set MWin_tVolumeType_Aux  1

# tkm_tSurfaceType
set tkm_tSurfaceType_Current   0
set tkm_tSurfaceType_Original  1
set tkm_tSurfaceType_Canonical 2

# tFunctionalVolume
set tFunctionalVolume_Overlay    0
set tFunctionalVolume_TimeCourse 1

# our global vars
set gOrientation 0
set gbLinkedCursor 1
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
set gsParcellationLabel ""
set gfwParcellationLabel ""
set gsFuncCoords ""
set gfwFuncCoords ""
set gfFuncValue 0
set gfwFuncValue ""
set gnZoomLevel 0
set gTool $DspA_tTool_Select
set gDisplayedVolume 0
set gb3DDisplay 0
set gbCursor 0
set gbSelection 0
set gbControlPoints 0
set gbFunctional 0
set gbParcellation 0
set gbAxes 0
set gbMaxIntProj 0
set gbMainSurface 0
set gbOriginalSurface 0
set gbCanonicalSurface 0
set gbDisplaySurfaceVertices 0
set gbInterpolateSurfaceVertices 0
set gnBrushRadius 1
set gBrushShape $DspA_tBrushShape_Circle
set gbBrush3D true
set gnBrushLow 0
set gnBrushHigh 0
set gnBrushNewValue 0
set gfVolumeColorScaleThresh 0
set gfVolumeColorScaleSquash 0
set gbVolumeDirty 0

foreach lMenuGroup {SurfaceLoading} {

    set glMenuGroups($lMenuGroup) {}
}

# ========================================================= UPDATES FROM MEDIT

proc UpdateLinkedCursorFlag { ibLinked } {

    global gbLinkedCursor
    set gbLinkedCursor $ibLinked
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

proc UpdateParcellationLabel { isLabel } {

    global gsParcellationLabel
    set gsParcellationLabel $isLabel
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
    global DspA_tDisplayFlag_ParcellationOverlay
    global DspA_tDisplayFlag_Axes DspA_tDisplayFlag_MaxIntProj
    global MWin_tVolumeType_Main MWin_tVolumeType_Aux

    global gDisplayedVolume gb3DDisplay gbCursor gbSelection gbControlPoints
    global gbFunctional gbMainSurface gbOriginalSurface gbCanonicalSurface
    global gbDisplaySurfaceVertices gbInterpolateSurfaceVertices
    global gbParcellation gbAxes gbMaxIntProj
    

    if { $DspA_tDisplayFlag_AuxVolume == $iFlag } {
  if { $ibValue == 0 } {
      set gDisplayedVolume $MWin_tVolumeType_Main
  } else {
      set gDisplayedVolume $MWin_tVolumeType_Aux
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
    if { $DspA_tDisplayFlag_ParcellationOverlay == $iFlag } {
  set gbParcellation $ibValue
    }
    if { $DspA_tDisplayFlag_Axes == $iFlag } {
  set gbAxes $ibValue
    }
    if { $DspA_tDisplayFlag_MaxIntProj == $iFlag } {
  set gbMaxIntProj $ibValue
    }
}

proc UpdateTool { iTool } {

    global gTool
    set gTool $iTool
}

proc UpdateBrush { inRadius iShape ib3D } {

    global gnBrushRadius gBrushShape gbBrush3D

    set gnBrushRadius $inRadius
    set gBrushShape   $iShape
    set gbBrush3D     $ib3D
}

proc UpdateBrushThreshold { inLow inHigh inNewValue } {

    global gnBrushLow gnBrushHigh gnBrushNewValue

    set gnBrushLow      $inLow
    set gnBrushHigh     $inHigh
    set gnBrushNewValue $inNewValue
}

proc UpdateVolumeColorScaleInfo { inThresh inSquash } {

    global gfVolumeColorScaleThresh gfVolumeColorScaleSquash 

    set gfVolumeColorScaleThresh $inThresh
    set gfVolumeColorScaleSquash $inSquash
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
  tkm_MakeEntry $fwMain "Enter a volume name:" sVolumeName 10

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { LoadVolume $sVolumeName } {}

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
  tkm_MakeEntry $fwMain "Enter a volume name:" sAuxVolumeName 10

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { LoadAuxVolume $sAuxVolumeName } {}

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
  tkm_MakeDirectorySelector $fwMain "Path to save this volume in:" \
    sVolumeDir

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { SaveVolumeAs $sVolumeName } {}

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
    { SaveLabel $sLabelName } {}

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

    set sLabelName ""
    set wwDialog .wwLoadLabelDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Label" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set fwButtons $wwDialog.fwButtons

  # prompt and entry field
  tkm_MakeFileSelector $fwMain "Enter a label name:" sLabelName
  
  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { LoadLabel $sLabelName } {}

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc LoadSurface { isSurfaceName } {

    global tkm_tSurfaceType_Current tkm_tSurfaceType_Original
    global tkm_tSurfaceType_Canonical
    global gLoadingSurface

    if { $tkm_tSurfaceType_Current == $gLoadingSurface } { 
  LoadMainSurface $isSurfaceName
    }

    if { $tkm_tSurfaceType_Original == $gLoadingSurface } {
  LoadOriginalSurface $isSurfaceName
    }
    
    if { $tkm_tSurfaceType_Canonical == $gLoadingSurface } {
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
  tkm_MakeFileSelector $fwMain "Enter a surface name:" sSurfaceName

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { LoadSurface $sSurfaceName } {}

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc FindVertex { inVertex } {

    global tkm_tSurfaceType_Current tkm_tSurfaceType_Original
    global tkm_tSurfaceType_Canonical
    global gFindingSurface

    if { $tkm_tSurfaceType_Current   == $gFindingSurface } {
  GotoMainVertex $inVertex
    }
    if { $tkm_tSurfaceType_Original  == $gFindingSurface } {
  GotoOriginalVertex $inVertex
    }
    if { $tkm_tSurfaceType_Canonical == $gFindingSurface } {
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
  tkm_MakeEntry $fwMain "Enter a vertex number:" nVertex 6
  
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
    { LoadFunctionalVolume $sPath $sStem }

  pack $fwPath $fwStem $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
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
    global DspA_tBrushShape_Square DspA_tBrushShape_Circle
    global gnBrushRadius gBrushShape gbBrush3D
    global gnSavedBrushRadius gSavedBrushShape gbSavedBrush3D

    set wwDialog .wwBrushInfoDlog

    set gnSavedBrushRadius $gnBrushRadius
    set gSavedBrushShape   $gBrushShape
    set gbSavedBrush3D     $gbBrush3D

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Brush Size" {-borderwidth 10}] } {
  
  set fwRadiusScale        $wwDialog.fwRadiusScale
  set fwShapeLabel         $wwDialog.fwShapeLabel
  set fwCircle             $wwDialog.fwCircle
  set fwSquare             $wwDialog.fwSquare
  set fw3DCheckbox         $wwDialog.fw3DCheckbox
  set fwButtons            $wwDialog.fwButtons
  
  # radius
  tkm_MakeSlider $fwRadiusScale "Radius" gnBrushRadius 1 100 100 "" 1
    
  # shape radio buttons
  tkm_MakeNormalLabel $fwShapeLabel "Shape"
  tkm_MakeRadioButton $fwCircle "Circle" \
    gBrushShape $DspA_tBrushShape_Circle ""
  tkm_MakeRadioButton $fwSquare "Square" \
    gBrushShape $DspA_tBrushShape_Square ""

  # 3d checkbox
  tkm_MakeCheckbox $fw3DCheckbox "3D" gbBrush3D ""
  
  # pack them in a column
  pack $fwRadiusScale $fwShapeLabel $fwCircle \
    $fwSquare $fw3DCheckbox             \
    -side top                           \
    -anchor w                           \
    -expand yes                         \
    -fill x
  
  # buttons. 
  tkm_MakeCancelApplyOKButtons $fwButtons $wwDialog \
    { SetBrush $gnBrushRadius $gBrushShape $gbBrush3D } \
    { SetBrush $gnSavedBrushRadius $gSavedBrushShape $gbSavedBrush3D }

  pack $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
   }
}

proc DoCustomBrushDlog {} {

    global gDialog
    global gnBrushLow gnBrushHigh gnBrushNewValue
    global gnSavedBrushLow gnSavedBrushHigh gnSavedBrushNewValue

    set wwDialog .wwCustomBrushDlog

    set gnSavedBrushLow      $gnBrushLow
    set gnSavedBrushHigh     $gnBrushHigh
    set gnSavedBrushNewValue $gnBrushNewValue

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Custom Brush" {-borderwidth 10}] } {

  set fwLowScale        $wwDialog.fwLowScale
  set fwHighScale        $wwDialog.fwHighScale
  set fwNewValueScale        $wwDialog.fwNewValueScale
  set fwButtons            $wwDialog.fwButtons
  
  # low, high, and new value sliders
  tkm_MakeSlider $fwLowScale "Low" gnBrushLow 0 255 200 "" 1
  tkm_MakeSlider $fwHighScale "High" gnBrushHigh 0 255 200 "" 1
  tkm_MakeSlider $fwNewValueScale "New Value" gnBrushNewValue \
    0 255 200 "" 1

  # pack them in a column
  pack $fwLowScale $fwHighScale $fwNewValueScale \
    -side top                           \
    -anchor w                           \
    -expand yes                         \
    -fill x
  
  # buttons. 
  tkm_MakeCancelApplyOKButtons $fwButtons $wwDialog \
    { SetBrushThreshold $gnBrushLow $gnBrushHigh $gnBrushNewValue } \
    {SetBrushThreshold $gnSavedBrushLow $gnSavedBrushHigh $gnSavedBrushNewValue  }

  pack $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
   }
}

proc DoVolumeColorScaleInfoDlog {} {

    global gDialog
    global gfVolumeColorScaleThresh gfVolumeColorScaleSquash 
    global gfSavedVolumeColorScaleThresh gfSavedVolumeColorScaleSquash 
 
    set wwDialog .wwVolumeColorScaleInfoDlog

    set gfSavedVolumeColorScaleSquash $gfVolumeColorScaleSquash
    set gfSavedVolumeColorScaleThresh $gfVolumeColorScaleThresh

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Brightness / Contrast" {-borderwidth 10}] } {

  set fwBrightness $wwDialog.fwBrightness
  set fwContrast   $wwDialog.fwContrast
  set fwButtons    $wwDialog.fwButtons
  
  # brightness and contrast sliders
  tkm_MakeSlider $fwBrightness "Brightness" gfVolumeColorScaleThresh \
    1 0 100 "" 1 0.01
  tkm_MakeSlider $fwContrast "Constrast" gfVolumeColorScaleSquash \
    0 20 100 "" 1

  # pack them in a column
  pack $fwBrightness $fwContrast \
    -side top                \
    -anchor w                \
    -expand yes              \
    -fill x
  
  # buttons
  tkm_MakeCancelApplyOKButtons $fwButtons $wwDialog \
    { SetVolumeColorScale $gfVolumeColorScaleThresh \
    $gfVolumeColorScaleSquash } \
    { SetVolumeColorScale $gfSavedVolumeColorScaleThresh \
    $gfSavedVolumeColorScaleSquash }

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
  tkm_MakeSlider $fwThresh "this value" nThreshold 0 255 200 "" 1

  #$fwBelow new value
  tkm_MakeSlider $fwNewValue "to this value" nNewValue 0 255 200 "" 1

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

proc DoRotateDlog {} {

    global gDialog
    
    set wwDialog .wwRotateDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Rotate" {-borderwidth 10}] } {

  set fwDegrees     $wwDialog.fwDegrees
  set fwDirection   $wwDialog.fwDirection
  set fwX           $wwDialog.fwX
  set fwY           $wwDialog.fwY
  set fwZ           $wwDialog.fwZ
  set fwButtons     $wwDialog.fwButtons

  set fDegrees 0
  set sDirection x

  # degrees
  tkm_MakeEntryWithIncDecButtons \
    $fwDegrees "Degrees" \
    fDegrees \
    "set fDegrees \$fDegrees" \
    "set fDegrees \[expr \$fDegrees - 0.5\]" \
    "set fDegrees \[expr \$fDegrees + 0.5\]" \
    5

  # direction radios
  tkm_MakeNormalLabel $fwDirection "Direction:"
  tkm_MakeRadioButton $fwX "X" sDirection x
  tkm_MakeRadioButton $fwY "Y" sDirection y
  tkm_MakeRadioButton $fwZ "Z" sDirection z

  # buttons. 
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { RotateBrain $fDegrees $sDirection }

  pack $fwDegrees $fwDirection $fwX $fwY $fwZ $fwButtons \
    -side top       \
    -anchor w       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoTranslateDlog {} {

    global gDialog
    
    set wwDialog .wwTranslateDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Translate" {-borderwidth 10}] } {

  set fwDistance    $wwDialog.fwDistance
  set fwDirection   $wwDialog.fwDirection
  set fwX           $wwDialog.fwX
  set fwY           $wwDialog.fwY
  set fwZ           $wwDialog.fwZ
  set fwButtons     $wwDialog.fwButtons

  set fDirection 0
  set sDirection x

  # degrees
  tkm_MakeEntryWithIncDecButtons \
    $fwDistance "Distance" \
    fDistance \
    "set fDistance \$fDistance" \
    "set fDistance \[expr \$fDistance - 0.5\]" \
    "set fDistance \[expr \$fDistance + 0.5\]" \
    5

  # direction radios
  tkm_MakeNormalLabel $fwDirection "Direction:"
  tkm_MakeRadioButton $fwX "X" sDirection x
  tkm_MakeRadioButton $fwY "Y" sDirection y
  tkm_MakeRadioButton $fwZ "Z" sDirection z

  # buttons. 
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { Translate $fDistance $sDirection }

  pack $fwDistance $fwDirection $fwX $fwY $fwZ $fwButtons \
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
    set wwDialog .wwLoadParcellationDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Parcellation" {-borderwidth 10}] } {

  set fwPath      $wwDialog.fwPath
  set fwColorFile $wwDialog.fwColorFile
  set fwButtons   $wwDialog.fwButtons

  # volume prompt and entry field
  tkm_MakeDirectorySelector $fwPath \
    "Enter a volume path" \
    sVolume

  # color prompt and entry field
  tkm_MakeFileSelector $fwColorFile \
    "Enter a color file name:" sColorFile

  # ok and cancel buttons. 
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { LoadParcellationVolume $sVolume $sColorFile }

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
    "Save RGB as:" sRGBName

  # ok and cancel buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { SaveRGB $sRGBName }

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
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

proc ShowParcellationLabel { ibShow } {

    global gfwVolValue gfwAuxVolValue gfwParcellationLabel

    if { $ibShow == 1 } {
  if {[catch { pack $gfwParcellationLabel  \
    -after $gfwAuxVolValue      \
    -side top                \
    -anchor w } sResult] } {
      if {[catch { pack $gfwParcellationLabel  \
        -after $gfwVolValue      \
        -side top                \
        -anchor w } sResult] } {
    puts "ShowParcellationLabel pack failed: $sResult"
    return;
      }
  }
    } else {
  pack forget $gfwParcellationLabel
    }
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

# ========================================================= BUILDING INTERFACE

proc CreateWindow { iwwTop } {

    global ksWindowName

    frame $iwwTop
    wm title . $ksWindowName
}

proc CreateMenuBar { ifwMenuBar } {

    global DspA_tOrientation_Sagittal DspA_tOrientation_Horizontal 
    global DspA_tOrientation_Coronal
    global DspA_tTool_Select DspA_tTool_Edit DspA_tTool_CtrlPts 
    global DspA_tTool_CustomEdit
    global gnVolX gnVolY gnVolZ
    global gb3D gbCursor gbSelection gbControlPoints gbFunctional
    global gbMainSurface gbOriginalSurface gbCanonicalSurface
    global gbDisplaySurfaceVertices gbInterpolateSurfaceVertices
    global gbParcellation gbAxes gbMaxIntProj
    global gTool
    global gDisplayMenu gnFunctionalOverlayMenuIndex gnTimeCourseMenuIndex

    frame $ifwMenuBar
    set mbwVolume      $ifwMenuBar.mbwVolume
    set mbwSelection   $ifwMenuBar.mbwSelection
    set mbwEdit        $ifwMenuBar.mbwEdit
    set mbwDisplay     $ifwMenuBar.mbwDisplay
    set mbwTools       $ifwMenuBar.mbwTools
    set mbwCtrlPts     $ifwMenuBar.mbwCtrlPts
    set mbwSurface     $ifwMenuBar.mbwSurface
    set mbwFunctional  $ifwMenuBar.mbwFunctional

    # volume menu button
    tkm_MakeMenu $mbwVolume "File" \
      { \
      \
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
      "Save Point" \
      SendCursor } \
      \
      { command \
      "Goto Point" \
      ReadCursor } \
      \
      { separator } \
      \
      { command \
      "Load Label..." \
      DoLoadLabelDlog } \
      \
      { command \
      "Save Label..." \
      DoSaveLabelDlog } \
      \
      { separator } \
      \
      { command \
      "Load Main Surface..." \
      "DoLoadSurfaceDlog $tkm_tSurfaceType_Current" } \
      \
      { command \
      "Load Original Surface..." \
      "DoLoadSurfaceDlog $tkm_tSurfaceType_Original" \
      tMenuGroup_SurfaceLoading } \
      \
      { command \
      "Load Pial Surface..." \
      "DoLoadSurfaceDlog $tkm_tSurfaceType_Canonical" \
      tMenuGroup_SurfaceLoading } \
      \
      { command \
      "Unload Surface" \
      UnloadSurface \
      SurfaceLoading } \
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
      { separator } \
      \
      { command \
      "Load Parcellation..." \
      DoLoadParcellationDlog } \
      \
      { separator } \
      \
      { command \
      "Save RGB..." \
      DoSaveRGBDlog } \
      \
      { separator } \
      \
      { command \
      "Quit" \
      QuitMedit } }
    
    # edit menu button
    tkm_MakeMenu $mbwEdit "Edit" \
      { \
      \
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
      "Clear Selection" \
      ClearSelection } \
      \
      { command \
      "Graph Selected Region" \
      GraphSelectedRegion } }

    # display menu
    tkm_MakeMenu $mbwDisplay "Display" \
      { { command \
      "Brightness / Contrast..." \
      DoVolumeColorScaleInfoDlog } \
      \
      { command \
      "Configure Functional Overlay..." \
      Overlay_DoConfigDlog \
      tMenuGroup_OverlayOptions } \
      \
      { command \
      "Configure Time Course Graph..." \
      TimeCourse_DoConfigDlog \
      tMenuGroup_TimeCourseOptions } \
       \
      { separator } \
       \
      { radio  \
      "Main Volume" \
      "SetDisplayFlag $DspA_tDisplayFlag_AuxVolume \
      $MWin_tVolumeType_Main"  \
      gDisplayedVolume  \
      0 } \
       \
      { radio  \
      "Aux Volume" \
      "SetDisplayFlag $DspA_tDisplayFlag_AuxVolume \
      $MWin_tVolumeType_Aux"  \
      gDisplayedVolume  \
      1 } \
       \
      { separator } \
       \
      { check  \
      "3D Multiple Views" \
      "SetDisplayConfig \$gb3D" \
      gb3D } \
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
      gbAxes } \
      \
      { check \
      "Maximum Intensity Projection" \
      "SetDisplayFlag $DspA_tDisplayFlag_MaxIntProj \
      \$gbMaxIntProj" \
      gbMaxIntProj } \
      \
      { check \
      "Selection / Label" \
      "SetDisplayFlag $DspA_tDisplayFlag_Selection \$gbSelection" \
      gbSelection } \
      \
      { check \
      "Control Points" \
      "SetDisplayFlag $DspA_tDisplayFlag_ControlPoints \
      \$gbControlPoints" \
      gbControlPoints } \
      \
      { check \
      "Functional Overlay" \
      "SetDisplayFlag $DspA_tDisplayFlag_FunctionalOverlay \
      \$gbFunctional" \
      gbFunctional } \
      \
      { check \
      "Parcellation Overlay" \
      "SetDisplayFlag $DspA_tDisplayFlag_ParcellationOverlay \
      \$gbParcellation" \
      gbParcellation } \
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
      gbDisplaySurfaceVertices } \
      \
      { check \
      "Interpolate Surface Vertices" \
      "SetDisplayFlag $DspA_tDisplayFlag_InterpolateSurfaceVertices \
      \$gbInterpolateSurfaceVertices"  \
      gbInterpolateSurfaceVertices }}

    # tools menu buttonx
    tkm_MakeMenu $mbwTools "Tools" \
      { { radio \
      "Select Voxels" \
      "SetTool $DspA_tTool_Select" \
      gTool \
      0 } \
      \
      { radio \
      "Edit Voxels" \
      "SetTool $DspA_tTool_Edit" \
      gTool \
      1 } \
      \
      { radio \
      "Select Ctrl Pts" \
      "SetTool $DspA_tTool_CtrlPts" \
      gTool \
      2 } \
      \
      { radio \
      "Edit Voxels w/ Custom Brush" \
      "SetTool $DspA_tTool_CustomEdit" \
      gTool \
      3 } \
      \
      { separator } \
      \
      { command \
      "Brush Size..." \
      DoBrushInfoDlog } \
      \
      { command \
      "Custom Brush..." \
      DoCustomBrushDlog } \
      \
      { separator } \
      \
      { command \
      "Threshold..." \
      DoThresholdDlog } \
      \
      { command \
      "Mirror" \
      Mirror } \
      \
      { command \
      "Translate..." \
      DoTranslateDlog } \
      \
      { command \
      "Rotate..." \
      DoRotateDlog } }
    
    # ctrl pts menu
    tkm_MakeMenu $mbwCtrlPts "Control Points" \
      { { command \
      "New" \
      { NewControlPoint $gnVolX $gnVolY $gnVolZ } } \
      \
      { command \
      "Select None" \
      DeselectAllControlPoints } \
      \
      { command \
      "Delete Selected" \
      DeleteSelectedControlPoints } \
      \
      { command \
      "Save" \
      WriteControlPointFile } }

    # surface menu
    tkm_MakeMenu $mbwSurface "Surface" \
      { { command  \
      "Show Nearest Main Vertex" \
      ShowNearestMainVertex }\
      \
      { command \
      "Show Nearest Original Vertex" \
      ShowNearestOriginalVertex } \
      \
      { command \
      "Show Nearest Pial Vertex" \
      ShowNearestCanonicalVertex } \
      \
      { command \
      "Find Main Vertex..." \
      { DoFindVertexDlog $tkm_tSurfaceType_Current } } \
      \
      { command \
      "Find Original Vertex..." \
      { DoFindVertexDlog $tkm_tSurfaceType_Original } } \
      \
      { command \
      "Find Pial Vertex..." \
      { DoFindVertexDlog $tkm_tSurfaceType_Canonical } } }
     
    pack $mbwVolume $mbwEdit $mbwDisplay $mbwTools \
      $mbwCtrlPts $mbwSurface  \
      -side left
}

proc CreateCursorFrame { ifwCursor } {

    global gfwVolCoords gfwRASCoords gfwTalCoords gfwVolValue
    global gfwAuxVolValue gfwFuncValue
    global gbLinkedCursor gsVolCoords gsRASCoords 
    global gsTalCoords gsVolName gnVolValue  gsFuncCoords
    global gsAuxVolValue gnAuxVolValue gfwFuncCoords gfFuncValue
    global gsParcellationLabel gfwParcellationLabel

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

    set gfwParcellationLabel $ifwCursor.fwParcellationLabel

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

    # the parcellation label
    tkm_MakeActiveLabel $gfwParcellationLabel "Parcellation: " gsParcellationLabel

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

    global DspA_tOrientation_Sagittal DspA_tOrientation_Horizontal 
    global DspA_tOrientation_Coronal
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
      "Sagittal (R -- L)" gOrientation $DspA_tOrientation_Sagittal \
      { SetOrientation $gOrientation }

    # sagittal slice slider
    tkm_MakeSlider $fwSagittalScale "" gnVolX 0 255 200 \
      { SetCursor $gnVolX $gnVolY $gnVolZ }       \
      1

    # horizontal orientaiton radio button
    tkm_MakeRadioButton $fwHorizontalButton \
      "Horizontal (Sup -- Inf)" gOrientation \
      $DspA_tOrientation_Horizontal \
      { SetOrientation $gOrientation }

    # horiontal slice slider
    tkm_MakeSlider $fwHorizontalScale "" gnVolY 0 255 200 \
      { SetCursor $gnVolX $gnVolY $gnVolZ }       \
      1

    # coronal orientaiton radio button
    tkm_MakeRadioButton $fwCoronalButton \
      "Coronal (Post -- Ant)" gOrientation \
      $DspA_tOrientation_Coronal \
      { SetOrientation $gOrientation }

    # coronal slice slider
    tkm_MakeSlider $fwCoronalScale "" gnVolZ 0 255 200 \
      { SetCursor $gnVolX $gnVolY $gnVolZ }       \
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

proc CreateToolbarFrame { ifwToolBar } {

    frame $ifwToolBar

    set fwButtons   $ifwToolBar.fwButtons

    tkm_MakeButtons $fwButtons \
      { {"Save Point" SendCursor}  \
      {"Goto Point" ReadCursor} }

    pack $fwButtons   \
      -side top \
      -fill x   \
      -expand yes
}

# = ====================================================================== MAIN

fixcolors

# build the window
set wwTop        .w
set fwMenuBar    $wwTop.fwMenuBar
set fwLeft       $wwTop.fwLeft
set fwRight      $wwTop.fwRight
set fwCursor     $fwLeft.fwCursor
set fwToolbar    $fwLeft.fwToolbar

CreateWindow         $wwTop

frame $fwLeft

CreateMenuBar        $fwMenuBar
CreateCursorFrame    $fwCursor
CreateToolbarFrame   $fwToolbar
CreateDisplayFrame   $fwRight

# pack the window
pack $fwMenuBar      \
  -side top    \
  -expand true \
  -fill x      \
  -anchor w

pack $fwCursor $fwToolbar \
  -side top         \
  -expand true      \
  -fill x           \
  -anchor w

pack $fwLeft $fwRight \
  -side left    \
  -padx 3       \
  -pady 3       \
  -expand true  \
  -fill x       \
  -fill y       \
  -anchor nw

pack $wwTop

proc CsurfInterface {} {

#    .w.fwMenuBar.mbwVolume.mwVolume delete 15 17
#    .w.fwMenuBar.mbwEdit.mwEdit delete 7
#    .w.fwMenuBar.mbwTools.mwTools delete 3
#    .w.fwMenuBar.mbwDisplay.mwDisplay delete 13
#    .w.fwMenuBar.mbwDisplay.mwDisplay delete 17 18
#    pack forget .w.fwMenuBar.mbwCtrlPts    
#    pack forget .w.fwMenuBar.mbwSurface

}


proc ErrorDlog { isMsg } {

    global gwwTop

    tk_messageBox -type ok \
      -icon error \
      -message $isMsg \
      -title "Error" \
      -parent $gwwTop
}

