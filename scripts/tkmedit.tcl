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
set kLabelFont   $pfont
set kNormalFont  $mfont
set kHighlightBgColor white

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
    global MWin_tVolumeType_Main MWin_tVolumeType_Aux

    global gDisplayedVolume gb3DDisplay gbCursor gbSelection gbControlPoints
    global gbFunctional gbMainSurface gbOriginalSurface gbCanonicalSurface
    global gbDisplaySurfaceVertices gbInterpolateSurfaceVertices
    global gbParcellation
    

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

proc DoLoadAuxVolumeDlog {} {

    global kLabelFont kNormalFont
    global gDialog

    set sVolumeName ""
    
    set wwDialog .wwLoadAuxVolumeDialog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Aux Volume" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set lwPrompt  $fwMain.lwPrompts
  set ewName    $fwMain.ewName

  set fwButtons $wwDialog.fwButtons
  set bwOK      $fwButtons.bwOK
  set bwCancel  $fwButtons.bwCancel

  # prompt and entry field
  frame $fwMain
  label $lwPrompt -text "Enter a volume name:" \
    -font $kLabelFont
  entry $ewName -textvariable sVolumeName \
    -font $kNormalFont
  pack $lwPrompt         \
    -side top      \
    -anchor w
  pack $ewName           \
    -side top      \
    -anchor w      \
    -expand yes    \
    -fill x

  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { LoadAuxVolume $sVolumeName; \
    Dialog_Close .wwLoadAuxVolumeDialog }
  button $bwCancel -text "Cancel" \
    -font $kNormalFont      \
    -command { Dialog_Close .wwLoadAuxVolumeDialog }
  pack $bwOK $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

  bind $wwDialog <Return> { .wwLoadAuxVolumeDialog.fwButtons.bwOK flash; .wwLoadAuxVolumeDialog.fwButtons.bwOK invoke }
   bind $wwDialog <Escape> { .wwLoadAuxVolumeDialog.fwButtons.bwCancel flash; .wwLoadAuxVolumeDialog.fwButtons.bwCancel invoke }
   }
}

proc DoSaveVolumeAsDlog {} {

    global kLabelFont kNormalFont
    global gDialog

    set sVolumeName ""
    set wwDialog .wwSaveVolumeAsDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save Volume As..." {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set lwPrompt  $fwMain.lwPrompt
  set ewName    $fwMain.ewName

  set fwButtons $wwDialog.fwButtons
  set bwOK      $fwButtons.bwOK
  set bwCancel  $fwButtons.bwCancel

  # prompt and entry field
  frame $fwMain
  label $lwPrompt -text "Path to save this volume in:" \
    -font $kLabelFont
  entry $ewName -textvariable sVolumeName \
    -font $kNormalFont
  pack $lwPrompt         \
    -side top      \
    -anchor w
  pack $ewName           \
    -side top      \
    -anchor w      \
    -expand yes    \
    -fill x

  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { SaveVolumeAs $sVolumeName; \
    Dialog_Close .wwSaveVolumeAsDlog }
  button $bwCancel -text "Cancel" \
    -font $kNormalFont      \
    -command { Dialog_Close .wwSaveVolumeAsDlog }
  pack $bwOK $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

  bind $wwDialog <Return> { .wwSaveVolumeAsDlog.fwButtons.bwOK flash; .wwSaveVolumeAsDlog.fwButtons.bwOK invoke }
  bind $wwDialog <Escape> { .wwSaveLabelDlog.fwButtons.bwCancel flash; .wwSaveLabelDlog.fwButtons.bwCancel invoke }
    }
}

proc DoSaveLabelDlog {} {

    global kLabelFont kNormalFont
    global gDialog

    set sLabelName ""
    set wwDialog .wwSaveLabelDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save Label" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set lwPrompt  $fwMain.lwPrompt
  set ewName    $fwMain.ewName

  set fwButtons $wwDialog.fwButtons
  set bwOK      $fwButtons.bwOK
  set bwCancel  $fwButtons.bwCancel

  # prompt and entry field
  frame $fwMain
  label $lwPrompt -text "Save this label as:" \
    -font $kLabelFont
  entry $ewName -textvariable sLabelName \
    -font $kNormalFont
  pack $lwPrompt         \
    -side top      \
    -anchor w
  pack $ewName           \
    -side top      \
    -anchor w      \
    -expand yes    \
    -fill x

  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { SaveLabel $sLabelName; \
    Dialog_Close .wwSaveLabelDlog }
  button $bwCancel -text "Cancel" \
    -font $kNormalFont      \
    -command { Dialog_Close .wwSaveLabelDlog }
  pack $bwOK $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

  bind $wwDialog <Return> { .wwSaveLabelDlog.fwButtons.bwOK flash; .wwSaveLabelDlog.fwButtons.bwOK invoke }
  bind $wwDialog <Escape> { .wwSaveLabelDlog.fwButtons.bwCancel flash; .wwSaveLabelDlog.fwButtons.bwCancel invoke }
    }
}

proc DoLoadLabelDlog {} {

    global kLabelFont kNormalFont
    global gDialog

    set sLabelName ""
    set wwDialog .wwLoadLabelDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Label" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set lwPrompt  $fwMain.lwPrompt
  set ewName    $fwMain.ewName

  set fwButtons $wwDialog.fwButtons
  set bwOK      $fwButtons.bwOK
  set bwCancel  $fwButtons.bwCancel

  # prompt and entry field
  frame $fwMain
  label $lwPrompt -text "Enter a label name:" \
    -font $kLabelFont
  entry $ewName -textvariable sLabelName \
    -font $kNormalFont
  pack $lwPrompt         \
    -side top      \
    -anchor w
  pack $ewName           \
    -side top      \
    -anchor w      \
    -expand yes    \
    -fill x

  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { LoadLabel $sLabelName; \
    Dialog_Close .wwLoadLabelDlog }
  button $bwCancel -text "Cancel" \
    -font $kNormalFont      \
    -command { Dialog_Close .wwLoadLabelDlog }
  pack $bwOK $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

  bind $wwDialog <Return> { .wwLoadLabelDlog.fwButtons.bwOK flash; .wwLoadLabelDlog.fwButtons.bwOK invoke }
  bind $wwDialog <Escape> { .wwLoadLabelDlog.fwButtons.bwCancel flash; .wwLoadLabelDlog.fwButtons.bwCancel invoke }
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

    global kLabelFont kNormalFont
    global gDialog
    global gLoadingSurface
 
    set gLoadingSurface $iSurface
    set sSurfaceName ""
    set wwDialog .wwLoadSurfaceDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Surface" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set lwPrompt  $fwMain.lwPrompt
  set ewName    $fwMain.ewName

  set fwButtons $wwDialog.fwButtons
  set bwOK      $fwButtons.bwOK
  set bwCancel  $fwButtons.bwCancel

  # prompt and entry field
  frame $fwMain
  label $lwPrompt -text "Enter a surface name:" \
    -font $kLabelFont
  entry $ewName -textvariable sSurfaceName \
    -font $kNormalFont
  pack $lwPrompt         \
    -side top      \
    -anchor w
  pack $ewName           \
    -side top      \
    -anchor w      \
    -expand yes    \
    -fill x

  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { LoadSurface $sSurfaceName; \
    Dialog_Close .wwLoadSurfaceDlog }
  button $bwCancel -text "Cancel" \
    -font $kNormalFont      \
    -command { Dialog_Close .wwLoadSurfaceDlog }
  pack $bwOK $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

  bind $wwDialog <Return> { .wwLoadSurfaceDlog.fwButtons.bwOK flash; .wwLoadSurfaceDlog.fwButtons.bwOK invoke }
  bind $wwDialog <Escape> { .wwLoadSurfaceDlog.fwButtons.bwCancel flash; .wwLoadSurfaceDlog.fwButtons.bwCancel invoke }
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

    global kLabelFont kNormalFont
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
  frame $fwMain
  label $lwPrompt -text "Enter a vertex number:" \
    -font $kLabelFont
  entry $ewName -textvariable nVertex \
    -font $kNormalFont
  pack $lwPrompt         \
    -side top      \
    -anchor w
  pack $ewName           \
    -side top      \
    -anchor w      \
    -expand yes    \
    -fill x

  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { FindVertex $nVertex; \
    Dialog_Close .wwFindVertexDlog }
  button $bwCancel -text "Cancel" \
    -font $kNormalFont      \
    -command { Dialog_Close .wwFindVertexDlog }
  pack $bwOK $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

  bind $wwDialog <Return> { .wwFindVertexDlog.fwButtons.bwOK flash; .wwFindVertexDlog.fwButtons.bwOK invoke }
  bind $wwDialog <Escape> { .wwFindVertexDlog.fwButtons.bwCancel flash; .wwFindVertexDlog.fwButtons.bwCancel invoke }
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

    global kLabelFont kNormalFont
    global tFunctionalVolume_Overlay tFunctionalVolume_TimeCourse
    global gDialog
    global gLoadingVolume

    set gLoadingVolume $iVolume
    set sPath ""
    set sStem ""
    set wwDialog .wwLoadFunctionalDlog

    # this doesn't work the second+ time the dlog is loaded
#    if { $iVolume == $tFunctionalVolume_Overlay } {
#  set sWindowName "Load Overlay Data"
#    } else {
#  set sWindowName "Load Time Course Data"
#    }
    set sWindowName "Load Functional Volume"

    # try to create the dlog...
    if { [Dialog_Create $wwDialog $sWindowName {-borderwidth 10}] } {

#  wm title $wwDialog $sWindowName

  set fwMain    $wwDialog.fwMain
  set lwPrompt  $fwMain.lwPrompt
  set ewPath    $fwMain.ewPath
  set ewStem    $fwMain.ewStem

  set fwButtons $wwDialog.fwButtons
  set bwOK      $fwButtons.bwOK
  set bwCancel  $fwButtons.bwCancel

  # prompt and entry fields
  frame $fwMain
  label $lwPrompt -text "Enter a path and stem:" \
    -font $kLabelFont
  entry $ewPath -textvariable sPath \
    -font $kNormalFont
  entry $ewStem -textvariable sStem \
    -font $kNormalFont        \
    -width 8
  pack $lwPrompt         \
    -side top      \
    -anchor w
  pack $ewPath           \
    -side left     \
    -anchor w      \
    -expand yes    \
    -fill x
  pack $ewStem           \
    -side left     \
    -anchor e      \
    -expand no

  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { LoadFunctionalVolume $sPath $sStem; \
    Dialog_Close .wwLoadFunctionalDlog }
  button $bwCancel -text "Cancel" \
    -font $kNormalFont      \
    -command { Dialog_Close .wwLoadFunctionalDlog }
  pack $bwOK $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

   bind $wwDialog <Return> { .wwLoadFunctionalDlog.fwButtons.bwOK flash; .wwLoadFunctionalDlog.fwButtons.bwOK invoke }
   bind $wwDialog <Escape> { .wwFindVertexDlog.fwButtons.bwCancel flash; .wwFindVertexDlog.fwButtons.bwCancel invoke }
   }


}

proc DoSaveDlog {} {

    global kLabelFont kNormalFont
    global gDialog

    set wwDialog .wwSaveDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save Volume" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set lwPrompt  $fwMain.lwPrompt

  set fwButtons $wwDialog.fwButtons
  set bwOK      $fwButtons.bwOK
  set bwCancel  $fwButtons.bwCancel

  # prompt and entry field
  frame $fwMain
  label $lwPrompt -text "Are you sure you wish to save changes to the volume?" \
    -font $kLabelFont
  pack $lwPrompt         \
    -side top      \
    -anchor w

  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { SaveVolume; \
    Dialog_Close .wwSaveDlog; }
  button $bwCancel -text "Cancel" \
    -font $kNormalFont      \
    -command { Dialog_Close .wwSaveDlog }
  pack $bwOK $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

   bind $wwDialog <Return> { .wwSaveDlog.fwButtons.bwOK flash; .wwSaveDlog.fwButtons.bwOK invoke }
   bind $wwDialog <Escape> { .wwSaveDlog.fwButtons.bwCancel flash; .wwSaveDlog.fwButtons.bwCancel invoke }
    }
}

proc DoBrushInfoDlog {} {

    global kLabelFont kNormalFont kHighlightBgColor
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

  set fwLabel              $wwDialog.fwLabel
  set lwLabel              $fwLabel.lwLabel
  
  set fwRadiusLabel        $wwDialog.fwRadiusLabel
  set lwRadius             $fwRadiusLabel.lwRadius
    
  set fwRadiusScale        $wwDialog.fwRadiusScale
  set swRadius             $fwRadiusScale.swRadius
  set ewRadius             $fwRadiusScale.ewRadius
  
  set fwShapeRadios        $wwDialog.fwShapeRadios
  set rwCircle             $fwShapeRadios.rwCircle
  set rwSquare             $fwShapeRadios.rwSquare
  
  set fw3DCheckbox         $wwDialog.fw3DCheckbox
  set cw3D                 $fw3DCheckbox.cw3D
    
  set fwButtons            $wwDialog.fwButtons
  set bwOK                 $fwButtons.bwOK
  set bwApply              $fwButtons.bwApply
  set bwCancel             $fwButtons.bwCancel
  
  # the label that goes at the top of the frame
  frame $fwLabel
  label $lwLabel -text "Brush" \
    -font $kLabelFont
  pack $lwLabel -side left \
    -anchor w
  
  # Radius label
  frame $fwRadiusLabel
  label $lwRadius -text "Radius" \
    -font $kNormalFont
  pack $lwRadius -side left \
    -anchor w
  
  # radius slider
  frame $fwRadiusScale
  scale $swRadius -orient horizontal          \
    -variable gnBrushRadius             \
    -from 1                             \
    -to 100                             \
    -showvalue false 
  entry $ewRadius -textvariable gnBrushRadius \
    -width 3
  pack $swRadius      \
    -side left  \
    -anchor w   \
    -fill x
  pack $ewRadius     \
    -side left \
    -anchor w  \
    -expand no
  
  # sahpe radio buttons
  frame $fwShapeRadios
  radiobutton $rwCircle -text "Circle"  \
    -variable gBrushShape         \
    -font $kNormalFont   \
    -highlightbackground $kHighlightBgColor \
    -relief flat                  \
    -value $DspA_tBrushShape_Circle
  radiobutton $rwSquare -text "Square " \
    -variable gBrushShape         \
    -font $kNormalFont   \
    -highlightbackground $kHighlightBgColor \
    -relief flat                  \
    -value $DspA_tBrushShape_Square
  pack $rwCircle $rwSquare              \
    -side top                     \
    -anchor w
  
  # 3d checkbox
  frame $fw3DCheckbox
  checkbutton $cw3D -text "3D"    \
    -variable gbBrush3D     \
    -font $kNormalFont      \
    -highlightbackground $kHighlightBgColor
  pack $cw3D        \
    -side top \
    -anchor w
  
  # pack them in a column
  pack $fwLabel $fwRadiusLabel $fwRadiusScale \
    $fwShapeRadios $fw3DCheckbox        \
    -side top                           \
    -anchor w                           \
    -expand yes                         \
    -fill x
  
  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { SetBrush $gnBrushRadius $gBrushShape $gbBrush3D; \
    Dialog_Close .wwBrushInfoDlog }
  button $bwApply -text "Apply"     \
    -font $kNormalFont  \
    -command { SetBrush $gnBrushRadius $gBrushShape $gbBrush3D; }
  button $bwCancel -text "Cancel" \
    -font $kNormalFont      \
    -command { SetBrush $gnSavedBrushRadius $gSavedBrushShape $gbSavedBrush3D; \
    Dialog_Close .wwBrushInfoDlog }
  pack $bwOK $bwApply $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5

  pack $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

    bind $wwDialog <Return> { .wwBrushInfoDlog.fwButtons.bwOK flash; .wwBrushInfoDlog.fwButtons.bwOK invoke }
    bind $wwDialog <Escape> { .wwBrushInfoDlog.fwButtons.bwCancel flash; .wwBrushInfoDlog.fwButtons.bwCancel invoke }
   }
}

proc DoCustomBrushDlog {} {

    global kLabelFont kNormalFont kHighlightBgColor
    global gDialog
    global gnBrushLow gnBrushHigh gnBrushNewValue
    global gnSavedBrushLow gnSavedBrushHigh gnSavedBrushNewValue

    set wwDialog .wwCustomBrushDlog

    set gnSavedBrushLow      $gnBrushLow
    set gnSavedBrushHigh     $gnBrushHigh
    set gnSavedBrushNewValue $gnBrushNewValue

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Custom Brush" {-borderwidth 10}] } {

  set fwLabel              $wwDialog.fwLabel
  set lwLabel              $fwLabel.lwLabel
  
  set fwLowLabel        $wwDialog.fwLowLabel
  set lwLow             $fwLowLabel.lwLow

  set fwLowScale        $wwDialog.fwLowScale
  set swLow             $fwLowScale.swLow
  set ewLow             $fwLowScale.ewLow
  
  set fwHighLabel        $wwDialog.fwHighLabel
  set lwHigh             $fwHighLabel.lwHigh

  set fwHighScale        $wwDialog.fwHighScale
  set swHigh             $fwHighScale.swHigh
  set ewHigh             $fwHighScale.ewHigh
  
  set fwNewValueLabel        $wwDialog.fwNewValueLabel
  set lwNewValue             $fwNewValueLabel.lwNewValue

  set fwNewValueScale        $wwDialog.fwNewValueScale
  set swNewValue             $fwNewValueScale.swNewValue
  set ewNewValue             $fwNewValueScale.ewNewValue
  
  set fwButtons            $wwDialog.fwButtons
  set bwOK                 $fwButtons.bwOK
  set bwApply              $fwButtons.bwApply
  set bwCancel             $fwButtons.bwCancel
  
  # the label that goes at the top of the frame
  frame $fwLabel
  label $lwLabel -text "Custom Brush" \
    -font $kLabelFont
  pack $lwLabel -side left \
    -anchor w
  
  # label
  frame $fwLowLabel
  label $lwLow -text "Low" \
    -font $kNormalFont
  pack $lwLow -side left \
    -anchor w
  
  # scale
  frame $fwLowScale
  scale $swLow -orient horizontal       \
    -variable gnBrushLow          \
    -from 0                             \
    -to 255                             \
    -showvalue false 
  entry $ewLow -textvariable gnBrushLow \
    -width 3
  pack $swLow   \
    -side left  \
    -anchor w   \
    -expand yes \
    -fill x
  pack $ewLow  \
    -side left \
    -anchor w  \
    -expand no

  # label
  frame $fwHighLabel
  label $lwHigh -text "High" \
    -font $kNormalFont
  pack $lwHigh -side left \
    -anchor w
  
  # scale
  frame $fwHighScale
  scale $swHigh -orient horizontal       \
    -variable gnBrushHigh          \
    -from 0                             \
    -to 255                             \
    -showvalue false 
  entry $ewHigh -textvariable gnBrushHigh \
    -width 3
  pack $swHigh   \
    -side left  \
    -anchor w   \
    -expand yes \
    -fill x
  pack $ewHigh  \
    -side left \
    -anchor w  \
    -expand no

  # label
  frame $fwNewValueLabel
  label $lwNewValue -text "New Value" \
    -font $kNormalFont
  pack $lwNewValue -side left \
    -anchor w
  
  # scale
  frame $fwNewValueScale
  scale $swNewValue -orient horizontal       \
    -variable gnBrushNewValue          \
    -from 0                             \
    -to 255                             \
    -showvalue false 
  entry $ewNewValue -textvariable gnBrushNewValue \
    -width 3
  pack $swNewValue   \
    -side left  \
    -anchor w   \
    -expand yes \
    -fill x
  pack $ewNewValue  \
    -side left \
    -anchor w  \
    -expand no

  # pack them in a column
  pack $fwLabel $fwLowLabel $fwLowScale \
    $fwHighLabel $fwHighScale \
    $fwNewValueLabel $fwNewValueScale \
    -side top                           \
    -anchor w                           \
    -expand yes                         \
    -fill x
  
  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { SetBrushThreshold $gnBrushLow $gnBrushHigh $gnBrushNewValue; \
    Dialog_Close .wwCustomBrushDlog }
  button $bwApply -text "Apply"     \
    -font $kNormalFont  \
    -command { SetBrushThreshold $gnBrushLow $gnBrushHigh $gnBrushNewValue; }
  button $bwCancel -text "Cancel" \
    -font $kNormalFont      \
    -command { SetBrushThreshold $gnSavedBrushLow $gnSavedBrushHigh $gnSavedBrushNewValue; \
    Dialog_Close .wwCustomBrushDlog }
  pack $bwOK $bwApply $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5

  pack $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

    bind $wwDialog <Return> { .wwCustomBrushDlog.fwButtons.bwOK flash; .wwCustomBrushDlog.fwButtons.bwOK invoke }
    bind $wwDialog <space> { .wwCustomBrushDlog.fwButtons.bwApply flash; .wwCustomBrushDlog.fwButtons.bwApply invoke }
    bind $wwDialog <Escape> { .wwCustomBrushDlog.fwButtons.bwCancel flash; .wwCustomBrushDlog.fwButtons.bwCancel invoke }
   }
}

proc DoVolumeColorScaleInfoDlog {} {

    global kLabelFont kNormalFont
    global gDialog
    global gfVolumeColorScaleThresh gfVolumeColorScaleSquash 
 
    set wwDialog .wwVolumeColorScaleInfoDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Brightness / Contrast" {-borderwidth 10}] } {

  set fwTop                $wwDialog.fwTop

  set fwLeft               $fwTop.fwLeft
  set lwThresh             $fwLeft.lwThresh
  set lwSquash             $fwLeft.lwSquash
  
  set fwRight              $fwTop.fwRight
  set fwScales             $fwRight.fwScales
  set swSquash             $fwScales.swSquash
  set swThresh             $fwScales.swThresh

  set fwEntries            $fwRight.fwEntries
  set ewSquash             $fwEntries.ewSquash
  set ewThresh             $fwEntries.ewThresh
  
  set fwButtons            $wwDialog.fwButtons
  set bwOK                 $fwButtons.bwOK
  
  frame $fwTop
  
  # labels
  frame $fwLeft
  label $lwThresh -text "Brightness" \
    -font $kNormalFont
  label $lwSquash -text "Contrast" \
    -font $kNormalFont
  pack $lwThresh $lwSquash \
    -side top        \
    -pady 5          \
    -anchor w

  frame $fwRight

  # scales
  frame $fwScales
  scale $swThresh -orient horizontal         \
    -command { SetVolumeColorScale $gfVolumeColorScaleThresh $gfVolumeColorScaleSquash } \
    -variable gfVolumeColorScaleThresh \
    -from 1                            \
    -to 0                              \
    -length 100                        \
    -resolution 0.01                   \
    -showvalue false
  scale $swSquash -orient horizontal         \
    -command { SetVolumeColorScale $gfVolumeColorScaleThresh $gfVolumeColorScaleSquash } \
    -variable gfVolumeColorScaleSquash \
    -from 0                            \
    -to 20                             \
    -length 100                        \
    -resolution 1                      \
    -showvalue false
  pack $swThresh $swSquash \
    -side top        \
    -pady 5          \
    -expand yes      \
    -anchor w

  # entries
  frame $fwEntries
  entry $ewThresh                                \
    -textvariable gfVolumeColorScaleThresh \
    -font $kNormalFont                     \
    -width 5
  bind $ewThresh <Return> { SetVolumeColorScale $gfVolumeColorScaleThresh $gfVolumeColorScaleSquash }
  entry $ewSquash                                \
    -textvariable gfVolumeColorScaleSquash \
    -font $kNormalFont                     \
    -width 5
  bind $ewSquash <Return> { SetVolumeColorScale $gfVolumeColorScaleThresh $gfVolumeColorScaleSquash }
  pack $ewThresh $ewSquash \
    -side top        \
    -pady 5          \
    -expand yes      \
    -anchor w
  
  # pack scales and entries
  pack $fwScales $fwEntries  \
    -side left         \
    -padx 5            \
    -anchor w          \
    -expand yes        \
    -fill x

  # pack left and right
  pack $fwLeft $fwRight    \
    -side left       \
    -padx 5          \
    -expand yes      \
    -anchor w

  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { Dialog_Close .wwVolumeColorScaleInfoDlog }
  pack $bwOK              \
    -side right

  pack $fwTop $fwButtons  \
    -side top    \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

    bind $wwDialog <Return> { .wwVolumeColorScaleInfoDlog.fwButtons.bwOK flash; .wwVolumeColorScaleInfoDlog.fwButtons.bwOK invoke }
    }
}

proc DoThresholdDlog {} {

    global kLabelFont kNormalFont kHighlightBgColor
    global gDialog
    
    set wwDialog .wwThresholdDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Threshold" {-borderwidth 10}] } {

  set fwTop                $wwDialog.fwTop

  set fwThresh      $wwDialog.fwThresh
  set lwThresh      $fwThresh.lwThresh
  set ewThresh      $fwThresh.ewThresh
  
  set fwDirection   $wwDialog.fwDirection
  set rwAbove       $fwDirection.rwAbove
  set rwBelow       $fwDirection.rwBelow

  set fwNewValue    $wwDialog.fwNewValue
  set lwNewValue    $fwNewValue.lwNewValue
  set ewNewValue    $fwNewValue.ewNewValue

  set fwButtons     $wwDialog.fwButtons
  set bwCancel      $fwButtons.bwCancel
  set bwApply       $fwButtons.bwApply
  set bwOK          $fwButtons.bwOK
  
  # threshold value
  frame $fwThresh
  label $lwThresh -text "Threshold:" \
    -font $kNormalFont
  entry $ewThresh -textvariable nThreshold \
    -font $kNormalFont
  pack $lwThresh     \
    -side left      \
    -anchor w
  pack  $ewThresh     \
    -side right      \
    -anchor e
  
  #direction radios
  frame $fwDirection
  radiobutton $rwBelow -text "Below"              \
    -variable bAbove                        \
    -font $kNormalFont                      \
    -highlightbackground $kHighlightBgColor \
    -relief flat                            \
    -value 0
  radiobutton $rwAbove -text "Above"              \
    -variable bAbove                        \
    -font $kNormalFont                      \
    -highlightbackground $kHighlightBgColor \
    -relief flat                            \
    -value 1
  pack $rwBelow $rwAbove \
    -side left     \
    -anchor w

  # new value
  frame $fwNewValue
  label $lwNewValue -text "New Value:" \
    -font $kNormalFont
  entry $ewNewValue -textvariable nNewValue \
    -font $kNormalFont
  pack $lwNewValue     \
    -side left      \
    -anchor w
  pack  $ewNewValue     \
    -side right      \
    -anchor e

  # buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { ThresholdVolume $nThreshold $bAbove $nNewValue; \
    Dialog_Close .wwThresholdDlog }
  button $bwApply -text "Apply"  \
    -font $kNormalFont     \
    -command { ThresholdVolume $nThreshold $bAbove $nNewValue }
  button $bwCancel -text "Close" \
    -font $kNormalFont      \
    -command { Dialog_Close .wwThresholdDlog }
  pack $bwOK $bwApply $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5
  
  pack $fwThresh $fwDirection $fwNewValue $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

    bind $wwDialog <Return> { .wwThresholdDlog.fwButtons.bwOK flash; .wwThresholdDlog.fwButtons.bwOK invoke }
    bind $wwDialog <space> { .wwThresholdDlog.fwButtons.bwApply flash; .wwThresholdDlog.fwButtons.bwApply invoke }
   bind $wwDialog <Escape> { .wwThresholdDlog.fwButtons.bwCancel flash; .wwThresholdDlog.fwButtons.bwCancel invoke }    }
}

proc DoLoadParcellationDlog {} {

    global kLabelFont kNormalFont
    global gDialog

    set sVolume ""
    set sColorFile ""
    set wwDialog .wwLoadParcellationDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Parcellation" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set lwVolume  $fwMain.lwVolume
  set ewVolume  $fwMain.ewVolume
  set lwColor   $fwMain.lwColor
  set ewColor   $fwMain.ewColor

  set fwButtons $wwDialog.fwButtons
  set bwOK      $fwButtons.bwOK
  set bwCancel  $fwButtons.bwCancel

  # volume prompt and entry field
  frame $fwMain
  label $lwVolume -text "Enter a volume path with prefix (i.e. /COR-):" \
    -font $kLabelFont
  entry $ewVolume -textvariable sVolume \
    -font $kNormalFont
  pack $lwVolume         \
    -side top      \
    -anchor w
  pack $ewVolume         \
    -side top      \
    -anchor w      \
    -expand yes    \
    -fill x

  # color prompt and entry field
  label $lwColor -text "Enter a color file name:" \
    -font $kLabelFont
  entry $ewColor -textvariable sColorFile \
    -font $kNormalFont
  pack $lwColor          \
    -side top      \
    -anchor w
  pack $ewColor          \
    -side top      \
    -anchor w      \
    -expand yes    \
    -fill x

  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { LoadParcellationVolume $sVolume $sColorFile; \
    Dialog_Close .wwLoadParcellationDlog }
  button $bwCancel -text "Cancel" \
    -font $kNormalFont      \
    -command { Dialog_Close .wwLoadParcellationDlog }
  pack $bwOK $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

    bind $wwDialog <Return> { .wwLoadParcellationDlog.fwButtons.bwOK invoke }
    bind $wwDialog <Escape> { .wwLoadParcellationDlog.fwButtons.bwCancel invoke }
    }
}

proc DoSaveRGBDlog {} {

    global kLabelFont kNormalFont
    global gDialog

    set sRGBName ""
    set wwDialog .wwSaveRGBDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save RGB" {-borderwidth 10}] } {

  set fwMain    $wwDialog.fwMain
  set lwPrompt  $fwMain.lwPrompt
  set ewName    $fwMain.ewName

  set fwButtons $wwDialog.fwButtons
  set bwOK      $fwButtons.bwOK
  set bwCancel  $fwButtons.bwCancel

  # prompt and entry field
  frame $fwMain
  label $lwPrompt -text "Save RGB as:" \
    -font $kLabelFont
  entry $ewName -textvariable sRGBName \
    -font $kNormalFont
  pack $lwPrompt         \
    -side top      \
    -anchor w
  pack $ewName           \
    -side top      \
    -anchor w      \
    -expand yes    \
    -fill x

  # ok and cancel buttons. the ok button will actually call the
  # important function here.
  frame $fwButtons
  button $bwOK -text "OK"     \
    -font $kNormalFont  \
    -command { Dialog_Close .wwSaveRGBDlog; \
               SaveRGB $sRGBName; }
  button $bwCancel -text "Cancel" \
    -font $kNormalFont      \
    -command { Dialog_Close .wwSaveRGBDlog }
  pack $bwOK $bwCancel \
    -side right  \
    -padx 5      \
    -pady 5

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5

    bind $wwDialog <Return> { .wwSaveRGBDlog.fwButtons.bwOK flash; .wwSaveRGBDlog.fwButtons.bwOK invoke }
    bind $wwDialog <Escape> { .wwSaveRGBDlog.fwButtons.bwCancel flash; .wwSaveRGBDlog.fwButtons.bwCancel invoke }
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
  pack forget $gfwTalCoords
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
  pack forget $gfwTalCoords
    }
}

proc ShowOverlayOptions { ibShow } {
    global gDisplayMenu gnFunctionalOverlayMenuIndex

    #  enable or disable menu item
    if { $ibShow == 0 } {
  $gDisplayMenu entryconfigure $gnFunctionalOverlayMenuIndex \
    -state disabled
    } else {
  $gDisplayMenu entryconfigure $gnFunctionalOverlayMenuIndex \
    -state normal
    }
}

proc ShowTimeCourseOptions { ibShow } {
    global gDisplayMenu gnTimeCourseMenuIndex

    # enable or disable menu item
    if { $ibShow == 0 } {
  $gDisplayMenu entryconfigure $gnTimeCourseMenuIndex \
    -state disabled
    } else {
  $gDisplayMenu entryconfigure $gnTimeCourseMenuIndex \
    -state normal
    }
}

# ========================================================= BUILDING INTERFACE

proc CreateWindow { iwwTop } {

    global ksWindowName

    frame $iwwTop
    wm title . $ksWindowName
}

proc CreateMenuBar { ifwMenuBar } {

    global kLabelFont kNormalFont
    global DspA_tOrientation_Sagittal DspA_tOrientation_Horizontal 
    global DspA_tOrientation_Coronal
    global DspA_tTool_Select DspA_tTool_Edit DspA_tTool_CtrlPts 
    global DspA_tTool_CustomEdit
    global gnVolX gnVolY gnVolZ
    global gb3D gbCursor gbSelection gbControlPoints gbFunctional
    global gbMainSurface gbOriginalSurface gbCanonicalSurface
    global gbDisplaySurfaceVertices gbInterpolateSurfaceVertices
    global gbParcellation
    global gTool
    global gDisplayMenu gnFunctionalOverlayMenuIndex gnTimeCourseMenuIndex


    frame $ifwMenuBar

    set mbwVolume      $ifwMenuBar.mbwVolume
    set mwVolume       $mbwVolume.mwVolume

    set mbwSelection   $ifwMenuBar.mbwSelection
    set mwSelection    $mbwSelection.mwSelection

    set mbwEdit        $ifwMenuBar.mbwEdit
    set mwEdit         $mbwEdit.mwEdit

    set mbwDisplay     $ifwMenuBar.mbwDisplay
    set mwDisplay      $mbwDisplay.mwDisplay

    set mbwTools       $ifwMenuBar.mbwTools
    set mwTools        $mbwTools.mwTools

    set mbwCtrlPts     $ifwMenuBar.mbwCtrlPts
    set mwCtrlPts      $mbwCtrlPts.mwCtrlPts

    set mbwSurface     $ifwMenuBar.mbwSurface
    set mwSurface      $mbwSurface.mwSurface
 
    set mbwFunctional  $ifwMenuBar.mbwFunctional
    set mwFunctional   $mbwFunctional.mwFunctional

    # volume menu button
    menubutton $mbwVolume -text "File" \
      -menu $mwVolume              \
      -underline 0                 \
      -font $kNormalFont
    
    # volume menu items
    menu $mwVolume
    $mwVolume add command -label "Load Aux Volume..." \
      -command { DoLoadAuxVolumeDlog }          \
      -underline 1                              \
      -font $kNormalFont
    $mwVolume add command -label "Save Volume" \
      -command { DoSaveDlog }            \
      -underline 0                       \
      -font $kNormalFont
    $mwVolume add command -label "Save Volume As..." \
      -command { DoSaveVolumeAsDlog }          \
      -font $kNormalFont
    $mwVolume add separator
    $mwVolume add command -label "Save Point" \
      -command { SendCursor }                         \
      -underline 20                                   \
      -font $kNormalFont
    $mwVolume add command -label "Goto Point" \
      -command { ReadCursor }                           \
      -underline 0                                      \
      -font $kNormalFont
    $mwVolume add separator
    $mwVolume add command -label "Load Label..." \
      -command { DoLoadLabelDlog }            \
      -underline 0                            \
      -font $kNormalFont
    $mwVolume add command -label "Save Label..." \
      -command { DoSaveLabelDlog }            \
      -underline 2                            \
      -font $kNormalFont
    $mwVolume add separator
    $mwVolume add command -label "Load Surface..." \
      -command { DoLoadSurfaceDlog $tkm_tSurfaceType_Current } \
      -underline 5                            \
      -font $kNormalFont
    $mwVolume add command -label "Unload Surface" \
      -command { UnloadSurface }             \
      -underline 0                           \
      -font $kNormalFont
    $mwVolume add separator
    $mwVolume add command -label "Load Overlay Data..." \
      -command { DoLoadFunctionalDlog $tFunctionalVolume_Overlay } \
      -underline 7                                    \
      -font $kNormalFont                              
    $mwVolume add command -label "Load Time Course Data..." \
      -command { DoLoadFunctionalDlog $tFunctionalVolume_TimeCourse } \
      -underline 5                                        \
      -font $kNormalFont
    $mwVolume add separator
    $mwVolume add command -label "Load Parcellation..." \
      -command { DoLoadParcellationDlog }         \
      -font $kNormalFont
    $mwVolume add separator
    $mwVolume add command -label "Save RGB..." \
      -command { DoSaveRGBDlog }         \
      -font $kNormalFont
    $mwVolume add separator
    $mwVolume add command -label "Quit" \
      -command { QuitMedit }      \
      -underline 0                \
      -font $kNormalFont

    # edit menu button
    menubutton $mbwEdit -text "Edit" \
      -menu $mwEdit            \
      -underline 0             \
      -font $kNormalFont
    
    # edit menu items
    menu $mwEdit
    $mwEdit add command -label "Undo Last Edit" \
      -command { UndoLastEdit }           \
      -font $kNormalFont
    $mwEdit add separator
    $mwEdit add command -label "Take Snapshot of Volume" \
      -command { SnapshotVolume }                  \
      -font $kNormalFont
    $mwEdit add command -label "Restore Volume from Snapshot" \
      -command { RestoreVolumeFromSnapshot }            \
      -font $kNormalFont
    $mwEdit add separator
    $mwEdit add command -label "Clear Selection" \
      -command { ClearSelection }          \
      -font $kNormalFont
    $mwEdit add command -label "Graph Selected Region" \
      -command { GraphSelectedRegion }                 \
      -underline 4                                     \
      -font $kNormalFont
  
    # display menu button
    menubutton $mbwDisplay -text "Display" \
      -menu $mwDisplay         \
      -underline 0             \
      -font $kNormalFont
    
    # display items
    menu $mwDisplay
    $mwDisplay add command -label "Brightness / Contrast..." \
      -command { DoVolumeColorScaleInfoDlog }       \
      -font $kNormalFont
    $mwDisplay add command -label "Configure Functional Overlay..." \
      -command { Overlay_DoConfigDlog }       \
      -font $kNormalFont                 \
      -state disabled
    $mwDisplay add command -label "Configure Time Course Graph..." \
      -command { TimeCourse_DoConfigDlog }       \
      -font $kNormalFont                 \
      -state disabled
    # save the location of this menu item so we can enable it later
    set gDisplayMenu                 $mwDisplay
    set gnFunctionalOverlayMenuIndex 2
    set gnTimeCourseMenuIndex 3
    $mwDisplay add separator
    $mwDisplay add radio -label "Main Volume" \
      -command { SetDisplayFlag $DspA_tDisplayFlag_AuxVolume $MWin_tVolumeType_Main } \
      -variable gDisplayedVolume        \
      -font $kNormalFont                \
      -value 0
    $mwDisplay add radio -label "Aux Volume" \
      -command { SetDisplayFlag $DspA_tDisplayFlag_AuxVolume $MWin_tVolumeType_Aux } \
      -variable gDisplayedVolume       \
      -font $kNormalFont               \
      -value 1
    $mwDisplay add separator
    $mwDisplay add check -label "3D Multiple Views" \
      -command { SetDisplayConfig $gb3D }     \
      -underline 0                            \
      -variable gb3D                          \
      -font $kNormalFont
    $mwDisplay add check -label "Cursor" \
      -command { SetDisplayFlag $DspA_tDisplayFlag_Cursor $gbCursor } \
      -underline 0                 \
      -variable gbCursor           \
      -font $kNormalFont 
    $mwDisplay add check -label "Selection / Label" \
      -command { SetDisplayFlag $DspA_tDisplayFlag_Selection $gbSelection } \
      -underline 12                           \
      -variable gbSelection                   \
      -font $kNormalFont
    $mwDisplay add check -label "Control Points" \
      -command { SetDisplayFlag $DspA_tDisplayFlag_ControlPoints $gbControlPoints } \
      -underline 3                         \
      -variable gbControlPoints            \
      -font $kNormalFont
    $mwDisplay add check -label "Functional Overlay" \
      -command { SetDisplayFlag $DspA_tDisplayFlag_FunctionalOverlay $gbFunctional } \
      -variable gbFunctional                   \
      -underline 0                             \
      -font $kNormalFont
    $mwDisplay add check -label "Parcellation Overlay" \
      -command { SetDisplayFlag $DspA_tDisplayFlag_ParcellationOverlay $gbParcellation } \
      -variable gbParcellation                 \
      -underline 2                             \
      -font $kNormalFont
    $mwDisplay add check -label "Main Surface" \
      -command { SetDisplayFlag $DspA_tDisplayFlag_MainSurface $gbMainSurface } \
      -underline 0                       \
      -variable gbMainSurface            \
      -font $kNormalFont
    $mwDisplay add check -label "Original Surface" \
      -command { SetDisplayFlag $DspA_tDisplayFlag_OriginalSurface $gbOriginalSurface } \
      -underline 0                           \
      -variable gbOriginalSurface            \
      -font $kNormalFont
    $mwDisplay add check -label "Pial Surface" \
      -command { SetDisplayFlag $DspA_tDisplayFlag_CanonicalSurface $gbCanonicalSurface } \
      -underline 0                       \
      -variable gbCanonicalSurface       \
      -font $kNormalFont
    $mwDisplay add check -label "Surface Vertices" \
      -command { SetDisplayFlag $DspA_tDisplayFlag_DisplaySurfaceVertices $gbDisplaySurfaceVertices } \
      -underline 8                           \
      -variable gbDisplaySurfaceVertices     \
      -font $kNormalFont
    $mwDisplay add check -label "Interpolate Surface Vertices" \
      -command { SetDisplayFlag $DspA_tDisplayFlag_InterpolateSurfaceVertices $gbInterpolateSurfaceVertices } \
      -underline 0                                       \
      -variable gbInterpolateSurfaceVertices             \
      -font $kNormalFont

    # tools menu buttonx
    menubutton $mbwTools -text "Tools" \
      -menu $mwTools             \
      -underline 0               \
      -font $kNormalFont
    
    # tool menu items
    menu $mwTools
    $mwTools add radio -label "Select Voxels"       \
      -command { SetTool $DspA_tTool_Select } \
      -variable gTool                         \
      -font $kNormalFont                      \
      -value $DspA_tTool_Select
    $mwTools add radio -label "Edit Voxels"        \
      -command { SetTool $DspA_tTool_Edit }  \
      -variable gTool                        \
      -font $kNormalFont                     \
      -value $DspA_tTool_Edit
    $mwTools add radio -label "Select Ctrl Pts"      \
      -command { SetTool $DspA_tTool_CtrlPts } \
      -variable gTool                          \
      -font $kNormalFont                       \
      -value $DspA_tTool_CtrlPts
    $mwTools add radio -label "Edit Voxels w/ Custom Brush" \
      -command { SetTool $DspA_tTool_CustomEdit }     \
      -variable gTool                                 \
      -font $kNormalFont                              \
      -value $DspA_tTool_CustomEdit
    $mwTools add separator
    $mwTools add command -label "Brush Size..." \
      -command { DoBrushInfoDlog }       \
      -font $kNormalFont
    $mwTools add command -label "Custom Brush..." \
      -command { DoCustomBrushDlog }       \
      -font $kNormalFont
    $mwTools add separator
    $mwTools add command -label "Threshold..." \
      -command { DoThresholdDlog }       \
      -font $kNormalFont

    # ctrl pts menu button
    menubutton $mbwCtrlPts -text "Control Points" \
      -menu $mwCtrlPts                      \
      -underline 0                          \
      -font $kNormalFont
    
    # ctrl pts menu items
    menu $mwCtrlPts
    $mwCtrlPts add command -label "New"                          \
      -command { NewControlPoint $gnVolX $gnVolY $gnVolZ } \
      -font $kNormalFont
    $mwCtrlPts add command -label "Select None"   \
      -command { DeselectAllControlPoints } \
      -underline 8                          \
      -font $kNormalFont
    $mwCtrlPts add command -label "Delete Selected"  \
      -command { DeleteSelectedControlPoints } \
      -underline 0                             \
      -font $kNormalFont
    $mwCtrlPts add command -label "Save"       \
      -command { WriteControlPointFile } \
      -font $kNormalFont
     
    # surface menu button
    menubutton $mbwSurface -text "Surface" \
      -menu $mwSurface               \
      -underline 1                   \
      -font $kNormalFont
    
    # surface menu items
    menu $mwSurface
    $mwSurface add command -label "Show Nearest Main Vertex" \
      -command { ShowNearestMainVertex }               \
      -font $kNormalFont
    $mwSurface add command -label "Show Nearest Original Vertex" \
      -command { ShowNearestOriginalVertex }               \
      -font $kNormalFont
    $mwSurface add command -label "Show Nearest Pial Vertex" \
      -command { ShowNearestCanonicalVertex }               \
      -font $kNormalFont
    $mwSurface add command -label "Find Main Vertex..."             \
      -command { DoFindVertexDlog $tkm_tSurfaceType_Current } \
      -font $kNormalFont
    $mwSurface add command -label "Find Original Vertex..."          \
      -command { DoFindVertexDlog $tkm_tSurfaceType_Original } \
      -font $kNormalFont
    $mwSurface add command -label "Find Pial Vertex..."          \
      -command { DoFindVertexDlog $tkm_tSurfaceType_Canonical } \
      -font $kNormalFont
     
    pack $mbwVolume $mbwEdit $mbwDisplay $mbwTools \
      $mbwCtrlPts $mbwSurface  \
      -side left
}

proc CreateCursorFrame { ifwCursor } {

    global kLabelFont kNormalFont kHighlightBgColor
    global gfwVolCoords gfwRASCoords gfwTalCoords gfwVolValue
    global gfwAuxVolValue gfwFuncValue
    global gbLinkedCursor gsVolCoords gsRASCoords 
    global gsTalCoords gsVolName gnVolValue  gsFuncCoords
    global gsAuxVolValue gnAuxVolValue gfwFuncCoords gfFuncValue

    frame $ifwCursor

    set fwLabel             $ifwCursor.fwLabel
    set lwLabel             $fwLabel.lwLabel
    
    set fwLinkCheckbox      $ifwCursor.fwLinkCheckbox
    set cwLink              $fwLinkCheckbox.cwLink
    
    set gfwVolCoords        $ifwCursor.fwVolCoords
    set ewVolCoordsLabel    $gfwVolCoords.ewVolCoordsLabel
    set ewVolCoords         $gfwVolCoords.ewVolCoords

    set gfwRASCoords        $ifwCursor.fwRASCoords
    set ewRASCoordsLabel    $gfwRASCoords.ewRASCoordsLabel
    set ewRASCoords         $gfwRASCoords.ewRASCoords

    set gfwTalCoords        $ifwCursor.fwTalCoords
    set ewTalCoordsLabel    $gfwTalCoords.ewTalCoordsLabel
    set ewTalCoords         $gfwTalCoords.ewTalCoords

    set gfwVolValue         $ifwCursor.fwVolValue
    set ewVolValueLabel     $gfwVolValue.ewVolValueLabel
    set ewVolValue          $gfwVolValue.ewVolValue

    set gfwAuxVolValue      $ifwCursor.fwAuxVolValue
    set ewAuxVolValueLabel  $gfwAuxVolValue.ewAuxVolValueLabel
    set ewAuxVolValue       $gfwAuxVolValue.ewAuxVolValue

    set gfwFuncCoords       $ifwCursor.fwFuncCoords

    set gfwFuncValue        $ifwCursor.fwFuncValue
    set ewFuncValueLabel    $gfwFuncValue.ewFuncValueLabel
    set ewFuncValue         $gfwFuncValue.ewFuncValue

    # the label that goes at the top of the frame
    frame $fwLabel
    label $lwLabel -text "Cursor" \
      -font $kLabelFont
    pack $lwLabel -side left \
      -anchor w

    # link checkbox
    frame $fwLinkCheckbox
    checkbutton $cwLink -text "Linked Cursors"    \
      -command { SetLinkedCursorFlag $gbLinkedCursor } \
      -variable gbLinkedCursor          \
      -font $kNormalFont      \
      -highlightbackground $kHighlightBgColor
    pack $cwLink -side left

    # the volume coords
    frame $gfwVolCoords
    label $ewVolCoordsLabel -text "Volume" \
      -font $kNormalFont
    entry $ewVolCoords                 \
      -textvariable gsVolCoords  \
      -width 26                  \
      -font $kNormalFont         \
      -state disabled            \
      -highlightbackground $kHighlightBgColor \
      -relief flat 
    pack $ewVolCoordsLabel $ewVolCoords -side left \
      -anchor w

    # the RAS coords
    frame $gfwRASCoords
    label $ewRASCoordsLabel -text "RAS" \
      -font $kNormalFont
    entry $ewRASCoords                 \
      -textvariable gsRASCoords  \
      -width 26                  \
      -font $kNormalFont         \
      -state disabled            \
      -highlightbackground $kHighlightBgColor \
      -relief flat 
    pack $ewRASCoordsLabel $ewRASCoords -side left \
      -anchor w

    # the Talairach coords
    frame $gfwTalCoords
    label $ewTalCoordsLabel -text "Talairach" \
      -font $kNormalFont
    entry $ewTalCoords                 \
      -textvariable gsTalCoords  \
      -width 26                  \
      -font $kNormalFont         \
      -state disabled            \
      -highlightbackground $kHighlightBgColor \
      -relief flat 
    pack $ewTalCoordsLabel $ewTalCoords -side left \
      -anchor w

     # the volume value
    frame $gfwVolValue
    entry $ewVolValueLabel             \
      -textvariable gsVolName    \
      -width 0                   \
      -font $kNormalFont         \
      -state disabled            \
      -highlightbackground $kHighlightBgColor \
      -relief flat
    entry $ewVolValue                  \
      -textvariable gnVolValue   \
      -width 5                   \
      -font $kNormalFont         \
      -state disabled            \
      -highlightbackground $kHighlightBgColor \
      -relief flat 
    pack $ewVolValueLabel $ewVolValue -side left \
      -anchor w

     # the aux volume value
    frame $gfwAuxVolValue
    entry $ewAuxVolValueLabel          \
      -textvariable gsAuxVolName \
      -width 0                   \
      -font $kNormalFont         \
      -state disabled            \
      -highlightbackground $kHighlightBgColor \
      -relief flat
    entry $ewAuxVolValue                     \
      -textvariable gnAuxVolValue   \
      -width 5                      \
      -font $kNormalFont            \
      -state disabled               \
      -highlightbackground $kHighlightBgColor    \
      -relief flat 
    pack $ewAuxVolValueLabel $ewAuxVolValue -side left \
      -anchor w

    # the Funcairach coords
    tkm_MakeActiveLabel $gfwFuncCoords \
      "Overlay" gsFuncCoords

    # the functional value
    frame $gfwFuncValue
    label $ewFuncValueLabel -text "Overlay value" \
      -font $kNormalFont
    entry $ewFuncValue                 \
      -textvariable gfFuncValue  \
      -width 5                   \
      -font $kNormalFont         \
      -state disabled            \
      -highlightbackground $kHighlightBgColor \
      -relief flat 
    pack $ewFuncValueLabel $ewFuncValue -side left \
      -anchor w

    # pack the subframes in a column. don't pack vol coords, aux value,
    # or func value. these can be explicity shown later.
    pack $fwLabel $fwLinkCheckbox $gfwRASCoords \
      $gfwTalCoords $gfwVolValue          \
      -side top                           \
      -anchor w

}

proc CreateDisplayFrame { ifwDisplay } {

    global kLabelFont kNormalFont kHighlightBgColor
    global DspA_tOrientation_Sagittal DspA_tOrientation_Horizontal 
    global DspA_tOrientation_Coronal
    global gOrientation gnVolX gnVolY gnVolZ gnZoomLevel
    global gb3D gbCursor gbSelection gbControlPoints gbFunctional
    global gbMainSurface gbOriginalSurface gbCanonicalSurface
    global gbDisplaySurfaceVertices gbInterpolateSurfaceVertices

    frame $ifwDisplay

    set fwLabel            $ifwDisplay.fwLabel
    set lwLabel            $fwLabel.lwLabel

    set fwSagittalButton   $ifwDisplay.fwSagittalButton
    set rwSagittal         $fwSagittalButton.rwSagittal

    set fwSagittalScale    $ifwDisplay.fwSagittalScale
    set swSagittal         $fwSagittalScale.swSagittal
    set ewSagittal         $fwSagittalScale.ewSagittal
    
    set fwHorizontalButton $ifwDisplay.fwHorizontalButton
    set rwHorizontal       $fwHorizontalButton.rwHorizontal

    set fwHorizontalScale  $ifwDisplay.fwHorizontalScale
    set swHorizontal       $fwHorizontalScale.swHorizontal
    set ewHorizontal       $fwHorizontalScale.ewHorizontal
    
    set fwCoronalButton    $ifwDisplay.fwCoronalButton
    set rwCoronal          $fwCoronalButton.rwCoronal

    set fwCoronalScale     $ifwDisplay.fwCoronalScale
    set swCoronal          $fwCoronalScale.swCoronal
    set ewCoronal          $fwCoronalScale.ewCoronal
    
    set fwZoomLabel        $ifwDisplay.fwZoomLabel
    set lwZoom             $fwZoomLabel.lwZoom
    
    set fwZoomScale        $ifwDisplay.fwZoomScale
    set swZoom             $fwZoomScale.swZoom
    set ewZoom             $fwZoomScale.ewZoom

    # the label that goes at the top of the frame
    frame $fwLabel
    label $lwLabel -text "Display" \
      -font $kLabelFont
    pack $lwLabel -side left \
      -anchor w

    # sagittal orientation radio button
    frame $fwSagittalButton
    radiobutton $rwSagittal -text "Sagittal" \
      -command { SetOrientation $gOrientation } \
      -font $kNormalFont \
      -highlightbackground $kHighlightBgColor \
      -variable gOrientation           \
      -relief flat                     \
      -value $DspA_tOrientation_Sagittal
    pack $rwSagittal -side left \
      -anchor w

    # sagittal slice slider
    frame $fwSagittalScale
    scale $swSagittal -orient horizontal \
      -variable gnVolX             \
      -from 0                      \
      -to 255                      \
      -length 200                  \
      -showvalue false
    bind $swSagittal <ButtonRelease> { SetCursor $gnVolX $gnVolY $gnVolZ }
    bind $swSagittal <B1-Motion>     { SetCursor $gnVolX $gnVolY $gnVolZ }
    entry $ewSagittal -textvariable gnVolX \
      -width 3                       \
      -selectbackground green        \
      -insertbackground black
    bind $ewSagittal <Return> { SetCursor $gnVolX $gnVolY $gnVolZ }
    pack $swSagittal    \
      -side left  \
      -anchor w   \
      -expand yes \
      -fill x
    pack $ewSagittal   \
      -side left \
      -anchor w  \
      -expand no

    # horizontal orientaiton radio button
    frame $fwHorizontalButton
    radiobutton $rwHorizontal -text "Horizontal" \
      -command { SetOrientation $gOrientation } \
      -font $kNormalFont   \
      -highlightbackground $kHighlightBgColor \
      -variable gOrientation               \
      -relief flat                         \
      -value $DspA_tOrientation_Horizontal
    pack $rwHorizontal -side left \
      -anchor w

    # horiontal slice slider
    frame $fwHorizontalScale
    scale $swHorizontal -orient horizontal \
      -variable gnVolY               \
      -from 0                        \
      -to 255                        \
      -showvalue false
    bind $swHorizontal <ButtonRelease> { SetCursor $gnVolX $gnVolY $gnVolZ }
    bind $swHorizontal <B1-Motion>     { SetCursor $gnVolX $gnVolY $gnVolZ }
    entry $ewHorizontal -textvariable gnVolY  \
      -selectbackground green           \
      -insertbackground black           \
      -width 3
    bind $ewHorizontal <Return> { SetCursor $gnVolX $gnVolY $gnVolZ }
    pack $swHorizontal    \
      -side left  \
      -anchor w   \
      -expand yes \
      -fill x
    pack $ewHorizontal   \
      -side left \
      -anchor w  \
      -expand no

    # coronal orientation radio button
    frame $fwCoronalButton
    radiobutton $rwCoronal -text "Coronal" \
      -command { SetOrientation $gOrientation } \
      -variable gOrientation         \
      -font $kNormalFont   \
      -highlightbackground $kHighlightBgColor \
      -relief flat                   \
      -value $DspA_tOrientation_Coronal
    pack $rwCoronal -side left \
      -anchor w

    # coronal slice slider
    frame $fwCoronalScale
    scale $swCoronal -orient horizontal \
      -variable gnVolZ            \
      -from 0                     \
      -to 255                     \
      -showvalue false
    bind $swCoronal <ButtonRelease> { SetCursor $gnVolX $gnVolY $gnVolZ }
    bind $swCoronal <B1-Motion>     { SetCursor $gnVolX $gnVolY $gnVolZ }
    entry $ewCoronal -textvariable gnVolZ \
      -selectbackground green       \
      -insertbackground black       \
      -width 3
    bind $ewCoronal <Return> { SetCursor $gnVolX $gnVolY $gnVolZ }
    pack $swCoronal     \
      -side left  \
      -anchor w   \
      -expand yes \
      -fill x
    pack $ewCoronal    \
      -side left \
      -anchor w  \
      -expand no

    # zoom label
    frame $fwZoomLabel
    label $lwZoom -text "Zoom" \
      -font $kNormalFont
    pack $lwZoom -side left \
      -anchor w

    # zoom slider
    frame $fwZoomScale
    scale $swZoom -orient horizontal \
      -command { SetZoomLevel $gnZoomLevel } \
      -variable gnZoomLevel    \
      -from 1                  \
      -to 16                   \
      -showvalue false
    bind $swZoom <ButtonRelease> { SetZoomLevel $gnZoomLevel }
    bind $swZoom <B1-Motion>     { SetZoomLevel $gnZoomLevel }
    entry $ewZoom -textvariable gnZoomLevel \
      -selectbackground green         \
      -insertbackground black         \
      -width 3
    bind $ewZoom <Return> { SetZoomLevel $gnZoomLevel }
    pack $swZoom    \
      -side left  \
      -anchor w   \
      -expand yes \
      -fill x
    pack $ewZoom   \
      -side left \
      -anchor w  \
      -expand no

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

    global kLabelFont kNormalFont

    frame $ifwToolBar

    set fwButtons   $ifwToolBar.fwButtons
    set bwSend      $fwButtons.bwSend
    set bwRead      $fwButtons.bwRead

    frame $fwButtons
    button $bwSend -text "Save Point" \
      -font $kNormalFont          \
      -command { SendCursor }
    button $bwRead -text "Goto Point" \
      -font $kNormalFont         \
      -command { ReadCursor }

    pack $bwSend $bwRead \
      -side left

    pack $fwButtons   \
      -side top \
      -fill y   \
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

    .w.fwMenuBar.mbwVolume.mwVolume delete 15 17
    .w.fwMenuBar.mbwEdit.mwEdit delete 7
    .w.fwMenuBar.mbwTools.mwTools delete 3
    .w.fwMenuBar.mbwDisplay.mwDisplay delete 13
    .w.fwMenuBar.mbwDisplay.mwDisplay delete 17 18
    pack forget .w.fwMenuBar.mbwCtrlPts    
    pack forget .w.fwMenuBar.mbwSurface

}







proc ErrorDlog { isMsg } {

    global gwwTop

    tk_messageBox -type ok \
      -icon error \
      -message $isMsg \
      -title "Error" \
      -parent $gwwTop
}