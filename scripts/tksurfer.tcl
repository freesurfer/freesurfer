#! /usr/bin/tixwish

package require BLT;

source $env(MRI_DIR)/lib/tcl/tkm_common.tcl

foreach sSourceFileName { tkm_wrappers.tcl } {

    set lPath [list "." "$env(MRI_DIR)/lib/tcl"]
    set bFound 0

    foreach sPath $lPath {
  
  if { $bFound == 0 } {
      set sFullFileName [ file join $sPath $sSourceFileName ]
      set nErr [catch { source $sFullFileName } sResult]
      if { $nErr == 0 } {
    dputs "Reading $sFullFileName"
    set bFound 1;
      }
  }
    }

    if { $bFound == 0 } {
  dputs "Couldn't load $sSourceFileName: Not found in $lPath"
    }
}

# ================================================================== CONSTANTS

set ksWindowName "TkSurfer Tools"
set ksImageDir   "$env(MRI_DIR)/lib/images/"

# =========================================================== LINKED VARIABLES

array set gaFileLocations {
    surf { insurf outsurf curv sulc patch annot \
      area origcoords ellcoords vrmlsurf targsurf targcurv }
    fs { fs fm }
    scripts { script }
    label { label }
    bem { dip dec }
    rgb { rgb }
    . { val named_rgbdir num_rgbdir }
}

array set gaFileNameSubDirs {
    kFileName_Surface   "surf"
    kFileName_FieldSign "fs"
    kFileName_Script    "scripts"
    kFileName_Label     "label"
    kFileName_BEM       "bem"
    kFileName_RGB       "rgb"
    kFileName_Home      "."
}

### current transforms
set xrot 0
set yrot 0
set zrot 0
set xtrans 0
set ytrans 0
set scalepercent 100 

### lights: can't be written back (floats in gl struct)
set light0 0.4
set light1 0.0
set light2 0.6
set light3 0.2
set offset 0.25

### events
set userok 0
set blinkflag FALSE
# set initdelay $blinkdelay

### misc defaults
set surfcolor 1
#set colscale 0
#set fthresh 0.1
set shrinksteps 5
set smoothsteps 5
set smoothtype val

# default transform values
set kanDefaultTransform(rotate) 90
set kanDefaultTransform(translate) 10
set kanDefaultTransform(scale) 110

#set home $env(HOME)
#set session ./
#set subject anders

# ====================================================== LINKED VAR MANAGEMENT

# holds local values of all the linked variables
set gaLinkedVars(light0) 0
set gaLinkedVars(light1) 0
set gaLinkedVars(light2) 0 
set gaLinkedVars(light3) 0
set gaLinkedVars(offset) 0
set gaLinkedVars(colscale) 0
set gaLinkedVars(truncphaseflag) 0
set gaLinkedVars(invphaseflag) 0
set gaLinkedVars(revphaseflag) 0
set gaLinkedVars(complexvalflag) 0
set gaLinkedVars(currentvaluefield) 0
set gaLinkedVars(fthresh) 0
set gaLinkedVars(fmid) 0
set gaLinkedVars(fslope) 0
set gaLinkedVars(fslope) 0
set gaLinkedVars(fnumconditions) 0
set gaLinkedVars(fnumtimepoints) 0
set gaLinkedVars(ftimepoint) 0
set gaLinkedVars(fcondition) 0
set gaLinkedVars(fmin) 0
set gaLinkedVars(fmax) 0
set gaLinkedVars(cslope) 0
set gaLinkedVars(cmid) 0
set gaLinkedVars(angle_cycles) 0
set gaLinkedVars(angle_offset) 0
set gaLinkedVars(sulcflag) 0
set gaLinkedVars(surfcolor) 0
set gaLinkedVars(vertexset) 0
set gaLinkedVars(overlayflag) 0
set gaLinkedVars(funcmin) 0
set gaLinkedVars(funcmax) 0
set gaLinkedVars(scalebarflag) 0
set gaLinkedVars(colscalebarflag) 0
set gaLinkedVars(verticesflag) 0
set gaLinkedVars(cmid) 0
set gaLinkedVars(dipavg) 0
set gaLinkedVars(mouseoverflag) 0
set gaLinkedVars(redrawlockflag) 0
set gaLinkedVars(timeresolution) 0
set gaLinkedVars(numprestimpoints) 0
set gaLinkedVars(colortablename) ""

    
# groups of variables that get sent to c code together
array set gaLinkedVarGroups {
    scene { light0 light1 light2 light3 offset }
    overlay { colscale truncphaseflag invphaseflag revphaseflag \
      complexvalflag foffset fthresh fmid fslope fmin fmax \
      fnumtimepoints fnumconditions ftimepoint fcondition}
    curvature { cslope cmid }
    phase { angle_offset angle_cycles }
    inflate { sulcflag }
    view { surfcolor vertexset overlayflag scalebarflag colscalebarflag \
      verticesflag currentvaluefield }
    cvavg { cmid dipavg }
    mouseover { mouseoverflag }
    all { light0 light1 light2 light3 offset colscale truncphaseflag invphaseflag revphaseflag complexvalflag fthresh foffset fmid fslope cslope cmid angle_offset angle_cycles sulcflag surfcolor vertexset overlayflag scalebarflag colscalebarflag verticesflag cmid dipavg mouseoverflag colortablename }
    redrawlock { redrawlockflag }
    graph { timeresolution numprestimpoints }
    label { colortablename }
}

proc SendLinkedVarGroup { iGroup } {
    global gaLinkedVarGroups gaLinkedVars
    set lVars $gaLinkedVarGroups($iGroup)
    foreach var $lVars {
  catch {
    upvar #0 $var varToUpdate;
    set varToUpdate $gaLinkedVars($var)
      }
    }
}

proc UpdateLinkedVarGroup { iGroup } {
    global gaLinkedVarGroups gaLinkedVars
    set lVars $gaLinkedVarGroups($iGroup)
    foreach var $lVars {
  catch {
    upvar #0 $var varToUpdate
    set gaLinkedVars($var) $varToUpdate
      }
    }
}

proc PrintLinkedVarGroup { iGroup } {
    global gaLinkedVarGroups gaLinkedVars
    set lVars $gaLinkedVarGroups($iGroup)
    foreach var $lVars {
  puts "$var=$gaLinkedVars($var)"
    }
}

proc PrintVar { iVar nIndex op } {
    puts "$iVar\($nIndex\) $op"
}

proc UpdateLabel { inSet inLabelIndex isLabel } {
    global glLabel gsaLabelContents
    if { $inSet == 0 } {
  set labelSet cursor
    } else {
  set labelSet mouseover
    }
    set gsaLabelContents([lindex $glLabel $inLabelIndex],value,$labelSet) $isLabel
}

proc UpdateUndoItemLabel { isLabel } {

    .w.fwMenuBar.mbwEdit.mw entryconfigure 1 -label $isLabel
}

proc UpdateValueLabelName { inValueIndex isName } {
    
    global gaScalarValueID gsaLabelContents gaSwapFieldInfo
    if { [info exists gaScalarValueID($inValueIndex,label)] == 0 } {
        puts "UpdateValueLabelName: $inValueIndex invalid"
        return
    }

    set label $gaScalarValueID($inValueIndex,label)

    # set the label contents name.
    set gsaLabelContents($label,name) $isName

    # set the swap field info name.
    set gaSwapFieldInfo($label,label) $isName

    # view->information menu
    .w.fwMenuBar.mbwView.mw.cmw2 entryconfigure [expr 8 + $inValueIndex] \
      -label $isName
    # view->overlay menu
    .w.fwMenuBar.mbwView.mw.cmw7 entryconfigure [expr 1 + $inValueIndex] \
      -label $isName
}

proc SwapValueLabelNames { inValueIndexA inValueIndexB } {

    global gaScalarValueID gsaLabelContents gaSwapFieldInfo
    if { [info exists gaScalarValueID($inValueIndexA,label)] == 0 } {
        puts "UpdateValueLabelName: $inValueIndex invalid"
        return
    }
    if { [info exists gaScalarValueID($inValueIndexB,label)] == 0 } {
        puts "UpdateValueLabelName: $inValueIndex invalid"
        return
    }

    set labelA $gaScalarValueID($inValueIndexA,label)
    set labelB $gaScalarValueID($inValueIndexB,label)

    # get old values
    set sLabelContentsA $gsaLabelContents($labelA,name)
    set sLabelContentsB $gsaLabelContents($labelB,name)

    # set the label contents name.
    set gsaLabelContents($labelA,name) $sLabelContentsB
    set gsaLabelContents($labelB,name) $sLabelContentsA

    # set the swap field info name.
    set gaSwapFieldInfo($labelA,label) $sLabelContentsB
    set gaSwapFieldInfo($labelB,label) $sLabelContentsA
}

# ==================================================================== GLOBALS

set gNextTransform(rotate,x) 0
set gNextTransform(rotate,y) 0
set gNextTransform(rotate,z) 0
set gNextTransform(rotate,deg) 0
set gNextTransform(translate,x) 0
set gNextTransform(translate,y) 0
set gNextTransform(translate,z) 0
set gNextTransform(translate,dist) 0
set gNextTransform(scale) 0
set gNextTransform(scale,amt) 100

# labels
set glLabel { \
  kLabel_VertexIndex \
  kLabel_Distance \
  kLabel_Coords_RAS \
  kLabel_Coords_MniTal \
  kLabel_Coords_Tal \
  kLabel_Coords_Index \
  kLabel_Coords_Normal \
  kLabel_Coords_Sphere_XYZ \
  kLabel_Coords_Sphere_RT \
  kLabel_Curvature \
  kLabel_Fieldsign \
  kLabel_Val \
  kLabel_Val2 \
  kLabel_ValBak \
  kLabel_Val2Bak \
  kLabel_ValStat \
  kLabel_Amplitude \
  kLabel_Angle \
  kLabel_Degree \
  kLabel_Annotation \
  kLabel_MRIValue \
  kLabel_Parcellation_Name }
foreach label $glLabel {
    set gfwaLabel($label,cursor) ""
    set gfwaLabel($label,mouseover) ""
}

set gsaLabelContents(kLabel_VertexIndex,name)       "Vertex Index"
set gsaLabelContents(kLabel_Distance,name)          "Distance"
set gsaLabelContents(kLabel_Coords_RAS,name)        "Vertex RAS"
set gsaLabelContents(kLabel_Coords_MniTal,name)     "Vertex MNI Talairach"
set gsaLabelContents(kLabel_Coords_Tal,name)        "Vertex Talairach"
set gsaLabelContents(kLabel_Coords_Index,name)      "MRI Index"
set gsaLabelContents(kLabel_Coords_Normal,name)     "Vertex Normal"
set gsaLabelContents(kLabel_Coords_Sphere_XYZ,name) "Spherical X, Y, Z"
set gsaLabelContents(kLabel_Coords_Sphere_RT,name)  "Spherical Rho, Theta"
set gsaLabelContents(kLabel_Curvature,name)         "Curvature"
set gsaLabelContents(kLabel_Fieldsign,name)         "Field Sign"
set gsaLabelContents(kLabel_Val,name)               "Overlay Layer 1"
set gsaLabelContents(kLabel_Val2,name)              "Overlay Layer 2"
set gsaLabelContents(kLabel_ValBak,name)            "Overlay Layer 3"
set gsaLabelContents(kLabel_Val2Bak,name)           "Overlay Layer 4"
set gsaLabelContents(kLabel_ValStat,name)           "Overlay Layer 5"
set gsaLabelContents(kLabel_Amplitude,name)         "Amplitude"
set gsaLabelContents(kLabel_Angle,name)             "Angle"
set gsaLabelContents(kLabel_Degree,name)            "Degree"
set gsaLabelContents(kLabel_Annotation,name)        "Annotation"
set gsaLabelContents(kLabel_MRIValue,name)          "MRI Value"
set gsaLabelContents(kLabel_Parcellation_Name,name) "Parcellation: "

foreach label $glLabel {
    set gsaLabelContents($label,value,cursor)    "none"
    set gsaLabelContents($label,value,mouseover) "none"
}

set glSwapField { \
  kField_Curv \
  kField_CurvBak \
  kField_Val \
  kField_Val2 \
  kField_ValBak \
  kField_Val2Bak \
  kField_Stat \
  kField_ImagVal }
set gaSwapFieldInfo(kField_Curv,label)     "curv"
set gaSwapFieldInfo(kField_CurvBak,label)  "curvbak"
set gaSwapFieldInfo(kField_Val,label)      "val"
set gaSwapFieldInfo(kField_Val2,label)     "val2"
set gaSwapFieldInfo(kField_ValBak,label)   "valbak"
set gaSwapFieldInfo(kField_Val2Bak,label)  "val2bak"
set gaSwapFieldInfo(kField_Stat,label)     "stat"
set gaSwapFieldInfo(kField_ImagVal,label)  "imag_val"
set n 0
foreach field $glSwapField {
    set gaSwapFieldInfo($field,index) $n
    incr n
}

# scalar values -> labels
set gaScalarValueID(kLabel_Val,index) 0
set gaScalarValueID(0,label) kLabel_Val
set gaScalarValueID(kLabel_Val2,index) 1
set gaScalarValueID(1,label) kLabel_Val2
set gaScalarValueID(kLabel_ValBak,index) 2
set gaScalarValueID(2,label) kLabel_ValBak
set gaScalarValueID(kLabel_Val2Bak,index) 3
set gaScalarValueID(3,label) kLabel_Val2Bak
set gaScalarValueID(kLabel_ValStat,index) 4
set gaScalarValueID(4,label) kLabel_ValStat

# tool bar frames
set gfwaToolBar(main)  ""

# ========================================================= BUILDING INTERFACE


proc DoFileDlog { which } {
    global tDlogSpecs
    tkm_DoFileDlog $tDlogSpecs($which)
}

proc DoSwapSurfaceFieldsDlog {} {

    global gDialog
    global glSwapField gaSwapFieldInfo

    set wwDialog .wwSwapSurfaceFieldsDlog

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Swap Surface Fields" {-borderwidth 10}] } {

  set fwMain             $wwDialog.fwMain
  set fwSwapField1       $fwMain.fwSwapField1
  set fwSwapField2       $fwMain.fwSwapField2
  set fwButtons          $wwDialog.fwButtons

  frame $fwMain

  tixOptionMenu $fwSwapField1 -label "Swap" \
    -variable swapField1 \
    -options {
      label.anchor e
      label.width 5
      menubutton.width 8
  }
  
  tixOptionMenu $fwSwapField2 -label "with" \
    -variable swapField2 \
    -options {
      label.anchor e
      label.width 5
      menubutton.width 8
  }

  foreach field $glSwapField {
      $fwSwapField1 add command $gaSwapFieldInfo($field,index) \
        -label $gaSwapFieldInfo($field,label)
      $fwSwapField2 add command $gaSwapFieldInfo($field,index) \
        -label $gaSwapFieldInfo($field,label)
  }
  
  # buttons.
  tkm_MakeApplyCloseButtons $fwButtons $wwDialog \
    { swap_vertex_fields $swapField1 $swapField2 } {}

  pack $fwSwapField1 $fwSwapField2 \
    -side left

  pack $fwMain $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoConfigLightingDlog {} {

    global gDialog
    global gaLinkedVars

    set wwDialog .wwConfigLightingDlog

    UpdateLinkedVarGroup scene

    if { [Dialog_Create $wwDialog "Configure Lighting" {-borderwidth 10}] } {

  set fwMain             $wwDialog.fwMain
  set fwLights           $wwDialog.fwLights
  set fwBrightness       $wwDialog.fwBrightness
  set fwButtons          $wwDialog.fwButtons

  frame $fwMain

  # sliders for lights 0 - 3
  tkm_MakeSliders $fwLights { \
    { {"Light 0" ""} gaLinkedVars(light0) 0 1 200 {} 1 0.1} \
    { {"Light 1" ""} gaLinkedVars(light1) 0 1 200 {} 1 0.1} \
    { {"Light 2" ""} gaLinkedVars(light2) 0 1 200 {} 1 0.1} \
    { {"Light 3" ""} gaLinkedVars(light3) 0 1 200 {} 1 0.1} }

  # slider for brightness offset
  tkm_MakeSliders $fwBrightness { \
    { {"Brightness" ""} gaLinkedVars(offset) 0 1 100 {} 1 0.1} }

  # buttons.
  tkm_MakeApplyCloseButtons $fwButtons $wwDialog \
    { SendLinkedVarGroup scene; \
    do_lighting_model $gaLinkedVars(light0) $gaLinkedVars(light1) $gaLinkedVars(light2) $gaLinkedVars(light3) $offset; \
    UpdateAndRedraw } {}

  pack $fwMain $fwLights $fwBrightness $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoConfigOverlayDisplayDlog {} {

    global gDialog
    global gaLinkedVars

    set wwDialog .wwConfigOverlayDisplayDlog

    UpdateLinkedVarGroup overlay

    if { [Dialog_Create $wwDialog "Configure Overlay Display" {-borderwidth 10}] } {

  set fwMain             $wwDialog.fwMain
  set fwTimePoint        $wwDialog.fwTimePoint
  set fwCondition        $wwDialog.fwCondition
  set fwColorScale       $wwDialog.fwColorScale
  set fwFlags            $wwDialog.fwFlags
  set lfwThreshold       $wwDialog.lfwThreshold
  set fwButtons          $wwDialog.fwButtons

  frame $fwMain

  set nMaxCondition [expr $gaLinkedVars(fnumconditions) - 1]
  if { $nMaxCondition < 0 } {
      set nMaxCondition 0
  }
  set nMaxTimePoint [expr $gaLinkedVars(fnumtimepoints) - 1]
  if { $nMaxTimePoint < 0 } {
      set nMaxTimePoint 0
  }
  
  tkm_MakeEntryWithIncDecButtons \
    $fwTimePoint "Time Point (0-$nMaxTimePoint)" \
    gaLinkedVars(ftimepoint) \
    {} 1
  
  tkm_MakeEntryWithIncDecButtons \
    $fwCondition "Condition (0-$nMaxCondition)" \
    gaLinkedVars(fcondition) \
    {} 1

  # color scale
  tkm_MakeRadioButtons $fwColorScale y "Color Scale" \
    gaLinkedVars(colscale) { \
    { text "Color Wheel (Complex)" 0 {} } \
    { text "RYGB Wheel (Complex)" 8 {} } \
    { text "Two Condition Green Red (Complex)" 4 {} } \
    { text "Green to Red (Signed)" 7 {} } \
    { text "Heat Scale (Stat, Positive)" 1 {} } \
    { text "Blue to Red (Signed)" 6 {} } \
    { text "Not Here (Signed)" 9 {} } }

  # checkboxes for truncate, inverse, reverse phase, complex values
  tkm_MakeCheckboxes $fwFlags y { \
    { text "Truncate" gaLinkedVars(truncphaseflag) {} } \
    { text "Inverse" gaLinkedVars(invphaseflag) {} } \
    { text "Reverse" gaLinkedVars(revphaseflag) {} } \
    { text "Complex" gaLinkedVars(complexvalflag) {} } }

  # sliders for thresh, mid, field for slope
  tixLabelFrame $lfwThreshold \
    -label "Threshold" \
    -labelside acrosstop \
    -options { label.padX 5 }

  set fwThresholdSub     [$lfwThreshold subwidget frame]
  set fwThresholdSliders $fwThresholdSub.fwThresholdSliders
  set fwThresholdSlope   $fwThresholdSub.fwThresholdSlope

  tkm_MakeSliders $fwThresholdSliders [list \
    [list {"Threshold offset"} gaLinkedVars(foffset) \
    -10000 10000 100 {} 1 0.25] \
    [list {"Threshold minimum"} gaLinkedVars(fthresh) \
    $gaLinkedVars(fmin) $gaLinkedVars(fmax) 100 {} 1 0.25] \
    [list {"Threshold midpoint"} gaLinkedVars(fmid) \
    $gaLinkedVars(fmin) $gaLinkedVars(fmax) 100 {} 1 0.25]]
  tkm_MakeEntry $fwThresholdSlope "Threshold slope" \
    gaLinkedVars(fslope) 6

  pack $fwThresholdSliders $fwThresholdSlope \
    -side top \
    -anchor w

  # buttons.
  tkm_MakeApplyCloseButtons $fwButtons $wwDialog \
    { SendLinkedVarGroup overlay;  \
      sclv_set_current_timepoint $gaLinkedVars(ftimepoint) $gaLinkedVars(fcondition); \
      UpdateAndRedraw; } {}

  pack $fwMain $fwTimePoint $fwCondition $fwColorScale \
    $fwFlags $lfwThreshold $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc UpdateOverlayDlogInfo {} {

    global gaLinkedVars
    
    # change the time point and condition range labels
    catch { \
  set nMaxCondition [expr $gaLinkedVars(fnumconditions) - 1]
  if { $nMaxCondition < 0 } {
      set nMaxCondition 0
  }
  set nMaxTimePoint [expr $gaLinkedVars(fnumtimepoints) - 1]
  if { $nMaxTimePoint < 0 } {
      set nMaxTimePoint 0
  }
  .wwConfigOverlayDisplayDlog.fwTimePoint.control config \
    -label "Time Point (0-$nMaxTimePoint)"
  .wwConfigOverlayDisplayDlog.fwCondition.control config \
    -label "Condition (0-$nMaxCondition)"
    }
  
    # change the range of the threshold sliders.
    catch {
  [.wwConfigOverlayDisplayDlog.lfwThreshold subwidget frame].fwThresholdSliders.sw0 config \
    -from $gaLinkedVars(fmin)
  [.wwConfigOverlayDisplayDlog.lfwThreshold subwidget frame].fwThresholdSliders.sw0 config \
    -to $gaLinkedVars(fmax)
    }

    catch {
  [.wwConfigOverlayDisplayDlog.lfwThreshold subwidget frame].fwThresholdSliders.sw1 config \
    -from $gaLinkedVars(fmin)
  [.wwConfigOverlayDisplayDlog.lfwThreshold subwidget frame].fwThresholdSliders.sw1 config \
    -to $gaLinkedVars(fmax)
    }
}

proc DoConfigCurvatureDisplayDlog {} {

    global gDialog
    global gaLinkedVars
    UpdateLinkedVarGroup curvature

    set wwDialog .wwConfigCurvatureDisplayDlog

    if { [Dialog_Create $wwDialog "Configure Curvature Display" {-borderwidth 10}] } {

  set fwMain             $wwDialog.fwMain
  set lfwThreshold       $wwDialog.lfwThreshold
  set fwButtons          $wwDialog.fwButtons

  frame $fwMain

  # fields for slope and midpoint
  tixLabelFrame $lfwThreshold \
    -label "Threshold" \
    -labelside acrosstop \
    -options { label.padX 5 }

  set fwThresholdSub     [$lfwThreshold subwidget frame]
  set fwThresholdSliders $fwThresholdSub.fwThresholdSliders
  set fwThresholdSlope   $fwThresholdSub.fwThresholdSlope

  tkm_MakeSliders $fwThresholdSliders [list \
    [list {"Threshold midpoint"} gaLinkedVars(cmid) \
    -1 1 100 {} 1 0.05 ]]
  tkm_MakeEntry $fwThresholdSlope "Threshold slope" \
    gaLinkedVars(cslope) 6

  pack $fwThresholdSliders $fwThresholdSlope \
    -side top \
    -anchor w

  # buttons.
  tkm_MakeApplyCloseButtons $fwButtons $wwDialog \
    { SendLinkedVarGroup curvature; UpdateAndRedraw } {}

  pack $fwMain $lfwThreshold $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoConfigPhaseEncodedDataDisplayDlog {} {

    global gDialog
    global gaLinkedVars
    UpdateLinkedVarGroup phase

    set wwDialog .wwConfigPhaseEncodedDataDisplayDlog

    if { [Dialog_Create $wwDialog "Configure Phase Encoded Data Display" {-borderwidth 10}] } {

  set fwMain             $wwDialog.fwMain
  set fwCycles           $wwDialog.fwCycles
  set fwOffset           $wwDialog.fwOffset
  set fwButtons          $wwDialog.fwButtons

  frame $fwMain
  
  # fields for angle cycles and angle offset
  tkm_MakeEntry $fwCycles "Angle Cycles: " gaLinkedVars(angle_cycles) 6 
  tkm_MakeEntry $fwOffset "Angle Offset: " gaLinkedVars(angle_offset) 6 

  # buttons.
  tkm_MakeApplyCloseButtons $fwButtons $wwDialog \
    { SendLinkedVarGroup phase; UpdateAndRedraw } {}

  pack $fwMain $fwCycles $fwOffset $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoLoadOverlayDlog {} {

    global gDialog gaLinkedVars
    global gaScalarValueID gsaLabelContents

    set wwDialog .wwLoadOverlayDlog

    set knWidth 400

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Load Overlay" {-borderwidth 10}] } {

  set fwFile             $wwDialog.fwFile
  set fwFileNote         $wwDialog.fwFileNote
  set fwField            $wwDialog.fwField
  set fwFieldNote        $wwDialog.fwFieldNote
  set fwButtons          $wwDialog.fwButtons

  set sFileName ""
  tkm_MakeFileSelector $fwFile "Load Overlay:" sFileName \
    [list ExpandFileName "" kFileName_Surface]

  tkm_MakeSmallLabel $fwFileNote "The file name of the values" 400

  tixOptionMenu $fwField -label "Into Field:" \
    -variable nFieldIndex \
    -options {
      label.anchor e
      label.width 5
      menubutton.width 8
  }
  
  tkm_MakeSmallLabel $fwFieldNote "The layer to load the values into" 400

  set nIndex 0
  while { [info exists gaScalarValueID($nIndex,label)] } {
      $fwField add command $nIndex \
        -label $gsaLabelContents($gaScalarValueID($nIndex,label),name)
      incr nIndex
  }
  
  # buttons.
        tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    {set val [ExpandFileName $sFileName kFileName_Surface]; \
    DoLoadValueFile $nFieldIndex }

  pack $fwFile $fwFileNote $fwField $fwFieldNote $fwButtons \
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
  } $wwDialog $knWidth $wwDialog] 
    }
}

proc DoLoadValueFile { inField } {

    global val

    # if ends in bfloat, pass to DoSpecifyStemAndRegistration
    set sExtension [file extension $val]
    if { $sExtension == ".bfloat" || $sExtension == ".bshort" } {
  DoSpecifyStemAndRegistration $inField

    } else {
  # else pass to normal function
  sclv_read_binary_values $inField
  UpdateAndRedraw
    }
}

proc DoSpecifyStemAndRegistration { inField } {

    global val sPath sStem sRegistration
    global gDialog gaLinkedVars

    set wwDialog .wwDoSpecifyStemAndRegistration

    set knWidth 400
    set sPath [file dirname $val]
    set sStem [lindex [split [file rootname [file tail $val]] _] 0]
    set sRegistration ""

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Specify Registration" {-borderwidth 10}] } {

  set fwStem             $wwDialog.fwStem
  set fwStemNote         $wwDialog.fwStemNote
  set fwReg              $wwDialog.fwReg
  set fwRegNote          $wwDialog.fwRegNote
  set fwButtons          $wwDialog.fwButtons


  tkm_MakeEntry $fwStem "Stem:" sStem

  tkm_MakeSmallLabel $fwStemNote "The stem of the volume" 400

  tkm_MakeFileSelector $fwReg "Registration file:" sRegistration

  tkm_MakeSmallLabel $fwRegNote "The file name of the registration file to load. Leave blank to use register.dat in the same directory." 

  # buttons.
        tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    "sclv_read_bfile_values $inField \$sPath \$sStem \$sRegistration; UpdateAndRedraw"

  pack $fwStem $fwStemNote $fwReg $fwRegNote $fwButtons \
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
  } $wwDialog $knWidth $wwDialog] 
    }
}

proc DoSaveValuesAsDlog {} {

    global gDialog
    global gaScalarValueID gsaLabelContents

    set wwDialog .wwSaveValuesAs

    # try to create the dlog...
    if { [Dialog_Create $wwDialog "Save Values As" {-borderwidth 10}] } {

  set fwFile             $wwDialog.fwFile
  set fwFileNote         $wwDialog.fwFileNote
  set fwField            $wwDialog.fwField
  set fwFieldNote        $wwDialog.fwFieldNote
  set fwButtons          $wwDialog.fwButtons

  set sFileName ""
  tkm_MakeFileSelector $fwFile "Save Values:" sFileName {}

  tkm_MakeSmallLabel $fwFileNote "The file name of the values file to create" 400

  tixOptionMenu $fwField -label "From Field:" \
    -variable nFieldIndex \
    -options {
      label.anchor e
      label.width 5
      menubutton.width 8
  }
  
  tkm_MakeSmallLabel $fwFieldNote "The layer to save" 400

  set nIndex 0
  while { [info exists gaScalarValueID($nIndex,label)] } {
      $fwField add command $nIndex \
        -label $gsaLabelContents($gaScalarValueID($nIndex,label),name)
      incr nIndex
  }
  
  # buttons.
        tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    {set val [ExpandFileName $sFileName kFileName_Surface]; \
    sclv_write_binary_values $nFieldIndex}

  pack $fwFile $fwFileNote $fwField $fwFieldNote $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoSmoothOverlayDlog {} {

    global gaScalarValueID gsaLabelContents
    global gDialog

    set wwDialog .wwSmoothOverlayDlog

    if { [Dialog_Create $wwDialog "Smooth Overlay" {-borderwidth 10}] } {

  set fwMain             $wwDialog.fwMain
  set fwSteps            $wwDialog.fwSteps
  set fwTarget           $wwDialog.fwTarget
  set fwButtons          $wwDialog.fwButtons

  frame $fwMain

  # field for number of steps
  tkm_MakeEntry $fwSteps "Number of Steps: " nSteps 6 

  # target scalar field
  tixOptionMenu $fwTarget -label "Target Field:" \
    -variable nFieldIndex \
    -options {
      label.anchor e
      label.width 5
      menubutton.width 8
  }
  
  set nIndex 0
  while { [info exists gaScalarValueID($nIndex,label)] } {
      $fwTarget add command $nIndex \
        -label $gsaLabelContents($gaScalarValueID($nIndex,label),name)
      incr nIndex
  }

  # buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { DoSmoothOverlay $nSteps $nFieldIndex; UpdateAndRedraw } {}

  pack $fwMain $fwSteps $fwTarget $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoSmoothCurvatureDlog {} {

    global gDialog

    set wwDialog .wwSmoothCurvDlog

    if { [Dialog_Create $wwDialog "Smooth Curvature" {-borderwidth 10}] } {

  set fwMain             $wwDialog.fwMain
  set fwSteps            $wwDialog.fwSteps
  set fwButtons          $wwDialog.fwButtons

  frame $fwMain

  # field for number of steps
  tkm_MakeEntry $fwSteps "Number of Steps: " nSteps 6 

  # buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { DoSmoothCurvature $nSteps; UpdateAndRedraw } {}

  pack $fwMain $fwSteps $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoInflateDlog {} {

    global gDialog
    global gaLinkedVars
    UpdateLinkedVarGroup inflate

    set wwDialog .wwInflateDlog

    if { [Dialog_Create $wwDialog "Inflate" {-borderwidth 10}] } {

  set fwMain             $wwDialog.fwMain
  set fwSteps            $wwDialog.fwSteps
  set fwSulc             $wwDialog.fwSulc
  set fwButtons          $wwDialog.fwButtons

  frame $fwMain

  # field for steps
  tkm_MakeEntry $fwSteps "Number of Steps: " nSteps 6 

  # cb for sulc flag
  tkm_MakeCheckboxes $fwSulc y { \
    {text "Sulc Sum" gaLinkedVars(sulcflag) {} } }

  # buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { SendLinkedVarGroup inflate; \
    DoInflate $nSteps; UpdateAndRedraw } {}

  pack $fwMain $fwSteps $fwSulc $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoDecimationDlog {} {

    global gDialog

    set wwDialog .wwDecimationDlog

    if { [Dialog_Create $wwDialog "Write Decimation" {-borderwidth 10}] } {

  set fwMain             $wwDialog.fwMain
  set fwFileName         $wwDialog.fwFileName
  set fwSpacing          $wwDialog.fwSpacing
  set fwButtons          $wwDialog.fwButtons

  frame $fwMain

  # make file name selector
  set sFileName "hi"
  tkm_MakeFileSelector $fwFileName \
    "Write Decimation File:" sFileName

  # field for spacing
  tkm_MakeEntry $fwSpacing "Spacing: " fSpacing 6 

  # buttons.
  tkm_MakeCancelOKButtons $fwButtons $wwDialog \
    { tkm_UpdateFileSelectorVariable $fwFileName; \
    DoDecimation $sFileName $fSpacing; UpdateAndRedraw }

  pack $fwMain $fwFileName $fwSpacing $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoSendToSubjectDlog {} {

    global gDialog

    set wwDialog .wwSendToSubjectDlog

    if { [Dialog_Create $wwDialog "Send To Subject" {-borderwidth 10}] } {

  set fwMain             $wwDialog.fwMain
  set fwSubject          $wwDialog.fwSubject
  set fwButtons          $wwDialog.fwButtons
  
  frame $fwMain
  
  # field for subject name
  tkm_MakeEntry $fwSubject "Subject: " sSubject 20
  
  # buttons.
  tkm_MakeApplyCloseButtons $fwButtons $wwDialog \
    { send_to_subject $sSubject } {}
  
  pack $fwMain $fwSubject $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc DoCustomFillDlog {} {

    global gDialog

    set wwDialog .wwCustomFillDlog

    if { [Dialog_Create $wwDialog "Custom Fill" {-borderwidth 10}] } {
  
  set fwMain             $wwDialog.fwMain
  set fwText             $wwDialog.fwText
  set fwFlags            $wwDialog.fwFlags
  set fwButtons          $wwDialog.fwButtons
  
  frame $fwMain
  
  # hint text
  tkm_MakeBigLabel $fwText "Fill:"

  # cbs for flags
  tkm_MakeCheckboxes $fwFlags y { \
    {text "Up to and including boundaries" bNoBoundary {} } \
    {text "Up to other labels" bNoLabel {} } \
    {text "Up to and including different curvature" bNoCmid {} } \
    {text "Up to functional values below threshold" bNoFThresh {} } \
      }
  
  # buttons.
  tkm_MakeApplyCloseButtons $fwButtons $wwDialog \
    { fill_flood_from_cursor $bNoBoundary $bNoLabel $bNoCmid $bNoFThresh; UpdateAndRedraw } {} "Fill"
  
  pack $fwMain $fwText $fwFlags $fwButtons \
    -side top       \
    -expand yes     \
    -fill x         \
    -padx 5         \
    -pady 5
    }
}

proc CreateWindow { iwwTop } {
    global ksWindowName
    frame $iwwTop
    wm title . $ksWindowName
    wm withdraw .
}

proc CreateMenuBar { ifwMenuBar } {

    global gaLinkedVars
    global gbShowToolBar gbShowLabel
    UpdateLinkedVarGroup view

    set mbwFile   $ifwMenuBar.mbwFile
    set mbwEdit   $ifwMenuBar.mbwEdit
    set mbwView   $ifwMenuBar.mbwView
    set mbwTools  $ifwMenuBar.mbwTools

    frame $ifwMenuBar -border 2 -relief raised

    # file menu button
    tkm_MakeMenu $mbwFile "File" { \
      {command \
      "Load Surface..." \
      {DoFileDlog LoadSurface} } \
      \
      { cascade "Load Surface Configuration..." { \
      { command \
      "Main Vertices" \
      {DoFileDlog LoadMainSurface} } \
      \
      { command \
      "Inflated Vertices" \
      {DoFileDlog LoadInflatedSurface} } \
      \
      { command \
      "White Vertices" \
      {DoFileDlog LoadWhiteSurface} } \
      \
      { command \
      "Pial Vertices" \
      {DoFileDlog LoadPialSurface} } \
      \
      { command \
      "Original Vertices" \
      {DoFileDlog LoadOriginalSurface} } } } \
      \
      {command \
      "Save Surface" \
      {} } \
      \
      {command \
      "Save Surface As..." \
      {DoFileDlog SaveSurfaceAs} } \
      \
      { separator } \
      \
      \
      {command "Load Overlay..." \
      {DoLoadOverlayDlog}} \
      \
      {command "Save Overlay As..." \
      {DoSaveValuesAsDlog} \
      mg_OverlayLoaded } \
      \
      {command \
      "Load Time Course..." \
      {DoFileDlog LoadTimeCourse} } \
      \
      { separator } \
      \
      {cascade "Curvature" { \
      {command "Load Curvature..." \
      {DoFileDlog LoadCurvature}} \
      \
      {command "Save Curvature" \
      {CheckFileAndDoCmd $curv write_binary_curv} \
      mg_CurvatureLoaded } \
      \
      {command "Save Curvature As..." \
      {DoFileDlog SaveCurvatureAs} \
      mg_CurvatureLoaded } \
    }   } \
      \
      {cascade "Patch" { \
      {command "Load Patch..." \
      {DoFileDlog LoadPatch}} \
      \
      {command "Save Patch" \
      {CheckFileAndDoCmd $patch write_binary_patch} \
      mg_PatchLoaded } \
      \
      {command "Save Patch As..." \
      {DoFileDlog SavePatchAs} \
      mg_PatchLoaded } \
    }   } \
      \
      {cascade "Label" { \
      { command "Load Color Table..." \
      { DoFileDlog LoadColorTable } } \
      \
      {command "Load Label..." \
      {DoFileDlog LoadLabel}} \
      \
      {command "Save Selected Label..." \
      {DoFileDlog SaveLabelAs} \
      mg_LabelLoaded } \
      \
      {command "Import Annotation..." \
      {DoFileDlog ImportAnnotation} } \
      \
      {command "Export Annotation..." \
      {DoFileDlog Export} \
      mg_LabelLoaded } \
      }   } \
      \
      {cascade "Field Sign" { \
      {command "Load Field Sign..." \
      {DoFileDlog LoadFieldSign}} \
      \
      {command "Save Field Sign" \
      {CheckFileAndDoCmd $fs write_fieldsign} \
      mg_FieldSignLoaded } \
      \
      {command "Save Field Sign As..." \
      {DoFileDlog SaveFieldSignAs} \
      mg_FieldSignLoaded } \
    }   } \
      {cascade "Field Mask" { \
      {command "Load Field Mask..." \
      {DoFileDlog LoadFieldMask}} \
      \
      {command "Save Field Mask" \
      {CheckFileAndDoCmd $fm write_fsmask} \
      mg_FieldMaskLoaded } \
      \
      {command "Save Field Mask As..." \
      {DoFileDlog SaveFieldMaskAs} \
      mg_FieldMaskLoaded } \
    }   } 
      { separator } \
      \
          {command \
      "Quit" \
      {exit} } }

    # edit menu 
    tkm_MakeMenu $mbwEdit "Edit" { \
      { command \
      "Nothing to Undo" \
      undo_last_action } \
      \
      { separator } \
      \
      { command \
      "Unmark All Vertices" \
      { clear_all_vertex_marks; UpdateAndRedraw } } \
  }
    
    # view menu
    tkm_MakeMenu $mbwView "View" { \
      { cascade \
      "Tool Bars" { \
      { check \
      "Main" \
      "ShowToolBar cut $gbShowToolBar(main)" \
      gbShowToolBar(main) } }} \
      { cascade \
      "Information" { \
      { check \
      "Vertex Index" \
      "ShowLabel kLabel_VertexIndex $gbShowLabel(kLabel_VertexIndex)"\
      gbShowLabel(kLabel_VertexIndex) } \
      { check \
      "Distance" \
      "ShowLabel kLabel_Distance $gbShowLabel(kLabel_Distance)"\
      gbShowLabel(kLabel_Distance) } \
      { check \
      "Vertex RAS" \
      "ShowLabel kLabel_Coords_RAS $gbShowLabel(kLabel_Coords_RAS)"\
      gbShowLabel(kLabel_Coords_RAS) } \
      { check \
      "Vertex MNI Talairach" \
           "ShowLabel kLabel_Coords_MniTal $gbShowLabel(kLabel_Coords_MniTal)"\
      gbShowLabel(kLabel_Coords_MniTal) } \
      { check \
      "Vertex Talairach" \
      "ShowLabel kLabel_Coords_Tal $gbShowLabel(kLabel_Coords_Tal)"\
      gbShowLabel(kLabel_Coords_Tal) } \
      { check \
      "MRI Index" \
      "ShowLabel kLabel_Coords_Index $gbShowLabel(kLabel_Coords_Index)"\
      gbShowLabel(kLabel_Coords_Index) } \
      { check \
      "Vertex Normal" \
    "ShowLabel kLabel_Coords_Normal $gbShowLabel(kLabel_Coords_Normal)"\
      gbShowLabel(kLabel_Coords_Normal) } \
      { check \
      "Spherical X, Y, Z" \
           "ShowLabel kLabel_Coords_Sphere_XYZ $gbShowLabel(kLabel_Coords_Sphere_XYZ)"\
      gbShowLabel(kLabel_Coords_Sphere_XYZ) } \
      { check \
      "Spherical Rho, Theta" \
           "ShowLabel kLabel_Coords_Sphere_RT $gbShowLabel(kLabel_Coords_Sphere_RT)"\
      gbShowLabel(kLabel_Coords_Sphere_RT) } \
      { check \
      "Curvature" \
      "ShowLabel kLabel_Curvature $gbShowLabel(kLabel_Curvature)"\
      gbShowLabel(kLabel_Curvature) \
      mg_CurvatureLoaded } \
      { check \
      "Field Sign" \
      "ShowLabel kLabel_Fieldsign $gbShowLabel(kLabel_Fieldsign)"\
      gbShowLabel(kLabel_Fieldsign) } \
      { check \
      "Overlay Layer 1" \
      "ShowLabel kLabel_Val $gbShowLabel(kLabel_Val)"\
      gbShowLabel(kLabel_Val) \
      mg_OverlayLoaded } \
      { check \
      "Overlay Layer 2" \
      "ShowLabel kLabel_Val2 $gbShowLabel(kLabel_Val2)"\
      gbShowLabel(kLabel_Val2) \
      mg_OverlayLoaded } \
      { check \
      "Overlay Layer 3" \
      "ShowLabel kLabel_ValBak $gbShowLabel(kLabel_ValBak)"\
      gbShowLabel(kLabel_ValBak) \
      mg_OverlayLoaded } \
      { check \
      "Overlay Layer 4" \
      "ShowLabel kLabel_Val2Bak $gbShowLabel(kLabel_Val2Bak)"\
      gbShowLabel(kLabel_Val2Bak) \
      mg_OverlayLoaded } \
      { check \
      "Overlay Layer 5" \
      "ShowLabel kLabel_ValStat $gbShowLabel(kLabel_ValStat)"\
      gbShowLabel(kLabel_ValStat) \
      mg_OverlayLoaded } \
      { check \
      "Amplitude" \
      "ShowLabel kLabel_Amplitude $gbShowLabel(kLabel_Amplitude)"\
      gbShowLabel(kLabel_Amplitude) } \
      { check \
      "Angle" \
      "ShowLabel kLabel_Angle $gbShowLabel(kLabel_Angle)"\
      gbShowLabel(kLabel_Angle) } \
      { check \
      "Degree" \
      "ShowLabel kLabel_Degree $gbShowLabel(kLabel_Degree)"\
      gbShowLabel(kLabel_Degree) } \
      { check \
      "Annotation" \
      "ShowLabel kLabel_Annotation $gbShowLabel(kLabel_Annotation)"\
      gbShowLabel(kLabel_Annotation) } \
      { check \
      "MRI Value" \
      "ShowLabel kLabel_MRIValue $gbShowLabel(kLabel_MRIValue)"\
      gbShowLabel(kLabel_MRIValue) } \
      { check \
      "Parcellation" \
      "ShowLabel kLabel_Parcellation_Name $gbShowLabel(kLabel_Parcellation_Name)"\
      gbShowLabel(kLabel_Parcellation_Name) } }}
      \
      { separator } \
      \
      { cascade "Configure..." { \
      { command "Lighting..." \
      {DoConfigLightingDlog} } \
      \
      { command "Overlay..." \
      {DoConfigOverlayDisplayDlog} \
       mg_OverlayLoaded  } \
      \
      { command "Time Course..." \
      {Graph_DoConfig} \
      mg_TimeCourseLoaded } \
      \
      { command "Curvature Display..." \
      {DoConfigCurvatureDisplayDlog} \
      mg_CurvatureLoaded } \
      \
      { command "Phase Encoded Data Display..." \
      {DoConfigPhaseEncodedDataDisplayDlog} } }}\
      \
      { separator } \
      \
      {cascade "Surface Configuration" { \
      { radio "Main" \
      { set_current_vertex_set $gaLinkedVars(vertexset); \
      UpdateLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(vertexset) \
      0 } \
      \
      { radio "Inflated" \
      { set_current_vertex_set $gaLinkedVars(vertexset); \
      UpdateLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(vertexset) \
      1 \
      mg_InflatedVSetLoaded } \
      \
      { radio "White" \
      { set_current_vertex_set $gaLinkedVars(vertexset); \
      UpdateLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(vertexset) \
      2 \
      mg_WhiteVSetLoaded } \
      \
      { radio "Pial" \
      { set_current_vertex_set $gaLinkedVars(vertexset); \
      UpdateLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(vertexset) \
      3 \
      mg_PialVSetLoaded } \
      \
      { radio "Original" \
      { set_current_vertex_set $gaLinkedVars(vertexset); \
      UpdateLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(vertexset) \
      4 \
      mg_OriginalVSetLoaded } } } \
      \
      {cascade "Overlay Layer" { \
      { radio "Overlay Layer 1" \
      { sclv_set_current_field 0; \
      SendLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(currentvaluefield) \
      0 \
      mg_OverlayLoaded } \
      \
      { radio "Overlay Layer 2" \
      { sclv_set_current_field 1; \
      SendLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(currentvaluefield) \
      1 \
      mg_OverlayLoaded } \
      \
      { radio "Overlay Layer 3" \
      { sclv_set_current_field 2; \
      SendLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(currentvaluefield) \
      2 \
      mg_OverlayLoaded } \
      \
      { radio "Overlay Layer 4" \
      { sclv_set_current_field 3; \
      SendLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(currentvaluefield) \
      3 \
      mg_OverlayLoaded } \
      \
      { radio "Overlay Layer 5" \
      { sclv_set_current_field 4; \
      SendLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(currentvaluefield) \
      4 \
      mg_OverlayLoaded } } } \
      \
      { separator } \
      \
      { check  \
      "Overlay" \
      { SendLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(overlayflag) } \
      \
      { check  \
      "Scale Bar" \
      { SendLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(scalebarflag) } \
      \
      { check  \
      "Color Scale Bar" \
      { SendLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(colscalebarflag) } \
      \
      { check  \
      "Wireframe Overlay" \
      { SendLinkedVarGroup view; UpdateAndRedraw } \
      gaLinkedVars(verticesflag) } }

    # tools menu
    tkm_MakeMenu $mbwTools "Tools" { \
      { command "Save Point" \
      DoSavePoint } \
      \
      { command "Goto Saved Point" \
      DoGotoPoint } \
      \
      { command "Send to Subject..." \
      { DoSendToSubjectDlog } } \
      \
      { separator } \
      \
      { command "Run Script..." \
      { DoFileDlog RunScript } } \
      \
      { separator } \
      \
      { cascade "Labels" { \
      { command "New Label from Marked Vertices" \
      { labl_new_from_marked_vertices; UpdateAndRedraw } } \
      \
      { command "Mark Seleted Label" \
      { labl_mark_vertices $gnSelectedLabel; UpdateAndRedraw } \
      mg_LabelLoaded } \
      \
      { command "Delete Selected Label" \
      { labl_remove $gnSelectedLabel; UpdateAndRedraw } \
      mg_LabelLoaded } } } \
      \
      { cascade "Cut" { \
      { command "Cut Line" \
      { cut_line 0; UpdateAndRedraw } } \
      \
      { command "Cut Closed Line" \
      { cut_line 1; UpdateAndRedraw } } \
      \
      { command "Cut Plane" \
      { cut_plane; UpdateAndRedraw } } \
      \
      { command "Clear Cuts" \
      { restore_ripflags 2; UpdateAndRedraw } } } } \
      \
      { cascade "Graph" { \
      { command "Graph Marked Vertices Avg" \
      { func_select_marked_vertices; func_graph_timecourse_selection} \
      mg_TimeCourseLoaded } \
      \
      { command "Graph Label Avg" \
      { func_select_label; func_graph_timecourse_selection } \
      mg_TimeCourseLoaded } \
      \
      { command "Save Graph to Postscript File" \
      { DoFileDlog SaveGraphToPS } \
      mg_TimeCourseLoaded } } } \
      \
      { cascade "Fill" { \
      { command "Make Fill Boundary" \
      { fbnd_new_line_from_marked_vertices; \
      clear_all_vertex_marks; UpdateAndRedraw } } \
      \
      { command "Delete Selected Boundary" \
      { fbnd_remove_selected_boundary } } \
      \
      { command "Custom Fill..." \
      { DoCustomFillDlog; } } \
      \
      { command "Fill Stats" \
      { floodfill_marked_patch 1; UpdateAndRedraw } \
      mg_OverlayLoaded } \
      \
      { command "Fill Stats" \
      { floodfill_marked_patch 1; UpdateAndRedraw } \
      mg_OverlayLoaded } \
      \
      { command "Fill Curvature" \
      { floodfill_marked_patch 2; UpdateAndRedraw } \
      mg_CurvatureLoaded } } } \
      \
      { cascade "Surface" { \
      { command "Smooth Curvature..." \
      { DoSmoothCurvatureDlog } \
      mg_CurvatureLoaded } \
      \
      { command "Clear Curvature" \
      { clear_curvature } \
      mg_CurvatureLoaded } \
      \
      { command "Smooth Overlay..." \
      { DoSmoothOverlayDlog } \
      mg_OverlayLoaded } \
      \
      { command "Inflate..." \
      { DoInflateDlog } } \
      \
      { command "Swap Surface Fields..." \
      { DoSwapSurfaceFieldsDlog } } \
      \
      { command "Write Decimation..." \
      { DoDecimationDlog } } \
      \
      {command "Write Dipoles..." \
      {DoFileDlog SaveDipolesAs}} \
      \
      { command "Average Background Midpoint" \
      { UpdateLinkedVarGroup cvavg; \
      set gaLinkedVars(cmid) $gaLinkedVars(dipavg); \
      SendLinkedVarGroup cvavg; UpdateAndRedraw } } } } \
      \
      { separator } \
      \
      { command "Save RGB As..." \
      { DoFileDlog SaveRGBAs } } \
      \
      { command "Make Frame" \
      { save_rgb_cmp_frame } } }
      

    pack $mbwFile $mbwEdit $mbwView $mbwTools \
      -side left
}

proc CreateNavigationArea { ifwTop } {

    global gNextTransform

    ResetTransform

    frame $ifwTop
    set fwNoteBook $ifwTop.fwNoteBook
    set fwButtons  $ifwTop.fwButtons

    tixNoteBook $fwNoteBook

    $fwNoteBook add nav1 -label "Historical"
    CreateOldNavigationArea [$fwNoteBook subwidget nav1].fw
    pack [$fwNoteBook subwidget nav1].fw

    $fwNoteBook add nav2 -label "Sensible"
    CreateSimpleNavigationArea [$fwNoteBook subwidget nav2].fw
    pack [$fwNoteBook subwidget nav2].fw

    frame $fwButtons -relief raised -border 2
    tkm_MakeButtons $fwButtons.restore { \
      { text "Restore View" { RestoreView } "" } }
    tkm_MakeButtons $fwButtons.redraw { \
      { text "Redraw View" { UpdateAndRedraw } "" } }

    bind $fwButtons.redraw.bw0 <B2-ButtonRelease> [list UpdateLockButton $fwButtons.redraw.bw0 gaLinkedVars(redrawlockflag)]

    pack $fwButtons.restore $fwButtons.redraw \
      -side left \
      -expand yes \
      -fill x

    pack $fwNoteBook $fwButtons \
      -side top \
      -expand yes \
      -fill x
}

proc CreateSmallNavigationArea { ifwTop } {

    global gaNavigationPane

    ResetTransform

    frame $ifwTop
    set fwNavigation $ifwTop.fwNavigation
    set fwCompressed $ifwTop.fwCompressed
    set fwOld        $ifwTop.fwOld
    set fwButtons    $ifwTop.fwButtons

    CreateSimpleNavigationArea     $fwNavigation
    CreateOldNavigationArea        $fwOld
    CreateCompressedNavigationArea $fwCompressed

    pack $fwCompressed \
      -side top
}

proc CreateSimpleNavigationArea { ifwTop } {

    global gNextTransform

    set knSliderWidth 150
    set knSliderHeight 100

    frame $ifwTop
    
    tkm_MakeBigLabel $ifwTop.rot_label "Rotate"
    tkm_MakeButtons $ifwTop.rot_n {{image icon_arrow_rot_x_pos \
      { rotate_brain_x -$gNextTransform(rotate,degrees); \
      UpdateAndRedraw } ""}}
    tkm_MakeButtons $ifwTop.rot_ne {{image icon_arrow_rot_z_neg \
      { rotate_brain_z $gNextTransform(rotate,degrees); \
      UpdateAndRedraw } ""}}
    tkm_MakeButtons $ifwTop.rot_e {{image icon_arrow_rot_y_neg \
      { rotate_brain_y $gNextTransform(rotate,degrees); \
      UpdateAndRedraw } ""}}
    tkm_MakeButtons $ifwTop.rot_s {{image icon_arrow_rot_x_neg \
      { rotate_brain_x $gNextTransform(rotate,degrees); \
      UpdateAndRedraw } ""}}
    tkm_MakeButtons $ifwTop.rot_w {{image icon_arrow_rot_y_pos \
      { rotate_brain_y -$gNextTransform(rotate,degrees); \
      UpdateAndRedraw } ""}}
    tkm_MakeButtons $ifwTop.rot_nw {{image icon_arrow_rot_z_pos \
      { rotate_brain_z -$gNextTransform(rotate,degrees); \
      UpdateAndRedraw } ""}}
    tkm_MakeSliders $ifwTop.rot_deg [list \
      [list {"Deg"} gNextTransform(rotate,degrees) \
      0 180.0 [expr $knSliderWidth / 2] {} 1 1 horizontal]]
    set nBaseCol 0
    grid $ifwTop.rot_label -column [expr 0 + $nBaseCol] -row 0 -columnspan 3
    grid $ifwTop.rot_n -column [expr 1 + $nBaseCol] -row 1
    grid $ifwTop.rot_ne -column [expr 2 + $nBaseCol] -row 1
    grid $ifwTop.rot_e -column [expr 2 + $nBaseCol] -row 2
    grid $ifwTop.rot_s -column [expr 1 + $nBaseCol] -row 3
    grid $ifwTop.rot_w -column [expr 0 + $nBaseCol] -row 2
    grid $ifwTop.rot_nw -column [expr 0 + $nBaseCol] -row 1
    grid $ifwTop.rot_deg -column [expr 0 + $nBaseCol] -row 4 -columnspan 3


    grid [frame $ifwTop.space1 -width 10] -column 3 -row 0

    tkm_MakeBigLabel $ifwTop.trans_label "Translate"
    tkm_MakeButtons $ifwTop.trans_n {{image icon_arrow_up \
      { translate_brain_y $gNextTransform(translate,dist); \
      UpdateAndRedraw } ""}}
    tkm_MakeButtons $ifwTop.trans_e {{image icon_arrow_right \
      { translate_brain_x $gNextTransform(translate,dist); \
      UpdateAndRedraw } ""}}
    tkm_MakeButtons $ifwTop.trans_s {{image icon_arrow_down \
      { translate_brain_y -$gNextTransform(translate,dist); \
      UpdateAndRedraw } ""}}
    tkm_MakeButtons $ifwTop.trans_w {{image icon_arrow_left \
      { translate_brain_x -$gNextTransform(translate,dist); \
      UpdateAndRedraw } ""}}
    tkm_MakeSliders $ifwTop.trans_mm [list \
      [list {"mm"} gNextTransform(translate,dist) \
      0 100.0 [expr $knSliderWidth / 2] {} 1 1 horizontal]]
    set nBaseCol 4
    grid $ifwTop.trans_label -column [expr 0 + $nBaseCol] -row 0 -columnspan 3
    grid $ifwTop.trans_n -column [expr 1 + $nBaseCol] -row 1
    grid $ifwTop.trans_e -column [expr 2 + $nBaseCol] -row 2
    grid $ifwTop.trans_s -column [expr 1 + $nBaseCol] -row 3
    grid $ifwTop.trans_w -column [expr 0 + $nBaseCol] -row 2
    grid $ifwTop.trans_mm -column [expr 0 + $nBaseCol] -row 4 -columnspan 3

    grid [frame $ifwTop.space2 -width 10] -column 7 -row 0

    tkm_MakeBigLabel $ifwTop.scl_label "Scale"
    tkm_MakeButtons $ifwTop.scl_s {{image icon_zoom_out \
      { scale_brain [expr 1.0 / [expr $gNextTransform(scale,amt)/100.0]]; \
      UpdateAndRedraw } ""}}
    tkm_MakeButtons $ifwTop.scl_b {{image icon_zoom_in \
      { scale_brain [expr $gNextTransform(scale,amt)/100.0]; \
      UpdateAndRedraw } ""}}
    tkm_MakeSliders $ifwTop.scl_per [list \
      [list {"%"} gNextTransform(scale,amt) \
      0 200.0 [expr $knSliderWidth / 2] {} 1 1 horizontal]]
    set nBaseCol 8
    grid $ifwTop.scl_label -column [expr 0 + $nBaseCol] -row 0 -columnspan 3
    grid $ifwTop.scl_s -column [expr 0 + $nBaseCol] -row 2
    grid $ifwTop.scl_b -column [expr 2 + $nBaseCol] -row 2
    grid $ifwTop.scl_per -column [expr 0 + $nBaseCol] -row 4 -columnspan 3
}

proc CreateCompressedNavigationArea { ifwTop } {

    global gNextTransform

    set knSliderWidth 150
    set knSliderHeight 100

    frame $ifwTop
    
    set fwRotation   $ifwTop.fwRotation
    set fwRotButtons $fwRotation.fwRotButtons
    set fwRotSlider  $fwRotation.fwRotSlider

    frame $fwRotation -relief raised -border 2

    tkm_MakeButtons $fwRotButtons { \
      {image icon_arrow_rot_x_pos \
      { rotate_brain_x -$gNextTransform(rotate,degrees); \
      UpdateAndRedraw } ""} \
      \
      {image icon_arrow_rot_z_neg \
      { rotate_brain_z $gNextTransform(rotate,degrees); \
      UpdateAndRedraw } ""} \
      \
      {image icon_arrow_rot_y_neg \
      { rotate_brain_y $gNextTransform(rotate,degrees); \
      UpdateAndRedraw } ""} \
      \
      {image icon_arrow_rot_x_neg \
      { rotate_brain_x $gNextTransform(rotate,degrees); \
      UpdateAndRedraw } ""} \
      \
      {image icon_arrow_rot_y_pos \
      { rotate_brain_y -$gNextTransform(rotate,degrees); \
      UpdateAndRedraw } ""} \
      \
      {image icon_arrow_rot_z_pos \
      { rotate_brain_z -$gNextTransform(rotate,degrees); \
      UpdateAndRedraw } ""} }

    tkm_MakeSliders $fwRotSlider [list \
      [list {"Deg"} gNextTransform(rotate,degrees) \
      0 180.0 [expr $knSliderWidth / 2] {} 1 1 horizontal]]

    pack $fwRotButtons $fwRotSlider -side top


    set fwTranslation   $ifwTop.fwTranslation
    set fwTransButtons  $fwTranslation.fwTransButtons
    set fwTransSlider   $fwTranslation.fwTransSlider

    frame $fwTranslation -relief raised -border 2

    tkm_MakeButtons $fwTransButtons { \
      {image icon_arrow_left \
      { translate_brain_x -$gNextTransform(translate,dist); \
      UpdateAndRedraw } ""} \
      \
      {image icon_arrow_down \
      { translate_brain_y -$gNextTransform(translate,dist); \
      UpdateAndRedraw } ""} \
      \
      {image icon_arrow_up \
      { translate_brain_y $gNextTransform(translate,dist); \
      UpdateAndRedraw } ""} \
      \
      {image icon_arrow_right \
      { translate_brain_x $gNextTransform(translate,dist); \
      UpdateAndRedraw } ""} }

    tkm_MakeSliders $fwTransSlider [list \
      [list {"mm"} gNextTransform(translate,dist) \
      0 100.0 [expr $knSliderWidth / 2] {} 1 1 horizontal]]

    pack $fwTransButtons $fwTransSlider -side top


    set fwScale       $ifwTop.fwScale
    set fwSclButtons  $fwScale.fwSclButtons
    set fwSclSlider   $fwScale.fwSclSlider

    frame $fwScale -relief raised -border 2

    tkm_MakeButtons $fwSclButtons { \
      {image icon_zoom_out \
      { scale_brain [expr 1.0 / [expr $gNextTransform(scale,amt)/100.0]]; \
      UpdateAndRedraw } ""} \
      \
      {image icon_zoom_in \
      { scale_brain [expr $gNextTransform(scale,amt)/100.0]; \
      UpdateAndRedraw } ""} }

    tkm_MakeSliders $fwSclSlider [list \
      [list {"%"} gNextTransform(scale,amt) \
      0 200.0 [expr $knSliderWidth / 2] {} 1 1 horizontal]]

    pack $fwSclButtons $fwSclSlider -side top

    pack $fwRotation $fwTranslation $fwScale -side left
}

proc CreateOldNavigationArea { ifwTop } {

    set knSliderWidth 150
    set knSliderHeight 100

    frame $ifwTop
    
    frame [set fwTransform   $ifwTop.fwTransform]
    frame [set fwRotate      $fwTransform.fwRotate]
    frame [set fwTop         $fwRotate.fwTop]
    frame [set fwBottom      $fwRotate.fwBottom]
    frame [set fwBottomRight $fwRotate.fwBottomRight]
    
    tkm_MakeBigLabel $fwTop.label "Rotate (deg)"
    tkm_MakeSliders $fwTop.y [list \
      [list {"Y"} gNextTransform(rotate,y) \
      180.0 -180.0 $knSliderWidth {} 0 1 horizontal]]
    tkm_MakeSliders $fwBottom.x [list \
    [list {"X"} gNextTransform(rotate,x) \
    180.0 -180.0 $knSliderHeight {} 0 1 vertical]]
    tkm_MakeEntry $fwBottomRight.x "X" gNextTransform(rotate,x) 5
    tkm_MakeEntry $fwBottomRight.y "Y" gNextTransform(rotate,y) 5
    tkm_MakeEntry $fwBottomRight.z "Z" gNextTransform(rotate,z) 5

    pack $fwTop.label $fwTop.y \
      -side top
    pack $fwBottom.x
    pack $fwBottomRight.x $fwBottomRight.y $fwBottomRight.z \
      -side top

    pack $fwTop \
      -side top
    pack $fwBottom \
      -side left
    pack $fwBottomRight \
      -side left \
      -expand yes \
      -fill x
  
    pack $fwRotate  \
      -side left \
      -padx 5 \
      -expand yes \
      -fill x

    frame [set fwTranslate   $fwTransform.fwTranslate]
    frame [set fwTop         $fwTranslate.fwTop]
    frame [set fwBottom      $fwTranslate.fwBottom]
    frame [set fwBottomRight $fwTranslate.fwBottomRight]
    
    tkm_MakeBigLabel $fwTop.label "Translate (mm)"
    tkm_MakeSliders $fwTop.x [list \
      [list {"X"} gNextTransform(translate,x) \
      -100.0 100.0 $knSliderWidth {} 0 1 horizontal]]
    tkm_MakeSliders $fwBottom.y [list \
    [list {"Y"} gNextTransform(translate,y) \
    100.0 -100.0 $knSliderHeight {} 0 1 vertical]]
    tkm_MakeEntry $fwBottomRight.x "X" gNextTransform(translate,x) 5
    tkm_MakeEntry $fwBottomRight.y "Y" gNextTransform(translate,y) 5

    pack $fwTop.label $fwTop.x \
      -side top
    pack $fwBottom.y
    pack $fwBottomRight.x $fwBottomRight.y \
      -side top

    pack $fwTop \
      -side top
    pack $fwBottom \
      -side left
    pack $fwBottomRight \
      -side left \
      -expand yes \
      -fill x
  
    pack $fwTranslate \
      -side left \
      -padx 5 \
      -expand yes \
      -fill x

    frame [set fwScale       $fwTransform.fwScale]
    frame [set fwTop         $fwScale.fwTop]
    frame [set fwLeft        $fwScale.fwLeft]
    frame [set fwRight       $fwScale.fwRight]
    
    tkm_MakeBigLabel $fwTop.label "Scale (%)"

    tkm_MakeSliders $fwLeft.scale [list \
    [list {""} gNextTransform(scale) \
    0.0 200.0 $knSliderHeight {} 0 1 vertical]]
    tkm_MakeEntry $fwRight.scale "" gNextTransform(scale) 5

    pack $fwTop.label \
      -side top
    pack $fwLeft.scale $fwRight.scale \
      -side left

    pack $fwTop \
      -side top 
    pack $fwLeft \
      -side left \
      -expand yes \
      -fill y 
    pack $fwRight \
      -side left \
      -expand yes \
      -fill x
  
    pack $fwScale \
      -side left \
      -expand yes \
      -fill y \
      -padx 5

    frame [set fwButton $ifwTop.fwButton]
    tkm_MakeButtons $fwButton.transform {{ text "Transform" \
      { DoTransform; ResetTransform; UpdateAndRedraw} }}
    pack $fwButton.transform

    pack $fwTransform $fwButton \
      -side top \
      -expand yes \
      -fill x \
      -pady 5
    pack $fwTransform
}

proc CreateCursorFrame { ifwTop } {

    set fwLabel             $ifwTop.fwMainLabel
    set fwLinkCheckbox      $ifwTop.fwLinkCheckbox
    set fwLabels            $ifwTop.fwLabels

    frame $ifwTop

    # the label that goes at the top of the frame
    tkm_MakeBigLabel $fwLabel "Cursor"

    # make the labels
    CreateLabelFrame $fwLabels cursor

    # pack the subframes in a column. 
    pack $fwLabel \
      -side top \
      -anchor w

    pack $fwLabels \
      -side top \
      -anchor s

}

proc CreateMouseoverFrame { ifwTop } {

    global gaLinkedVars
    UpdateLinkedVarGroup mouseover

    set fwTop               $ifwTop.fwTop
    set fwLabel             $fwTop.fwMainLabel
    set fwOn                $fwTop.fwOn
    set fwLabels            $ifwTop.fwLabels

    frame $ifwTop
    frame $fwTop

    # the label that goes at the top of the frame
    tkm_MakeBigLabel $fwLabel "Mouse"

    # make the labels
    CreateLabelFrame $fwLabels mouseover

    # pack the subframes in a column. 
    pack $fwLabel  \
      -side left \
      -anchor w

    pack $fwTop \
      -expand y -fill x \
      -side top \
      -anchor w

    pack $fwLabels \
      -side bottom \
      -anchor s
}

proc CreateLabelFrame { ifwTop iSet } {

    global glLabel gfwaLabel gsaLabelContents

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

  # if it's a value label, make it a normal entry (editable). else 
  # make it an active label (uneditable).
  if { $label == "kLabel_Val" ||         \
    $label == "kLabel_Val2" ||     \
    $label == "kLabel_ValBak" ||   \
    $label == "kLabel_Val2Bak" ||  \
    $label == "kLabel_ValStat" } { 
      tkm_MakeEntry $fwLabel "" gsaLabelContents($label,name) 14
  } else {
      tkm_MakeActiveLabel $fwLabel "" gsaLabelContents($label,name) 14
  }

  # active leabel for the contents (uneditable).
  tkm_MakeActiveLabel $fwValue "" gsaLabelContents($label,value,$iSet) 18
  pack $fwLabel $fwValue \
    -side left \
    -anchor w
    }

}

proc ShowLabel { isLabel ibShow } {

    global gbShowLabel
    PackLabel $isLabel cursor $ibShow 
    PackLabel $isLabel mouseover $ibShow
    set gbShowLabel($isLabel) $ibShow
}

proc ShowValueLabel { inValueIndex ibShow } {
    
    global gaScalarValueID gsaLabelContents
    if { [info exists gaScalarValueID($inValueIndex,label)] == 0 } {
  puts "ShowValueLabel: $inValueIndex invalid"
  return
    }

    ShowLabel $gaScalarValueID($inValueIndex,label) $ibShow
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

proc CreateToolBar { ifwToolBar } {

    global gfwaToolBar

    frame $ifwToolBar

    # main toolbar
    set gfwaToolBar(main)   $ifwToolBar.fwMainBar
    set fwCut              $gfwaToolBar(main).fwCut
    set fwPoint            $gfwaToolBar(main).fwPoint
    set fwView             $gfwaToolBar(main).fwView
    set fwFill             $gfwaToolBar(main).fwFill
    
    frame $gfwaToolBar(main) -border 2 -relief raised
    
    tkm_MakeButtons $fwCut { \
      { image icon_cut_line { cut_line 0; UpdateAndRedraw } \
      "Cut Line" } \
      { image icon_cut_closed_line { cut_line 1; UpdateAndRedraw } \
      "Cut Closed Line" } \
      { image icon_cut_plane { cut_plane; UpdateAndRedraw } \
      "Cut Plane" } \
      { image icon_cut_area { floodfill_marked_patch 0; \
      UpdateAndRedraw } "Fill Uncut Area" } \
      { image icon_cut_clear { restore_ripflags 2; UpdateAndRedraw } \
      "Clear Cuts" } }
    
    tkm_MakeButtons $fwPoint { \
      { image icon_cursor_save { DoSavePoint } "Save Point" } \
      { image icon_cursor_goto { DoGotoPoint } "Goto Saved Point" } }
    
    tkm_MakeButtons $fwView { \
      { image icon_home { RestoreView } "Restore View" } \
      { image icon_redraw { UpdateAndRedraw } "Redraw View" } }
    
    tkm_MakeButtons $fwFill { \
      { image icon_draw_line { fbnd_new_line_from_marked_vertices; \
      UpdateAndRedraw } "Make Fill Boundary" } \
      { image icon_draw_line_closed { close_marked_vertices; \
      fbnd_new_line_from_marked_vertices; \
      UpdateAndRedraw } "Make Closed Fill Boundary" } \
      { image icon_fill { DoCustomFillDlog } "Custom Fill" } \
      { image icon_erase_line { fbnd_remove_selected_boundary; } \
      "Remove Selected Boundary" } }
    
    bind $fwView.bw1 <B2-ButtonRelease> [list UpdateLockButton $fwView.bw1 gaLinkedVars(redrawlockflag)]
    
    pack $fwCut $fwPoint $fwView $fwFill \
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

proc MoveToolWindow { inX inY } {
    wm geometry . +$inX+$inY
    wm deiconify .
}
# ====================================================================== GRAPH

# constants
set ksGraphWindowName "Time Course"
set knLineWidth(active) 4
set knLineWidth(inactive) 2
set knMinWidthPerTick 80

# gGraphSetting($dataSet,{visible,condition,label})
# gConditionData($condition,{points,errors})

set glAllColors {Red Green Blue Purple Brown Pink Gray LightBlue Yellow Orange}
set gnMaxColors [llength $glAllColors]
set nCondition 0
foreach dataSet $glAllColors {
    if { 0 == $nCondition } {
  set gGraphSetting($dataSet,visible)  0
    } else {
  set gGraphSetting($dataSet,visible)  1
    }
    set gGraphSetting($dataSet,condition) $nCondition
    set gGraphSetting($dataSet,label) "Condition $nCondition"
    incr nCondition
}
set glGraphColors [lrange $glAllColors 0 end] ;# the ones in use
set gbErrorBars                     0
set gbTimeCourseOffset              0
set gbShowTimeCourseOffsetOptions   0
set gnMaxNumErrorBars               0
set gnNumDataSets                0
set gnNumTimeCourseConditions    0
set gsTimeCourseLocation ""
set gsTimeCourseDataName ""
set gbAutoRangeGraph 1
set gbPreStimOffset 0
set gnFixedAxesSize(x1) 0
set gnFixedAxesSize(x2) 0
set gnFixedAxesSize(y1) 0
set gnFixedAxesSize(y2) 0

set gbFunctionalWindowOpen 0

proc CreateGraphWindow { iwwTop } {
    global gwwGraphWindow ksGraphWindowName

    set gwwGraphWindow $iwwTop
    toplevel $iwwTop

    wm title $iwwTop $ksGraphWindowName
    wm geometry $iwwTop 600x400
}

proc CreateGraphFrame { ifwGraph } {

    global gwGraph
    global knLineWidth

    frame $ifwGraph

    set gwGraph      $ifwGraph.gwGraph
    set fwNotes      $ifwGraph.fwNotes

    blt::graph $gwGraph -title "Time Course" \
      -plotbackground white
    
    tkm_MakeNormalLabel $fwNotes "Click-1 to set time point. Click-2 and drag to zoom in. Click-3 to unzoom."

    pack $gwGraph        \
      -side top    \
      -fill both   \
      -expand true

    pack $fwNotes -side top

    pack $ifwGraph       \
      -side top    \
      -padx 3      \
      -pady 3      \
      -expand true \
      -fill both


    $gwGraph legend bind all <Enter> { 
  $gwGraph element configure [$gwGraph legend get current] \
    -linewidth $knLineWidth(active)
  $gwGraph legend activate [$gwGraph legend get current]
    }
    
    
    $gwGraph legend bind all <Leave> { 
  $gwGraph element configure [$gwGraph legend get current] \
    -linewidth $knLineWidth(inactive)
  $gwGraph legend deactivate [$gwGraph legend get current]
    }
    
    bind $gwGraph <ButtonPress-2> { Graph_RegionStart %W %x %y }
    bind $gwGraph <B2-Motion> { Graph_RegionMotion %W %x %y }
    bind $gwGraph <ButtonRelease-2> { Graph_RegionEnd %W %x %y }
    bind $gwGraph <ButtonRelease-3> { Graph_Unzoom %W }
}

proc Graph_ShowWindow {} {
    global gwwGraphWindow
    wm deiconify $gwwGraphWindow
}

proc Graph_HideWindow {} {
    global gwwGraphWindow
    wm withdraw $gwwGraphWindow
}

proc Graph_BeginData {} {
    global gConditionData gnNumDataSets
    for { set nCond 0 } { $nCond < $gnNumDataSets } { incr nCond } {
  set gConditionData($nCond,points) {}
  set gConditionData($nCond,errors) {}
    }
}

proc Graph_EndData {} {
    UpdateLinkedVarGroup graph
    Graph_Draw
}

proc Graph_ClearGraph {} {
    global gwGraph
    set lElements [$gwGraph element names *]
    foreach element $lElements {
  $gwGraph element delete $element
    }
}

proc Graph_SetPointsData { inCondition ilPoints } {
    global gConditionData gnNumDataSets
    # save the data. update the number of data sets. 
    set gConditionData($inCondition,points) $ilPoints
  set gnNumDataSets [expr $inCondition + 1]
}

proc Graph_SetErrorData { inCondition ilErrors } {
    global gConditionData
    # save the data
    set gConditionData($inCondition,errors) $ilErrors
}

proc Graph_Draw {} {

    global knMinWidthPerTick
    global gwGraph gnFixedAxesSize
    global gConditionData gGraphSetting glGraphColors
    global gbErrorBars gbAutoRangeGraph
    global gnMaxNumErrorBars gnNumDataSets
    global gaLinkedVars gbPreStimOffset

    Graph_ClearGraph

    # if no data, return
    if {$gnNumDataSets == 0} { return; }
    
    # if there is only one condition, show it (0 is hidden by default)
    if { $gnNumDataSets == 1 } {
  set gGraphSetting(Red,visible) 1
    }
    
    foreach dataSet $glGraphColors {
  
  # skip if not visible.
  if { $gGraphSetting($dataSet,visible) == 0 } { continue; }
  
  # get the condition index for this color
  set nCondition $gGraphSetting($dataSet,condition)
  
  # get the data. continue if we couldn't get it. get its length.
  if { [catch { set lGraphData $gConditionData($nCondition,points) } \
    sResult]} { continue; }
  set nLength [llength $lGraphData]
  if { $nLength <= 0 } { continue; }
  
  # try to size the spacing of the ticks on the x axis appropriatly.
  # get the x range of this data and find out how many points we have.
  # then get the width of the graph and divide it by the number
  # of points. if < knMinWidthPerTick pixels, set the step size to
  # the width divided by knMinWidthPerTick. otherwise set it to half
  # time res (half because the minor tick goes in between each
  # major tick).
  set nNumPoints [expr $nLength / 2]
  set nWidth [$gwGraph cget -width]
  set nWidthPerTick [expr $nWidth / $nNumPoints]
  if { $nWidthPerTick < $knMinWidthPerTick } {
      set nNumMarks [expr $nWidth / $knMinWidthPerTick]
      set nWidthPerTick [expr $nNumPoints / $nNumMarks]
      $gwGraph axis configure x -stepsize $nWidthPerTick
  } else {
      set nWidthPerTick [expr $gaLinkedVars(timeresolution) / 2]
      if { $nWidthPerTick < 1 } {
    set nWidthPerTick 1
      }
      $gwGraph axis configure x -stepsize $nWidthPerTick
  }

  # if we're subtracting the prestim avg..
  if { $gbPreStimOffset && $gaLinkedVars(numprestimpoints) > 0 } {

      # get the sum of all the points before the stim.
      set fPreStimSum 0
      set nNumPreStimPoints $gaLinkedVars(numprestimpoints)
      for { set nTP 0 } { $nTP < $nNumPreStimPoints } { incr nTP } {
    set nIndex [expr [expr $nTP * 2] + 1];
    set fPreStimSum [expr double($fPreStimSum) + double([lindex $lGraphData $nIndex])];
      }

      # find the avg.
      set fPreStimAvg [expr double($fPreStimSum) / double($nNumPreStimPoints)]
      # subtract from all points.
      for { set nTP 0 } { $nTP < $nNumPoints } { incr nTP } {
    set nIndex [expr [expr $nTP * 2] + 1];
    set fOldValue [lindex $lGraphData $nIndex]
    set fNewValue [expr double($fOldValue) - double($fPreStimAvg)]
    set lGraphData [lreplace $lGraphData $nIndex $nIndex $fNewValue]
      }
  }


  # graph the data
  $gwGraph element create line$dataSet \
    -data $lGraphData \
    -color $dataSet \
    -label $gGraphSetting($dataSet,label) \
    -pixels 2 \
    -linewidth 2
  
  # if we're drawing error bars...
  if { 1 == $gbErrorBars } {
      
      # get the data for this condition.
      if { [catch { set lErrors $gConditionData($nCondition,errors) } \
        sResult] } { continue; }
      
      # get the num of errors. save the highest number
      set nLength [llength $lErrors]
      if { $nLength > $gnMaxNumErrorBars } { 
    set gnMaxNumErrorBars $nLength
      }

      # for each value...
      for {set nErrorIndex 0} \
        {$nErrorIndex < $nLength} \
        {incr nErrorIndex} {
    
    # get error amt
    set nError [lindex $lErrors $nErrorIndex]
    
    # if 0, continue
    if { $nError == 0 } { continue; }
    
    # get the index of the data in the graph data list
    set nGraphIndex [expr $nErrorIndex * 2]
    
    # get the x/y coords at this point on the graph
    set nX [lindex $lGraphData $nGraphIndex]
    set nY [lindex $lGraphData [expr $nGraphIndex + 1]];
    
    # draw a graph line from the top to the bottom
    $gwGraph element create error$dataSet$nErrorIndex \
      -data [list $nX [expr $nY - $nError] \
      $nX [expr $nY + $nError] ] \
      -color $dataSet \
      -symbol splus \
      -label "" \
      -pixels 5
      }
  }
    }
    
    # draw some axes
    $gwGraph marker create line \
      -coords [list 0 -Inf 0 Inf] \
      -name xAxis
    $gwGraph marker create line \
      -coords [list -Inf 0 Inf 0] \
      -name yAxis
    
    # if autorange is off, set the min and max of the axes to the saved
    # amounts.
    if { $gbAutoRangeGraph == 0 } {
  $gwGraph axis configure y -min $gnFixedAxesSize(y1) \
    -max $gnFixedAxesSize(y2)
    }
}

proc Graph_SaveToPS { isFileName } {
    global gwGraph
    catch {$gwGraph postscript output $isFileName}
}

proc Graph_UpdateSize { } {
    global gbAutoRangeGraph gwGraph gnFixedAxesSize
    if { $gbAutoRangeGraph == 0 } {
  set gnFixedAxesSize(y1) [lindex [$gwGraph axis limits y] 0]
  set gnFixedAxesSize(y2) [lindex [$gwGraph axis limits y] 1]
    } else {
  $gwGraph axis configure x y -min {} -max {}
    }
}

# zoom callbacks
proc Graph_Zoom { iwGraph inX1 inY1 inX2 inY2 } {

    if { $inX1 < $inX2 } {
  $iwGraph axis configure x -min $inX1 -max $inX2
    } elseif { $inX1 > $inX2 } {
  $iwGraph axis configure x -min $inX2 -max $inX1
    }
    if { $inY1 < $inY2 } {
  $iwGraph axis configure y -min $inY1 -max $inY2
    } elseif { $inY1 > $inY2 } {
  $iwGraph axis configure y -min $inY2 -max $inY1
    }
}
proc Graph_Unzoom { iwGraph } {
    $iwGraph axis configure x y -min {} -max {}
}
proc Graph_RegionStart { iwGraph inX inY } {
    global gnRegionStart
    $iwGraph marker create line -coords { } -name zoomBox \
      -dashes dash -xor yes
    set gnRegionStart(x) [$iwGraph axis invtransform x $inX]
    set gnRegionStart(y) [$iwGraph axis invtransform y $inY]
}
proc Graph_RegionMotion { iwGraph inX inY } {
    global gnRegionStart
    set nX [$iwGraph axis invtransform x $inX]
    set nY [$iwGraph axis invtransform y $inY]
    $iwGraph marker configure zoomBox -coords [list \
      $gnRegionStart(x) $gnRegionStart(y) \
      $gnRegionStart(x) $nY $nX $nY \
      $nX $gnRegionStart(y) \
      $gnRegionStart(x) $gnRegionStart(y)]
}
proc Graph_RegionEnd { iwGraph inX inY } {
    global gnRegionStart
    $iwGraph marker delete zoomBox
    set nX [$iwGraph axis invtransform x $inX]
    set nY [$iwGraph axis invtransform y $inY]
    Graph_Zoom $iwGraph $gnRegionStart(x) $gnRegionStart(y) $nX $nY
}

proc Graph_SetTestData { inNumConditions inNumTimePoints } {

    global gaLinkedVars

    Graph_ShowWindow
    Graph_BeginData
    
    set min [expr 0 - [expr $gaLinkedVars(numprestimpoints) * $gaLinkedVars(timeresolution)]]
    set max [expr [expr $inNumTimePoints * $gaLinkedVars(timeresolution)] + $min]
    for { set cn 0 } { $cn < $inNumConditions } { incr cn } {
  set lData {}
  for { set tp 0 } { $tp < $inNumTimePoints } { incr tp } {
      set second [expr [expr $tp * $gaLinkedVars(timeresolution)] - [expr $gaLinkedVars(numprestimpoints) * $gaLinkedVars(timeresolution)]];
      lappend lData $second [expr [expr $tp + 2] + $cn]
  }
  Graph_SetPointsData $cn $lData
    }
    
    Graph_EndData
}

proc Graph_SetNumConditions { inNumConditions } {

    global gnNumTimeCourseConditions 
    global glAllColors gnMaxColors glGraphColors gGraphSetting

    set gnNumTimeCourseConditions $inNumConditions

    # set the graph colors to the first num_conditions possible colors
    set nLastColorToUse [expr $gnNumTimeCourseConditions - 1]
    if { $nLastColorToUse > [expr $gnMaxColors - 1] } {
  set nLastColorToUse [expr $gnMaxColors - 1]
    }
    set glGraphColors [lrange $glAllColors 0 $nLastColorToUse]
    
    # make sure none of our graph setting conditions are invalid
    foreach dataSet $glGraphColors {
  if { $gGraphSetting($dataSet,condition) >= $gnNumTimeCourseConditions } {
      set gGraphSetting($dataSet,condition) 0
      set gGraphSetting($dataSet,visible)   0
  }
    }

}

proc Graph_DoConfigDlog {} {
    
    global gDialog gbShowTimeCourseOffsetOptions
    global gbErrorBars gbTimeCourseOffset
    global gGraphSetting gnNumTimeCourseConditions glGraphColors 
    global gbPreStimOffset
    global gaLinkedVars
    
    set wwDialog .wwTimeCourseConfigDlog
    
    if { [Dialog_Create $wwDialog "Configure Graph" {-borderwidth 10}] } {
  
  set nMaxCondition [expr $gnNumTimeCourseConditions - 1]
  
  set fwConditions          $wwDialog.fwConditions
  set lfwDisplay            $wwDialog.lfwDisplay
  set fwButtons             $wwDialog.fwButtons
  
  frame $fwConditions
  
  tkm_MakeBigLabel $fwConditions.fwNameLabel "Label"
  tkm_MakeBigLabel $fwConditions.fwVisibleLabel "Visible"
  tkm_MakeBigLabel $fwConditions.fwColorLabel "Color"
  tkm_MakeBigLabel $fwConditions.fwConditionLabel "Condition Shown"
  grid $fwConditions.fwNameLabel -column 0 -row 0 -padx 5
  grid configure $fwConditions.fwNameLabel -sticky w
  grid $fwConditions.fwColorLabel -column 1 -row 0 -padx 5
  grid configure $fwConditions.fwColorLabel -sticky w
  grid $fwConditions.fwVisibleLabel -column 2 -row 0 -padx 5
  grid configure $fwConditions.fwVisibleLabel -sticky w
  grid $fwConditions.fwConditionLabel -column 3 -row 0 -padx 5
  grid configure $fwConditions.fwConditionLabel -sticky w
  
  set nRow 1
  foreach dataSet $glGraphColors {
      
      set fw $fwConditions
      
      # entry for the name
      tkm_MakeEntry $fw.fwEntry$dataSet "" gGraphSetting($dataSet,label) 20
      
      # make a list of the args, then make a list of that list. 
      tkm_MakeNormalLabel $fw.fwLabel$dataSet "$dataSet"
      tkm_MakeCheckboxes $fw.cbVisible$dataSet y \
        [list [list text "" gGraphSetting($dataSet,visible) "Graph_Draw"]]
      
      # this goes on the right, an entry for the condition this
      # color is displaying.
      tkm_MakeEntryWithIncDecButtons \
        $fw.fwControl$dataSet \
        "Condition (0-$nMaxCondition)" \
        gGraphSetting($dataSet,condition) \
        "Graph_Draw" \
        1
      
      grid $fw.fwEntry$dataSet -column 0 -row $nRow -padx 5
      grid $fw.fwLabel$dataSet -column 1 -row $nRow -padx 5
      grid configure $fw.fwLabel$dataSet -sticky w -padx 5
      grid $fw.cbVisible$dataSet -column 2 -row $nRow -padx 5
      grid $fw.fwControl$dataSet -column 3 -row $nRow -padx 5

      incr nRow
  }

  tixLabelFrame $lfwDisplay \
    -label "Display" \
    -labelside acrosstop \
    -options { label.padX 5 }

  set fwDisplaySub          [$lfwDisplay subwidget frame]
  set fwOptions             $fwDisplaySub.fwOptions
  set fwPreStimPoints       $fwDisplaySub.fwPreStimPoints
  set fwTimeRes             $fwDisplaySub.fwTimeRes

  set lOffsetOptions {}
  if { $gbShowTimeCourseOffsetOptions == 1 } {
      set lOffsetOptions [list text "Show percent change" \
        gbTimeCourseOffset \
        {set gbTimeCourseOffset $gbTimeCourseOffset} ]
  } 

  tkm_MakeCheckboxes $fwOptions v [list \
    { text "Show error bars" gbErrorBars \
    {set gbErrorBars $gbErrorBars} } \
    { text "Automatically size graph" gbAutoRangeGraph \
    {set gbAutoRangeGraph $gbAutoRangeGraph } } \
    $lOffsetOptions \
    { text "Subtract pre-stim average" gbPreStimOffset \
    {set gbPreStimOffset $gbPreStimOffset} } ]

  tkm_MakeActiveLabel $fwPreStimPoints \
    "Number of pre-stim points: " gaLinkedVars(numprestimpoints) 5
  tkm_MakeActiveLabel $fwTimeRes \
    "Time resolution: " gaLinkedVars(timeresolution) 5


  tkm_MakeApplyCloseButtons $fwButtons $wwDialog {Graph_UpdateSize;Graph_Draw;}


  pack $fwConditions \
    $lfwDisplay $fwOptions \
    $fwPreStimPoints $fwTimeRes $fwButtons \
    -side top \
    -anchor w \
    -expand yes \
    -fill x
    }
}

proc Graph_ShowOffsetOptions { ibShow } {
    global gbShowTimeCourseOffsetOptions
    puts "Graph_ShowOffsetOptions $ibShow"
    set gbShowTimeCourseOffsetOptions $ibShow
}

# =============================================================== LABEL WINDOW

set ksLabelListWindowName "Labels"
# the number of labels we know about
set gnNumLabels 0
# the strucutre names we know about
set glStructures {}
# structure option widget
set gowStructures ""
set gsColorTableFileName ""
# currently selected label
set gnSelectedLabel 0

proc LblLst_CreateWindow { iwwTop } {

    global gwwLabelListWindow ksLabelListWindowName

    # creates the window and sets its title and initial size
    set gwwLabelListWindow $iwwTop
    toplevel $iwwTop
    wm title $iwwTop $ksLabelListWindowName
    wm geometry $iwwTop 450x325
    wm minsize $iwwTop 450 324
}

proc LblLst_CreateLabelList { ifwList } {

    global glwLabelList gnNumLabels
    global gaLabelInfo
    global glStructures gowStructures
    global gsColorTableFileName
    global gaLinkedVars

    set fwTop             $ifwList
    set fwColorTable      $fwTop.fwColorTable
    set fwLabelFrame      $fwTop.fwLabelFrame
    set fwLabelList       $fwLabelFrame.fwLabelList
    set fwLabelInfo       $fwLabelFrame.fwLabelInfo

    frame $fwTop
    
    frame $fwColorTable
    set fwColorTableFileName $fwColorTable.fwColorTableFileName
    set fwColorTableLoad     $fwColorTable.fwColorTableLoad

    # make the file selector widget and a load button
    set sColorTableFileName ""
    tkm_MakeFileSelector $fwColorTableFileName \
      "Color Table:" gaLinkedVars(colortablename)
    tkm_MakeButtons $fwColorTableLoad [list \
      [list text "Load" {labl_load_color_table $gaLinkedVars(colortablename); UpdateLinkedVarGroup label} ] ]

    pack $fwColorTableFileName \
      -side left \
      -fill x \
      -expand yes
    pack $fwColorTableLoad \
      -side left \
      -anchor e

    # this is the list box of labels
    frame $fwLabelFrame -relief raised -border 2

    tixScrolledListBox $fwLabelList \
      -options { listbox.selectmode single } \
      -browsecmd LblLst_SelectHilitedLabel
    set glwLabelList [$fwLabelList subwidget listbox]
    

    # the info area for the selected label
    frame $fwLabelInfo
    set fwName     $fwLabelInfo.fwName
    set fwVisible  $fwLabelInfo.fwVisible
    set fwLabel    $fwLabelInfo.fwLabel
    set fwInfoButtons  $fwLabelInfo.fwButtons

    tkm_MakeEntry $fwName "Name" gaLabelInfo(name)
    tkm_MakeCheckboxes $fwVisible y [list \
      [list text "Visible" gaLabelInfo(visible)] ]
    tixOptionMenu $fwLabel -label "Structure" \
      -variable gaLabelInfo(structure)
    set gowStructures $fwLabel
    foreach entry $glStructures {
  $gowStructures add command $entry -label $entry
    }

    tkm_MakeButtons $fwInfoButtons [list \
      [list text "Apply" { LblLst_SendCurrentInfo; redraw } \
      "Apply these settings to the selected label" ] \
      [list text "Auto Name" \
      { labl_set_name_from_table $gnSelectedLabel; redraw } \
      "Set the label name from the structure" ] \
      [list text "Mark" { labl_mark_vertices $gnSelectedLabel; redraw } \
      "Mark vertices in this label" ] \
      [list text "Quick Save" \
      { labl_save $gnSelectedLabel [ExpandFileName $gaLabelInfo(name) kFileName_Surface]} \
      "Save this label to label file" ] \
      [list text "Delete" { labl_remove $gnSelectedLabel; redraw } \
      "Delete this label from the list" ] ] y
   
    pack $fwName $fwVisible $fwLabel \
      -side top \
      -expand yes \
      -fill x
    pack  $fwInfoButtons \
      -side top 

    pack $fwLabelList \
      -side left \
      -expand yes \
      -fill both \
      -anchor n
    pack  $fwLabelInfo \
      -side left \
      -expand yes \
      -fill x \
      -anchor n

    pack $fwColorTable \
      -side top \
      -fill x \
      -expand yes \
      -pady 5

    pack $fwLabelFrame \
      -side top \
      -fill both \
      -expand yes

    set gnNumLabels 0
}


proc LblLst_ShowWindow {} {
    global gwwLabelListWindow
    wm deiconify $gwwLabelListWindow
}

proc LblLst_HideWindow {} {
    global gwwLabelListWindow
    wm withdraw $gwwLabelListWindow
}

proc LblLst_UpdateCurrentInfo { isName inStructure ibVisible } {

    global gaLabelInfo
    global glStructures

    # sets the info for the curernt label
    set gaLabelInfo(name) $isName
    set gaLabelInfo(structure) [lindex $glStructures $inStructure]
    set gaLabelInfo(visible) $ibVisible
}

proc LblLst_UpdateInfo { inIndex isName inStructure ibVisible } {

    global gaLabelInfo
    global glStructures
    global gnSelectedLabel
    global glwLabelList

    # delete the list entry in the list box and reinsert it with the
    # new name
    $glwLabelList delete $inIndex
    $glwLabelList insert $inIndex $isName

    # if this is the label currently selected, update the info area too.
    if { $inIndex == $gnSelectedLabel } {
  set gaLabelInfo(name) $isName
  set gaLabelInfo(structure) [lindex $glStructures $inStructure]
  set gaLabelInfo(visible) $ibVisible
    }
}

proc LblLst_SelectHilitedLabel {} {

    global glwLabelList

    # find the hilighted label in the list box and select it
    set nSelection [$glwLabelList curselection]
    if {$nSelection != ""} {
  LblLst_SelectLabel $nSelection
    }
}

proc LblLst_SelectLabel { inIndex } {

    global gnSelectedLabel

    # select this label. this should in turn send us an update 
    # of its information.
    set gnSelectedLabel $inIndex
    labl_select $inIndex
}

proc LblLst_SendCurrentInfo {} {

    global gnSelectedLabel
    global glwLabelList
    global gaLabelInfo
    global gowStructures
    
    # send the contents of the label info.
    set nSelection $gnSelectedLabel
    set nStructure [lsearch [$gowStructures entries] $gaLabelInfo(structure)]
    labl_set_info $nSelection $gaLabelInfo(name) \
      $nStructure $gaLabelInfo(visible)

    # delete the list entry in the list box and reinsert it with the
    # new name
    $glwLabelList delete $nSelection
    $glwLabelList insert $nSelection $gaLabelInfo(name)
}

proc LblLst_AddLabel { isName } {

    global glwLabelList
    global gnNumLabels

    # add a label entry to the end of the list box
    $glwLabelList insert end $isName
    incr gnNumLabels
}

proc LblLst_RemoveLabel { inIndex } {

    global glwLabelList

    # delete the list entry in the list box 
    $glwLabelList delete $inIndex
}

proc LblLst_SetStructures { ilStructures } {

    global glStructures gowStructures

    # set our list of structures.
    set glStructures $ilStructures

    # if we have the widget, delete all the entries and reinsert them.
    if { $gowStructures != "" } {

  foreach entry [$gowStructures entries] {
      $gowStructures delete $entry
  }
  foreach entry $glStructures {
      $gowStructures add command $entry -label $entry
  }
    }
}

# ======================================================== COMPATIBILITY LAYER

proc setfile { iVarName isFileName } {
    upvar $iVarName localvar
    set localvar [ExpandFileName $isFileName]
}

proc save_rgb_named { isName } {

    global rgb named_rgbdir rel_rgbname

    # cp arg->global for uplevel
    set rel_rgbname $name

    # echo cfunc; update abbrev
    uplevel {setfile rgb $named_rgbdir/$rel_rgbname}

    # cfunc
    save_rgb_named_orig $rel_rgbname
}

# ======================================================================= MISC 

proc prompt {} {
 puts "% "
}

proc LoadSurface { isFileName } {
    global insurf hemi ext
    set insurf [ExpandFileName $isFileName kFileName_Surface]
    #read_binary_surf
    UpdateAndRedraw
    set hemi [file rootname [file tail $insurf]]
    set ext [string trimleft [file tail $insurf] $hemi.]
}

proc ExpandFileName { isFileName {iFileType ""} } {

    global session home subject
    global gaFileNameSubDirs

    # look at the first char
    set sFirstChar [string range $isFileName 0 0]
    switch $sFirstChar {
  "/" {
      set sExpandedFileName $isFileName
  }
  "~" {
      if { [string range $isFileName 1 1] == "/" } {
    set sTail [string range $isFileName 2 end]
    set sSubject $subject
    set sSubDir [file dirname $sTail] ;# path portion of file name
    set sFileName [file tail $sTail]  ;# only file name 
    if { $sSubDir == "." } {
        set sExpandedFileName \
          $home/$sSubject/$sFileName 
    } else {
        set sExpandedFileName \
          $home/$sSubject/$sSubDir/$sFileName
    }
      } else {
    set sTail [string range $isFileName 1 end]
    set sSubject [file dirname $sTail]
    if { $sSubject == "." } {
        set sSubject [file tail $sTail]
        set sExpandedFileName \
          $home/$sSubject
    } else {
        set sFileName [file tail $sTail]
        set sExpandedFileName \
          $home/$sSubject/$sFileName
    }
      }
  }
  "*" {
      set sExpandedFileName \
        $session/[string range $isFileName 2 end]
  } 
  "#" {
      if { [info exists env(CSURF_DIR)] } {
    set sExpandedFileName \
           $env(CSURF_DIR)/lib/tcl/[string range $isFileName 2 end]
      } else {
    set sExpandedFileName $isFileName
      }
  }
  default {
      # see if we have a file name type and can find a default
      # subdirectory.
      set sSubDir ""
      if { $iFileType != "" } {
    if { [info exists gaFileNameSubDirs($iFileType)] } {
        set sSubDir $gaFileNameSubDirs($iFileType)
    }
      }
      if { $sSubDir == "" } {
    puts "No expansion found!!!!"
    set sExpandedFileName $isFileName
      } elseif { $sSubDir == "." } {
    set sExpandedFileName \
      $home/$subject/$isFileName 
      } else {
    set sExpandedFileName \
      $home/$subject/$sSubDir/$isFileName
      }
  }
    }
    
    return $sExpandedFileName
}

proc UpdateLockButton { ibwButton ibFlag } {

    upvar #0 $ibFlag bFlag
    if { $bFlag == 0 } {
  set bFlag 1
  $ibwButton config -relief sunken
    } else {
  set bFlag 0
  $ibwButton config -relief raised
    }

    SendLinkedVarGroup redrawlock
}

proc UpdateAndRedraw {} {
    # using after so that dlogs obscuring the graphics window can 
    # close before the redraw.
    after 500 { redraw; }
}

proc RestoreView {} {
    global flag2d
    ResetTransform
    make_lateral_view
    if {$flag2d} {
  restore_zero_position
  rotate_brain_x -90
    }
    UpdateAndRedraw
}

proc DoTransform { } {
    global gNextTransform
    rotate_brain_x $gNextTransform(rotate,x)
    rotate_brain_y $gNextTransform(rotate,y)
    rotate_brain_z $gNextTransform(rotate,z)
    translate_brain_x $gNextTransform(translate,x)
    translate_brain_y $gNextTransform(translate,y)
    scale_brain [expr $gNextTransform(scale)/100.0]
}

proc ResetTransform { } {
    global gNextTransform kanDefaultTransform
    set gNextTransform(rotate,degrees) $kanDefaultTransform(rotate)
    set gNextTransform(translate,dist) $kanDefaultTransform(translate)
    set gNextTransform(scale,amt) $kanDefaultTransform(scale)
}

proc DoSmoothOverlay { inSteps inField } {

    sclv_smooth $inSteps $inField
}

proc DoSmoothCurvature { inSteps } {

    smooth_curv $inSteps
}

proc DoInflate { inSteps } {
    shrink $inSteps
}

proc DoDecimation { isFileName ifSpacing } { 
    if { $ifSpacing > 0.999} {
  subsample_dist $ifSpacing
    } else {
  subsample_orient $ifSpacing
    }
    
    # write decimation
    set dec [ExpandFileName $isFileName kFileName_BEM]
    write_binary_decimation
}

proc DoSavePoint {} {
    find_orig_vertex_coordinates
}

proc DoGotoPoint {} {
    select_orig_vertex_coordinates
}


proc CreateImages {} {

    global ksImageDir

    foreach image_name { icon_cursor_goto icon_cursor_save \
      icon_arrow_up icon_arrow_down icon_arrow_left icon_arrow_right \
      icon_arrow_rot_z_pos icon_arrow_rot_z_neg \
      icon_zoom_in icon_zoom_out \
      icon_arrow_rot_y_neg icon_arrow_rot_y_pos \
      icon_arrow_rot_x_neg icon_arrow_rot_x_pos \
      icon_cut_area icon_cut_closed_line icon_cut_line \
      icon_cut_plane icon_cut_clear \
      icon_draw_line icon_draw_line_closed icon_fill icon_erase_line \
      icon_surface_main icon_surface_original icon_surface_pial \
      icon_home icon_redraw } {

  if { [catch {image create photo  $image_name -file \
    [ file join $ksImageDir $image_name.gif ]} sResult] != 0 } {
      dputs "Error loading $image_name:"
      dputs $sResult
  }
    }
}


proc SetKeyBindings {} {

    # redraw
    bind . <Alt-r> { UpdateAndRedraw }

    # movement 
    bind . <Alt-Up>    { rotate_brain_x 18.0; UpdateAndRedraw }
    bind . <Alt-Down>  { rotate_brain_x -18.0; UpdateAndRedraw }
    bind . <Alt-Right> { rotate_brain_y -18.0; UpdateAndRedraw }
    bind . <Alt-Left>  { rotate_brain_y 18.0; UpdateAndRedraw }
    bind . <Alt-Shift-Left>  { rotate_brain_y 180.0; UpdateAndRedraw }
    bind . <Alt-Shift-Right>  { rotate_brain_y -180.0; UpdateAndRedraw }
    bind . <Alt-KP_Next>  { rotate_brain_z 10.0; UpdateAndRedraw }
    bind . <Alt-KP_Prior> { rotate_brain_z -10.0; UpdateAndRedraw }
    bind . <Alt-KP_Right>   { translate_brain_x 10.0; UpdateAndRedraw }
    bind . <Alt-KP_Left>    { translate_brain_x -10.0; UpdateAndRedraw }
    bind . <Alt-KP_Up>      { translate_brain_y 10.0; UpdateAndRedraw }
    bind . <Alt-KP_Down>    { translate_brain_y -10.0; UpdateAndRedraw }
    bind . <Alt-braceleft>    { scale_brain 0.75; UpdateAndRedraw }
    bind . <Alt-braceright>   { scale_brain 1.25; UpdateAndRedraw }
    bind . <Alt-bracketleft>  { scale_brain 0.95; UpdateAndRedraw }
    bind . <Alt-bracketright> { scale_brain 1.05; UpdateAndRedraw }
}

set tDlogSpecs(LoadSurface) [list \
  -title "Load Surface" \
  -prompt1 "Load Surface:" \
  -default1 [list ExpandFileName "" kFileName_Surface] \
  -note1 "The file name of the surface" \
  -okCmd {LoadSurface %s1} ]
set tDlogSpecs(SaveSurfaceAs) [list \
  -title "Save Surface As" \
  -prompt1 "Save Surface:" \
  -default1 [list ExpandFileName "" kFileName_Surface] \
  -note1 "The file name of the surface to save" \
  -okCmd {set outsurf [ExpandFileName %s1 kFileName_Surface]; \
  CheckFileAndDoCmd $outsurf write_binary_surface;} } ]

set tDlogSpecs(LoadMainSurface) [list \
  -title "Load Main Vertices" \
  -prompt1 "Load Main Vertices:" \
  -default1 [list ExpandFileName "" kFileName_Surface] \
  -note1 "The file name of the surface to load into main vertices" \
  -okCmd {set filename [ExpandFileName %s1 kFileName_Surface]; \
  read_surface_vertex_set 0 $filename} ]
set tDlogSpecs(LoadInflatedSurface) [list \
  -title "Load Inflated Vertices" \
  -prompt1 "Load Inflated Vertices:" \
  -default1 [list ExpandFileName "" kFileName_Surface] \
  -note1 "The file name of the surface to load into inflated vertices" \
  -okCmd {set filename [ExpandFileName %s1 kFileName_Surface]; \
  read_surface_vertex_set 1 $filename} ]
set tDlogSpecs(LoadWhiteSurface) [list \
  -title "Load White Vertices" \
  -prompt1 "Load White Vertices:" \
  -default1 [list ExpandFileName "" kFileName_Surface] \
  -note1 "The file name of the surface to load into white vertices" \
  -okCmd {set filename [ExpandFileName %s1 kFileName_Surface]; \
  read_surface_vertex_set 2 $filename} ]
set tDlogSpecs(LoadPialSurface) [list \
  -title "Load Pial Vertices" \
  -prompt1 "Load Pial Vertices:" \
  -default1 [list ExpandFileName "" kFileName_Surface] \
  -note1 "The file name of the surface to load into pial vertices" \
  -okCmd {set filename [ExpandFileName %s1 kFileName_Surface]; \
  read_surface_vertex_set 3 $filename} ]
set tDlogSpecs(LoadOriginalSurface) [list \
  -title "Load Original Vertices" \
  -prompt1 "Load Original Vertices:" \
  -default1 [list ExpandFileName "" kFileName_Surface] \
  -note1 "The file name of the surface to load into original vertices" \
  -okCmd {set filename [ExpandFileName %s1 kFileName_Surface]; \
  read_surface_vertex_set 4 $filename} ]

set tDlogSpecs(LoadTimeCourse) [list \
  -title "Load Time Course" \
  -prompt1 "Load Volume:" \
  -type1 dir \
  -note1 "The directory containing the binary volume to load" \
  -default1 [list ExpandFileName "" kFileName_Home] \
  -prompt2 "Stem:" \
  -type2 text \
  -note2 "The stem of the binary volume" \
  -prompt3 "Registration File:" \
  -note3 "The file name of the registration file to load. Leave blank to use register.dat in the same directory." \
  -default3 [list ExpandFileName "" kFileName_Home] \
  -okCmd {func_load_timecourse %s1 %s2 %s3;}]

set tDlogSpecs(LoadCurvature) [list \
  -title "Load Curvature" \
  -prompt1 "Load Curvature:" \
  -default1 [list ExpandFileName "" kFileName_Surface] \
  -note1 "The file name of the curvature data" \
  -okCmd {set curv [ExpandFileName %s1 kFileName_Surface]; \
  read_binary_curv; set curvflag 1; RestoreView; } ]
set tDlogSpecs(SaveCurvatureAs) [list \
  -title "Save Curvature As" \
  -prompt1 "Save Curvature:" \
  -default1 [list ExpandFileName "" kFileName_Surface] \
  -note1 "The file name of the curvature data to save" \
  -okCmd {set curv [ExpandFileName %s1 kFileName_Surface]; \
  CheckFileAndDoCmd $curv write_binary_curv} ]

set tDlogSpecs(LoadPatch) [list \
  -title "Load Patch" \
  -prompt1 "Load Patch:" \
  -default1 [list ExpandFileName "" kFileName_Surface] \
  -note1 "The file name of the patch data" \
  -okCmd {set patch [ExpandFileName %s1 kFileName_Surface]; \
  read_binary_patch; RestoreView; UpdateAndRedraw; } ]
set tDlogSpecs(SavePatchAs) [list \
  -title "Save Patch As" \
  -prompt1 "Save Patch:" \
  -default1 [list ExpandFileName "" kFileName_Surface] \
  -note1 "The file name of the patch data to save" \
  -okCmd {set patch [ExpandFileName %s1 kFileName_Surface]; \
  CheckFileAndDoCmd $patch write_binary_patch} ]

set tDlogSpecs(LoadColorTable) [list \
  -title "Load Color Table" \
  -prompt1 "Load Color Table:" \
  -default1 [list ExpandFileName "" kFileName_Label] \
  -note1 "The file name of the color table" \
  -okCmd {labl_load_color_table [ExpandFileName %s1 kFileName_Label]; \
  UpdateLinkedVarGroup label} ]
set tDlogSpecs(LoadLabel) [list \
  -title "Load Label" \
  -prompt1 "Load Label:" \
  -default1 [list ExpandFileName "" kFileName_Label] \
  -note1 "The file name of the label data" \
  -okCmd {labl_load [ExpandFileName %s1 kFileName_Label] }]
set tDlogSpecs(SaveLabelAs) [list \
  -title "Save Selected Label" \
  -prompt1 "Save Selected Label:" \
  -default1 [list ExpandFileName "" kFileName_Label] \
  -note1 "The file name of the label data to save" \
  -okCmd {labl_save $gnSelectedLabel [ExpandFileName %s1 kFileName_Label] }]
set tDlogSpecs(ImportAnnotation) [list \
  -title "Import Annotaion" \
  -prompt1 "Import Annotation:" \
  -default1 [list ExpandFileName "" kFileName_Label] \
  -note1 "The file name of the annotaion" \
  -okCmd {labl_import_annotation [ExpandFileName %s1 kFileName_Label]; \
  UpdateAndRedraw;} ]
set tDlogSpecs(ExportAnnotation) [list \
  -title "Export Annotaion" \
  -prompt1 "Export Annotation:" \
  -default1 [list ExpandFileName "" kFileName_Label] \
  -note1 "The file name of the annotaion to save" \
  -okCmd {labl_export_annotation [ExpandFileName %s1 kFileName_Label]} ]


set tDlogSpecs(SaveDipolesAs) [list \
  -title "Save Dipoles As" \
  -prompt1 "Save Dipoles As:" \
  -default1 [list ExpandFileName "" kFileName_BEM] \
  -note1 "The file name of the dipoles data to save" \
  -okCmd {set dip [ExpandFileName %s1 kFileName_BEM]; \
  CheckFileAndDoCmd $dip write_binary_dipoles;} ]

set tDlogSpecs(LoadFieldSign) [list \
  -title "Load Field Sign" \
  -prompt1 "Load Field Sign:" \
  -default1 [list ExpandFileName "" kFileName_BEM] \
  -note1 "The file name of the field sign data" \
  -okCmd {set fs [ExpandFileName %s1 kFileName_BEM]; \
  read_fieldsign; RestoreView;} ]
set tDlogSpecs(SaveFieldSignAs) [list \
  -title "Save Field Sign As" \
  -prompt1 "Save Field Sign:" \
  -default1 [list ExpandFileName "" kFileName_FieldSign] \
  -note1 "The file name of the field sign data to save" \
  -okCmd {set fs [ExpandFileName %s1 kFileName_FieldSign]; \
  CheckFileAndDoCmd $fs write_fieldsign} ]

set tDlogSpecs(LoadFieldMask) [list \
  -title "Load Field Mask" \
  -prompt1 "Load Field Mask:" \
  -default1 [list ExpandFileName "" kFileName_FieldSign] \
  -note1 "The file name of the field mask data" \
  -okCmd {set fm [ExpandFileName %s1 kFileName_FieldSign]; \
  read_fsmask; RestoreView;} ]
set tDlogSpecs(SaveFieldMaskAs) [list \
  -title "Save Field Mask As" \
  -prompt1 "Save Field Mask:" \
  -default1 [list ExpandFileName "" kFileName_FieldSign] \
  -note1 "The file name of the field mask data to save" \
  -okCmd {set fm [ExpandFileName %s1 kFileName_FieldSign]; \
  CheckFileAndDoCmd $fm write_fsmask} ]

set tDlogSpecs(RunScript) [list \
  -title "Run Script" \
  -prompt1 "Run Script:" \
  -default1 [list ExpandFileName "" kFileName_Script] \
  -note1 "The file name of the TCL script to run" \
  -okCmd {source [ExpandFileName %s1 kFileName_Script]} ]

set tDlogSpecs(SaveGraphToPS) [list \
  -title "Save Time Course" \
  -prompt1 "Save Time Course As:" \
  -note1 "The file name of the PostScript file to create" \
  -default1 [list ExpandFileName "" kFileName_Home] \
  -okCmd {Graph_SaveToPS [ExpandFileName %s1 kFileName_Home]} ]

set tDlogSpecs(SaveRGBAs) [list \
  -title "Save RGB As" \
  -prompt1 "Save RGB:" \
  -default1 [list ExpandFileName "" kFileName_RGB] \
  -note1 "The file name of the RGB file to save" \
  -okCmd {set rgb [ExpandFileName %s1 kFileName_RGB]; save_rgb} ]


proc CheckFileAndDoCmd { iFile iFunction } {

    global gDialog
    
    if { [file exists $iFile] == 0 } {
  $iFunction
    } else {

  set wwDialog .wwOKReplaceDlog
  if { [Dialog_Create $wwDialog "File Exists" {-borderwidth 10}] } {
      
      set fwWarning          $wwDialog.fwWarning
      set fwMessage          $wwDialog.fwMessage
      set fwButtons          $wwDialog.fwButtons
      
      tkm_MakeBigLabel $fwWarning "$iFile exists."
      tkm_MakeNormalLabel $fwMessage "Okay to overwrite?"
      
      # buttons.
      tkm_MakeCancelOKButtons $fwButtons $wwDialog [list $iFunction]
      
      pack $fwWarning $fwMessage $fwButtons \
        -side top       \
        -expand yes     \
        -fill x         \
        -padx 5         \
        -pady 5
  }
    }
}

# ======================================================================= MAIN

CreateImages

set wwTop        .w
set fwMenuBar    $wwTop.fwMenuBar
set fwToolBar    $wwTop.fwToolBar
set fwNavigation $wwTop.fwNavigation
set fwCursor     $wwTop.fwCursor
set fwMouseover  $wwTop.fwMouseover

CreateWindow         $wwTop
CreateMenuBar        $fwMenuBar
CreateToolBar        $fwToolBar
#CreateNavigationArea $fwNavigation
CreateSmallNavigationArea $fwNavigation
CreateCursorFrame    $fwCursor
CreateMouseoverFrame $fwMouseover

pack $fwMenuBar $fwToolBar $fwNavigation \
  -side top \
  -expand true \
  -fill x

pack $fwCursor $fwMouseover \
  -side left \
  -padx 5 \
  -expand true \
  -fill x \
  -fill y \
  -anchor s

pack $wwTop

ShowToolBar main 1

# set up graph window 
CreateGraphWindow .wwGraph
CreateGraphFrame .wwGraph.gwGraph
Graph_UpdateSize
Graph_HideWindow

# set up label list window
LblLst_CreateWindow .wwLabelList
LblLst_CreateLabelList .wwLabelList.lwLabel
pack .wwLabelList.lwLabel \
  -fill both \
  -expand yes
LblLst_HideWindow

# wacky bindings
SetKeyBindings

# enable default labels
ShowLabel kLabel_VertexIndex 1
ShowLabel kLabel_Coords_RAS 1
ShowLabel kLabel_Coords_Tal 1

# make sure window is shoowing
MoveToolWindow 0 0

labl_load_color_table $env(CSURF_DIR)/surface_labels.txt

# we did it!
dputs "Successfully parsed tksurfer.tcl"



#Graph_SetTestData 3 10
