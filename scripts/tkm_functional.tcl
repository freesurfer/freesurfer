#! /usr/bin/wish4.0

# our graph code
source $env(MRI_DIR)/lib/tcl/tkm_graph.tcl

source $env(MRI_DIR)/lib/tcl/wrappers.tcl

# hungarian notation notes

# prefix  type                   library (if applicable)
# ww      window widget          tk
# fw      frame widget           tk
# sw      selector widget        tix
# cw      control widget         tix
# lw      label widget           tk
# lfw     label frame widget     tix
# vw      canvas widget          tk
# gw      graph widget           emu_graph
# ew      entry widget           tk
# lew     label entry widget     tix
# rw      radio button widget    tk
# cbw     checkbutton widget    tk
# bw      button widget          tk
# dw      dialog box widget      tix
# nw      notebook widget        tix
# cbxw    combobox widget       tix
# l       list                   tcl
# n       integer / number
# s       string
# f       float
# g       global

set ksWindowName "Functional Data"
set kLabelFont   $pfont
set kNormalFont  $mfont

# unique identifiers for names
set ksGraphDataName(0) cond0
set ksGraphDataName(1) cond1
set ksGraphDataName(2) cond2
set ksGraphDataName(3) cond3
set ksGraphDataName(4) cond4

# color for conditions
set ksGraphDataColor(0) red
set ksGraphDataColor(1) blue
set ksGraphDataColor(2) green
set ksGraphDataColor(3) purple
set ksGraphDataColor(4) yellow

# the vmar/hmark stuff likes its colors with braces around them
set ksGraphDataColor2(0) {red}
set ksGraphDataColor2(1) {blue}
set ksGraphDataColor2(2) {green}
set ksGraphDataColor2(3) {purple}
set ksGraphDataColor2(4) {yellow}

# our global vars
set gnIsDisplay              1
set gnTimeResolution         0
set gnNumPreStimPoints       0
set gfCurrentTimeCourseValue 0
set gfCurrentOverlayValue    0
set gnTimePoint              0
set gnTimeSecond             0
set gnCondition              0
set gsTimeCoursePath         ""
set gsTimeCourseStem         ""
set gsOverlayPath            ""
set gsOverlayStem            ""
set gfThresholdMin           0
set gfThresholdMid           0
set gfThresholdMax           0
set glGraphData(0)           {}
set glGraphData(1)           {}
set glGraphData(2)           {}
set glGraphData(3)           {}
set glGraphData(4)           {}

set gsGraphName      "timeCourse"
set gnNumberErrorBars 0

# ============================================================== DRAWING GRAPH

proc ClearGraphData {} {

    global ksGraphDataName
    global gsGraphName gcwGraph gnNumberErrorBars

    $gsGraphName clearmark currentTimePointMark

    # for every condition...
    for {set nIndex 0} {$nIndex < 5} {incr nIndex} {
  
  # clear its data and error bars
  $gsGraphName data $ksGraphDataName($nIndex) -coords {}

  for {set nBarIndex 0}                     \
    {$nBarIndex < $gnNumberErrorBars} \
    {incr nBarIndex} {
      set sTag errorBars$nIndex-$nBarIndex
      $gcwGraph delete -withtag $sTag
      }
    }
}

proc SetGraphData { inDataSet ilPoints } {

    global ksGraphDataName ksGraphDataColor ksGraphDataColor2 
    global gsGraphName gcwGraph 
    global glGraphData gnCondition gnTimeSecond gsGraphDataName

    # save our graph data
    set glGraphData($inDataSet) $ilPoints

    # set the points. use the color for the this condition.
    $gsGraphName data $ksGraphDataName($inDataSet)    \
      -colour $ksGraphDataColor($inDataSet) \
      -points 1                             \
      -lines 1                              \
      -coords $ilPoints

    # draw some axes
    $gsGraphName vmark 0 yAxis {black}
    $gsGraphName hmark 0 xAxis {black}
    
    # create the time bar at the current second in the color for the
    # current condition
    $gsGraphName vmark $gnTimeSecond currentTimePointMark \
      $ksGraphDataColor2($gnCondition)

    RecalcCurrentTimeCourseValue
}

proc SetGraphErrorBars { inDataSet ilErrors } {

    global ksGraphDataName ksGraphDataColor ksGraphDataColor2
    global gsGraphName gcwGraph 
    global gnCondition gnTimeSecond glGraphData
    global gnNumberErrorBars

    # for each error...
    set nLength [llength $ilErrors]
    set gnNumberErrorBars $nLength
    for {set nIndex 0} {$nIndex < $nLength} {incr nIndex} {

  # get the error amount.
  set nError [lindex $ilErrors $nIndex]

  # get the x/y coords for the corresponding point.
  set nPointIndex [expr $nIndex * 2]
  set nGraphX [lindex $glGraphData($inDataSet) $nPointIndex]
  set nGraphY [lindex $glGraphData($inDataSet) [expr $nPointIndex + 1]]

  # convert it to canvas coords
  set nCanvasX [$gsGraphName x2canvas $nGraphX]
  set nCanvasY [$gsGraphName y2canvas $nGraphY]

  # draw a vertical line at the point, extending above and beyond
  # the point the amounto f the error. use the right color.
  set sTag errorBars$inDataSet-$nIndex
  $gcwGraph delete -withtag $sTag
  $gcwGraph create line                          \
    $nCanvasX [expr $nCanvasY - $nError]   \
    $nCanvasX [expr $nCanvasY + $nError]   \
    -fill $ksGraphDataColor($inDataSet)    \
    -tag $sTag
    }
}

# called when the user clicks in the canvas.
proc HandleGraphClick { inX inY } {

    global gsGraphName

    # get the x coord of the graph and send it to the c code.
    set nX [$gsGraphName canvas2x $inX]
    SetCurrentFunctionalTimeSecond $nX
}

# ============================================================== UPDATING VARS

proc UpdateTimeCoursePath { isPath } {

    global gsTimeCoursePath

    set gsTimeCoursePath $isPath
}

proc UpdateTimeCourseStem { isStem } {

    global gsTimeCourseStem

    set gsTimeCourseStem $isStem
}

# called by the c code when the time resolution is updated
proc UpdateTimeResolution { inTimeResolution } {

    global gnTimeResolution
    
    set gnTimeResolution $inTimeResolution
}

# called by the c code when the number of prestim points is changed
proc UpdateNumPreStimPoints { inNumPreStimPoints } {

    global gnNumPreStimPoints

    set gnNumPreStimPoints $inNumPreStimPoints
}

proc RecalcCurrentTimeCourseValue {} {

    global glGraphData
    global gnTimePoint gnCondition gfCurrentTimeCourseValue

    # get the index of the y value of the point we're interested in
    set nIndex [expr [expr $gnTimePoint * 2] + 1]

    # pull the value out of our graph data
    set gfCurrentTimeCourseValue [lindex $glGraphData($gnCondition) $nIndex]
}

proc UpdateOverlayPath { isPath } {

    global gsOverlayPath

    set gsOverlayPath $isPath
}

proc UpdateOverlayStem { isStem } {

    global gsOverlayStem

    set gsOverlayStem $isStem
}

# called by the c code when the display flag is updated
proc UpdateDisplayStatus { inIsDisplay } {

    global gnIsDisplay

    set gnIsDisplay $inIsDisplay
}

# called by the c code when the time point is updated.
proc UpdateTimePoint { inTimePoint inTimeSecond } {

    global gnTimePoint gnTimeSecond gsGraphName

    # set our globals
    set gnTimePoint  $inTimePoint
    set gnTimeSecond $inTimeSecond

    # move the current time bar.
    $gsGraphName movevmark currentTimePointMark $gnTimeSecond

    RecalcCurrentTimeCourseValue
}

# called by the c code when the condition is updated
proc UpdateCondition { inCondition } {

    global gnCondition gnTimeSecond gsGraphName
    global ksGraphDataColor2

    # set our global
    set gnCondition $inCondition

    # clear the existing time bar and make a new one with the color
    # of the new condition.
    $gsGraphName clearmark currentTimePointMark
    $gsGraphName vmark $gnTimeSecond currentTimePointMark \
      $ksGraphDataColor2($gnCondition)

    RecalcCurrentTimeCourseValue
}

proc UpdateCurrentOverlayValue { ifValue } {

    global gfCurrentOverlayValue

    set gfCurrentOverlayValue $ifValue
}

# called by the c code when the threshold min is updated
proc UpdateThresholdMin { ifMin } {

    global gfThresholdMin

    # set our global
    set gfThresholdMin $ifMin
}

# called by the c code when the threshold mid is updated
proc UpdateThresholdMid { ifMid } {

    global gfThresholdMid

    # set our global
    set gfThresholdMid $ifMid
}

# called by the c code when the threshold max is updated
proc UpdateThresholdMax { ifMax } {

    global gfThresholdMax

    # set our global
    set gfThresholdMax $ifMax
}

# ========================================================= BUILDING INTERFACE

proc CreateWindow { iwwTop } {

    global ksWindowName

    toplevel $iwwTop
    wm title $iwwTop $ksWindowName
}

proc CreateGraphFrame { ifwGraph } {

    global gcwGraph
    global kLabelFont
    global gsGraphName

    frame $ifwGraph

    set lwTimeCourse $ifwGraph.lwTimeCourse
    set gcwGraph     $ifwGraph.gcwGraph

    label $lwTimeCourse -text "Time Course Graph" \
      -anchor w \
      -font $kLabelFont

    canvas $gcwGraph -width 400 -height 300
    
    emu_graph $gsGraphName -canvas $gcwGraph -width 330 -height 250
    
    pack $lwTimeCourse $gcwGraph -side top \
      -fill both       \
      -expand true

    pack $ifwGraph -side top \
      -padx 3          \
      -pady 3          \
      -expand true     \
      -fill x

    # clicking in the canvas calls HandleGraphClick with the coords clicked
    bind $gcwGraph <ButtonRelease-1> {
  HandleGraphClick %x %y
    }
}

proc CreateTimeCourseFrame { ifwTimeCourse } {

    global kLabelFont kNormalFont
    global gnTimeResolution gnNumPreStimPoints
    global gsTimeCoursePath gsTimeCourseStem gfCurrentTimeCourseValue

    frame $ifwTimeCourse 
    
    set fwLabel $ifwTimeCourse.fwLabel
    set lwLabel $fwLabel.lwLabel

    set fwFileName $ifwTimeCourse.fwFileName
    set lwPath     $fwFileName.lwPath
    set ewPath     $fwFileName.ewPath
    set lwStem     $fwFileName.lwStem
    set ewStem     $fwFileName.ewStem

    set fwTimePoints  $ifwTimeCourse.fwTimePoints
    set lwTimeRes     $fwTimePoints.lwTimeRes
    set ewTimeRes     $fwTimePoints.ewTimeRes
    set lwPreStim     $fwTimePoints.lwPreStim
    set ewPreStim     $fwTimePoints.ewPreStim

    set fwCurrentValue $ifwTimeCourse.fwCurrentValue
    set lwCurrentValue $fwCurrentValue.lwCurrentValue
    set lwActualValue  $fwCurrentValue.lwActualValue


    # the label that goes at the top of the frame
    frame $fwLabel
    label $lwLabel -text "Time Course Volume" \
      -font $kLabelFont
    pack $lwLabel -anchor w

    # the path and stem entries
    frame $fwFileName
    label $lwPath -text "Path" \
      -font $kNormalFont
    entry $ewPath -textvariable gsTimeCoursePath \
      -bd 0
    label $lwStem -text "Stem" \
      -font $kNormalFont
    entry $ewStem -textvariable gsTimeCourseStem \
      -width 3                             \
      -bd 0

    # pack them in a row, with the path entry filling
    pack $lwPath -side left
    pack $ewPath -side left \
      -fill x         \
      -expand true
    pack $lwStem $ewStem -side left

    # time res and prestim entries
    frame $fwTimePoints
    label $lwTimeRes -text "Time Resolution: " \
      -font $kNormalFont                 \
      -anchor e
    entry $ewTimeRes -textvariable gnTimeResolution \
      -width 3                                \
      -bd 0
    label $lwPreStim -text "# Pre-stim Time Pts: " \
      -font $kNormalFont                     \
      -anchor e
    entry $ewPreStim -textvariable gnNumPreStimPoints \
      -width 3                                  \
      -bd 0

    # pack them in a row, to the left
    pack $lwTimeRes $ewTimeRes  $lwPreStim $ewPreStim -side left

    # make the label that will contain the actual current value
    frame $fwCurrentValue
    label $lwCurrentValue -text "Current value:" \
      -font $kNormalFont                   \
      -anchor w
    label $lwActualValue -textvariable gfCurrentTimeCourseValue \
      -font $kNormalFont                   \
      -anchor w
    pack $lwCurrentValue $lwActualValue -side left

    # pack the subframes in a column
    pack $fwLabel $fwFileName $fwTimePoints $fwCurrentValue -side top \
      -expand true                                              \
      -fill x                                                   \
      -anchor w

    # pack the time course frame
    pack $ifwTimeCourse -side top \
      -padx 3               \
      -pady 3               \
      -expand true          \
      -fill x               \
      -anchor w

    # pressing return in the fields calls the c set functions
    bind $ewPath <Return> {
  LoadTimeCourseData $gsTimeCoursePath $gsTimeCourseStem
    }
    bind $ewStem <Return> {
  LoadTimeCourseData $gsTimeCoursePath $gsTimeCourseStem
    }
    bind $ewTimeRes <Return> {
  SetCurrentFunctionalTimeResolution $gnTimeResolution
    }
    bind $ewPreStim <Return> {
  SetCurrentFunctionalNumPreStimPoints $gnNumPreStimPoints
    }
}

proc CreateOverlayFrame { ifwOverlay } {

    global kLabelFont kNormalFont
    global gnIsDisplay gnTimePoint gnCondition
    global gsOverlayPath gsOverlayStem gfCurrentOverlayValue

    frame $ifwOverlay
    
    set fwLabel $ifwOverlay.fwLabel
    set lwLabel $fwLabel.lwLabel

    set fwVisible  $ifwOverlay.fwVisible
    set cbwVisible $fwVisible.cbwVisible

    set fwFileName $ifwOverlay.fwFileName
    set lwPath     $fwFileName.lwPath
    set ewPath     $fwFileName.ewPath
    set lwStem     $fwFileName.lwStem
    set ewStem     $fwFileName.ewStem

    set fwCurrentPlane  $ifwOverlay.fwCurrentPlane
    set lwTimePoint     $fwCurrentPlane.lwTimePoint
    set ewTimePoint     $fwCurrentPlane.ewTimePoint
    set bwIncTimePoint  $fwCurrentPlane.vwIncTimePoint
    set bwDecTimePoint  $fwCurrentPlane.vwDecTimePoint
    set lwCondition     $fwCurrentPlane.lwCondition
    set ewCondition     $fwCurrentPlane.ewCondition
    set bwIncCondition  $fwCurrentPlane.vwIncCondition
    set bwDecCondition  $fwCurrentPlane.vwDecCondition

    set fwCurrentValue $ifwOverlay.fwCurrentValue
    set lwCurrentValue $fwCurrentValue.lwCurrentValue
    set lwActualValue  $fwCurrentValue.lwActualValue

    set fwThreshold     $ifwOverlay.fwThreshold
    set lwThreshold     $fwThreshold.lwThreshold
    set lwMin           $fwThreshold.lwMin
    set ewMin           $fwThreshold.ewMin
    set lwMid           $fwThreshold.lwMid
    set ewMid           $fwThreshold.ewMid
    set lwMax           $fwThreshold.lwMax
    set ewMax           $fwThreshold.ewMax

    # the label that goes at the top of the frame
    frame $fwLabel
    label $lwLabel -text "Overlay Volume" \
      -font $kLabelFont
    pack $lwLabel -anchor w

    # visible checkbutton
    frame $fwVisible -bd 0
    checkbutton $cbwVisible -text "Visible" \
      -variable gnIsDisplay          \
      -font $kNormalFont             \
      -highlightbackground white
    pack $cbwVisible -side left

    # the path and stem entries
    frame $fwFileName
    label $lwPath -text "Path" \
      -font $kNormalFont
    entry $ewPath -textvariable gsOverlayPath \
      -bd 0
    label $lwStem -text "Stem" \
      -font $kNormalFont
    entry $ewStem -textvariable gsOverlayStem \
      -width 3                          \
      -bd 0

    # pack them in a row, with the path entry filling
    pack $lwPath -side left
    pack $ewPath -side left \
      -fill x         \
      -expand true
    pack $lwStem $ewStem -side left

    # cur time point and condition entries with their buttons
    frame $fwCurrentPlane
    label $lwTimePoint -text "Time Point: " \
      -font $kNormalFont \
      -anchor e
    entry $ewTimePoint -textvariable gnTimePoint \
      -width 3                             \
      -bd 0
    button $bwDecTimePoint -text "-" \
      -font $kNormalFont       \
      -padx 2                  \
      -pady 0                  \
      -bd 1
    button $bwIncTimePoint -text "+" \
      -font $kNormalFont       \
      -padx 2                  \
      -pady 0                  \
      -bd 1
    label $lwCondition -text "Condition: " \
      -font $kNormalFont \
      -anchor e
    entry $ewCondition -textvariable gnCondition \
      -width 3                              \
      -bd 0
    button $bwDecCondition -text "-" \
      -font $kNormalFont       \
      -padx 2                  \
      -pady 0                  \
      -bd 1
    button $bwIncCondition -text "+" \
      -font $kNormalFont       \
      -padx 2                  \
      -pady 0                  \
      -bd 1

    # pack them in a row, to the left
    pack $lwTimePoint $ewTimePoint          \
      $bwDecTimePoint $bwIncTimePoint \
      $lwCondition $ewCondition       \
      $bwDecCondition $bwIncCondition -side left \
      -padx 2

    # make the label that will contain the actual current value
    frame $fwCurrentValue
    label $lwCurrentValue -text "Current value:" \
      -font $kNormalFont                   \
      -anchor w
    label $lwActualValue -textvariable gfCurrentOverlayValue \
      -font $kNormalFont                   \
      -anchor w
    pack $lwCurrentValue $lwActualValue -side left

    # threshold min, mid, and max
    frame $fwThreshold
    label $lwThreshold -text "Color scale thresholds: " \
      -font $kNormalFont                          \
      -anchor e
    label $lwMin -text "Min: " \
      -font $kNormalFont \
      -anchor e
    entry $ewMin -textvariable gfThresholdMin \
      -width 6                          \
      -bd 0
    label $lwMid -text "Mid: " \
      -font $kNormalFont \
      -anchor e
    entry $ewMid -textvariable gfThresholdMid \
      -width 6                          \
      -bd 0
    label $lwMax -text "Max: " \
      -font $kNormalFont \
      -anchor e
    entry $ewMax -textvariable gfThresholdMax \
      -width 6                          \
      -bd 0
    
    # pack them in a row, to the left
    pack $lwThreshold $lwMin $ewMin $lwMid $ewMid $lwMax $ewMax -side left \
      -padx 2

    # pack the subframes in a column
    pack $fwLabel $fwVisible $fwFileName     \
      $fwCurrentPlane $fwCurrentValue  \
      $fwThreshold  -side top          \
      -expand true                     \
      -fill x                          \
      -anchor w

    # pack the time course frame
    pack $ifwOverlay -side top  \
      -padx 3             \
      -pady 3             \
      -expand true        \
      -fill x             \
      -anchor w

    # clicking visible changes the value and calls teh c code func
    bind $ewPath <Return> {
  LoadOverlayData $gsOverlayPath $gsOverlayStem
    }
    bind $ewStem <Return> {
  LoadOverlayData $gsOverlayPath $gsOverlayStem
    }
    bind $cbwVisible <ButtonRelease-1> {
  SetCurrentFunctionalDisplayStatus $gnIsDisplay
    }
    bind $ewTimePoint <Return> {
  SetCurrentFunctionalTimePoint $gnTimePoint
    }
    bind $bwDecTimePoint <ButtonRelease-1> {
  SetCurrentFunctionalTimePoint [expr $gnTimePoint - 1]
    }
    bind $bwIncTimePoint <ButtonRelease-1> {
  SetCurrentFunctionalTimePoint [expr $gnTimePoint + 1]
    }
    bind $ewCondition <Return> {
  SetCurrentFunctionalCondition $gnCondition
    }
    bind $bwDecCondition <ButtonRelease-1> {
  SetCurrentFunctionalCondition [expr $gnCondition - 1]
    }
    bind $bwIncCondition <ButtonRelease-1> {
  SetCurrentFunctionalCondition [expr $gnCondition + 1]
    }
    bind $ewMin <Return> {
  SetCurrentFunctionalThresholdMin $gfThresholdMin
    }
    bind $ewMid <Return> {
  SetCurrentFunctionalThresholdMid $gfThresholdMid
    }
    bind $ewMax <Return> {
  SetCurrentFunctionalThresholdMax $gfThresholdMax
    }
}

# = ====================================================================== MAIN

# build the window
set wwTop        .w
set fwGraph      $wwTop.fwGraph
set fwTimeCourse $wwTop.fwTimeCourse
set fwOverlay    $wwTop.fwOverlay

CreateWindow          $wwTop
CreateGraphFrame      $fwGraph
CreateTimeCourseFrame $fwTimeCourse
CreateOverlayFrame    $fwOverlay

set gfwGraph $fwGraph
