#! /usr/bin/tixwish

# our graph code
source $env(MRI_DIR)/lib/tcl/tkm_graph.tcl
source $env(MRI_DIR)/lib/tcl/tkm_dialog.tcl
source $env(MRI_DIR)/lib/tcl/tkm_wrappers.tcl

# hungarian notation notes:
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
# cbw     checkbutton widget     tk
# bw      button widget          tk
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
set ksWindowName "Functional Data"
set kLabelFont   $pfont
set kNormalFont  $mfont
set gsGraphName      "timeCourse"

# FunV_tDisplayFlag
set FunV_tDisplayFlag_Ol_TruncateNegative 0
set FunV_tDisplayFlag_Ol_TruncatePositive 1
set FunV_tDisplayFlag_Ol_ReversePhase     2
set FunV_tDisplayFlag_Ol_OffsetValues     3
set FunV_tDisplayFlag_TC_GraphWindowOpen  4
set FunV_tDisplayFlag_TC_OffsetValues     5
set FunV_tDisplayFlag_TC_PreStimOffset    6

# our global vars
set gnTimePoint      0
set gnTimeSecond     0
set gnCondition      0
set gbTruncateNegative 0
set gbTruncatePositive 0
set gbReverse        0
set gbOverlayOffset  0
set gbShowOverlayOffsetOptions   0
set gfThreshold(min)   0
set gfThreshold(mid)   0
set gfThreshold(slope) 0
set gnOverlayNumTimePoints 0
set gnOverlayNumConditions 0
set gsOverlayDataName ""

set glAllColors {Red Green Blue Purple Yellow Orange Brown Pink Gray LightBlue}
set gnMaxColors [llength $glAllColors]
set nCondition 0
foreach dataSet $glAllColors {
    set gGraphSetting($dataSet,visible)  1
    set gGraphSetting($dataSet,condition) $nCondition
    incr nCondition
    set ksGraphDataColor($dataSet) $dataSet
    set ksGraphDataColor2($dataSet) "{$dataSet}"
}
set glGraphColors [lrange $glAllColors 0 end]
set gbErrorBars                     0
set gbTimeCourseOffset              0
set gbShowTimeCourseOffsetOptions   0
set gnPreStimPoints                 0
set gnTimeResolution                0
set gnMaxNumErrorBars               0
set gnNumDataSets                0
set gnNumTimeCourseConditions    0
set gsTimeCourseLocation ""
set gsTimeCourseDataName ""
set gbAutoRangeGraph 1
set gbPreStimOffset 0

set gbFunctionalWindowOpen 0

# ============================================================== DRAWING GRAPH

proc TimeCourse_BeginDrawingGraph {} {

    TimeCourse_ClearData
}

proc TimeCourse_ClearData {} {

    global ksGraphDataColor
    global gsGraphName gcwGraph gnMaxNumErrorBars gbAutoRangeGraph
    global gcwGraph gsGraphName glGraphColors

    $gsGraphName configure -autorange $gbAutoRangeGraph

    # get size of graph
    set nWidth [winfo width $gcwGraph]
    set nHeight [winfo height $gcwGraph]

    # just recreate the graph. 
    # set size of graph widget. leave room for axes.
    $gsGraphName configure \
      -width [expr $nWidth - 70] \
      -height [expr $nHeight - 50] \
      -nticks_x [expr $nWidth / 30] \
      -nticks_y [expr $nHeight / 30]

    if { [catch {$gsGraphName clearmark currentTimePointMark} sResult] } {}

    # for every condition...
    foreach dataSet $glGraphColors {

  # clear its data and error bars
  $gsGraphName data $dataSet                  \
    -colour $ksGraphDataColor($dataSet) \
    -points 0                           \
    -lines 0                            \
    -coords {0 0 1 1}

  # for each possible error here..
  for {set nErrorIndex 0} {$nErrorIndex < $gnMaxNumErrorBars}\
    {incr nErrorIndex} {
      
      # build a tag for this bar and delete the error bars
      set sTag errorBars$dataSet-$nErrorIndex
      $gcwGraph delete -withtag $sTag
  }
    }
}

proc TimeCourse_UpdateGraphData { inDataSet ilPoints } {

    global glGraphData gnNumDataSets

    # save our graph data
    set glGraphData($inDataSet) $ilPoints

    if {[expr $inDataSet + 1] > $gnNumDataSets} {
  set gnNumDataSets [expr $inDataSet + 1]
    }

}

proc TimeCourse_UpdateErrorData { inDataSet ilErrors } {

    global glErrorData

    set glErrorData($inDataSet) $ilErrors
}

proc TimeCourse_EndDrawingGraph {} {

    TimeCourse_DrawGraph
}

proc TimeCourse_DrawGraph {} {

    global gsGraphName ksGraphDataColor ksGraphDataColor2 knGraphColorIndex
    global glErrorData glGraphData gGraphSetting glGraphColors
    global gnCondition gnTimeSecond
    global gsGraphName gcwGraph gbErrorBars gbAutoRangeGraph
    global gnMaxNumErrorBars gnNumDataSets gnNumTimeCourseConditions

    # if no data, return
    if {$gnNumDataSets == 0} {
  return;
    }
    
    # find the x mins and maxes for our values.
    set nXMin 100000
    set nXMax -100000
    set nYMin 100000
    set nYMax -100000
    foreach dataSet $glGraphColors {
  
  if { $gGraphSetting($dataSet,visible) == 1 } {

      # get the condition index for this color
      set nCondition $gGraphSetting($dataSet,condition)
      
      # make sure we have this data
      if { [catch { set lGraphData $glGraphData($nCondition) } s]} {
    continue;
      }

      # get the number of values
      set nDataLength [llength $lGraphData]
      set nDataPointsLength [expr $nDataLength / 2]

      # go thru values and check for mins and maxes.
      for {set nIndex 0} {$nIndex < $nDataPointsLength} {incr nIndex} {

    # get func value
    set nGraphIndex [expr [expr $nIndex * 2] + 1]
    set nValue [lindex $lGraphData $nGraphIndex];

    # now compare to y min/max
    if { $nValue < $nYMin } { set nYMin $nValue }
    if { $nValue > $nYMax } { set nYMax $nValue }
      }

      # if we are displayinge rror data...
      if { $gbErrorBars == 1 } {
    
    # get the errors for this condition
    if { [catch { set lErrors $glErrorData($nCondition) } s] } {
        continue;
    }  
    
    # get the number of values. make sure it's the same as the
    # number of data values.
    set nErrorLength [llength $lErrors]
    if { $nDataPointsLength != $nErrorLength } {
        continue;
    }
    
    # scan thru. use error length since we want to ignore 
    # the x point values in the graph data.
    for {set nIndex 0} {$nIndex < $nErrorLength} {incr nIndex} {
        
        # get func value
        set nGraphIndex [expr [expr $nIndex * 2] + 1]
        set nValue [lindex $lGraphData $nGraphIndex];

        # get error amt and add/subtract to get a low/hi value
        set nError [lindex $lErrors $nIndex];
        set nLow [expr $nValue - $nError]
        set nHigh [expr $nValue + $nError]

        # now compare to y min/max
        if { $nLow < $nYMin } { set nYMin $nLow }
        if { $nHigh > $nYMax } { set nYMax $nHigh }
    }
      }
  }

  # the x min and max is just the first and last x value in any
  # data set
  set nLength [llength $glGraphData(0)]
  set nXMin [lindex $glGraphData(0) 0]
  set nXMax [lindex $glGraphData(0) [expr $nLength - 2]]
    }
    set nPadding [expr [expr $nYMax - $nYMin] / 15]

    puts "y min $nYMin y max $nYMax padding $nPadding -> [expr $nYMin - $nPadding] [expr $nYMax + $nPadding]"

    # now size the graph, but only if they haven't selected autorange.
    if { $gbAutoRangeGraph == 1 } {
  $gsGraphName configure -autorange 0 \
    -xmin $nXMin \
    -xmax $nXMax \
    -ymin [expr $nYMin - $nPadding] \
    -ymax [expr $nYMax + $nPadding]
    } else {
      $gsGraphName configure -autorange 0
    }

    foreach dataSet $glGraphColors {
  
  # if not visible, draw some dummy data in this color. make
  # sure to use the already established min and max as the dummy
  # data so we don't mess up our range. 
  if { $gGraphSetting($dataSet,visible) == 0 } {
#      set lDummyVals "{$nXMin $nYMin $nXMax $nYMax}"
      set lDummyVals {0 0 1 1}
      $gsGraphName data $dataSet \
        -colour $ksGraphDataColor($dataSet) \
        -points  0 \
        -lines  0 \
        -coords $lDummyVals
      continue;
  }
      
  # get the condition index for this color
  set nCondition $gGraphSetting($dataSet,condition)
  
  # make sure we have this data
  if { [catch { set lGraphData $glGraphData($nCondition) } sResult]} {
      # no data
      continue;
  }
  set nLength [llength $lGraphData]
  if { $nLength <= 0 } {
      continue;
  }

  # graph the data
  $gsGraphName data $dataSet                  \
    -colour $ksGraphDataColor($dataSet) \
    -points 1                           \
    -lines 1                            \
    -coords $lGraphData

  # get the errors for this condition
  if { [catch { set lErrors $glErrorData($nCondition) } sResult] } {
      # no error data
      continue;
  }
  
  # for each error...
  set nLength [llength $lErrors]
  if { $nLength > $gnMaxNumErrorBars } { 
      set gnMaxNumErrorBars $nLength
  }
  for {set nErrorIndex 0} {$nErrorIndex < $nLength} {incr nErrorIndex} {
      
      # build a tag for this bar.
      set sTag errorBars$dataSet-$nErrorIndex

      # actually, if they are hidden, just delete this tag.
      if { $gbErrorBars == 0 } {
    $gcwGraph delete -withtag $sTag
    continue;
      }
  
      # actually, if this color is not visible, just delete this tag.
      if { $gGraphSetting($dataSet,visible) == 0 } {
    $gcwGraph delete -withtag $sTag
    continue;
      }
  
      # get error amt
      set nError [lindex $lErrors $nErrorIndex]
      
      # get the x/y coords at this point on the graph
      set nPointIndex [expr $nErrorIndex * 2]
      set nGraphX [lindex $lGraphData $nPointIndex]
      set nGraphY [lindex $lGraphData [expr $nPointIndex + 1]]
      
      # get the y bounds, with errors, then cap to graph
      set nGraphY1 [expr $nGraphY - $nError]
      set nGraphY2 [expr $nGraphY + $nError]
      if { $nGraphY1 < [$gsGraphName cget ymin] } {
    set nGraphY1 [$gsGraphName cget ymin]
      }
      if { $nGraphY2 > [$gsGraphName cget ymax] } {
    set nGraphY2 [$gsGraphName cget ymax]
      }
      
      # convert it to canvas coords
      set nCanvasX  [$gsGraphName x2canvas $nGraphX]
      set nCanvasY1 [$gsGraphName y2canvas $nGraphY1]
      set nCanvasY2 [$gsGraphName y2canvas $nGraphY2]
      
      # draw a vertical line at this point, extending above and below
      # it, using the right color.
      $gcwGraph delete -withtag $sTag
      $gcwGraph create line \
        [expr $nCanvasX - 5] $nCanvasY1 \
        [expr $nCanvasX + 5] $nCanvasY1 \
        $nCanvasX $nCanvasY1 \
        $nCanvasX $nCanvasY2 \
        [expr $nCanvasX - 5] $nCanvasY2 \
        [expr $nCanvasX + 5] $nCanvasY2 \
        -fill $ksGraphDataColor($dataSet) \
        -tag $sTag
  }
    }
    
    # draw some axes
    $gsGraphName vmark 0 yAxis black
    $gsGraphName hmark 0 xAxis black

    TimeCourse_DrawCurrentTimePoint
}

proc TimeCourse_DrawCurrentTimePoint {} {

    global gsGraphName gcwGraph gnTimeSecond

    # delete the old one
    $gcwGraph delete currentTimePoint

    # draw a dashed line at current time point
    set nX [$gsGraphName x2canvas $gnTimeSecond]
    set nYTop [$gsGraphName y2canvas [$gsGraphName cget ymax]]
    set nYBottom [$gsGraphName y2canvas [$gsGraphName cget ymin]]
    set nYStep [expr [expr $nYTop - $nYBottom] / 10]
    set nYCurrent $nYBottom
    for {set nStep 0} {$nStep < 5} {incr nStep } {
  $gcwGraph create line \
    $nX $nYCurrent \
    $nX [expr $nYCurrent + $nYStep] \
    -fill red \
    -tag currentTimePoint
  incr nYCurrent [expr $nYStep * 2]
    }
}

# called when the user clicks in the canvas.
proc HandleGraphClick { inX inY } {

    global gsGraphName

    # get the x coord of the graph and send it to the c code.
    set nSecond [$gsGraphName canvas2x $inX]
}

# ============================================================== UPDATING VARS

proc ShowFunctionalWindow {} {

    global gwwTop gbFunctionalWindowOpen
    wm deiconify $gwwTop
    set gbFunctionalWindowOpen 1
}

proc HideFunctionalWindow {} {

    global gwwTop gbFunctionalWindowOpen
    wm withdraw $gwwTop
    set gbFunctionalWindowOpen 0
}

proc Overlay_DoConfigDlog {} { 

    global gDialog
    global FunV_tDisplayFlag_Ol_TruncateNegative 
    global FunV_tDisplayFlag_Ol_TruncatePositive 
    global FunV_tDisplayFlag_Ol_ReversePhase
    global gnOverlayNumTimePoints gnOverlayNumConditions
    global gnTimePoint gnCondition gbOverlayOffset gbShowOverlayOffsetOptions
    global gbTruncateNegative gbReverse gbTruncatePositive
    global gfThreshold

    set wwDialog .wwOverlayConfigDlog

    if { [Dialog_Create $wwDialog "Configure Functional Overlay" {-borderwidth 10}] } {

  set fwTimePoint       $wwDialog.fwTimePoint
  set fwCondition       $wwDialog.fwCondition
  set fwTruncate        $wwDialog.fwTruncate
  set fwTruncateNeg     $wwDialog.fwTruncateNeg
  set fwReverse         $wwDialog.fwReverse
  set fwOffset          $wwDialog.fwOffset
  set fwThresholdMin    $wwDialog.fwThresholdMin
  set fwThresholdMid    $wwDialog.fwThresholdMid
  set fwThresholdSlope  $wwDialog.fwThresholdSlope
  set fwButtons         $wwDialog.fwButtons

  Overlay_SaveConfiguration;

  set nMaxCondition [expr $gnOverlayNumConditions - 1]
  set nMaxTimePoint [expr $gnOverlayNumTimePoints - 1]

  tkm_MakeEntryWithIncDecButtons \
    $fwTimePoint "Time Point (0-$nMaxTimePoint)" \
    gnTimePoint \
    {} 1

  tkm_MakeEntryWithIncDecButtons \
    $fwCondition "Condition (0-$nMaxCondition)" \
    gnCondition \
    {} 1

  tkm_MakeCheckbox $fwTruncate "Truncate negative values" \
    gbTruncateNegative \
    "set gbTruncateNegative \$gbTruncateNegative"

  tkm_MakeCheckbox $fwTruncateNeg "Truncate positive values" \
    gbTruncatePositive \
    "set gbTruncatePositive \$gbTruncatePositive"

  tkm_MakeCheckbox $fwReverse "Reverse values" gbReverse \
    "set gbReverse \$gbReverse"

  if { $gbShowOverlayOffsetOptions == 1 } {
      tkm_MakeCheckbox $fwOffset "Show percent change" gbOverlayOffset \
        "set gbOverlayOffset \$gbOverlayOffset"
  } else {
      frame $fwOffset
  }

  tkm_MakeEntry $fwThresholdMin "Threshold minimum" gfThreshold(min) 6
  tkm_MakeEntry $fwThresholdMid "Threshold midpoint" gfThreshold(mid) 6
  tkm_MakeEntry $fwThresholdSlope "Threshold slope" gfThreshold(slope) 6

  tkm_MakeCancelApplyOKButtons $fwButtons $wwDialog \
    { Overlay_SetConfiguration; } \
    { Overlay_RestoreConfiguration; }
  
  pack $fwTimePoint $fwCondition $fwTruncate $fwTruncateNeg \
    $fwReverse $fwOffset \
    $fwThresholdMin $fwThresholdMid $fwThresholdSlope \
    $fwButtons \
    -side top \
    -anchor w \
    -expand yes \
    -fill x
    }
}

proc TimeCourse_DoConfigDlog {} {

    global gDialog gbShowTimeCourseOffsetOptions
    global gbErrorBars gnPreStimPoints gnTimeResolution gbTimeCourseOffset
    global gGraphSetting gnNumTimeCourseConditions glGraphColors 
    global gbPreStimOffset

    set wwDialog .wwTimeCourseConfigDlog

    if { [Dialog_Create $wwDialog "Configure Time Course" {-borderwidth 10}] } {

  set nMaxCondition [expr $gnNumTimeCourseConditions - 1]

  set fwConditions          $wwDialog.fwConditions
  set fwLabels              $fwConditions.fwLabels
  set fwVisibleLabel        $fwLabels.fwVisibleLabel
  set fwConditionLabel      $fwLabels.fwConditionLabel
  set fwErrorBars           $wwDialog.fwErrorBars
  set fwAutoRange           $wwDialog.fwAutoRange
  set fwOffset              $wwDialog.fwOffset
  set fwPreStimOffset       $wwDialog.fwPreStimOffset
  set fwPreStimPoints       $wwDialog.fwPreStimPoints
  set fwTimeRes             $wwDialog.fwTimeRes
  set fwButtons             $wwDialog.fwButtons

  TimeCourse_SaveConfiguration;

  frame $fwConditions

  frame $fwLabels
  tkm_MakeBigLabel $fwVisibleLabel "Visible"
  pack $fwVisibleLabel \
    -side left \
    -anchor w
  tkm_MakeBigLabel $fwConditionLabel "Condition Shown"
  pack $fwConditionLabel \
    -side right \
    -anchor e
  pack $fwLabels \
    -side top \
    -expand yes \
    -fill x

  foreach dataSet $glGraphColors {

      set fw $fwConditions.fwCondition$dataSet
      frame $fw
      tkm_MakeCheckbox $fw.cbVisibile \
        "" gGraphSetting($dataSet,visible) \
        "TimeCourse_SetGraphSetting $dataSet visible \$gGraphSetting($dataSet,visible)"
      tkm_MakeNormalLabel $fw.lwLabel "$dataSet"
      tkm_MakeEntryWithIncDecButtons \
        $fw.fwControl \
        "Condition (0-$nMaxCondition)" \
        gGraphSetting($dataSet,condition) \
        "TimeCourse_SetGraphSetting $dataSet condition" \
        1
      pack $fw.cbVisibile $fw.lwLabel \
        -side left \
        -anchor w
      pack $fw.fwControl \
        -side right
      pack $fw \
        -expand yes \
        -fill x
  }


  tkm_MakeCheckbox \
    $fwErrorBars "Show error bars" gbErrorBars \
    {set gbErrorBars $gbErrorBars}

  tkm_MakeCheckbox \
    $fwAutoRange "Automatically size graph" gbAutoRangeGraph \
    {set gbAutoRangeGraph $gbAutoRangeGraph}

  if { $gbShowTimeCourseOffsetOptions == 1 } {
      tkm_MakeCheckbox \
        $fwOffset "Show percent change" gbTimeCourseOffset \
        {set gbTimeCourseOffset $gbTimeCourseOffset}
  } else {
      frame $fwOffset
  }

  tkm_MakeCheckbox \
    $fwPreStimOffset "Subtract pre-stim average" gbPreStimOffset \
    {set gbPreStimOffset $gbPreStimOffset}

  tkm_MakeEntryWithIncDecButtons \
    $fwPreStimPoints "Number of pre-stim points" gnPreStimPoints \
    {} 1

  tkm_MakeActiveLabel \
    $fwTimeRes "Time resolution" gnTimeResolution

  tkm_MakeCancelApplyOKButtons $fwButtons $wwDialog \
    { TimeCourse_SetConfiguration; TimeCourse_DrawGraph } \
    { TimeCourse_RestoreConfiguration; TimeCourse_DrawGraph }

  pack $fwConditions $fwErrorBars $fwAutoRange $fwOffset \
    $fwPreStimOffset $fwPreStimPoints $fwTimeRes $fwButtons \
    -side top \
    -anchor w \
    -expand yes \
    -fill x
    }
}


# ================================================== MANAGING INTERNAL STUFF

proc TimeCourse_SetGraphSetting { iColor iType iValue } {
    global gGraphSetting gnNumTimeCourseConditions
    
    if { $iType == "condition" } {
  if { $iValue >= $gnNumTimeCourseConditions } {
      ErrorDlog "Invalid condition for $iColor. Must be between 0 and $gnNumTimeCourseConditions."
      return;
  }
  if { $iValue < 0 } {
      ErrorDlog "Invalid condition for $iColor. Must be between 0 and $gnNumTimeCourseConditions."
      return;
  }
    }

    set gGraphSetting($iColor,$iType) $iValue
}

proc Overlay_SaveConfiguration {} {

    global gnTimePoint gnCondition gbTruncateNegative gbReverse gfThreshold 
    global gbOverlayOffset gbTruncatePositive
    global gnSavedTimePoint gnSavedCondition gbSavedTruncateNegative 
    global gbSavedTruncatePositive
    global gbSavedReverse gfSavedThreshold gbSavedOverlayOffset
    
    set gnSavedTimePoint $gnTimePoint
    set gnSavedCondition $gnCondition
    set gbSavedTruncateNegative $gbTruncateNegative
    set gbSavedTruncatePositive $gbTruncatePositive
    set gbSavedReverse $gbReverse
    set gbSavedOverlayOffset $gbOverlayOffset
    foreach entry {min mid slope} {
  set gfSavedThreshold($entry) $gfThreshold($entry)
    }


}

proc Overlay_RestoreConfiguration {} {

    global gnTimePoint gnCondition gbTruncateNegative gbReverse gfThreshold
    global gbOverlayOffset gbTruncatePositive
    global gnSavedTimePoint gnSavedCondition gbSavedTruncateNegative
    global gbSavedTruncatePositive
    global gbSavedReverse gfSavedThreshold gbSavedOverlayOffset
    
    set gnTimePoint $gnSavedTimePoint
    set gnCondition $gnSavedCondition
    set gbTruncateNegative $gbSavedTruncateNegative
    set gbTruncatePositive $gbSavedTruncatePositive
    set gbReverse $gbSavedReverse
    set gbTruncateNegative $gbSavedTruncateNegative
    set gbOverlayOffset $gbSavedOverlayOffset
    foreach entry {min mid slope} {
  set gfThreshold($entry) $gfSavedThreshold($entry)
    }

    Overlay_SetConfiguration
}

proc Overlay_SetConfiguration {} {

    global gnTimePoint gnCondition gbTruncateNegative gbReverse gfThreshold
    global gbOverlayOffset gbTruncatePositive
    global FunV_tDisplayFlag_Ol_TruncateNegative
    global FunV_tDisplayFlag_Ol_TruncatePositive
    global FunV_tDisplayFlag_Ol_ReversePhase
    global FunV_tDisplayFlag_Ol_OffsetValues

    Overlay_SetTimePoint $gnTimePoint
    Overlay_SetCondition $gnCondition
    Overlay_SetDisplayFlag $FunV_tDisplayFlag_Ol_TruncateNegative $gbTruncateNegative
    Overlay_SetDisplayFlag $FunV_tDisplayFlag_Ol_ReversePhase $gbReverse
    Overlay_SetDisplayFlag $FunV_tDisplayFlag_Ol_TruncatePositive $gbTruncatePositive
    Overlay_SetDisplayFlag $FunV_tDisplayFlag_Ol_OffsetValues $gbOverlayOffset
    Overlay_SetThreshold $gfThreshold(min) $gfThreshold(mid) $gfThreshold(slope)
}

proc TimeCourse_SaveConfiguration {} {

    global gbErrorBars gnPreStimPoints gnTimeResolution gGraphSetting
    global gbTimeCourseOffset gbPreStimOffset
    global gbSavedErrorBars gnSavedPreStimPoints gnSavedTimeResolution
    global gSavedGraphSetting glGraphColors gbSavedTimeCourseOffset
    global gbSavedPreStimOffset

    set gbSavedErrorBars $gbErrorBars
    set gnSavedPreStimPoints $gnPreStimPoints
    set gnSavedTimeResolution $gnTimeResolution
    set gbSavedTimeCourseOffset $gbTimeCourseOffset
    set gbSavedPreStimOffset $gbPreStimOffset
    foreach color $glGraphColors {
  foreach entry {visible condition} {
      set gSavedGraphSetting($color,$entry) $gGraphSetting($color,$entry)
  }
    }
}

proc TimeCourse_RestoreConfiguration {} {

    global gbErrorBars gnPreStimPoints gnTimeResolution gGraphSetting
    global gbTimeCourseOffset gbPreStimOffset
    global gbSavedErrorBars gnSavedPreStimPoints gnSavedTimeResolution
    global gSavedGraphSetting glGraphColors gbSavedTimeCourseOffset
    global gbSavedPreStimOffset

    set gbErrorBars $gbSavedErrorBars
    set gnPreStimPoints $gnSavedPreStimPoints
    set gnTimeResolution $gnSavedTimeResolution
    set gbTimeCourseOffset $gbSavedTimeCourseOffset
    set gbPreStimOffset $gbSavedPreStimOffset
    foreach color $glGraphColors {
  foreach entry {visible condition} {
      set gGraphSetting($color,$entry) $gSavedGraphSetting($color,$entry)
  }
    }

    TimeCourse_SetConfiguration;
}

proc TimeCourse_SetConfiguration {} {

    global gbErrorBars gnPreStimPoints gnTimeResolution gbTimeCourseOffset
    global gbPreStimOffset
    global FunV_tDisplayFlag_TC_OffsetValues FunV_tDisplayFlag_TC_PreStimOffset

    TimeCourse_SetErrorBarsFlag $gbErrorBars
    TimeCourse_SetNumPreStimPoints $gnPreStimPoints
    TimeCourse_SetTimeResolution $gnTimeResolution
    TimeCourse_SetDisplayFlag $FunV_tDisplayFlag_TC_OffsetValues $gbTimeCourseOffset
    TimeCourse_SetDisplayFlag $FunV_tDisplayFlag_TC_PreStimOffset $gbPreStimOffset
}
# =======================================================================


# ============================================================= C CALLBACKS

proc Overlay_UpdateNumTimePoints { inNumPts } { 
    global gnOverlayNumTimePoints
    set gnOverlayNumTimePoints $inNumPts
}

proc Overlay_UpdateNumConditions { inNumConditions } { 
    global gnOverlayNumConditions
    set gnOverlayNumConditions $inNumConditions
}

proc Overlay_UpdateDisplayFlag { iFlag ibNewValue } { 
    
    global FunV_tDisplayFlag_Ol_ReversePhase 
    global FunV_tDisplayFlag_Ol_TruncateNegative
    global FunV_tDisplayFlag_Ol_TruncatePositive

    global gbTruncateNegative gbReverse gbTruncatePositive

    if { $FunV_tDisplayFlag_Ol_ReversePhase == $iFlag } {
  set gbReverse $ibNewValue
    }
    if { $FunV_tDisplayFlag_Ol_TruncateNegative == $iFlag } {
  set gbTruncateNegative $ibNewValue
    }
    if { $FunV_tDisplayFlag_Ol_TruncatePositive == $iFlag } {
  set gbTruncatePositive $ibNewValue
    }
}

proc Overlay_UpdateDataName { isName } { 
    global gsOverlayDataName
    set gsOverlayDataName $isName
} 

proc Overlay_UpdateTimePoint { inTimePoint inTimeSecond } {

    global gnTimePoint gnTimeSecond gbFunctionalWindowOpen
    set gnTimePoint  $inTimePoint
    set gnTimeSecond $inTimeSecond

    if { $gbFunctionalWindowOpen == 1 } {
  TimeCourse_DrawCurrentTimePoint
    }
}

proc Overlay_UpdateCondition { inCondition } {
    global gnCondition 
    set gnCondition $inCondition
}

proc Overlay_UpdateThreshold { ifMin ifMid ifSlope } {
    global gfThreshold
    set gfThreshold(min) $ifMin
    set gfThreshold(mid) $ifMid
    set gfThreshold(slope) $ifSlope
}

proc Overlay_ShowOffsetOptions { ibShow } {

    global gbShowOverlayOffsetOptions
    set gbShowOverlayOffsetOptions $ibShow
}

proc TimeCourse_UpdateNumConditions { inNumConditions } { 
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

proc TimeCourse_UpdateNumPreStimPoints { inNumPts } { 
    global gnPreStimPoints
    set gnPreStimPoints $inNumPts
 }

proc TimeCourse_UpdateTimeResolution { inTimeRes } { 
    global gnTimeResolution
    set gnTimeResolution $inTimeRes
 }

proc TimeCourse_UpdateDisplayFlag { iFlag ibNewValue } { 

    global FunV_tDisplayFlag_TC_GraphWindowOpen
    
    if { $FunV_tDisplayFlag_TC_GraphWindowOpen == $iFlag } {
  if { $ibNewValue == 0 } {
      HideFunctionalWindow
  } else {
      ShowFunctionalWindow
  }
    }
}

proc TimeCourse_UpdateDataName { isName } { 
    global gwwTop gsTimeCourseDataName
    wm title $gwwTop $isName
    set gsTimeCourseDataName $isName
}

proc TimeCourse_UpdateLocationName { isName } { 
    global gsTimeCourseLocation
    set gsTimeCourseLocation $isName
}

proc TimeCourse_SetErrorBarsFlag { ibNewValue } {
    global gbErrorBars
    set gbErrorBars $ibNewValue
}

proc TimeCourse_ShowOffsetOptions { ibShow } {

    global gbShowTimeCourseOffsetOptions
    set gbShowTimeCourseOffsetOptions $ibShow
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
    global gsGraphName gsTimeCourseDataName gsTimeCourseLocation

    frame $ifwGraph

    set fwDataName   $ifwGraph.fwDataName
    set gcwGraph     $ifwGraph.gcwGraph
    set fwLocation   $ifwGraph.fwLocation

    tkm_MakeActiveLabel $fwDataName "" gsTimeCourseDataName

    canvas $gcwGraph -width 400 -height 300
    
    bind $gcwGraph <Configure> { ResizeGraphFrame }

    emu_graph $gsGraphName -canvas $gcwGraph -width 330 -height 250
    
    tkm_MakeActiveLabel $fwLocation "" gsTimeCourseLocation
    
    pack $fwDataName    \
      -side top     \
      -expand false \
      -anchor w

    pack $gcwGraph \
      -side top            \
      -fill both           \
      -expand true

    pack $ifwGraph       \
      -side top    \
      -padx 3     \
      -pady 3      \
      -expand true \
      -fill both

    pack $fwLocation   \
      -side top     \
      -expand false \
      -anchor w

    # clicking in the canvas calls HandleGraphClick with the coords clicked
    bind $gcwGraph <ButtonRelease-1> {
  HandleGraphClick %x %y
    }
}

proc ResizeGraphFrame {} {

    global gcwGraph gsGraphName gwwTop

    # redraw graph
    TimeCourse_ClearData
    TimeCourse_DrawGraph
}

# = ====================================================================== MAIN

# build the window
set gwwTop       .wwFunctional
set fwGraph      $gwwTop.fwGraph

CreateWindow          $gwwTop
CreateGraphFrame      $fwGraph

set gfwGraph $fwGraph

HideFunctionalWindow

#TimeCourse_DoConfigDlog
#Overlay_DoConfigDlog
#TestData

proc TestData {} {

    TimeCourse_UpdateGraphData 0 {1 1 2 2 3 3 4 4 5 5}
    TimeCourse_UpdateGraphData 1 {1 8 2 9 3 6 4 8 5 9}
    TimeCourse_UpdateErrorData 1 {.9 .8 .7 .6 .5}
    TimeCourse_UpdateGraphData 2 {1 -3 2 0 3 -4 4 -8 5 8}
    TimeCourse_UpdateGraphData 3 {1 2 2 9 3 -9 4 8 5 10}
    TimeCourse_UpdateGraphData 4 {1 -4 2 8 3 -9 4 8 5 1}
    TimeCourse_UpdateGraphData 6 {1 3 2 8 3 9 4 9 5 -2.3}
    TimeCourse_UpdateGraphData 7 {1 8 2 -1 3 -8 4 5 5 -1}
    TimeCourse_UpdateNumConditions 7
}





