#! /usr/bin/wish4.0

source $env(MRI_DIR)/lib/tcl/wrappers.tcl

# our graph code
source $env(MRI_DIR)/lib/tcl/tkm_graph.tcl

# unique identifiers for names
set kGraphDataName(0) cond0
set kGraphDataName(1) cond1
set kGraphDataName(2) cond2
set kGraphDataName(3) cond3
set kGraphDataName(4) cond4

# color for conditions
set kGraphDataColor(0) red
set kGraphDataColor(1) blue
set kGraphDataColor(2) green
set kGraphDataColor(3) purple
set kGraphDataColor(4) yellow

# the vmar/hmark stuff likes its colors with braces around them
set kGraphDataColor2(0) {red}
set kGraphDataColor2(1) {blue}
set kGraphDataColor2(2) {green}
set kGraphDataColor2(3) {purple}
set kGraphDataColor2(4) {yellow}

# set the points and errors for a particular condition
proc SetGraphData { inDataSet inPoints } {

    global kGraphDataName kGraphDataColor kGraphDataColor2 gGraph gCanvas 
    global gCondition gTimeSecond gGraphData

    # save our graph data
    set gGraphData($inDataSet) $inPoints

    # set the points. use the color for the this condition.
    gGraph data $kGraphDataName($inDataSet) \
      -colour $kGraphDataColor($inDataSet) \
      -points 1 -lines 1 \
      -coords $inPoints

    # draw some axes
    gGraph vmark 0 yAxis {black}
    gGraph hmark 0 xAxis {black}
    
    # create the time bar at the current second in the color for the
    # current condition
    gGraph vmark $gTimeSecond gCurrentTimePointMark \
      $kGraphDataColor2($gCondition)
}

proc SetGraphErrorBars { inDataSet inErrors } {

    global kGraphDataName kGraphDataColor kGraphDataColor2 gGraph gCanvas 
    global gCondition gTimeSecond gGraphData

    # for each error...
    set theLength [llength $inErrors]
    for {set theIndex 0} {$theIndex < $theLength} {incr theIndex} {

  # get the error amount.
  set theError [lindex $inErrors $theIndex]

  # get the x/y coords for the corresponding point.
  set theGraphX [lindex $gGraphData($inDataSet) [expr $theIndex * 2]]
  set theGraphY [lindex $gGraphData($inDataSet) \
    [expr [expr $theIndex * 2] + 1]]

  # convert it to canvas coords
  set theCanvasX [gGraph x2canvas $theGraphX]
  set theCanvasY [gGraph y2canvas $theGraphY]

  # draw a vertical line at the point, extending above and beyond
  # the point the amounto f the error. use the right color.
  set theTag errorBars$inDataSet-$theIndex
  $gCanvas delete -withtag $theTag
  $gCanvas create line $theCanvasX [expr $theCanvasY - $theError] \
    $theCanvasX [expr $theCanvasY + $theError] \
    -fill $kGraphDataColor($inDataSet) \
    -tag $theTag
    }
}

# called when the user clicks in the canvas.
proc HandleGraphClick { inX inY } {

    global gGraph

    # get the x coord of the graph and send it to the c code.
    set theX [gGraph canvas2x $inX]
    SetCurrentFunctionalTimeSecond $theX
}

# called by the c code when the display flag is updated
proc UpdateDisplayStatus { inIsDisplay } {

    global gIsDisplay

    set gIsDisplay $inIsDisplay
}

# called by the c code when the time resolution is updated
proc UpdateTimeResolution { inTimeResolution } {

    global gTimeResolution
    
    set gTimeResolution $inTimeResolution
}

# called by the c code when the number of prestim points is changed
proc UpdateNumPreStimPoints { inNumPreStimPoints } {

    global gNumPreStimPoints

    set gNumPreStimPoints $inNumPreStimPoints
}

# called by the c code when the time point is updated.
proc UpdateTimePoint { inTimePoint inTimeSecond } {

    global gGraph gTimePoint gTimeSecond

    # set our globals
    set gTimePoint $inTimePoint
    set gTimeSecond $inTimeSecond

    # move the current time bar.
    gGraph movevmark gCurrentTimePointMark $gTimeSecond
}

# called by the c code when the condition is updated
proc UpdateCondition { inCondition } {

    global gCondition gTimeSecond kGraphDataColor2

    # set our global
    set gCondition $inCondition

    # clear the existing time bar and make a new one with the color
    # of the new condition.
    gGraph clearmark gCurrentTimePointMark
    gGraph vmark $gTimeSecond gCurrentTimePointMark \
      $kGraphDataColor2($gCondition)
}

proc UpdatePath { inPath } {

    global gPath

    set gPath $inPath
}

proc UpdateStem { inStem } {

    global gStem

    set gStem $inStem
}

# called by the c code when the threshold min is updated
proc UpdateThresholdMin { inMin } {

    global gThresholdMin

    # set our global
    set gThresholdMin $inMin
}

# called by the c code when the threshold mid is updated
proc UpdateThresholdMid { inMid } {

    global gThresholdMid

    # set our global
    set gThresholdMid $inMid
}

# called by the c code when the threshold max is updated
proc UpdateThresholdMax { inMax } {

    global gThresholdMax

    # set our global
    set gThresholdMax $inMax
}


# our global vars
set gIsDisplay 1
set gTimeResolution 0
set gNumPreStimPoints 0
set gTimePoint 0
set gTimeSecond 0
set gCondition 0
set gPath ""
set gStem ""
set gThresholdMin 0
set gThresholdMid 0
set gThresholdMax 0
set gGraphData(0) {}
set gGraphData(1) {}
set gGraphData(2) {}
set gGraphData(3) {}
set gGraphData(4) {}
set gGraphData(5) {}

# use var names for our widget heirarchy
set theWindow .functionalWindow
 set gCanvas $theWindow.graphCanvas
 set theGraphKeyFrame $theWindow.graphKeyFrame
  set theGraphKeyLabel1 $theGraphKeyFrame.label1
  set theGraphKeyLabel2 $theGraphKeyFrame.label2
  set theGraphKeyLabel3 $theGraphKeyFrame.label3
  set theGraphKeyLabel4 $theGraphKeyFrame.label4
  set theGraphKeyLabel5 $theGraphKeyFrame.label5
  set theGraphKeyLabel6 $theGraphKeyFrame.label6
 set theParametersFrame $theWindow.paramtersFrames
  set theDisplayFrame $theParametersFrame.displayFrame
   set theDisplayCheckbutton $theDisplayFrame.checkbuton
  set theDirectoryFrame $theParametersFrame.directoryFrame
   set thePathLabel $theDirectoryFrame.pathLabel
   set thePathField $theDirectoryFrame.pathField
   set theStemLabel $theDirectoryFrame.stemLabel
   set theStemField $theDirectoryFrame.stemField
  set theTimeConfigFrame $theParametersFrame.timeConfigFrame
   set theTimeResolutionLabel $theTimeConfigFrame.timeResolutionLabel
   set theTimeResolutionField $theTimeConfigFrame.timeResolutionField
   set thePreStimLabel $theTimeConfigFrame.preStimLabel
   set thePreStimField $theTimeConfigFrame.preStimField
  set theTimePointFrame $theParametersFrame.timePointFrame
   set theTimePointLabel $theTimePointFrame.label
   set theTimePointField $theTimePointFrame.field
   set theTimePointIncrement $theTimePointFrame.increment
   set theTimePointDecrement $theTimePointFrame.decrement
  set theConditionFrame $theParametersFrame.conditionFrame
   set theConditionLabel $theConditionFrame.label
   set theConditionField $theConditionFrame.field
   set theConditionIncrement $theConditionFrame.increment
   set theConditionDecrement $theConditionFrame.decrement
  set theThresholdFrame $theParametersFrame.thresholdFrame
   set theThresholdMinLabel $theThresholdFrame.minLabel
   set theThresholdMinField $theThresholdFrame.minField
   set theThresholdMidLabel $theThresholdFrame.midLabel
   set theThresholdMidField $theThresholdFrame.midField
   set theThresholdMaxLabel $theThresholdFrame.maxLabel
   set theThresholdMaxField $theThresholdFrame.maxField

# create top level window and set its name
toplevel $theWindow
wm title $theWindow "Functional Data"

# canvas for the graph
canvas $gCanvas -width 400 -height 300
pack $gCanvas

# create graph
emu_graph gGraph -canvas $gCanvas
bind $gCanvas <ButtonRelease-1> { 
    HandleGraphClick %x %y
 }

# graph key frame is a bunch of labels that display the color of each
# condition. secret: you can click on the label to select that condition.
frame $theGraphKeyFrame
pack $theGraphKeyFrame -anchor w

label $theGraphKeyLabel1 -text "Graph key:"
pack $theGraphKeyLabel1 -side left

label $theGraphKeyLabel2 -text "Cond 0" -foreground $kGraphDataColor(0)
bind $theGraphKeyLabel2 <ButtonRelease-1> {
    SetCurrentFunctionalCondition 0
}
pack $theGraphKeyLabel2 -side left

label $theGraphKeyLabel3 -text "Cond 1" -foreground $kGraphDataColor(1)
bind $theGraphKeyLabel3 <ButtonRelease-1> {
    SetCurrentFunctionalCondition 1
}
pack $theGraphKeyLabel3 -side left

label $theGraphKeyLabel4 -text "Cond 2" -foreground $kGraphDataColor(2)
bind $theGraphKeyLabel4 <ButtonRelease-1> {
    SetCurrentFunctionalCondition 2
}
pack $theGraphKeyLabel4 -side left

label $theGraphKeyLabel5 -text "Cond 3" -foreground $kGraphDataColor(3)
bind $theGraphKeyLabel5 <ButtonRelease-1> {
    SetCurrentFunctionalCondition 3
}
pack $theGraphKeyLabel5 -side left

label $theGraphKeyLabel6 -text "Cond 4" -foreground $kGraphDataColor(4)
bind $theGraphKeyLabel6 <ButtonRelease-1> {
    SetCurrentFunctionalCondition 4
}
pack $theGraphKeyLabel6 -side left

# parameter frame is a frame for all parameters.
frame $theParametersFrame
pack $theParametersFrame -anchor w

# checkbutton for display flag
frame $theDisplayFrame
pack $theDisplayFrame -anchor w

checkbutton $theDisplayCheckbutton -variable gIsDisplay -text "Display Overlay"
bind $theDisplayCheckbutton <ButtonRelease-1> {
    global gIsDisplay
    SetCurrentFunctionalDisplayStatus $gIsDisplay
}
pack $theDisplayCheckbutton -side left

# field for the path
frame $theDirectoryFrame
pack $theDirectoryFrame -fill x -expand true

label $thePathLabel -text "Directory:"
pack $thePathLabel -side left

entry $thePathField -textvariable gPath
pack $thePathField -side left -fill x -expand true

# field for stem
label $theStemLabel -text "Stem:"
pack $theStemLabel -side left

entry $theStemField -textvariable gStem -width 3
pack $theStemField -side left

# use to change the time res and number of pre stim values to affect
# the display of the xaxis in the graph.
frame $theTimeConfigFrame
pack $theTimeConfigFrame -anchor w

label $theTimeResolutionLabel -text "Time resolution (TR):"
pack $theTimeResolutionLabel -side left

entry $theTimeResolutionField -width 2 -textvariable gTimeResolution
bind $theTimeResolutionField <Return> {
    SetCurrentFunctionalTimeResolution $gTimeResolution
}
pack $theTimeResolutionField -side left

label $thePreStimLabel -text "Number of pre-stim points:"
pack $thePreStimLabel -side left

entry $thePreStimField -width 2 -textvariable gNumPreStimPoints
bind $thePreStimField <Return> {
    SetCurrentFunctionalNumPreStimPoints $gNumPreStimPoints
}
pack $thePreStimField -side left

# use to change the current time point. must press return to activate the
# change. the new value is sent to the c code. the c code then calls
# the script's UpdateXXX function which sets the global variable with the
# new value. that way, the c code can override the value if there is an error.
frame $theTimePointFrame
pack $theTimePointFrame -anchor w

label $theTimePointLabel -text "Current time point:"
pack $theTimePointLabel -side left

entry $theTimePointField -width 2 -textvariable gTimePoint
bind $theTimePointField <Return> {
    SetCurrentFunctionalTimePoint $gTimePoint
}
pack $theTimePointField -side left

# handy little plus/minus buttons for incrementing and decrementing the
# variable.
button $theTimePointDecrement -text -
bind $theTimePointDecrement <ButtonRelease-1> {
    global gTimePoint
    SetCurrentFunctionalTimePoint [expr $gTimePoint - 1]
}
pack $theTimePointDecrement -side left -padx 5

button $theTimePointIncrement -text +
bind $theTimePointIncrement <ButtonRelease-1> {
    global gTimePoint
    SetCurrentFunctionalTimePoint [expr $gTimePoint + 1]
}
pack $theTimePointIncrement -side left -padx 5

# entry field for the condition.
frame $theConditionFrame
pack $theConditionFrame -anchor w

label $theConditionLabel -text "Current condition:"
pack $theConditionLabel -side left

entry $theConditionField -width 2 -textvariable gCondition
bind $theConditionField <Return> {
    SetCurrentFunctionalCondition $gCondition
}
pack $theConditionField -side left

# dec and inc buttons for the condition
button $theConditionDecrement -text -
bind $theConditionDecrement <ButtonRelease-1> {
    global gCondition
    SetCurrentFunctionalCondition [expr $gCondition - 1]
}
pack $theConditionDecrement -side left -padx 5

button $theConditionIncrement -text +
bind $theConditionIncrement <ButtonRelease-1> {
    global gCondition
    SetCurrentFunctionalCondition [expr $gCondition + 1]
}
pack $theConditionIncrement -side left -padx 5

# entry field for the threshold min.
frame $theThresholdFrame
pack $theThresholdFrame -anchor w

label $theThresholdMinLabel -text "Overlay threshold: min:"
pack $theThresholdMinLabel -side left

entry $theThresholdMinField -width 5 -textvariable gThresholdMin
bind $theThresholdMinField <Return> {
    SetCurrentFunctionalThresholdMin $gThresholdMin
}
pack $theThresholdMinField -side left

label $theThresholdMidLabel -text "mid:"
pack $theThresholdMidLabel -side left

# entry field for the threshold mid.
entry $theThresholdMidField -width 5 -textvariable gThresholdMid
bind $theThresholdMidField <Return> {
    SetCurrentFunctionalThresholdMid $gThresholdMid
}
pack $theThresholdMidField -side left

label $theThresholdMaxLabel -text "max:"
pack $theThresholdMaxLabel -side left

# entry field for the threshold max.
entry $theThresholdMaxField -width 5 -textvariable gThresholdMax
bind $theThresholdMaxField <Return> {
    SetCurrentFunctionalThresholdMax $gThresholdMax
}
pack $theThresholdMaxField -side left



# draw some inital data
#set theDummyData {0 0 1 1}
#set theDummyErrors {0 0}
#SetGraphData 0 $theDummyData $theDummyError

#set theTestData0 {-6 3 -4 -2.2 -2 0.333 0 4.001 2 6 4 -3.2 6 4.67 8 10.40 \
#    10 14.9901 12 12.1 14 24.1 16 9.2 18 2.09 20 -9.19 22 10.99}
#set theTestData1 {-6 -5 -4 34.2 -2 12.333 0 12.001 2 2 4 10 6 34.67 8 90.40 \
#    10 47.9901 12 89.1 14 90.1 16 30.2 18 35.09 20 14.19 22 13.99}
#set theTestData2 {-6 -4 -4 12.2 -2 10.333 0 45.001 2 4 4 20 6 25.67 8 13.40 \
#    10 24.9901 12 40.1 14 45.1 16 43.2 18 39.09 20 14.19 22 13.99}
#set theErrors {10 6.5 8.9 12 9.4 7.9 9.1 5.3 7.9 4.5 13.2 16.4 9.1 8.9 9.55}

#SetGraphData 0 $theTestData0 $theErrors
#SetGraphData 1 $theTestData1 $theErrors
#SetGraphData 2 $theTestData2 $theErrors
