#! /usr/pubsw/bin/tixwish

# $Id: TclChartWindow.tcl,v 1.10 2006/07/13 17:46:18 kteich Exp $

package require Tix;
package require BLT;


# Look for tkUtils.tcl.
set sDefaultScriptsDir ""
catch { set sDefaultScriptsDir "$env(FREESURFER_HOME)/lib/tcl" }
set sUtilsDir ""
catch { set sUtilsDir "$env(TKUTILS_SCRIPTS_DIR)" }

set bLoaded 0
foreach sPath [list $sUtilsDir "." "../scripts" $sDefaultScriptsDir] {
    if { !$bLoaded } {
	set sFullFileName [ file join $sPath tkUtils.tcl ]
	if { [file readable $sFullFileName] } {
	    puts "Reading $sFullFileName"
	    source $sFullFileName
	    set bLoaded 1
	}
    }
}
if { !$bLoaded } { exit }


# constant values for stuff
set kValid(lMarkers) {square circle diamond plus cross splus scross triangle}
set kValid(lColors) {red blue green yellow black purple orange pink brown}

# gData
#   lID
#   $ID
#     title
#     xAxisLabel
#     yAxisLabel
#     showLegend
#     lGroups
#     $nGroup
#       lPoints - list of arrays, fields are (x, y, label)
#       connected
#       label
#       red
#       green
#       blue
#     lXAxisMarkers - list of markers, fields are (value, label, red, green, blue)

# gChart
#   $ID
#     state
#       info
#       pointsChanged
#       hiElement
#       window
#         geometry

# gWidgets
#   $ID
#     wwTop
#     gwChart
#     lwInfo
#     windowBuilt

# Builds the main window. 
proc Chart_BuildWindow { iID } {
    global gWidgets gData gChart

    set wwTop         .chart-$iID
    set gwChart        $wwTop.gwChart
    set lwInfo         $wwTop.lwInfo
    set fwReport       $wwTop.fwReport
    set bwReport       $fwReport.bwReport
    set bwPicture      $fwReport.bwPicture

    # Fill out default values.
    if { ![info exists gData($iID,title)] } {
	set gData($iID,title) ""
    }
    if { ![info exists gData($iID,xAxisLabel)] } {
	set gData($iID,xAxisLabel) ""
    }
    if { ![info exists gData($iID,yAxisLabel)] } {
	set gData($iID,yAxisLabel) ""
    }
    if { ![info exists gData($iID,showLegend)] } {
	set gData($iID,showLegend) false
    }
    if { ![info exists gData($iID,lGroups)] } {
	set gData($iID,lGroups) {}
    }

    # Make the to window and set its title.
    toplevel $wwTop -height 500 -width 500
    wm title $wwTop $gData($iID,title)

    # Make the graph.
    blt::graph $gwChart \
	-title $gData($iID,title) \
	-plotbackground white \
	-relief raised -border 2

    # Bind our callbacks.
    $gwChart legend bind all <Enter> [list Chart_CBLegendEnter $iID %W]
    $gwChart legend bind all <Leave> [list Chart_CBLegendLeave $iID %W]
    $gwChart legend bind all <ButtonPress-1> [list Chart_CBLegendClick $iID %W]
    bind $gwChart <Motion> [list Chart_CBGraphMotion $iID %W %x %y]
    bind $gwChart <Destroy> [list Chart_CBCloseWindow $iID] 

    # Set the y axis label to the measurement name.
    $gwChart axis configure y -title $gData($iID,yAxisLabel)

    # Make the info label.
    set gChart($iID,state,info) ""
    tkuMakeActiveLabel $lwInfo \
	-variable gChart($iID,state,info)

    # Make the report button.
    frame $fwReport
    button $bwReport \
	-text "Generate Report" \
	-command "Chart_DoGenerateReportDlog $iID"
    button $bwPicture \
	-text "Save Picture" \
	-command "Chart_SaveAsPostscript $iID"
    pack $bwReport $bwPicture \
	-side right

    # Place everything in the window.
    grid $gwChart       -column 0 -row 0 -columnspan 3 -sticky news
    grid $lwInfo        -column 0 -row 1 -sticky nwe
    grid $fwReport      -column 0 -row 2 -sticky nwe
    grid columnconfigure $wwTop 0 -weight 1
    grid rowconfigure $wwTop 0 -weight 1
    grid rowconfigure $wwTop 1 -weight 0

    # Set the names in the gWidgets array.
    set gWidgets($iID,wwTop)          $wwTop
    set gWidgets($iID,gwChart)        $gwChart
    set gWidgets($iID,lwInfo)         $lwInfo

    # Create the pen for our active element.
    $gwChart pen create activeElement \
	-symbol circle -color red -pixels 0.2i -fill ""

    # Note that the window has been built.
    set gWidgets($iID,windowBuilt) 1

    # Insert the ID.
    lappend gData(lID) $iID
}

# This plots the current data on the graph. It is fast enough that it
# can be called any time the data is changed to completely redraw it
# from scratch.
proc Chart_PlotData { iID } {
    global gWidgets gChart gData

    # Don't plot if the window isn't built or we don't have data.
    if { ![info exists gWidgets($iID,windowBuilt)] ||
	 !$gWidgets($iID,windowBuilt) } {
	return
    }

    set gw $gWidgets($iID,gwChart)

    # Update the window title.
    wm title $gWidgets($iID,wwTop) $gData($iID,title)

    # Set the graph title.
    $gw configure -title $gData($iID,title)

    # Set the axis titles.
    $gw axis configure x -title $gData($iID,xAxisLabel)
    $gw axis configure y -title $gData($iID,yAxisLabel)

    # Remove all the elements and markers from the graph.
    set lElements [$gw element names *]
    foreach element $lElements {
	$gw element delete $element
    }
    set lMarkers [$gw marker names *]
    foreach marker $lMarkers {
	$gw marker delete $marker
    }
    
    # If we have no group, return.
    if { ![info exists gData($iID,lGroups)] || 
	 [llength $gData($iID,lGroups)] == 0 } {
	return
    }

    # For each group...
    foreach nGroup $gData($iID,lGroups) {
	
	set color blue
	if { [info exists gData($iID,$nGroup,red)] } {
	    set r $gData($iID,$nGroup,red)
	    set g $gData($iID,$nGroup,green)
	    set b $gData($iID,$nGroup,blue)
	    set color [format "#%.2x%.2x%.2x" $r $g $b]
	}
	
	# Set options based on our draw data.
	if { [info exists gData($iID,$nGroup,connected)] &&
	     $gData($iID,$nGroup,connected) } {
	    
	    set lPoints {}
	    for { set nPoint 0 } \
		{ $nPoint < [llength $gData($iID,$nGroup,lPoints)] } \
		{ incr nPoint } {
		    
		    array set point \
			[lindex $gData($iID,$nGroup,lPoints) $nPoint]

		    lappend lPoints $point(x) $point(y)
		}    

	    set sLabel ""
	    if { [info exists gData($iID,$nGroup,label)] } {
		set sLabel $gData($iID,$nGroup,label)
	    }

	    $gw element create group$nGroup \
		-data $lPoints \
		-label $sLabel \
		-color $color \
		-linewidth 1 -outlinewidth 1 \
		-activepen activeElement

	} else {
	    
	    # Create an element for each point.
	    for { set nPoint 0 } \
		{ $nPoint < [llength $gData($iID,$nGroup,lPoints)] } \
		{ incr nPoint } {
		    
		    array set point \
			[lindex $gData($iID,$nGroup,lPoints) $nPoint]
		    
		    $gw element create group$nGroup-point$nPoint \
			-data [list $point(x) $point(y)] \
			-label $point(label) \
			-color $color \
			-linewidth 0 -outlinewidth 1 \
			-activepen activeElement
		}
	}
    }
    
    # Draw the markers.
    if { [info exists gData($iID,lXAxisMarkers)] } {
	for { set nMarker 0 } \
	    { $nMarker < [llength $gData($iID,lXAxisMarkers)] } \
	    { incr nMarker } {

		array set xMarker \
		    [lindex  $gData($iID,lXAxisMarkers) $nMarker]
		
		$gw marker create line \
		    -coords [list $xMarker(value) -Inf $xMarker(value) Inf] \
		    -name marker$nMarker \
		    -outline [format "#%.2x%.2x%.2x" $xMarker(red) \
				$xMarker(green) $xMarker(blue)]

		$gw marker create text \
		    -text $xMarker(label) \
		    -background "" -fill "" \
		    -rotate 270 \
		    -anchor nw \
		    -coords [list $xMarker(value) Inf] \
		    -name textMarker$nMarker \
		    -outline [format "#%.2x%.2x%.2x" $xMarker(red) \
				$xMarker(green) $xMarker(blue)]
	    }
    }
    
    # Show or hide the legend.
    if { $gData($iID,showLegend) } {
	$gw legend configure -hide false
    } else {
 	$gw legend configure -hide true
    }
    
    set gChart($iID,state,pointsChanged) 0
}

# Hilight/UnhilightElement works on an element by name (which could be
# a subject or class, depending on viewing mode). It will
# select/unselect the element name in the legend and change the
# drawing pen of the element in the graph, which if activated draws it
# with a red circle around it.
proc Chart_HilightElement { iID iElement } {
    global gWidgets
    $gWidgets($iID,gwChart) legend activate $iElement
    $gWidgets($iID,gwChart) element activate $iElement
}

proc Chart_UnhilightElement { iID iElement } {
    global gWidgets
    $gWidgets($iID,gwChart) legend deactivate $iElement
    $gWidgets($iID,gwChart) element deactivate $iElement
}

# Focus/Unfocus is called to 'mouseover' an element. It
# Hilight/Unhilights an element and puts or removes the subject name
# in a text marker in the graph.
proc Chart_UnfocusElement { iID } {
    global gChart gWidgets
    
    # If we have a focused element, unhighlight it, set the
    # highlighted element name to null, and delete the hover text
    # marker.
    if { [info exists gChart($iID,state,hiElement)] && \
 	     "$gChart($iID,state,hiElement)" != "" } {
 	Chart_UnhilightElement $iID $gChart($iID,state,hiElement)
 	set gChart($iID,state,hiElement) ""
 	$gWidgets($iID,gwChart) marker delete hover
    }
}

proc Chart_FocusElement { iID iElement inSubjInClass iX iY } {
    global gChart gWidgets gData
    
    # Don't focus on error bars.
    if { [string match error* $iElement] } {
 	return
    }
    
    # Set the highlighted element name and highlight the element.
    set gChart($iID,state,hiElement) $iElement
    Chart_HilightElement $iID $gChart($iID,state,hiElement)
    
    $gWidgets($iID,gwChart) marker create text \
 	-name hover \
	-text [$gWidgets($iID,gwChart) element cget $iElement -label] \
	-anchor nw \
 	-coords [list $iX $iY]
}


# Finds the element under the mouse.
proc Chart_FindMousedElement { iID iX iY } {
    global gWidgets
    set bFound [$gWidgets($iID,gwChart) element closest $iX $iY aFound -halo 10]
    if { $bFound } {
 	return [list $aFound(name) $aFound(index) $aFound(x) $aFound(y)]
    }
    return ""
}



# Our callbacks.
proc Chart_CBCloseWindow { iID } {
    global gWidgets
    set gWidgets($iID,windowBuilt) 0

    # Tells the C code to delete the associated C object.
    DeleteTclChartWindow $iID
}

proc Chart_CBLegendEnter { iID igw } {
    Chart_HilightElement $iID [$igw legend get current]
}

proc Chart_CBLegendLeave { iID igw } {
    Chart_UnhilightElement $iID [$igw legend get current]
}

proc Chart_CBLegendClick { iID igw } {
}

proc Chart_CBGraphMotion { iID igw iX iY } {
    Chart_UnfocusElement $iID
    set lResult [Chart_FindMousedElement $iID $iX $iY]
    set element [lindex $lResult 0]
    if { "$element" != "" } { 
 	set index [lindex $lResult 1]
 	set x [lindex $lResult 2]
 	set y [lindex $lResult 3]
 	Chart_FocusElement $iID $element $index $x $y
    }
}

proc Chart_DoGenerateReportDlog { iID } {

    tkuDoFileDlog -title "Generate Chart Report" \
	-prompt1 "Include group label column" \
	-defaultvalue1 1 \
	-type1 checkbox \
	\
	-prompt2 "Include point label column" \
	-defaultvalue2 1 \
	-type2 checkbox \
	\
	-prompt3 "Include X value column" \
	-defaultvalue3 1 \
	-type3 checkbox \
	\
	-prompt4 "Include Y value column" \
	-defaultvalue4 1 \
	-type4 checkbox \
	\
	-prompt5 "Save file: " \
	-okCmd "GenerateChartReport $iID %s1 %s2 %s3 %s4 %s5"
}

proc Chart_SaveAsPostscript { iID } {
    global gWidgets

    tkuDoFileDlog -title "Save Chart Picture" \
	-prompt1 "Save Postscript Picture As: " \
	-okCmd "$gWidgets($iID,gwChart) postscript output %s1"

}

# ============================================================ PUBLIC


# Call once before anything else to initialize the data structures.
proc Chart_Init {} {
    global gWidgets gbLibLoaded gData
    set gData(lID) {}
}

# Create the window with this ID. This should be called once per
# window.
proc Chart_NewWindow { iID } {
    global gData gWidgets
    if { ![info exists gWidgets($iID,windowBuilt)] ||
 	 !$gWidgets($iID,windowBuilt) } {
 	Chart_BuildWindow $iID
    }
}

# Close the window.
proc Chart_CloseWindow { iID } {
    global gData gWidgets
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }
    wm withdraw $gWidgets($iID,wwTop)
}

# Show or hide the window.
proc Chart_ShowWindow { iID } {
    global gData gWidgets
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }
    wm deiconify $gWidgets($iID,wwTop)
    if { [info exists gWidgets($iID,state,window,geometry)] } {
 	wm geometry $gWidgets($iID,wwTop) $gWidgets($iID,state,window,geometry)
    }
    Chart_PlotData $iID
}

proc Chart_HideWindow { iID } {
    global gData gWidgets
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }
    if { [info exists gWidgets($iID,wwTop)] } {
 	set gWidgets($iID,state,window,geometry) \
 	    [wm geometry $gWidgets($iID,wwTop)]
 	wm inconify $gWidgets($iID,wwTop)
    }
}

# Set the window title.
proc Chart_SetWindowTitle { iID isTitle } {
    global gData
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }
    set gData($iID,title) $isTitle
    Chart_PlotData $iID
}

# Set the label under the x axis.
proc Chart_SetXAxisLabel { iID isLabel } {
    global gData
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }
    set gData($iID,xAxisLabel) $isLabel
    Chart_PlotData $iID
}

# Set the label to the side of the y axis.
proc Chart_SetYAxisLabel { iID isLabel } {
    global gData
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }
    set gData($iID,yAxisLabel) $isLabel
    Chart_PlotData $iID
}

# Show or hide the legend.
proc Chart_SetShowLegend { iID ibShowLegend } {
    global gData
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }
    set gData($iID,showLegend) $ibShowLegend
    Chart_PlotData $iID
}

# Set the info string displayed under the graph.
proc Chart_SetInfo { iID isInfo } {
    global gData gChart 
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }
    set gChart($iID,state,info) $isInfo
}

proc Chart_ClearData { iID } {
    global gData
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }

    foreach nGroup $gData($iID,lGroups) {
	set gData($iID,$nGroup,lPoints) {}
    }

    set gData($iID,lGroups) {}
    set gData($iID,lXAxisMarkers) {}
}

# This function expects data to come in as a list of array lists. Each
# element in ilPoints should be a list of array label/value pairs, e.g.
#  [list x 10.4 y 20.98 label "Hello World"]
proc Chart_SetPointData { iID inGroup ilPoints } {
    global gData
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }

    # Add this group number to the list of groups if it's not already
    # there.
    if { ![info exists gData($iID,lGroups)] ||
	 [lsearch $gData($iID,lGroups) $inGroup] == -1 } { 
	lappend gData($iID,lGroups) $inGroup
    }

    # Set the point data.
    set gData($iID,$inGroup,lPoints) $ilPoints

    Chart_PlotData $iID
}

# This function expects data to come in as a list of array lists. Each
# element in ilMarkers should be a list of array label/value pairs,
# e.g. [list value 5.0 label "Hello" red 255 green 0 blue 128]
proc Chart_SetXAxisMarkers { iID ilMarkers } {
    global gData
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }

    # Set the marker data.
    set gData($iID,lXAxisMarkers) $ilMarkers

    Chart_PlotData $iID
}

# Set display information about the groups.
proc Chart_SetGroupConnected { iID inGroup ibConnected } {
    global gData
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }

    # Set the value.
    set gData($iID,$inGroup,connected) $ibConnected
}

proc Chart_SetGroupLabel { iID inGroup isLabel } {
    global gData
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }

    # Set the value.
    set gData($iID,$inGroup,label) $isLabel
}

proc Chart_SetGroupColor { iID inGroup iRed iGreen iBlue } {
    global gData
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }

    # Set the values.
    set gData($iID,$inGroup,red) $iRed
    set gData($iID,$inGroup,green) $iGreen
    set gData($iID,$inGroup,blue) $iBlue
}

# Save the current plot graphic to a postscript file.
proc Chart_SaveToPostscript { iID ifnPS } {
    global gData gWidgets 
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }
    set err [catch {$gWidgets($iID,gwChart) postscript output $ifnPS} sResult]
    if { $err } {
	puts "Could not save postscript file: $sResult"
    }
}
