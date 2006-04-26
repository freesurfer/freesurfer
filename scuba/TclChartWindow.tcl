#! /usr/pubsw/bin/tixwish

# $Id: TclChartWindow.tcl,v 1.4 2006/04/26 17:57:11 kteich Exp $

package require Tix;
package require BLT;

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


# Also look for tkUtils.tcl.
set sDefaultScriptsDir ""
catch { set sDefaultScriptsDir "$env(FREESURFER_HOME)/lib/tcl" }
set sUtilsDir ""
catch { set sUtilsDir "$env(TKUTILS_SCRIPTS_DIR)" }

set fnUtils \
    [FindFile tkUtils.tcl \
	 [list $sUtilsDir "." "../scripts" $sDefaultScriptsDir]]
if { [string compare $fnUtils ""] == 0 } { exit }
source $fnUtils


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
#     lPoints

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
    if { ![info exists gData($iID,lPoints)] } {
	set gData($iID,lPoints) {}
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

    # Place everythingin the window.
    grid $gwChart       -column 0 -row 0 -columnspan 3 -sticky news
    grid $lwInfo        -column 0 -row 1 -sticky nwe
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
    
    # If we have no points, return.
    if { ![info exists gData($iID,lPoints)] || 
	 [llength $gData($iID,lPoints)] == 0 } {
	return
    }

    # Create an element for each point.
    for { set nPoint 0 } { $nPoint < [llength $gData($iID,lPoints)] } { incr nPoint } {
	
	array set point [lindex $gData($iID,lPoints) $nPoint]
	
	$gw element create point$nPoint \
	    -data [list $point(x) $point(y)] \
	    -label $point(label) \
	    -linewidth 0 -outlinewidth 1 \
	    -activepen activeElement
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

# This function expects data to come in as a list of array lists. Each
# element in ilPoints should be a list of array label/value pairs, e.g.
#  [list x 10.4 y 20.98 label "Hello World"]
proc Chart_SetPointData { iID ilPoints } {
    global gData
    if { [lsearch $gData(lID) $iID] == -1 } { puts "ID not found"; return }
    set gData($iID,lPoints) $ilPoints
    Chart_PlotData $iID
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
