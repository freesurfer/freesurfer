#! /usr/bin/tixwish

##
## test.tcl
##
## CVS Revision Info:
##    $Author: nicks $
##    $Date: 2007/01/05 00:20:54 $
##    $Revision: 1.10 $
##
## Copyright (C) 2002-2007,
## The General Hospital Corporation (Boston, MA). 
## All rights reserved.
##
## Distribution, usage and copying of this software is covered under the
## terms found in the License Agreement file named 'COPYING' found in the
## FreeSurfer source code root directory, and duplicated here:
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
##
## General inquiries: freesurfer@nmr.mgh.harvard.edu
## Bug reports: analysis-bugs@nmr.mgh.harvard.edu
##

source "../scripts/fsgdfPlot.tcl"

set fnTestDataDir ""
if { [info exists env(FSDEV_TEST_DATA)] } {
    set fnTestDataDir $env(FSDEV_TEST_DATA)
} else {
    set fnTestDataDir ./test_data
}

set glID {}

FsgdfPlot_Init
set ID [FsgdfPlot_Read "$fnTestDataDir/fsgdf/y-lh.fsgd"]
if { $ID < 0 } { puts "!!!! FsgdfPlot_ParseHeader failed, ID was $ID." ; exit }
FsgdfPlot_ShowWindow $ID

lappend glID $ID

# set ID [FsgdfPlot_Read "$fnTestDataDir/fsgdf/novar.fsgd"]
# if { $ID < 0 } { puts "!!!! FsgdfPlot_ParseHeader failed." ; exit }
# FsgdfPlot_ShowWindow $ID

lappend glID $ID


proc ScheduleRedraw {} {
    global gRedraws glID
    after 500 {

	foreach ID $glID {
	    set vno [expr round(rand() * 10000)]
	    FsgdfPlot_SetPoint $ID $vno 0 0
	    FsgdfPlot_SetInfo $ID "vno $vno"

	    FsgdfPlot_SaveToTable $ID test.table
	    FsgdfPlot_SaveToPostscript $ID test.ps
	}

	incr gRedraws -1
	if { $gRedraws > 0 } [list ScheduleRedraw]
    }
}

set gRedraws 10
ScheduleRedraw


if { 0 } {
after 5000 {
    FsgdfPlot_BeginPointList $gID
    for { set vno 1000 } { $vno < 10000 } { incr vno } {
	FsgdfPlot_AddPoint $gID $vno 0 0
    }
    FsgdfPlot_EndPointList $gID
    FsgdfPlot_SetInfo $gID "vnos 1000 - 10000"

}
}

wm withdraw .
