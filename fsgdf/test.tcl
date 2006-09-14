#! /usr/bin/tixwish
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
