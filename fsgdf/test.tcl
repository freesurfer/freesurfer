#! /usr/bin/tixwish
source "../scripts/fsgdfPlot.tcl"

set fnTestDataDir ""
if { [info exists env(FSDEV_TEST_DATA)] } {
    set fnTestDataDir $env(FSDEV_TEST_DATA)
} else {
    set fnTestDataDir /space/lyon/1/fsdev/test_data
}

FsgdfPlot_Init
set err [FsgdfPlot_ParseHeader "$fnTestDataDir/fsgdf/y-lh.fsgd"]
if { $err } { puts "!!!! FsgdfPlot_ParseHeader failed." ; exit }
FsgdfPlot_ShowWindow

#FsgdfPlot_SaveToTable test.table
#FsgdfPlot_SaveToPostscript test.ps


proc ScheduleRedraw {} {
    global gRedraws
    after 500 {
	set vno [expr round(rand() * 10000)]
	FsgdfPlot_SetPoint $vno 0 0
	FsgdfPlot_SetInfo "vno $vno"
	incr gRedraws -1
	if { $gRedraws > 0 } { ScheduleRedraw }
    }
}

set gRedraws 10
ScheduleRedraw




if { 0 } {
after 5000 {
    FsgdfPlot_BeginPointList
    for { set vno 1000 } { $vno < 10000 } { incr vno } {
	FsgdfPlot_AddPoint $vno 0 0
    }
    FsgdfPlot_EndPointList
    FsgdfPlot_SetInfo "vnos 1000 - 10000"

}
}

wm withdraw .