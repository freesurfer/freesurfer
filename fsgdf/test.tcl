#! /usr/bin/tixwish
source "fsgdfPlot.tcl"

FsgdfPlot_Init
FsgdfPlot_ParseHeader "/home/kteich/fsgdf_data/y-lh.fsgd"
FsgdfPlot_ShowWindow
FsgdfPlot_SetPoint 10000 0 0
FsgdfPlot_SetInfo "vno 10000"

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