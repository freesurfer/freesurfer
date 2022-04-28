#
# $Id: 14Point.fs,v 1.1.1.1.2.1 2001/11/04 05:24:15 idiscovery Exp $
#
proc tixSetFontset {} {
    global tixOption tcl_platform

    switch -- $tcl_platform(platform) "windows" {
	# This should be Tahoma for Win2000/XP
	set font "MS Sans Serif"
	set fixedfont "Courier New"
    } unix {
	set font "helvetica"
	set fixedfont "courier"
    }

    set tixOption(font)         [list $font -14]
    set tixOption(bold_font)    [list $font -14 bold]
    set tixOption(menu_font)    [list $font -14]
    set tixOption(italic_font)  [list $font -14 bold italic]
    set tixOption(fixed_font)   [list $fixedfont -14]
    set tixOption(border1)      1
}
