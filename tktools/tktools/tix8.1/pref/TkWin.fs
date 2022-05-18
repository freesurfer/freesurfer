#
# $Id: TkWin.fs,v 1.1 2000/10/12 01:41:26 idiscovery Exp $
#

proc tixSetFontset {} {

    global tixOption

    set tixOption(font)         "windows-message"
    set tixOption(bold_font)    "windows-status"
    set tixOption(menu_font)    "windows-menu"
    set tixOption(italic_font)  "windows-message"
    set tixOption(fixed_font)   "systemfixed"
    set tixOption(border1)      1
}

