#
# $Id: Old12Pt.fs,v 1.1.1.1 2000/05/17 11:08:48 idiscovery Exp $
#

proc tixSetFontset:12Point {} {
    global tixOption

    set tixOption(font)         -*-helvetica-medium-r-normal-*-12-*-*-*-*-*-*-*
    set tixOption(bold_font)    -*-helvetica-bold-r-normal-*-12-*-*-*-*-*-*-*
    set tixOption(menu_font)    -*-helvetica-bold-o-normal-*-12-*-*-*-*-*-*-*
    set tixOption(italic_font)  -*-helvetica-bold-o-normal-*-12-*-*-*-*-*-*-*
    set tixOption(fixed_font) -*-courier-medium-r-*-*-12-*-*-*-*-*-*-*
    set tixOption(border1)      1
}
