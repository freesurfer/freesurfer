# -*-mode: tcl; fill-column: 75; tab-width: 8; coding: iso-latin-1-unix -*-
#
#	$Id: MkManag.tcl,v 1.2.2.1 2001/11/04 05:10:08 idiscovery Exp $
#
# MkManag.tcl --
#
#	This file implements the "Manager" page in the widget demo
#
#	This file has not been properly documented. It is NOT intended
#	to be used as an introductory demo program about Tix
#	programming. For such demos, please see the files in the
#	demos/samples directory or go to the "Samples" page in the
#	"widget demo"
#
#
# Copyright (c) 1996, Expert Interface Technologies
#
# See the file "license.terms" for information on usage and redistribution
# of this file, and for a DISCLAIMER OF ALL WARRANTIES.
#

proc MkManager {nb page} {
    set w [$nb subwidget $page]

    set name [tixOptionName $w]
    option add *$name*TixLabelFrame*label.padX 4

    tixLabelFrame $w.pane -label "tixPanedWindow"
    tixLabelFrame $w.note -label "tixNoteBook"

    MkPanedWindow [$w.pane subwidget frame]
    MkNoteBook    [$w.note subwidget frame]

    tixForm $w.pane -top 0 -left 0   -right $w.note -bottom -1
    tixForm $w.note -top 0  -right -1 -bottom -1
}

proc MkPanedWindow {w} {
    set name [tixOptionName $w]

    message $w.msg \
	-relief flat -width 240 -anchor n\
	-text {The PanedWindow widget allows the user to interactively\
manipulate the\
sizes of several panes. The panes can be arranged either vertically or\
horizontally.}

    label $w.group -text "Newsgroup: comp.lang.tcl"

    tixPanedWindow $w.pane

    set p1 [$w.pane add list -min 70 -size 100]
    set p2 [$w.pane add text -min 70]

    tixScrolledListBox $p1.list
    $p1.list subwidget listbox config -font [tix option get fixed_font]

    tixScrolledText    $p2.text
    $p2.text subwidget text    config -font [tix option get fixed_font]

    $p1.list subwidget listbox insert end \
	"  12324 Re: TK is good for your health" \
	"+ 12325 Re: TK is good for your health" \
	"+ 12326 Re: Tix is even better for your health (Was: TK is good...)" \
	"  12327 Re: Tix is even better for your health (Was: TK is good...)" \
	"+ 12328 Re: Tix is even better for your health (Was: TK is good...)" \
	"  12329 Re: Tix is even better for your health (Was: TK is good...)" \
	"+ 12330 Re: Tix is even better for your health (Was: TK is good...)"

    $p2.text subwidget text config -wrap none -bg \
	[$p1.list subwidget listbox cget -bg]
    $p2.text subwidget text insert end {
Mon, 19 Jun 1995 11:39:52        comp.lang.tcl              Thread   34 of  220
Lines 353       A new way to put text and bitmaps together iNo responses
ioi@xpi.com                  Ioi K. Lam at Expert Interface Technologies

Hi,

I have implemented a new image type called "compound". It allows you
to glue together a bunch of bitmaps, images and text strings together
to form a bigger image. Then you can use this image with widgets that
support the -image option. This way you can display very fancy stuffs
in your GUI. For example, you can display a text string string
together with a bitmap, at the same time, inside a TK button widget.


You can also you is in other places such as putting fancy bitmap+text
in menus, tabs of tixNoteBook widgets, etc. This feature will be
included in the next release of Tix (4.0b1). Count on it to make jazzy
interfaces!}

    pack $p1.list -expand yes -fill both -padx 4 -pady 6
    pack $p2.text -expand yes -fill both -padx 4 -pady 6

    pack $w.msg   -side top -padx 3 -pady 3 -fill both
    pack $w.group -side top -padx 3 -pady 3 -fill both
    pack $w.pane  -side top -padx 3 -pady 3 -expand yes -fill both
}

proc MkNoteBook  {w} {

    message $w.msg \
	-relief flat -width 240 -anchor n\
	-text {The NoteBook widget allows you to lay out a complex\
interface into individual pages.}

    # We use these options to set the sizes of the subwidgets inside the
    # notebook, so that they are well-aligned on the screen.
    #
    set name [tixOptionName $w]
    option add *$name*TixControl*entry.width 10
    option add *$name*TixControl*label.width 18
    option add *$name*TixControl*label.anchor e
    option add *$name*TixNoteBook*tagPadX 8


    tixNoteBook $w.nb -ipadx 6 -ipady 6

    # Create the two tabs on the notebook. The -underline option
    # puts a underline on the first character of the labels of the tabs.
    # Keyboard accelerators will be defined automatically according
    # to the underlined character.	
    #
    $w.nb add hard_disk -label "Hard Disk" -underline 8
    $w.nb add network   -label "Network"   -underline 0
   
    # Create the first page
    #
    set f [$w.nb subwidget hard_disk]

    # the frame for the buttons that are present in all the pages
    #
    frame $f.common
    pack $f.common -side right -padx 2 -pady 2 -fill y
    CreateCommonButtons $w $f.common


    # Create the controls that only belong to this page
    #
    tixControl $f.a -value 12   -label "Access Time: "
    tixControl $f.w -value 400  -label "Write Throughput: "
    tixControl $f.r -value 400  -label "Read Throughput: "
    tixControl $f.c -value 1021 -label "Capacity: "
    pack $f.a $f.w $f.r $f.c  -side top -padx 20 -pady 2
    
    # Create the second page	
    #
    set f [$w.nb subwidget network]

    # the frame for the buttons that are present in all the pages
    #
    frame $f.common
    pack $f.common -side right -padx 2 -pady 2 -fill y

    tixControl $f.a -value 12   -label "Access Time: "
    tixControl $f.w -value 400  -label "Write Throughput: "
    tixControl $f.r -value 400  -label "Read Throughput: "
    tixControl $f.c -value 1021 -label "Capacity: "
    tixControl $f.u -value 10   -label "Users: "
    
    CreateCommonButtons $w $f.common

    pack $f.a $f.w $f.r $f.c $f.u -side top -padx 20 -pady 2
    pack $w.msg  -side top -padx 3 -pady 3 -fill both
    pack $w.nb   -expand yes -fill both -padx 5 -pady 5 -side top
 }

proc CreateCommonButtons {w f} {
    button $f.ok     -text OK     -width 6
    button $f.cancel -text Cancel -width 6

    pack $f.ok $f.cancel -side top -padx 2 -pady 2
}
