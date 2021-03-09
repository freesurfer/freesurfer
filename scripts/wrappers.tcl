#! /usr/bin/tclsh

##
## wrappers.tcl
##
##
## Copyright (c) 1996-1999 Martin Sereno and Anders Dale
## Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
##
## Terms and conditions for use, reproduction, distribution and contribution
## are found in the 'FreeSurfer Software License Agreement' contained
## in the file 'LICENSE' found in the FreeSurfer distribution, and here:
##
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
##
## Reporting: freesurfer@nmr.mgh.harvard.edu
##
############################################################################
############################################################################
### standard widget wrappers

if ![info exists sfont] { set sfont \
  -b&h-lucidatypewriter-medium-r-normal-sans-10-100-75-75-m-60-iso8859-1 }
if ![info exists ffont] { set ffont \
  -b&h-lucidatypewriter-medium-r-normal-sans-12-120-75-75-m-70-iso8859-1 }
if ![info exists ffontb] { set ffontb \
  -b&h-lucidatypewriter-bold-r-normal-sans-12-120-75-75-m-70-iso8859-1 }
if ![info exists ffontbb] { set ffontbb \
  -*-helvetica-bold-o-normal-*-15-*-*-*-*-*-*-* }
if ![info exists mfont] { set mfont \
  -*-helvetica-medium-r-normal-*-13-*-*-*-*-*-*-* }  ;# 12bold,13med*,14med
if ![info exists pfont] { set pfont \
  -*-helvetica-bold-r-normal-*-14-*-*-*-*-*-*-* }

if { [info exists subjtmpdir] && ![file writable $subjtmpdir] } {
  set subjtmpdir /tmp/surfer.tmp.$env(USER)   ;# $env(USER)-[pid]
  if { ![file exists $subjtmpdir] } { exec mkdir $subjtmpdir } ;# else exists
}
if ![info exists promptflag] { set promptflag 0 }
if ![info exists program] { set program unknown }

### colors
#set tk_strictMotif 1
set bgcol       #c5c5bd
set actbgcol    #ececec    ;# default=#ececec(gray)
set hilbgcol    #c5c5bd
set hilcol      #333333
set selcol      green      ;# blue, green
set selbgcol    lightgreen ;# lightblue, lightgreen
set entbgcol    #dfdfdf    ;# #dfdfdf(offwhite)
set runbgcol    #b365ff    ;# darkor, #982dff(pr), #b365ff(lpr)
set runactbgcol #cc9aff    ;# orang, #b365ff(lpr), #cc9aff(llpr)
set menubg      $entbgcol
set menuactbg   #b2c2b2    ;# #d7afff(lllpr), selbgcol, #85a385, #b2c2b2(grygrn)
set menuchk     green      ;# green, #006000(dkgrn), green

proc labval { frame label defvalue tab } {
  global ffont ffontb
  set f [frame $frame.$label]
  pack $f -fill x -padx 2
  label $f.la -text $label: -width $tab -anchor e -font $ffont
  pack $f.la -side left
  label $f.lb -text $defvalue -font $ffontb
  pack $f.lb -side left
}

proc edlabval { frame label defvalue rwflag tab width {rowcol col} } {
  global ffont ffontb
  set f [frame $frame.$label]
  #pack $f -fill x -padx 2   ;# old
  if { $rowcol == "row" } { set side left }
  if { $rowcol == "col" } { set side top  }
  pack $f -side $side -fill x -padx 2
  if { ![string match none* $label] } {
    label $f.la -text $label: -width $tab -anchor e -font $ffont
    pack $f.la -side left
  }
  entry $f.e -font $ffont -width $width ;# -textvariable textvariable
  $f.e config -selectbackground green -insertbackground black
  $f.e insert 0 $defvalue
  pack $f.e -side left -fill x
  if { $rwflag == "rw" } {
    button $f.br -text READ -font $ffontb -pady 0
    pack $f.br -side left
    button $f.bw -text WRITE -font $ffontb -pady 0
    pack $f.bw -side right
  }
  if {$rwflag == "w"} {
    button $f.bw -text WRITE -font $ffontb -pady 0
    pack $f.bw -side right
  }
  if {$rwflag == "r"} {
    button $f.br -text READ -font $ffontb -pady 0
    pack $f.br -side left
  }
  if {$rwflag == "s"} {
    button $f.bs -text SET -font $ffontb -pady 0
    pack $f.bs -side right
  }
  if {$rwflag == "n"} {
    return
  }
}

proc radios { frame title buttonlabel variable value titlewidth rowcol } {
  global ffont ffontb
  if { $rowcol == "row" } { set side left }
  if { $rowcol == "col" } { set side top  }
  set f [frame $frame.a$buttonlabel]   ;# lower case
  pack $f -side $side -fill x -padx 2
  if { [string length $title] > 0 } {
    label $f.la -text $title: -width $titlewidth -anchor e -font $ffontb
    pack $f.la -side $side
  }
  radiobutton $f.ra -text $buttonlabel -variable $variable -value $value \
     -font $ffont
  pack $f.ra -side $side -anchor w
}

proc checks { frame title buttonlabel variable rowcol } {
  global ffont ffontb
  if { $rowcol == "row" } { set side left }
  if { $rowcol == "col" } { set side top }
  set f [frame $frame.a$buttonlabel]   ;# lower case
  pack $f -side $side -fill x -padx 2  ;# should toggle w/rowcol, also anchor
  if { [string length $title] > 0 } {
    label $f.la -text $title: -width 12 -anchor e -font $ffontb
    pack $f.la -side $side
  }
  checkbutton $f.ck -text $buttonlabel -variable $variable -font $ffont
  if { $rowcol == "col" } { $f.ck config -anchor w }
  pack $f.ck -side $side -fill x
}

proc buttons { frame buttonlabel function rowcol {pady 2} {padx 11} } {
  global ffontb
  if { $rowcol == "row" } { set side left }
  if { $rowcol == "col" } { set side top }
  set f [frame $frame.a$buttonlabel]  ;# lower case
  pack $f -side $side -fill x
  button $f.bu -text $buttonlabel -font $ffontb -command $function \
     -pady $pady -padx $padx
  pack $f.bu -side $side
}

proc okreplace { absfilename {altmsg none} {altretlabel none} {midbut none} } {
  global userok pfont ffontbb bgcol entbgcol
  set file [file tail $absfilename]
  set dir [file dirname $absfilename]
  set f [toplevel .okreplace -borderwidth 10]
  positionpopup $f
  wm protocol $f WM_DELETE_WINDOW killbox
  label $f.la -text REPLACE -font $ffontbb -bd 2 -relief groove -padx 7 -pady 7 
  if {$altmsg == "none"} {
    set msgtext "The file: $file in $dir exists.  Replace it?"
  } else {
    set msgtext $altmsg
  }
  #message $f.msg -text $msgtext -font $pfont -width 280 -padx 50 -pady 20
  message $f.msg -text $msgtext -font $pfont -width 370 -padx 30 -pady 20
  set b [frame $f.buttons -borderwidth 10]
  pack $f.la -side top -anchor w
  pack $f.msg $f.buttons -side top -fill x
  button $b.cancel -text "Cancel" -font $pfont -command {set userok 0}
  if {$midbut != "none"} {
    button $b.mid -text "$midbut" -font $pfont -command {set userok 2}
  }
  if {$altretlabel == "none"} {
    set buttext Replace
  } else {
    set buttext $altretlabel
  }
  button $b.replace -text "$buttext" -font $pfont -command {set userok 1}
  pack $b.replace -side right -padx 5 -anchor e
  if {$midbut != "none"} {pack $b.mid -side right -padx 5 -anchor e}
  pack $b.cancel -side right -padx 5 -anchor e
  bind .okreplace <Control-c> {set userok 0}
  bind .okreplace <Return> {set userok 1}
  foreach win "$f $f.msg $f.buttons" { $win configure -background $bgcol }
  foreach win "$f.la $b.cancel $b.replace" { $win configure -bg $entbgcol }
  if {$midbut != "none"} {$b.mid configure -bg $entbgcol}
  tkwait visibility $f
  wm protocol . WM_DELETE_WINDOW killbox
  raise $f
  focus $f
  grab $f
  tkwait variable userok
  grab release $f
  destroy $f
  wm protocol . WM_DELETE_WINDOW testclose
  return $userok
}

proc okclose { absfilename } {
  global userok pfont ffontbb bgcol entbgcol
  set file [file tail $absfilename]
  set dir [file dirname $absfilename]
  set f [toplevel .okclose -borderwidth 10]
  positionpopup $f
  wm protocol $f WM_DELETE_WINDOW killbox
  label $f.la -text CLOSE -font $ffontbb -bd 2 -relief groove -padx 7 -pady 7 
  set msgtext "Save changes to $file?"
  #message $f.msg -text $msgtext -font $pfont -width 280 -padx 50 -pady 20
  message $f.msg -text $msgtext -font $pfont -width 370 -padx 30 -pady 20
  set b [frame $f.buttons -borderwidth 10]
  pack $f.la -side top -anchor w
  pack $f.msg $f.buttons -side top -fill x
  button $b.cancel -text "Cancel" -font $pfont -command {set userok 0}
  button $b.dontsave -text "Don't Save" -font $pfont -command {set userok 1}
  button $b.save -text "Save" -font $pfont -command {set userok 2}
  pack $b.save $b.dontsave $b.cancel -side right -padx 5 -anchor e
  bind $f <Control-c> {set userok 0}
  bind $f <Return> {set userok 2}
  foreach win "$f $f.msg $f.buttons" { $win config -background $bgcol }
  foreach win "$f.la $b.cancel $b.dontsave $b.save" {$win config -bg $entbgcol}
  tkwait visibility $f
  wm protocol . WM_DELETE_WINDOW killbox      ;# rootwin kill is cancel
  raise $f            ;# 4dwm focus follows mouse (in rootwins killbox)
  focus $f
  grab $f
  tkwait variable userok
  grab release $f
  destroy $f
  wm protocol . WM_DELETE_WINDOW testclose    ;# re-arm
  return $userok
}

proc confirmalert { errortext } {
  global userok pfont ffontbb bgcol entbgcol
  set nm .confirm
  while {[info commands $nm] == "$nm"} {set nm ${nm}x}
  set f [toplevel $nm -borderwidth 10]
  positionpopup $f
  wm protocol $f WM_DELETE_WINDOW killbox
  label $f.la -text ALERT -font $ffontbb -bd 2 -relief groove -padx 7 -pady 7 
  #message $f.msg -text $errortext -font $pfont -width 280 -padx 50 -pady 20
  message $f.msg -text $errortext -font $pfont -width 370 -padx 30 -pady 20
  set b [frame $f.buttons -borderwidth 10]
  pack $f.la -side top -anchor w
  pack $f.msg $f.buttons -side top -fill x
  button $b.cancel -text "OK" -font $pfont -command {set userok 0}
  pack $b.cancel -side right -padx 5 -anchor e
  bind $f <Control-c> {set userok 0}
  bind $f <Return> {set userok 0}
  foreach win "$f $f.msg $b" { $win config -background $bgcol }
  foreach win "$f.la $b.cancel" { $win config -background $entbgcol }
  tkwait visibility $f
  wm protocol . WM_DELETE_WINDOW killbox
  raise $f
  focus $f
  #grab $f   ;# errors if another app grabs
  tkwait variable userok
  #grab release $f
  destroy $f
  wm protocol . WM_DELETE_WINDOW testclose
  return $userok
}

proc helpwin { helpfile {cols 50} {rows 10} {location std} } {
  global userok env ffont ffontbb pfont program bgcol entbgcol selbgcol
  set nm .help
  while {[info commands $nm] == "$nm"} {set nm ${nm}x}
  if {$location == "std"} {
    set helpfile $env(FREESURFER_HOME)/lib/help/$program/$helpfile
  } else {
    set helpfile $location/$helpfile
  }
  if [catch {open $helpfile}] { return }
  set f [toplevel $nm -borderwidth 10]
  positionpopup $f
  wm protocol $f WM_DELETE_WINDOW killbox
  wm resizable $f 0 1
  label $f.la -text HELP -font $ffontbb -bd 2 -relief groove -padx 7 -pady 7 
  pack $f.la -side top -anchor w

  frame $f.t
  set log [text $f.t.log -width $cols -height $rows -font $ffont \
    -borderwidth 2 -relief raised -setgrid true -selectbackground $selbgcol \
    -yscrollcommand "$f.t.scroll set"]
  scrollbar $f.t.scroll -command "$f.t.log yview"
  pack $f.t.scroll -side right -fill y
  pack $f.t.log -side left -expand true -fill y
  pack $f.t -side top -expand true -fill y
  $f.t.log config -state normal
  set id [open $helpfile]
  $f.t.log insert end [read $id]
  close $id

  set b [frame $f.but -borderwidth 1]
  pack $b -side right
  button $b.ok -text "OK" -font $pfont -command "destroy $f"
  pack $b.ok -side right -padx 5 -anchor e
  foreach win "$f $f.but $f.t.scroll" { $win config -background $bgcol }
  foreach win "$f.la $b.ok $f.t.log" { $win config -background $entbgcol }
  bind $nm <Return> "destroy $f"
  focus $f
}

proc killbox { } {
  global userok
  set userok 0
}

proc testreplace { file function } {
  if { [file exists $file] } {
    if { [okreplace $file] } {
      $function
    } else {
      return
    }
  } else {
    $function
  }
}

proc prompt { } {
  global promptflag
  if { $promptflag } { puts -nonewline "% " }
  flush stdout
}

proc tmpcontrols { title variablelist } {
  global ffontb userok
  set frame .tmpcontrols_popup
  set f [toplevel $frame -borderwidth 10]
  wm title $f $title
  positionpopup $f
  wm protocol $f WM_DELETE_WINDOW killbox
  label $f.la -text $title -font $ffontb
  pack $f.la -side top
  set longest 0
  foreach var $variablelist {
    if {[string length $var] > $longest} {set longest [string length $var]}
  }
  foreach var $variablelist {
    edlabval $f $var 0 n [expr $longest+2] 5
    $f.$var.e config -textvariable $var
  }
  buttons $f "RUN SCRIPT" { set userok 1 } col
  buttons $f "CANCEL"     { set userok 0 } col
  focus $f
  grab $f
  tkwait variable userok
  grab release $f
  destroy $f
  return $userok
}

proc controls { title variablelist buttonname buttoncommand } {
  global ffontb
  set frame .a"$title"_popup
  set f [toplevel $frame -borderwidth 10]
  wm title $f $title
  positionpopup $f
  label $f.la -text $title -font $ffontb
  pack $f.la -side top
  set longest 0
  foreach var $variablelist {
    if {[string length $var] > $longest} {set longest [string length $var]}
  }
  foreach var $variablelist {
    edlabval $f $var 0 n [expr $longest+2] 10
    $f.$var.e config -textvariable $var
  }
  buttons $f $buttonname $buttoncommand col
  buttons $f "CANCEL" "destroy $f" col
}

proc fixcolors { {skipwin .mbar} } {
  global bgcol actbgcol hilbgcol hilcol selcol selbgcol entbgcol
  foreach win [info commands .*] {
    if [string match $skipwin* $win] { continue }
    if [string match *.shell $win] { continue }
    $win config -background $bgcol -highlightbackground $hilbgcol \
       -highlightcolor $hilcol 
  }
  foreach win [info commands *.ck] {$win config -selectcolor $selcol}
  foreach win [info commands *.ra] {$win config -selectcolor $selcol}
  foreach win "[info commands *.e] [info commands *.entry]" {
    $win config -background $entbgcol -selectbackground $selbgcol
  }
 foreach win "[info commands *.a*.bu]" {$win config -activebackground $actbgcol}
}

proc positionpopup { win {size {}} } {
  set winx [expr [winfo rootx .]+45]   ;# rel to root
  set winy [expr [winfo rooty .]+25]
  wm geometry $win $size+${winx}+${winy}
}

proc menusetup { mbar } {
  global Menu menubg mfont
  frame $mbar -borderwidth 1 -relief raised -background $menubg
  pack $mbar -side top -anchor w -fill x
  set Menu(mbar) $mbar
  set Menu(uid) 0
}

proc menuhead { label side {accelindex -1}} {
  global Menu menubg menuactbg mfont
  if [info exists Menu(menu,$label)] { return }
  set name $Menu(mbar).mb$Menu(uid)
  set mname $name.me
  incr Menu(uid)
  set mb [menubutton $name -text $label -menu $mname \
         -background $menubg -activebackground $menuactbg -font $mfont]
  if {$accelindex >= 0} {$mb config -underline $accelindex}
  pack $mb -side $side -padx 5
  set menu \
    [menu $mname -tearoff 0 -background $menubg -activebackground $menuactbg]
  set Menu(menu,$label) $menu
}

proc menucmd { mname label cmd {accelindex -1}} {
  global Menu menubg menuactbg mfont
  if [catch {set Menu(menu,$mname)} menu] { return }
  $menu add command -label $label -command $cmd \
     -background $menubg -activebackground $menuactbg -font $mfont
  if {$accelindex >= 0} {$menu config -underline $accelindex}
}

proc menucheck { mname label var {cmd {}} } {
  global Menu menubg menuactbg mfont menuchk
  if [catch {set Menu(menu,$mname)} menu] { return }
  $menu add check -label $label -command $cmd -variable $var -font $mfont \
       -background $menubg -activebackground $menuactbg -selectcolor $menuchk
}
 
proc menusepar { mname } {
  global Menu
  if [catch {set Menu(menu,$mname)} menu] { return }
  $menu add separator
}

proc menubind { what sequence mname label } {
  global Menu
  if [catch {set Menu(menu,$mname)} menu] { return }
  if [catch {$menu index $label} index] { return }
  set cmd [$menu entrycget $index -command]
  bind $what $sequence $cmd
  $menu entryconfigure $index -accelerator $sequence
}

proc menudisable { mname label } {
  global Menu
  if [catch {set Menu(menu,$mname)} menu] { return }
  if [catch {$menu index $label} index] { return }
  $menu entryconfigure $index -state disabled
}

proc menuenable { mname label } {
  global Menu
  if [catch {set Menu(menu,$mname)} menu] { return }
  if [catch {$menu index $label} index] { return }
  $menu entryconfigure $index -state normal
}
