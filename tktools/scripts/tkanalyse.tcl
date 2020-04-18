#! /usr/bin/tclsh

##
## tkanalyse.tcl
##
##
## Copyright (C) 1996 Martin Sereno and Anders Dale
## Copyright (C) 2002-2011, CorTechs Labs, Inc. (La Jolla, CA) and
## The General Hospital Corporation (Boston, MA).
##
## Terms and conditions for use, reproduction, distribution and contribution
## are found in the 'FreeSurfer/CorTechs Software License Agreement' contained
## in the file 'license.cortechs.txt' found in the FreeSurfer distribution,
## and here:
##
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferCorTechsLicense
##
## Reporting: freesurfer@nmr.mgh.harvard.edu
##
############################################################################

set program tkanalyse
set ffont -b&h-lucidatypewriter-medium-r-normal-sans-12-120-75-75-m-70-iso8859-1
set ffontb -b&h-lucidatypewriter-bold-r-normal-sans-12-120-75-75-m-70-iso8859-1
set ffontbb -*-helvetica-bold-*-*-*-15-*-*-*-*-*-*-*
set sfont -b&h-lucidatypewriter-medium-*-*-*-10-*-*-*-*-*-*-*

set sclenx 140
set sclenfx 180
set scleny 127

### analyse vars set at startup
# $home $session $subject $registerdat $analysedat

### literal file list,subdirs for setfile
set rgbdir rgb
set rgbfiles { rgb }

### default abbreviations
set rgbabbrev /tmp/tkanalyse.rgb

### events
set userok 0
set lastslice $curslice

### slice (slider!) limits
set structmax [expr $nslices-1]
set functmax [expr $nperslice-1]

### dispmode's (4=phase defunct: same as pphase thresh=0)
set struct  1
set funct   2 
set pval    3
set KS      5
set pphase  6

### source standard widget wrapper code
source $env(FREESURFER_HOME)/lib/tcl/wrappers.tcl

############################################################################
proc setfile { varName value } {  ;# makes "varName"abbrev if it doesn't exist
  upvar $varName localvar
  upvar ${varName}abbrev localabbrev    ;# setfile rgb ~/tmp/myrgb
  global home subject session
  global rgbdir
  global rgbfiles
  ### set subdir using name of var and premade lists
  if { [lsearch -exact $rgbfiles $varName] >= 0 } {
    set subdir $rgbdir
  } else {
    puts "bad arg: don't know how to set file type: $varName"
    puts "  setfile {$rgbfiles} value"
    prompt
    return
  }
  ### do expansions and make fullname
  if { [string range $value 0 0] == "/"} {
    set fullname $value
  } elseif { [string range $value 0 0] == "~"} {
    if { [string length $value] > 1 && [string range $value 1 1] != "/" } {
      puts "tkmedit: can't reset subject home dir without restarting"
      prompt
      return
    }
    if { [string length $value] < 3 } {
      puts "tkmedit: no filename specified for setfile"
      prompt
      return
    }
    set tildegone [string range $value 2 end]
    set subdir [file dirname $tildegone]      ;# overwrite, may be multilevel
    set filename [file tail $tildegone]
    if { $subdir == "." } {
      set fullname $home/$subject/$filename
    } else {
      set fullname $home/$subject/$subdir/$filename
    }
  } elseif { [string range $value 0 0] == "*"} {
    set stargone [string range $value 2 end]
    set fullname $session/$stargone
  } else {  ;# relative (guess session vs. subjects)
    if { $subdir == "." } {
      set fullname $session/$value
    } else {
      set fullname $session/$subdir/$value
    }
  }
  set localvar $fullname
  #puts $fullname

  ### attempt to re-abbrev (first ~, then *, else set absolute)
  set homename $home/$subject
  set homelen [string length $homename]
  set endhome [incr homelen -1]
  set begintail [incr homelen 1]
  if { $homename == [string range $fullname 0 $endhome] } {
    set localabbrev ~[string range $fullname $begintail end]
    return
  }
  set sessionlen [string length $session]
  set endsession [incr sessionlen -1]
  set begintail [incr sessionlen 1]
  if { $session == [string range $fullname 0 $endsession] } {
    set localabbrev *[string range $fullname $begintail end]
    return
  }
  set localabbrev $fullname
}

proc tosingle { } {
  global doublebufferflag openglwindowflag
  if {!$openglwindowflag} {
    puts "tksurfer: ### tosingle failed: no gl window open"
    prompt
    return
  }
  set doublebufferflag FALSE
  singlebuffer
  gconfig
  redraw
}

proc todouble { } {
  global doublebufferflag openglwindowflag
  if {!$openglwindowflag} {
    puts "tksurfer: ### todouble failed: no gl window open"
    prompt
    return
  }
  set doublebufferflag TRUE
  doublebuffer
  gconfig
  redraw
}

proc fixbrightbounds { varName index op } {
  global dispmode
  global struct funct pval KS pphase
  set f .right.sc.bri.bot.left.y
  switch $op {
    w {
      switch $dispmode \
 $struct {$f config -from 0.13 -to 0.0 -var sfact  -tick 0.13 -res 0.01} \
 $funct  {$f config -from 0.13 -to 0.0 -var ffact  -tick 0.13 -res 0.01} \
 $pval   {$f config -from 25.0 -to 0.0 -var pfact  -tick 25.0 -res 0.5}  \
 $KS     {$f config -from 1200 -to 0   -var KSfact -tick 1200 -res 25}   \
 $pphase {$f config -from 25.0 -to 0.0 -var pfact  -tick 25.0 -res 0.5}
    }
  }  
}

proc fixthreshbounds { varName index op } {
  global dispmode
  global pval KS pphase
  set f .right.sc.thr.bot.left.y
  switch $op {
    w {
      switch $dispmode \
 $pval   {$f config -from 25.0 -to 0.0 -var pthresh  -tick 25.0 -res 0.1}  \
 $KS     {$f config -from 2.5  -to 0.0 -var KSthresh -tick 2.5  -res 0.02} \
 $pphase {$f config -from 25.0 -to 0.0 -var pthresh  -tick 25.0 -res 0.1}
    }
  }
}

proc testclose { } {
  if {0} {   ;# C editedflag
    set resp [okclose somefile]
    if {$resp > 1} {      }
    if {$resp > 0} { exit }
  } else {
    exit
  }
}

############################################################################
wm title . "tkanalyse ($subject--[file tail [exec pwd]])"
wm geometry . 447x323+117+0
wm protocol . WM_DELETE_WINDOW testclose
wm resizable . 0 0

frame .head
pack .head -side top
  frame .head.left
  pack .head.left -side left
  frame .head.right
  pack .head.right -side left

frame .left -borderwidth 1
pack .left -side left

  frame .left.view
  pack .left.view -side top -fill y
  foreach v { struct funct } {   ;# two identical slice panels
    frame .left.view.$v -borderwidth 2 -relief groove
    pack .left.view.$v -side top -fill x
      frame .left.view.$v.top
      pack .left.view.$v.top -side top
      frame .left.view.$v.bot
      pack .left.view.$v.bot -side top
  }

  frame .left.mode -borderwidth 2 -relief groove
  pack .left.mode -side left

  frame .left.spacer
  pack .left.spacer -side left

  frame .left.lr
  pack .left.lr -side left
    frame .left.lr.plot -borderwidth 2 -relief groove
    pack .left.lr.plot -side top
    frame .left.lr.squ
    pack .left.lr.squ -side top

frame .right -borderwidth 1
pack .right -side left

  frame .right.point -borderwidth 2 -relief groove
  pack .right.point -side top -fill x

  frame .right.compute -borderwidth 2 -relief groove
  pack .right.compute -side top -fill x
    frame .right.compute.top
    pack .right.compute.top -side top
    frame .right.compute.mid
    pack .right.compute.mid -side top
    frame .right.compute.bot
    pack .right.compute.bot -side top
  
  frame .right.sc
  pack .right.sc -side top

    frame .right.sc.bri -borderwidth 2 -relief groove
    pack .right.sc.bri -side left
      frame .right.sc.bri.top
      pack .right.sc.bri.top -side top
      frame .right.sc.bri.bot
      pack .right.sc.bri.bot -side top
        frame .right.sc.bri.bot.left
        pack .right.sc.bri.bot.left -side left
        frame .right.sc.bri.bot.right
        pack .right.sc.bri.bot.right -side left

    frame .right.sc.thr -borderwidth 2 -relief groove
    pack .right.sc.thr -side left
      frame .right.sc.thr.top
      pack .right.sc.thr.top -side top
      frame .right.sc.thr.bot
      pack .right.sc.thr.bot -side top
        frame .right.sc.thr.bot.left
        pack .right.sc.thr.bot.left -side left
        frame .right.sc.thr.bot.right
        pack .right.sc.thr.bot.right -side left

  frame .right.rgb -borderwidth 2 -relief groove
  pack .right.rgb -side top

############################################################################
### title 
#set f .head.left
#buttons $f "POP GL" { pop_gl_window } row 0 5
set f .head.right
edlabval $f "scan" [exec pwd] n 6 55
$f.scan.e config -font $ffontb -state disabled
$f.scan.e xview end

### STRUCT image: button
set f .left.view.struct.top
buttons $f "STRUCT IMAGE" { if {$dispmode == $funct} {set dispmode $struct}; \
                            goto_slice $curslice; redraw; \
                            focus .left.view.struct.bot.sc } row
edlabval $f "slice" 0 n 6 3
$f.slice.e config -textvariable curslice -font $sfont
bind $f.slice.e <Return> { if {$dispmode == $funct} {set dispmode $struct}; \
                           goto_slice $curslice; redraw; update idletasks }
### STRUCT image: scale
set f .left.view.struct.bot
scale $f.sc -from 0 -to $structmax -length $sclenx -variable curslice \
   -orient horizontal -tickinterval $structmax -showvalue false -font $sfont \
   -width 11 -resolution 1
pack $f.sc -side top
bind $f.sc <ButtonRelease> {if {$dispmode == $funct} {set dispmode $struct}; \
                            if {$curslice != $lastslice} {goto_slice $curslice;\
                           redraw; update idletasks; set lastslice $curslice };\
                           focus .left.view.struct.bot.sc }
bind $f.sc <B1-Motion> { if {$dispmode == $funct} {set dispmode $struct}; \
                         if {$curslice != $lastslice} {goto_slice $curslice; \
                           redraw; update idletasks; set lastslice $curslice } }
### FUNCT images: button
set f .left.view.funct.top
buttons $f "FUNCT IMAGES" { set dispmode $funct; f_update; redraw; \
                            focus .left.view.funct.bot.sc  } row
edlabval $f "image" 0 n 6 3
$f.image.e config -textvariable fimc -font $sfont
bind $f.image.e <Return> { set dispmode $funct; f_update; \
                           goto_slice $fimc; redraw; update idletasks }
### FUNCT images: scale
set f .left.view.funct.bot
scale $f.sc -from 0 -to $functmax -length $sclenfx -variable fimc \
   -orient horizontal -tickinterval $functmax -showvalue false -font $sfont \
   -width 11 -resolution 1
pack $f.sc -side top
bind $f.sc <ButtonRelease> { set dispmode $funct; f_update; \
                             goto_slice $fimc; redraw; update idletasks; \
                             focus .left.view.funct.bot.sc }
bind $f.sc <B1-Motion> { set dispmode $funct; f_update; \
                         goto_slice $fimc; redraw; update idletasks}

### mode field
set f .left.mode
label $f.la -text "MODE" -font $ffontb -pady 1
pack $f.la -side top
radios $f  ""  "struct"  dispmode  $struct  7  col
radios $f  ""  "funct"   dispmode  $funct   7  col
radios $f  ""  "pval"    dispmode  $pval    7  col
radios $f  ""  "pphase"  dispmode  $pphase  7  col
radios $f  ""  "KS"      dispmode  $KS      7  col
bind $f.astruct.ra <ButtonRelease> { redraw }
bind $f.afunct.ra <ButtonRelease>  { f_update; set dispmode $funct; redraw }
bind $f.apval.ra <ButtonRelease>   { redraw }
bind $f.apphase.ra <ButtonRelease> { redraw }
bind $f.aKS.ra <ButtonRelease>     { redraw }

### spacer
set f .left.spacer
label $f.la -text ""
pack $f.la -side left

### plot field
set f .left.lr.plot
label $f.la -text "Plot Pixel" -font $ffontb -pady 0
pack $f.la -side top
buttons $f "RAW"        { f_update; plot_timeseries } col 0 11
buttons $f "FOURIER"    { f_update; plot_fourier }    col 0 11
buttons $f "KS"         { f_update; plot_KS }         col 0 11
buttons $f "COND STATS" { f_update; plot_condition_stats } col 0 11

### misc entries
set f .left.lr.squ
edlabval $f "fsquash" 0 n 8 4
$f.fsquash.e config -textvariable fsquash
bind $f.fsquash.e <Return> { set_scale; redraw }
edlabval $f "fthresh" 0 n 8 4
$f.fthresh.e config -textvariable fthresh
bind $f.fthresh.e <Return> { set_scale; redraw }

### point panel
set f .right.point
buttons $f "POP GL" { pop_gl_window } row 0 4
#"$f.aPOP GL.bu" config -font $ffontbb
label $f.la -text " " -font $ffont  ;# space
pack $f.la -side left
#buttons $f "SAVE PNT" { write_point } row 2 8
buttons $f "SAVE PNT" { write_point; \
                catch { send tksurfer select_orig_vertex_coordinates } } row 2 8
buttons $f "GOTO PNT" { goto_point; redraw; set curslice $curslice } row 2 8

### compute field
set f .right.compute.top
label $f.la -text "COMPUTE SLICE  " -font $ffontb -pady 1
pack $f.la -side left
checks $f "" "rmtrend" trendflag row
set f .right.compute.mid
buttons $f "FOURIER"  { f_update; compute_fourier_image; \
                                if {$dispmode < $pval} {set dispmode $pval}; \
                                                               redraw } row 2 7
buttons $f "KS"       { f_update; compute_KS_image; \
                                              set dispmode $KS; redraw} row 2 7
buttons $f "COND"     { f_update; compute_condition_stats; redraw} row 2 7

buttons $f "COLLAPSE" { f_update; collapse_f_timecourses; redraw }  row 2 7

set f .right.compute.bot
buttons $f "SMOOTH RAW"     { f_update; smooth_f_images; redraw }       row 2 9
buttons $f "SMOOTH FOURIER" { f_update; smooth_fourier_images; redraw } row 2 9

### brightness panel
set f .right.sc.bri.top
label $f.la -text "BRIGHTNESS" -font $ffontb -pady 1
pack $f.la -side top
set f .right.sc.bri.bot.left
scale $f.y -from 25.0 -to 0.0 -length $scleny -variable pfact \
   -orient vertical -tickinterval 25.0 -showvalue false -font $sfont \
   -width 11 -resolution 0.5
pack $f.y -side left
bind $f.y <ButtonRelease> { redraw }
trace variable dispmode w fixbrightbounds
# entries
set f .right.sc.bri.bot.right
edlabval $f "s"  0 n 3 5
edlabval $f "f"  0 n 3 5
edlabval $f "p"  0 n 3 5
edlabval $f "ks" 0 n 3 5
$f.s.e config -textvariable sfact -font $sfont
$f.f.e config -textvariable ffact -font $sfont
$f.p.e config -textvariable pfact -font $sfont
$f.ks.e config -textvariable KSfact -font $sfont
$f.s.la config -font $sfont
$f.f.la config -font $sfont
$f.p.la config -font $sfont
$f.ks.la config -font $sfont -text "KS:"
bind $f.s.e <Return> { redraw }
bind $f.f.e <Return> { redraw }
bind $f.p.e <Return> { redraw }
bind $f.ks.e <Return> { redraw }

### threshold panel
set f .right.sc.thr.top
label $f.la -text "THRESHOLD" -font $ffontb -pady 1
pack $f.la -side top
# scale
set f .right.sc.thr.bot.left
scale $f.y -from 25.0 -to 0.0 -length $scleny -variable pthresh \
   -orient vertical -tickinterval 25.0 -showvalue false -font $sfont \
   -width 11 -resolution 0.1
pack $f.y -side left
bind $f.y <ButtonRelease> { redraw }
trace variable dispmode w fixthreshbounds
# entries
set f .right.sc.thr.bot.right
label $f.la1 -font $ffont  ;# space
pack $f.la1 -side top
label $f.la2 -font $ffont  ;# space
pack $f.la2 -side top
edlabval $f "p"  0 n 3 5
edlabval $f "ks" 0 n 3 5
$f.p.e config -textvariable pthresh -font $sfont
$f.ks.e config -textvariable KSthresh -font $sfont
$f.p.la config -font $sfont
$f.ks.la config -font $sfont -text "KS:"
bind $f.p.e <Return> { redraw }
bind $f.ks.e <Return> { redraw }

### save rgb button,field
set f .right.rgb
edlabval $f  "rgb"   $rgb   w 4 19
$f.rgb.e config -textvariable rgbabbrev
$f.rgb.bw config -command { if {$doublebufferflag} {tosingle}; \
                            setfile rgb [.right.rgb.rgb.e get]; \
                            testreplace $rgb save_rgb }

############################################################################
### shortcut key bindings--to find keysyms: bind . <KeyPress> { puts %K }
# commands
bind . <Alt-f> { .left.lr.plot.aTIME.bu invoke }
bind . <Alt-F> { .left.lr.plot.aSPECT.bu invoke }
bind . <Alt-k> { .left.lr.plot.aKS.bu invoke }
bind . <Alt-K> { .right.compute.top.aKS.bu invoke }
bind . <Alt-s> { .right.compute.top.aFOURIER.bu invoke }
bind . <Alt-r> { ".right.point.aSAVE PNT.bu" invoke }
bind . <Alt-w> { ".right.point.aGOTO PNT.bu" invoke }
bind . <Alt-R> { read_f_images; redraw }
bind . <Alt-a> { ".right.compute.top.aSMOOTH FOURIER.bu" invoke }
bind . <Alt-A> { ".right.compute.top.aSMOOTH RAW.bu" invoke }
bind . <Alt-T> { ".right.compute.top.aCOLLAPSE.bu" invoke }
# planes
bind . <Alt-x> {  }
bind . <Alt-y> {  }
bind . <Alt-z> {  }
# flags
bind . <Alt-Key-M> { set trendflag 1; redraw }
bind . <Alt-m> { set trendflag 0; redraw }
# slices
bind . <Alt-Right> {if {$dispmode == $funct} {focus .left.view.funct.bot.sc \
                      } else { focus .left.view.struct.bot.sc};up_slice;redraw }
bind . <Alt-Left> {if {$dispmode == $funct} {focus .left.view.funct.bot.sc \
                    } else { focus .left.view.struct.bot.sc};down_slice;redraw }
# modes
bind . <Alt-Key-1> { set dispmode $struct; redraw }
bind . <Alt-Key-2> { f_update; set dispmode $funct; redraw }
bind . <Alt-Key-3> { set dispmode $pval;   redraw }
bind . <Alt-Key-5> { set dispmode $KS;     redraw }
bind . <Alt-Key-6> { set dispmode $pphase; redraw }
# contrast
bind . <Alt-asterisk> { set fsquash [expr $fsquash * 1.1]; set_scale; redraw }
bind . <Alt-slash> { set fsquash [expr $fsquash / 1.1]; set_scale; redraw }
bind . <Alt-plus> { set fthresh [expr $fthresh + 0.05]; set_scale; redraw }
bind . <Alt-minus> { set fthresh [expr $fthresh - 0.05]; set_scale; redraw }
bind . <Alt-Up> { up_bright; redraw }
bind . <Alt-Down> { down_bright; redraw }
# thresh
bind . <Alt-braceleft>    { down_thresh; redraw }
bind . <Alt-braceright>   { up_thresh; redraw }
#bind . <Alt-bracketleft>  {  }
#bind . <Alt-bracketright> {  }
# alternate thresh
#bind . <Alt-KP_Right> {  }
#bind . <Alt-KP_Left>  {  }
bind . <Alt-KP_Up>    { up_thresh; update idletasks; redraw }
bind . <Alt-KP_Down>  { down_thresh; update idletasks; redraw }

############################################################################
### right-click help
#bind .draw.main.aREDRAW.bu <B3-ButtonRelease> { helpwin redraw }

############################################################################
puts "tkanalyse.tcl: startup done"

