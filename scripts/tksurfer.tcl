#! /usr/local/bin/tclsh7.4
############################################################################
#  Copyright (c) 1999 CorTechs - All rights reserved
############################################################################
set program tksurfer
### fonts from wrappers.tcl

### slider sizing
set sclenx 120  ;# (125)
set scleny 90   ;# (100)

### some tksurfer file vars set at startup
# $home $session $script $subject $hemi $ext $insurf
# $outsurf $targsurf $curv $targcurv $sulc $patch $annot $area $fs $fm $dip $dec
 
### literal subdirs for setfile
set scriptsdir scripts
set surfdir surf
set fsdir fs
set labeldir label
set bemdir bem
set rgbdir rgb

### file lists for setfile
set sessionfiles { val }
set surffiles \
   { insurf outsurf curv sulc patch annot area origcoords ellcoords vrmlsurf }
set fsfiles { fs fm } 
set scriptsfiles { script } 
set labelfiles { label } 
set bemfiles { dip dec } 
set targfiles { targsurf targcurv } 
set rgbfiles { rgb }
set rgbdirs { named_rgbdir num_rgbdir }   ;# not files but same expansions

### current transforms
set xrot 0
set yrot 0
set zrot 0
set xtrans 0
set ytrans 0
set scalepercent 100 

### lights: can't be written back (floats in gl struct)
set light0 0.4
set light1 0.0
set light2 0.6
set light3 0.2
#set offset 0.25

### events
set userok 0
set blinkflag FALSE
set initdelay $blinkdelay

### misc defaults
set surfcolor 1
#set colscale 0
#set fthresh 0.1
set shrinksteps 5
set smoothsteps 5
set smoothtype val

### source standard widget wrapper code
source $env(CSURF_DIR)/lib/tcl/wrappers.tcl

############################################################################
proc setsession { name } {
  global session home subject
  if { [string range $name 0 0] == "~" } {
    if { [string length $name] == 1 } {
      set session $home/$subject
      return
    }
    if { [string length $name] < 3 || [string range $name 1 1] != "/"} {
      puts "bad tilde name for setsession"
      prompt
      return
    }
    set tildegone [string range $name 2 end]
    set subdir [file dirname $tildegone]
    set sessiontail [file tail $tildegone]
    if { $subdir == "." } {
      set session $home/$subject/$sessiontail
    } else {
      set session $home/$subject/$subdir/$sessiontail
    }
  } elseif { [string range $name 0 0] == "/" } {
    set session $name
  } else {
    puts "invalid name for setsession: must start with ~/ or /"
    prompt
    return
  }
}

proc setfile { varName value } {  ;# makes "varName"abbrev if it doesn't exist!
  upvar $varName localvar                 ;# setfile val */4/sig3-lh.w
  upvar ${varName}abbrev localabbrev      ;# setfile outsurf rh.new
  global home subject                     ;# setfile val /tmp/data.w
  global env                              ;# setfile script local.tcl
  global session                          ;# setfile named_rgbdir ~/tmp/myrgb
  global sessionfiles                     ;# setfile script #/builtin.tcl
  global surfdir fsdir scriptsdir bemdir rgbdir labeldir
  global surffiles fsfiles scriptsfiles bemfiles rgbfiles labelfiles targfiles
  global rgbdirs
  ### set subdir using name of var and premade lists
  set onebrainflag TRUE
  if { [lsearch -exact $surffiles $varName] >= 0 } {
    set subdir $surfdir
  } elseif { [lsearch -exact $fsfiles $varName] >= 0 } {
    set subdir $fsdir
  } elseif { [lsearch -exact $scriptsfiles $varName] >= 0 } {
    set subdir $scriptsdir
  } elseif { [lsearch -exact $labelfiles $varName] >= 0 } {
    set subdir $labeldir
  } elseif { [lsearch -exact $bemfiles $varName] >= 0 } {
    set subdir $bemdir
  } elseif { [lsearch -exact $rgbfiles $varName] >= 0 } {
    set subdir $rgbdir
  } elseif { [lsearch -exact $sessionfiles $varName] >= 0 } {
    set subdir .
  } elseif { [lsearch -exact $rgbdirs $varName] >= 0 } {
    set subdir .
  } elseif { [lsearch -exact $targfiles $varName] >= 0 } {
    set subdir $surfdir
    set onebrainflag FALSE
  } else {
    puts "bad arg: don't know how to set file type: $varName"
    puts "  setfile {$surffiles} value"
    puts "  setfile {$labelfiles} value"
    puts "  setfile {$fsfiles} value"
    puts "  setfile {$sessionfiles} value"
    puts "  setfile {$bemfiles} value"
    puts "  setfile {$targfiles} value (~name OK)"
    puts "  setfile {$rgbfiles} value"
    puts "  setfile {$rgbdirs} value"
    puts "  setfile {$scriptsfiles} value"
    prompt
    return
  } 
  ### do expansions and make fullname
  if { [string range $value 0 0] == "/"} {
    set fullname $value
  } elseif { [string range $value 0 0] == "~"} {
    if { [string range $value 1 1] != "/" } {
      if {$onebrainflag} {
        puts "other subjects only allowed for targfiles: $targfiles" 
        prompt
        return
      }
      set tildegone [string range $value 1 end]
      set tmpsubject [file dirname $tildegone]
      set subdir [file dirname [file tail $tildegone]]
      set filename [file tail [file tail $tildegone]]
    } else {
      set tildegone [string range $value 2 end]
      set tmpsubject $subject
      set subdir [file dirname $tildegone]    ;# overwrite, may be multilevel
      set filename [file tail $tildegone]
    }
    if { $subdir == "." } { 
      set fullname $home/$tmpsubject/$filename
    } else {
      set fullname $home/$tmpsubject/$subdir/$filename
    }
  } elseif { [string range $value 0 0] == "*"} {
    set stargone [string range $value 2 end]
    set fullname $session/$stargone
  } elseif { [string range $value 0 0] == "#" } {
    set poundgone [string range $value 2 end]
    set fullname $env(CSURF_DIR)/lib/tcl/$poundgone
  } else {  ;# relative (guess session vs. subjects vs. pwd)
    if {$subdir == $surfdir} {
      set fullname $home/$subject/$subdir/$value
    } elseif {$subdir == $bemdir} {
      set fullname $home/$subject/$subdir/$value
    } elseif {$subdir == $scriptsdir} {
      set fullname [exec pwd]/$value
    } elseif { $subdir == "." } {
      set fullname $session/$value
    } else {
      set fullname $session/$subdir/$value
    }
  }
  set localvar $fullname
  #puts $fullname

  ### attempt to re-abbrev (first ~, then *, else #, else set absolute)
  if {$varName != "script"} {
    set homename $home/$subject
    set homelen [string length $homename]
    set endhome [incr homelen -1]
    set begintail [incr homelen 1]
    if { "$homename" == "[string range $fullname 0 $endhome]" } {
      set localabbrev ~[string range $fullname $begintail end]
      return
    }
    set sessionlen [string length $session]
    set endsession [incr sessionlen -1]
    set begintail [incr sessionlen 1]
    if { "$session" == "[string range $fullname 0 $endsession]" } {
      set localabbrev *[string range $fullname $begintail end]
      return
    }
  } else {
    set scripthome $env(CSURF_DIR)/lib/tcl
    set scripthomelen [string length $scripthome]
    set endscript [incr scripthomelen -1]
    set begintail [incr scripthomelen 1]
    if { "$scripthome" == "[string range $fullname 0 $endscript]" } {
      set localabbrev #[string range $fullname $begintail end]
      return
    }
  }
  set localabbrev $fullname
}

proc transform { } {
  global xrot yrot zrot xtrans ytrans scalepercent
  rotate_brain_x $xrot
  rotate_brain_y $yrot
  rotate_brain_z $zrot
  translate_brain_x $xtrans
  translate_brain_y $ytrans
  scale_brain [expr $scalepercent/100.0]
}

proc resettransform { } {
  global xrot yrot zrot xtrans ytrans scalepercent
  set xrot 0; set yrot 0; set zrot 0
  set xtrans 0; set ytrans 0
  set scalepercent 100
}

proc mriforce { updown } {
  global stressthresh mmid whitemid graymid
  if { $updown == "up" } {
    set stressthresh [expr $stressthresh * 1.1]
    set mmid [expr $mmid + 1.0]
    set whitemid [expr $whitemid + 1.0]
    set graymid [expr $graymid + 1.0]
  } elseif { $updown == "down" } {
    set stressthresh [expr $stressthresh / 1.1]
    set mmid [expr $mmid - 1.0]
    set whitemid [expr $whitemid - 1.0]
    set graymid [expr $graymid - 1.0]
  }
}

proc smooth { steps type } {
  if { $type == "curv" } { 
    smooth_curv $steps
  }
  if { $type == "val" } {
    smooth_val $steps
  }
  if { $type == "sparse" } {
    smooth_val_sparse $steps
  }
}

proc redrawbutton {  } {
  global light0 light1 light2 light3 offset
  transform
  resize_window 0
  do_lighting_model $light0 $light1 $light2 $light3 $offset
  redraw
  resettransform
  prompt
}

proc restore {  } {
  global flag2d
  resettransform
  make_lateral_view
  if {$flag2d} {
    restore_zero_position
    rotate_brain_x -90
  }
}

proc restore_second {  } {
  global flag2d
  resettransform
  make_lateral_view_second
  if {$flag2d} {
    puts "no flat second surface yet"
    prompt
  }
}

proc fixdecname { varName index op } {
  global dec decabbrev dip_spacing hemi
  # trace execs in context of invoking code (global here so no uplevel setfile)
  setfile dec $decabbrev                                ;# read curr
  if {$dip_spacing > 0.999} {
    set dectail [format "%s-%d.dec" $hemi $dip_spacing] ;# fix tail one way
  } else {
    set dectail [format "%s-%3.2f.dec" $hemi $dip_spacing] ;# fix tail other way
  }  
  set dec [file dirname $dec]/$dectail                 ;# mod abs name
  setfile dec $dec                                     ;# update abbrev
}

proc fixcolscale { varName index op } {    ;# TODO: floss set*color
  global colscale complexvalflag
  if {$colscale == 0} { set complexvalflag 1 }  ;# eccen/polar wheeel
  if {$colscale == 1} { set complexvalflag 1 }  ;# heat, but ignores complex
  if {$colscale == 2} { set complexvalflag 1 }  ;# cyan/blue--both work
  if {$colscale == 3} { set complexvalflag 1 }  ;# blue/green/red--both work
  if {$colscale == 4} { set complexvalflag 1 }  ;# red/green two condition
  if {$colscale == 5} { set complexvalflag 1 }  ;# gray--both work
  if {$colscale == 6} { set complexvalflag 0 }  ;# red/blue w/white ends
  if {$colscale == 7} { set complexvalflag 0 }  ;# grn/red w/mag/cyan near 0
  if {$colscale == 8} { set complexvalflag 1 }  ;# RYGB wheel
  if {$colscale == 9} { set complexvalflag 0 }  ;# fieldsign (special signed)
}

proc save_rgb_named { name } {   ;# save_rgb_named wrapper for interface update
  global rgb named_rgbdir rel_rgbname
  set rel_rgbname $name                            ;# cp arg->global for uplevel
  uplevel {setfile rgb $named_rgbdir/$rel_rgbname} ;# echo cfunc; update abbrev
  save_rgb_named_orig $rel_rgbname                 ;# cfunc
}

proc dontsave_rgb_named { name } {  ;# just interface update
  global rgb named_rgbdir rel_rgbname
  set rel_rgbname $name
  uplevel {setfile rgb $named_rgbdir/$rel_rgbname}
}

proc findsendto { } {
  global fulltkmedit fulltkregister
  set fulltkmedit ""
  set fulltkregister ""
  # BSD: [lrange [exec ps -uxww | grep /tkmedit] 10 10]
  catch { set fulltkmedit [lrange [exec ps -af | grep /tkmedit] 7 7] }
  catch { set fulltkregister [lrange [exec ps -af | grep /tkregister] 7 7] }
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

proc macro { } {
  global tksurferinterface
  move_window 8 32
  wm withdraw .
  pack forget .mright .draw.mesh.me.aMESH.acm .shrink.cut.mini .files.val.aPUSH
  pack .subhead.hemisurf -side left
  pack .subhead.hemisurf.home -side top
  pack .subhead.hemisurf.hemisurf -side top
  pack .subhead.script   -side left
  pack .subhead.script -after .subhead.hemisurf
  pack .subhead.script.script -side top
  pack .subhead.script.aLS
  pack .files.outsurf    -before .files.curv -fill x
  pack .files.dipoles    -after .files.label -fill x
  pack .files.decimated  -after .files.dipoles
  pack .files.fieldsign  -after .files.decimated
  pack .files.fsmask     -after .files.fieldsign
  pack .files.val        -after .files.fsmask
  .files.val.la config -width 10
  pack .files.session(*) -after .files.val
  pack .draw             -after .files
  pack .draw.main.aFRAME -after .draw.main.aRESTORE -side left
  pack .draw.main.aFRAME.bu
  pack .draw.curvswitch  -after .draw.main
  pack .draw.valswitch   -after .draw.curvswitch
  pack .draw.mesh        -after .draw.valswitch
  pack .draw.switch2d    -after .draw.mesh
  pack .draw.mesh.sb     -after .draw.mesh.b
  pack .display          -after .xform
  pack .display.ri
  pack .shrink           -after .display
  pack .shrink.smooth    -after .shrink.cut
  pack .shrink.head      -after .shrink.smooth
  pack .shrink.main      -after .shrink.head
  pack .shrink.cut.aFLAT
  pack .compute          -after .shrink
  pack .draw.mesh.me     -before .draw.mesh.r
  .draw.valswitch.aw1.la config -text Scl: -width 4
  pack .draw.valswitch.aw2 -after .draw.valswitch.aw1 -side left
  pack .draw.valswitch.abi -after .draw.valswitch.aw2 -side left
  pack .draw.valswitch.aGR -after .draw.valswitch.abi -side left
  .draw.curvswitch.anone.la config -text BgSurfCol: -width 10
  pack .draw.curvswitch.aarea .draw.curvswitch.ashear -side left
  {.draw.main.aSEND PNT.bu} config -text "SEND PNT"
  {.draw.main.aGOTO PNT.bu} config -text "GOTO PNT"
  pack .shrink.cut.aROI -before .shrink.cut.aREGION
  pack .shrink.cut.aUNDO -after .shrink.cut.aLINE
  pack .shrink.cut.aINIT -after .shrink.cut.aUNDO
  .shrink.cut.aREGION.bu config -text "REGION"
  .shrink.cut.aPLANE.bu config -text "PLANE"
  .shrink.cut.aLINE.bu config -text "LINE"
  .shrink.cut.aROI.bu config -text "ROI"
  .shrink.cut.aUNDO.bu config -text "UNDO"
  .shrink.cut.aINIT.bu config -text "INIT"
  .shrink.smooth.steps.la config -text steps: -width 6
  pack .shrink.smooth.steps.asparse
  .subhead.hemisurf.home.e config -width 12
  .subhead.hemisurf.hemisurf.e config -width 12
  pack forget .draw.main.acurv
  .files.fieldsign.br config -command { setfile fs [.files.fieldsign.e get];\
                                                      read_fieldsign }
  .files.fieldsign.bw config -command { setfile fs [.files.fieldsign.e get]; \
                                   testreplace $fs  write_fieldsign }
  wm geometry . --8+32  ;# 356x985
  after 500 { wm deiconify . }
  set tksurferinterface macro
}

proc mini+ { } {
  global tksurferinterface
  move_window 8 32
  wm withdraw .
  pack forget .mright .draw.mesh.me.aMESH.acm .shrink.cut.mini .files.val.aPUSH
  pack .subhead.hemisurf -side left
  pack .subhead.hemisurf.home -side top
  pack .subhead.hemisurf.hemisurf -side top
  pack .subhead.script   -side left
  pack .subhead.script -after .subhead.hemisurf
  pack .subhead.script.script -side top
  pack .subhead.script.aLS
  pack .files.outsurf    -before .files.curv -fill x
  pack .files.fieldsign  -after .files.label
  pack .files.fsmask     -after .files.fieldsign
  pack .files.val        -after .files.fsmask
  .files.val.la config -width 10
  pack .files.session(*) -after .files.val
  pack .draw             -after .files
  pack .draw.main.aFRAME -after .draw.main.aRESTORE -side left
  pack .draw.main.aFRAME.bu
  pack .draw.curvswitch  -after .draw.main
  pack .draw.valswitch   -after .draw.curvswitch
  pack .display          -after .xform
  pack .display.ri
  pack .shrink           -after .display
  pack .shrink.smooth    -after .shrink.cut
  pack .shrink.head      -after .shrink.smooth
  pack .shrink.cut.aFLAT
  pack .draw.mesh.me     -before .draw.mesh.r
  pack .draw.mesh.sb     -after .draw.mesh.b
  .draw.valswitch.aw1.la config -text Scl: -width 4
  pack .draw.valswitch.aw2 -after .draw.valswitch.aw1 -side left
  pack .draw.valswitch.abi -after .draw.valswitch.aw2 -side left
  pack .draw.valswitch.aGR -after .draw.valswitch.abi -side left
  .draw.curvswitch.anone.la config -text BgSurfCol: -width 10
  pack .draw.curvswitch.aarea .draw.curvswitch.ashear -side left
  {.draw.main.aSEND PNT.bu} config -text "SEND PNT"
  {.draw.main.aGOTO PNT.bu} config -text "GOTO PNT"
  pack .shrink.cut.aROI -before .shrink.cut.aREGION
  pack .shrink.cut.aUNDO -after .shrink.cut.aLINE
  pack .shrink.cut.aINIT -after .shrink.cut.aUNDO
  .shrink.cut.aREGION.bu config -text "REGION"
  .shrink.cut.aLINE.bu config -text "LINE"
  .shrink.cut.aROI.bu config -text "ROI"
  .shrink.cut.aUNDO.bu config -text "UNDO"
  .shrink.cut.aINIT.bu config -text "INIT"
  .shrink.smooth.steps.la config -text steps: -width 6
  pack .shrink.smooth.steps.asparse
  .subhead.hemisurf.home.e config -width 12
  .subhead.hemisurf.hemisurf.e config -width 12
  set winlist { .files.dipoles .files.decimated .draw.switch2d \
                .shrink.main .compute .draw.main.acurv}
  foreach win $winlist { pack forget $win }
  .files.fieldsign.br config -command { setfile fs [.files.fieldsign.e get];\
                                                      read_fieldsign }
  .files.fieldsign.bw config -command { setfile fs [.files.fieldsign.e get]; \
                                   testreplace $fs  write_fieldsign }
  wm geometry . --8+32  ;# 356x728
  after 500 { wm deiconify . }
  set tksurferinterface mini+
}

proc mini { } {
  global tksurferinterface
  move_window 658 32
  wm withdraw .
  set winlist {.subhead.hemisurf .subhead.script .files.outsurf \
               .files.fsmask .files.dipoles .files.decimated .files.session(*) \
               .draw.curvswitch .draw.switch2d .draw.mesh.sb .draw.main.acurv \
               .draw.main.aFRAME .draw.main.aFRAME.bu .display.ri \
               .draw.valswitch.aw2 .draw.valswitch.abi .draw.valswitch.aGR \
               .draw.curvswitch.aarea .draw.curvswitch.ashear \
               .shrink.smooth .shrink.head .shrink.main .shrink.cut.aFLAT \
               .shrink.smooth.steps.asparse .compute }
  foreach win $winlist { pack forget $win }
  if {0} { ;# surface stuff
    pack .subhead.hemisurf
    pack .subhead.hemisurf.home -side left
    pack .subhead.hemisurf.hemisurf -side left
    .subhead.hemisurf.home.e config -width 16
    .subhead.hemisurf.hemisurf.e config -width 18
    .subhead.hemisurf.home.e xview end
    .subhead.hemisurf.hemisurf.e xview end
  } else { ;# script stuff
    pack forget .subhead.hemisurf
    pack .subhead.script
    pack .subhead.script.script -side right
    pack forget .subhead.script.aLS
    .subhead.script.script.e config -width 17
  }
  # TODO: rm/putback -bd, change scripts entry to all-tcl combobox

  pack .files.fieldsign  -after .files.label
  pack .files.val        -after .files.fieldsign
  ## rearrange
  pack .mright -before .subhead -side right
  .draw config -bd 2 -relief groove
  pack .draw -in .mright 
  pack .draw.curvswitch
  pack .draw.valswitch -after .draw.curvswitch
  pack .draw.mesh -after .draw.valswitch
  pack .draw.mesh.me    -after .draw.mesh.b
  pack .draw.mesh.me.aMESH.acm
  pack .display -after .draw
  pack .shrink -after .display
  pack .shrink.smooth
  pack .shrink.cut.mini -before .shrink.cut.aREGION -side bottom
  raise .draw
  raise .display
  raise .shrink
  pack .shrink.cut.aROI .shrink.cut.aUNDO .shrink.cut.aINIT \
     -in .shrink.cut.mini -side left
  ## squish
  {.draw.main.aSEND PNT.bu} config -text SEND
  {.draw.main.aGOTO PNT.bu} config -text GOTO
  .draw.valswitch.aw1.la config -text ColScl -width 6
  .draw.curvswitch.anone.la config -text BackgroundColor -width 15
  #.shrink.cut.aREGION.bu config -text "CUT REGION"
  .shrink.cut.aREGION.bu config -text "CUTAREA"
  .shrink.cut.aPLANE.bu config -text "CUTPLANE"
  .shrink.cut.aLINE.bu config -text "CUTLINE"
  .shrink.cut.aROI.bu config -text "FILL"
  .shrink.cut.aUNDO.bu config -text "UNDO CUT"
  .shrink.cut.aINIT.bu config -text "RECONNECTALL"
  .shrink.smooth.steps.la config -text "" -width 0
  ## auto do fsmask--reads fieldsign suffix but assumes fsmask suffix is .fm
  .files.fieldsign.br config -command { setfile fs [.files.fieldsign.e get];\
                                                      read_fieldsign; \
                     setfile fm [string range [.files.fieldsign.e get] 0 \
                    [expr [string length [.files.fieldsign.e get]] - 4]].fm; \
                                                      read_fsmask }
  .files.fieldsign.bw config -command { setfile fs [.files.fieldsign.e get]; \
                     setfile fm [string range [.files.fieldsign.e get] 0 \
                    [expr [string length [.files.fieldsign.e get]] - 4]].fm; \
                                   testreplace $fm  write_fsmask; \
                                   testreplace $fs  write_fieldsign }
  pack .files.val.aPUSH -before .files.val.la -side left
  .files.val.la config -width 4
  wm geometry . +658+664
  after 500 { wm deiconify . }
  set tksurferinterface mini
}

proc csurf { } {
  global tksurferinterface
  move_window 658 32
  set smoothtype val
  puts "using csurf interface"
  wm withdraw .
  set winlist {.subhead.hemisurf .subhead.script .files.outsurf \
               .files.fsmask .files.dipoles .files.decimated .files.val \
               .files.session(*) .draw.valswitch .files.val.aPUSH \
               .draw.curvswitch .draw.switch2d .draw.mesh.sb .draw.main.acurv \
               .draw.main.aFRAME .draw.main.aFRAME.bu .display.ri \
               .draw.valswitch.aw2 .draw.valswitch.abi .draw.valswitch.aGR \
               .draw.curvswitch.aarea .draw.curvswitch.ashear \
               .shrink.smooth .shrink.head .shrink.main .shrink.cut.aFLAT \
               .shrink.smooth.steps.asparse .compute .files.fieldsign \
               .display.le.acomplexval.ck .display.le.arevphase.ck \
               .draw.mesh
               .shrink.smooth.steps.acurv.ra .shrink.smooth.steps.aval.ra}
  foreach win $winlist { pack forget $win }
  if {0} { ;# surface stuff
    pack .subhead.hemisurf
    pack .subhead.hemisurf.home -side left
    pack .subhead.hemisurf.hemisurf -side left
    .subhead.hemisurf.home.e config -width 16
    .subhead.hemisurf.hemisurf.e config -width 18
    .subhead.hemisurf.home.e xview end
    .subhead.hemisurf.hemisurf.e xview end
  } else { ;# script stuff
    pack forget .subhead.hemisurf
#    pack .subhead.script
#    pack .subhead.script.script -side right
#    pack forget .subhead.script.aLS
#    .subhead.script.script.e config -width 17
  }
  # TODO: rm/putback -bd, change scripts entry to all-tcl combobox

#  pack .files.fieldsign  -after .files.label
#  pack .files.val        -after .files.label
  ## rearrange
  pack .mright -before .subhead -side right
  .draw config -bd 2 -relief groove
  pack .draw -in .mright 
  pack .draw.curvswitch
#  pack .draw.valswitch -after .draw.curvswitch
#  pack .draw.mesh -after .draw.valswitch
#  pack .draw.mesh.me    -after .draw.mesh.b
  pack .draw.mesh.me.aMESH.acm
  pack .display -after .draw
  pack .shrink -after .display 
  pack .shrink.smooth -anchor w
  pack .shrink.cut.mini -before .shrink.cut.aREGION -side bottom -anchor w
  raise .draw
  raise .display
  raise .shrink
  pack .shrink.cut.aROI .shrink.cut.aUNDO .shrink.cut.aINIT \
     -in .shrink.cut.mini -side left  -anchor w
  ## squish
  {.draw.main.aSEND PNT.bu} config -text "SAVE PT"
  {.draw.main.aGOTO PNT.bu} config -text "GOTO PT"
#  .draw.valswitch.aw1.la config -text ColScl -width 6
  .draw.curvswitch.anone.la config -text BackgroundColor -width 15
  #.shrink.cut.aREGION.bu config -text "CUT REGION"
  .shrink.cut.aREGION.bu config -text "CUTAREA"
  .shrink.cut.aPLANE.bu config -text "CUTPLANE"
  .shrink.cut.aLINE.bu config -text "CUTLINE"
  .shrink.cut.aROI.bu config -text "FILL STATS"
  .shrink.cut.aUNDO.bu config -text "UNDO CUT"
  .shrink.cut.aINIT.bu config -text "RECONNECTALL"
  .shrink.smooth.steps.la config -text "" -width 0
  ## auto do fsmask--reads fieldsign suffix but assumes fsmask suffix is .fm
#  .files.fieldsign.br config -command { setfile fs [.files.fieldsign.e get];\
#                                                      read_fieldsign; \
#                     setfile fm [string range [.files.fieldsign.e get] 0 \
#                    [expr [string length [.files.fieldsign.e get]] - 4]].fm; \
#                                                      read_fsmask }
#  .files.fieldsign.bw config -command { setfile fs [.files.fieldsign.e get]; \
#                     setfile fm [string range [.files.fieldsign.e get] 0 \
#                    [expr [string length [.files.fieldsign.e get]] - 4]].fm; \
#                                   testreplace $fm  write_fsmask; \
#                                   testreplace $fs  write_fieldsign }
#  pack .files.val.aPUSH -before .files.val.la -side left
#  .files.val.la config -width 4
  wm geometry . +658+664
  after 500 { wm deiconify . }
  set tksurferinterface csurf
}


proc micro { } {
  global tksurferinterface
  move_window 658 32
  wm withdraw .
  pack forget .mright
  set winlist {.subhead.hemisurf .subhead.script .files.outsurf \
               .files.fieldsign .files.fsmask .files.dipoles .files.decimated \
               .files.val .files.session(*) .draw.curvswitch .draw.valswitch \
               .draw.switch2d .draw.main.aFRAME .draw.main.aFRAME.bu .display \
               .shrink.smooth .shrink.head .shrink.main .shrink.cut.aFLAT \
               .compute .draw.mesh.me.aMESH.acm .shrink.cut.mini \
               .files.val.aPUSH}
  foreach win $winlist { pack forget $win }
  pack .subhead.hemisurf
  pack .subhead.hemisurf.home -side left
  pack .subhead.hemisurf.hemisurf -side left
  pack .draw -after .files
  pack .draw.main.acurv  -after .draw.main.aRESTORE -side left
  pack .draw.main.acurv.ck
  pack .draw.mesh.me     -after .draw.mesh.b   ;# avoid missed REDRAW click
  pack .draw.mesh.sb     -after .draw.mesh.me
  pack .shrink -after .xform
  .subhead.hemisurf.home.e config -width 16
  .subhead.hemisurf.hemisurf.e config -width 18
  .subhead.hemisurf.home.e xview end
  .subhead.hemisurf.hemisurf.e xview end
  {.draw.main.aSEND PNT.bu} config -text "SEND PNT"
  {.draw.main.aGOTO PNT.bu} config -text "GOTO PNT"
  pack .shrink.cut.aROI -before .shrink.cut.aREGION
  pack .shrink.cut.aUNDO -after .shrink.cut.aLINE
  pack .shrink.cut.aINIT -after .shrink.cut.aUNDO
  .shrink.cut.aROI.bu config -text "ROI"
  .shrink.cut.aREGION.bu config -text "CUTREGION"
  .shrink.cut.aPLANE.bu config -text "PLANE"
  .shrink.cut.aLINE.bu config -text "CUTLINE"
  .shrink.cut.aUNDO.bu config -text "UNDO"
  .shrink.cut.aINIT.bu config -text "INIT"
  wm geometry . +658+664
  after 500 { wm deiconify . }
  set tksurferinterface micro
}

############################################################################
wm title . "surfer tools" ;# "tksurfer ($subject--$hemi.$ext)"
wm protocol . WM_DELETE_WINDOW testclose
wm resizable . 0 0

frame .subhead
pack .subhead -side top
  frame .subhead.hemisurf
  pack .subhead.hemisurf -side left
  frame .subhead.script -borderwidth 2 -relief groove
  pack .subhead.script -side right -fill x -expand true

frame .files
pack .files -side top

frame .draw
pack .draw -side top
  frame .draw.main
  pack .draw.main -side top ;# -fill x
  frame .draw.mainsub
  pack .draw.mainsub -side top -fill x
  frame .draw.curvswitch
  pack .draw.curvswitch -side top -fill x
  frame .draw.valswitch
  pack .draw.valswitch -side top -fill x
  frame .draw.mesh
  pack .draw.mesh -side top -fill x
    frame .draw.mesh.me
    pack .draw.mesh.me -side left -fill x
    frame .draw.mesh.r
    pack .draw.mesh.r -side left -fill x
    frame .draw.mesh.g
    pack .draw.mesh.g -side left -fill x
    frame .draw.mesh.b
    pack .draw.mesh.b -side left -fill x
    frame .draw.mesh.sb
    pack .draw.mesh.sb -side left -fill x
  frame .draw.switch2d
  pack .draw.switch2d -side top -fill x

frame .xform
pack .xform -side top
  # rotate
  frame .xform.rotate -borderwidth 2 -relief groove
  pack .xform.rotate -side left
    frame .xform.rotate.top
    pack .xform.rotate.top -side top
    frame .xform.rotate.bot
    pack .xform.rotate.bot -side top
      frame .xform.rotate.bot.la
      pack .xform.rotate.bot.la -side right -anchor center
  # translate
  frame .xform.translate -borderwidth 2 -relief groove
  pack .xform.translate -side left
    frame .xform.translate.top
    pack .xform.translate.top -side top
    frame .xform.translate.bot
    pack .xform.translate.bot -side top
      frame .xform.translate.bot.la
      pack .xform.translate.bot.la -side right -anchor center
  # scale
  frame .xform.scale -borderwidth 2 -relief groove
  pack .xform.scale -side left -anchor n -fill both -expand true
    frame .xform.scale.top
    pack .xform.scale.top -side top
    frame .xform.scale.bot
    pack .xform.scale.bot -side top

frame .display ;# -borderwidth 2 -relief groove
pack .display -side top
  frame .display.le
  pack .display.le -side left
  frame .display.mi
  pack .display.mi -side left
  frame .display.ri
  pack .display.ri -side left

frame .shrink -borderwidth 2 -relief groove
pack .shrink
  frame .shrink.cut
  pack .shrink.cut
    frame .shrink.cut.mini
    pack .shrink.cut.mini -side bottom
  frame .shrink.smooth ;# -borderwidth 2 -relief groove
  pack .shrink.smooth
  frame .shrink.head
  pack .shrink.head -side top
  frame .shrink.main
  pack .shrink.main -side top
    frame .shrink.main.le
    pack .shrink.main.le -side left
    frame .shrink.main.mi
    pack .shrink.main.mi -side left
    frame .shrink.main.ri
    pack .shrink.main.ri -side left

frame .compute ;# -borderwidth 2 -relief groove
pack .compute -side top
  frame .compute.le
  pack .compute.le -side left
  frame .compute.lm
  pack .compute.lm -side left
  frame .compute.rm
  pack .compute.rm -side left
  frame .compute.ri
  pack .compute.ri -side left

frame .mright -borderwidth 2 -relief groove
pack .mright -side right

############################################################################
### title
set f .subhead.hemisurf
edlabval $f "home" $home/$subject n 5 12
.subhead.hemisurf.home.e config -font $ffontb -state disabled
edlabval $f "hemisurf" $hemi.$ext n 5 12
#.subhead.hemisurf.hemisurf.e config -font $ffontb -state disabled

### script box
set f .subhead.script
edlabval $f "script" $script n 7 19
$f.script.e config -textvariable scriptabbrev
setfile script $script       ;# reset scriptabbrev (to show cmdline or default)
buttons $f RUN { setfile script [.subhead.script.script.e get]; \
                                  catch { source $script }; prompt } row 0 6
#was: jot -p -520,0,-987,0
buttons $f EDIT { setfile script [.subhead.script.script.e get]; \
                              exec nedit $script & } row 0 6
buttons $f READENV { source $env(CSURF_DIR)/lib/tcl/readenv.tcl; \
                   setfile script $env(CSURF_DIR)/lib/tcl/readenv.tcl } row 0 6
buttons $f LS \
     { puts ""; puts [exec ls -FC $env(CSURF_DIR)/lib/tcl]; prompt } row 0 6
bind $f.script.e <Return> "$f.aRUN.bu invoke"

### editable file list region
set f .files
if [info exists env(tksurferinterface)] {
    if {$env(tksurferinterface) == "csurf"} {
        #edlabval $f "insurf"     default  r  6 20
        edlabval $f "outsurf"    default   w 6 20
        edlabval $f "curv"       default  rw 6 20
        edlabval $f "patch"      default  rw 6 20
        edlabval $f "label"      default  rw 6 20
        edlabval $f "dipoles"    default   w 6 20
        edlabval $f "decimated"  default   w 6 19
        edlabval $f "fieldsign"  default  rw 6 20
        edlabval $f "fsmask"     default  rw 6 20
        edlabval $f "val"        default  rw 6 20
        edlabval $f "session(*)" $session  s 13 27
        edlabval $f "rgb"        default   w 6 27
    } else {
        #edlabval $f "insurf"     default  r  10 20
        edlabval $f "outsurf"    default   w 10 20
        edlabval $f "curv"       default  rw 10 20
        edlabval $f "patch"      default  rw 10 20
        edlabval $f "label"      default  rw 10 20
        edlabval $f "dipoles"    default   w 10 20
        edlabval $f "decimated"  default   w 10 19
        edlabval $f "fieldsign"  default  rw 10 20
        edlabval $f "fsmask"     default  rw 10 20
        edlabval $f "val"        default  rw 10 20
        edlabval $f "session(*)" $session  s 13 27
        edlabval $f "rgb"        default   w 10 27
    }
} else {
    #edlabval $f "insurf"     default  r  10 20
    edlabval $f "outsurf"    default   w 10 20
    edlabval $f "curv"       default  rw 10 20
    edlabval $f "patch"      default  rw 10 20
    edlabval $f "label"      default  rw 10 20
    edlabval $f "dipoles"    default   w 10 20
    edlabval $f "decimated"  default   w 10 19
    edlabval $f "fieldsign"  default  rw 10 20
    edlabval $f "fsmask"     default  rw 10 20
    edlabval $f "val"        default  rw 10 20
    edlabval $f "session(*)" $session  s 13 27
    edlabval $f "rgb"        default   w 10 27
}

### setfile creates abbrev var from abs var set by surfer at startup
setfile insurf $insurf
setfile outsurf $outsurf
setfile curv $curv
setfile patch $patch
setfile origcoords $origcoords
setfile ellcoords $ellcoords
setfile label $label
setfile dip $dip
setfile dec $dec
setfile fs $fs
setfile fm $fm
setfile val $val
setfile rgb $rgb

### entries use just-created abbrevs
#$f.insurf.e config -textvariable insurfabbrev
$f.outsurf.e config -textvariable outsurfabbrev
$f.curv.e config -textvariable curvabbrev
$f.patch.e config -textvariable patchabbrev
$f.label.e config -textvariable labelabbrev
$f.dipoles.e config -textvariable dipabbrev
$f.decimated.e config -textvariable decabbrev
$f.fieldsign.e config -textvariable fsabbrev
$f.fsmask.e config -textvariable fmabbrev
$f.val.e config -textvariable valabbrev
$f.rgb.e config -textvariable rgbabbrev

### just mini
set f .files
buttons $f.val PUSH { shift_values; set complexvalflag 1 } row 0 4
pack $f.val.aPUSH -before $f.val.la -side left
pack forget $f.val.aPUSH

### one non-standard (std now commented out)
.subhead.hemisurf.hemisurf.e config -textvariable insurfabbrev
bind .subhead.hemisurf.hemisurf.e <Return> \
   { setfile insurf [.subhead.hemisurf.hemisurf.e get]; \
                         read_binary_surf; redrawbutton; \
                         set hemi [file rootname [file tail $insurf]]; \
                         set ext [string trimleft [file tail $insurf] $hemi.]; \
                         wm title . "tksurfer ($subject--$hemi.$ext)" }

### read entry fields at button press (plus 3 non-standard configs)
#$f.insurf.br config -command {setfile insurf [.files.insurf.e get]; \
                                                      read_binary_surf }
$f.outsurf.bw config -command {setfile outsurf [.files.outsurf.e get]; \
                               testreplace $outsurf write_binary_surface }

$f.curv.br config -command { setfile curv [.files.curv.e get]; \
                                           read_binary_curv; set curvflag 1 }
$f.curv.bw config -command { setfile curv [.files.curv.e get]; \
                                 testreplace $curv  write_binary_curv }

$f.patch.br config -command { setfile patch [.files.patch.e get]; \
                                                      read_binary_patch;restore}
$f.patch.bw config -command { setfile patch [.files.patch.e get]; \
                                testreplace $patch  write_binary_patch }

$f.fieldsign.br config -command { setfile fs [.files.fieldsign.e get];\
                                                      read_fieldsign }
$f.fieldsign.bw config -command { setfile fs [.files.fieldsign.e get]; \
                                   testreplace $fs  write_fieldsign }

$f.fsmask.br config -command { setfile fm [.files.fsmask.e get]; \
                                                      read_fsmask }
$f.fsmask.bw config -command { setfile fm [.files.fsmask.e get]; \
                                   testreplace $fm  write_fsmask }

$f.label.br config -command { setfile label [.files.label.e get]; \
                          read_and_color_labeled_vertices $meshr $meshg $meshb }

$f.label.bw config -command { setfile label [.files.label.e get]; \
                                testreplace $label  write_labeled_vertices }

$f.dipoles.bw config -command { setfile dip [.files.dipoles.e get]; \
                                  testreplace $dip  write_binary_dipoles }

$f.decimated.bw config -command { setfile dec [.files.decimated.e get]; \
                                  if {$dip_spacing > 0.999} { \
                                    subsample_dist $dip_spacing \
                                  } else { \
                                    subsample_orient $dip_spacing \
                                  }; \
                                  testreplace $dec  write_binary_decimation }

edlabval $f.decimated "sp" $dip_spacing  n 3 4
$f.decimated.sp.e config -textvariable dip_spacing
trace variable dip_spacing w fixdecname

$f.val.br config -command { setfile val [.files.val.e get]; \
             read_binary_values; fix_nonzero_vals; set overlayflag 1 }
$f.val.bw config -command { setfile val [.files.val.e get]; \
                                  testreplace $val  write_binary_values }

$f.session(*).e config -textvariable session
$f.session(*).bs config -command {setsession [.files.session(*).e get]}

set rel_rgbname [file tail $rgb]  ;# dir resettable; global for uplevel setfile
$f.rgb.bw config -command { setfile rgb [.files.rgb.e get]; \
                                  testreplace $rgb  save_rgb; \
                                  set rel_rgbname [file tail $rgb] }
### main buttons
set f .draw.main
buttons $f REDRAW { redrawbutton } row 0 4
#$f.aREDRAW.bu config -font $ffontbb
buttons $f RESTORE { restore } row 2 4
buttons $f "FRAME" { save_rgb_cmp_frame } row 2 4
#buttons $f "SAVE NUMRGB" { save_rgb_num } row    ;# works too
checks $f "" "curv" curvflag row
.draw.main.acurv.ck config -command { if {$curvflag} {set surfcolor 1} \
                                      else           {set surfcolor 0} }
pack forget .draw.main.acurv.ck
buttons $f "SEND PNT" { find_orig_vertex_coordinates; findsendto; \
   catch { send tkanalyse { goto_point; redraw; set curslice $curslice } }; \
   catch { send $fulltkmedit { goto_point; redraw; unzoomcoords $plane; \
                          sendupdate; sendgoto } }; \
   catch { send $fulltkregister { goto_point; set updateflag TRUE; \
                             unzoomcoords $plane } } } row 1 4
buttons $f "GOTO PNT" { select_orig_vertex_coordinates } row 1 4

### for blink doublebuffer (now in tmp proc)
#set f .draw.mainsub
#buttons $f REDRAW2 { redraw_second } row 0 4
#$f.aREDRAW2.bu config -font $ffontbb
#buttons $f RESTORE2 { restore_second } row 2 4
#label $f.la -text " " -font $ffont  ;# space
#pack $f.la -side left
#buttons $f SWAPBUFFERS { swap_buffers } row 2 4
#bind $f.aSWAPBUFFERS.bu <ButtonPress-1> { set blinkflag TRUE }
#bind $f.aSWAPBUFFERS.bu <ButtonRelease-1> { set blinkflag FALSE; \
                                      set blinkdelay $initdelay }
### rotate: label
set f .xform.rotate.top
label $f.la -text "ROTATE (deg)" -font $ffontb -pady 0
pack $f.la -side top
### horiz
scale $f.y -from 180 -to -180 -length $sclenx -variable yrot \
    -orient horizontal -tickinterval 90 -showvalue false -font $sfont -width 11
pack $f.y -side top
#bind $f.y <ButtonRelease-1> { redrawbutton }
### vertical
set f .xform.rotate.bot
scale $f.x -from 180 -to -180 -length $scleny -variable xrot \
    -orient vertical -tickinterval 90 -showvalue false -font $sfont -width 11
pack $f.x -side left ;# -anchor e
#bind $f.x <ButtonRelease-1> { redrawbutton }
### entries
set f .xform.rotate.bot.la
edlabval $f "yrot" 0 n 5 4
edlabval $f "xrot" 0 n 5 4
label $f.space -text "" -font $sfont  ;# a little space
pack $f.space -side top
edlabval $f "zrot" 0 n 5 4
$f.yrot.e config -textvariable yrot -font $sfont
$f.xrot.e config -textvariable xrot -font $sfont
$f.zrot.e config -textvariable zrot -font $sfont
$f.yrot.la config -font $sfont
$f.xrot.la config -font $sfont
$f.zrot.la config -font $sfont
bind $f.xrot.e <Return> { redrawbutton }
bind $f.yrot.e <Return> { redrawbutton }
bind $f.zrot.e <Return> { redrawbutton }

### translate: label
set f .xform.translate.top
label $f.la -text "TRANSLATE (%)" -font $ffontb -pady 0
pack $f.la -side top
### horiz
scale $f.x -from -100 -to 100 -length $sclenx -variable xtrans \
   -orient horizontal -tickinterval 50 -showvalue false -width 11 -font $sfont
pack $f.x -side top
#bind $f.x <ButtonRelease-1> { redrawbutton }
### vertical
set f .xform.translate.bot
scale $f.y -from 100 -to -100 -length $scleny -variable ytrans \
    -orient vertical -tickinterval 50 -showvalue false -width 11 -font $sfont
pack $f.y -side left
#bind $f.y <ButtonRelease-1> { redrawbutton }
### entries
set f .xform.translate.bot.la
edlabval $f "x" 0 n 2 4
edlabval $f "y" 0 n 2 4
$f.x.e config -textvariable xtrans -font $sfont
$f.y.e config -textvariable ytrans -font $sfont
$f.x.la config -font $sfont
$f.y.la config -font $sfont
bind $f.x.e <Return> { redrawbutton }
bind $f.y.e <Return> { redrawbutton }

### scale: label
set f .xform.scale.top
label $f.la -text "SCALE (%)" -font $ffontb -pady 0
pack $f.la -side top
### vertical
scale $f.xyz -from 300 -to 1 -length [expr $scleny+15] -variable scalepercent \
   -orient vertical -tickinterval 100 -showvalue false -font $sfont -width 11
pack $f.xyz -side top -expand false
#bind $f.xyz <ButtonRelease-1> { redrawbutton }
### entry
set f .xform.scale.bot
edlabval $f "%" 0 n 2 4
$f.%.e config -textvariable scalepercent -font $sfont
$f.%.la config -font $sfont
bind $f.%.e <Return> { redrawbutton }

### BGcolor radiorow
set f .draw.curvswitch
radios $f "BgSurfCol" "none" surfcolor  0 10 row
radios $f ""          "curv" surfcolor  1 10 row
radios $f ""          "area" surfcolor  2 10 row
radios $f ""          "shear" surfcolor 3 10 row

### colscale radiorow
set f .draw.valswitch
# complex
radios $f "Scl" "w1"  colscale 0 4 row;  $f.aw1.ra config -padx 0
radios $f ""    "w2"  colscale 8 4 row;  $f.aw2.ra config -padx 0
radios $f ""    "bi"  colscale 4 4 row;  $f.abi.ra config -padx 0
# positive
radios $f ""    "ht"  colscale 1 4 row;  $f.aht.ra config -padx 0
#radios $f ""   "br"  colscale 2 4 row
#radios $f ""   "bgr" colscale 3 4 row
#radios $f ""   "gr"  colscale 5 4 row
# signed
radios $f ""    "BRw" colscale 6 4 row;  $f.aBRw.ra config -padx 0
radios $f ""    "GR"  colscale 7 4 row;  $f.aGR.ra config -padx 0
radios $f ""    "fs"  colscale 9 4 row;  $f.afs.ra config -padx 0
trace variable colscale w fixcolscale

### mesh, scalebars
set f .draw.mesh
buttons $f.me "MESH" \
  {if {$verticesflag} { set verticesflag 0; \
    .draw.mesh.me.aMESH.bu config -relief raised } \
  else {set verticesflag 1; \
    .draw.mesh.me.aMESH.bu config -relief sunken}; \
  redraw } row 0 7
edlabval $f.r "r" $meshr n 2 3
edlabval $f.g "g" $meshg n 2 3
edlabval $f.b "b" $meshb n 2 3
$f.r.r.e config -textvariable meshr
$f.g.g.e config -textvariable meshg
$f.b.b.e config -textvariable meshb
label $f.sb.la -text " " -font $ffont  ;# space
pack $f.sb.la -side left
checks $f.sb ""  "cm"        scalebarflag row
checks $f.sb ""  "colbar"  colscalebarflag row
buttons .draw.mesh.me.aMESH "cm" \
  {if {$scalebarflag} { set scalebarflag 0; set colscalebarflag 0; \
    .draw.mesh.me.aMESH.acm.bu config -relief raised } \
  else {set scalebarflag 1; set colscalebarflag 1; \
    .draw.mesh.me.aMESH.acm.bu config -relief sunken}; \
  redraw } row 0 2
pack forget .draw.mesh.me.aMESH.acm

### 2d vectors checkrow
set f .draw.switch2d
checks $f "2d Vectors" "shear" shearvecflag row
checks $f ""           "norm" normvecflag row
checks $f ""           "move" movevecflag row
buttons $f "CVAVG" { set cmid $dipavg } row 1

### display panel left
set f .display.le
checks $f "" "truncate stats" truncphaseflag col
checks $f "" "invert stats" invphaseflag col
checks $f "" "revphase" revphaseflag col
if [info exists env(tksurferinterface)] {
    if {$env(tksurferinterface) == "csurf"} {
    checks $f "" "scale bar" scalebarflag col
    checks $f "" "color bar" colscalebarflag col
    }
} else {
  checks $f "" "scale bar" scalebarflag col
  checks $f "" "color bar" colscalebarflag col
}


checks $f "" "overlay" overlayflag col
checks $f "" "complexval" complexvalflag col
if [info exists env(tksurferinterface)] {
    if {$env(tksurferinterface) != "csurf"} {
        edlabval $f "fthresh" $fthresh n 8 4
        $f.fthresh.e config -textvariable fthresh
    }
}
### display panel middle
set f .display.mi
if [info exists env(tksurferinterface)] {
    if {$env(tksurferinterface) == "csurf"} {
        edlabval $f "fthresh" $fthresh n 12 4
        $f.fthresh.e config -textvariable fthresh
    }
} else {
    edlabval $f "fthresh" $fthresh n 12 4
    $f.fthresh.e config -textvariable fthresh
}

#edlabval $f "fcurv"   $fcurv   n 12 4
edlabval $f "fmid"    $fmid    n 12 4
edlabval $f "fslope"  $fslope  n 12 4
edlabval $f "cslope"  $cslope  n 12 4
edlabval $f "cmid"    $cmid    n 12 4
if [info exists env(tksurferinterface)] {
    if {$env(tksurferinterface) != "csurf"} {
        edlabval $f "anglecycles" $angle_cycles n 12 4
        edlabval $f "angleoffset" $angle_offset n 12 4
        $f.anglecycles.e config -textvariable angle_cycles
        $f.angleoffset.e config -textvariable angle_offset
    }
}
$f.fslope.e config -textvariable fslope 
#$f.fcurv.e config -textvariable fcurv
$f.fmid.e config -textvariable fmid
$f.cslope.e config -textvariable cslope 
$f.cmid.e config -textvariable cmid
### display panel right
set f .display.ri
edlabval $f "light0" $light0 n 9 4
edlabval $f "light1" $light1 n 9 4
edlabval $f "light2" $light2 n 9 4
edlabval $f "light3" $light3 n 9 4
edlabval $f "offset" $offset n 9 4
edlabval $f "blufact" $blufact n 9 4
$f.light0.e config -textvariable light0
$f.light1.e config -textvariable light1
$f.light2.e config -textvariable light2
$f.light3.e config -textvariable light3
$f.offset.e config -textvariable offset
$f.blufact.e config -textvariable blufact

### cut surface buttons
set f .shrink.cut
buttons $f "ROI" { floodfill_marked_patch 1; redrawbutton } row 1 2
buttons $f "REGION" { cut_line 1; redrawbutton } row 1 2
buttons $f "FILLAREA" { floodfill_marked_patch 0; redrawbutton } row 1 2
buttons $f "PLANE" { cut_plane; redrawbutton } row 1 2
buttons $f "LINE" { cut_line 0; redrawbutton } row 1 2
buttons $f "UNDO" { restore_ripflags 1; redrawbutton } row 1 2
buttons $f "INIT" { restore_ripflags 2; redrawbutton } row 1 2
#buttons $f "INIT" { restore_ripflags 2; read_binary_surf; read_binary_curv; \
                    set flag2d 0; restore; redrawbutton } row 1 2
buttons $f "FLAT" { flatten; set flag2d 1; restore; redrawbutton } row 1 2
### shrink smooth row
set f .shrink.smooth
buttons $f SMOOTH { smooth $smoothsteps $smoothtype } row 1
edlabval $f "steps" 1 n 6 4
$f.steps.e config -textvariable smoothsteps
radios $f.steps "" "curv" smoothtype curv     12 row
radios $f.steps "" "val" smoothtype val       12 row
radios $f.steps "" "sparse" smoothtype sparse 12 row
### shrink panel head 
set f .shrink.head
buttons $f SHRINK { shrink $shrinksteps } row 1
edlabval $f "steps" 1 n 6 4 
$f.steps.e config -textvariable shrinksteps 
label $f.steps.la2 -text " MRIforce:" -font $ffont
pack $f.steps.la2 -side left 
buttons $f.steps UP { mriforce up } row 1 6
buttons $f.steps DOWN { mriforce down } row 1 6
checkbutton $f.steps.ck -variable MRIflag -font $ffont
pack $f.steps.ck -side left 
### shrink panel main left
set f .shrink.main.le
edlabval $f "wt"  $wt  n 6 5
edlabval $f "wa"  $wa  n 6 5
edlabval $f "ws"  $ws  n 6 5
edlabval $f "wn"  $wn  n 6 5
edlabval $f "wsh" $wsh n 6 5
edlabval $f "wbn" $wbn n 6 5
#edlabval $f "wc"  $wc n 6 5
$f.wt.e  config -textvariable wt
$f.wa.e  config -textvariable wa
$f.ws.e  config -textvariable ws
$f.wn.e  config -textvariable wn
$f.wsh.e config -textvariable wsh
$f.wbn.e config -textvariable wbn
#$f.wc.e  config -textvariable wc
### shrink panel main middle
set f .shrink.main.mi
#checks $f "" "ncthrflag" ncthreshflag  col
edlabval $f "ncthr"   $ncthresh n 7 5 
checkbutton $f.ncthr.ck -variable ncthreshflag -font $ffont
pack $f.ncthr.ck -side left
buttons $f "NORM AREA" { normalize_area } col 1
checks $f "" "area"      areaflag      col
#checks $f "" "expand"   expandflag    col
checks $f "" "momentum"  momentumflag  col
edlabval $f "update"  $update   n 8 5
edlabval $f "decay"   $decay    n 8 5
$f.ncthr.e config -textvariable ncthresh
$f.update.e config -textvariable update
$f.decay.e config -textvariable decay
### center override
#$f.ancthrflag.ck config -anchor center
$f.aarea.ck config -anchor center
$f.amomentum.ck config -anchor center

### shrink panel main right
set f .shrink.main.ri
#checks $f "" "MRI data on"   MRIflag  col
#{.shrink.main.ri.aMRI data on.ck} config -anchor center
edlabval $f  "mstrength"    $mstrength    n 10 5
edlabval $f  "mslope"       $mslope       n 10 5
edlabval $f  "mmid"         $mmid         n 10 5
edlabval $f  "whitemid"     $whitemid     n 10 5
edlabval $f  "graymid"      $graymid      n 10 5
edlabval $f  "cthk"         $cthk         n 10 5
$f.mstrength.e config -textvariable mstrength
$f.mslope.e config -textvariable mslope
$f.mmid.e config -textvariable mmid 
$f.whitemid.e config -textvariable whitemid
$f.graymid.e config -textvariable graymid
$f.cthk.e config -textvariable cthk

### compute panel left
buttons .compute.le "COMP CURV" { compute_curvature } col 0
buttons .compute.le "CLEAR CURV"   { clear_curvature }   col 0 5
### compute panel left middle
checks .compute.lm "" "sulc sum" sulcflag col
buttons .compute.lm "COMP THK" { compute_cortical_thickness } col 0
### compute panel right middle
buttons .compute.rm "COMP FS"  { compute_fieldsign }  col 0
buttons .compute.rm "COMP CMF" { compute_CMF }       col 0
### compute panel left
buttons .compute.ri "pushval" { shift_values; set complexvalflag 1 }   col 0
buttons .compute.ri "swapcmplx"  { swap_values}    col 0 5
#buttons .compute.ri "big swap"  { big_swap_values } col 0

############################################################################
### shortcut key bindings
bind . <Alt-a> { set ncthreshflag [expr !$ncthreshflag] }
bind . <Alt-A> { set areaflag [expr !$areaflag] }
bind . <Alt-B> { .files.patch.bw invoke }
bind . <Alt-d> { set shrinksteps 100; update idletasks; shrink $shrinksteps }
bind . <Alt-D> { set shrinksteps 1000; update idletasks; shrink $shrinksteps }
bind . <Alt-c> { compute_curvature }
bind . <Alt-C> { clear_curvature }
bind . <Alt-f> { ".draw.main.aSEND PNT.bu" invoke }
bind . <Alt-F> { set momentumflag [expr !$momentumflag] }
bind . <Alt-G> { testreplace $area write_binary_areas }
bind . <Alt-h> { help; prompt }
bind . <Alt-K> { setfile curv [.files.curv.e get]; read_binary_curv }
bind . <Alt-l> { set smoothsteps 5; set smoothtype curv; update idletasks; \
                                            smooth $smoothsteps $smoothtype }
bind . <Alt-m> { set smoothsteps 5; set smoothtype sparse; update idletasks; \
                                            smooth $smoothsteps $smoothtype }
bind . <Alt-Key-M> { set MRIflag [expr !$MRIflag] }  ;# need "Key" just for M??
bind . <Alt-N> { .files.dipoles.bw invoke; .files.decimated.bw invoke }
bind . <Alt-o> { set overlayflag [expr !$overlayflag] }
bind . <Alt-Q> { .files.val.bw invoke }
bind . <Alt-r> { redrawbutton }
#bind . <Alt-R> { setfile sulc [.files.sulc.e get]; read_binary_sulc }
bind . <Alt-s> { set shrinksteps 1; update idletasks; shrink $shrinksteps }
bind . <Alt-S> { set shrinksteps 10; update idletasks; shrink $shrinksteps }
bind . <Alt-U> { normalize_area }
bind . <Alt-v> { set verticesflag [expr !$verticesflag] }
bind . <Alt-w> { .files.outsurf.bw invoke }
bind . <Alt-T> { save_rgb }
bind . <Alt-Key-0> { restore_zero_position }
bind . <Alt-Key-1> { restore_initial_position }
bind . <Alt-Key-2> { make_lateral_view }
bind . <Alt-slash>    { mriforce down }
bind . <Alt-asterisk> { mriforce up }
bind . <Alt-Up>    { set xrot [expr $xrot+18.0]}
bind . <Alt-Down>  { set xrot [expr $xrot-18.0]}
bind . <Alt-Right> { set yrot [expr $yrot-18.0]}
bind . <Alt-Left>  { set yrot [expr $yrot+18.0]}
bind . <Alt-Shift-Left>  { set yrot [expr $yrot+180.0]}
bind . <Alt-Shift-Right>  { set yrot [expr $yrot-180.0]}
bind . <Alt-KP_Next>  { set zrot [expr $zrot+10.0]}
bind . <Alt-KP_Prior> { set zrot [expr $zrot-10.0]}
bind . <Alt-KP_Right>   { set xtrans [expr $xtrans+10.0]}
bind . <Alt-KP_Left>    { set xtrans [expr $xtrans-10.0]}
bind . <Alt-KP_Up>      { set ytrans [expr $ytrans+10.0]}
bind . <Alt-KP_Down>    { set ytrans [expr $ytrans-10.0]}
bind . <Alt-parenleft>  { set cthk [expr $cthk-0.1]}
bind . <Alt-parenright> { set cthk [expr $cthk+0.1]}
bind . <Alt-braceleft>    { set scalepercent [expr $scalepercent/1.25]}
bind . <Alt-braceright>   { set scalepercent [expr $scalepercent*1.25]}
bind . <Alt-bracketleft>  { set scalepercent [expr $scalepercent/1.05]}
bind . <Alt-bracketright> { set scalepercent [expr $scalepercent*1.05]}
bind . <Control-F1> { micro }
bind . <Control-F2> { mini  }
bind . <Control-F3> { mini+ }
bind . <Control-F4> { macro }
bind . <Control-F5> { csurf }

############################################################################
### right-click help
bind .draw.main.aREDRAW.bu <ButtonRelease-3> { helpwin redraw }
bind .draw.main.aRESTORE.bu <ButtonRelease-3> { helpwin restore }
bind .files.decimated.sp.e  <ButtonRelease-3> { helpwin dip_spacing }
bind ".draw.main.aSEND PNT.bu" <ButtonRelease-3> { helpwin save_pnt }
bind ".draw.main.aGOTO PNT.bu" <ButtonRelease-3> { helpwin goto_pnt }
bind .display.le.arevphase.ck <ButtonRelease-3> { helpwin revphase }
bind ".display.le.ainvert stats.ck" <ButtonRelease-3> { helpwin invphase }
bind ".display.le.atruncate stats.ck" <ButtonRelease-3> { helpwin truncphase }
bind .display.ri.light0.e <ButtonRelease-3> { helpwin light0 }
bind .display.ri.light1.e <ButtonRelease-3> { helpwin light1 }
bind .display.ri.light2.e <ButtonRelease-3> { helpwin light2 }
bind .display.ri.light3.e <ButtonRelease-3> { helpwin light3 }
bind .shrink.cut.aREGION.bu <ButtonRelease-3> { helpwin cutregion }
bind .shrink.cut.aFILLAREA.bu <ButtonRelease-3> { helpwin fill }
bind .shrink.cut.aPLANE.bu <ButtonRelease-3> { helpwin cutplane }
bind .shrink.cut.aLINE.bu <ButtonRelease-3> { helpwin cutline }
bind .shrink.cut.aFLAT.bu <ButtonRelease-3> { helpwin flatten }
bind .shrink.cut.aUNDO.bu <ButtonRelease-3> { helpwin undo }
bind .shrink.cut.aINIT.bu <ButtonRelease-3> { helpwin init }
bind .shrink.cut.aROI.bu <ButtonRelease-3> { helpwin roi }
bind .subhead.hemisurf.home.e <ButtonRelease-3> { helpwin home }
bind .subhead.script.script.e <ButtonRelease-3> { helpwin script }
bind .subhead.hemisurf.hemisurf.e <ButtonRelease-3> { helpwin hemisurf }
bind .subhead.script.aRUN.bu <ButtonRelease-3> { helpwin run }
bind .subhead.script.aEDIT.bu <ButtonRelease-3> { helpwin edit }
bind .subhead.script.aREADENV.bu <ButtonRelease-3> { helpwin readenv }
bind .subhead.script.aLS.bu <ButtonRelease-3> { helpwin ls }
bind .files.outsurf.e <ButtonRelease-3> { helpwin outsurf }
bind .files.outsurf.bw <ButtonRelease-3> { helpwin outsurf_write }
bind .files.curv.e <ButtonRelease-3> { helpwin curv }
bind .files.curv.br <ButtonRelease-3> { helpwin curv_read }
bind .files.curv.bw <ButtonRelease-3> { helpwin curv_write }
bind .files.patch.la <ButtonRelease-3> { helpwin cutinstructions }
bind .files.val.e <ButtonRelease-3> { helpwin val }
bind .files.val.br <ButtonRelease-3> { helpwin val_read }
bind .files.val.bw <ButtonRelease-3> { helpwin val_write }
bind .files.rgb.e <ButtonRelease-3> { helpwin rgb }
bind .files.rgb.bw <ButtonRelease-3> { helpwin rgb_write }
bind .files.session(*).e <ButtonRelease-3> { helpwin session }
bind .files.session(*).bs <ButtonRelease-3> { helpwin session_set }
bind .draw.curvswitch.anone.la <ButtonRelease-3> { helpwin curv_bgcol }
bind .draw.curvswitch.anone.ra <ButtonRelease-3> { helpwin curv_none }
bind .draw.curvswitch.acurv.ra <ButtonRelease-3> { helpwin curv_curv }
bind .draw.valswitch.aw1.la <ButtonRelease-3> { helpwin val_scale }
bind .draw.valswitch.aw1.ra <ButtonRelease-3> { helpwin val_w1 }
bind .draw.valswitch.aw2.ra <ButtonRelease-3> { helpwin val_w2 }
bind .draw.valswitch.abi.ra <ButtonRelease-3> { helpwin val_bi }
bind .draw.valswitch.aht.ra <ButtonRelease-3> { helpwin val_ht }
bind .draw.valswitch.aBRw.ra <ButtonRelease-3> { helpwin val_BRw }
bind .draw.valswitch.aGR.ra <ButtonRelease-3> { helpwin val_GR }
bind .draw.valswitch.afs.ra <ButtonRelease-3> { helpwin val_fs }

### help todo
bind .files.patch.e <ButtonRelease-3> { helpwin patch }
bind .files.patch.br <ButtonRelease-3> { helpwin patch_read }
bind .files.patch.bw <ButtonRelease-3> { helpwin patch_write }
bind .files.label.e <ButtonRelease-3> { helpwin label }
bind .files.label.bw <ButtonRelease-3> { helpwin label_write }
bind .files.dipoles.e <ButtonRelease-3> { helpwin dipoles }
bind .files.dipoles.bw <ButtonRelease-3> { helpwin dipoles_write }
bind .files.decimated.e <ButtonRelease-3> { helpwin decimated }
bind .files.decimated.bw <ButtonRelease-3> { helpwin decimated_write }
bind .files.fieldsign.e <ButtonRelease-3> { helpwin fieldsign }
bind .files.fieldsign.br <ButtonRelease-3> { helpwin fieldsign_read }
bind .files.fieldsign.bw <ButtonRelease-3> { helpwin fieldsign_write }
bind .files.fsmask.e <ButtonRelease-3> { helpwin fsmask }
bind .files.fsmask.br <ButtonRelease-3> { helpwin fsmask_read}
bind .files.fsmask.bw <ButtonRelease-3> { helpwin fsmask_write }
bind .draw.curvswitch.aarea.ra <ButtonRelease-3> { helpwin curv_area }
bind .draw.curvswitch.ashear.ra <ButtonRelease-3> { helpwin curv_shear }
bind .draw.main.aFRAME.bu <ButtonRelease-3> { helpwin frame }
bind .draw.mesh.me.aMESH.bu <ButtonRelease-3> { helpwin mesh }
bind .draw.mesh.r.r.e <ButtonRelease-3> { helpwin mesh_r }
bind .draw.mesh.g.g.e <ButtonRelease-3> { helpwin mesh_g }
bind .draw.mesh.b.b.e <ButtonRelease-3> { helpwin mesh_b }
bind .draw.mesh.sb.acm.ck <ButtonRelease-3> { helpwin cm }
bind .draw.mesh.sb.acolbar.ck <ButtonRelease-3> { helpwin colbar }
bind .draw.switch2d.ashear.ck <ButtonRelease-3> { helpwin 2d_shear }
bind .draw.switch2d.anorm.ck <ButtonRelease-3> { helpwin 2d_norm }
bind .draw.switch2d.amove.ck <ButtonRelease-3> { helpwin 2d_move }
bind .draw.switch2d.aCVAVG.bu <ButtonRelease-3> { helpwin cvavg }
bind .display.le.aoverlay.ck <ButtonRelease-3> { helpwin overlay }
bind .display.le.acomplexval.ck <ButtonRelease-3> { helpwin complexval }

if [info exists env(tksurferinterface)] {
    if {$env(tksurferinterface) != "csurf"} {
        bind .display.le.fthresh.e <ButtonRelease-3> { helpwin fthresh }
    } else {
        bind .display.mi.fthresh.e <ButtonRelease-3> { helpwin fthresh }
    }
} else {
    bind .display.mi.fthresh.e <ButtonRelease-3> { helpwin fthresh }
}


bind .display.mi.fslope.e <ButtonRelease-3> { helpwin fslope }
bind .display.mi.fmid.e <ButtonRelease-3> { helpwin fmid }
bind .display.mi.cslope.e <ButtonRelease-3> { helpwin cslope }
bind .display.mi.cmid.e <ButtonRelease-3> { helpwin cmid }
if [info exists env(tksurferinterface)] {
    if {$env(tksurferinterface) != "csurf"} {
      bind .display.mi.anglecycles.e <ButtonRelease-3> { helpwin anglecycles }
      bind .display.mi.angleoffset.e <ButtonRelease-3> { helpwin angleoffset }
      bind .display.mi.angleoffset.e <ButtonRelease-3> { helpwin angleoffset }
  }
}
bind .display.ri.offset.e <ButtonRelease-3> { helpwin offset }
bind .display.ri.blufact.e <ButtonRelease-3> { helpwin blufact }
bind .shrink.smooth.aSMOOTH.bu <ButtonRelease-3> { helpwin smooth }
bind .shrink.smooth.steps.e <ButtonRelease-3> { helpwin smooth_steps }
bind .shrink.smooth.steps.acurv.ra <ButtonRelease-3> { helpwin smooth_curv }
bind .shrink.smooth.steps.aval.ra <ButtonRelease-3> { helpwin smooth_val }
bind .shrink.smooth.steps.asparse.ra <ButtonRelease-3> { helpwin smooth_sparse}
bind .shrink.head.aSHRINK.bu <ButtonRelease-3> { helpwin shrink }
bind .shrink.head.steps.e <ButtonRelease-3> { helpwin shrink_steps }
bind .shrink.head.steps.aUP.bu <ButtonRelease-3> { helpwin shrink_up }
bind .shrink.head.steps.aDOWN.bu <ButtonRelease-3> { helpwin shrink_down }
bind .shrink.head.steps.ck <ButtonRelease-3> { helpwin shrink_steps_ck }
bind .shrink.main.le.wt.e <ButtonRelease-3> { helpwin shrink_wt }
bind .shrink.main.le.wa.e <ButtonRelease-3> { helpwin shrink_wa }
bind .shrink.main.le.ws.e <ButtonRelease-3> { helpwin shrink_ws }
bind .shrink.main.le.wn.e <ButtonRelease-3> { helpwin shrink_wn }
bind .shrink.main.le.wsh.e <ButtonRelease-3> { helpwin shrink_wsh }
bind .shrink.main.le.wbn.e <ButtonRelease-3> { helpwin shrink_wbn }
bind .shrink.main.mi.ncthr.e <ButtonRelease-3> { helpwin shrink_ncthr }
bind .shrink.main.mi.ncthr.ck <ButtonRelease-3> { helpwin shrink_ncthr_ck }
bind ".shrink.main.mi.aNORM AREA.bu" <ButtonRelease-3> { helpwin norm_area }
bind .shrink.main.mi.aarea.ck <ButtonRelease-3> { helpwin shrink_area }
bind .shrink.main.mi.amomentum.ck <ButtonRelease-3> { helpwin shrink_momentum }
bind .shrink.main.mi.update.e <ButtonRelease-3> { helpwin shrink_update }
bind .shrink.main.mi.decay.e <ButtonRelease-3> { helpwin shrink_decay }
bind .shrink.main.ri.mstrength.e <ButtonRelease-3> { helpwin shrink_mstrength }
bind .shrink.main.ri.mslope.e <ButtonRelease-3> { helpwin shrink_mslope }
bind .shrink.main.ri.mmid.e <ButtonRelease-3> { helpwin shrink_mmid }
bind .shrink.main.ri.whitemid.e <ButtonRelease-3> { helpwin shrink_whitemid }
bind .shrink.main.ri.graymid.e <ButtonRelease-3> { helpwin shrink_graymid }
bind .shrink.main.ri.cthk.e <ButtonRelease-3> { helpwin shrink_cthk }
bind ".compute.le.aCOMP CURV.bu" <ButtonRelease-3> { helpwin comp_curv }
bind ".compute.le.aCLEAR CURV.bu" <ButtonRelease-3> { helpwin clear_curv }
bind ".compute.lm.asulc sum.ck" <ButtonRelease-3> { helpwin sulc_sum }
bind ".compute.lm.aCOMP THK.bu" <ButtonRelease-3> { helpwin comp_thk }
bind ".compute.rm.aCOMP FS.bu" <ButtonRelease-3> { helpwin comp_fs }
bind ".compute.rm.aCOMP CMF.bu" <ButtonRelease-3> { helpwin comp_cmf }
bind .compute.ri.apushval.bu <ButtonRelease-3> { helpwin pushval }
bind .compute.ri.aswapcmplx.bu <ButtonRelease-3> { helpwin swapcmplx }

############################################################################
puts "tksurfer.tcl: startup done"

# N.B.: w/command line script, these execute before gl window up
fixcolors
if [info exists env(tksurferinterface)] {
  if {$env(tksurferinterface) == "macro"} { macro; update }
  if {$env(tksurferinterface) == "mini+"} { mini+; update }
  if {$env(tksurferinterface) == "mini"}  { mini;  update }
  if {$env(tksurferinterface) == "micro"} { micro; update }
  if {$env(tksurferinterface) == "csurf"} { csurf; update }
} else {
  csurf
  puts "tksurfer.tcl: default csurf interface (change: csurf,macro,mini+,mini,micro)"
  puts "tksurfer.tcl: or: setenv tksurferinterface {csurf,macro,mini+,mini,micro}"
}

### extras: fonts tested
#set ffont -b&h-lucidatypewriter-medium-*-*-*-10-*-*-*-*-*-*-*
#set ffontb -b&h-lucidatypewriter-bold-*-*-*-10-*-*-*-*-*-*-*
#set ffont -*-helvetica-med-*-*-*-10-*-*-*-*-*-*-*
#set ffontb -*-helvetica-bold-*-*-*-10-*-*-*-*-*-*-*
#set ffont -b&h-lucidabright-medium-i-normal--10-100-75-75-p-57-iso8859-1
#set ffontb -b&h-lucidabright-demibold-i-normal--10-100-75-75-p-59-iso8859-1
#set ffont -sgi-screen-medium-r-normal-*-10-100-72-72-m-60-iso8859-1
#set ffontb -sgi-screen-bold-r-normal-*-10-100-72-72-m-70-iso8859-1
#set ffont -sgi-screen-medium-r-normal-*-11-110-72-72-m-*-ascii
#set ffontb -sgi-screen-bold-r-normal-*-11-110-72-72-m-*-ascii
#set ffont -sgi-screen-medium-r-normal-*-11-110-72-72-m-70-iso8859-1
#set ffontb -sgi-screen-bold-r-normal-*-11-110-72-72-m-80-iso8859-1
#set ffont -sgi-screen-medium-r-normal-*-10-100-72-72-m-*-ascii
#set ffontb -sgi-screen-bold-r-normal-*-10-100-72-72-m-*-ascii

