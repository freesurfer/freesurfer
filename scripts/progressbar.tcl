#
# Progressbar Widget written in pure tcl
#
# @(#)progressbar.tcl v1.3 00/04/28 (c) 2000 Alexander Schoepe
#
# Progressbar is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# Progressbar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Progressbar; see the file COPYING.  If not, write to
# the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Special thanks to Bryan Oakley his mclistbox was a very good example
# how to build a widget - GREAT SOURCE!
# Futhermore special thanks to Joerg Mehring who ask me just miniutes
# before releasing "Is there a variable parameter?". So the release had
# to wait some more time and a variable parameter has been added.
#
# Always remember there is NO SUPPORT! eMail: tcl@sesam.com
#
# ###############################  USAGE  #################################
#
#    NAME
#       progressbar - Create and manipulate progressbar widgets
#    SYNOPSIS
#       progressbar pathName ?options?
#    STANDARD OPTIONS
#       -borderwidth or -bd, borderWidth, BorderWidth
#       -relief, relief, Relief
#    WIDGET-SPECIFIC OPTIONS
#       -background, background, Background
#       -width, width, Width
#       -color, color, Color
#       -percent, percent, Percent
#       -shape, shape, Shape
#       -variable, variable, Variable
#       -width, width, Width
#    WIDGET COMMAND
#       pathName cget option
#       pathName configure ?option? ?value option value ...?
#
# ##############################  XSAMPLE  ################################
#
#     package require progressbar 1.3
#
#     pack [set w [::progressbar::progressbar .pb]]
#     for {set percent 0} {$percent <= 100} {incr percent} {
#       $w configure -percent $percent
#       update
#       after 100
#     }
#     destroy $w
#
#     set percent 0
#     pack [::progressbar::progressbar .pb -variable percent]
#     for {} {$percent <= 100} {incr percent} {
#       update
#       after 100
#     }
#     destroy $w
#
# ############################  NO WARRANTY  ##############################
#
# BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
# FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW.  EXCEPT WHEN
# OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
# PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED
# OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE ENTIRE RISK AS
# TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU.  SHOULD THE
# PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING,
# REPAIR OR CORRECTION.
#
# #######################   KEINE GEWÄHRLEISTUNG  #########################
#
# Da das Programm ohne jegliche Kosten lizenziert wird, besteht
# keinerlei Gewährleistung für das Programm, soweit dies gesetzlich zulässig
# ist. Sofern nicht anderweitig schriftlich bestätigt, stellen die
# Copyright-Inhaber und/oder Dritte das Programm so zur Verfügung, "wie es
# ist", ohne irgendeine Gewährleistung, weder ausdrücklich noch implizit,
# einschließlich - aber nicht begrenzt auf - Marktreife oder Verwendbarkeit
# für einen bestimmten Zweck. Das volle Risiko bezüglich Qualität und
# Leistungsfähigkeit des Programms liegt bei Ihnen. Sollte sich das Programm
# als fehlerhaft herausstellen, liegen die Kosten für notwendigen Service,
# Reparatur oder Korrektur bei Ihnen.
#
# #######################  LIMITATION DE GARANTIE  ########################
#
# Parce que l'utilisation de ce Programme est libre et gratuite, aucune
# garantie n'est fournie, comme le permet la loi. Sauf mention écrite,
# les détenteurs du copyright et/ou les tiers fournissent le Programme en
# l'état, sans aucune sorte de garantie explicite ou implicite, y compris
# les garanties de commercialisation ou d'adaptation dans un but
# particulier. Vous assumez tous les risques quant à la qualité et aux
# effets du Programme. Si le Programme est défectueux, Vous assumez le
# coût de tous les services, corrections ou réparations nécessaires.
#
# ##############################  HISTORY  ################################
#
#  1.1  first release
#  1.2  bug fixes removing variable trace
#  1.3  another bug fix removing variable trace
#


package require Tk 8.0
package provide progressbar 1.3


namespace eval ::progressbar {
  # this is the public interface
  namespace export progressbar

  # these contain references to available options
  variable widgetOptions

  # these contain references to available commands
  variable widgetCommands

  # these contain references to available options for shape option
  variable widgetShapes

  # these contain references to global variables
  variable widgetGlobals

  set widgetGlobals(debug) 0
}


# ::progressbar::Init --
#
#     Initialize the global (well, namespace) variables. This should
#     only be called once, immediately prior to creating the first
#     instance of the widget
#
# Results:
#
#     All state variables are set to their default values; all of 
#     the option database entries will exist.
#
# Returns:
# 
#     empty string

proc ::progressbar::Init {} {
  variable widgetOptions
  variable widgetCommands
  variable widgetGlobals
  variable widgetShapes

  if {$widgetGlobals(debug)} {
    puts stderr "pb_Init"
  }

  # here we match up command line options with option database names
  # and classes. As it turns out, this is a handy reference of all of the
  # available options. Note that if an item has a value with only one
  # item (like -bd, for example) it is a synonym and the value is the
  # actual item.

  array set widgetOptions [list \
    -background          {background          Background} \
    -bd                  -borderwidth \
    -bg                  -background \
    -borderwidth         {borderWidth         BorderWidth} \
    -color               {color               Color} \
    -pc                  -percent \
    -percent             {percent             Percent} \
    -relief              {relief              Relief} \
    -shape               {shape               Shape} \
    -variable            {variable            Variable} \
    -width               {width               Width} \
  ] 

  # this defines the valid widget commands. It's important to
  # list them here; we use this list to validate commands and
  # expand abbreviations.

  set widgetCommands [list \
    cget \
    configure \
  ]

  # this defines the valid shape options. It's important to
  # list them here; we use this list to validate options and
  # expand abbreviations.

  set widgetShapes [list \
    3D \
    3d \
    flat \
  ]

  set widgetGlobals(toDraw) {
    rect #bdbdbd es0 {[expr $mark +3] 2 [expr $width -2] 11} {-outline ""}
    line #525252 es1 {[expr $mark +1] 2 [expr $mark +1] 11} {}
    line #8c8c8c es2 {[expr $mark +2] 11 [expr $mark +2] 2 \
      [expr $width -4] 2} {}
    line #8c8c8c es3 {[expr $mark +3] 11 [expr $width -3] 11 \
      [expr $width -3] 3} {}
    line $rgb(0) pb0 {4 11 [expr $mark -1] 11 [expr $mark -1] 3} {}
    line $rgb(1) pb1 {3 11 3 10 [expr $mark -2] 10 [expr $mark -2] 2 \
      [expr $mark -1] 2 4 2} {}
    line $rgb(2) pb2 {3 2 2 2 2 11 2 10 3 10 3 9 [expr $mark -3] 9 \
      [expr $mark -3] 3 [expr $mark -2] 3 4 3} {}
    line $rgb(3) pb3 {3 3 3 9 3 8 [expr $mark -3] 8 [expr $mark -3] 4 4 4} {}
    line $rgb(4) pb4 {3 4 3 8 3 7 [expr $mark -3] 7 [expr $mark -3] 5 4 5} {}
    line $rgb(5) pb5 {3 5 3 7 3 6 [expr $mark -3] 6} {}
    line #000000 mrk {$mark 1 $mark 12} {}
    line #adadad fr0 {0 12 0 0 [expr $width -1] 0} {}
    line #ffffff fr1 {1 13 [expr $width -1] 13 [expr $width -1] 1} {}
    line #000000 fr2 {1 1 [expr $width -2] 1 [expr $width -2] 12 1 12 1 1} {}
  }

  set widgetGlobals(@blue0)\
    {#000052 #0031ce #3163ff #639cff #9cceff #efefef}
  set widgetGlobals(@blue1)\
    {#000021 #00639c #009cce #00ceff #63ffff #ceffff}
  set widgetGlobals(@blue2)\
    {#000052 #31319c #6363ce #9c9cff #ceceff #efefef}
  set widgetGlobals(@blue3)	\
    {#21214a #52527b #63639c #8484bd #b5b5ef #ceceff}
  set widgetGlobals(@blue4)\
    {#29396b #4a6b9c #6384b5 #739cd6 #94b5ef #adceff}
  set widgetGlobals(@green0)	\
    {#003131 #08736b #318c94 #5abdad #63dece #ceffef}
  set widgetGlobals(@green1)\
    {#001000 #003100 #316331 #639c63 #9cce9c #ceffce}
  set widgetGlobals(@green2)\
    {#002100 #006331 #319c63 #31ce63 #63ff9c #ceffce}
  set widgetGlobals(@green3)\
    {#003131 #316363 #427b7b #639c9c #9ccece #bdefef}
  set widgetGlobals(@yellow0)\
    {#101010 #636300 #9c9c00 #cece00 #ffff00 #ffff9c}
  set widgetGlobals(@yellow1)\
    {#8c7321 #cead39 #e7c642 #f7de63 #f7de63 #ffffe7}
  set widgetGlobals(@red0)\
    {#420000 #9c0000 #ce3131 #ff6363 #ff9c9c #ffcece}
  set widgetGlobals(@red1)\
    {#210000 #9c3100 #ce6331 #ff9c63 #ffce9c #ffffce}
  set widgetGlobals(@magenta0)\
    {#210000 #630063 #9c319c #ce63ce #ff9cff #ffceff}
  set widgetGlobals(@brown0)\
    {#210000 #633100 #9c6331 #ce9c63 #efb573 #ffdeb5}
  set widgetGlobals(@brown1)\
    {#310000 #7b4242 #9c6363 #ce9c9c #efcece #ffdede}
  set widgetGlobals(@gray0)\
    {#212121 #525252 #737373 #adadad #cecece #efefef}

  # this initializes the option database. Kinda gross, but it works
  # (I think).
  set tmpWidget ".__tmp__"

  # steal some options from frame widgets; we only want a subset
  # so we'll use a slightly different method. No harm in *not*
  # adding in the one or two that we don't use... :-)
  label $tmpWidget
  foreach option [list Background Relief] {
    set values [$tmpWidget configure -[string tolower $option]]
    option add *Progressbar.$option [lindex $values 3]
  }
  destroy $tmpWidget

  # these are unique to us...
  option add *Progressbar.borderWidth	5		widgetDefault
  option add *Progressbar.color		@blue0		widgetDefault
  option add *Progressbar.percent	0		widgetDefault
  option add *Progressbar.shape		3D		widgetDefault
  option add *Progressbar.variable	{}		widgetDefault
  option add *Progressbar.width		180		widgetDefault

  # define the class bindings
  # this allows us to clean up some things when we go away
  bind Progressbar <Destroy> [list ::progressbar::DestroyHandler %W]
}


# ::progressbar::progressbar --
#
#     This is the command that gets exported. It creates a new
#     progressbar widget.
#
# Arguments:
#
#     w        path of new widget to create
#     args     additional option/value pairs (eg: -background white, etc.)
#
# Results:
#
#     It creates the widget and sets up all of the default bindings
#
# Returns:
#
#     The name of the newly create widget

proc ::progressbar::progressbar {args} {
  variable widgetOptions
  variable widgetGlobals

  if {$widgetGlobals(debug)} {
    puts stderr "pb_progressbar '$args'"
  }

  # perform a one time initialization
  if {![info exists widgetOptions]} {
    Init
  }

  # make sure we at least have a widget name
  if {[llength $args] == 0} {
    error "wrong # args: should be \"progressbar pathName ?options?\""
  }

  # ... and make sure a widget doesn't already exist by that name
  if {[winfo exists [lindex $args 0]]} {
    error "window name \"[lindex $args 0]\" already exists"
  }

  # and check that all of the args are valid
  foreach {name value} [lrange $args 1 end] {
    Canonize [lindex $args 0] option $name
  }

  # build it...
  set w [eval Build $args]

  # and we are done!
  return $w
}


# ::progressbar::Build --
#
#    This does all of the work necessary to create the basic
#    progressbar. 
#
# Arguments:
#
#    w        widget name
#    args     additional option/value pairs
#
# Results:
#
#    Creates a new widget with the given name. Also creates a new
#    namespace patterened after the widget name, as a child namespace
#    to ::progressbar
#
# Returns:
#
#    the name of the widget

proc ::progressbar::Build {w args} {
  variable widgetOptions
  variable widgetGlobals

  if {$widgetGlobals(debug)} {
    puts stderr "pb_Build '$w' '$args'"
  }

  # create the namespace for this instance, and define a few
  # variables
  namespace eval ::progressbar::$w {
    variable options
    variable widgets
  }

  # this gives us access to the namespace variables within
  # this proc
  upvar ::progressbar::${w}::widgets widgets
  upvar ::progressbar::${w}::options options

  # this is our widget -- a frame of class Progressbar. Naturally,
  # it will contain other widgets. We create it here because
  # we need it to be able to set our default options.
  set widgets(this) [frame $w -class Progressbar]

  # this defines all of the default options. We get the
  # values from the option database. Note that if an array
  # value is a list of length one it is an alias to another
  # option, so we just ignore it
  foreach name [array names widgetOptions] {
    if {[llength $widgetOptions($name)] == 1} continue
    set optName  [lindex $widgetOptions($name) 0]
    set optClass [lindex $widgetOptions($name) 1]
    set options($name) [option get $w $optName $optClass]
    if {$widgetGlobals(debug) > 1} {
      puts stderr "pb_Build:Opt '$w' '$optName' '$optClass' '$options($name)'"
    }
  }

  # now apply any of the options supplied on the command
  # line. This may overwrite our defaults, which is OK
  if {[llength $args] > 0} {
    array set options $args
  }
  
  # this will only set the name of canvas's widget, we will
  # later create the canvas in our drawing procedure.
  set widgets(canvas) $w.pb

  # we will later rename the frame's widget proc to be our
  # own custom widget proc. We need to keep track of this
  # new name, so we'll define and store it here...
  set widgets(frame) ::progressbar::${w}::$w

  # this moves the original frame widget proc into our
  # namespace and gives it a handy name
  rename ::$w $widgets(frame)

  # now, create our widget proc. Obviously (?) it goes in
  # the global namespace. All progressbar widgets will actually
  # share the same widget proc to cut down on the amount of
  # bloat. 
  proc ::$w {command args} \
    "eval ::progressbar::WidgetProc {$w} \$command \$args"

  # ok, the thing exists... let's do a bit more configuration. 
  if {[catch "Configure $widgets(this) [array get options]" error]} {
    catch {destroy $w}
  }

  return $w
}


# ::progressbar::WidgetProc --
#
#    This gets uses as the widgetproc for an progressbar widget. 
#    Notice where the widget is created and you'll see that the
#    actual widget proc merely evals this proc with all of the
#    arguments intact.
#
#    Note that some widget commands are defined "inline" (ie:
#    within this proc), and some do most of their work in 
#    separate procs. This is merely because sometimes it was
#    easier to do it one way or the other.
#
#    w         widget pathname
#    command   widget subcommand
#    args      additional arguments; varies with the subcommand
#
# Results:
#
#    Performs the requested widget command

proc ::progressbar::WidgetProc {w command args} {
  variable widgetOptions
  variable widgetGlobals

  if {$widgetGlobals(debug)} {
    puts stderr "pb_WidgetProc '$w' '$command' '$args'"
  }

  upvar ::progressbar::${w}::widgets   widgets
  upvar ::progressbar::${w}::options   options

  set command [Canonize $w command $command]

  set result ""

  # here we go. Error checking be damned!
  switch $command {
    cget {
      if {[llength $args] != 1} {
	error "wrong # args: should be $w cget option"
      }
      set opt [Canonize $w option [lindex $args 0]]
      set result $options($opt)
    }

    configure {
      set result [eval Configure {$w} $args]
    }
  }
  return $result
}


# ::progressbar::HumanizeList --
#
#    Returns a human-readable form of a list by separating items
#    by columns, but separating the last two elements with "or"
#    (eg: foo, bar or baz)
#
# Arguments:
#
#    list    a valid tcl list
#
# Results:
#
#    A string which as all of the elements joined with ", " or 
#    the word " or "

proc ::progressbar::HumanizeList {list} {
  variable widgetGlobals

  if {$widgetGlobals(debug)} {
    puts stderr "pb_HumanizeList $list"
  }

  if {[llength $list] == 1} {
    return [lindex $list 0]
  } else {
    set list [lsort $list]
    set secondToLast [expr {[llength $list] -2}]
    set most [lrange $list 0 $secondToLast]
    set last [lindex $list end]

    return "[join $most {, }] or $last"
  }
}


# ::progressbar::Canonize --
#
#    takes a (possibly abbreviated) option or command name and either 
#    returns the canonical name or an error
#
# Arguments:
#
#    w        widget pathname
#    object   type of object to canonize; must be one of "command",
#             "option", "column" or "column option".
#    opt      the option (or command) to be canonized
#
# Returns:
#
#    Returns either the canonical form of an option or command,
#    or raises an error if the option or command is unknown or
#    ambiguous.

proc ::progressbar::Canonize {w object opt} {
  variable widgetOptions
  variable widgetCommands
  variable widgetGlobals
  variable widgetShapes

  if {$widgetGlobals(debug)} {
    puts stderr "pb_Canonize '$w' '$object' '$opt'"
  }

  switch $object {
    command {
      if {[lsearch -exact $widgetCommands $opt] >= 0} {
	return $opt
      }

      # command names aren't stored in an array, and there
      # isn't a way to get all the matches in a list, so
      # we'll stuff the columns in a temporary array so
      # we can use [array names]
      set list $widgetCommands
      foreach element $list {
	set tmp($element) ""
      }
      set matches [array names tmp ${opt}*]
    }

    option {
      if {[info exists widgetOptions($opt)] \
	  && [llength $widgetOptions($opt)] == 3} {
	return $opt
      }
      set list [array names widgetOptions]
      set matches [array names widgetOptions ${opt}*]
    }

    shape {
      if {[lsearch -exact $widgetShapes $opt] >= 0} {
	return $opt
      }

      # same procedure as command
      set list $widgetShapes
      foreach element $list {
	set tmp($element) ""
      }
      set matches [array names tmp ${opt}*]
    }
  }
  if {[llength $matches] == 0} {
    set choices [HumanizeList $list]
    error "unknown $object \"$opt\"; must be one of $choices"
  } elseif {[llength $matches] == 1} {
    # deal with option aliases
    set opt [lindex $matches 0]
    switch $object {
      option {
	if {[llength $widgetOptions($opt)] == 1} {
	  set opt $widgetOptions($opt)
	}
      }
    }
    return $opt
  } else {
      set choices [HumanizeList $list]
      error "ambiguous $object \"$opt\"; must be one of $choices"
  }
}


# ::progressbar::RGBs --
#
#    Calculates RGB colors
#
# Arguments:
#
#    color  basic color as rgb or name
#
# Returns:
#    
#    A list of 6 calculated RGBs values and the original value.
 
proc ::progressbar::RGBs {color} {
  variable widgetGlobals

  if {$widgetGlobals(debug)} {
    puts stderr "pb_RGB '$color'"
  }

  # get rgb values of given color
  set color [winfo rgb . $color]

  set R [expr int([lindex $color 0] / 256)]
  set G [expr int([lindex $color 1] / 256)]
  set B [expr int([lindex $color 2] / 256)]

  set rgb {}
  foreach factor {0.13 0.32 0.45 0.68 0.8 0.93} {
    set r [expr int($R * $factor)]
    set g [expr int($G * $factor)]
    set b [expr int($B * $factor)]
    lappend rgb [format "#%02x%02x%02x" $r $g $b]
  }
  lappend rgb [format "#%02x%02x%02x" $R $G $B]

  return $rgb
}


# ::progressbar::Configure --
#
#    Implements the "configure" widget subcommand
#
# Arguments:
#
#    w      widget pathname
#    args   zero or more option/value pairs (or a single option)
#
# Results:
#    
#    Performs typcial "configure" type requests on the widget
 
proc ::progressbar::Configure {w args} {
  variable widgetOptions
  variable widgetGlobals

  if {$widgetGlobals(debug)} {
    puts stderr "pb_Configure '$w' '$args'"
  }

  upvar ${w}::widgets widgets
  upvar ${w}::options options
  
  if {[llength $args] == 0} {
    # hmmm. User must be wanting all configuration information
    # note that if the value of an array element is of length
    # one it is an alias, which needs to be handled slightly
    # differently
    set results {}
    foreach opt [lsort [array names widgetOptions]] {
      if {[llength $widgetOptions($opt)] == 1} {
	set alias $widgetOptions($opt)
	set optName $widgetOptions($alias)
	lappend results [list $opt $optName]
      } else {
	set optName  [lindex $widgetOptions($opt) 0]
	set optClass [lindex $widgetOptions($opt) 1]
	set default [option get $w $optName $optClass]
	lappend results [list $opt $optName $optClass $default $options($opt)]
      }
    }
    return $results
  }
  
  # one argument means we are looking for configuration
  # information on a single option
  if {[llength $args] == 1} {
    set opt [Canonize $w option [lindex $args 0]]
    set optName  [lindex $widgetOptions($opt) 0]
    set optClass [lindex $widgetOptions($opt) 1]
    set default [option get $w $optName $optClass]
    set results [list $opt $optName $optClass $default $options($opt)]
    return $results
  }

  # if we have an odd number of values, bail. 
  if {[expr {[llength $args]%2}] == 1} {
    # hmmm. An odd number of elements in args
    error "value for \"[lindex $args end]\" missing"
  }
  
  # Great. An even number of options. Let's make sure they 
  # are all valid before we do anything. Note that Canonize
  # will generate an error if it finds a bogus option; otherwise
  # it returns the canonical option name
  foreach {name value} $args {
    set name [Canonize $w option $name]
    set opts($name) $value
  }

  # process all of the configuration options
  foreach option [array names opts] {
    set newValue $opts($option)
    if {[info exists options($option)]} {
      set oldValue $options($option)
    }

    if {$widgetGlobals(debug) > 2} {
      puts stderr "pb_Configure:Opt '$option' n='$newValue' o='$oldValue'"
    }
    switch -- $option {
      -background  -
      -borderwidth -
      -relief      {
	if {[winfo exists $widgets(this)]} {
	  $widgets(frame) configure $option $newValue
	  set options($option) [$widgets(frame) cget $option]
	}
      }
      -color {
        switch -- $newValue {
	  @blue0    -
	  @blue1    -
	  @blue2    -
	  @blue3    -
	  @blue4    -
	  @green0   -
	  @green1   -
	  @green2   -
	  @green3   -
	  @yellow0  -
	  @yellow1  -
	  @red0     -
	  @red1     -
	  @magenta0 -
	  @brown0   -
	  @brown1   -
	  @gray0    {
	    set options(rgb) $widgetGlobals($newValue)
	  }
	  @* {
	    set options(rgb) $widgetGlobals(@saphir)
	  }
	  default {
	    set options(rgb) [RGBs $newValue]
	  }
	}
	set options(rgbHasChanged) 1
      }
      -percent {
	set options($option) $newValue
      }
      -shape {
	set options($option) [Canonize $w shape $newValue]
	set options(rgbHasChanged) 1
      }
      -variable {
	# hmmm .. are there any traces left? Yes! Destroy!
	if {[info procs Trace($w)] != ""} {
	  uplevel 3 trace vdelete $oldValue wu ::progressbar::Trace($w)
	  unset widgetGlobals($w)
	  rename Trace($w) {}
	}
	if {$newValue != ""} {
	  # there is a new variable to trace. build a new proc to trace it.
	  proc ::progressbar::Trace($w) {name1 name2 op} "
	    variable widgetGlobals

	    if {\$widgetGlobals(debug)} {
	      puts stderr \"pb_Trace($w) '\$name1' '\$name2' '\$op'\"
	    }
	    switch -- \$op {
	      w {
		if {\$name2 != \"\"} {
		  upvar 1 \${name1}(\$name2) var
		  catch {$w configure -percent \$var}
		} else {
		  upvar 1 \$name1 var
		  catch {$w configure -percent \$var}
		}
	      }
	      u {
		if {\[info procs Trace($w)\] != \"\"} { \
		  unset widgetGlobals($w); \
		  rename Trace($w) {}; \
		}
	      }
	    }
	  "
	  # install trace proc for variable
	  uplevel 3 trace variable $newValue wu ::progressbar::Trace($w)
	}
	set options($option) $newValue
	set widgetGlobals($w) $newValue
      }
      -width {
	if {$newValue < 20} {
	  error "a -width of less than 20 is not supported."
	}
	if {[winfo exists $widgets(canvas)]} {
	  $widgets(canvas) configure $option $newValue
	  set options($option) [$widgets(canvas) cget $option]
	} else {
          set options($option) $newValue
	}
      }
    }
  }

  Draw $w
}


# ::progressbar::DestroyHandler {w} --
# 
#    Cleans up after a progressbar widget is destroyed
#
# Arguments:
#
#    w    widget pathname
#
# Results:
#
#    The namespace that was created for the widget is deleted,
#    the widget proc and variable traces are removed.

proc ::progressbar::DestroyHandler {w} {
  variable widgetGlobals

  if {$widgetGlobals(debug)} {
    puts stderr "pb_DestroyHandler '$w'"
  }

  # hmmm .. are there any traces left? Yes! Destroy!
  if {[info procs Trace($w)] != ""} {
    uplevel 1 trace vdelete $widgetGlobals($w) wu ::progressbar::Trace($w)
    unset widgetGlobals($w)
    rename Trace($w) {}
  }

  # if the widget actually being destroyed is of class Progressbar,
  # crush the namespace and kill the proc. Get it? Crush. Kill. 
  # Destroy. Heh. Danger Will Robinson! Oh, man! I'm so funny it
  # brings tears to my eyes.
  if {[string compare [winfo class $w] "Progressbar"] == 0} {
    namespace delete ::progressbar::$w
    rename $w {}
  }
}


# ::progressbar::Draw --
#
#    Implements the draw subroutine
#
# Arguments:
#
#    w      widget pathname
#
# Results:
#    
#    Performs the drawing of progressbar

proc ::progressbar::Draw {w} {
  variable widgetGlobals

  if {$widgetGlobals(debug) > 2} {
    puts stderr "pb_Draw '$w'"
  }

  upvar ${w}::widgets widgets
  upvar ${w}::options options

  set width   $options(-width)
  set percent $options(-percent)

  if {$options(-shape) == "flat"} {
    set minDisplay 0
    if {[llength $options(rgb)] == 7} {
      set rgb(0) [lindex $options(rgb) 6]
    } else {
      set rgb(0) [lindex $options(rgb) 2]
    }
    set rgb(1) $rgb(0)
    set rgb(2) $rgb(0)
    set rgb(3) $rgb(0)
    set rgb(4) $rgb(0)
    set rgb(5) $rgb(0)
  } else {
    set minDisplay 7
    set rgb(0) [lindex $options(rgb) 0]
    set rgb(1) [lindex $options(rgb) 1]
    set rgb(2) [lindex $options(rgb) 2]
    set rgb(3) [lindex $options(rgb) 3]
    set rgb(4) [lindex $options(rgb) 4]
    set rgb(5) [lindex $options(rgb) 5]
  }

  if {$percent < 0} {
    set percent 0
  } elseif {$percent > 100} {
    set percent 100
  }
  if {$percent == 0} {
    set mark $minDisplay
  } else {
    set mark [expr (($width - $minDisplay) / 100.0 * $percent) + $minDisplay]
  }

  if {![winfo exists $widgets(canvas)]} {
    canvas $widgets(canvas) -width $width -height 14 -bd 0 -highlightthickness 0
    pack $widgets(canvas) -side left -anchor nw -fill both

    foreach {type color tag coords opts} $widgetGlobals(toDraw) {
      eval $widgets(canvas) create $type $coords -fill $color -tag t$tag $opts
    }

    set options(rgbHasChanged) 0
    # nothing more to do
    return
  }

  foreach {type color tag coords opts} $widgetGlobals(toDraw) {
    eval $widgets(canvas) coords t$tag $coords
    if {$options(rgbHasChanged)} {
      eval $widgets(canvas) itemconfigure t$tag -fill $color
    }
  }
  set options(rgbHasChanged) 0
}

