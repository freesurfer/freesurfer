##-----------------------------------------------------------------------
##                   
##        This file is part of the EMU Speech Database System
##           (c) Copyright 1996, 1997
##       Speech Hearing and Language Research Centre,
##         Macquarie University, Sydney, Australia.
##
##  Permission to use, copy, modify, distribute this software and its    
##  documentation for research, educational and individual use only, is  
##  hereby granted without fee, subject to the following conditions:     
##   1. The code must retain the above copyright notice, this list of    
##      conditions and the following disclaimer.                         
##   2. Any modifications must be clearly marked as such.                
##   3. Original authors' names are not deleted.                         
##  This software may not be used for commercial purposes without        
##  specific prior written permission from the authors.                  
##                                  
##  MACQUARIE UNIVERSITY AND THE CONTRIBUTORS TO THIS WORK        
##  DISCLAIM ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING      
##  ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT   
##  SHALL MACQUARIE UNIVERSITY NOR THE CONTRIBUTORS BE LIABLE     
##  FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES    
##  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN   
##  AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION,          
##  ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF       
##  THIS SOFTWARE. THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, AND THE
##  AUTHORS AND DISTRIBUTORS HAVE NO OBLIGATION TO PROVIDE MAINTENANCE,
##  SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
##
##-----------------------------------------------------------------------
## $Id: tkm_graph.tcl,v 1.1 2000/03/23 17:21:30 inverse Exp $ 


#------------------------------------------------------------------
# 
# emu_graph widget 
# 
#-------------------------------------------------------------------


## legal options 
set emu_graph(options) { 
    width height xref yref ticklen axistextoffset
    nticks_x nticks_y font autorange canvas xmin xmax ymin ymax
}

## default values for some options
set emu_graph(default) { \
    -width           300\
    -height          200\
    -xref            50\
    -yref            30\
    -ticklen         5\
    -axistextoffset  5\
    -nticks_x         5\
    -nticks_y         5\
    -font            fixed \
    -autorange       1\
}

set emu_graph(dataoptions) {
    points lines colour coords mask maskthresh trackdata trackcolumn redraw
}

set emu_graph(datadefault) {
    -points          0 \
    -lines           1\
    -colour          red\
    -trackcolumn     0 \
    -redraw          1
}

# here we will record all the graph names as they are made
set emu_graph(graphs) {}

proc emu_graph args {
    global emu_graph

    set graph [lindex $args 0]
    lappend emu_graph(graphs) $graph
    
    ## remove any existing config info under this name
    foreach key [array names emu_graph] {
        if [string match "$graph,*" $key] {
            unset emu_graph($key)
        }
    }

    ## prepend the default options to args, they can then be overridden if they
    ## appear later in the args
    set args [concat $graph 0 $emu_graph(default) [lrange $args 1 end]] 

    ## now parse the options
    set restargs [eval "emu_graph:internal_configure $args"]
    
    # shouldn't be any more args
    if { $restargs != {} } {
  error "Usage: emu_graph graph \[options\] ($restargs)"
    }
    set emu_graph($graph,datasets) {}
    
    # define the widget command
    proc $graph args "\
  global emu_graph; \
  set method emu_graph..\[lindex \$args 0\];\
  if {\[info procs \$method\] != {}} {\
      eval \[concat \$method $graph \[lrange \$args 1 end\]\];\
  } else { \
      error \"no method \[lindex \$args 0\] for emu_graph,\n options are \[emu_graph:methods\]\"\
        }"
}

## find the names of all methods emu_graph..*, just giving the * bit
proc emu_graph:methods {} {
    set mm [info procs "emu_graph..*"]
    set result {}
    foreach meth $mm {
        regexp "emu_graph..(.*)" $meth ignore realmeth
        lappend result $realmeth
    }
    return $result
}

proc emu_graph..data args {

    global emu_graph

    set graph [lindex $args 0]
    set tag [lindex $args 1]

    if {[llength $tag]>1 || [string match "-*" $tag]} {
  error "Usage: graph data tag \[options\] x y"
    }

    set args [concat $graph $tag $emu_graph(datadefault) [lrange $args 2 end]] 

    ## now parse the options
    set restargs [eval "emu_graph:internal_configure $args"]
    
    if { [llength $restargs] != 0 } {
  error "Usage: graph data tag \[options\]"
    }

    ## append tag only if not already exists, Mark Koennecke
    set mark_list $emu_graph($graph,datasets)
    if { [lsearch -exact $mark_list $tag] < 0 } {    
       set emu_graph($graph,datasets) [lappend emu_graph($graph,datasets) $tag]
    }
    set datalength 0 
    ## if we have data as coords then verify that each element is a pair 
    ## and remember the length for later
    if { [info exists emu_graph($graph,$tag,coords)] } {
  set ncoords [llength $emu_graph($graph,$tag,coords)]
  if { int($ncoords/2)*2 != $ncoords } {
      set emu_graph($graph,$tag,coords) {}
      error "bad data format in emu_graph $graph, data $tag\n -- length of coords must be even, it was $ncoords"
  }
  set datalength [expr $ncoords/2]
    }
    ## if we have data as trackdata, remember it's length
    if { [info exists emu_graph($graph,$tag,trackdata)] } {
  set datalength [$emu_graph($graph,$tag,trackdata) length]
    }

    # if there's a mask, chech that there's also a maskthresh and 
    # that the length of the mask is the same as the data
    if { $datalength != 0 && [info exists emu_graph($graph,$tag,mask)] } {
  if { ![info exists emu_graph($graph,$tag,maskthresh)] } {
      error "No threshold supplied with masking vector in emu_graph, use -maskthresh N"
  }
  if { [llength $emu_graph($graph,$tag,mask)] != $datalength } {
      error "Mask vector and coords have different lengths ([llength $emu_graph($graph,$tag,$mask)] and  $datalength)"
  }
    }
    if {$datalength != 0 && $emu_graph($graph,$tag,redraw)} {
  emu_graph..redraw $graph
    }
}

## make an image the backdrop of the graph, fit the graph axes around the
## image -- used for displaying a spectrogram image under formant plots
proc emu_graph..image {graph image xmin xmax ymin ymax} {
    global emu_graph
    ## if we're doing this then the image dictates the axis ranges
    set emu_graph($graph,autorange) 0
    set emu_graph($graph,xmin) $xmin
    set emu_graph($graph,xmax) $xmax
    set emu_graph($graph,ymin) $ymin
    set emu_graph($graph,ymax) $ymax

    set emu_graph($graph,width) [image width $image]
    set emu_graph($graph,height) [image height $image]

    set emu_graph($graph,image) $image
    
    emu_graph..redraw $graph
}



proc emu_graph..configure args {
    set newargs [concat [lindex $args 0] 0 [lrange $args 1 end]]
    eval "emu_graph:internal_configure $newargs"
}

proc emu_graph:internal_configure args {
    ## rest of args is a list of option value pairs, set emu_graph($canvas,option)
    ## to value for each option, if args remain after last option (ie
    ## something not beginning with a - or after a --, they are returned

    global emu_graph
    
    set graph [lindex $args 0]
    set datatag [lindex $args 1]
    set args [lrange $args 2 end]
    
    if {![emu_graph..is_graph $graph]} {
        error "$graph is not an emu_graph"
    }

    # if we're setting options for a data set we modify $graph
    # to include the tag to make the array entry 
    # emu_graph($graph,$tag,$option)
    if { $datatag != 0 } {
  set graph "$graph,$datatag"    
  set validoptions $emu_graph(dataoptions)
    } else {
  set validoptions $emu_graph(options)
    }  

    
    set curropt ""
    for {set i 0} { $i<[llength $args] } { incr i 2 } {
        if { [lindex $args $i] == "--" } {
            # terminating switch, return rest of args
            return [lrange $args [expr $i+1] end]
        } elseif { [regexp -- "-(.*)" [lindex $args $i] ignore curropt] } {
            # have a switch, get value and set option
            # but check first that it's kosher
            if { [lsearch $validoptions $curropt] >= 0 } {
                if { $i+1<[llength $args] } {
                    set emu_graph($graph,$curropt) [lindex $args [expr $i+1]]
                }
            } else {
                error "Bad option -$curropt to emu_graph\n\tuse one of $validoptions"
            }
        } else {
            ## options have run out, return rest of args
            return [lrange $args $i end]
        }
    }
}

proc emu_graph..redraw {graph} {
    global emu_graph
    if {![emu_graph..is_graph $graph]} {
        error "$graph is not an emu_graph"
    }
    # draw it if we have a canvas
    if {[info exists emu_graph($graph,canvas)]} {
  $emu_graph($graph,canvas) delete withtag graph$graph
        emu_graph:axes $graph
        emu_graph:plot_data $graph
    } else {
        error "$graph isn't associated with a canvas, use the -canvas option"
    }
}


proc emu_graph..is_graph {graph} {
    global emu_graph
    return [expr [lsearch $emu_graph(graphs) $graph] >= 0]
}

proc emu_graph:auto_range {graph} {

    global emu_graph

    if {![emu_graph..is_graph $graph]} {
        error "$graph is not an emu_graph"
    }

    ## we only autorange if the option is on or if there is no range set
    if { $emu_graph($graph,autorange) ||
         !([info exists emu_graph($graph,xmin)] &&
           [info exists emu_graph($graph,xmax)] &&
           [info exists emu_graph($graph,ymin)] &&
           [info exists emu_graph($graph,ymax)]) } {


        set xyrange {{1e19 -1e19} {1e19 -1e19}}
        ## look at each dataset, find max/min for all
        foreach tag $emu_graph($graph,datasets) {

      if { [info exists emu_graph($graph,$tag,mask)] } {
    set mask $emu_graph($graph,$tag,mask)
    set maskthresh $emu_graph($graph,$tag,maskthresh)
      } else {
    set mask 0
    set maskthresh 0
      }

      if { [info exists emu_graph($graph,$tag,trackdata)] } {
    ## get ranges from the data
    set yrange [$emu_graph($graph,$tag,trackdata) \
        range $emu_graph($graph,$tag,trackcolumn)]
    ## xrange is start and end times
    set xrange [list [$emu_graph($graph,$tag,trackdata) start]\
        [$emu_graph($graph,$tag,trackdata) end]]
    set xyrange [list $xrange $yrange]
      } elseif { [info exists emu_graph($graph,$tag,coords)] } {
    set xyrange [emu_graph:maxrange $xyrange \
         [emu_graph:range\
              $emu_graph($graph,$tag,coords)\
              $mask $maskthresh]]
      }
        }


  set xrange [lindex $xyrange 0]
  set yrange [lindex $xyrange 1]
  
        set xextra 0
        set yextra 0
  
        set emu_graph($graph,xmin) [expr [lindex $xrange 0]-$xextra]
        set emu_graph($graph,xmax) [expr [lindex $xrange 1]+$xextra]
        set emu_graph($graph,ymin) [expr [lindex $yrange 0]-$yextra]
        set emu_graph($graph,ymax) [expr [lindex $yrange 1]+$yextra]
    }
}

proc emu_graph:plot_data {graph} {

    global emu_graph

    if {![emu_graph..is_graph $graph]} { 
        error "$graph is not an emu_graph"
    }

    set canvas $emu_graph($graph,canvas)

    if {[info exists emu_graph($graph,image)]} {
  $canvas create image \
      [emu_graph..x2canvas $graph $emu_graph($graph,xmin)]  \
      [emu_graph..y2canvas $graph $emu_graph($graph,ymax)] \
      -anchor nw -image $emu_graph($graph,image) \
      -tags [list graph$graph image$graph]
    }

    foreach tag $emu_graph($graph,datasets) {
        # plot the points, first delete any old ones
        $canvas delete -withtag $tag 

        set tags [list graph$graph data$graph $tag]

  if { [info exists emu_graph($graph,$tag,trackdata)] } {
      set coords \
    [$emu_graph($graph,$tag,trackdata) coords\
         $emu_graph($graph,$tag,trackcolumn)\
         $emu_graph($graph,xmin) $emu_graph($graph,xmax)\
         $emu_graph($graph,xfactor) $emu_graph($graph,xref)\
         $emu_graph($graph,ymin) $emu_graph($graph,ymax)\
         $emu_graph($graph,yfactor) $emu_graph($graph,yref)]
  } elseif { [info exists emu_graph($graph,$tag,coords)] } {
      set coords \
    [emu_graph:scale_points $graph $emu_graph($graph,$tag,coords)]
  } else {
      set coords {}
  }
  
  # we may have a masking vector
  if { [info exists emu_graph($graph,$tag,mask)] } {
      set mask $emu_graph($graph,$tag,mask)
      set maskthresh $emu_graph($graph,$tag,maskthresh)
  } else {
      set mask 0
  }
  
  if { $emu_graph($graph,$tag,lines) } {
      ## create the line as a single line item
      eval "$canvas create line $coords -fill $emu_graph($graph,$tag,colour) -tag {$tags}"
  }

        for {set i 0} {$i < [llength $coords]-1} {incr i 2} {
      ## we'll draw the point if were either not masking or if
      ## the mask value is over the threshold
      if { $mask == 0 || \
         [lindex $mask [expr $i/2]] >= $maskthresh } {
    set point [lrange $coords $i [expr $i+1]]
    if { [emu_graph:point_in_bounds $graph $point] } {
        
        if { $emu_graph($graph,$tag,points) } {
      set ox [lindex $point 0]
      set oy [lindex $point 1]
      $canvas create oval \
          [expr $ox-2] [expr $oy-2]\
          [expr $ox+2] [expr $oy+2]\
          -fill $emu_graph($graph,$tag,colour) \
          -outline $emu_graph($graph,$tag,colour) -width 0 \
          -tag $tags
        }
    }
      }
  }
    }
}
                   
#
# check whether point is in bounds, where point is already scaled to canvas coords
#
proc emu_graph:point_in_bounds {graph point} {
    global emu_graph
    set x [lindex $point 0]
    set y [lindex $point 1]
 
    if { $x >= $emu_graph($graph,xref) && 
   $x <= $emu_graph($graph,xref)+$emu_graph($graph,width)  &&
   $y <= $emu_graph($graph,yref)+$emu_graph($graph,height) && 
   $y >= $emu_graph($graph,yref) } {
  return 1 
    } else {
  return 0
    }
}


proc emu_graph:scale_factor {graph} {

    global emu_graph

    if {![emu_graph..is_graph $graph]} {
        error "$graph is not an emu_graph"
    }

    ## calculate scale factors for x and y
    set width  $emu_graph($graph,width)
    set height $emu_graph($graph,height)
    set xdelta [expr double($emu_graph($graph,xmax) - $emu_graph($graph,xmin))]
    set ydelta [expr double($emu_graph($graph,ymax) - $emu_graph($graph,ymin))]
    if {$xdelta == 0} { set xdelta 0.001 }
    if {$ydelta == 0} { set ydelta 0.001 }

    set xfactor [expr double($width)/$xdelta]
    set yfactor [expr double($height)/$ydelta]

    set emu_graph($graph,xfactor) $xfactor
    set emu_graph($graph,yfactor) $yfactor

}

proc emu_graph:axes {graph} {
    # generate axes for a plot
    global emu_graph

    if {![emu_graph..is_graph $graph]} {
        error "$graph is not an emu_graph"
    }

    ## make sure we have the current scale factors etc
    emu_graph:auto_range $graph
    emu_graph:scale_factor $graph

    set x_min $emu_graph($graph,xmin)
    set x_max $emu_graph($graph,xmax)
    set y_min $emu_graph($graph,ymin)  
    set y_max $emu_graph($graph,ymax)

    set y_min_c [emu_graph..y2canvas $graph $y_min]
    set y_max_c [emu_graph..y2canvas $graph $y_max]
    set x_min_c [emu_graph..x2canvas $graph $x_min]
    set x_max_c [emu_graph..x2canvas $graph $x_max]

    # parameters affecting axis drawing
    set ticklen        $emu_graph($graph,ticklen)
    set axistextoffset $emu_graph($graph,axistextoffset)
    set nticks_x       $emu_graph($graph,nticks_x)
    set nticks_y       $emu_graph($graph,nticks_y)
    set graphfont      $emu_graph($graph,font)

    set canvas         $emu_graph($graph,canvas)

    # clear up any existing axes
    $canvas delete -withtag axis

    $canvas create rect $x_min_c $y_min_c $x_max_c $y_max_c\
        -outline black -tag [list graph$graph axis]
                                          
    # y-pos of tick end points and of axis tick labels
    set ticky [expr $y_min_c-$ticklen]
    set texty [expr $y_min_c+$axistextoffset]
    # put ticks and numbers on the axis 
    # starting at next nice number above x_min
    set nicex_min [emu_graph:nicenum $x_min 1]
    set delta_x [emu_graph:nicenum [expr double($x_max-$x_min)/$nticks_x] 1]

    for {set t $nicex_min} {$t <= $x_max} {set t [expr $t+$delta_x]} {
  ## it may be that $f is one tick below y_min, don't draw it if it is
  ## this is because of a problem with rounding down in nicenum
  if {$t >= $x_min} {
      set x [emu_graph..x2canvas $graph $t]
      $canvas create line $x $y_min_c $x $ticky \
    -tag [list graph$graph axis]
      $canvas create line $x $y_max_c $x [expr $y_max_c+$ticklen]\
    -tag [list graph$graph axis]
      # and add the label
      $canvas create text [emu_graph..x2canvas $graph $t] $texty \
    -text $t -font $graphfont -tag [list graph$graph axis]\
    -anchor n
  }
    }

    # now the y axis
    set tickx1   [expr [emu_graph..x2canvas $graph $x_min]+$ticklen]
    set tickx2   [expr [emu_graph..x2canvas $graph $x_max]-$ticklen]
    set textx    [expr [emu_graph..x2canvas $graph $x_min]-$axistextoffset]

    set nicey_min [emu_graph:nicenum $y_min 1]   
    set delta_y [emu_graph:nicenum [expr double($y_max-$y_min)/$nticks_y] 1]

    for { set f $nicey_min } { $f <= $y_max } { set f [expr $f + $delta_y] } {
  ## it may be that $f is one tick below y_min, don't draw it if it is
  ## this is because of a problem with rounding down in nicenum
  if {$f >= $y_min} {
      set y [emu_graph..y2canvas $graph $f]
      $canvas create line [emu_graph..x2canvas $graph $x_min]\
    $y $tickx1 $y -tag [list graph$graph axis]
      $canvas create line [emu_graph..x2canvas $graph $x_max]\
    $y $tickx2 $y -tag [list graph$graph axis]
      # and add the label
      $canvas create text $textx $y -text $f -anchor e \
    -tag [list graph$graph axis] -font $graphfont
        }
    }
}

# scale_points with inlined scaling, Mark Koennecke
proc emu_graph:scale_points {graph coords} {
    global emu_graph

    if {![emu_graph..is_graph $graph]} {
        error "$graph is not an emu_graph"
    }

    set result {}
    for {set i 0} {$i < [llength $coords]-1} {incr i 2} {
  lappend result [expr int(([lindex $coords $i] - $emu_graph($graph,xmin)) \
                * $emu_graph($graph,xfactor) + $emu_graph($graph,xref))]

  lappend result [expr int(($emu_graph($graph,ymax) - \
         [lindex $coords [expr $i+1]]) \
                * $emu_graph($graph,yfactor) + $emu_graph($graph,yref))]
    }
    return $result
}

proc emu_graph..bbox {graph} {
    global emu_graph
    return [$emu_graph($graph,canvas) bbox graph$graph]
}

proc emu_graph..cget {graph option} { 
    global emu_graph
    # remove leading - if present
    if { [regexp -- "-(.*)" $option ignore optname] } {
  set option optname
    }
    # now look for it in the options store
    if {[info exists emu_graph($graph,$option)] } {
  return $emu_graph($graph,$option)
    } else {
  error "emu_graph has no configuration option $option"
    }
}


proc emu_graph..y2canvas {graph y} {
    global emu_graph

    if {![emu_graph..is_graph $graph]} {
        error "$graph is not an emu_graph"
    }

    return [expr int(($emu_graph($graph,ymax) - $y) \
                * $emu_graph($graph,yfactor) + $emu_graph($graph,yref))]
}

proc emu_graph..canvas2y {graph cy} {
    global emu_graph

    if {![emu_graph..is_graph $graph]} {
        error "$graph is not an emu_graph"
    }

    return [expr $emu_graph($graph,ymax) - \
          ($cy - $emu_graph($graph,yref))/$emu_graph($graph,yfactor)]
}

proc emu_graph..canvas2x {graph cx} {
    global emu_graph

    if {![emu_graph..is_graph $graph]} {
        error "$graph is not an emu_graph"
    }

    return [expr $emu_graph($graph,xmin) + \
          double($cx - $emu_graph($graph,xref))/double($emu_graph($graph,xfactor))]
}

proc emu_graph..x2canvas {graph x} {
    global emu_graph

    if {![emu_graph..is_graph $graph]} {
        error "$graph is not an emu_graph"
    }
    return  [expr int(($x - $emu_graph($graph,xmin)) \
                * $emu_graph($graph,xfactor) + $emu_graph($graph,xref))]
}


## find the ranges of a list of coordinates {{x y} {x' y'} {x'' y''}...}
## returns two ranges {{xmin xmax} {ymin ymax}}
proc emu_graph:range {list {mask 0} {maskthresh 0}} {
    set xmax -10e64
    set xmin 10e64
    set ymax -10e64
    set ymin 10e64
    for {set i 0} {$i < [llength $list]-1} {incr i 2} {
  set x [lindex $list $i]
  set y [lindex $list [expr $i+1]]

  if { $mask == 0 || \
     [lindex $mask [expr $i/2]] >= $maskthresh } {
  
      if {$y > $ymax} {
    set ymax $y
      }

      if {$y < $ymin} {
    set ymin $y
      }
  }
  # don't worry about the mask for x -- we still want to line up with
  # other plots 
  if {$x>$xmax} {
      set xmax $x
  }
  
  if {$x < $xmin} {
      set xmin $x
  }
    }
    return [list [list $xmin $xmax] [list $ymin $ymax]] 
}


## find the maxima of the sets of ranges a and b which are both {{xmin xmax} {ymin ymax}}
proc emu_graph:maxrange {a b} {
    return [list [emu_graph:maxrange-aux [lindex $a 0] [lindex $b 0]]\
    [emu_graph:maxrange-aux [lindex $a 1] [lindex $b 1]]]
}


## find the maxima of the ranges a and b which are both {min max} pairs
proc emu_graph:maxrange-aux {a b} {
    # get the smallest minimum
    if {[lindex $a 0] < [lindex $b 0]} {
        set first [lindex $a 0]
    } else {
        set first [lindex $b 0]
    }
    # and the largest maximum
    if {[lindex $a 1] > [lindex $b 1]} {
        set second [lindex $a 1]
    } else {
        set second [lindex $b 1]
    }
    return [list $first $second]
}

             
## translated from C-code in Blt, who got it from:
##      Taken from Paul Heckbert's "Nice Numbers for Graph Labels" in
##      Graphics Gems (pp 61-63).  Finds a "nice" number approximately
##      equal to x.

proc emu_graph:nicenum {x floor} {

    if {$x == 0} { return 0 }
    
    set negative 0
    if {$x < 0} {
        set x [expr -$x]
        set negative 1
    }

    set exponX [expr floor(log10($x))]
    set fractX [expr $x/pow(10,$exponX)]; # between 1 and 10
    if {$floor} {
        if {$fractX < 1.5} {
            set nf 1.0
        } elseif {$fractX < 3.0} {
            set nf 2.0
        } elseif {$fractX < 7.0} {
            set nf 5.0
        } else {                
         set nf 10.0
        }
    } elseif {$fractX <= 1.0} {
        set nf 1.0
    } elseif {$fractX <= 2.0} {
        set nf 2.0
    } elseif {$fractX <= 5.0} {
        set nf 5.0
    } else {
        set nf 10.0
    }
    if { $negative } {
        return [expr -$nf * pow(10,$exponX)]
    } else {
        return [expr $nf * pow(10,$exponX)]
    }
}                

#
# put a vertical or horizontal mark on the graph 
#
proc emu_graph..vmark {graph x tag {color red}} {
    global emu_graph
    if { $x >= $emu_graph($graph,xmin) && $x <= $emu_graph($graph,xmax) } {
  set tags [list graph$graph vmark $tag]
  # if there's already an item with this tag then delete it
  $emu_graph($graph,canvas) delete $tag
  set cx [emu_graph..x2canvas $graph $x]
  $emu_graph($graph,canvas) create line \
      $cx [emu_graph..y2canvas $graph $emu_graph($graph,ymin)]\
      $cx [emu_graph..y2canvas $graph $emu_graph($graph,ymax)]\
      -fill $color -tags $tags
    }
}

proc emu_graph..hmark {graph y tag {color red}} {
    global emu_graph
    if { $y >= $emu_graph($graph,ymin) && $y <= $emu_graph($graph,ymax) } {
  # if there's already an item with this tag then delete it
  $emu_graph($graph,canvas) delete $tag
  set tags [list graph$graph vmark $tag]
  set cy [emu_graph..y2canvas $graph $y]
  $emu_graph($graph,canvas) create line \
      [emu_graph..x2canvas $graph $emu_graph($graph,xmin)] $cy\
      [emu_graph..x2canvas $graph $emu_graph($graph,xmax)] $cy\
      -fill $color -tags $tags
    }
}

proc emu_graph..clearmark {graph tag} {
    global emu_graph
    $emu_graph($graph,canvas) delete $tag
}


proc emu_graph..movevmark {graph tag newx} {
    global emu_graph
    set cx [emu_graph..x2canvas $graph $newx]
    $emu_graph($graph,canvas) coords $tag \
      $cx [emu_graph..y2canvas $graph $emu_graph($graph,ymin)]\
      $cx [emu_graph..y2canvas $graph $emu_graph($graph,ymax)]
}

proc emu_graph..movehmark {graph tag newy} {
    global emu_graph
    set cy [emu_graph..y2canvas $graph $newy]
    $emu_graph($graph,canvas) coords $tag \
      [emu_graph..x2canvas $graph $emu_graph($graph,xmin)] $cy\
      [emu_graph..x2canvas $graph $emu_graph($graph,xmax)] $cy\
}

## if we're running inside the tcl-plugin start up a demo graph
if [info exists embed_args] {

set testdata(1) {5.2427 1293.93 5.2477 1331.18 5.2527 1389.83 5.2577 1452.29 5.2627 1483.41 5.2677 547.065 5.2727 726.348 5.2777 1479.08 5.2827 1479.08 5.2877 1295.13 5.2927 1259.15 5.2977 349.753 5.3027 422.846 5.3077 434.979 5.3127 458.625 5.3177 386.405 5.3227 383.799 5.3277 406.285 5.3327 301.165 5.3377 313.201 5.3427 341.098 5.3477 283.14 5.3527 261.498 5.3577 193.357 5.3627 1345.98 5.3677 113.406 5.3727 288.299 5.3777 332.471 5.3827 652.12 5.3877 769.216 5.3927 863.813 5.3977 1034.27 5.4027 952.184 5.4077 980.757 5.4127 1050.41 5.4177 1075.02 5.4227 1162.34 5.4277 1361.41 5.4327 787.691 5.4377 787.691 5.4427 1166.31 5.4477 1342.4 5.4527 149.062 5.4577 200.061 5.4627 79.2487 5.4677 1251.47 5.4727 1254.47 5.4777 1264.44 5.4827 1402.37 5.4877 102.372 5.4927 125.9 5.4977 141.364 5.5027 169.658 5.5077 175.506 5.5127 192.736 5.5177 202.981 5.5227 205.917 5.5277 213.05 5.5327 231.236 5.5377 236.086 5.5427 242.616 5.5477 165.794 5.5527 194.984 5.5577 200.887 5.5627 193.183 5.5677 194.956 5.5727 173.402 5.5777 165.963 5.5827 160.621 5.5877 155.256 5.5927 153.418 5.5977 149.943 5.6027 147.725 5.6077 147.008 5.6127 141.461 5.6177 146.052 5.6227 140.273 5.6277 126.636 5.6327 65.9178 5.6377 95.3058 5.6427 115.762 5.6477 108.024 5.6527 119.894 5.6577 61.8934 5.6627 61.2774 5.6677 94.3914 5.6727 161.267 5.6777 167.213 5.6827 181.157 5.6877 191.344 5.6927 200.091 5.6977 209.809 5.7027 213.563 5.7077 222.195 5.7127 229.789 5.7177 234.246 5.7227 239.496 5.7277 246.445 5.7327 253.146 5.7377 259.44 5.7427 265.127 5.7477 270.522 5.7527 274.769 5.7577 277.723 5.7627 279.279 5.7677 276.289 5.7727 275.203 5.7777 275.299 5.7827 273.669 5.7877 273.254 5.7927 273.9 5.7977 277.771 5.8027 284.016 5.8077 286.347 5.8127 286.552 5.8177 283.044 5.8227 265.953 5.8277 256.408 5.8327 252.82 5.8377 261.414 5.8427 263.411 5.8477 270.692 5.8527 279.679 5.8577 284.105 5.8627 287.694 5.8677 285.845 5.8727 257.974 5.8777 251.776 5.8827 247.367 5.8877 252.394 5.8927 256.593 5.8977 266.3 5.9027 259.87 5.9077 266.837 5.9127 263.266 5.9177 275.026 5.9227 274.145 5.9277 278.991 5.9327 281.76 5.9377 287.249 5.9427 289.784 5.9477 283.262 5.9527 290.556 5.9577 282.54 5.9627 292.311 5.9677 286.609 5.9727 276.07 5.9777 271.818 5.9827 243.69 5.9877 241.797 5.9927 245.383 5.9977 187.553 6.0027 180.693 6.0077 177.431 6.0127 273.668 6.0177 386.014 6.0227 331.685 6.0277 266.005 6.0327 215.005 6.0377 136.092 6.0427 161.138 6.0477 169.616 6.0527 139.869 6.0577 126.257 6.0627 134.615 6.0677 97.8916 6.0727 57.0648 6.0777 285.589 6.0827 285.589 6.0877 285.589 6.0927 341.395 6.0977 306.123 6.1027 492.25 6.1077 550.71 6.1127 525.322 6.1177 290.034 6.1227 135.066 6.1277 135.252 6.1327 155.995 6.1377 163.653 6.1427 608.521 6.1477 702.777 6.1527 65.9631 6.1577 217.884 6.1627 236.254 6.1677 132.745 6.1727 204.82 6.1777 452.779 6.1827 401.706 6.1877 294.19 6.1927 298.304 6.1977 268.792 6.2027 160.209 6.2077 1497.77 6.2127 226.696 6.2177 261.783 6.2227 203.92 6.2277 235.108 6.2327 245.154 6.2377 1489.82 6.2427 1489.82 6.2477 229.219 6.2527 299.341 6.2577 328.306 6.2627 348.67 6.2677 343.632 6.2727 343.792 6.2777 302 6.2827 285.527 6.2877 252.514 6.2927 186.773 6.2977 196.208 6.3027 218.415 6.3077 193.701 6.3127 53.904 6.3177 177.946 6.3227 115.312 6.3277 1134.81 6.3327 1072.45 6.3377 1085.3 6.3427 1309.02 6.3477 1117.63 6.3527 163.301 6.3577 216.819 6.3627 250.833 6.3677 149.651 6.3727 1348.28 6.3777 1018.11 6.3827 982.854 6.3877 1495.12 6.3927 1495.12 6.3977 671.37 6.4027 663.733 6.4077 675.005 6.4127 70.7996 6.4177 70.7996 6.4227 364.594 6.4277 692.625 6.4327 883.546 6.4377 924.274 6.4427 1430.86 6.4477 1350.97 6.4527 1197.4 6.4577 130.285 6.4627 220.82 6.4677 155.536 6.4727 118.564 6.4777 76.9127 6.4827 1148.64 6.4877 145.15 6.4927 188.478 6.4977 162.321 6.5027 1073.48 6.5077 1349.41 6.5127 163.626 6.5177 173.534 6.5227 166.885 6.5277 177.888 6.5327 1489.04 6.5377 1489.04 6.5427 1324.77 6.5477 1199.66 6.5527 129.896 6.5577 156.591 6.5627 179.702 6.5677 190.091 6.5727 176.552 6.5777 893.741 6.5827 1086.26 6.5877 1208.79 6.5927 110.355 6.5977 165.744 6.6027 146.021 6.6077 1359.5 6.6127 1209.28 6.6177 1205.59 6.6227 50.0568 6.6277 92.5615 6.6327 98.198 6.6377 59.4589 6.6427 1312.96 6.6477 52.0352 6.6527 127.656 6.6577 1136.44 6.6627 1193.87 6.6677 1324.67 6.6727 124.253 6.6777 145.293 6.6827 166.718 }

set testdata(2) {5.2427 2155.64 5.2477 2241.15 5.2527 2209.37 5.2577 2152.84 5.2627 2341.81 5.2677 1531.98 5.2727 1635.13 5.2777 1614.01 5.2827 2163.22 5.2877 2069.39 5.2927 2036.02 5.2977 1341.43 5.3027 1485.32 5.3077 1519.54 5.3127 1466.85 5.3177 1420.68 5.3227 1446.84 5.3277 1404.27 5.3327 1407.09 5.3377 1508.31 5.3427 1551.89 5.3477 1580.92 5.3527 2051.87 5.3577 2009.86 5.3627 2004.56 5.3677 1513.93 5.3727 1574.87 5.3777 1558.08 5.3827 1585.99 5.3877 1645.55 5.3927 1716.24 5.3977 1850.07 5.4027 1882.44 5.4077 1970.22 5.4127 2040.86 5.4177 2028.98 5.4227 2055.94 5.4277 2356.96 5.4327 1504.82 5.4377 1900.94 5.4427 2888.75 5.4477 2941.43 5.4527 1337.75 5.4577 1352.66 5.4627 1309.07 5.4677 2083.66 5.4727 2120.8 5.4777 2126.94 5.4827 2112.8 5.4877 2074.21 5.4927 2040.18 5.4977 1837.34 5.5027 2360.21 5.5077 2253.5 5.5127 2370.49 5.5177 2414.08 5.5227 2449.26 5.5277 2476.48 5.5327 2496.99 5.5377 2456.35 5.5427 2484.49 5.5477 2530.51 5.5527 2517.11 5.5577 2461.54 5.5627 2467.12 5.5677 2536.55 5.5727 2569.76 5.5777 2566.48 5.5827 2584.98 5.5877 2615.32 5.5927 2633.47 5.5977 2628.86 5.6027 2647.55 5.6077 2653.26 5.6127 2527.41 5.6177 2452.12 5.6227 2427.92 5.6277 2405.02 5.6327 2412.69 5.6377 2429.29 5.6427 2441.29 5.6477 2461.38 5.6527 2464.8 5.6577 2412.82 5.6627 2412.42 5.6677 2428.7 5.6727 2431.71 5.6777 2383.84 5.6827 2348.3 5.6877 2319.1 5.6927 2302.26 5.6977 2278.65 5.7027 1902.82 5.7077 1777.49 5.7127 1754.36 5.7177 1766.45 5.7227 1748.45 5.7277 1756.38 5.7327 1757.26 5.7377 1746.65 5.7427 1738.14 5.7477 1730 5.7527 1715.38 5.7577 1694.73 5.7627 1686.83 5.7677 1677.98 5.7727 1674.04 5.7777 1669.32 5.7827 1664.19 5.7877 1656.83 5.7927 1647.58 5.7977 1650.31 5.8027 1645.9 5.8077 1634.69 5.8127 1616.63 5.8177 1601.82 5.8227 1599.21 5.8277 1591.5 5.8327 1586.6 5.8377 1582.47 5.8427 1568.21 5.8477 1560.35 5.8527 1569.81 5.8577 1577.35 5.8627 1579.55 5.8677 1574.93 5.8727 1559.16 5.8777 1554.09 5.8827 1557.95 5.8877 1555.56 5.8927 1543.04 5.8977 1537.07 5.9027 1548.85 5.9077 1554.14 5.9127 1560.53 5.9177 1567.53 5.9227 1581.01 5.9277 1587.47 5.9327 1591.32 5.9377 1596.55 5.9427 1597.62 5.9477 1605.31 5.9527 1614.26 5.9577 1621.91 5.9627 1614.74 5.9677 1619.75 5.9727 1616.01 5.9777 1607.64 5.9827 1585.59 5.9877 1575.05 5.9927 1591.74 5.9977 1617.42 6.0027 1674.97 6.0077 1407.35 6.0127 1405.27 6.0177 1456.01 6.0227 1506.22 6.0277 1578.06 6.0327 1616.38 6.0377 1617.34 6.0427 1602.66 6.0477 1611.36 6.0527 1613.02 6.0577 1612.22 6.0627 1613.75 6.0677 1593.92 6.0727 1579.04 6.0777 1554.53 6.0827 1570.99 6.0877 1598.9 6.0927 1600.57 6.0977 1591.8 6.1027 1584.85 6.1077 1559.88 6.1127 1557.95 6.1177 1494.91 6.1227 1302.22 6.1277 1295.18 6.1327 1354.68 6.1377 1493.09 6.1427 1556.57 6.1477 1547.87 6.1527 1507.87 6.1577 1460.95 6.1627 1483.76 6.1677 1510.53 6.1727 1534.02 6.1777 1546.94 6.1827 1578.17 6.1877 1542.53 6.1927 1502.37 6.1977 1489.4 6.2027 1488.69 6.2077 2311.07 6.2127 1473.84 6.2177 1457.43 6.2227 1446.58 6.2277 1438.76 6.2327 1476.91 6.2377 1512.78 6.2427 2480.61 6.2477 1497.51 6.2527 1506.22 6.2577 1462.15 6.2627 1482.93 6.2677 1396.97 6.2727 1395.89 6.2777 1417.18 6.2827 2157.26 6.2877 1289.23 6.2927 1392.59 6.2977 1452.16 6.3027 1453.7 6.3077 1412.13 6.3127 1406.89 6.3177 1401.58 6.3227 1133.01 6.3277 2107.78 6.3327 1696.01 6.3377 2544.26 6.3427 2622.38 6.3477 2791.33 6.3527 2249.94 6.3577 1377.65 6.3627 1380.64 6.3677 1554.22 6.3727 2242.07 6.3777 2090.77 6.3827 1853.09 6.3877 1505.49 6.3927 2455.56 6.3977 2050.95 6.4027 1908.49 6.4077 1834.92 6.4127 1560.19 6.4177 1481.89 6.4227 1521.75 6.4277 1547.3 6.4327 1499.24 6.4377 1443.66 6.4427 2775.07 6.4477 2765.38 6.4527 2698.63 6.4577 2077.95 6.4627 1316.87 6.4677 1463.7 6.4727 1391.62 6.4777 1323.92 6.4827 2118.03 6.4877 1256.1 6.4927 1252.03 6.4977 1189.28 6.5027 2046.83 6.5077 2012.69 6.5127 1617.16 6.5177 1605.77 6.5227 1474.4 6.5277 1494.03 6.5327 1528.2 6.5377 2371.11 6.5427 2342.87 6.5477 2462.19 6.5527 1199.02 6.5577 1338.16 6.5627 1848.08 6.5677 1417.69 6.5727 1341.52 6.5777 1912.57 6.5827 2063.28 6.5877 2130.67 6.5927 1340.41 6.5977 1399.35 6.6027 1504.75 6.6077 2122.97 6.6127 2058.16 6.6177 2077.44 6.6227 2108.84 6.6277 2192.75 6.6327 1457 6.6377 1273.54 6.6427 2226.83 6.6477 1219.01 6.6527 1110.91 6.6577 2123.42 6.6627 1936.94 6.6677 2854.69 6.6727 2221.59 6.6777 1328.31 6.6827 1438.65}

    canvas .c -width 500 -height 400
    pack .c 
    emu_graph foo -canvas .c -width 400 -height 300
    foo data d1 -colour red -points 1 -lines 0 \
  -coords $testdata(1)
    foo data d2 -colour blue -points 0 -lines 1 \
  -coords $testdata(2)
} 

proc tdtest {} {
    emutemplate T andosl
    T trackdata tdata jc:c:mscjc001 rms 0 0 
 
    emu_graph foo -canvas .c
    foo data d1 -colour red -points 0 -lines 1 -trackdata tdata
#    foo data d2 -colour blue -points 0 -lines 1 -coords [tdata coords]
}

