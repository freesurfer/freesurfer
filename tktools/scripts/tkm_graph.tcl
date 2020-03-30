##
## tkm_graph.tcl - functions for use with BLT graph
##
##
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


package require BLT;

proc grf_New { iwGraph isTitle } {

    blt::graph $iwGraph -title $isTitle \
	    -plotbackground white
 
    $iwGraph legend bind all <Enter> { 
	$iwGraph element configure [$iwGraph legend get current] \
		    -linewidth $knLineWidth(active)
	$iwGraph legend activate [$iwGraph legend get current]
    }
    
    $iwGraph legend bind all <Leave> { 
	$iwGraph element configure [$iwGraph legend get current] \
		-linewidth $knLineWidth(inactive)
	$iwGraph legend deactivate [$iwGraph legend get current]
    }

    bind $iwGraph <ButtonPress-2> { grf_RegionStart %W %x %y }
    bind $iwGraph <B2-Motion> { grf_RegionMotion %W %x %y }
    bind $iwGraph <ButtonRelease-2> { grf_RegionEnd %W %x %y }
    bind $iwGraph <ButtonRelease-3> { grf_Unzoom %W }
}

proc grf_Zoom { iwGraph inX1 inY1 inX2 inY2 } {

    if { $inX1 < $inX2 } {
	$iwGraph axis configure x -min $inX1 -max $inX2
    } elseif { $inX1 > $inX2 } {
	$iwGraph axis configure x -min $inX2 -max $inX1
    }
    if { $inY1 < $inY2 } {
	$iwGraph axis configure y -min $inY1 -max $inY2
    } elseif { $inY1 > $inY2 } {
	$iwGraph axis configure y -min $inY2 -max $inY1
    }
}

proc grf_Unzoom { iwGraph } {
    $iwGraph axis configure x y -min {} -max {}
}

proc grf_RegionStart { iwGraph inX inY } {
    global gnRegionStart
    $iwGraph marker create line -coords { } -name zoomBox \
	    -dashes dash -xor yes
    set gnRegionStart(x) [$iwGraph axis invtransform x $inX]
    set gnRegionStart(y) [$iwGraph axis invtransform y $inY]
}

proc grf_RegionMotion { iwGraph inX inY } {
    global gnRegionStart
    set nX [$iwGraph axis invtransform x $inX]
    set nY [$iwGraph axis invtransform y $inY]
    $iwGraph marker configure zoomBox -coords [list \
	    $gnRegionStart(x) $gnRegionStart(y) \
	    $gnRegionStart(x) $nY $nX $nY \
	    $nX $gnRegionStart(y) \
	    $gnRegionStart(x) $gnRegionStart(y)]
}

proc grf_RegionEnd { iwGraph inX inY } {
    global gnRegionStart
    $iwGraph marker delete zoomBox
    set nX [$iwGraph axis invtransform x $inX]
    set nY [$iwGraph axis invtransform y $inY]
    Zoom $iwGraph $gnRegionStart(x) $gnRegionStart(y) $nX $nY
}
