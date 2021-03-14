##
## make_movie.tcl
##
##
## Copyright © 2021
## The General Hospital Corporation (Boston, MA). 
## All rights reserved.
##
## Distribution, usage and copying of this software is covered under the
## terms found in the License Agreement file named 'COPYING' found in the
## FreeSurfer source code root directory, and duplicated here:
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
##
## General inquiries: freesurfer@nmr.mgh.harvard.edu
## Bug reports: analysis-bugs@nmr.mgh.harvard.edu
##

set err [catch {

		::tkcon_tcl_puts "Starting..."

		set prefix I13_BWH_in_vivo_angio_nl_dist

    # Load an LUT file for the segementation to use. Get its ID. Set
    # its label to something good.
    set lutID [MakeNewColorLUT]
		::tkcon_tcl_puts "Loading hires LUT..."
#    SetColorLUTFileName $lutID [file join $env(FREESURFER_HOME) jeans_labels.txt]
    SetColorLUTFileName $lutID [file join $env(FREESURFER_HOME) FreeSurferColorLUT.txt]
    SetColorLUTLabel $lutID "Hires LUT"
		UpdateLUTList

		::tkcon_tcl_puts "Loading lowres LUT..."
    set lowreslutID [MakeNewColorLUT]
    SetColorLUTFileName $lowreslutID [file join $env(FREESURFER_HOME) tkmeditColorsCMA]
    SetColorLUTLabel $lowreslutID "Lowres LUT"
		UpdateLUTList
    
    # Load base level. The LoadVolume commmand makes a new layer in
    # it, and automatically adds it to the current view at the first
    # available draw level, 0. Give it a label.
    set fnBaseVolume /space/rock/4/users/hires/recons/I13_BWH_in_vivo/mri/norm.mgz
		::tkcon_tcl_puts "Loading $fnBaseVolume..."
    set baseLayerID [LoadVolume $fnBaseVolume 1 [GetMainFrameID]]
		::tkcon_tcl_puts "Done."
    SetLayerLabel $baseLayerID "Base volume"


    # Load lowres seg level. The LoadVolume commmand makes a new layer in
    # it, and automatically adds it to the current view at the first
    # available draw level, 0. Give it a label.
    set fnSegVol /space/rock/4/users/hires/recons/I13_BWH_in_vivo/mri/aseg.mgz
		::tkcon_tcl_puts "Loading $fnSegVol..."
    set lowresSegLayerID [LoadVolume $fnSegVol 1 [GetMainFrameID]]
    SetLayerLabel $lowresSegLayerID "lowres segmentation"
		::tkcon_tcl_puts "Done."
    set lowrescolID [Get2DMRILayerVolumeCollection $lowresSegLayerID]
    Set2DMRILayerColorMapMethod $lowresSegLayerID lut
    SetLayerOpacity $lowresSegLayerID 0.3
    Set2DMRILayerDrawZeroClear $lowresSegLayerID true
		Set2DMRILayerColorLUT $lowresSegLayerID $lowreslutID
    SetLayerLabel $lowresSegLayerID "Lowres Segmentation Overlay"


    # Load the first segmentation volume, making a layer for it and
    # adding to the next available draw level, 1. Save the layer
    # ID. Also get the data collection ID as we will use that later to
    # change the volume data to which the collection is pointing.
    set nVol 0
#    set fnSegVol [format "%s_%04d.mgz" $prefix $nVol]
    set fnSegVol [format "%s_%04d.mgz" $prefix $nVol]
		::tkcon_tcl_puts "Loading $fnSegVol..."
    set layerID [LoadVolume $fnSegVol 1 [GetMainFrameID]]
		::tkcon_tcl_puts "Done."
    set colID [Get2DMRILayerVolumeCollection $layerID]


    # Set up the layer as a segmentation volume: lut color map, draw
    # zero values clear, and assign the LUT we got earlier. Give it a
    # label.
    Set2DMRILayerColorMapMethod $layerID lut
    Set2DMRILayerDrawZeroClear $layerID true
    Set2DMRILayerColorLUT $layerID $lutID
    SetLayerLabel $layerID "Segmentation Overlay"


    # Set a location in the view.
    SetViewRASCenter 0 28 22 -40
    SetViewRASCenter 0 25.88 11.00 -35.62
    SetViewRASCenter 0 31.38 29.00 -40.38
    SetViewRASCenter 0 27.75 36.50 -15.24
    SetViewZoomLevel 0 8
		SetViewInPlane 0 y
		set gaView(current,inPlane) y
		UpdateFrame [GetMainFrameID]


    # For each volume...
    for { set nVol 0 } { $nVol <= 13 } { incr nVol 1 } {
        
        # Generate a volume file name using our base name and the number.
        set fnSegVol [format "%s_%04d.mgz" $prefix $nVol]
        # Set the volume file name in the collection to this new
        # volume and tell the collection to load the volume.
        SetVolumeCollectionFileName $colID $fnSegVol
				::tkcon_tcl_puts "Loading $fnSegVol..."
        SetStatusBarText "Loading $fnSegVol..."
        LoadVolumeFromFileName $colID
				::tkcon_tcl_puts "Done."
        # Force a window update.
        UpdateFrame [GetMainFrameID]
        # Take a screen shot.
        set fnCapture [format "%s-%04d.tiff" $prefix $nVol]
				::tkcon_tcl_puts "writing image to $fnCapture..."
        CaptureFrameToFile [GetMainFrameID] $fnCapture
    }


} sResult]

# Check for errors.
if { 0 != $err } { 
    ::tkcon_tcl_puts "Script failed: $sResult"
} else { 
    ::tkcon_tcl_puts "Script complete." 
}
