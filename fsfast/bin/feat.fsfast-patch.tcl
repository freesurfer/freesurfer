# feat.fsfast-patch.tcl 
# 
# This is a modified version of FSL's FEAT TCL program. 
# All restrictions of the FSL LICENSE below apply.
#
# $Id: feat.fsfast-patch.tcl,v 1.3 2007/09/16 23:17:42 greve Exp $

#{{{ copyright and setup 

#   FEAT - FMRI Expert Analysis Tool
#
#   Stephen Smith, FMRIB Image Analysis Group
#
#   Copyright (C) 1999-2006 University of Oxford
#
#   Part of FSL - FMRIB's Software Library
#   http://www.fmrib.ox.ac.uk/fsl
#   fsl@fmrib.ox.ac.uk
#   
#   Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
#   Imaging of the Brain), Department of Clinical Neurology, Oxford
#   University, Oxford, UK
#   
#   
#   LICENCE
#   
#   FMRIB Software Library, Release 3.3 (c) 2006, The University of
#   Oxford (the "Software")
#   
#   The Software remains the property of the University of Oxford ("the
#   University").
#   
#   The Software is distributed "AS IS" under this Licence solely for
#   non-commercial use in the hope that it will be useful, but in order
#   that the University as a charitable foundation protects its assets for
#   the benefit of its educational and research purposes, the University
#   makes clear that no condition is made or to be implied, nor is any
#   warranty given or to be implied, as to the accuracy of the Software,
#   or that it will be suitable for any particular purpose or for use
#   under any specific conditions. Furthermore, the University disclaims
#   all responsibility for the use which is made of the Software. It
#   further disclaims any liability for the outcomes arising from using
#   the Software.
#   
#   The Licensee agrees to indemnify the University and hold the
#   University harmless from and against any and all claims, damages and
#   liabilities asserted by third parties (including claims for
#   negligence) which arise directly or indirectly from the use of the
#   Software or the sale of any products based on the Software.
#   
#   No part of the Software may be reproduced, modified, transmitted or
#   transferred in any form or by any means, electronic or mechanical,
#   without the express permission of the University. The permission of
#   the University is not required if the said reproduction, modification,
#   transmission or transference is done without financial return, the
#   conditions of this Licence are imposed upon the receiver of the
#   product, and all original and amended source code is included in any
#   transmitted product. You may be held legally responsible for any
#   copyright infringement that is caused or encouraged by your failure to
#   abide by these terms and conditions.
#   
#   You are not permitted under this Licence to use this Software
#   commercially. Use for which any financial return is received shall be
#   defined as commercial use, and includes (1) integration of all or part
#   of the source code or the Software into a product for sale or license
#   by or on behalf of Licensee to third parties or (2) use of the
#   Software or any derivative of it for research with the final aim of
#   developing software products for sale or license to a third party or
#   (3) use of the Software or any derivative of it for research with the
#   final aim of developing non-software products for sale or license to a
#   third party, or (4) use of the Software to provide any service to an
#   external organisation for which payment is received. If you are
#   interested in using the Software commercially, please contact Isis
#   Innovation Limited ("Isis"), the technology transfer company of the
#   University, to negotiate a licence. Contact details are:
#   innovation@isis.ox.ac.uk quoting reference DE/1112.

source ${FSLDIR}/tcl/fslstart.tcl

set VARS(history) {}

#}}}

##### FEAT GUI #####
#{{{ feat5 GUI

proc feat5 { w } {
    global fmri PXHOME FSLDIR USER feat_files unwarp_files unwarp_files_mag initial_highres_files highres_files VARS argc argv PWD gui_ext HOME
 
    #{{{ main window

feat5:setupdefaults

toplevel $w

wm title      $w "FEAT - FMRI Expert Analysis Tool - v$fmri(version)"
wm iconname   $w "FEAT [ expr int($fmri(version)) ]"
wm iconbitmap $w @${FSLDIR}/tcl/fmrib.xbm

tixBalloon $w.bhelp -initwait 3000 -options { background #b0ffb0 wrapLength 800 }

#}}}
    #{{{ mode

frame $w.mode

tixOptionMenu $w.mode.level -variable fmri(level)
$w.mode.level add command 1 -label "First-level analysis"
$w.mode.level add command 2 -label "Higher-level analysis"

$w.bhelp bind $w.mode.level -msg "Use \"First-level analysis\" for analysing each session's data -
i.e. the time-series analysis of the raw 4D FMRI data.

Use \"Higher-level analysis\" for combining first-level analyses. You
can use this hierarchically - for example at second-level to analyse
across several sessions and then at third-level to analyse across
several subjects."

tixOptionMenu $w.mode.analysis -variable fmri(analysis)
$w.mode.analysis add command 7 -label "Full analysis"
$w.mode.analysis add command 1 -label "Pre-stats"
$w.mode.analysis add command 3 -label "Pre-stats + Stats"
$w.mode.analysis add command 2 -label "                     Stats"
$w.mode.analysis add command 6 -label "                     Stats + Post-stats"
$w.mode.analysis add command 4 -label "                                  Post-stats"
$w.mode.analysis add command 0 -label "Registration only"

$w.bhelp bind $w.mode.analysis -msg "You can run a full analysis - Pre-Stats; Stats; Post-stats;
Registration - or a (sensible) subset of these options.

If you select \"Post-stats\" or \"Registration only\", you will need
to select a FEAT directory (or directories) instead of starting with
4D image data; the results already produced in those FEAT directories
will then be used as appropriate.

Note that if you want to run only \"Post-stats\", you must select the
FEAT directory/directories before editing the contrasts or
thresholding parameters, as these will get reset on selection of the
FEAT directory/directories."

pack $w.mode.level $w.mode.analysis -in $w.mode -side left -anchor w

#}}}
    #{{{ main notebook

tixNoteBook $w.nb
$w.nb add misc      -label "Misc"
$w.nb add data      -label "Data"
$w.nb add filtering -label "Pre-stats"
$w.nb add stats     -label "Stats"
$w.nb add poststats -label "Post-stats"
$w.nb add reg       -label "Registration"

#{{{ Misc

set fmri(miscf) [ $w.nb subwidget misc ]

#{{{ balloon help

checkbutton $w.help -variable fmri(help_yn) -text "Balloon help" -command "feat5:updatehelp $w" -justify right

$w.bhelp bind $w.help -msg "And don't expect this message to appear whilst you've turned this
option off!"

feat5:updatehelp $w

#}}}
#{{{ featwatcher

checkbutton $w.featwatcher -variable fmri(featwatcher_yn) -text "Featwatcher" -justify right

$w.bhelp bind $w.featwatcher -msg "Start the Featwatcher GUI whenever FEAT is run?"

#}}}
#{{{ delay

tixControl $w.delay -label "Delay before starting (hours) " \
	-variable fmri(delay) -step 1 -min 0 -selectmode immediate

$w.bhelp bind $w.delay -msg "If you are using FEAT to carry out multiple analyses, you might want
to do this overnight, to reduce the load on your computer. You can use
this delay feature to tell FEAT how long to wait before starting the
analyses." 

#}}}
#{{{ brain threshold

tixControl $w.brain_thresh -label "Brain/background threshold, % " \
	-variable fmri(brain_thresh) -step 1 -min 0 -max 100 -selectmode immediate

$w.bhelp bind $w.brain_thresh -msg "This is automatically calculated, as a % of the maximum input image
intensity. It is used in intensity normalisation, brain mask
generation and various other places in the FEAT analysis."

#}}}
#{{{ design efficiency

tixLabelFrame $w.contrastest -label "Design efficiency"
set contrastestsub [ $w.contrastest subwidget frame ]

tixControl $w.contrastest.noise -label "Noise level % " \
    -variable fmri(noise) -step .25 -min 0.0001 -max 100 -selectmode immediate
grid $w.contrastest.noise -in $contrastestsub -column 0 -row 0 -padx 3 -pady 3

tixControl $w.contrastest.noisear -label "Temporal smoothness " \
    -variable fmri(noisear) -step .1 -min -0.99 -max 0.99 -selectmode immediate
grid $w.contrastest.noisear -in $contrastestsub -column 1 -row 0 -padx 3 -pady 3

button $w.contrastest.estnoise -text "Estimate from data" -command "feat5:estnoise"
grid $w.contrastest.estnoise -in $contrastestsub -column 0 -row 1 -padx 3 -pady 3

tixControl $w.contrastest.critical_z -label "Z threshold " \
	-variable fmri(critical_z) -step 1 -min 0.0001 -max 100 -selectmode immediate
$w.bhelp bind $w.contrastest.critical_z -msg "This is the Z value used to determine what level of activation would
be statistically significant, to be used only in the design
efficiency calculation. Increasing this will result in higher
estimates of required effect."

grid $w.contrastest.critical_z -in $contrastestsub -column 1 -row 1 -padx 3 -pady 3

$w.bhelp bind $w.contrastest -msg "The \"Noise level %\" and \"Temporal smoothness\" together
characterise the noise in the data, to be used only in the design
efficiency estimation.

The \"Noise level %\" is the standard deviation (over time) for a
typical voxel, expressed as a percentage of the baseline signal level.

The \"Temporal smoothness\" is the smoothness coefficient in a simple
AR(1) autocorrelation model (much simpler than that actually used in
the FILM timeseries analysis but good enough for the efficiency
calculation here).

If you want to get a rough estimate of this noise level and temporal
smoothness from your actual input data, press the \"Estimate from
data\" button (after you have told FEAT where your input data
is). This takes about 30-60 seconds to estimate. This applies just the
spatial and temporal filtering (i.e., no motion correction) that you
have specified in the \"Pre-stats\" section, and gives a reasonable
approximation of the noise characteristics that will remain in the
fully preprocessed data, once FEAT has run."

#}}}
#{{{ relative filenames for custom timings and highres images

checkbutton $w.relative -variable fmri(relative_yn) -text "Use relative filenames for multiple analyses" -command "feat5:updatereg $w"

$w.bhelp bind $w.relative -msg "blah...."

#}}}
#{{{ new directory if re-thresholding

set tmpval $fmri(newdir_yn)
tixOptionMenu $w.newdir_yn -variable fmri(newdir_yn)
$w.newdir_yn add command 0 -label "Overwrite original post-stats results"
$w.newdir_yn add command 1 -label "Copy original FEAT directory for new Post-stats / Registration"
set fmri(newdir_yn) $tmpval

$w.bhelp bind $w.newdir_yn -msg "If you are just re-running post-stats or registration, you can either
choose to overwrite the original post-stats and registration results
or create a complete copy of the original FEAT directory, with the new
results in it."

#}}}
#{{{ cleanup first-level standard-space data

checkbutton $w.sscleanup -variable fmri(sscleanup_yn) -text "Cleanup first-level standard-space images" -justify right

$w.bhelp bind $w.sscleanup -msg "When you run a higher-level analysis, the first thing that happens is
that first-level images are transformed into standard-space (in
<firstlevel>.feat/reg_standard subdirectories) for feeding into the
higher-level analysis. This takes up quite a lot of disk space, so if
you want to save disk space, turn this option on and these these
upsampled images will get deleted once they have been fed into the
higher-level analysis. However, generating them can take quite a lot
of time, so if you want to run several higher-level analyses, all
using the same first-level FEAT directories, then leave this option
turned off."

#}}}

pack $w.help $w.featwatcher $w.delay $w.brain_thresh $w.contrastest -in $fmri(miscf) -anchor w -side top -padx 5 -pady 1

#}}}
#{{{ Data

set fmri(dataf) [ $w.nb subwidget data ]

#{{{ input type for higher-level

tixOptionMenu $w.inputtype -variable fmri(inputtype)
$w.inputtype add command 1 -label "Inputs are lower-level FEAT directories"
$w.inputtype add command 2 -label "Inputs are 3D cope images from FEAT directories"

$w.inputtype configure  -command "feat5:updateselect $w"

$w.bhelp bind $w.inputtype -msg "Select the kind of input you want to feed into the higher-level FEAT
analysis.

If you choose to select \"FEAT directories\", then the higher-level design
will get applied across the selected FEAT directories; each
lower-level FEAT directory forms a \"time-point\" in the higher-level
model. For example, each lower-level FEAT directory represents a
single session in a multiple-session higher-level analysis or a single
subject in a multiple-subject analysis. If the lower-level FEAT
directories contain more than one contrast (cope), then the
higher-level analysis is run separately for each one; in this case,
the higher-level \".gfeat\" directory will end up contain more than
one \".feat\" directory, one for each lower-level contrast. This
option requires all lower-level FEAT directories to include the same
set of contrasts.

If you choose to select \"3D cope images from FEAT directories\", then
you explicitly control which cope corresponds to which \"time-point\"
to be fed into the higher-level analysis. For example, if you have only
one lower-level FEAT directory, containing multiple contrasts (copes),
and you want to carry out a higher-level analysis across these
contrasts, then this is the correct option. With this option (as with
the previous one), the chosen copes will automatically get transformed
into standard space if necesary."

#}}}

#{{{ multiple analyses

frame $w.multiple

tixControl $w.multiple.number -label "Number of analyses " \
	-variable fmri(multiple) -step 1 -min 1 -selectmode immediate \
	-options { entry.width 3 }

button $w.multiple.setup -text "Select 4D data" -command "feat5:multiple_select $w 0 \"Select input data\" "

pack $w.multiple.number $w.multiple.setup -in $w.multiple -side left -padx 5

$w.bhelp bind $w.multiple -msg "Set the filename of the 4D input image (e.g. /home/sibelius/func.hdr).
You can setup FEAT to process many input images, one after another, as
long as they all require exactly the same analysis. Each one will
generate its own FEAT directory, the name of which is based on the
input data's filename.
Alternatively, if you are running either just \"Post-stats\" or
\"Registration only\", or running \"Higher-level analysis\", the
selection of 4D data changes to the selection of FEAT directories.
Note that in this case you should select the FEAT directories before
setting up anything else in FEAT (such as changing the
thresholds). This is because quite a lot of FEAT settings are loaded
from the first selected FEAT directory, possibly over-writing any
settings which you wish to change!"

$w.multiple.number configure -command "feat5:updateselect $w"

#}}}

#{{{ output directory

FSLFileEntry $w.outputdir \
	-variable fmri(outputdir) \
	-directory $PWD \
	-label " Output directory" \
	-title "Name the output directory" \
	-width 35 \
	-filterhist VARS(history)

$w.bhelp bind $w.outputdir -msg "If this is left blank, the output FEAT directory name is derived from
the input data name. (For higher-level analysis, the output name is
derived from the first lower-level FEAT directory selected as input.)

If, however, you wish to explicitly choose the output FEAT directory
name, for example, so that you can include in the name a hint about
the particular analysis that was carried out, you can set this
here.

This output directory naming behaviour is modified if you are setting
up multiple first-level analyses, where you are selecting multiple
input data sets and will end up with multiple output FEAT
directories. In this case, whatever you enter here will be used and
appended to what would have been the default output directory name if
you had entered nothing. For example, if you are setting up 3 analyses
with input data names \"/home/neo/fmri1.hdr\", \"/home/neo/fmri2.hdr\"
and \"/home/neo/fmri3.hdr\", and set the output name to \"analysisA\",
the output directories will end up as
\"/home/neo/fmri1_analysisA.feat\" etc."

#}}}

frame $w.datamain
#{{{ npts & ndelete

frame $w.nptsndelete

#{{{ npts

set fmri(npts) 0

tixControl $w.npts -label "Total volumes " -variable fmri(npts) \
	-step 1 -min 0 -selectmode immediate -options { entry.width 4 }

$w.bhelp bind $w.npts -msg "The number of FMRI volumes in the time series, including any initial
volumes that you wish to delete. This will get set automatically once
valid input data has been selected.

Alternatively you can set this number by hand before selecting data so
that you can setup and view a model without having any data, for
experimental planning purposes etc."

#}}}
#{{{ ndelete

tixControl $w.ndelete -label "       Delete volumes " -variable fmri(ndelete) \
	-step 1 -min 0 -selectmode immediate -options { entry.width 3 }

$w.bhelp bind $w.ndelete -msg "The number of initial FMRI volumes to delete before any further
processing. You should have decided on this number when the scans were
acquired.  Typically your experiment would have begun after these
initial scans (sometimes called \"dummy scans\"). These should be the
volumes that are not wanted because steady-state imaging has not yet
been reached - typically two or three volumes. These volumes are
deleted as soon as FEAT is started, so all 4D data sets produced by
FEAT will not contain the deleted volumes.

Note that \"Delete volumes\" should not be used to correct for the
time lag between stimulation and the measured response - this is
corrected for in the design matrix by convolving the input stimulation
waveform with a blurring-and-delaying haemodynamic response function.

Most importantly, remember when setting up the design matrix, that the
timings in the design matrix start at t=0 seconds, and this
corresponds to the start of the first image taken after the deleted
scans. In other words, the design matrix starts AFTER the \"deleted
scans\" have been deleted."

#}}}

pack $w.npts $w.ndelete -in $w.nptsndelete -side left

#}}}
#{{{ TR & highpass

frame $w.trparadigm_hp

#{{{ TR

tixControl $w.tr -label "TR (s) " -variable fmri(tr) \
	-step 0.25 -min 0.001 -selectmode immediate -options { entry.width 4 }
$w.bhelp bind $w.tr -msg "The time (in seconds) between scanning successive FMRI volumes."

#}}}
#{{{ High pass

tixControl $w.paradigm_hp -label "     High pass filter cutoff (s) " \
	-variable fmri(paradigm_hp) -step 5 -min 1 -selectmode immediate -options { entry.width 5 }

$w.bhelp bind $w.paradigm_hp -msg "The high pass frequency cutoff point (seconds), that is, the longest
temporal period that you will allow.

A sensible setting in the case of an rArA or rArBrArB type block
design is the (r+A) or (r+A+r+B) total cycle time.

For event-related designs the rule is not so simple, but in general
the cutoff can typically be reduced at least to 50s.

This value is setup here rather than in Pre-stats because it also
affects the generation of the model; the same high pass filtering is
applied to the model as to the data, to get the best possible match
between the model and data."

#}}}

pack $w.tr $w.paradigm_hp -in $w.trparadigm_hp -side left

#}}}

pack $w.nptsndelete $w.trparadigm_hp -in $w.datamain -side top -padx 5 -pady 3 -anchor w

pack $w.multiple $w.datamain -in $fmri(dataf) -anchor w -side top

#{{{ FSL logo

set graphpic [ image create photo -file ${FSLDIR}/tcl/fsl-logo-tiny.ppm ]
button $w.logo -image $graphpic -command "FmribWebHelp file: ${FSLDIR}/doc/index.html" -borderwidth 0
pack $w.logo -in $fmri(dataf) -anchor e -side bottom -padx 5 -pady 5

#}}}

#}}}
#{{{ Pre-statistics processing

set fmri(filteringf) [ $w.nb subwidget filtering ]

#{{{ motion correction

frame $w.mc
set fmri(mcf) $w.mc

set tmpval $fmri(mc)
tixOptionMenu $w.mc.menu -label "Motion correction: " -variable fmri(mc)
$w.mc.menu add command 0 -label "None"
$w.mc.menu add command 1 -label "MCFLIRT"
set fmri(mc) $tmpval

pack $w.mc.menu -in $fmri(mcf) -side top -side left

$w.bhelp bind $w.mc -msg "You will normally want to apply motion correction; this attempts to
remove the effect of subject head motion during the
experiment. MCFLIRT uses FLIRT (FMRIB's Linear Registration Tool)
tuned to the problem of FMRI motion correction, applying rigid-body
transformations.

Note that there is no \"spin history\" (aka \"correction for movement\")
option with MCFLIRT. This is because this is still a poorly understood
correction method which is under further investigation."

#}}}
#{{{ B0 unwarping

frame $w.unwarpf
set fmri(unwarpf) $w.unwarpf

label $fmri(unwarpf).label -text "B0 unwarping"

checkbutton $fmri(unwarpf).yn -variable fmri(regunwarp_yn) -command "feat5:updateprestats $w 0"

tixLabelFrame $fmri(unwarpf).lf -label "B0 unwarping"
set fmri(unwarpff) [ $fmri(unwarpf).lf subwidget frame ]

set unwarp_files(1) "~"
set unwarp_files_mag(1) "~"

FSLFileEntry $fmri(unwarpff).unwarpsingle \
	-variable unwarp_files(1) \
	-pattern "IMAGE" \
	-directory $PWD \
	-label "Fieldmap       " \
	-title "Select the B0 fieldmap image" \
	-width 35 \
	-filterhist VARS(history) \
	-command "feat5:multiple_check $w 1 1 0"

button $fmri(unwarpff).unwarpmultiple -text "Select the B0 fieldmap images" \
	-command "feat5:multiple_select $w 1 \"Select the B0 fieldmap images\" "

FSLFileEntry $fmri(unwarpff).unwarpmagsingle \
	-variable unwarp_files_mag(1) \
	-pattern "IMAGE" \
	-directory $PWD \
	-label "Fieldmap mag" \
	-title "Select the B0 fieldmap magnitude image" \
	-width 35 \
	-filterhist VARS(history) \
	-command "feat5:multiple_check $w 2 1 0"

button $fmri(unwarpff).unwarpmagmultiple -text "Select the B0 fieldmap magnitude images" \
	-command "feat5:multiple_select $w 2 \"Select the B0 fieldmap magnitude images\" "

frame $fmri(unwarpff).opts1

tixControl $fmri(unwarpff).opts1.dwell -label "EPI dwell time (ms)" -variable fmri(dwell) \
	-step 0.1 -options { entry.width 5 }

tixControl $fmri(unwarpff).opts1.te -label "  EPI TE (ms)" -variable fmri(te) \
	-step 1 -options { entry.width 5 }

frame $fmri(unwarpff).opts2

tixControl $fmri(unwarpff).opts2.signallossthresh -label "  % Signal loss threshold" -variable fmri(signallossthresh) \
	-step 1 -selectmode immediate -options { entry.width 3 }

tixOptionMenu $fmri(unwarpff).opts2.unwarp_dir -label "Unwarp direction " -variable fmri(unwarp_dir)
$fmri(unwarpff).opts2.unwarp_dir add command x  -label "x"
$fmri(unwarpff).opts2.unwarp_dir add command x- -label "-x"
$fmri(unwarpff).opts2.unwarp_dir add command y  -label "y"
$fmri(unwarpff).opts2.unwarp_dir add command y- -label "-y"
$fmri(unwarpff).opts2.unwarp_dir add command z  -label "z"
$fmri(unwarpff).opts2.unwarp_dir add command z- -label "-z"
set fmri(unwarp_dir) y-

pack $fmri(unwarpff).opts1.dwell $fmri(unwarpff).opts1.te -in $fmri(unwarpff).opts1 -side left -anchor s
pack $fmri(unwarpff).opts2.unwarp_dir $fmri(unwarpff).opts2.signallossthresh -in $fmri(unwarpff).opts2 -side left -anchor s

pack $fmri(unwarpff).unwarpsingle $fmri(unwarpff).opts1 $fmri(unwarpff).opts2 -in $fmri(unwarpff) -anchor w -side top -pady 2 -padx 3

pack $fmri(unwarpf).label $fmri(unwarpf).yn -in $fmri(unwarpf) -side left

$w.bhelp bind $fmri(unwarpf) -msg "B0 unwarping is carried out using FUGUE. Here you need to enter the B0 fieldmap images which should be created before you run FEAT and usually require site/scanner/sequence specific processing.  See the PRELUDE/FUGUE documentation for more information on creating these images.  The two images that are required are (1) a fieldmap image which must have units of rad/s, and (2) a registered magnitude image (this is usually a standard magnitude-only reconstructed image from the fieldmap sequence data).

Next you need to enter the \"EPI Dwell time\" in milliseconds (this is the time between echoes in successive k-space lines - often also known as the echo spacing for EPI) and \"EPI TE\" (echo time) also in milliseconds. These two values relate to your FMRI EPI data, not the fieldmap data.

You also need to specify the \"Unwarp direction\", which is the phase-encoding direction of your FMRI EPI data. The sign of this direction will depend on both the sign of the phase encode blips in the EPI sequence and on the sign of the fieldmap.  As it can be difficult to predict this sign when using a particular site/scanner/sequence for the first time, it is usual to try both positive and negative values in turn and see which gives better undistortion (the wrong sign will increase the amount of distortion rather than decrease it).

Finally, you need to specify a \"% Signal loss threshold\". This determines where the signal loss in the EPI is too great for registration to get a good match between the EPI data and other images. Areas where the % signal loss in the EPI exceeds this threshold will get masked out of the registration process between the EPI and the fieldmap and structural images.

If you are running both motion correction and B0 unwarping, the motion correction resampling does not get applied at the same time as the motion estimation; instead the motion correction gets applied simultaneously with the application of the B0 unwarping, in order to minimise interpolation-related image blurring.

Once you have run FEAT you should definitely check the unwarping report (click on the mini-movie, shown in the main FEAT report page, that flicks between distorted and undistorted versions of example_func). In particular you should check that it looks like the unwarping has occurred in the correct direction (and change the unwarp direction and/or sign if it is not)."

#}}}
#{{{ slice timing correction

frame $w.st
set fmri(stf) $w.st

FSLFileEntry $w.st_file \
	-variable fmri(st_file) \
	-pattern "*" \
	-directory $PWD \
	-label "" \
	-title "Select a slice order/timings file" \
	-width 20 \
	-filterhist VARS(history)

set tmpval $fmri(st)
tixOptionMenu $w.st.menu -label "Slice timing correction: " -variable fmri(st) -command "feat5:updateprestats $w"
$w.st.menu add command 0 -label "None"
$w.st.menu add command 1 -label "Regular up (0, 1, 2 ... n-1)"
$w.st.menu add command 2 -label "Regular down (n-1, n-2 ... 0)"
$w.st.menu add command 5 -label "Interleaved (0, 2, 4 ... 1, 3, 5 ... )"
$w.st.menu add command 3 -label "Use slice order file"
$w.st.menu add command 4 -label "Use slice timings file"
set fmri(st) $tmpval

pack $w.st.menu -in $fmri(stf) -side top -side left

$w.bhelp bind $w.st -msg "Slice timing correction corrects each voxel's time-series for the fact
that later processing assumes that all slices were acquired exactly
half-way through the relevant volume's acquisition time (TR), whereas
in fact each slice is taken at slightly different times.

Slice timing correction works by using (Hanning-windowed) sinc
interpolation to shift each time-series by an appropriate fraction of
a TR relative to the middle of the TR period. It is necessary to know
in what order the slices were acquired and set the appropriate option
here.

If slices were acquired from the bottom of the brain to the top select 
\"Regular up\".  If slices were acquired from the top of the brain
to the bottom select \"Regular down\".

If the slices were acquired with interleaved order (0, 2, 4 ... 1, 3,
5 ...) then choose the \"Interleaved\" option.

If slices were not acquired in regular order you will need to  use a
slice order file or a slice timings file. If a slice order file is to
be used, create a text file with a single number on each line, where
the first line states which slice was acquired first, the second line
states which slice was acquired second, etc. The first slice is
numbered 1 not 0.

If a slice timings file is to be used, put one value (ie for each
slice) on each line of a text file. The units are in TRs, with 0.5
corresponding to no shift. Therefore a sensible range of values will
be between 0 and 1."

#}}}
#{{{ spin history

#frame $w.sh
#
#label $w.sh.label -text "Adjustment for movement"
#
#checkbutton $w.sh.yn -variable fmri(sh_yn)
#
#pack $w.sh.label $w.sh.yn -in $w.sh -padx 5 -side left
#
#$w.bhelp bind $w.sh -msg "blah"

#}}}
#{{{ BET brain extraction

frame $w.bet

label $w.bet.label -text "BET brain extraction"

checkbutton $w.bet.yn -variable fmri(bet_yn)

pack $w.bet.label $w.bet.yn -in $w.bet -side left

$w.bhelp bind $w.bet -msg "This uses BET brain extraction to create a brain mask from the first
volume in the FMRI data. This is normally better than simple
intensity-based thresholding for getting rid of unwanted voxels in
FMRI data. Note that here, BET is setup to run in a quite liberal way so that
there is very little danger of removing valid brain voxels.

If the field-of-view of the image (in any direction) is less than 30mm
then BET is turned off by default.

Note that, with respect to any structural image(s) used in FEAT
registration, you need to have already run BET on those before running
FEAT."

#}}}
#{{{ spatial filtering

tixControl $w.smooth -label "Spatial smoothing FWHM (mm) " \
	-variable fmri(smooth) -step 1 -min 0 -selectmode immediate

$w.bhelp bind $w.smooth -msg "This determines the extent of the spatial smoothing, carried out on
each volume of the FMRI data set separately. This is intended to
reduce noise without reducing valid activation; this is successful as
long as the underlying activation area is larger than the extent of
the smoothing. Thus if you are looking for very small activation areas
then you should maybe reduce smoothing from the default of 5mm, and if
you are looking for larger areas, you can increase it, maybe to 10 or
even 15mm.

To turn off spatial smoothing simply set FWHM to 0."

#}}}
#{{{ intensity normalization

frame $w.norm

label $w.norm.label -text "Intensity normalization"

checkbutton $w.norm.yn -variable fmri(norm_yn)

pack $w.norm.label $w.norm.yn -in $w.norm -side left

$w.bhelp bind $w.norm -msg "This forces every FMRI volume to have the same mean intensity. For
each volume it calculates the mean intensity and then scales the
intensity across the whole volume so that the global mean becomes a
preset constant. This step is normally discouraged - hence is turned
off by default. When this step is not carried out, the whole 4D data
set is still normalised by a single scaling factor (\"grand mean
scaling\") - each volume is scaled by the same amount. This is so that
higher-level analyses are valid."

#}}}
#{{{ temporal filtering

frame $w.temp

label $w.temp.label -text "Temporal filtering    "

label $w.temp.pslabel -text "Perfusion subtraction"
checkbutton $w.temp.ps_yn -variable fmri(perfsub_yn) -command "feat5:updateperfusion $w"

set tmpval $fmri(tagfirst)
tixOptionMenu $w.temp.tcmenu -variable fmri(tagfirst)
$w.temp.tcmenu add command 1 -label "First timepoint is tag"
$w.temp.tcmenu add command 0 -label "First timepoint is control"
set fmri(tagfirst) $tmpval

label $w.temp.hplabel -text "Highpass"
checkbutton $w.temp.hp_yn -variable fmri(temphp_yn)

label $w.temp.lplabel -text "Lowpass"
checkbutton $w.temp.lp_yn -variable fmri(templp_yn)

pack $w.temp.label $w.temp.pslabel $w.temp.ps_yn $w.temp.hplabel $w.temp.hp_yn $w.temp.lplabel $w.temp.lp_yn -in $w.temp -side top -side left

$w.bhelp bind $w.temp -msg "\"Perfusion subtraction\" is a pre-processing step for perfusion FMRI
(as opposed to normal BOLD FMRI) data. It subtracts even from odd
timepoints in order to convert tag-control alternating timepoints into
a perfusion-only signal. If you are setting up a full perfusion model
(where you model the full alternating tag/control timeseries in the
design matrix) then you should NOT use this option. The subtraction
results in a temporal shift of the sampled signal to half a TR
earlier; hence you should ideally shift your model forwards in time by
half a TR, for example by reducing custom timings by half a TR or by
increasing the model shape phase by half a TR. When you select this
option, FILM prewhitening is turned off (because it is not
well-matched to the autocorrelation resulting from the subtraction
filter) and instead the varcope and degrees-of-freedom are corrected
after running FILM in OLS mode. See the \"Perfusion\" section of the
manual for more information.


\"Highpass\" temporal filtering uses a local fit of a straight line
(Gaussian-weighted within the line to give a smooth response) to
remove low frequency artefacts. This is preferable to sharp rolloff
FIR-based filtering as it does not introduce autocorrelations into the
data.

\"Lowpass\" temporal filtering reduces high frequency noise by Gaussian
smoothing (sigma=2.8s), but also reduces the strength of the signal of
interest, particularly for single-event experiments. It is not
generally considered to be helpful, so is turned off by default.

By default, the temporal filtering that is applied to the data will also be
applied to the model."

#}}}
#{{{ melodic

frame $w.melodic

label $w.melodic.label -text "MELODIC ICA data exploration"

checkbutton $w.melodic.yn -variable fmri(melodic_yn)

pack $w.melodic.label $w.melodic.yn -in $w.melodic -side top -side left

$w.bhelp bind $w.melodic -msg "This runs MELODIC, the ICA (Independent Component Analysis) tool in
FSL. We recommend that you run this, in order to gain insight into
unexpected artefacts or activation in your data.

You can even use this MELODIC output to \"de-noise\" your data; see
the FEAT manual for information on how to do this."

#}}}

feat5:updateprestats $w 0
pack $fmri(mcf) $fmri(unwarpf) $fmri(stf) $w.bet $w.smooth $w.norm $w.temp $w.melodic -in $fmri(filteringf) -anchor w -pady 1 -padx 5

#}}}
#{{{ Stats

set fmri(statsf) [ $w.nb subwidget stats ]

checkbutton $w.prewhiten -text "Use FILM prewhitening" -variable fmri(prewhiten_yn)
$w.bhelp bind $w.prewhiten -msg "For normal first-level time series analysis you should use
prewhitening to make the statistics valid and maximally efficient. For
other data - for example, very long TR (>30s) FMRI data, PET data or
data with very few time points (<50) - this should be turned off."

set fmri(motionevs) 0
set tmpval $fmri(motionevs)
tixOptionMenu $w.motionevs -label "Add motion parameters to model " -variable fmri(motionevs)
$w.motionevs add command 0 -label "No"
$w.motionevs add command 1 -label "Yes, orthogonalise rest of model wrt motion"
#$w.motionevs add command 2 -label "Yes, orthogonalise motion wrt rest of model"
$w.motionevs add command 3 -label "Yes, don't orthogonalise"
set fmri(motionevs) $tmpval
$w.bhelp bind $w.motionevs -msg "You may want to include the head motion parameters (as estimated by
MCFLIRT motion correction in the Pre-stats processing) as confound EVs
in your model. This can sometimes help remove the residual effects of
motion that are still left in the data even after motion correction.

This is not strongly recommended as there is still much to learn about
residual motion effects; simply adding such confound EVs is quite a
simplistic solution. We would recommend instead turning on MELODIC in
the FEAT Pre-stats and using ICA-based denoising as a better
alternative to removing residual motion effects (see the FEAT manual
for more information on that). However, if you do wish to include
motion parameters in your model then select one of the \"Yes\"
options. If you do this, then once the motion correction has been run,
the translation and rotation parameters are added as extra confound
EVs in the model.

If you do wish to include motion confounds, we recommend using the
\"Yes, orthogonalise rest of model wrt motion\" option. This is the
most conservative approach; all of the signal that is correlated both
to the original EVs and to the motion parameters will be declared as
being artefact, so you are protecting against false positives.

If instead you select the \"Yes, don't orthogonalise\" option, then
signal which correlates with both original EVs and motion parameters
will be modelled by both, in proportion to exactly how well it
correlates with each. This is less conservative than the above option,
in terms of possible false positives."

button $w.wizard -width 20 -text "Model setup wizard" -command "feat5:wizard $w"
$w.bhelp bind $w.wizard -msg "This lets you easily setup simple common experimental designs.

At first level, the options are regular rest-A-rest-A... or
rest-A-rest-B-rest-A-rest-B... designs (block or single-event) for
normal BOLD FMRI, or a rest-A-rest-A... design for full modelling of
perfusion FMRI data.

At second level, the options are one-group t-test, two-group-unpaired
and two-group-paired t-tests.

If you need to further adjust the resulting setup, use \"Model setup
wizard\" first, then press the \"Full model setup\" button."

button $w.model -width 20 -text "Full model setup" -command "feat5:setup_model $w"
$w.bhelp bind $w.model -msg "This allows complete control of the model-based analysis to be used."

set fmri(w_model) 0

set tmpval $fmri(mixed_yn)
tixOptionMenu $w.mixed -variable fmri(mixed_yn)
$w.mixed add command 3 -label "Fixed effects"
$w.mixed add command 0 -label "Mixed effects: Simple OLS"
$w.mixed add command 2 -label "Mixed effects: FLAME 1"
$w.mixed add command 1 -label "Mixed effects: FLAME 1+2"
set fmri(mixed_yn) $tmpval

$w.bhelp bind $w.mixed -msg "The main choice here is between fixed effects (FE) and mixed effects (ME) higher-level modelling. FE modelling is more \"sensitive\" to activation than ME, but is restricted in the inferences that can be made from its results; because FE ignores cross-session/subject variance, reported activation is with respect to the group of sessions or subjects present, and not representative of the wider population. ME does model the session/subject variability, and it therefore allows inference to be made about the wider population from which the sessions/subjects were drawn.

The FE option implements a standard weighted fixed effects model.  No random effects variances are modelled or estimated. The FE error variances are the variances (varcopes) from the previous level. Weighting is introduced by allowing these variances to be unequal (heteroscedastic). Degrees-of-freedom are calculated by summing the effective degrees-of-freedom for each input from the previous level and subtracting the number of higher-level regressors.

We now discuss the different ME options.

OLS (ordinary least squares) is a fast estimation technique which ignores all lower-level variance estimation and applies a very simple higher-level model. This is the least accurate of the ME options.

For the most accurate estimation of higher-level activation you should use FLAME (FMRIB's Local Analysis of Mixed Effects) modelling and estimation. This is a sophisticated two-stage process using Bayesian modelling and estimation (for example it allows separate modelling of the variance in different subject groups, and forces random effects variance to be non-negative).

The first stage of FLAME is significantly more accurate than OLS, and nearly as fast. The second stage of FLAME increases accuracy slightly over the first stage, but is quite a lot slower (typically 45-200 minutes). It takes all voxels which FLAME stage 1 shows to be near threshold and carries out a full MCMC-based analysis at these points, to get the most accurate estimate of activation.

We generally recommend using \"FLAME 1\", as it is MUCH faster than running both stages, and nearly as accurate. 

If you are carrying out a mid-level analysis (e.g., cross-sessions) and will be feeding this into an even higher-level analysis (e.g., cross-subjects), then you should definitely use the \"FLAME 1\" option, as it is not possible for FLAME to know in advance of the highest-level analysis what voxels will ultimately be near threshold. Also, given that the second stage is only run on certain voxels, you should definitely only run stage 1 if you are going to use inference methods (such as mixture modelling) that analyse the unthresholded statistic histogram.

If you do decide to run \"FLAME 1+2\" and the FEAT logs indicate a large difference between the stage 1 and stage 2 estimations (or, for example, the final thresholded zstat image looks \"speckled\"), this is an indication that your data is highly non-Gaussian (e.g., has one or more strong outlier subjects, or has two clearly different groups of subjects being modelled as a single group). In such a case, stage 1 estimation is quite inaccurate (OLS even more so), hence the larger-than-normal difference between stages 1 and 2. The only really good solution is to investigate in your data what is going on - for example, to find the bad outlier."

#}}}
#{{{ Post-Stats

set fmri(poststatsf) [ $w.nb subwidget poststats ]

#{{{ edit contrasts

button $w.modelcon -text "Edit contrasts" -command "feat5:setup_model $w"

$w.bhelp bind $w.modelcon -msg "This allows setup of contrasts and F-tests, to be run on a previous
analysis."

#}}}
#{{{ pre-thresholding masking

FSLFileEntry $fmri(poststatsf).threshmask \
	-variable fmri(threshmask) \
	-pattern "IMAGE" \
	-directory "$PWD" \
	-label "Pre-threshold masking" \
	-title "Select mask" \
	-width 30 \
	-filterhist VARS(history)

$w.bhelp bind $fmri(poststatsf).threshmask -msg "If you choose a mask for \"Pre-threshold masking\" then all stats
images will be masked by the chosen mask before thresholding. There
are two reasons why you might want to do this. The first is that you
might want to constrain your search for activation to a particular
area. The second is that in doing so, you are reducing the number of
voxels tested and therefore will make any
multiple-comparison-correction in the thresholding less stringent.

The mask image chosen does not have to be a binary mask - for example,
it can be a thresholded stats image from a previous analysis (in the
same space as the data to be analysed here); only voxels containing
zero in the \"mask\" image will get zeroed in this masking process."

#}}}
#{{{ thresholding

tixLabelFrame $w.thresh -label "Thresholding"
set fmri(lfthresh) [ $w.thresh subwidget frame ]

set tmpval $fmri(thresh)
tixOptionMenu $w.thresh.menu -variable fmri(thresh)
$w.thresh.menu add command 0 -label "None"
$w.thresh.menu add command 1 -label "Uncorrected"
$w.thresh.menu add command 2 -label "Voxel"
$w.thresh.menu add command 3 -label "Cluster"
set fmri(thresh) $tmpval

tixControl $w.uncorrected -label "Uncorrected voxel P threshold" \
	-variable fmri(prob_thresh) -step 0.005 -min 0 -max 1 -selectmode immediate
tixControl $w.voxel -label "Corrected voxel P threshold" \
	-variable fmri(prob_thresh) -step 0.005 -min 0 -max 1 -selectmode immediate
tixControl $w.prob_thresh -label " Cluster P threshold" \
	-variable fmri(prob_thresh) -step 0.005 -min 0 -max 1 -selectmode immediate
tixControl $w.z_thresh -label "Z threshold" \
	-variable fmri(z_thresh) -step 0.1 -min 0 -selectmode immediate

pack $w.thresh.menu -in $fmri(lfthresh) -side top -padx 5 -side left

$w.bhelp bind $w.thresh -msg "After carrying out the initial statistical test, the resulting Z
statistic image is then normally thresholded to show which voxels or
clusters of voxels are activated at a particular significance level.

If \"Cluster\" thresholding is selected, a Z statistic threshold is
used to define contiguous clusters.  Then each cluster's estimated
significance level (from GRF-theory) is compared with the cluster
probability threshold. Significant clusters are then used to mask the
original Z statistic image for later production of colour blobs. This
method of thresholding is an alternative to \"Voxel\"-based
correction, and is normally more sensitive to activation. You may well
want to increase the cluster creation \"Z threshold\" if you have high
levels of activation.

If \"Voxel\" thresholding is selected, GRF-theory-based maximum height
thresholding is carried out, with thresholding at the level set, using
one-tailed testing. This test is less overly-conservative than
Bonferroni correction.

You can also choose to simply threshold the uncorrected Z statistic
values, or apply no thresholding at all."

#}}}
#{{{ contrast masking

button $w.conmask -text "Contrast masking" -command "feat5:setup_conmask $w"

set fmri(conmask_help) "Setup the masking of contrasts by other contrasts; after thresholding
of all contrasts has taken place you can further \"threshold\" a given
Z statistic image by masking it with non-zeroed voxels from other
contrasts.

This means that of the voxels which passed thresholding in the
contrast (or F-test) of interest, only those which also survived
thresholding in the other contrasts (or F-tests) are kept.

As a further option, the generated masks can be derived from all
positive Z statistic voxels in the mask contrasts rather than all
voxels that survived thresholding."

$w.bhelp bind $w.conmask -msg $fmri(conmask_help) 

#}}}
#{{{ rendering

tixLabelFrame $w.render -label "Rendering"
set fmri(lfrendering) [ $w.render subwidget frame ]
set fmri(lfrenderingtop) [ frame $fmri(lfrendering).top ]

#{{{ Z display min and max

set tmpvalzdisplay $fmri(zdisplay)
tixOptionMenu $w.zmaxmenu -variable fmri(zdisplay)
$w.zmaxmenu add command 0 -label "Use actual Z min/max"
$w.zmaxmenu add command 1 -label "Use preset Z min/max"

tixControl $w.zmin -label "Min" -variable fmri(zmin) -step 1 -min 0 -selectmode immediate -options { entry.width 4 }
tixControl $w.zmax -label "Max" -variable fmri(zmax) -step 1 -min 0 -selectmode immediate -options { entry.width 4 }

$w.bhelp bind $w.zmaxmenu -msg "The Z statistic range selected for rendering is automatically
calculated by default, to run from red (minimum Z statistic after
thresholding) to yellow (maximum Z statistic). If more than one colour
rendered image is to be produced (i.e., when multiple constrasts are
created) then the overall range of Z values is automatically found
from all of the Z statistic images, for consistent Z statistic
colour-coding.

If multiple analyses are to be carried out, \"Use preset Z min/max\"
should be chosen, and the min/max values set by hand. Again, this
ensures consistency of Z statistic colour-coding - if several
experiments are to be reported side-by-side, colours will refer to the
same Z statistic values in each picture. When using this option, you
should choose a conservatively wide range for the min and max (e.g.,
min=1, max=15), to make sure that you do not carry out unintentional
thresholding via colour rendering."

#}}}
#{{{ render type

set tmpvalrendertype $fmri(rendertype)
tixOptionMenu $w.rendertype -variable fmri(rendertype)
$w.rendertype add command 0 -label "Solid blobs"
$w.rendertype add command 1 -label "Transparent blobs"

$w.bhelp bind $w.rendertype -msg "With \"Solid colours\" you don't see any sign of the background images
within the colour blobs; with \"Transparent colours\" you will see
through the colour blobs to the background intensity"

#}}}

pack $w.zmaxmenu $w.rendertype -in $fmri(lfrenderingtop) -side left -anchor n -padx 5 -pady 5
pack $fmri(lfrenderingtop) -in $fmri(lfrendering) -anchor w

#}}}

pack $fmri(poststatsf).threshmask -in $fmri(poststatsf) -side top -anchor w -padx 5 -pady 5
pack $w.thresh $w.render          -in $fmri(poststatsf) -side top -anchor w

$w.thresh.menu configure -command "feat5:updatepoststats $w"
$w.zmaxmenu configure -command "feat5:updatepoststats $w"
set fmri(zdisplay) $tmpvalzdisplay
set fmri(rendertype) $tmpvalrendertype

#{{{ background image for group stats

set tmpval $fmri(bgimage)
tixOptionMenu $w.bgimage -variable fmri(bgimage) -label " Background image"
$w.bgimage add command 1 -label "Mean highres"
$w.bgimage add command 2 -label "First highres"
$w.bgimage add command 3 -label "Mean functional"
$w.bgimage add command 4 -label "First functional"
$w.bgimage add command 5 -label "Standard space template"
set fmri(bgimage) $tmpval

$w.bhelp bind $w.bgimage -msg "With \"Higher-level analysis\" you can select what image will be used
as the background image for the activation colour overlays. The
default of \"Mean highres\" is probably the best for relating
activation to underlying structure. For a sharper underlying image,
(but one which is not so representative of the group of subjects), you
can instead choose to use the highres image from the first selected
subject.

You can alternatively choose to use the original lowres functional
data for the overlays, or the standard-space template image."

#}}}

#}}}
#{{{ Registration

set fmri(regf) [ $w.nb subwidget reg ]

set fmri(regreduce_dof) 0

#{{{ high res

frame $fmri(regf).initial_highres

checkbutton $fmri(regf).initial_highres.yn -variable fmri(reginitial_highres_yn) -command "feat5:updatereg_hr_init $w"

label $fmri(regf).initial_highres.label -text "Initial structural image"

tixLabelFrame $fmri(regf).initial_highres.lf -label "Initial structural image"
set fmri(initial_highresf) [ $fmri(regf).initial_highres.lf subwidget frame ]

set initial_highres_files(1) "~"

FSLFileEntry $fmri(initial_highresf).initial_highressingle \
	-variable initial_highres_files(1) \
	-pattern "IMAGE" \
	-directory $PWD \
	-label "" \
	-title "Select initial structural image" \
	-width 45 \
	-filterhist VARS(history) \
	-command "feat5:multiple_check $w 3 1 0"

button $fmri(initial_highresf).initial_highresmultiple -text "Select initial structural images" \
	-command "feat5:multiple_select $w 3 \"Select initial structural images\" "

frame $fmri(initial_highresf).opts

set tmpval $fmri(reginitial_highres_search)
tixOptionMenu $fmri(initial_highresf).opts.search -variable fmri(reginitial_highres_search) -label "  Linear "
$fmri(initial_highresf).opts.search add command 0 -label "No search"
$fmri(initial_highresf).opts.search add command 90 -label "Normal search"
$fmri(initial_highresf).opts.search add command 180 -label "Full search"
set fmri(reginitial_highres_search) $tmpval

set tmpval $fmri(reginitial_highres_dof)
tixOptionMenu $fmri(initial_highresf).opts.dof -variable fmri(reginitial_highres_dof)
$fmri(initial_highresf).opts.dof add command 3 -label "3 DOF (translation-only)"
$fmri(initial_highresf).opts.dof add command 6 -label "6 DOF"
$fmri(initial_highresf).opts.dof add command 7 -label "7 DOF"
$fmri(initial_highresf).opts.dof add command 9 -label "9 DOF"
$fmri(initial_highresf).opts.dof add command 12 -label "12 DOF"
set fmri(reginitial_highres_dof) $tmpval

label $fmri(initial_highresf).opts.nonlinear_label -text "     Nonlinear"
checkbutton $fmri(initial_highresf).opts.nonlinear_yn -variable fmri(reginitial_highres_nonlinear_yn)

#pack $fmri(initial_highresf).opts.search $fmri(initial_highresf).opts.dof $fmri(initial_highresf).opts.nonlinear_label \
#	$fmri(initial_highresf).opts.nonlinear_yn -in $fmri(initial_highresf).opts -side left
pack $fmri(initial_highresf).opts.search $fmri(initial_highresf).opts.dof -in $fmri(initial_highresf).opts -side left

pack $fmri(initial_highresf).initial_highressingle $fmri(initial_highresf).opts -in $fmri(initial_highresf) -anchor w -side top -pady 2 -padx 3
pack $fmri(regf).initial_highres.yn $fmri(regf).initial_highres.label -in $fmri(regf).initial_highres -side left

$w.bhelp bind $fmri(regf).initial_highres -msg "This is the initial high resolution structural image which the low
resolution functional data will be registered to, and this in turn
will be registered to the main highres image. It only makes sense to
have this initial highres image if a main highres image is also
specified and used in the registration. 

One example of an initial highres structural image might be a
medium-quality structural scan taken during a day's scanning, if a
higher-quality image has been previously taken for the subject. A
second example might be a full-brain image with the same MR sequence
as the functional data, useful if the actual functional data is only
partial-brain. It is strongly recommended that this image have
non-brain structures already removed, for example by using BET.

If the field-of-view of the functional data (in any direction) is less
than 120mm, then the registration of the functional data will by
default have a reduced degree-of-freedom, for registration stability.

If you are attempting to register partial field-of-view functional
data to a whole-brain image then \"3 DOF\" is recommended - in this
case only translations are allowed.

If the orientation of any image is different from any other image it
may be necessary to change the search to \"Full search\"."

#}}}
#{{{ high res2

frame $fmri(regf).highres

checkbutton $fmri(regf).highres.yn -variable fmri(reghighres_yn) -command "feat5:updatereg_hr $w"

label $fmri(regf).highres.label -text "Main structural image"

tixLabelFrame $fmri(regf).highres.lf -label "Main structural image"
set fmri(highresf) [ $fmri(regf).highres.lf subwidget frame ]

set highres_files(1) "~"

FSLFileEntry $fmri(highresf).highressingle \
	-variable highres_files(1) \
	-pattern "IMAGE" \
	-directory $PWD \
	-label "" \
	-title "Select main structural image" \
	-width 45 \
	-filterhist VARS(history) \
	-command "feat5:multiple_check $w 4 1 0"

button $fmri(highresf).highresmultiple -text "Select main structural images" \
	-command "feat5:multiple_select $w 4 \"Select main structural images\" "

frame $fmri(highresf).opts

set tmpval $fmri(reghighres_search)
tixOptionMenu $fmri(highresf).opts.search -variable fmri(reghighres_search) -label "  Linear "
$fmri(highresf).opts.search add command 0 -label "No search"
$fmri(highresf).opts.search add command 90 -label "Normal search"
$fmri(highresf).opts.search add command 180 -label "Full search"
set fmri(reghighres_search) $tmpval

set tmpval $fmri(reghighres_dof)
tixOptionMenu $fmri(highresf).opts.dof -variable fmri(reghighres_dof)
$fmri(highresf).opts.dof add command 3 -label "3 DOF (translation-only)"
$fmri(highresf).opts.dof add command 6 -label "6 DOF"
$fmri(highresf).opts.dof add command 7 -label "7 DOF"
$fmri(highresf).opts.dof add command 9 -label "9 DOF"
$fmri(highresf).opts.dof add command 12 -label "12 DOF"
set fmri(reghighres_dof) $tmpval

label $fmri(highresf).opts.nonlinear_label -text "     Nonlinear"
checkbutton $fmri(highresf).opts.nonlinear_yn -variable fmri(reghighres_nonlinear_yn)

#pack $fmri(highresf).opts.search $fmri(highresf).opts.dof $fmri(highresf).opts.nonlinear_label \
#	$fmri(highresf).opts.nonlinear_yn -in $fmri(highresf).opts -side left
pack $fmri(highresf).opts.search $fmri(highresf).opts.dof -in $fmri(highresf).opts -side left

pack $fmri(highresf).highressingle $fmri(highresf).opts -in $fmri(highresf) -anchor w -side top -pady 2 -padx 3
pack $fmri(regf).highres.yn $fmri(regf).highres.label -in $fmri(regf).highres -side left

$w.bhelp bind $fmri(regf).highres -msg "This is the main high resolution structural image which the low resolution
functional data will be registered to (optionally via the \"initial
structural image\"), and this in turn will be registered to the
standard brain. It is strongly recommended that this image have
non-brain structures already removed, for example by using BET.

If the field-of-view of the functional data (in any direction) is less
than 120mm, then the registration of the functional data will by
default have a reduced degree-of-freedom, for registration stability.

If you are attempting to register partial field-of-view functional
data to a whole-brain image then \"3 DOF\" is recommended - in this
case only translations are allowed.

If the orientation of any image is different from any other image it
may be necessary to change the search to \"Full search\"."

#}}}
#{{{ standard

frame $fmri(regf).standard

checkbutton $fmri(regf).standard.yn -variable fmri(regstandard_yn) -command "feat5:updatereg $w"

label $fmri(regf).standard.label -text "Standard space"

tixLabelFrame $fmri(regf).standard.lf -label "Standard space"
set fmri(standardf) [ $fmri(regf).standard.lf subwidget frame ]

FSLFileEntry  $fmri(standardf).standardentry \
	-variable fmri(regstandard) \
	-pattern "IMAGE" \
	-directory "${FSLDIR}/etc/standard/" \
	-label "" \
	-title "Select the standard image" \
	-width 45 \
	-filterhist VARS(history)

frame $fmri(standardf).opts

set tmpval $fmri(regstandard_search)
tixOptionMenu $fmri(standardf).opts.search -variable fmri(regstandard_search) -label "  Linear "
$fmri(standardf).opts.search add command 0 -label "No search"
$fmri(standardf).opts.search add command 90 -label "Normal search"
$fmri(standardf).opts.search add command 180 -label "Full search"
set fmri(regstandard_search) $tmpval

set tmpval $fmri(regstandard_dof)
tixOptionMenu $fmri(standardf).opts.dof -variable fmri(regstandard_dof)
$fmri(standardf).opts.dof add command 3 -label "3 DOF (translation-only)"
$fmri(standardf).opts.dof add command 6 -label "6 DOF"
$fmri(standardf).opts.dof add command 7 -label "7 DOF"
$fmri(standardf).opts.dof add command 9 -label "9 DOF"
$fmri(standardf).opts.dof add command 12 -label "12 DOF"
set fmri(regstandard_dof) $tmpval

label $fmri(standardf).opts.nonlinear_label -text "     Nonlinear"
checkbutton $fmri(standardf).opts.nonlinear_yn -variable fmri(regstandard_nonlinear_yn)

#pack $fmri(standardf).opts.search $fmri(standardf).opts.dof $fmri(standardf).opts.nonlinear_label \
#	$fmri(standardf).opts.nonlinear_yn -in $fmri(standardf).opts -side left
pack $fmri(standardf).opts.search $fmri(standardf).opts.dof -in $fmri(standardf).opts -side left

pack $fmri(standardf).standardentry $fmri(standardf).opts -in $fmri(standardf) -anchor w -side top -pady 2 -padx 3
pack $fmri(regf).standard.yn $fmri(regf).standard.label -in $fmri(regf).standard -side left

$w.bhelp bind $fmri(regf).standard -msg "This is the standard (reference) image;
it should be an image already in MNI152/Talairach/etc. space,
ideally with the non-brain structures already removed.

If the field-of-view of the functional data (in any direction) is less
than 120mm, then the registration of the functional data will by
default have a reduced degree-of-freedom, for registration stability.

If you are attempting to register partial field-of-view functional
data to a whole-brain image then \"3 DOF\" is recommended - in this
case only translations are allowed.

If the orientation of any image is different from any other image it
may be necessary to change the search to \"Full search\"."

#}}}

pack $fmri(regf).initial_highres $fmri(regf).highres $fmri(regf).standard -in $fmri(regf) -side top -anchor w -pady 0

#}}}

set fmri(level) 1
set fmri(analysis) 7
$w.mode.level    configure -command "feat5:updatelevel $w"
$w.mode.analysis configure -command "feat5:updateanalysis $w"

set tmpval $fmri(paradigm_hp)
feat5:updatelevel $w 0
set fmri(paradigm_hp) $tmpval

$w.nb raise data

#}}}
    #{{{ button Frame

frame $w.btns
    
button $w.btns.apply -command "feat5:apply $w" -text "Go"

FSLFileDirSelectDialog $w.savedialog \
	-title  "Save FEAT setup" \
	-command "feat5:write $w 1 1 0"
$w.savedialog subwidget fsbox config -pattern "*.fsf" -filterhist VARS(history)
button $w.btns.save -command "$w.savedialog popup" -text "Save"

FSLFileDirSelectDialog $w.loaddialog \
	-title  "Load FEAT setup" \
	-command "feat5:load $w 1"
$w.loaddialog subwidget fsbox config -pattern "*.fsf" -filterhist VARS(history)
button $w.btns.load -command "$w.loaddialog popup" -text "Load"

button $w.btns.cancel -command "destroy $w" -text "Exit"

button $w.btns.help -command "FmribWebHelp file: ${FSLDIR}/doc/feat5/index.html" -text "Help"

#{{{ Utils

menubutton $w.btns.utils -text "Utils" -menu $w.btns.utils.menu -relief raised -bd 2

menu $w.btns.utils.menu

$w.btns.utils.menu add command -label "Make_flobs - create optimal basis set (of HRF convolution kernels)" -command { exec sh -c "${FSLDIR}/bin/Make_flobs$gui_ext" & }

$w.btns.utils.menu add command -label "Featquery - get FEAT stats from ROI mask or co-ordinates" -command { exec sh -c "${FSLDIR}/bin/Featquery$gui_ext" & }

$w.btns.utils.menu add command -label "High-res FEAT stats colour rendering" -command { exec sh -c "${FSLDIR}/bin/Renderhighres$gui_ext" & }

#}}}

pack $w.btns.apply $w.btns.save $w.btns.load $w.btns.cancel $w.btns.help $w.btns.utils -in $w.btns -side left -expand yes

#}}}

    pack $w.mode $w.nb -in $w -side top -anchor n -padx 10 -pady 10
    pack $w.btns -in $w -side bottom -fill x -padx 10 -pady 10

    #{{{ load fsf file

if { $argc > 0 } {

    set inputname [ lindex $argv 0 ]

    if { [ string first / $inputname ] != 0 && [ string first ~ $inputname ] != 0 } {
	set inputname ${PWD}/$inputname
    }

    if { [ string compare [ file extension $inputname ] .fsf ] == 0 } {

	if { [ file readable $inputname ] } {
	    puts "Loading FEAT setup file $inputname"
	    feat5:load $w 1 $inputname
	} else {
	    MxPause "setup file $inputname doesn't exist!"
	}
    }

}

#}}}
    #{{{ updates needed after the loading of settings

if { $fmri(perfsub_yn) } {
    pack $w.temp.tcmenu -in $w.temp -after $w.temp.ps_yn -side top -side left -padx 5
}

#}}}
}

#}}}

##### FEAT analysis procedures #####
#{{{ feat5:getconname

proc feat5:getconname { featdir contrastnumber } {
    source ${featdir}/design.fsf
    if { [ info exists fmri(conname_real.$contrastnumber) ] } {
	return $fmri(conname_real.$contrastnumber)
    } else {
	return ""
    }
}

#}}}
#{{{ feat5:connectivity

proc feat5:connectivity { image } {

    global FSLDIR

    set CONNECTIVITY "26"

#    maybe actually USE this test? but needs careful thought.......
#
#    set CONNECTIVITY "6"
#
#    if { [ expr abs([ exec sh -c "$FSLDIR/bin/avwval $image pixdim3" ]) ] < 3.1 } {
#	set CONNECTIVITY "26"
#    }

    return $CONNECTIVITY
}

#}}}
#{{{ feat5:poststats

proc feat5:poststats { rerunning stdspace } {
    global FSLDIR FD report logout comout fmri
    set ps ""
    set rs ""

    if { $fmri(thresh) == 0 } {
	return 0
    }

    set maskcomments "<p>"
    
    #{{{ setup raw stats list

cd ${FD}/stats
set rawstatslist [ remove_ext [ lsort -dictionary [ imglob -oneperimage zstat*.* ] ] [ lsort -dictionary [ imglob -oneperimage zfstat*.* ] ] ]
cd ${FD}

#}}}
    #{{{ setup standard-space/non-standard-space thingies

set STDOPT   ""
set STDEXT   ""
set SLICER   "-A"
set VOXorMM  ""

if { $stdspace != 0 } {
    set STDOPT   "-std"
    set STDEXT   "_std"
    set SLICER   "-S 2"
    set VOXorMM  "--mm"
}

#}}}
    #{{{ brain-mask Z stats and find smoothness

foreach rawstats $rawstatslist {

    if { ! $rerunning } {
	fsl:exec "$FSLDIR/bin/avwmaths stats/$rawstats -mas mask thresh_$rawstats"

	if { $fmri(threshmask) != "" && [ imtest $fmri(threshmask) ] } {
	    fsl:exec "$FSLDIR/bin/avwmaths thresh_$rawstats -mas $fmri(threshmask) thresh_$rawstats"
	    set maskcomments "<p>Z-stat images were masked with $fmri(threshmask) before thresholding.<br>"
	}
    }

    if { [ file exists stats/smoothness ] } {
	set fmri(DLH$rawstats)     [ exec sh -c " grep DLH    stats/smoothness | awk '{ print \$2 }'" ]
	set fmri(RESELS$rawstats)  [ exec sh -c " grep RESELS stats/smoothness | awk '{ print \$2 }'" ]
    } else {
	fsl:exec "$FSLDIR/bin/smoothest -m mask -z stats/$rawstats > stats/${rawstats}.smoothness"
	set fmri(DLH$rawstats)     [ exec sh -c " grep DLH    stats/${rawstats}.smoothness | awk '{ print \$2 }'" ]
	set fmri(RESELS$rawstats)  [ exec sh -c " grep RESELS stats/${rawstats}.smoothness | awk '{ print \$2 }'" ]
    }

    if { ! $rerunning } {
	set fmri(VOLUME$rawstats) [ exec sh -c " ${FSLDIR}/bin/avwstats thresh_$rawstats -V | awk '{ print \$1 }'" ]
	fsl:exec "echo $fmri(VOLUME$rawstats) > thresh_${rawstats}.vol"
    } else {
	if { ! [ info exists fmri(VOLUME$rawstats) ] } {
	    if { [ file exists thresh_${rawstats}.vol ] } {
		set fmri(VOLUME$rawstats) [ exec sh -c "cat thresh_${rawstats}.vol" ]
	    } else {
		set fmri(VOLUME$rawstats) [ exec sh -c "${FSLDIR}/bin/avwstats mask -V | awk '{ print \$1 }'" ]
	    }
	}
    }

    fsl:echo $logout "$rawstats: DLH=$fmri(DLH$rawstats) VOLUME=$fmri(VOLUME$rawstats) RESELS=$fmri(RESELS$rawstats)"

}

#}}}
    #{{{ thresholding and contrast masking

if { ! $rerunning } {

    set firsttime 1
    foreach rawstats $rawstatslist {
	set i [ string trimleft $rawstats "abcdefghijklmnopqrstuvwxyz_" ]

	if { $firsttime == 1 } {
	    set ps "$ps Z (Gaussianised T/F) statistic images were thresholded"
	}

	if { $fmri(thresh) < 3 } {
	    #{{{ voxel-based thresholding

if { $firsttime == 1 } {
    if { $fmri(thresh) == 1 } {
	set ps "$ps at P=$fmri(prob_thresh) (uncorrected)."
	set zthresh [ fsl:exec "${FSLDIR}/bin/ptoz $fmri(prob_thresh)" ]
    } else {
	set ps "$ps using GRF-theory-based maximum height thresholding with a (corrected) significance threshold of P=$fmri(prob_thresh) \[Worsley 1992]."
	set rs "$rs\[Worsley 1992\] K.J. Worsley, A.C. Evans, S. Marrett and P. Neelin. A three-dimensional statistical analysis
                for CBF activation studies in human brain. Journal of Cerebral Blood Flow and Metabolism 12(900-918) 1992.<br> 
                "
    }
}

if { $fmri(thresh) == 2 } {
    set zthresh [ fsl:exec "${FSLDIR}/bin/ptoz $fmri(prob_thresh) -g [ expr int ( $fmri(VOLUME$rawstats) / $fmri(RESELS$rawstats) ) ]" ]
}

fsl:exec "$FSLDIR/bin/avwmaths thresh_$rawstats -thr $zthresh thresh_$rawstats"

#}}}
	} else {
	    #{{{ cluster thresholding

if { $firsttime == 1 } {
    set ps "$ps using clusters determined by Z>$fmri(z_thresh) and a (corrected) cluster significance threshold of P=$fmri(prob_thresh) \[Worsley 1992]."
    set rs "$rs\[Worsley 1992\] K.J. Worsley, A.C. Evans, S. Marrett and P. Neelin. A three-dimensional statistical analysis
    for CBF activation studies in human brain. Journal of Cerebral Blood Flow and Metabolism 12(900-918) 1992.<br> 
    "
}

set COPE ""
if { [ string first "zfstat" $rawstats ] < 0 && [ imtest stats/cope${i} ] } {
    set COPE "-c stats/cope$i"
}

fsl:exec "$FSLDIR/bin/cluster -i thresh_$rawstats $COPE -t $fmri(z_thresh) -p $fmri(prob_thresh) -d $fmri(DLH$rawstats) --volume=$fmri(VOLUME$rawstats) --othresh=thresh_$rawstats -o cluster_mask_$rawstats --connectivity=[ feat5:connectivity thresh_$rawstats ] $VOXorMM --olmax=lmax_${rawstats}${STDEXT}.txt > cluster_${rawstats}${STDEXT}.txt"

fsl:exec "$FSLDIR/bin/cluster2html . cluster_$rawstats $STDOPT"

#}}}
	}

	set firsttime 0
    }

    #{{{ contrast masking

if { $fmri(conmask1_1) } {

    fsl:exec "mkdir conmask"

    foreach rawstats $rawstatslist {

	set i [ string trimleft $rawstats "abcdefghijklmnopqrstuvwxyz_" ]

	# check for being F-test
	set I $i
	if { [ string first "zstat" $rawstats ] == -1 } {
	    incr I $fmri(ncon_real)
	}

	set theinput thresh_$rawstats

	for { set C 1 } { $C <= [ expr $fmri(ncon_real) + $fmri(nftests_real) ] } { incr C } {
	    if { $C != $I } {

		set F ""
		set c $C
		if { $C > $fmri(ncon_real) } {
		    set F f
		    set c [ expr $C - $fmri(ncon_real) ]
		}

		if { $fmri(conmask${I}_$C) } {
		    
		    set themask thresh_z${F}stat$c
		    if { $fmri(conmask_zerothresh_yn) } {
			set themask stats/z${F}stat$c
		    }

		    fsl:echo $logout "Masking $theinput with $themask"

		    set maskcomments "$maskcomments
After all thresholding, $rawstats was masked with $themask.<br>"

		    if { [ imtest $themask ] } {
			fsl:exec "${FSLDIR}/bin/avwmaths $theinput -mas $themask conmask/thresh_$rawstats"
		    } else {
			fsl:exec "${FSLDIR}/bin/avwmaths $theinput -mul 0 conmask/thresh_$rawstats"
		    }

		    set theinput conmask/thresh_$rawstats
		}
		
	    }
	}
    }

    fsl:exec "/bin/mv -f conmask/* . ; rmdir conmask"

    #{{{ redo clustering

if { $fmri(thresh) == 3 } {
    foreach rawstats $rawstatslist {
	set i [ string trimleft $rawstats "abcdefghijklmnopqrstuvwxyz_" ]
	set COPE ""
	if { [ string first "zfstat" $rawstats ] < 0 && [ imtest stats/cope${i} ] } {
	    set COPE "-c stats/cope$i"
	}

	# we're not going to re-test cluster size so pthresh is set to 1000

	fsl:exec "$FSLDIR/bin/cluster -i thresh_$rawstats $COPE -t $fmri(z_thresh) -d $fmri(DLH$rawstats) --volume=$fmri(VOLUME$rawstats) --othresh=thresh_$rawstats -o cluster_mask_$rawstats --connectivity=[ feat5:connectivity thresh_$rawstats ] $VOXorMM --olmax=lmax_${rawstats}${STDEXT}.txt > cluster_${rawstats}${STDEXT}.txt"

	fsl:exec "$FSLDIR/bin/cluster2html . cluster_$rawstats $STDOPT"
    }
}

#}}}
}

#}}}
}

#}}}
    #{{{ re-run cluster for StdSpace

if { $rerunning && [ file exists reg/example_func2standard.mat ] && $fmri(thresh) == 3 } {

    set z_thresh    $fmri(z_thresh)
    set prob_thresh $fmri(prob_thresh)

    if { $fmri(analysis) == 0 } {
	set z_thresh    [ exec sh -c " grep 'set fmri(z_thresh)'    design.fsf | awk '{ print \$3 }'" ]
	set prob_thresh [ exec sh -c " grep 'set fmri(prob_thresh)' design.fsf | awk '{ print \$3 }'" ]
    }

    foreach rawstats $rawstatslist {

	set i [ string trimleft $rawstats "abcdefghijklmnopqrstuvwxyz_" ]

	set COPE ""
	if { [ string first "zfstat" $rawstats ] < 0 && [ imtest stats/cope${i} ] } {
	    set COPE "-c stats/cope$i"
	}

	fsl:exec "$FSLDIR/bin/cluster -i thresh_$rawstats ${COPE} -t $z_thresh -d $fmri(DLH$rawstats) --volume=$fmri(VOLUME$rawstats) -x reg/example_func2standard.mat --stdvol=reg/standard --mm --connectivity=[ feat5:connectivity thresh_$rawstats ] --olmax=lmax_${rawstats}_std.txt > cluster_${rawstats}_std.txt"
	fsl:exec "$FSLDIR/bin/cluster2html . cluster_${rawstats} -std"
    }
}

#}}}
    #{{{ rendering

if { ! $rerunning } {

    #{{{ Find group Z min and max

if { $fmri(zdisplay) == 0 } {

    set fmri(zmin) 100000
    set fmri(zmax) -100000

    foreach rawstats $rawstatslist {
    
	set zminmax [ fsl:exec "${FSLDIR}/bin/avwstats thresh_$rawstats -l 0.0001 -R" ]
	set zmin [ lindex $zminmax 0 ]
	set zmax [ lindex $zminmax 1 ]

	# test for non-empty thresh_zstats image
	if { $zmax > 0.0001 } {
	    if { $fmri(zmin) > $zmin } {
		set fmri(zmin) $zmin
	    }
	    if { $fmri(zmax) < $zmax } {
		set fmri(zmax) $zmax
	    }
	}
    }

    if { $fmri(zmin) > 99999 } {
	set fmri(zmin) 2.3
	set fmri(zmax) 8
    }
}

fsl:echo $logout "Rendering using zmin=$fmri(zmin) zmax=$fmri(zmax)"

#}}}

    set underlying example_func
    #    if { ! [ imtest $underlying ] } {
    #	set underlying mean_highres
    #    }

    set firsttime 1
    foreach rawstats $rawstatslist {
	
	set i [ string trimleft $rawstats "abcdefghijklmnopqrstuvwxyz_" ]
	
	if { [ string first "zstat" $rawstats ] < 0 || $fmri(conpic_real.$i) == 1 } {
	    #{{{ Rendering

set conname "$rawstats"
if { [ string first "zfstat" $rawstats ] < 0 } {
    set conname "$conname &nbsp;&nbsp;-&nbsp;&nbsp; C${i}"
    if { $fmri(conname_real.$i) != "" } {
	set conname "$conname ($fmri(conname_real.$i))"
    }
} else {
    set conname "$conname &nbsp;&nbsp;-&nbsp;&nbsp; F${i}"

    set start 1
    for { set c 1 } { $c <= $fmri(ncon_real) } { incr c 1 } {
	if { $fmri(ftest_real${i}.${c}) == 1 } {
	    if { $start == 1 } {
		set conname "$conname ("
		set start 0
	    } else {
		set conname "$conname & "
	    }
	    set conname "${conname}C$c"
	}
    }
    set conname "$conname)"
}

fsl:exec "$FSLDIR/bin/overlay $fmri(rendertype) 0 $underlying -a thresh_$rawstats $fmri(zmin) $fmri(zmax) rendered_thresh_$rawstats"
fsl:exec "${FSLDIR}/bin/slicer rendered_thresh_$rawstats $SLICER 750 rendered_thresh_${rawstats}.png"

if { $firsttime == 1 } {
    fsl:exec "/bin/cp ${FSLDIR}/etc/luts/ramp.gif .ramp.gif"
    feat5:report_insert poststatsps $ps
    feat5:report_insert poststatsrs $rs
    feat5:report_insert_start poststatspics
    puts $report "
<hr><b>Thresholded activation images</b>
&nbsp; &nbsp; &nbsp; &nbsp; 
[ expr int($fmri(zmin)*10)/10.0 ]
<IMG BORDER=0 SRC=\".ramp.gif\">
[ expr int($fmri(zmax)*10)/10.0 ]
$maskcomments
"
    set firsttime 0
}

if { $fmri(thresh) == 3 } {
    puts $report "<p>$conname<br>
    <a href=\"cluster_${rawstats}${STDEXT}.html\"><IMG BORDER=0 SRC=\"rendered_thresh_${rawstats}.png\"></a>
    "
} else {
    puts $report "<p>$conname<br>
    <IMG BORDER=0 SRC=\"rendered_thresh_${rawstats}.png\">
    "
}

#}}}
	}
    }

    feat5:report_insert_stop poststatspics
}

#}}}

}

#}}}
#{{{ feat5:flirt

proc feat5:flirt { in ref dof search interp existing_mats report init in_weighting } {

    global FSLDIR FSLCONVERT logout comout

    set out ${in}2$ref

    if { $existing_mats } {

	fsl:exec "${FSLDIR}/bin/flirt -ref $ref -in $in -out $out -applyxfm -init ${out}.mat -interp $interp"

    } else {

	if { $dof == 3 } {
	    set dof "6 -schedule ${FSLDIR}/etc/flirtsch/sch3Dtrans_3dof"
	}

	fsl:exec "${FSLDIR}/bin/flirt -ref $ref -in $in -out $out -omat ${out}.mat -cost corratio -dof $dof -searchrx -$search $search -searchry -$search $search -searchrz -$search $search -interp $interp $init $in_weighting" -r

    }

    fsl:exec "${FSLDIR}/bin/convert_xfm -inverse -omat ${ref}2${in}.mat ${out}.mat"

    fsl:exec "${FSLDIR}/bin/slicer $out $ref -s 2 -x 0.35 sla -x 0.45 slb -x 0.55 slc -x 0.65 sld -y 0.35 sle -y 0.45 slf -y 0.55 slg -y 0.65 slh -z 0.35 sli -z 0.45 slj -z 0.55 slk -z 0.65 sll ; $FSLCONVERT -colors 100 +append sla slb slc sld sle slf slg slh sli slj slk sll ${out}1.gif ; ${FSLDIR}/bin/slicer $ref $out -s 2 -x 0.35 sla -x 0.45 slb -x 0.55 slc -x 0.65 sld -y 0.35 sle -y 0.45 slf -y 0.55 slg -y 0.65 slh -z 0.35 sli -z 0.45 slj -z 0.55 slk -z 0.65 sll ; $FSLCONVERT -colors 100 +append sla slb slc sld sle slf slg slh sli slj slk sll ${out}2.gif ; $FSLCONVERT -colors 100 -append ${out}1.gif ${out}2.gif ${out}.gif; /bin/rm -f sla slb slc sld sle slf slg slh sli slj slk sll"

    #imrm $out

    if { $report != "" } {
	puts $report "<hr><p>Registration of $in to $ref"
	if { $in == "example_func_orig_distorted" } {
	    puts $report " (for comparison with unwarped example_func vs highres, shown above)"
	}
	puts $report "<br><br><a href=\"${out}.gif\"><IMG BORDER=0 SRC=\"${out}.gif\" WIDTH=2000></a><br>"
    }
}

#}}}
#{{{ feat5:find_std

proc feat5:find_std { featdir image } {
    global FSLDIR

    if {  [ file exists ${featdir}/design.lev ] } {
	if { [ imtest ${featdir}/$image ] } {
	    return ${featdir}/$image
	} else {
	    return 0
	}
    } else {
	if { [ imtest ${featdir}/reg_standard/$image ] } {
	    return ${featdir}/reg_standard/$image
	} elseif { $image == "standard" && [ imtest ${featdir}/reg/standard ] } {
	    return ${featdir}/reg/standard
	} else {
	    return 0
	}
    }
}

#}}}
#{{{ feat5:report_insert

proc feat5:report_insert { sectionlabel insertstring } {
    global report

    feat5:report_insert_start $sectionlabel
    puts $report "$insertstring"
    feat5:report_insert_stop $sectionlabel
}

proc feat5:report_insert_start { sectionlabel } {
    global report

    if { [ info exists report ] } {
	catch { close $report } errmsg
    }

    catch { exec sh -c "mv report.html tmpreport.html" } errmsg
    set iptr [ open tmpreport.html r ]
    set report [ open report.html w ]
    fconfigure $report -buffering none
    set foundit 0

    while { ! $foundit && [ gets $iptr line ] >= 0 } {
	if { [ regexp "<!--${sectionlabel}start-->" $line ] } {
	    set foundit 1
	}
	puts $report $line
    }

    close $iptr
}

proc feat5:report_insert_stop { sectionlabel } {
    global report

    set iptr [ open tmpreport.html r ]
    set foundit 0

    while { [ gets $iptr line ] >= 0 } {
	if { [ regexp "<!--${sectionlabel}stop-->" $line ] } {
	    set foundit 1
	}
	if { $foundit } {
	    puts $report $line
	}
    }
    
    close $iptr
    exec sh -c "rm -f tmpreport.html"
}

#}}}
#{{{ stringstrip

proc stringstrip { in ext } {

    set lengthin [ string length $in ]
    set lengthext [ string length $ext ]

    if { $lengthin > $lengthext } {

	if { [ string compare [ string range $in [ expr $lengthin - $lengthext ] $lengthin ] $ext ] == 0 } {
	    return [ string range $in 0 [ expr $lengthin - $lengthext - 1 ] ]
	} else { 
	    return $in
	}

    } else {
	return $in
    }
}

proc feat5:strip { in } {
    set in [ stringstrip $in .feat ]
    set in [ stringstrip $in .gfeat ]
    set in [ remove_ext $in ]
    return $in
}

#}}}
#{{{ feat5:proc

proc feat5:proc { fsfroot } {

    #{{{ setup and load fsf

global FSLDIR FSLCONVERT FSLGNUPLOT FSLSLASH PWD HOME HOSTNAME OSFLAVOUR logout fmri feat_files unwarp_files unwarp_files_mag initial_highres_files highres_files FD report ps rs logout comout gui_ext

if { [ string range $fsfroot 0 0 ] != "/" } {
    set fsfroot ${PWD}/$fsfroot
}

feat5:setupdefaults
feat5:load -1 1 ${fsfroot}.fsf

#}}}
    #{{{ prepare data for higher-level analyses

if { $fmri(level)>1 && $fmri(analysis)!=4 } {

    #{{{ if inputtype=copes, fix feat_files

if { $fmri(inputtype) == 2 } {
    for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
	set copes($session) $feat_files($session)
	set feat_files($session) [ file dirname [ file dirname $feat_files($session) ] ]
    }
}	

#}}}
    #{{{ setup gfeat directory

set thedate [ exec date ]

if { $fmri(outputdir) != "" } {
    set featname [ feat5:strip $fmri(outputdir) ].gfeat
} else {
    set featname [ feat5:strip [ remove_ext $feat_files(1) ] ].gfeat
}
set FD [ new_filename $featname ]
if { ! [ file writable [ file dirname $FD ] ] } {
    set FD [ new_filename ${HOME}/[ file tail $featname ] ]
}
fsl:exec "mkdir $FD" -n

cd $FD
set FD [ pwd ]
set gFD $FD
set logout ${FD}/report.log
set comout ${FD}/report.com

fsl:echo $logout "Started higher-level FEAT at $thedate on $HOSTNAME
FEAT output is ${FD}"

if { $fmri(featwatcher_yn) } {
    catch { exec sh -c "${FSLDIR}/bin/Featwatcher${gui_ext} $FD" >& /dev/null & } errmsg
} else {
    if { $OSFLAVOUR != "cygwin" } {
	set tailpid [ exec tail -f report.log & ]
    }
}

fsl:exec "/bin/cp ${fsfroot}.fsf design.fsf"
fsl:exec "${FSLDIR}/bin/feat_model design"

set fsfroot ${FD}/design

#}}}
    #{{{ write start of report web page

set greport [ open report.html "w" ]
fconfigure $greport -buffering none

puts $greport "<HTML>

<TITLE>Higher-Level FEAT</TITLE>

<BODY BACKGROUND=\"file:${FSLSLASH}${FSLDIR}/doc/images/fsl-bg.jpg\">

<hr><CENTER>
<H1>Higher-Level FEAT Index</H1>
${FD}/report.html<br>
$thedate
</CENTER>
"

puts $greport "<hr><H2>Input (lower-level) reports</H2>"

for { set i 1 } { $i <= $fmri(multiple) } { incr i 1 } {
    puts $greport "${i} <A HREF=\"$feat_files($i)/report.html\">$feat_files($i)/report.html</A><br>
    "
}

#}}}
    #{{{ run featregapply on each input feat directory

for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
    fsl:exec "${FSLDIR}/bin/featregapply $feat_files($session)"
}

#}}}
    #{{{ setup background image

if { [ file exists $feat_files(1)/design.lev ] } {
    set fmri(bgimage) 4
}

set bg_list ""

switch $fmri(bgimage) {
    1 {
	for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
	    set find_std [ feat5:find_std $feat_files($session) reg/highres ]
	    if { $find_std != 0 } {
		set bg_list "$bg_list $find_std"
	    }
	}
    }
    2 {
	for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
	    set find_std [ feat5:find_std $feat_files($session) reg/highres ]
	    if { $find_std != 0 } {
		set bg_list "$bg_list $find_std"
		set session $fmri(multiple)
	    }
	}
    }
    3 {
	for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
	    set find_std [ feat5:find_std $feat_files($session) example_func ]
	    if { $find_std != 0 } {
		set bg_list "$bg_list $find_std"
	    }
	}
    }
    4 {
	for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
	    set find_std [ feat5:find_std $feat_files($session) example_func ]
	    if { $find_std != 0 } {
		set bg_list "$bg_list $find_std"
		set session $fmri(multiple)
	    }
	}
    }
}

if { [ string length $bg_list ] == 0 } {
    fsl:exec "${FSLDIR}/bin/avwmaths [ feat5:find_std $feat_files(1) standard ] bg_image"
} else {
    fsl:exec "${FSLDIR}/bin/avwmerge -t bg_image $bg_list"
    if { [ llength $bg_list ] > 1 } {
	fsl:exec "${FSLDIR}/bin/ip bg_image bg_image 10 -i 1000"
	fsl:exec "${FSLDIR}/bin/avwmaths bg_image -Tmean bg_image"
    }
}

#}}}
    #{{{ setup mask image and inputreg report

set mask_list ""

for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
    set mask_list "$mask_list [ feat5:find_std $feat_files($session) mask ]"
}

fsl:exec "${FSLDIR}/bin/avwmerge -t mask $mask_list"

# make input reg report if at second level
if { ! [ file exists $feat_files(1)/design.lev ] } {

    #{{{ create inputreg images

fsl:exec "mkdir inputreg"

cd inputreg

fsl:exec "${FSLDIR}/bin/avwmaths_16SI ../mask -mul $fmri(multiple) -Tmean masksum"
fsl:exec "${FSLDIR}/bin/avwmaths_16SI masksum -thr $fmri(multiple) -add masksum masksum"
fsl:exec "$FSLDIR/bin/overlay 0 0 -c ../bg_image -a masksum 0.9 [ expr 2 * $fmri(multiple) ] masksum_overlay"
fsl:exec "${FSLDIR}/bin/slicer masksum_overlay -S 2 750 masksum_overlay.png"
imrm masksum_overlay

fsl:exec "${FSLDIR}/bin/avwmaths_16SI masksum -mul 0 maskunique"
for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
    fsl:exec "${FSLDIR}/bin/avwmaths_16SI [ feat5:find_std $feat_files($session) mask ] -mul -1 -add 1 -mul $session -add maskunique maskunique"
}
fsl:exec "${FSLDIR}/bin/avwmaths_16SI masksum -thr [ expr $fmri(multiple) - 1 ] -uthr [ expr $fmri(multiple) - 1 ] -bin -mul maskunique maskunique"
fsl:exec "$FSLDIR/bin/overlay 0 0 ../bg_image -a maskunique 0.9 $fmri(multiple) maskunique_overlay"
fsl:exec "${FSLDIR}/bin/slicer maskunique_overlay -S 2 750 maskunique_overlay.png"
imrm maskunique_overlay

#}}}
    #{{{ create inputreg webpage

fsl:exec "/bin/cp ${FSLDIR}/etc/luts/ramp.gif .ramp.gif"

set regreport [ open index.html "w" ]

puts $regreport "<HTML>

<TITLE>Input Registration Summary</TITLE>

<BODY BACKGROUND=\"file:${FSLSLASH}${FSLDIR}/doc/images/fsl-bg.jpg\">

<hr><CENTER>
<H1>Input Registration Summary</H1>
for higher-level FEAT analysis &nbsp; &nbsp; <a href=\"../report.html\">${FD}</a><br>
[ exec date ]
</CENTER>

<hr><p><font size=+1><b>Summaries of functional-to-standard registrations for all inputs</b></font><br><br>
"

for { set i 1 } { $i <= $fmri(multiple) } { incr i 1 } {
    puts $regreport "${i} <A HREF=\"$feat_files($i)/reg/index.html\">$feat_files($i)<br>
<IMG BORDER=0 SRC=\"$feat_files($i)/reg/example_func2standard.gif\"></A><br><br>
"
}

puts $regreport "<hr><p><font size=+1><b>Sum of all input masks after transformation to standard space
&nbsp; &nbsp; &nbsp; &nbsp; 1 <IMG BORDER=0 SRC=\".ramp.gif\"> $fmri(multiple)</b></font><br><br>
<IMG BORDER=0 SRC=\"masksum_overlay.png\">

<hr><p><font size=+1><b>Unique missing-mask voxels
&nbsp; &nbsp; &nbsp; &nbsp; 1 <IMG BORDER=0 SRC=\".ramp.gif\"> $fmri(multiple)</b></font><br><br>
This shows voxels where only one mask is missing, to enable easy identification of single gross registration problems. For detail, view image ${FD}/inputreg/maskunique<br><br>
<IMG BORDER=0 SRC=\"maskunique_overlay.png\">

<HR><FONT SIZE=1>This page produced automatically by FEAT - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL</A>.</FONT>
</BODY></HTML>
"

close $regreport

#}}}

    cd $FD

    puts $greport "<hr><H2><a href=\"inputreg/index.html\">Summary of lower-level registrations and masks</a></H2>"

}

puts $greport "<hr><H2>Higher-level reports</H2>"

fsl:exec "${FSLDIR}/bin/avwmaths mask -Tmin mask"

#}}}
    #{{{ setup mean_func image

set mean_list ""

for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
    set find_std [ feat5:find_std $feat_files($session) mean_func ]
    if { $find_std != 0 } {
	set mean_list "$mean_list $find_std"
    }
}

if { [ string length $mean_list ] == 0 } {
    fsl:exec "${FSLDIR}/bin/avwmaths_32R mask -bin -mul 10000 mean_func"
} else {
    fsl:exec "${FSLDIR}/bin/avwmerge -t mean_func $mean_list"
    fsl:exec "${FSLDIR}/bin/avwmaths mean_func -Tmean mean_func"
}

#}}}

    if { $fmri(inputtype) == 1 } {
	#{{{ create 4D inputs for "FEAT directories" option

for { set nci 1 } { $nci <=  $fmri(ncopeinputs) } { incr nci 1 } {   

    if { $fmri(copeinput.$nci) } {

	set contrastname($nci) "[ feat5:getconname $feat_files(1) $nci ]"

	set cope_list ""
	set mean_lcon 0
	for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
	    set find_std [ feat5:find_std $feat_files($session) stats/cope$nci ]
	    if { $find_std != 0 } {
		set cope_list "$cope_list $find_std"
	    } else {
		fsl:echo $logout "Error - not all input FEAT directories have the same set of contrasts! (currently looking in $feat_files($session) for cope${nci})"
	        exit 1
            }
	    if { [ file exists $feat_files($session)/design.con ] } {
		set awknum [ expr 1 + $nci ]
		set mean_lcon_incr [ exec sh -c "grep PPheights $feat_files($session)/design.con | awk '{ print \$$awknum }'" ]
		if { [ file exists $feat_files($session)/design.lcon ] } {
		    set mean_lcon_incr [ expr $mean_lcon_incr * [ exec sh -c "cat $feat_files($session)/design.lcon" ] ]
		}
		set mean_lcon [ expr $mean_lcon + $mean_lcon_incr ]
	    } else {
		fsl:echo $logout "Error - not all input FEAT directories have valid design.con contrast files."
		exit 1
	    }
        }

	fsl:exec "${FSLDIR}/bin/avwmerge -t cope$nci $cope_list"
	fsl:exec "${FSLDIR}/bin/avwmaths cope$nci -mas mask cope$nci"

	set mean_lcon [ expr $mean_lcon / $fmri(multiple) ]
	if { $mean_lcon < 0.01 } {
	    set mean_lcon 1
	}
	fsl:exec "echo -n '$mean_lcon ' >> design.lcon"

	set varcope_list ""
	for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
	    set varcope_list "$varcope_list [ feat5:find_std $feat_files($session) stats/varcope$nci ]"
	}
	fsl:exec "${FSLDIR}/bin/avwmerge -t varcope$nci $varcope_list"
	fsl:exec "${FSLDIR}/bin/avwmaths varcope$nci -mas mask varcope$nci"
	
	#{{{ setup t_dof

	set tdof_list ""
	for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
	    if { $fmri(mixed_yn) != 3 || [ file exists $feat_files($session)/design.lev ] } {
		set find_std [ feat5:find_std $feat_files($session) stats/tdof_t$nci ]
		if { $find_std != 0 } {
		    set tdof_list "$tdof_list $find_std"
		}
	    } else {
		set THEDOF [ exec sh -c "cat $feat_files($session)/stats/dof" ]
		fsl:exec "${FSLDIR}/bin/avwmaths $feat_files($session)/reg_standard/stats/cope$nci -mul 0 -add $THEDOF $feat_files($session)/reg_standard/stats/FEtdof_t$nci"
		set tdof_list "$tdof_list $feat_files($session)/reg_standard/stats/FEtdof_t$nci"
	    }
	}
	if { $tdof_list != "" } {
	    fsl:exec "${FSLDIR}/bin/avwmerge -t tdof_t$nci $tdof_list"
	    fsl:exec "${FSLDIR}/bin/avwmaths tdof_t$nci -mas mask tdof_t$nci"
	}

#}}}

    } else {
	fsl:exec "echo -n '1 ' >> design.lcon"
    }
}

if { $fmri(sscleanup_yn) } {
    for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
	fsl:exec "${FSLDIR}/bin/featregapply $feat_files($session) -c"
    }
}

for { set nci 1 } { $nci <=  $fmri(ncopeinputs) } { incr nci 1 } {   
    if { $fmri(copeinput.$nci) } {
	set feat_files($nci) ${FD}/cope$nci
    } else {
	set feat_files($nci) -1
    }
}
set fmri(multiple) $fmri(ncopeinputs)

#}}}
    } elseif { $fmri(inputtype) == 2 } {
	#{{{ create 4D inputs for "3D cope input" option

set cope_list ""
set varcope_list ""
set tdof_list ""
set mean_lcon 0

for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
    set copename [ file tail $copes($session) ]
    set copenumber [ remove_ext [ string trimleft $copename "cope" ] ]

    set cope_list "$cope_list [ feat5:find_std $feat_files($session) stats/$copename ]"
    set varcope_list "$varcope_list [ feat5:find_std $feat_files($session) stats/var$copename ]"

    if { [ file exists $feat_files($session)/design.con ] } {
	set awknum [ expr 1 + $copenumber ]
	set mean_lcon_incr [ exec sh -c "grep PPheights $feat_files($session)/design.con | awk '{ print \$$awknum }'" ]
	if { [ file exists $feat_files($session)/design.lcon ] } {
	    set mean_lcon_incr [ expr $mean_lcon_incr * [ exec sh -c "cat $feat_files($session)/design.lcon" ] ]
	}
	set mean_lcon [ expr $mean_lcon + $mean_lcon_incr ]
    }

    if { $fmri(mixed_yn) != 3 || [ file exists $feat_files(1)/design.lev ] } {
	set find_std [ feat5:find_std $feat_files($session) stats/tdof_t$copenumber ]
	if { $find_std != 0 } {
	    set tdof_list "$tdof_list $find_std"
	}
    } else {
	set THEDOF [ exec sh -c "cat $feat_files($session)/stats/dof" ]
	fsl:exec "${FSLDIR}/bin/avwmaths $feat_files($session)/reg_standard/stats/$copename -mul 0 -add $THEDOF $feat_files($session)/reg_standard/stats/FEtdof_t${copenumber}"
	set tdof_list "$tdof_list $feat_files($session)/reg_standard/stats/FEtdof_t${copenumber}"
    }
}

fsl:exec "${FSLDIR}/bin/avwmerge -t cope1 $cope_list"
fsl:exec "${FSLDIR}/bin/avwmaths cope1 -mas mask cope1"

set mean_lcon [ expr $mean_lcon / $fmri(multiple) ]
if { $mean_lcon < 0.01 } {
    set mean_lcon 1
}
fsl:exec "echo $mean_lcon > design.lcon"

fsl:exec "${FSLDIR}/bin/avwmerge -t varcope1 $varcope_list"
fsl:exec "${FSLDIR}/bin/avwmaths varcope1 -mas mask varcope1"

if { $tdof_list != "" } {
    fsl:exec "${FSLDIR}/bin/avwmerge -t tdof_t1 $tdof_list"
    fsl:exec "${FSLDIR}/bin/avwmaths tdof_t1 -mas mask tdof_t1"
}

if { $fmri(sscleanup_yn) } {
    for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
	fsl:exec "${FSLDIR}/bin/featregapply $feat_files($session) -c"
    }
}

set contrastname(1) ""
set fmri(multiple) 1
set feat_files(1) ${FD}/cope1

#}}}
    }

    #{{{ finish

fsl:echo $logout "Starting higher-level stats..."

after 2000

catch { exec kill -9 $tailpid } errmsg

#}}}
}

#}}}

    for { set session 1 } { $session <= $fmri(multiple) } { incr session 1 } {
	if { $fmri(level) == 1 || $feat_files($session) != -1 } {
	
	    #{{{ setup FEAT directory

set FD $feat_files($session)
set thedate [ exec date ]

if { $fmri(analysis) == 0 || $fmri(analysis) == 4 } {
    # 0 = No first-level analysis (registration and/or group stats only)
    # 4 = Contrasts, Thresholding, Rendering
    #{{{ copy old featdir if required, and backup old files inside

if { $fmri(newdir_yn) } {
    set OldFD $feat_files($session)
    
    set FD [ new_filename $OldFD ]
    
    if { ! [ file writable [ file dirname $FD ] ] } {
	set FD [ new_filename ${HOME}/[ file tail $OldFD ] ]
    }

    file copy $OldFD $FD
}

cd $FD

set BACKUPS [ new_filename old ]
fsl:exec "/bin/mkdir $BACKUPS ; mv report.log report.com design* $BACKUPS ; /bin/cp $BACKUPS/design.fsf $BACKUPS/design.lev ." -n

if { $fmri(analysis) == 4 } {
    fsl:exec "/bin/mv ?endered_thresh_* cluster_* lmax_* prob_mask_* thresh_* tsplot stats/cope* stats/neff* stats/tstat* stats/varcope* stats/zstat* stats/fstat* stats/zfstat* $BACKUPS" -n
}

#}}}
} else {
  # 7 : Full first-level analysis
  # 1 : Pre-Stats
  # 3 : Pre-Stats + Stats
  # 2 :             Stats
  # 6 :             Stats + Contrasts, Thresholding, Rendering
  #{{{ setup new featdir name

if { $fmri(level) == 1 && $fmri(outputdir) != "" } {
    if { $fmri(multiple) == 1 } {
	set featname [ feat5:strip $fmri(outputdir) ].feat
    } else {
	set featname [ file rootname $feat_files($session) ]_[ file rootname [ file tail $fmri(outputdir) ] ].feat
    }
} else {
    set featname [ feat5:strip [ remove_ext $feat_files($session) ] ].feat
}
#set FD [ new_filename $featname ]; ; # This adds the +'s to the feat dir if it exists
set FD $featname; # This does not create new file name
if { ! [ file writable [ file dirname $FD ] ] } {
    set FD [ new_filename ${HOME}/[ file tail $featname ] ]
}

fsl:exec "/bin/mkdir $FD" -n

set funcdata [ remove_ext $feat_files($session) ]

#}}}
}

cd $FD
set FD [ pwd ]
set feat_files($session) $FD
set logout ${FD}/report.log
set comout ${FD}/report.com

fsl:echo $logout "Started FEAT at $thedate on $HOSTNAME
FEAT output is ${FD}"

if { $fmri(featwatcher_yn) } {
    catch { exec sh -c "${FSLDIR}/bin/Featwatcher${gui_ext} $FD" >& /dev/null & } errmsg
} else {
    if { $OSFLAVOUR != "cygwin" } {
	set tailpid [ exec tail -f report.log & ]
    }
}

fsl:exec "/bin/cp ${fsfroot}.fsf design.fsf"

#{{{ process relative filenames

# remove this - it's defunct......
# if { $fmri(relative_yn) && $session > 1 } {

#     set fix_design [ open design.fsf "a" ]

#     set datarootdir [ file dirname $feat_files($session) ]

#     for { set i 1 } { $i <= $fmri(evs_orig) } { incr i 1 } {
# 	if { $fmri(shape$i) > 2 } {
# 	    set fmri(custom$i) ${datarootdir}/[ file tail $fmri(custom$i)]
# 	    puts $fix_design "
# # Custom EV file (EV $i)
# set fmri(custom$i) \"$fmri(custom$i)\""
#         }
#     }

#     close $fix_design
# }

#}}}

set motreg ""
if { $fmri(motionevs) > 0 && [ file exists mc/prefiltered_func_data_mcf.par ] } {
    set motreg "mc/prefiltered_func_data_mcf.par"
}
fsl:exec "${FSLDIR}/bin/feat_model design $motreg"
set fsfroot ${FD}/design

#{{{ setup higher-level stuff

if { $fmri(level)>1 } {
    if { $fmri(analysis)!=4 } {
	fsl:echo design.lev "$contrastname($session)"
	fsl:exec "cat ../design.lcon | awk '{ print \$$session }' > design.lcon"
    } else {
	set contrastname($session) [ fsl:exec "cat design.lev" ]
    }
}

#}}}

#}}}

	    if { $fmri(analysis) > 0 } {
		#{{{ write start of report file

set conname ""
if { $fmri(level) > 1 } {
    set conname " for Lower-level Contrast $session"
    if { $contrastname($session) != "" } {
	set conname "$conname ($contrastname($session))"
    }
}

# if this is a re-run of an existing featdir, find the original FEAT version
# (if it is older than 5.63 we'll restart the report.html from scratch as the
# poststats insertion stuff doesn't work well otherwise)
set OLDFEAT 0
set orig_feat_version $fmri(version)
if { [ file exists report.html ] } {
    set orig_feat_version [ fsl:exec "grep 'This page produced' report.html | awk '{ print \$8 }' -" ]
}
if { [ string is double $orig_feat_version ] && $orig_feat_version > 4 && $orig_feat_version < 5.63 } {
    set OLDFEAT 1
}

if { $OLDFEAT || $fmri(analysis)!=4 } {

    set report [ open report.html "w" ]
    fconfigure $report -buffering none

    puts $report "<HTML><TITLE>FEAT Report</TITLE><BODY BACKGROUND=\"file:${FSLSLASH}${FSLDIR}/doc/images/fsl-bg.jpg\">
<hr><CENTER><H1>FEAT Report$conname</H1>
${FD}/report.html<br>$thedate</CENTER>"

    set ps "<hr><b>Analysis methods</b>
<p>Analysis was carried out using FEAT (FMRI Expert Analysis Tool) Version $fmri(version), part of FSL (FMRIB's Software Library, www.fmrib.ox.ac.uk/fsl)."

    set rs "<p><b>References</b><br>
"
}

#}}}

		if { $fmri(level) == 1 } {
		    #{{{ pre-stats

if { $fmri(analysis) != 4 } {

    puts $report "<!--prestatsstart-->"

    #{{{ check npts, delete images, make example_func

# check npts
set total_volumes [ exec sh -c "${FSLDIR}/bin/avwnvols $funcdata 2> /dev/null" ]
fsl:echo $logout "Total original volumes = $total_volumes"
if { $total_volumes != $fmri(npts) } {
    fsl:echo $logout "Error - $funcdata has a different number of time points to that in FEAT setup"
    return 1
}

# delete images
if { $fmri(ndelete) > 0 } {
    fsl:echo $logout "Deleting $fmri(ndelete) volume(s) - BE WARNED for future analysis!"
    set total_volumes [ expr $total_volumes - $fmri(ndelete) ]
    fsl:exec "${FSLDIR}/bin/avwroi $funcdata prefiltered_func_data $fmri(ndelete) $total_volumes"
    set funcdata prefiltered_func_data
}

# choose halfway image and copy to example_func (unless alternative example_func setup)
set target_vol_number [ expr $total_volumes / 2 ]
if { [ imtest $fmri(alternative_example_func) ] } {
    fsl:exec "${FSLDIR}/bin/avwmaths $fmri(alternative_example_func) example_func"
} else {
    fsl:exec "${FSLDIR}/bin/avwroi $funcdata example_func $target_vol_number 1"
} 

#}}}

    if { $fmri(filtering_yn) } {

	set ps "$ps The following pre-statistics processing was applied"

	#{{{ motion correction

#mc: 0=none 1=MCFLIRT

if { $fmri(mc) != 0 } {
    set ps "$ps; motion correction using MCFLIRT \[Jenkinson 2002\]"
    set rs "$rs\[<a href=\"http://www.fmrib.ox.ac.uk/analysis/techrep/#TR02MJ1\">Jenkinson 2002</a>\] M. Jenkinson and P. Bannister and M. Brady and S. Smith. Improved optimisation for the robust and accurate linear registration and motion correction of brain images. NeuroImage 17:2(825-841) 2002.<br>
    "

    fsl:exec "${FSLDIR}/bin/mcflirt -in $funcdata -out prefiltered_func_data_mcf -mats -plots -refvol $target_vol_number -rmsrel -rmsabs" -r
    if { ! $fmri(regunwarp_yn) } {
	set funcdata prefiltered_func_data_mcf
    }

    fsl:exec "/bin/mkdir mc ; /bin/mv -f prefiltered_func_data_mcf.mat prefiltered_func_data_mcf.par prefiltered_func_data_mcf_abs.rms prefiltered_func_data_mcf_abs_mean.rms prefiltered_func_data_mcf_rel.rms prefiltered_func_data_mcf_rel_mean.rms mc"
    cd mc
    
    #{{{ make plots

fsl:exec "echo \"set term pbm color ; set size 1,0.3 ; set title 'MCFLIRT estimated rotations (radians)' ; plot 'prefiltered_func_data_mcf.par' using 1 title 'x' with lines , 'prefiltered_func_data_mcf.par' using 2 title 'y' with lines , 'prefiltered_func_data_mcf.par' using 3 title 'z' with lines\" | $FSLGNUPLOT > rot.ppm"
fsl:exec "$FSLCONVERT rot.ppm rot.gif"

fsl:exec "echo \"set term pbm color ; set size 1,0.3 ; set title 'MCFLIRT estimated translations (mm)' ; plot 'prefiltered_func_data_mcf.par' using 4 title 'x' with lines , 'prefiltered_func_data_mcf.par' using 5 title 'y' with lines , 'prefiltered_func_data_mcf.par' using 6 title 'z' with lines\" | $FSLGNUPLOT > trans.ppm"
fsl:exec "$FSLCONVERT trans.ppm trans.gif"

fsl:exec "echo \"set term pbm color ; set size 1,0.3 ; set title 'MCFLIRT estimated mean displacement (mm)' ; plot 'prefiltered_func_data_mcf_abs.rms' title 'absolute' with lines , 'prefiltered_func_data_mcf_rel.rms' title 'relative' with lines\" | $FSLGNUPLOT > disp.ppm"
fsl:exec "$FSLCONVERT disp.ppm disp.gif"

fsl:exec "/bin/rm -f rot.ppm trans.ppm disp.ppm"

#}}}
    #{{{ extract mean displacements

set mcchannel [ open prefiltered_func_data_mcf_abs_mean.rms "r" ]
gets $mcchannel line
scan $line "%f" absrms
set absrms [ expr int($absrms*100.0)/100.0 ]
close $mcchannel

set mcchannel [ open prefiltered_func_data_mcf_rel_mean.rms "r" ]
gets $mcchannel line
scan $line "%f" relrms
set relrms [ expr int($relrms*100.0)/100.0 ]
close $mcchannel

#}}}
    #{{{ web page report

set mcreport [ open index.html "w" ]

puts $mcreport "<HTML>

<TITLE>FEAT Motion Correction Report</TITLE>

<BODY BACKGROUND=\"file:${FSLSLASH}${FSLDIR}/doc/images/fsl-bg.jpg\">

<hr><CENTER>
<H1>FEAT Motion Correction Report</H1>
for FEAT session &nbsp; &nbsp; <a href=\"../report.html\">${FD}</a><br>
[ exec date ]
<hr>

<p><IMG BORDER=0 SRC=\"rot.gif\">
<p><IMG BORDER=0 SRC=\"trans.gif\">
<p><IMG BORDER=0 SRC=\"disp.gif\">

<p>Mean (across voxels) voxel displacements:<br>
absolute (each time point with respect to the reference image) = ${absrms}mm;<br>
relative (each time point with respect to the previous timepoint) = ${relrms}mm.

</CENTER>
<HR><FONT SIZE=1>This page produced automatically by FEAT - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL</A>.</FONT>

</BODY></HTML>
"

close $mcreport

#}}}
    #{{{ main web report output

#set pixdimx [ expr abs([ exec sh -c "$FSLDIR/bin/avwval example_func pixdim1" ]) ]
#set pixdimy [ expr abs([ exec sh -c "$FSLDIR/bin/avwval example_func pixdim2" ]) ]
#set pixdimz [ expr abs([ exec sh -c "$FSLDIR/bin/avwval example_func pixdim3" ]) ]
#set pixdim [ expr  pow($pixdimx * $pixdimy * $pixdimz,0.33333) ]
#
#set MC_TOLERANCE [ expr 0.5 * $pixdim / 5.0 ]

set MC_TOLERANCE 0.5

set mcwarning ""

if { $relrms > $MC_TOLERANCE } {
    set mcwarning " - warning - high levels of motion detected"
}

puts $report "<hr><a href=\"mc/index.html\">Motion correction report</a> (mean displacements: absolute=${absrms}mm, relative=${relrms}mm)$mcwarning"

#}}}

    cd $FD
}

#}}}
	#{{{ B0 unwarping

if { $fmri(regunwarp_yn) } {

    fsl:exec "/bin/mkdir unwarp"
    set OLDPWD [ pwd ]
    cd unwarp

    #{{{ start web page report

set unwarpreport [ open index.html "w" ]

puts $unwarpreport "<HTML>

<TITLE>FEAT / FUGUE Fieldmap Unwarping Report</TITLE>

<BODY BACKGROUND=\"file:${FSLSLASH}${FSLDIR}/doc/images/fsl-bg.jpg\">

<hr><CENTER>
<H1>FEAT / FUGUE Fieldmap Unwarping Report</H1>
for FEAT session &nbsp; &nbsp; <a href=\"../report.html\">${FD}</a><br>
[ exec date ]
</CENTER>
"

#}}}
    #{{{ do the unwarping calculations

    # FM = space of fieldmap
    # EF = space of example_func
    # UD = undistorted (in any space)
    # D  = distorted (in any space)

    # copy in unwarp input files into reg subdir
    fsl:exec "${FSLDIR}/bin/avwmaths ../example_func EF_D_example_func"
    fsl:exec "${FSLDIR}/bin/avwmaths $unwarp_files($session) FM_UD_fmap"
    fsl:exec "${FSLDIR}/bin/avwmaths $unwarp_files_mag($session) FM_UD_fmap_mag"

    # brain-extract fmap_mag and de-mean fmap
    if { $fmri(bet_yn) } {
	fsl:exec "${FSLDIR}/bin/bet FM_UD_fmap_mag FM_UD_fmap_mag_brain -m"
    } else {
	fsl:exec "${FSLDIR}/bin/avwmaths FM_UD_fmap_mag FM_UD_fmap_mag_brain"
	fsl:exec "${FSLDIR}/bin/avwmaths FM_UD_fmap_mag -bin FM_UD_fmap_mag_brain_mask"
    }
    fsl:exec "${FSLDIR}/bin/avwmaths FM_UD_fmap -sub [ fsl:exec "${FSLDIR}/bin/avwstats++ FM_UD_fmap -k FM_UD_fmap_mag_brain_mask -P 50" ] FM_UD_fmap"

    # create report picture of fmap overlaid onto whole-head mag image
    set fmapmin [ fsl:exec "${FSLDIR}/bin/avwstats++ FM_UD_fmap -R | awk '{ print \$1 }'" ]
    fsl:exec "${FSLDIR}/bin/avwmaths FM_UD_fmap -sub $fmapmin -add 10 -mas FM_UD_fmap_mag_brain_mask grot"
    set fmapminmax [ fsl:exec "${FSLDIR}/bin/avwstats++ grot -l 1 -p 0.1 -p 95" ]
    fsl:exec "${FSLDIR}/bin/overlay 0 0 FM_UD_fmap_mag -a grot $fmapminmax grot"
    fsl:exec "${FSLDIR}/bin/slicer grot -s 3 -x 0.35 sla -x 0.45 slb -x 0.55 slc -x 0.65 sld -y 0.35 sle -y 0.45 slf -y 0.55 slg -y 0.65 slh -z 0.35 sli -z 0.45 slj -z 0.55 slk -z 0.65 sll ; $FSLCONVERT -colors 100 +append sla slb slc sld sle slf slg slh sli slj slk sll fmap+mag.gif"
    puts $unwarpreport "<hr><p>Brain-masked B0 fieldmap in colour, overlaid on top of fieldmap magnitude image<br>
<a href=\"fmap+mag.gif\"><IMG BORDER=0 SRC=\"fmap+mag.gif\" WIDTH=1200></a><br><br>"

    # get a sigloss estimate and make a siglossed mag for forward warp
    set epi_te [ expr $fmri(te) / 1000.0 ]
    fsl:exec "${FSLDIR}/bin/sigloss -i FM_UD_fmap --te=$epi_te -m FM_UD_fmap_mag_brain_mask -s FM_UD_fmap_sigloss"
    set siglossthresh [ expr 1.0 - ( $fmri(signallossthresh) / 100.0 ) ]
    fsl:exec "${FSLDIR}/bin/avwmaths_32R FM_UD_fmap_sigloss -mul FM_UD_fmap_mag_brain FM_UD_fmap_mag_brain_siglossed"

    # make a warped version of FM_UD_fmap_mag to match with the EPI
    set dwell [ expr $fmri(dwell) / 1000.0 ]
    fsl:exec "${FSLDIR}/bin/fugue -i FM_UD_fmap_mag_brain_siglossed --loadfmap=FM_UD_fmap --mask=FM_UD_fmap_mag_brain_mask --dwell=$dwell -w FM_D_fmap_mag_brain_siglossed --nokspace --unwarpdir=$fmri(unwarp_dir)"
    fsl:exec "${FSLDIR}/bin/fugue -i FM_UD_fmap_sigloss             --loadfmap=FM_UD_fmap --mask=FM_UD_fmap_mag_brain_mask --dwell=$dwell -w FM_D_fmap_sigloss             --nokspace --unwarpdir=$fmri(unwarp_dir)"
    fsl:exec "${FSLDIR}/bin/avwmaths FM_D_fmap_sigloss -thr $siglossthresh FM_D_fmap_sigloss"
    fsl:exec "${FSLDIR}/bin/flirt -in EF_D_example_func -ref FM_D_fmap_mag_brain_siglossed -omat EF_2_FM.mat -o grot -dof 6 -refweight FM_D_fmap_sigloss"
    fsl:exec "${FSLDIR}/bin/convert_xfm -omat FM_2_EF.mat -inverse EF_2_FM.mat"

    # put fmap stuff into space of EF_D_example_func
    fsl:exec "${FSLDIR}/bin/flirt -in FM_UD_fmap                -ref EF_D_example_func -init FM_2_EF.mat -applyxfm -out EF_UD_fmap"
    fsl:exec "${FSLDIR}/bin/flirt -in FM_UD_fmap_mag_brain      -ref EF_D_example_func -init FM_2_EF.mat -applyxfm -out EF_UD_fmap_mag_brain"
    fsl:exec "${FSLDIR}/bin/flirt -in FM_UD_fmap_mag_brain_mask -ref EF_D_example_func -init FM_2_EF.mat -applyxfm -out EF_UD_fmap_mag_brain_mask"
    fsl:exec "${FSLDIR}/bin/flirt -in FM_UD_fmap_sigloss        -ref EF_D_example_func -init FM_2_EF.mat -applyxfm -out EF_UD_fmap_sigloss"
    fsl:exec "${FSLDIR}/bin/avwmaths_32R FM_UD_fmap_mag_brain_mask -thr 0.5 -bin FM_UD_fmap_mag_brain_mask"
    fsl:exec "${FSLDIR}/bin/avwmaths_32R EF_UD_fmap_sigloss -thr $siglossthresh EF_UD_fmap_sigloss"

    # create report pic for sigloss
    fsl:exec "${FSLDIR}/bin/overlay 1 0 FM_UD_fmap_mag_brain -a EF_UD_fmap_sigloss 0 1 grot"
    fsl:exec "${FSLDIR}/bin/slicer grot -s 3 -x 0.35 sla -x 0.45 slb -x 0.55 slc -x 0.65 sld -y 0.35 sle -y 0.45 slf -y 0.55 slg -y 0.65 slh -z 0.35 sli -z 0.45 slj -z 0.55 slk -z 0.65 sll ; $FSLCONVERT -colors 100 +append sla slb slc sld sle slf slg slh sli slj slk sll EF_UD_sigloss+mag.gif"
    puts $unwarpreport "Thresholded signal loss weighting image<br>
<a href=\"EF_UD_sigloss+mag.gif\"><IMG BORDER=0 SRC=\"EF_UD_sigloss+mag.gif\" WIDTH=1200></a><br><br>"

    # apply warp to EF_D_example_func and save unwarp-shiftmap then convert to unwarp-warpfield
    fsl:exec "${FSLDIR}/bin/fugue --loadfmap=EF_UD_fmap --dwell=$dwell --mask=EF_UD_fmap_mag_brain_mask -i EF_D_example_func -u EF_UD_example_func --unwarpdir=$fmri(unwarp_dir) --saveshift=EF_UD_shift"
    fsl:exec "${FSLDIR}/bin/convertwarp -s EF_UD_shift -o EF_UD_warp -r EF_D_example_func --shiftdir=$fmri(unwarp_dir)"

    # create report pic for shift extent
    set shiftminmax [ fsl:exec "${FSLDIR}/bin/avwstats++ EF_UD_shift -R -P 1 -P 99" ]
    set shiftminR [ format %.1f [ lindex $shiftminmax 0 ] ]
    set shiftmaxR [ format %.1f [ lindex $shiftminmax 1 ] ]
    set shiftminr [ expr [ lindex $shiftminmax 2 ] * -1.0 ]
    set shiftmaxr [ lindex $shiftminmax 3 ]
    fsl:exec "${FSLDIR}/bin/avwmaths EF_UD_shift -mul -1 grot"
    fsl:exec "${FSLDIR}/bin/overlay 1 0 FM_UD_fmap_mag_brain -a EF_UD_shift 0 $shiftmaxr grot 0 $shiftminr grot"
    fsl:exec "${FSLDIR}/bin/slicer grot -s 3 -x 0.35 sla -x 0.45 slb -x 0.55 slc -x 0.65 sld -y 0.35 sle -y 0.45 slf -y 0.55 slg -y 0.65 slh -z 0.35 sli -z 0.45 slj -z 0.55 slk -z 0.65 sll ; $FSLCONVERT -colors 100 +append sla slb slc sld sle slf slg slh sli slj slk sll EF_UD_shift+mag.gif"
    fsl:exec "/bin/cp ${FSLDIR}/etc/luts/ramp.gif .ramp.gif"
    fsl:exec "/bin/cp ${FSLDIR}/etc/luts/ramp2.gif .ramp2.gif"
puts $unwarpreport "Unwarping shift map, in voxels &nbsp;&nbsp;&nbsp; ${shiftminR} <IMG BORDER=0 SRC=\".ramp2.gif\"> 0 <IMG BORDER=0 SRC=\".ramp.gif\"> ${shiftmaxR}<br>
<a href=\"EF_UD_shift+mag.gif\"><IMG BORDER=0 SRC=\"EF_UD_shift+mag.gif\" WIDTH=1200></a><br><br>"

    # create report pics in EF space
    fsl:exec "${FSLDIR}/bin/slicer EF_D_example_func    -s 3 -x 0.35 sla -x 0.45 slb -x 0.55 slc -x 0.65 sld -y 0.35 sle -y 0.45 slf -y 0.55 slg -y 0.65 slh -z 0.35 sli -z 0.45 slj -z 0.55 slk -z 0.65 sll ; $FSLCONVERT -colors 100 +append sla slb slc sld sle slf slg slh sli slj slk sll EF_D_example_func.gif"
    fsl:exec "${FSLDIR}/bin/slicer EF_UD_example_func   -s 3 -x 0.35 sla -x 0.45 slb -x 0.55 slc -x 0.65 sld -y 0.35 sle -y 0.45 slf -y 0.55 slg -y 0.65 slh -z 0.35 sli -z 0.45 slj -z 0.55 slk -z 0.65 sll ; $FSLCONVERT -colors 100 +append sla slb slc sld sle slf slg slh sli slj slk sll EF_UD_example_func.gif"
    fsl:exec "${FSLDIR}/bin/slicer EF_UD_fmap_mag_brain -s 3 -x 0.35 sla -x 0.45 slb -x 0.55 slc -x 0.65 sld -y 0.35 sle -y 0.45 slf -y 0.55 slg -y 0.65 slh -z 0.35 sli -z 0.45 slj -z 0.55 slk -z 0.65 sll ; $FSLCONVERT -colors 100 +append sla slb slc sld sle slf slg slh sli slj slk sll EF_UD_fmap_mag_brain.gif"
    fsl:exec "${FSLDIR}/bin/gifmerge -50 -l0 EF_D_example_func.gif EF_UD_example_func.gif                          > EF_UD_movie2.gif"
    fsl:exec "${FSLDIR}/bin/gifmerge -50 -l0 EF_D_example_func.gif EF_UD_example_func.gif EF_UD_fmap_mag_brain.gif > EF_UD_movie3.gif ; /bin/rm -f sla slb slc sld sle slf slg slh sli slj slk sll grot*"
    puts $unwarpreport "Original distorted example_func (example_func_orig_distorted)<br>
<a href=\"EF_D_example_func.gif\"><IMG BORDER=0 SRC=\"EF_D_example_func.gif\" WIDTH=1200></a><br><br>
Undistorted example_func (example_func)<br>
<a href=\"EF_UD_example_func.gif\"><IMG BORDER=0 SRC=\"EF_UD_example_func.gif\" WIDTH=1200></a><br><br>
Non-distorted fieldmap magnitude brain-extracted image in space of example_func (unwarp/EF_UD_fmap_mag_brain)<br>
<a href=\"EF_UD_fmap_mag_brain.gif\"><IMG BORDER=0 SRC=\"EF_UD_fmap_mag_brain.gif\" WIDTH=1200></a><br><br>
Movie of all 3 images.<br>
<a href=\"EF_UD_movie3.gif\"><IMG BORDER=0 SRC=\"EF_UD_movie3.gif\" WIDTH=1200></a><br><br>
"

#}}}
    #{{{ end web page report

puts $unwarpreport "
<HR><FONT SIZE=1>This page produced automatically by FEAT - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL</A>.</FONT>
</BODY></HTML>
"

close $unwarpreport

#}}}
    #{{{ main FEAT web page stuff

set ps "$ps; fieldmap-based EPI unwarping using PRELUDE+FUGUE \[Jenkinson 2003, 2004\]"
set rs "$rs\[<a href=\"http://www.fmrib.ox.ac.uk/analysis/techrep/#TR01MJ1\">Jenkinson 2003</a>\] M. Jenkinson. A fast, automated, n-dimensional phase unwrapping algorithm. Magnetic Resonance in Medicine 49(1):193-197 2003.<br>
\[<a href=\"http://www.fmrib.ox.ac.uk/~mark/work/hbm2004.ps\">Jenkinson 2004</a>\] M. Jenkinson. Improving the registration of B0-disorted EPI images using calculated cost function weights. Tenth Int. Conf. on Functional Mapping of the Human Brain 2004.<br>
"

puts $report "<!--unwarpstart-->
<hr><a href=\"unwarp/index.html\">Fieldmap unwarping report<br><br>
Comparison of original (distorted) and unwarped example_func<br>
<IMG BORDER=0 SRC=\"unwarp/EF_UD_movie2.gif\" WIDTH=1000></a>
<!--unwarphrstart-->
<!--unwarphrstop-->
<!--unwarpstop-->"

#}}}
    #{{{ apply warping and motion correction to example_func and 4D data

cd $OLDPWD

immv example_func example_func_orig_distorted
fsl:exec "${FSLDIR}/bin/applywarp -i example_func_orig_distorted -o example_func -w unwarp/EF_UD_warp -r example_func_orig_distorted --abs --mask=unwarp/EF_UD_fmap_mag_brain_mask"

# now either apply unwarping one vol at a time (including applying individual mcflirt transforms at same time),
# or if mcflirt transforms don't exist, just apply warp to 4D $funcdata
if { [ file exists mc/prefiltered_func_data_mcf.mat/MAT_0000 ] } {
    for { set i 0 } { $i < $total_volumes } { incr i 1 } {
	set pad [format %04d $i]
	fsl:exec "${FSLDIR}/bin/avwroi $funcdata grot$pad $i 1"
	fsl:exec "${FSLDIR}/bin/applywarp -i grot$pad -o grot$pad --premat=mc/prefiltered_func_data_mcf.mat/MAT_$pad -w unwarp/EF_UD_warp -r example_func --abs --mask=unwarp/EF_UD_fmap_mag_brain_mask"
    }
    fsl:exec "${FSLDIR}/bin/avwmerge -t prefiltered_func_data_unwarp [ imglob -oneperimage grot* ]"
    fsl:exec "/bin/rm -f grot*"
} else {
    fsl:exec "${FSLDIR}/bin/applywarp -i $funcdata -o prefiltered_func_data_unwarp -w unwarp/EF_UD_warp -r example_func --abs --mask=unwarp/EF_UD_fmap_mag_brain_mask"
}

set funcdata prefiltered_func_data_unwarp

#}}}
    #{{{ finish warp-related stuff (NOT NEEDED)

# # to be run AFTER all the standard-space (etc) registration
#
# # generate warp field for total resampling
# fsl:exec "${FSLDIR}/bin/flirt -in EF_UD_shift -ref standard -out grot -applyxfm"
# # compensate for change in voxel size (as shiftmap units are voxels!)
# if { [ string first x $fmri(unwarp_dir) ] == 0 } {
#     set DIM 1
# } elseif { [ string first y $fmri(unwarp_dir) ] == 0 } {
#     set DIM 2
# } else {
#     set DIM 3
# }
# set VOXSIZERATIO [ expr abs( [ exec sh -c "$FSLDIR/bin/avwval example_func pixdim$DIM" ] / [ exec sh -c "$FSLDIR/bin/avwval standard pixdim$DIM" ] ) ]
# fsl:exec "${FSLDIR}/bin/avwmaths grot -mul $VOXSIZERATIO grot"
# # must convert the shiftmap and apply the postaff separately!
# fsl:exec "${FSLDIR}/bin/convertwarp -s grot -o grot -r standard --shiftdir=$fmri(unwarp_dir)"
# fsl:exec "${FSLDIR}/bin/convertwarp -w grot --postmat=${ef}2standard.mat -r standard -o SS_UD_example_func2standard_warp"

# # annoying mask generation to get rid of edge effects in apply warp
# fsl:exec "${FSLDIR}/bin/avwmaths SS_UD_example_func2standard_warp -mul 0 -add 1 grot"
# fsl:exec "${FSLDIR}/bin/convertwarp -w grot --postmat=${ef}2standard.mat -r standard -o grot"
# fsl:exec "${FSLDIR}/bin/avwmaths grot -Tmean -thr 0.9999 -bin SS_UD_warpmask"

# # example all-in-one unwarping+upsampling
# #fsl:exec "${FSLDIR}/bin/applywarp -i ../example_func -o SS_UD_example_func -w SS_UD_example_func2standard_warp -r standard --abs"
# #fsl:exec "${FSLDIR}/bin/avwmaths SS_UD_example_func -mas SS_UD_warpmask SS_UD_example_func"

#}}}
}

#}}}
	#{{{ slice timing correction

if { $fmri(st) > 0 } {

    set ps "$ps; slice-timing correction using Fourier-space time-series phase-shifting"

    set st_opts ""
    
    switch $fmri(st) {
	2 {
	    set st_opts "--down"
	}
	3 {
	    set st_opts "--ocustom=$fmri(st_file)"
	}
	4 {
	    set st_opts "--tcustom=$fmri(st_file)"
	}
	5 {
	    set st_opts "--odd"
	}
    }

    fsl:exec "${FSLDIR}/bin/slicetimer -i $funcdata --out=prefiltered_func_data_st -r $fmri(tr) $st_opts"
    set funcdata prefiltered_func_data_st
}

#}}}
	#{{{ BET

if { $fmri(bet_yn) } {

    if { [ info exists fmri(alternative_mask) ] && [ imtest $fmri(alternative_mask) ] } {
        puts "-------Using alternative mask---------"
	# Mask with external mask
        # First, copy over
	fsl:exec "${FSLDIR}/bin/avwmaths $fmri(alternative_mask) mask"
        # Now do the masking
	fsl:exec "avwmaths++ $funcdata $fmri(alternative_mask) prefiltered_func_data_bet"

    } else {
	set ps "$ps; non-brain removal using BET \[Smith 2002\]"
	set rs "$rs\[<a href=\"http://www.fmrib.ox.ac.uk/analysis/techrep/#TR00SMS2\">Smith 2002</a>\] S. Smith. Fast Robust Automated Brain Extraction. Human Brain Mapping 17:3(143-155) 2002.<br>
    "
	fsl:exec "${FSLDIR}/bin/betfunc $funcdata prefiltered_func_data_bet"
    }
    set funcdata prefiltered_func_data_bet
}

#}}}
	#{{{ filtering

set DATATYPE ""
if { [ exec sh -c "$FSLDIR/bin/avwval $funcdata datatype" ] != 4 } {
    set DATATYPE "_32R"
}

set thecommand "${FSLDIR}/bin/ip$DATATYPE ${funcdata} filtered_func_data $fmri(brain_thresh)"

#{{{ Spatial filtering

if { $fmri(smooth) > 0.01 } {
    set thecommand "$thecommand -s [ expr $fmri(smooth) / 2.355 ]"
    set ps "$ps; spatial smoothing using a Gaussian kernel of FWHM $fmri(smooth)mm"
}

#}}}
#{{{ Masking and global intensity normalization

set thecommand "$thecommand -m mask"

set global_mean 10000

if { $fmri(brain_thresh) == 0 } {
    fsl:echo $logout "
Warning - using low global_mean (1000) because thresholding set to 0.
This must be also the case for all other subjects in a group analysis.
"
    set global_mean 1000
}

if { $fmri(norm_yn)} {
    set thecommand "$thecommand -i $global_mean"
    set ps "$ps; global (volumetric) multiplicative mean intensity renormalization"
} else {
    set thecommand "$thecommand -I $global_mean"
    set ps "$ps; mean-based intensity normalisation of all volumes by the same factor"
}

#}}}
#{{{ Perfusion subtraction

if { $fmri(perfsub_yn) } {

    set ps "$ps; control-tag perfusion subtraction with half-TR sinc interpolation and no decimation"

    set absbrainthresh [ fsl:exec "$thecommand | grep THRESHOLD | awk '{ print \$2 }' -" ]

    set perfinvert ""
    if { ! $fmri(tagfirst) } {
	set perfinvert "-c"
    }

    fsl:exec "${FSLDIR}/bin/perfusion_subtract filtered_func_data filtered_func_data $perfinvert"

    set thecommand "${FSLDIR}/bin/ip$DATATYPE filtered_func_data filtered_func_data $fmri(brain_thresh)"
}

#}}}
#{{{ Temporal filtering

if { $fmri(temphp_yn) || $fmri(templp_yn) } {

    set hp_sigma_vol -1
    if { $fmri(temphp_yn) } {
	set hp_sigma_sec [ expr $fmri(paradigm_hp) / 2.0 ]
	set hp_sigma_vol [ expr $hp_sigma_sec / $fmri(tr) ]
	set ps "$ps; highpass temporal filtering (Gaussian-weighted least-squares straight line fitting, with sigma=${hp_sigma_sec}s)"
    }

    set lp_sigma_vol -1
    if { $fmri(templp_yn) } {
	set lp_sigma_sec 2.8
	set lp_sigma_vol [ expr $lp_sigma_sec / $fmri(tr) ]
	set ps "$ps; Gaussian lowpass temporal filtering HWHM ${lp_sigma_sec}s"
    }

    set thecommand "$thecommand -t $hp_sigma_vol $lp_sigma_vol"

}

#}}}

set ps "$ps."

set absbrainthresh [ fsl:exec "$thecommand | grep THRESHOLD | awk '{ print \$2 }' -" ]

#}}}
	#{{{ MELODIC

if { $fmri(melodic_yn) } {
    set ps "$ps ICA-based exploratory data analysis was carried out using MELODIC \[Beckmann 2004\], in order to investigate the possible presence of unexpected artefacts or activation."
    set rs "$rs\[<a href=\"http://www.fmrib.ox.ac.uk/analysis/techrep/#TR02CB1\">Beckmann 2004</a>\] C.F. Beckmann and S.M. Smith. Probabilistic Independent Component Analysis for Functional Magnetic Resonance Imaging. IEEE Trans. on Medical Imaging 23:2(137-152) 2004.<br>
    "

    puts $report "<hr><a href=\"filtered_func_data.ica/report/00index.html\">MELODIC data exploration report</a>"

    fsl:exec "${FSLDIR}/bin/melodic -i filtered_func_data -o filtered_func_data.ica -v --nobet --bgthreshold=1 --tr=$fmri(tr) -d 0 --mmthresh=\"0.5\" --report"
}

#}}}

    } else {
        fsl:exec "${FSLDIR}/bin/avwmaths $funcdata filtered_func_data"
    }

    #{{{ set TR in header of filtered_func_data if not already correct

set IMTR [ exec sh -c "$FSLDIR/bin/avwval filtered_func_data pixdim4" ]

if { [ expr abs($IMTR - $fmri(tr)) ] > 0.01 } {
    fsl:exec "${FSLDIR}/bin/avwhd -x filtered_func_data | sed 's/  dt = .*/  dt = '$fmri(tr)'/g' | ${FSLDIR}/bin/avwcreatehd - filtered_func_data"
}

#}}}
    #{{{ make mean_filtered_func_data

fsl:exec "${FSLDIR}/bin/avwmaths filtered_func_data -Tmean mean_func"

#}}}
    #{{{ make binary mask if it hasn't already been made

if { ! [ imtest mask ] } {
    set absbrainthresh [ fsl:exec "${FSLDIR}/bin/ip filtered_func_data filtered_func_data $fmri(brain_thresh) -m mask | grep THRESHOLD | awk '{ print \$2 }' -" ]
}

# Overwrite mask with external mask
if { [ info exists fmri(alternative_mask) ] && [ imtest $fmri(alternative_mask) ] } {
    puts "-------Using alternative mask---------"
    fsl:exec "${FSLDIR}/bin/avwmaths $fmri(alternative_mask) mask"
}

#}}}
    #{{{ set absbrainthresh if brain_thresh==0

if { $fmri(brain_thresh) == 0 } {
    set absbrainthresh [ fsl:exec "$FSLDIR/bin/avwstats filtered_func_data -R | awk '{ print \$1 }' -" ]
} 

#}}}

    fsl:exec "rm -rf prefiltered_func_data*"
    puts $report "<!--prestatsstop-->"
}

#}}}
		} else {
		    #{{{ make copies/moves of higher-level files

if { $fmri(analysis) != 4 } {
    imcp ../bg_image example_func
    imcp ../mean_func mean_func
    imcp ../mask mask
    immv ../[ file tail $funcdata ] filtered_func_data
    immv ../var[ file tail $funcdata ] var_filtered_func_data
    immv ../tdof_t[ string trimleft [ file tail $funcdata ] "cope" ] tdof_filtered_func_data
}

set absbrainthresh -1e10

#}}}
		}
		
		#{{{ stats

if { $fmri(stats_yn) || $fmri(analysis) == 4 } {

    if { $fmri(analysis) != 4 } {

	if { $fmri(level) == 1 } {

	    #{{{ FILM

#{{{ copy input timing files into feat directory

for { set evs 1 } { $evs <= $fmri(evs_orig) } { incr evs 1 } {

    if { $fmri(shape${evs}) == 2 || $fmri(shape${evs}) == 3 } {

	if { ! [ file exists custom_timing_files ] } {
	    fsl:exec "mkdir custom_timing_files"
	}

	fsl:exec "cp $fmri(custom${evs}) custom_timing_files/ev${evs}.txt"
    }
}

#}}}

set film_opts "-sa -ms 5 -sp $FSLDIR/bin/susan_smooth"
set film_text " with local autocorrelation correction"

if { $fmri(motionevs) > 0 && [ file exists mc/prefiltered_func_data_mcf.par ] } {
    catch { fsl:exec "$FSLDIR/bin/feat_model design mc/prefiltered_func_data_mcf.par" } ErrMsg
}

if { ! $fmri(prewhiten_yn) } {
    set film_opts "-noest" 
    set film_text ""
}

fsl:exec "$FSLDIR/bin/film_gls -rn stats $film_opts filtered_func_data design.mat $absbrainthresh" -r

if { ! [ imtest stats/pe1 ] } {
    fsl:echo report.log "Error: FILM did not complete - it probably ran out of memory"
    fsl:echo "" "Error: FILM did not complete - it probably ran out of memory"
    exit 1
}


set ps "$ps Time-series statistical analysis was carried out using FILM $film_text \[Woolrich 2001\]."
set rs "$rs\[<a href=\"http://www.fmrib.ox.ac.uk/analysis/techrep/#TR01MW1\">Woolrich 2001</a>\] M.W. Woolrich, B.D. Ripley, J.M. Brady and S.M. Smith. Temporal Autocorrelation in Univariate Linear Modelling of FMRI Data. NeuroImage 14:6(1370-1386) 2001.<br>
"

#{{{ correct stats for reduced DOF and altered whitened DM if perfusion subtraction was run

if { $fmri(perfsub_yn) } {
    set THEDOF [ expr ( $fmri(npts) / 2 ) - $fmri(evs_real) ]
    if { $THEDOF < 1 } {
	set THEDOF 1
    }
    fsl:exec "echo $THEDOF > stats/dof"

    fsl:exec "${FSLDIR}/bin/avwmaths stats/corrections -mul 2 stats/corrections"
}

#}}}

#}}}

	} else {

	    #{{{ FLAME

set DOFS ""
if { [ imtest tdof_filtered_func_data ] } {
    if { [ exec sh -c "${FSLDIR}/bin/avwnvols tdof_filtered_func_data 2> /dev/null" ] == $fmri(npts) } {
	set DOFS "--dvc=tdof_filtered_func_data"
    }
}

set FTESTS ""
if { [ file exists design.fts ] } {
    set FTESTS "--fc=design.fts"
}

set FLAME ""

set ps "$ps Higher-level analysis was carried out using"

if { $fmri(mixed_yn) == 3 } {
    set FLAME "--fe"
    set ps "$ps a fixed effects model, by forcing the random effects variance to zero in FLAME (FMRIB's Local Analysis of Mixed Effects) \[Beckmann 2003, Woolrich 2004\]."
} else {

    if { $fmri(mixed_yn) == 0 } {
	set FLAME "--jols"
	set ps "$ps OLS (ordinary least squares) simple mixed effects."
    } else {

	if { $fmri(mixed_yn) == 1 } {
	    if { $fmri(thresh) == 3 } {
		set zlt [ expr $fmri(z_thresh) - 0.05 ]
		set zut [ expr $fmri(z_thresh) + 0.35 ]
	    } else {
		set zlt 2
		set zut 20
	    }
	    set stageoneonly ""
	} else {
	    set zlt 100000
	    set zut 100000
	    set stageoneonly "stage 1 only (i.e., without the final MCMC-based stage) "
	}

	set FLAME "--ols --nj=10000 --bi=500 --se=1 --fm --zlt=$zlt --zut=$zut"

	set ps "$ps FLAME (FMRIB's Local Analysis of Mixed Effects) ${stageoneonly}\[Beckmann 2003, Woolrich 2004\]."
    }
}

if { $fmri(mixed_yn) != 0 } {
    set rs "$rs\[<a href=\"http://www.fmrib.ox.ac.uk/analysis/techrep/#TR01CB1\">Beckmann 2003</a>\] C. Beckmann, M. Jenkinson and S.M. Smith. General multi-level linear modelling for group analysis in FMRI. NeuroImage 20(1052-1063) 2003.<br>
               \[<a href=\"http://www.fmrib.ox.ac.uk/analysis/techrep/#TR03MW1\">Woolrich 2004</a>\] M.W. Woolrich, T.E.J Behrens, C.F. Beckmann, M. Jenkinson and S.M. Smith. Multi-level linear modelling for FMRI group analysis using Bayesian inference. NeuroImage 21:4(1732-1747) 2004<br>
               "
}

set NumPoints [ exec sh -c "grep NumPoints design.mat | awk '{ print \$2 }'" ]
set NumWaves  [ exec sh -c "grep NumWaves  design.mat | awk '{ print \$2 }'" ]

# if we have a LOT of subjects then run FLAME a slice at a time
if { $NumPoints < 40 } {

    fsl:exec "$FSLDIR/bin/flame --cope=filtered_func_data --vc=var_filtered_func_data $DOFS --mask=mask --ld=stats --dm=design.mat --cs=design.grp --tc=design.con $FTESTS $FLAME" -r

} else {

    set DIMX [ exec sh -c "$FSLDIR/bin/avwval filtered_func_data dim1" ]
    set DIMY [ exec sh -c "$FSLDIR/bin/avwval filtered_func_data dim2" ]
    set DIMZ [ exec sh -c "$FSLDIR/bin/avwval filtered_func_data dim3" ]
    for { set slice 0 } { $slice < $DIMZ } { incr slice 1 } {
	fsl:exec "$FSLDIR/bin/avwroi mask tmpmask 0 $DIMX 0 $DIMY $slice 1"
	fsl:exec "$FSLDIR/bin/avwroi filtered_func_data tmpcope 0 $DIMX 0 $DIMY $slice 1"
	fsl:exec "$FSLDIR/bin/avwroi var_filtered_func_data tmpvarcope 0 $DIMX 0 $DIMY $slice 1"
	if { [ imtest tdof_filtered_func_data ] } {
	    fsl:exec "$FSLDIR/bin/avwroi tdof_filtered_func_data tmptdof 0 $DIMX 0 $DIMY $slice 1"
	    set DOFS "--dvc=tmptdof"
	}
	fsl:exec "$FSLDIR/bin/flame --cope=tmpcope --vc=tmpvarcope $DOFS --mask=tmpmask --ld=stats$slice --dm=design.mat --cs=design.grp --tc=design.con $FTESTS $FLAME" -r
    }
    foreach f [ imglob -oneperimage stats0/* ] {
	set froot [ file tail $f ]
	fsl:exec "$FSLDIR/bin/avwmerge -z stats0/$froot [ lsort -dictionary [ imglob -oneperimage stats*/$froot ] ]"
    }
    fsl:exec "mv stats0 stats"
    fsl:exec "/bin/rm -rf stats?* tmp*"

}

fsl:exec "echo [ expr $NumPoints - $NumWaves ] > stats/dof"

fsl:exec "/bin/rm -f stats/zem* stats/zols* stats/mask*"
#fsl:exec "gzip -9 stats/mean* stats/fdof* stats/tstat* stats/fstat*"

#}}}

	}
	
	if { [ imtest stats/res4d ] } {
	    fsl:exec "$FSLDIR/bin/smoothest -d [ exec sh -c "cat stats/dof" ] -m mask -r stats/res4d > stats/smoothness"
	    if { $fmri(cleanup_residuals_yn) } {
		fsl:exec "rm -f stats/res4d*"
	    }
	}

    }

    if { $fmri(level) == 1 } {
	set FTESTS ""
	if { [ file exists design.fts ] } {
	    set FTESTS "-f design.fts"
	}
	fsl:exec "$FSLDIR/bin/contrast_mgr $FTESTS stats design.con"
    }

}

#}}}
		#{{{ close report

if { $OLDFEAT || $fmri(analysis)!=4 } {

puts $report "<!--poststatspicsstart-->
<!--poststatspicsstop-->
<!--poststatstsplotstart-->
<!--poststatstsplotstop-->
<!--regstart-->
<!--regstop-->
$ps
<!--poststatspsstart-->
<!--poststatspsstop-->
<!--regpsstart-->
<!--regpsstop-->
$rs
<!--poststatsrsstart-->
<!--poststatsrsstop-->
<!--regrsstart-->
<!--regrsstop-->"

if { [ file exists design.gif ] } {
    puts $report "<p><p><b>Design matrix</b><br><a href=\"design.mat\"><IMG BORDER=0 SRC=\"design.gif\"></a>"

    if { [ file exists design_cov.gif ] } {
	puts $report "<p><p><b>Covariance matrix</b><br><IMG BORDER=0 SRC=\"design_cov.gif\">"
    }
}

puts $report "<HR><FONT SIZE=1>This page produced automatically by FEAT $fmri(version) - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL</A>.</FONT></BODY></HTML>"
close $report
}

#}}}
		#{{{ poststats
		
if { $fmri(poststats_yn) } {
    feat5:poststats 0 [ expr $fmri(level) - 1 ]
}

#}}}
		#{{{ time series plots and report output for each contrast

if { $fmri(stats_yn) || $fmri(poststats_yn) } {

    set fmrifile filtered_func_data
    if { ! [ imtest $fmrifile ] } {
	set fmrifile [ string trimright [ file root ${FD} ] + ]
    }

    fsl:exec "mkdir tsplot"
    fsl:exec "${FSLDIR}/bin/tsplot . -f $fmrifile -o tsplot"

    catch { exec sh -c "cat tsplot/tsplot_index" } errmsg
    regsub -all "tsplot" $errmsg "tsplot/tsplot" errmsg

    feat5:report_insert poststatstsplot "<hr><b>Time series plots</b><p>
$errmsg"

}

#}}}
	    }

	    if { $fmri(level) == 1 } {
		#{{{ registration and re-thresholding

if { $fmri(reginitial_highres_yn) || $fmri(reghighres_yn) || $fmri(regstandard_yn) } {
    #{{{ setup varables and web page report

set existing_mats 0

if { [ file exists reg ] } {
    fsl:exec "mv reg $BACKUPS"
}
fsl:exec "/bin/mkdir reg"
cd reg

imcp ../example_func example_func

# test for weighting image from fieldmap unwarping to use with example_func
set ef_weighting_flag ""
if { [ imtest ../unwarp/EF_UD_fmap_sigloss ] } {
    set ef_weighting_flag "-inweight ../unwarp/EF_UD_fmap_sigloss"
}

# test for pre-unwarping example_func image (in order to create unwarping evaluation images)
set doefd 0
if { [ imtest ../example_func_orig_distorted ] } {
    set doefd 1
    imcp ../example_func_orig_distorted example_func_orig_distorted 
}

set reg_report [ open index.html "w" ]
fconfigure $reg_report -buffering none

puts $reg_report "<HTML>

<TITLE>FEAT Registration Report</TITLE>

<BODY BACKGROUND=\"file:${FSLSLASH}${FSLDIR}/doc/images/fsl-bg.jpg\">

<hr><CENTER>
<H1>FEAT Registration Report</H1>
for FEAT session &nbsp; &nbsp; <a href=\"../report.html\">${FD}</a><br>
[ exec date ]
</CENTER>
"

#}}}
    #{{{ setup initial transforms

set init_initial_highres ""
if { [ file exists $fmri(init_initial_highres) ] } {
    set init_initial_highres "-init $fmri(init_initial_highres)"
}

set init_highres ""
if { [ file exists $fmri(init_highres) ] } {
    set init_highres "-init $fmri(init_highres)"
}

set init_standard ""
if { [ file exists $fmri(init_standard) ] } {
    set init_standard "-init $fmri(init_standard)"
}

#}}}
    #{{{ setup flirt files

if { $fmri(reginitial_highres_yn) } {
    if { [ info exists initial_highres_files($session) ] } {
	fsl:exec "${FSLDIR}/bin/avwmaths [ remove_ext $initial_highres_files($session) ] initial_highres"
    } else { 
	if { ! [ imtest initial_highres ] } {
	    fsl:echo $logout "Warning - registration to initial_highres turned on but
no initial_highres image specified in setup file or in
FEAT directory! Will not register to initial_highres."
	    set fmri(reginitial_highres_yn) 0
	}
    }
}

if { $fmri(reghighres_yn) } {
    if { [ info exists highres_files($session) ] } {
	fsl:exec "${FSLDIR}/bin/avwmaths [ remove_ext $highres_files($session) ] highres"
    } else { 
	if { ! [ imtest highres ] } {
	    fsl:echo $logout "Warning - registration to highres turned on but
no highres image specified in setup file or in
FEAT directory! Will not register to highres."
	    set fmri(reghighres_yn) 0
	}
    }
}

if { $fmri(regstandard_yn) } {
    imcp [ remove_ext $fmri(regstandard) ] standard
}

#}}}
    #{{{ -> highres

if { $fmri(reghighres_yn) } {

    if { $fmri(reginitial_highres_yn) } {

	feat5:flirt example_func initial_highres $fmri(reginitial_highres_dof) $fmri(reginitial_highres_search) trilinear $existing_mats $reg_report $init_initial_highres $ef_weighting_flag
	if { $doefd} {
	    feat5:flirt example_func_orig_distorted initial_highres $fmri(reginitial_highres_dof) $fmri(reginitial_highres_search) trilinear $existing_mats "" $init_initial_highres ""
	}

	feat5:flirt initial_highres highres $fmri(reghighres_dof) $fmri(reghighres_search) trilinear $existing_mats $reg_report $init_highres ""

	fsl:exec "${FSLDIR}/bin/convert_xfm -omat example_func2highres.mat -concat initial_highres2highres.mat example_func2initial_highres.mat"

        feat5:flirt example_func highres 0 0 trilinear 1 $reg_report "" ""
	if { $doefd } {
	    feat5:flirt example_func_orig_distorted highres 0 0 trilinear 1 $reg_report "" ""
	}

    } else {

	feat5:flirt example_func highres $fmri(reghighres_dof) $fmri(reghighres_search) trilinear $existing_mats $reg_report $init_highres $ef_weighting_flag
	if { $doefd} {
	    feat5:flirt example_func_orig_distorted highres $fmri(reghighres_dof) $fmri(reghighres_search) trilinear $existing_mats $reg_report $init_highres ""
	}

    }
}

#}}}
    #{{{ -> standard

if { $fmri(regstandard_yn) } {

    if { $fmri(reghighres_yn) } {

	feat5:flirt highres standard $fmri(regstandard_dof) $fmri(regstandard_search) trilinear $existing_mats $reg_report $init_standard ""

	fsl:exec "${FSLDIR}/bin/convert_xfm -omat example_func2standard.mat -concat highres2standard.mat example_func2highres.mat"

        feat5:flirt example_func standard 0 0 trilinear 1 $reg_report "" ""
	if { $doefd} {
	    feat5:flirt example_func_orig_distorted standard 0 0 trilinear 1 "" "" ""
	}

    } else {

	feat5:flirt example_func standard $fmri(regstandard_dof) $fmri(regstandard_search) trilinear $existing_mats $reg_report $init_standard $ef_weighting_flag
	if { $doefd} {
	    feat5:flirt example_func_orig_distorted standard $fmri(regstandard_dof) $fmri(regstandard_search) trilinear $existing_mats "" $init_standard ""
	}

    }

    # prepare unwarping evaluation short summary image (example_func vs highres)
    if { $doefd && [ imtest example_func2highres ] } {
	fsl:echo .coord "51 57 40"
	set h [ fsl:exec "${FSLDIR}/bin/img2imgcoord -src standard -dest highres -xfm standard2highres.mat .coord | tail -1" ]
	fsl:exec "${FSLDIR}/bin/slicer example_func2highres highres -s 3 -x -[ expr round([ lindex $h 0 ]) ] sla -y -[ expr round([ lindex $h 1 ]) ] slb -z -[ expr round([ lindex $h 2 ]) ] slc"
	fsl:exec "${FSLDIR}/bin/slicer highres example_func2highres -s 3 -x -[ expr round([ lindex $h 0 ]) ] sld -y -[ expr round([ lindex $h 1 ]) ] sle -z -[ expr round([ lindex $h 2 ]) ] slf"
	fsl:exec "$FSLCONVERT -colors 100 +append sla sld slb sle slc slf example_func2highres3sl.gif"
    }

}

#}}}
    #{{{ timing and final outputs

puts $reg_report "<HR><FONT SIZE=1>This page produced automatically by FEAT - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL</A>.</FONT>

</BODY></HTML>
"

fsl:echo $logout "Finished registration at [ exec date ]"

close $reg_report

#}}}
    #{{{ put biblio stuff & unwarp pic & reg link etc. into original report

cd $FD

if { [ file exists reg/example_func2highres3sl.gif ] } {
    feat5:report_insert unwarphr "<br>Unwarped example_func vs. highres for evaluation of unwarping<br>
<a href=\"reg/index.html\"><IMG BORDER=0 SRC=\"reg/example_func2highres3sl.gif\" WIDTH=1000></a><br>"
}

feat5:report_insert_start reg
puts $report "<hr><a href=\"reg/index.html\">Registration report</a>"
if { [ file exists ${FD}/reg/example_func2standard1.gif ] } {
    puts $report "<br><br><a href=\"reg/index.html\"><IMG BORDER=0 SRC=\"reg/example_func2standard1.gif\" WIDTH=1000></a><br>"
}
feat5:report_insert_stop reg

feat5:report_insert regps "Registration to high resolution and/or standard images was carried out using FLIRT \[Jenkinson 2001, 2002\]."

# first blank regrs in report.html, to allow to check for pre-existence of Jenkinson 2002 ref.
feat5:report_insert regrs ""
set reg_rs "\[<a href=\"http://www.fmrib.ox.ac.uk/analysis/techrep/#TR00MJ2\">Jenkinson 2001</a>\] M. Jenkinson and S.M. Smith. A Global Optimisation Method for Robust Affine Registration of Brain Images. Medical Image Analysis 5:2(143-156) 2001.<br>"
if { [ catch { exec sh -c "grep \"Jenkinson 2002\" report.html" } errmsg ] != 0 } {
    set reg_rs "$reg_rs
\[<a href=\"http://www.fmrib.ox.ac.uk/analysis/techrep/#TR02MJ1\">Jenkinson 2002</a>\] M. Jenkinson, P. Bannister, M. Brady and S. Smith. Improved Optimisation for the Robust and Accurate Linear Registration and Motion Correction of Brain Images. NeuroImage 17:2(825-841) 2002.<br>"
}
feat5:report_insert regrs "$reg_rs"

#}}}
}

if { $fmri(poststats_yn) && [ file exists ${FD}/reg/example_func2standard.mat ] && [ file exists ${FD}/stats ] } {
    feat5:poststats 1 0
}

#}}}
	    }
    
	    #{{{ final logout report and finish up

if { $fmri(level) > 1 } {
    set conname ""
    if { $contrastname($session) != "" } {
	set conname "($contrastname($session))"
    }

    puts $greport "<a href=\"[ file tail $FD ]/report.html\">Lower-level contrast $session $conname</a><br>
"
}

fsl:echo $logout "
Finished FEAT at [ exec date ]
To view the FEAT report point your web browser at
${FD}/report.html
"

after 2000
catch { exec kill -9 $tailpid } errmsg

#}}}
	}
    }

    #{{{ final level 2 outputs

if { $fmri(level) > 1 } {

    if { [ file exists design.mat ] } {
	puts $greport "<hr><p><p><b>Design matrix</b><br><a href=\"design.mat\"><IMG BORDER=0 SRC=\"design.gif\"></a>
"
    }

    puts $greport "<HR><FONT SIZE=1>This page produced automatically by FEAT $fmri(version) - a part of <A HREF=\"http://www.fmrib.ox.ac.uk/fsl\">FSL</A>.</FONT>

</BODY></HTML>
"

    close $greport

    fsl:echo "" "Finished higher-level FEAT at [ exec date ]
To view the FEAT report point your web browser at
${gFD}/report.html
"

}

#}}}
}

#}}}

#{{{ call GUI 

if { ! [ info exists INGUI ] } {
    wm withdraw .
    feat5 .r
    tkwait window .r
}

#}}}
