/**
 * @brief Places surface based on an intensity input image. This is meant to provide
 * replacement functionality for mris_make_surfaces in a form that is easier to 
 * maintain.
 */
/*
 * Original Author: Douglas N Greve (but basically a rewrite of mris_make_surfaces by BF)
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

/*
This version of mris_place_surface yields identical output as
mris_make_surfaces under these command-line conditions for simple T1w
input for both cross and long.  The creation of the white and pial
surfaces are done separately using mris_make_surfaces below whereas in
recon-all they are done with one command line. The white surfaces will
be the same but the pial surfaces will be slightly different.
mris_make_surfaces is too complicated to figure out
why. mris_make_surfaces pial calculation simplifies a little when the
white is not also computed, and I was able to match it exactly with
mris_place_surface. This version of mris_place_surface serves as a
documentation for mris_make_surfaces and reference for changes to
mris_place_surface.

NOTE: as of 12/3/2019, this command lines probably do not work anymore

Cross-sectional analysis

setenv SUBJECTS_DIR /autofs/cluster/fsm/users/greve/subjects/trt.fsm030
set subject = dev.xli.fsm030.01
cd $SUBJECTS_DIR/$subject/surf

mris_make_surfaces -c -cortex 0 -output .mms -aseg ../mri/aseg.presurf -white white.preaparc -noaparc -whiteonly -mgz -T1 brain.finalsurfs dev.xli.fsm030.01 lh
mris_make_surfaces -c -cortex 0 -output .mms -aseg ../mri/aseg.presurf -orig_white white.preaparc -white white -whiteonly -mgz -T1 brain.finalsurfs dev.xli.fsm030.01 lh
mris_make_surfaces -nowhite -c -cortex 0 -output .mms -orig_white white -orig_pial white.preaparc -aseg ../mri/aseg.presurf -mgz -T1 brain.finalsurfs dev.xli.fsm030.01 lh

mris_make_surfaces -c -cortex 0 -output .mms -aseg ../mri/aseg.presurf -white white.preaparc -noaparc -whiteonly -mgz -T1 brain.finalsurfs dev.xli.fsm030.01 rh
mris_make_surfaces -c -cortex 0 -output .mms -aseg ../mri/aseg.presurf -orig_white white.preaparc -white white -whiteonly -mgz -T1 brain.finalsurfs dev.xli.fsm030.01 rh
mris_make_surfaces -nowhite -c -cortex 0 -output .mms -orig_white white -orig_pial white.preaparc -aseg ../mri/aseg.presurf -mgz -T1 brain.finalsurfs dev.xli.fsm030.01 rh

mris_place_surface --s dev.xli.fsm030.01 lh orig place.white.preaparc --white
mris_place_surface --s dev.xli.fsm030.01 lh white.preaparc place.white --white
mris_place_surface --s dev.xli.fsm030.01 lh white place.pial --pial --init lh.white.preaparc

mris_place_surface --s dev.xli.fsm030.01 rh orig place.white.preaparc --white
mris_place_surface --s dev.xli.fsm030.01 rh white.preaparc place.white --white
mris_place_surface --s dev.xli.fsm030.01 rh white place.pial --pial --init rh.white.preaparc

Longitudinal analysis - note that there is a bug in recon-all that causes mris_make_surface to not use the white_preaparc
when creating the white and pial surfaces

setenv SUBJECTS_DIR /autofs/cluster/fsm/users/greve/subjects/trt.fsm030
set subject = fsm030.01.long.base.fsm030
cd $SUBJECTS_DIR/$subject/surf

mris_make_surfaces -c -cortex 0 -output .mms -orig_white orig_white -orig orig_white -long -max 3.5 -aseg ../mri/aseg.presurf -white white.preaparc -whiteonly -mgz -noaparc -T1 brain.finalsurfs fsm030.01.long.base.fsm030 lh
mris_make_surfaces -c -cortex 0 -output .mms -orig_white orig_white -orig orig_white -long -max 3.5 -aseg ../mri/aseg.presurf -white white -whiteonly -mgz -T1 brain.finalsurfs fsm030.01.long.base.fsm030 lh
mris_make_surfaces -nowhite -c -cortex 0 -output .mms -orig orig_white -orig_white orig_white -orig_pial orig_pial -long -max 3.5 -aseg ../mri/aseg.presurf -mgz -T1 brain.finalsurfs fsm030.01.long.base.fsm030 lh

mris_make_surfaces -c -cortex 0 -output .mms -orig_white orig_white -orig orig_white -long -max 3.5 -aseg ../mri/aseg.presurf -white white.preaparc -whiteonly -mgz -noaparc -T1 brain.finalsurfs fsm030.01.long.base.fsm030 rh
mris_make_surfaces -c -cortex 0 -output .mms -orig_white orig_white -orig orig_white -long -max 3.5 -aseg ../mri/aseg.presurf -white white -whiteonly -mgz -T1 brain.finalsurfs fsm030.01.long.base.fsm030 rh
mris_make_surfaces -nowhite -c -cortex 0 -output .mms -orig orig_white -orig_white orig_white -orig_pial orig_pial -long -max 3.5 -aseg ../mri/aseg.presurf -mgz -T1 brain.finalsurfs fsm030.01.long.base.fsm030 rh

mris_place_surface --s fsm030.01.long.base.fsm030 lh orig_white white.preaparc.mps --white --long --max-thickness 3.5
mris_place_surface --use-aparc --s fsm030.01.long.base.fsm030 lh orig_white white.mps --white --long --max-thickness 3.5
mris_place_surface --use-aparc --long --max-thickness 3.5 --s fsm030.01.long.base.fsm030 lh orig_white pial.mps --pial --init lh.orig_pial

mris_place_surface --s fsm030.01.long.base.fsm030 rh orig_white white.preaparc.mps --white --long --max-thickness 3.5
mris_place_surface --use-aparc --s fsm030.01.long.base.fsm030 rh orig_white white.mps --white --long --max-thickness 3.5
mris_place_surface --use-aparc --long --max-thickness 3.5 --s fsm030.01.long.base.fsm030 rh orig_white pial.mps --pial --init rh.orig_pial

# To run with T2 (or FLAIR)
# The -nsigma_zzz does not appear to have an effect
mris_make_surfaces -c -cortex 0 -output .redo -orig_white white -orig_pial woT2.pial \
  -aseg ../mri/aseg.presurf -nowhite -mgz -T1 brain.finalsurfs -T2 ../mri/T2 \
  -nsigma_above 2 -nsigma_below 5 dev.hrdecremdn.996782 rh
# This does not give an exact replication, but it is close
mris_place_surface --adgws-in autodetstats.lh.dat 
  --seg ../mri/aseg.presurf.mgz --wm ../mri/wm.mgz --invol ../mri/brain.finalsurfs.mgz 
  --lh --i lh.woT2.pial --init lh.woT2.pial --o lh.T2.pial 
  --pial --nsmooth 0 --rip-label ../label/lh.cortex.label
  --aparc ../label/lh.aparc.annot --repulse-surf lh.white 
  --mmvol--mmvol ../mri/T2.mgz T2 --white-surf lh.white
# The role of the input (--i) vs init (--init) surface needs to be clarified

# To run with just an input surface and input volume
 mris_place_surface --adgws-in autodetstats.lh.dat --invol ../mri/brain.finalsurfs.mgz --lh --i ../surf/lh.orig --white --o lh.junk --no-intensity-proc --no-rip-midline

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <float.h>
#include <errno.h>

#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "surfgrad.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "label.h"
#include "annotation.h"
#include "cmdargs.h"
#include "cma.h"
#include "romp_support.h"
#include "mris_multimodal_refinement.h"

class RIP_MNGR{
public:
  int RipVertices(void);
  MRIS *surf;
  MRI *seg, *invol;
  const char *hemi;
  int RipFreeze = 1;
  int RipLesion = 0;
  int RipWMSA = 0;
  int nRipSegs= 0;
  int RipSegNo[100];
  int RipBG = 0;
  int RipMidline = 1;
  char *riplabelfile = NULL;
  char *ripsurffile=NULL;
  MRIS *ripsurf;
  char *aparcpath=NULL;
  double dmin = -2.0, dmax = +2.0, dstep = 0.5;
};

int MRISripBasalGanglia(MRIS *surf, MRI *seg, const double dmin, const double dmax, const double dstep);
int MRISripSegs(MRIS *surf, MRI *seg, const double dmin, const double dmax, const double dstep);
int MRISpinMedialWallToWhite(MRIS *surf, const LABEL *cortex);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);

struct utsname uts;
char *cmdline, cwd[2000];
int debug = 0, checkoptsonly = 0;

int main(int argc, char *argv[]) ;

const char *Progname = "mris_place_surfaces";

INTEGRATION_PARMS parms, old_parms ;
int lh_label = LH_LABEL ;
int rh_label = RH_LABEL ;
double mid_gray = 67.5;

int max_pial_averages = 16 ;
int min_pial_averages = 2 ;
int max_white_averages = 4 ;
int min_white_averages = 0 ;
float pial_sigma = 2.0f ;
float white_sigma = 2.0f ;
float max_cbv_dist = 5.0 ; // same as max_thickness in MMS
int vavgs = 5 ;
int nthreads = 1;
int nbrs = 2;
int nsmoothsurf = 0 ;

char *SUBJECTS_DIR;
char *insurfpath = NULL;
char *blendsurfpath = NULL;
double blendweight = 0;
char *outsurfpath = NULL;
char *involpath=NULL;
char *segvolpath=NULL;
char *wmvolpath=NULL;
char *aparcpath=NULL;
char *repulsesurfpath = NULL;
char *whitesurfpath = NULL;

// These are for when using subject 
char *subject = NULL, *insurfname = NULL, *outsurfname = NULL;
const char* hemi = NULL;
const char *involname="brain.finalsurfs.mgz", *segvolname="aseg.presurf.mgz",*wmvolname="wm.mgz",*aparcname="aparc";

char tmpstr[2000];
int err=0;
//int longitudinal = 0; // currently does nothing
int surftype = -1; //GRAY_WHITE; // GRAY_CSF
int UseAParc = 0;
char *adgwsinfile = NULL;
char *adgwsoutfile = NULL;
MRIS *surf;
char *ripflagout = NULL;
RIP_MNGR ripmngr;
LABEL *pinlabel = NULL;
int DoIntensityProc = 1;

double shrinkThresh = -1;

int mm_contrast_type = -1;
char *mmvolpath = NULL;
MRI *mmvol = NULL;
float T2_min_inside = 110 ;
float T2_max_inside = 300 ;
double T2_min_outside = 130;
double T2_max_outside = 300 ;
double max_outward_dist = 1 ; // for normal cases, for pathology might be much larger
double wm_weight = 3 ;    // threshold for finding gm outliers: thresh = (wt*wm + gm )/(wt+1)
double Ghisto_left_inside_peak_pct = 0.1; //0.01 for flair. typo?
double Ghisto_right_inside_peak_pct = 0.01;
double Ghisto_left_outside_peak_pct = 0.5;
double Ghisto_right_outside_peak_pct = 0.5;
int n_averages=0;
int UseMMRefine = 0;
AutoDetGWStats adgws;
char *coversegpath = NULL;
MRI *mri_cover_seg = NULL;

/*--------------------------------------------------*/
int main(int argc, char **argv) 
{
  int nargs, i, msec;
  double        spring_scale = 1;
  MRI *invol, *seg=NULL, *wm, *involCBV, *involPS;
  Timer timer ;
  char *cmdline2, cwd[2000];
  //char *field=NULL;

  nargs = handleVersionOption(argc, argv, "mris_place_surface");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);
  cmdline2 = argv2cmdline(argc,argv);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  Gdiag |= DIAG_SHOW ;

  memset(&parms, 0, sizeof(parms)) ;
  // don't let gradient use exterior information (slows things down)
  parms.fill_interior = 0 ;
  parms.projection = NO_PROJECTION ;
  parms.tol = 1e-4 ;
  parms.dt = 0.5f ;
  parms.base_dt = parms.dt ;

  parms.l_curv = 1.0 ;
  parms.l_intensity = 0.2 ;
  parms.l_tspring = 1.0f ;
  parms.l_nspring = 0.5f ;
  parms.l_spring = 0.0f ;
  parms.l_surf_repulse = 0.0 ;
  parms.l_spring_nzr = 0.0 ;
  parms.l_spring_nzr_len = 0.0 ;
  parms.l_hinge = 0;
  parms.l_tsmooth = 0;

  parms.niterations = 0 ;
  parms.write_iterations = 0 /*WRITE_ITERATIONS */;
  parms.integration_type = INTEGRATE_MOMENTUM ;
  parms.momentum = 0.0 /*0.8*/ ;
  parms.dt_increase = 1.0 /* DT_INCREASE */;
  parms.dt_decrease = 0.50 /* DT_DECREASE*/ ;
  parms.error_ratio = 50.0 /*ERROR_RATIO */;
  if(parms.momentum < 0.0) parms.momentum = 0.0 ;
  parms.niterations = 100;

  if(argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if(checkoptsonly) return(0);

  if(surftype == GRAY_WHITE){
    // White
    parms.l_repulse = 5.0 ;
    parms.l_surf_repulse = 0.0 ;
  }

  if(surftype == GRAY_CSF){
    // Pial
    parms.l_repulse = 0.0;
    parms.l_surf_repulse = 5.0;
  }

  // print out version of this program and mrisurf.c
  printf("%s\n",getVersion().c_str());
  printf("%s\n",getVersion().c_str());
  printf("\n");
  printf("cd %s\n",cwd);
  printf("setenv SUBJECTS_DIR %s\n",getenv("SUBJECTS_DIR"));
  printf("%s\n",cmdline2);
  printf("\n");
  fflush(stdout);

  if(adgwsinfile == NULL){
    // Note: in long stream orig = orig_white
    err = adgws.AutoDetectStats(subject, hemi);
    if(err) exit(1);
  }
  if(adgwsoutfile){
    err = adgws.Write(adgwsoutfile);
    if(err) exit(1);
  }

  printf("Reading in input surface %s\n",insurfpath);
  surf = MRISread(insurfpath);
  if(surf==NULL) exit(1);

  if(blendsurfpath){
    // This is done in long with blendweight=0.25
    // Save input coords into tmp
    MRISsaveVertexPositions(surf, TMP_VERTICES);
    // Read blend surface coords into the current coords
    err = MRISreadVertexPositions(surf, blendsurfpath);
    if(err) exit(1);
    // Average blend with current. 
    MRISblendXYZandTXYZ(surf, blendweight, 1.0-blendweight);
  }

  MRISedges(surf);
  MRIScorners(surf);
  MRIScomputeMetricProperties(surf);
  if(nbrs > 1) MRISsetNeighborhoodSizeAndDist(surf, nbrs) ;
  if(nsmoothsurf > 0) {
    printf("Smoothing surface with %d iterations\n",nsmoothsurf);
    // In mris_make_surface, this is not done when orig_white is specified, ie,
    // it is done when the orig surface is used for initiation (eg, when 
    // creating white.preaparc). Don't smooth for pial.
    MRISaverageVertexPositions(surf, nsmoothsurf) ;
  }
  else printf("Not smoothing input surface\n");

  MRIScomputeMetricProperties(surf);
  MRISstoreMetricProperties(surf) ;
  MRISremoveIntersections(surf); // done in mris_make_surfaces
  MRISsetVals(surf,-1) ;  /* clear white matter intensities */
  MRISsetVal2(surf, 0) ;   // will be marked for vertices near lesions
  MRISclearMark2s(surf) ;  
  MRISfaceMetric(surf,0);
  MRISedgeMetric(surf,0);
  MRIScornerMetric(surf,0);
  MRISprettyPrintSurfQualityStats(stdout, surf);

  if(whitesurfpath){
    // white coords used in t2/flair and to force pial to white in medial wall
    printf("Reading white surface coordinates from %s\n",whitesurfpath);
    MRISsaveVertexPositions(surf, TMP_VERTICES) ; 
    err = MRISreadVertexPositions(surf, whitesurfpath);
    if(err) exit(1);
    MRISsaveVertexPositions(surf, WHITE_VERTICES) ; // This is used for T2 pial
    MRISrestoreVertexPositions(surf, TMP_VERTICES);
  }
  else MRISsaveVertexPositions(surf, WHITE_VERTICES) ;

  if(repulsesurfpath){
    printf("Reading repulsion surface coordinates from %s\n",repulsesurfpath);
    MRISsaveVertexPositions(surf, TMP_VERTICES) ; 
    err = MRISreadVertexPositions(surf, repulsesurfpath);
    if(err) exit(1);
    MRISsaveVertexPositions(surf, ORIGINAL_VERTICES) ; // This is used for repulsion
    MRISrestoreVertexPositions(surf, TMP_VERTICES);
  }
  else MRISsaveVertexPositions(surf, ORIGINAL_VERTICES) ; // This is used for repulsion

  if(surftype == GRAY_CSF){
    MRISsaveVertexPositions(surf, PIAL_VERTICES) ; 
  }

  // This might not do anything
  MRISsaveVertexPositions(surf, TMP_VERTICES) ; 

  if(aparcpath) {
    printf("Reading in aparc %s\n",aparcpath);
    if (MRISreadAnnotation(surf, aparcpath) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read annotation",aparcpath) ;
  }
  else printf("Not reading in aparc\n");

  printf("Reading in input volume %s\n",involpath);
  invol = MRIread(involpath);
  if(invol==NULL) exit(1);

  if(DoIntensityProc){
    // =========== intensity volume preprocessing ==============
    // It would be nice to do this externally
    printf("Reading in wm volume %s\n",wmvolpath);
    wm = MRIread(wmvolpath);
    if(wm==NULL) exit(1);
    // Clip invol invol voxel intensity to 110 (if it is in the wmmask)
    MRIclipBrightWM(invol, wm) ;
    
    MRI *mri_labeled = MRIfindBrightNonWM(invol, wm) ;
    MRIfree(&wm);
    if(surftype == GRAY_WHITE){
      printf("Masking bright non-wm for white surface\n");
      // Replace bright and borderbright invol voxels with 0
      // Can this be done for pial as well?
      MRImask(invol, mri_labeled, invol, BRIGHT_LABEL, 0) ;
      MRImask(invol, mri_labeled, invol, BRIGHT_BORDER_LABEL, 0) ;
    }
    
    if(surftype == GRAY_CSF){
      // Modify the volume for pial surface
      printf("Masking bright non-wm for pial surface mid_gray = %g\n",adgws.MID_GRAY);
      // Replace brightborder voxels with MID_GRAY (can be done once?)
      // Why would you want to do this? The brightborder voxels are >= 100,
      // but this would make them look like GM and force the pial outside of them.
      // This would happen, eg, for a bright vessel in cortex
      MRImask(invol, mri_labeled, invol, BRIGHT_BORDER_LABEL, adgws.MID_GRAY);
      // Replace bright voxels with 255 (this gets changed below)
      // Not sure why this is needed except that it gets set to 0 below which
      // could look like a strong negative gradient
      // Create the volume used to computed the border values
      involCBV = MRImask(invol, mri_labeled, NULL, BRIGHT_LABEL, 255) ;
      /* From mris_make_surfaces: Replace bright stuff such as eye
	 sockets with 255.  Simply zeroing it out would make the border
	 always go through the sockets, and ignore subtle local minima
	 in intensity at the border of the sockets.  Will set to 0
	 after border values have been computed so that it doesn't mess
	 up gradients.  */
      // Replace bright voxels with 0 (mask them out)
      // This undoes some of the masking above (or vice versa)
      // Create the volume used to position the surface
      involPS = MRImask(invol, mri_labeled, NULL, BRIGHT_LABEL, 0);
    }
    else {
      // Use the same input for white surface
      involCBV = invol;
      involPS  = invol;
    }
    MRIfree(&mri_labeled);
    // ========= End intensity volume preproc ===================
  }
  else {
    involCBV = invol;
    involPS  = invol;
  }

  if(segvolpath){
    printf("Reading in seg volume %s\n",segvolpath);
    seg = MRIread(segvolpath);
    if(seg==NULL) exit(1);
    // Check dims
  }

  // Manage ripping of vertices
  ripmngr.hemi  = hemi;
  ripmngr.invol = invol;
  ripmngr.seg   = seg;
  ripmngr.surf  = surf;
  ripmngr.RipVertices();

  if(mmvolpath){
    printf("Reading in multimodal volume %s\n",mmvolpath);
    mmvol = MRIread(mmvolpath);
    if(mmvol==NULL) exit(1);
    // should check that the dimensions, resolution are the same
    parms.l_intensity = 0.000;
    parms.l_location  = 0.500;
    parms.l_repulse   = 0.025;
    parms.l_nspring   = 1.000;
    parms.l_curv      = 0.500;
    parms.check_tol = 1 ;
    wm_weight = 3 ;
    max_pial_averages = 4;
    max_outward_dist = 3 ;
    if(mm_contrast_type == CONTRAST_FLAIR) Ghisto_left_inside_peak_pct = 0.01; //0.01 for flair. typo?
    // Check dims
  }

  if(parms.l_hinge > 0 || parms.l_spring_nzr > 0){
    if(parms.l_spring_nzr){
      double  *edgestats = MRISedgeStats(surf, 0, NULL, NULL);
      parms.l_spring_nzr_len = edgestats[1];
      free(edgestats);
    }
  }

  // Print out to help determine what is what
  VERTEX *vv = &surf->vertices[(int)round(surf->nvertices/2.0)];
  printf("vertex %d: xyz = (%g,%g,%g) oxyz = (%g,%g,%g) wxzy = (%g,%g,%g) pxyz = (%g,%g,%g) \n",
	 (int)round(surf->nvertices/2.0),vv->x,vv->y,vv->z, vv->origx,vv->origy,vv->origz, 
	 vv->whitex,vv->whitey,vv->whitez,vv->pialx,vv->pialy,vv->pialz); fflush(stdout);

  double inside_hi=0, border_hi=0, border_low=0, outside_low=0, outside_hi=0,current_sigma=0;
  int n_min_averages=0;
  if(surftype == GRAY_WHITE){
    current_sigma = white_sigma ;
    if(n_averages == 0) n_averages = max_white_averages;     
    n_min_averages = min_white_averages; 
    inside_hi = adgws.white_inside_hi;
    border_hi = adgws.white_border_hi;
    border_low = adgws.white_border_low;
    outside_low = adgws.white_outside_low;
    outside_hi = adgws.white_outside_hi;
  }
  if(surftype == GRAY_CSF){
    current_sigma = pial_sigma ;
    if(n_averages == 0) n_averages = max_pial_averages;     
    n_min_averages = min_pial_averages; 
    inside_hi = adgws.pial_inside_hi;
    border_hi = adgws.pial_border_hi;
    border_low = adgws.pial_border_low;
    outside_low = adgws.pial_outside_low;
    outside_hi = adgws.pial_outside_hi;
  }

  timer.reset() ;
  printf("n_averages %d\n",n_averages);
  for (i = 0 ;  n_averages >= n_min_averages ; n_averages /= 2, current_sigma /= 2, i++) {

    printf("Iteration %d =========================================\n",i);
    printf("n_averages=%d, current_sigma=%g\n",n_averages,current_sigma); fflush(stdout);

    if(shrinkThresh >= 0 && i==1){
      printf("Shrinking big triangles %g\n",shrinkThresh);
      MRISshrinkFaces(surf, shrinkThresh, 0);
    }

    if(seg && surftype == GRAY_WHITE && ripmngr.RipMidline){
      // This is done each cycle with white.preaparc
      printf("Freezing midline and others\n");  fflush(stdout);
      ripmngr.RipVertices();
    }

    if(mri_cover_seg) {
      if(i == 0){
	MRI *mri_bin, *mri_tmp;
	mri_bin = MRIclone(involCBV, NULL) ;
	mri_tmp = MRIScoverSeg(surf, mri_bin, mri_cover_seg, surftype);
	MRIfree(&mri_bin) ; 
	MRIfree(&involCBV);
	involCBV = mri_tmp ;
      }
    }

    parms.sigma = current_sigma ;
    parms.n_averages = n_averages ;

    INTEGRATION_PARMS_copy(&old_parms, &parms) ;
    if(mmvol == NULL){
      // Compute the target intensity value (l_intensity)
      printf("Computing target border values \n");
      // The outputs are set in each vertex structure:
      //   v->val2 = current_sigma; // smoothing level along gradient used to find the target
      //   v->val  = max_mag_val; // intensity at target location
      //   v->d = max_mag_dist;   // dist to target along normal
      //   v->mean = max_mag;     // derivative at target intensity
      //   v->marked = 1;         // vertex has good data
      //   v->targx = v->x + v->nx * v->d; // same for y and z
      MRIScomputeBorderValues(surf, involCBV, NULL, inside_hi,border_hi,border_low,outside_low,outside_hi,
			      current_sigma, 2*max_cbv_dist, parms.fp, surftype, NULL, 0, parms.flags,seg,-1,-1) ;
      // Note: 3rd input (NULL) was "mri_smooth" in mris_make_surfaces, but
      // this was always a copy of the input (mri_T1 or invol); it is not used in CBV
      
      if(seg && surftype == GRAY_WHITE){
	printf("Finding expansion regions\n"); fflush(stdout);
	// Masks out the v->curv field of vertices with long distances (v->d)
	MRISfindExpansionRegions(surf) ;
      }

      if(vavgs > 0) {
	printf("Averaging target values for %d iterations...\n",vavgs) ;
	// MRIScomputeBorderValues() sets v->marked=1 for all unripped
	MRISaverageMarkedVals(surf, vavgs) ;
      }

      /* BF's note from MMS: There are frequently regions of gray
	whose intensity is fairly flat. We want to make sure the
	surface settles at the innermost edge of this region, so on
	the first pass, set the target intensities artificially high
	so that the surface will move all the way to white matter
	before moving outwards to seek the border (I know it's a hack,
	but it improves the surface in a few areas. The alternative is
	to explicitly put a gradient-seeking term in the cost
	functional instead of just using one to find the target
	intensities). */
      
    } 
    else {
      // Compute the target xyz coordinate (l_location)
      printf("Computing pial target locations using multimodal (%d)\n",mm_contrast_type); fflush(stdout);
      if(!UseMMRefine){
	MRIScomputePialTargetLocationsMultiModal(surf, mmvol, NULL, 0, 
					       mm_contrast_type, seg, 
					       T2_min_inside, T2_max_inside, 
					       T2_min_outside, T2_max_outside, 
					       max_outward_dist,
					       Ghisto_left_inside_peak_pct, Ghisto_right_inside_peak_pct, 
					       Ghisto_left_outside_peak_pct, Ghisto_right_outside_peak_pct, 
					       wm_weight, pial_sigma, invol) ;
      }
      else{
	printf("UseMMRefine\n");
	MRIS_MultimodalRefinement* refine = new MRIS_MultimodalRefinement();
	MRI* whiteMR  = MRIcopy(invol,NULL);
	MRI* vesselMR = MRIcopy(invol,NULL);
	refine->SegmentWM(invol,mmvol, whiteMR);
	refine->SegmentVessel(invol,mmvol, vesselMR);
	refine->SetStep(.4);
	refine->SetNumberOfSteps(12);
	refine->SetGradientSigma(.3);
	refine->SetSegmentation(seg);
	refine->FindMaximumGradient(mm_contrast_type == CONTRAST_T2);
	refine->addImage(invol);
	refine->addImage(mmvol);
	refine->SetWhiteMR(whiteMR);
	refine->SetVesselMR(vesselMR);
	refine->getTarget(surf); //, debugVertex);
	MRIfree(&whiteMR);
	MRIfree(&vesselMR);
	delete refine;
        MRIS* surftarget = MRISclone(surf);
	for(int v =0; v<surf->nvertices;v++)
	{
		surftarget->vertices[v].x =surf->vertices[v].targx;	
		surftarget->vertices[v].y =surf->vertices[v].targy;	
		surftarget->vertices[v].z =surf->vertices[v].targz;	
	}
  	 MRISwrite(surftarget,std::string(std::string(outsurfpath)+std::string(".target") ).c_str());
	MRISfree(&surftarget);
	 }		
    }


    // This appears to adjust the cost weights based on the iteration but in
    // practice, they never change because spring scale is 1
    parms.l_nspring *= spring_scale ; 
    parms.l_spring *= spring_scale ; 
    // This line with tspring being ajusted twice was probably originally a typo
    // but it has existed this way for a long time. It was changed after 
    // version 6, but I just changed it back for consistency. 
    parms.l_tspring *= spring_scale ;  parms.l_tspring *= spring_scale ;

    if(surftype == GRAY_CSF){
      // This had a bad effect on highres white and no effect on 1mm
      // Used on pial in mris_make_surfaces
      parms.l_tspring = MIN(1.0,parms.l_tspring) ;
    }
    parms.l_nspring = MIN(1.0, parms.l_nspring) ;
    parms.l_spring = MIN(1.0, parms.l_spring) ;
    printf("Positioning Surface: tspring = %g, nspring = %g, spring = %g, niters = %d ",
	   parms.l_tspring,parms.l_nspring,parms.l_spring,parms.niterations); 
    printf("l_repulse = %g, l_surf_repulse = %g, checktol = %d\n",parms.l_repulse,
	   parms.l_surf_repulse,parms.check_tol);fflush(stdout);

    printf("Positioning pial surface\n");fflush(stdout);
    // Note: 3rd input (invol) was "mri_smooth" in mris_make_surfaces, but
    // this was always a copy of the input (mri_T1 or invol)
    MRISpositionSurface(surf, involPS, involPS, &parms);

    //sprintf(tmpstr,"lh.place.postpos%02d",i);
    //MRISwrite(surf, tmpstr);

    old_parms.start_t = parms.start_t ;
    INTEGRATION_PARMS_copy(&parms, &old_parms) ;

    // MMS runs	MRISsetCroppedToZero(mris) with T2/FLAIR, but not clear it does anything

    if(!n_averages)
      break ; 

  } // end major loop placing the white surface using

  if(pinlabel){
    printf("Pinning medial wall to white surface\n");
    MRISpinMedialWallToWhite(surf, pinlabel);
  }

  MRISremoveIntersections(surf); //matches mris_make_surface

  printf("\n\n");
  printf("Writing output to %s\n",outsurfpath);
  err = MRISwrite(surf,outsurfpath);
  if(err){
    printf("ERROR: writing to %s\n",outsurfpath);
    exit(1);
  }

  if(ripflagout){
    printf("Writing ripflagout to %s\n",ripflagout);
    const char *field = "ripflag";    
    MRISwriteField(surf, &field, 1, ripflagout);
  }

  msec = timer.milliseconds() ;
  printf("#ET# mris_place_surface %5.2f minutes\n", (float)msec/(60*1000.0f));
  printf("#VMPC# mris_place_surfaces VmPeak  %d\n",GetVmPeak());
  printf("mris_place_surface done\n");

  return(0);

}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if(!strcasecmp(option, "--version")) print_version() ;
    else if(!strcasecmp(option, "--debug"))   debug = 1;
    else if(!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if(!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    //else if(!strcmp(option, "--long")) longitudinal = 1; currently does nothing
    else if(!strcmp(option, "--white")) surftype = GRAY_WHITE;
    else if(!strcmp(option, "--pial"))  surftype = GRAY_CSF;
    else if(!strcmp(option, "--use-aparc")) UseAParc = 1;
    else if(!strcmp(option, "--rip-wmsa"))    ripmngr.RipWMSA = 1;
    else if(!strcmp(option, "--no-rip-wmsa")) ripmngr.RipWMSA = 0;
    else if(!strcmp(option, "--rip-freeze"))    ripmngr.RipFreeze = 1;
    else if(!strcmp(option, "--no-rip-freeze")) ripmngr.RipFreeze = 0;
    else if(!strcmp(option, "--rip-lesion"))    ripmngr.RipLesion = 1;
    else if(!strcmp(option, "--no-rip-lesion"))    ripmngr.RipLesion = 0;
    else if(!strcmp(option, "--rip-bg"))    ripmngr.RipBG = 1;
    else if(!strcmp(option, "--no-rip-bg")) ripmngr.RipBG = 0;
    else if(!strcmp(option, "--rip-midline"))     ripmngr.RipMidline = 1;
    else if(!strcmp(option, "--no-rip-midline"))  ripmngr.RipMidline = 0;
    else if(!strcmp(option, "--no-intensity-proc"))  DoIntensityProc = 0;
    else if(!strcmp(option, "--first-peak-d1"))    CBVfindFirstPeakD1 = 1;
    else if(!strcmp(option, "--no-first-peak-d1")) CBVfindFirstPeakD1 = 0;
    else if(!strcmp(option, "--first-peak-d2"))    CBVfindFirstPeakD2 = 1;
    else if(!strcmp(option, "--no-first-peak-d2")) CBVfindFirstPeakD2 = 0;
    else if(!strcmp(option, "--lh"))  hemi = "lh";
    else if(!strcmp(option, "--rh"))   hemi = "rh";
    else if(!strcmp(option, "--rip-projection")){
      if(nargc < 3) CMDargNErr(option,3);
      sscanf(pargv[0],"%lf",&ripmngr.dmin);
      sscanf(pargv[1],"%lf",&ripmngr.dmax);
      sscanf(pargv[2],"%lf",&ripmngr.dstep);
      nargsused = 3;
    }
    else if(!strcmp(option, "--pin-medial-wall")){
      if(nargc < 1) CMDargNErr(option,1);
      pinlabel = LabelRead("",pargv[0]);
      if(pinlabel == NULL) exit(1);
      nargsused = 1;
    }
    else if(!strcmp(option, "--no-pin-medial-wall"))  pinlabel = NULL;
    else if(!strcasecmp(option, "--seg")){
      if(nargc < 1) CMDargNErr(option,1);
      segvolpath = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--wm")){
      if(nargc < 1) CMDargNErr(option,1);
      wmvolpath = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--invol")){
      if(nargc < 1) CMDargNErr(option,1);
      involpath = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--mmvol")){
      if(nargc < 2) CMDargNErr(option,2);
      mmvolpath = pargv[0];
      if(!stricmp(pargv[1],"t2"))         mm_contrast_type = CONTRAST_T2;
      else if(!stricmp(pargv[1],"flair")) mm_contrast_type = CONTRAST_FLAIR;
      else {
	printf("ERROR: mmvol must be either t2 or flair weighted\n");
	exit(1);
      }
      surftype = GRAY_CSF;
      nargsused = 2;
    } 
    else if(!strcasecmp(option, "--cover-seg") || !strcasecmp(option, "--cover_seg")){
      if(nargc < 1) CMDargNErr(option,1);
      coversegpath = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--mm-refine"))   UseMMRefine = 1;
    else if(!strcmp(option, "--i")){
      if(nargc < 1) CMDargNErr(option,1);
      insurfpath = pargv[0];
      parms.l_tspring = 0.3;
      parms.l_nspring = 0.3;
      nargsused = 1;
    }
    else if(!strcmp(option, "--blend-surf")){
      if(nargc < 2) CMDargNErr(option,2);
      sscanf(pargv[0],"%lf",&blendweight);
      blendsurfpath = pargv[1];
      nargsused = 2;
    }
    else if(!strcasecmp(option, "--rip-surf")){
      if(nargc < 1) CMDargNErr(option,1);
      ripmngr.ripsurffile = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--rip-label")){
      if(nargc < 1) CMDargNErr(option,1);
      ripmngr.riplabelfile = pargv[0];
      ripmngr.RipMidline = 0;
      nargsused = 1;
    } 
    else if(!strcmp(option, "--o")){
      if(nargc < 1) CMDargNErr(option,1);
      outsurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--repulse-surf")){
      if(nargc < 1) CMDargNErr(option,1);
      repulsesurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--white-surf")){
      if(nargc < 1) CMDargNErr(option,1);
      whitesurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--aparc")){
      if(nargc < 1) CMDargNErr(option,1);
      aparcpath = pargv[0];
      ripmngr.aparcpath = aparcpath;
      UseAParc = 1;
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--max-cbv-dist")){
      // Limit on the distance that CBV will search along the normal (inside and out)
      // This is called max_thickness in MMS. Usually set to 3.5mm for long
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%f",&max_cbv_dist);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--adgws") || !strcasecmp(option, "--adgws-in")){
      if(nargc < 1) CMDargNErr(option,1);
      adgwsinfile = pargv[0];
      err = adgws.Read(adgwsinfile);
      if(err){
	printf("ERROR: reading %s\n",adgwsinfile);
	exit(1);
      }
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--adgws-out")){
      if(nargc < 1) CMDargNErr(option,1);
      adgwsoutfile = pargv[0];
      nargsused = 1;
    } 

    // The white/pial border/inside/outside hi/low must be specified
    // AFTER --agws-in
    else if(!strcmp(option, "--white_border_hi")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.white_border_hi = atof(pargv[0]);
      nargsused = 1;
    }
    else if(!strcmp(option, "--white_border_low")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.white_border_low = atof(pargv[0]);
      nargsused = 1;
    }
    else if(!strcmp(option, "--white_outside_low")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.white_outside_low = atof(pargv[0]);
      nargsused = 1;
    }
    else if(!strcmp(option, "--white_inside_hi")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.white_inside_hi = atof(pargv[0]);
      nargsused = 1;
    }
    else if(!strcmp(option, "--white_outside_hi")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.white_outside_hi = atof(pargv[0]);
      nargsused = 1;
    }
    else if(!strcmp(option, "--pial_border_hi")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.pial_border_hi = atof(pargv[0]);
      nargsused = 1;
    }
    else if(!strcmp(option, "--pial_border_low")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.pial_border_low = atof(pargv[0]);
      nargsused = 1;
    }
    else if(!strcmp(option, "--pial_outside_low")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.pial_outside_low = atof(pargv[0]);
      nargsused = 1;
    }
    else if(!strcmp(option, "--pial_inside_hi")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.pial_inside_hi = atof(pargv[0]);
      nargsused = 1;
    }
    else if(!strcmp(option, "--pial_outside_hi")) {
      if(nargc < 1) CMDargNErr(option,1);
      adgws.pial_outside_hi = atof(pargv[0]);
      nargsused = 1;
    }


    else if(!strcmp(option, "--nsmooth")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nsmoothsurf);
      nargsused = 1;
    }
    else if(!strcmp(option, "--mm_min_inside")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%f",&T2_min_inside);
      nargsused = 1;
    }
    else if(!strcmp(option, "--mm_max_inside")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%f",&T2_max_inside);
      nargsused = 1;
    }
    else if(!strcmp(option, "--mm_min_outside")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&T2_min_outside);
      nargsused = 1;
    }
    else if(!strcmp(option, "--mm_max_outside")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&T2_max_outside);
      nargsused = 1;
    }
    // ======== Cost function weights ================
    else if (!stricmp(option, "--intensity")) {
      parms.l_intensity = atof(pargv[0]) ;
      printf("l_intensity = %2.3f\n", parms.l_intensity);
      nargsused = 1;
    }
    else if (!stricmp(option, "--hinge")) {
      parms.l_hinge = atof(pargv[0]) ;
      printf("l_hinge = %2.3f\n", parms.l_hinge);
      nargsused = 1;
    }
    else if (!stricmp(option, "--spring_nzr")) {
      parms.l_spring_nzr = atof(pargv[0]) ;
      printf("l_spring_nzr = %2.3f\n", parms.l_spring_nzr) ;
      nargsused = 1;
    }
    else if (!stricmp(option, "--spring")) {
      parms.l_spring = atof(pargv[0]) ;
      printf("l_spring = %2.3f\n", parms.l_spring) ;
      nargsused = 1;
    }
    else if (!stricmp(option, "--tspring")) {
      parms.l_tspring = atof(pargv[0]) ;
      printf("l_tspring = %2.3f\n", parms.l_tspring) ;
      nargsused = 1;
    }
    else if (!stricmp(option, "--nspring")) {
      parms.l_nspring = atof(pargv[0]) ;
      printf("l_nspring = %2.3f\n", parms.l_nspring) ;
      nargsused = 1;
    }
    else if (!stricmp(option, "--curv")){
      parms.l_curv = atof(pargv[0]) ;
      printf("l_curv = %2.3f\n", parms.l_curv) ;
      nargsused = 1;
    }
    else if (!stricmp(option, "--surf-repulse")){
      sscanf(pargv[0],"%f",&parms.l_surf_repulse);
      printf("l_surf_repulse = %2.3f\n", parms.l_surf_repulse);
      nargsused = 1;
    }
    else if (!stricmp(option, "--repulse")){
      sscanf(pargv[0],"%f",&parms.l_repulse);
      printf("l_repulse = %2.3f\n", parms.l_repulse) ;
      nargsused = 1;
    }
    // ======== End Cost function weights ================
    else if (!stricmp(option, "--n_averages")){
      sscanf(pargv[0],"%d",&n_averages);
      nargsused = 1;
    }
    else if (!stricmp(option, "--shrink")){
      sscanf(pargv[0],"%lf",&shrinkThresh);
      nargsused = 1;
    }
    else if(!strcmp(option, "--debug-vertex")){
      if(nargc < 1) CMDargNErr(option,1);
      Gdiag_no = atoi(pargv[0]) ;
      printf("Gdiag_no set to %d\n",Gdiag_no);
      nargsused = 1;
    }
    else if(!strcmp(option, "--ripflag-out")){
      // output ripflag overlay
      if(nargc < 1) CMDargNErr(option,1);
      ripflagout = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--sd")){
      if(nargc < 1) CMDargNErr(option,1);
      printf("using %s as SUBJECTS_DIR...\n", pargv[0]) ;
      setenv("SUBJECTS_DIR",pargv[0],1);
      nargsused = 1;
    }
    else if(!strcmp(option, "--s")){
      // This is meant to make it easier to run by populating the input paths
      // with files known from the subjectsdir. Might not work anymore
      if(nargc < 4) CMDargNErr(option,4);
      subject = pargv[0];
      hemi = pargv[1];
      insurfname = pargv[2];
      outsurfname = pargv[3];
      nargsused = 4;
      // Use aparc if input is white.preaparc (output is white) or
      // input is white (output is pial)
      if(strcmp(insurfname,"white.preaparc")==0 || strcmp(insurfname,"white")==0) UseAParc=1;
    }
    else if(!strcasecmp(option, "--segvolname")){
      // Applies to --s
      if(nargc < 1) CMDargNErr(option,1);
      segvolname = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--involname")){
      // Applies to --s
      if(nargc < 1) CMDargNErr(option,1);
      involname = pargv[0];
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--thickness")){
      // This appears to give the same result as mris_make_surfaces
      // except in a few vertices "near" the edge of the ripped
      // region. The way the thickness calc works is that it searchs
      // within a nbhd_size hop radius for the closest vertex on the
      // other surface excluding ripped vertices. nbhd_size=20 is
      // pretty big (about 6mm FWHM), so a ripped vertex might have an
      // influence pretty far away.  However, most of the time, the
      // closest vertex is within a few hops, so you don't usually see
      // effects far away, but they certainly can be there.
      if(nargc < 5) {
	printf("ERROR: usage --thickness white pial nbhd_size(20) maxthickness(5) out\n");
	exit(1);
      }
      surf = MRISread(pargv[0]);
      if(surf==NULL) exit(1);
      MRIScomputeMetricProperties(surf); // probably don't need to do this
      MRISsaveVertexPositions(surf, ORIGINAL_VERTICES) ;
      err = MRISreadVertexPositions(surf, pargv[1]);
      if(err) exit(1);
      MRIScomputeMetricProperties(surf);
      int nbhd_size;
      sscanf(pargv[2],"%d",&nbhd_size);
      float max_thickness;
      sscanf(pargv[3],"%f",&max_thickness);
      MRISmeasureCorticalThickness(surf, nbhd_size, max_thickness);
      err = MRISwriteCurvature(surf, pargv[4]);
      if(err) exit(1);
      exit(0);
      nargsused = 4;
    } 
    else if(!strcasecmp(option, "--curv-map")){
      // This gives the same result as mris_make_surfaces for
      // white. For pial, there is a band around the ripped vertices
      // where it is different (but the same in the center). The
      // center is the same because, in MMS, the curv would have been
      // computed from the white and the pial is the same as the white
      // in the ripped vertices. This function differs from MMS
      // because the ripped vertices are included in the curv
      // calculation for vertices near the ripped area (and vice
      // versa). This has implications for the surface-based
      // registration (with white.preaparc).
      if(nargc < 4) {
	printf("ERROR: usage --curv-map surf nbrs(2) curvature_avgs(10) out\n");
	exit(1);
      }
      surf = MRISread(pargv[0]);
      if(surf==NULL) exit(1);
      int curvature_avgs;
      sscanf(pargv[1],"%d",&nbrs);
      sscanf(pargv[2],"%d",&curvature_avgs);
      printf("insurf  %s, nbrs %d, curvature_avgs %d\n",pargv[0],nbrs,curvature_avgs);
      if(nbrs > 1) MRISsetNeighborhoodSizeAndDist(surf, nbrs) ;
      MRIScomputeMetricProperties(surf);
      MRIScomputeSecondFundamentalForm(surf) ;
      MRISuseMeanCurvature(surf) ;
      MRISaverageCurvatures(surf, curvature_avgs) ;
      err = MRISwriteCurvature(surf, pargv[3]);
      if(err) exit(1);
      exit(0);
      nargsused = 3;
    } 
    else if(!strcasecmp(option, "--area-map")){
      // This appears to give the same result as mris_make_surfaces
      // for white. For pial, it is the same except in vertices of the
      // ripped region. In MRISwriteArea(), the area of non-ripped
      // vertices is copied into the curv, then curv is saved. In MMS,
      // the curv field already has something in it after pial
      // placement, so it gets copied too. This does not happen for
      // white because the surface is unripped prior to this
      // computation.
      if(nargc < 2) {
	printf("ERROR: usage --area surf out\n");
	exit(1);
      }
      surf = MRISread(pargv[0]);
      if(surf==NULL) exit(1);
      MRIScomputeMetricProperties(surf);
      err = MRISwriteArea(surf, pargv[1]);
      if(err) exit(1);
      exit(0);
      nargsused = 2;
    } 
    else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nthreads);
      #ifdef _OPENMP
      omp_set_num_threads(nthreads);
      #endif
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--max-threads")){
      nthreads = 1;
      #ifdef _OPENMP
      nthreads = omp_get_max_threads();
      omp_set_num_threads(nthreads);
      #endif
    } 
    else if(!strcasecmp(option, "--max-threads-1") || !strcasecmp(option, "--max-threads-minus-1")){
      nthreads = 1;
      #ifdef _OPENMP
      nthreads = omp_get_max_threads()-1;
      if(nthreads < 0) nthreads = 1;
      omp_set_num_threads(nthreads);
      #endif
    } 
    else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* --------------------------------------------- */
static void check_options(void) {
  if(surftype == -1){
    printf("ERROR: must specify surface type --white or --pial\n");
    exit(1);
  }
  if(insurfpath == NULL && subject == NULL){
    printf("ERROR: no input surface set\n");
    exit(1);
  }
  if(outsurfpath == NULL && subject == NULL){
    printf("ERROR: no output surface set\n");
    exit(1);
  }
  if(insurfpath != NULL && subject != NULL){
    printf("ERROR: cannot use both --i and --s\n");
    exit(1);
  }
  if(outsurfpath != NULL && subject != NULL){
    printf("ERROR: cannot use both --o and --s\n");
    exit(1);
  }

  // When using --s  (might not work anymore)
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(insurfpath == NULL && subject != NULL){
    sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,insurfname);
    insurfpath = strcpyalloc(tmpstr);
    // Turn off surface smoothing unless input is orig (mris_make_surfaces)
    if(strcmp(insurfname,"orig")!=0) nsmoothsurf = 0;
    if(UseAParc){
      sprintf(tmpstr,"%s/%s/label/%s.%s.annot",SUBJECTS_DIR,subject,hemi,aparcname);
      aparcpath = strcpyalloc(tmpstr);
    }
    sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,outsurfname);
    outsurfpath = strcpyalloc(tmpstr);
    sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,involname);
    involpath = strcpyalloc(tmpstr);
    sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,segvolname);
    segvolpath = strcpyalloc(tmpstr);
    sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,wmvolname);
    wmvolpath = strcpyalloc(tmpstr);
  }

  if(mmvol){
    if(surftype != GRAY_CSF){
      printf("ERROR: surface type must be pial with multimodal input\n");
      exit(1);
    }
    if(whitesurfpath == NULL){
      printf("ERROR: need to specify --white-surf with multimodal input\n");
    }
  }

  if(surftype == GRAY_CSF && repulsesurfpath == NULL){
    printf("ERROR: must supply a repulse surface when placing pial surface\n");
    exit(1);
  }

  if(pinlabel &&  whitesurfpath == NULL){
    printf("ERROR: must spec --white-surf with --pin-medial-wall\n");
    exit(1);
  }

  if(coversegpath){
    printf("Reading in %s",coversegpath);
    mri_cover_seg = MRIread(coversegpath);
    if(!mri_cover_seg) exit(1);
  }

  return;
}



#include "mris_place_surface.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mris_place_surface_help_xml,mris_place_surface_help_xml_len);
}

static void
print_help(void)
{
  print_usage() ;
  exit(1) ;
}
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  return;
}

int MRISripSegs(MRIS *surf, MRI *seg, const int *SegNo, const int nSegNos, 
		const double dmin, const double dmax, const double dstep)
{
  int vno, nripped=0, k;

  printf("Starting MRISripSegs() d = (%g %g %g) segnos: ",dmin,dmax,dstep);
  for(k=0; k < nSegNos; k++) printf("%d ",SegNo[k]);
  printf("\n");

  for(vno=0; vno < surf->nvertices; vno++){
    VERTEX *v;
    int segid;
    double   xv, yv, zv, xs, ys, zs, d, val ;

    v = &(surf->vertices[vno]);
    if(v->ripflag)  continue ;

    for (d = dmin ; d <= dmax ; d += dstep) {
      xs = v->x + d*v->nx ;
      ys = v->y + d*v->ny ;
      zs = v->z + d*v->nz ;

      // Sample the aseg at this distance
      MRISsurfaceRASToVoxelCached(surf, seg, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeType(seg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
      segid = nint(val) ;
      for(k=0; k < nSegNos; k++){
	if(segid == SegNo[k]){
	  //if(!IS_WMSA(segid)) continue;
	  //printf("Ripping vertex %d at %g %g %g  depth=%g seg=%d\n",vno,round(xv),round(yv),round(zv),d,segid);
	  v->ripflag = 1;
	  nripped ++;
	  break;
	}
      } // k
      if(v->ripflag) break; // no need to continue along normal
    } // dist

  } // vertex

  printf("MRISripSegs(): %g %g %g ripped %d\n",dmin,dmax,dstep,nripped);
  return(nripped);
}


int MRISripBasalGanglia(MRIS *surf, MRI *seg, const double dmin, const double dmax, const double dstep)
{
  int vno, nripped=0;
  int indices[100], nindices=0, n;

  if(! surf->ct){
    printf("ERROR: MRISripPutamenNucAcc(): surface must have annotation\n");
    return(-1);
  }

  nindices=0;
  CTABfindName(surf->ct, "medialorbitofrontal", &indices[nindices++]);
  CTABfindName(surf->ct, "rostralanteriorcingulate", &indices[nindices++]);
  CTABfindName(surf->ct, "insula", &indices[nindices++]);

  for(vno=0; vno < surf->nvertices; vno++){
    VERTEX *v;
    int segid, hit, index;
    double   xv, yv, zv, xs, ys, zs, d, val ;

    v = &(surf->vertices[vno]);
    if(v->ripflag)  continue ;

    CTABfindAnnotation(surf->ct, v->annotation, &index);
    hit = 0;
    for(n=0; n < nindices; n++){
      if(index == indices[n]){
	hit = 1;
	break;
      }
    }
    if(! hit) continue;

    for (d = dmin ; d <= dmax ; d += dstep) {
      xs = v->x + d*v->nx ;
      ys = v->y + d*v->ny ;
      zs = v->z + d*v->nz ;

      // Sample the aseg at this distance
      MRISsurfaceRASToVoxelCached(surf, seg, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeType(seg, xv, yv, zv, &val, SAMPLE_NEAREST) ;
      segid = nint(val) ;

      // Add external and extreme capsules?
      if(segid != Left_Putamen && segid != Right_Putamen &&
	 segid != Left_Caudate && segid != Right_Caudate &&
	 segid != Left_Claustrum && segid != Right_Claustrum &&
	 segid != Left_Accumbens_area && segid != Right_Accumbens_area) continue;

      v->ripflag = 1;
      nripped ++;
      break;
    }
  }

  printf("MRISripBasalGanglia(): %g %g %g ripped %d\n",dmin,dmax,dstep,nripped);
  return(nripped);
}


int RIP_MNGR::RipVertices(void)
{

  if(RipFreeze){
    printf("Ripping frozen voxels\n");
    RipSegNo[nRipSegs++] = 247;
  }
  if(RipWMSA){
    printf("Ripping WMSA voxels\n");
    RipSegNo[nRipSegs++] = 77;
    RipSegNo[nRipSegs++] = 78;
    RipSegNo[nRipSegs++] = 79;
  }
  if(RipLesion){
    printf("Ripping Lesion voxels\n");
    RipSegNo[nRipSegs++] = 25;
    RipSegNo[nRipSegs++] = 57;
  }

  if(riplabelfile){
    printf("Ripping vertices not in label %s\n",riplabelfile);
    LABEL *riplabel = LabelRead("",riplabelfile);
    if(riplabel == NULL) exit(1);
    MRISripNotLabel(surf,riplabel);
    LabelFree(&riplabel);
  }

  int ripsurfneeded = 0;
  if(RipMidline || RipBG || nRipSegs){
    ripsurfneeded = 1;
    if(seg == NULL){
      printf("ERROR: need seg for ripping\n");
      exit(1);
    }
  }

  if(ripsurffile && !ripsurfneeded){
    printf("INFO: ripsurf specified but not needed, so ignoring\n");
    return(0);
  }

  if(ripsurffile){
    printf("Reading in ripping surface %s\n",ripsurffile);
    ripsurf = MRISread(ripsurffile);
    if(ripsurf==NULL) exit(1);
    if(ripsurf->nvertices != surf->nvertices){
      printf("ERROR: ripsurf dim mismatch %d %d\n",ripsurf->nvertices,surf->nvertices);
      exit(1);
    }
    if(aparcpath) {
      printf("Reading in aparc %s for ripsurf\n",aparcpath);
      if(MRISreadAnnotation(ripsurf, aparcpath) != NO_ERROR)
	ErrorExit(ERROR_NOFILE, "%s: could not read annotation",aparcpath) ;
    }
  }
  else {
    if(ripsurfneeded){
      printf("INFO: rip surface needed but not specified, so using input surface\n");
    }
    ripsurf = surf;
  }

  // Merging of surface and volume. This volume-based ripping
  // really needs to be done with a white-type surface (ie, not pial)
  if(seg && RipMidline){
    // probably want to use white for this
    printf("Freezing midline and others\n");  fflush(stdout);
    MRISripMidline(ripsurf, seg, invol, hemi, surftype, 0) ;
  }
  if(RipBG){
    // probably want to use white for this
    printf("Ripping BG\n");
    MRISripBasalGanglia(ripsurf, seg, dmin,dmax,dstep);
  }
  if(nRipSegs){
    // probably want to use white for this
    printf("Ripping WMSA\n");
    MRISripSegs(ripsurf, seg, RipSegNo, nRipSegs, dmin,dmax,dstep);
  }

  if(ripsurffile){
    // Copy ripflags back into the input surface
    int vno;
    for(vno=0; vno < surf->nvertices; vno++){
      if(ripsurf->vertices[vno].ripflag)  
	surf->vertices[vno].ripflag = 1;
    }
    MRISfree(&ripsurf);
  }

  return(0);
}

/*!
\fn int MRISpinMedialWallToWhite(MRIS *surf, const LABEL *cortex)
\brief Set the current coordinates to that of white{xyz} in vertices that
are not represented in the label. The label is named "cortex" but it
could be anything. This function is used after the pial surface is placed
without freezing the medial wall next to hippocampus and amyg. This function
can then be used to force the pial surface back to the white surface in all
of the medial wall.
*/
int MRISpinMedialWallToWhite(MRIS *surf, const LABEL *cortex)
{
  int vno,i;
  int *InLabel;
  VERTEX *v;

  InLabel = (int*) calloc(sizeof(int),surf->nvertices);
  for (i = 0; i < cortex->n_points; i++)
    InLabel[cortex->lv[i].vno] = 1;

  for(vno = 0; vno < surf->nvertices; vno++){
    if(InLabel[vno]) continue;
    v = &(surf->vertices[vno]);
    v->x = v->whitex;
    v->y = v->whitey;
    v->z = v->whitez;
  }

  free(InLabel);
  return(0);
}




