/**
 * @file  mri_place_surface.c
 * @brief Places surface based on an intensity input image. This is meant to provide
 * replacement functionality for mris_make_surfaces in a form that is easier to 
 * maintain.
 */
/*
 * Original Author: Douglas N Greve (but basically a rewrite of mris_make_surfaces by BF)
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2017/02/15 21:04:18 $
 *    $Revision: 1.246 $
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
#include "romp_support.h"

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

static char vcid[] =
"$Id: mri_glmfit.c,v 1.246 2017/02/15 21:04:18 greve Exp $";
const char *Progname = "mri_glmfit";

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
// Note "max_thickness" here is really a limit on the distance that CBV will
// search along the normal (inside and out)
float max_thickness = 5.0 ;
int vavgs = 5 ;
int nthreads = 1;
int nbrs = 2;
int nsmoothsurf = 5 ;

char *SUBJECTS_DIR;
char *insurfpath = NULL;
char *outsurfpath = NULL;
char *involpath=NULL;
char *segvolpath=NULL;
char *wmvolpath=NULL;
char *aparcpath=NULL;
char *initsurfpath=NULL;

char *subject = NULL,*hemi = NULL,*insurfname = NULL, *outsurfname = NULL;
char *involname="brain.finalsurfs.mgz", *segvolname="aseg.presurf.mgz",*wmvolname="wm.mgz",*aparcname="aparc";

char tmpstr[2000];
int err=0;
int longitudinal = 0;
int surftype = -1; //GRAY_WHITE; // GRAY_CSF
int UseAParc = 0;

/*--------------------------------------------------*/
int main(int argc, char **argv) 
{
  int nargs, i, msec;
  double        spring_scale = 1;
  MRIS *surf;
  MRI *invol, *seg, *wm;
  Timer timer ;
  char *cmdline2, cwd[2000];
  //char *field;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
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
  printf("%s\n",vcid);
  printf("%s\n",MRISurfSrcVersion());
  printf("\n");
  printf("cd %s\n",cwd);
  printf("setenv SUBJECTS_DIR %s\n",getenv("SUBJECTS_DIR"));
  printf("%s\n",cmdline2);
  printf("\n");
  fflush(stdout);

  // Note: in long stream orig = orig_white
  AutoDetGWStats adgws;
  err = adgws.AutoDetectStats(subject, hemi);
  if(err) exit(1);
  double inside_hi=0, border_hi=0, border_low=0, outside_low=0, outside_hi=0,current_sigma=0;
  int n_averages=0, n_min_averages=0;
  if(surftype == GRAY_WHITE){
    current_sigma = white_sigma ;
    n_averages = max_white_averages;     
    n_min_averages = min_white_averages; 
    inside_hi = adgws.white_inside_hi;
    border_hi = adgws.white_border_hi;
    border_low = adgws.white_border_low;
    outside_low = adgws.white_outside_low;
    outside_hi = adgws.white_outside_hi;
  }
  if(surftype == GRAY_CSF){
    current_sigma = pial_sigma ;
    n_averages = max_pial_averages;     
    n_min_averages = min_pial_averages; 
    inside_hi = adgws.pial_inside_hi;
    border_hi = adgws.pial_border_hi;
    border_low = adgws.pial_border_low;
    outside_low = adgws.pial_outside_low;
    outside_hi = adgws.pial_outside_hi;
  }

  printf("Reading in input surface %s\n",insurfpath);
  surf = MRISread(insurfpath);
  if(surf==NULL) exit(1);
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
  MRISsaveVertexPositions(surf, ORIGINAL_VERTICES) ;
  MRISsaveVertexPositions(surf, WHITE_VERTICES) ;
  MRISsetVals(surf,-1) ;  /* clear white matter intensities */
  MRISfaceMetric(surf,0);
  MRISedgeMetric(surf,0);
  MRIScornerMetric(surf,0);
  MRISprettyPrintSurfQualityStats(stdout, surf);

  if(aparcpath) {
    printf("Reading in aparc %s\n",aparcpath);
    if (MRISreadAnnotation(surf, aparcpath) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read annotation",aparcpath) ;
  }
  else printf("Not reading in aparc\n");

  printf("Reading in input volume %s\n",involpath);
  invol = MRIread(involpath);
  if(invol==NULL) exit(1);

  //  =========== intensity volume preprocessing ==============
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
    MRImask(invol, mri_labeled, invol, BRIGHT_LABEL, 0) ;
    MRImask(invol, mri_labeled, invol, BRIGHT_BORDER_LABEL, 0) ;
  }

  printf("Reading in seg volume %s\n",segvolpath);
  seg = MRIread(segvolpath);
  if(seg==NULL) exit(1);

  // First interaction of surface and volume
  if(seg && surftype == GRAY_CSF){
    printf("Freezing midline and others\n");  fflush(stdout);
    MRISripMidline(surf, seg, invol, hemi, surftype, 0) ;
  }

  if(initsurfpath){
    if(longitudinal){
      // Save initial surface location (white) into TMP_VERTICES (v->tx, v->ty, v->tz)
      // This will be used for the repulsion cost
      MRISsaveVertexPositions(surf, TMP_VERTICES);
    }
    // Now read in true init surface
    printf("Reading vertex coordinates from %s\n",initsurfpath);
    err = MRISreadVertexPositions(surf, initsurfpath);
    if(err) exit(1);
    MRISremoveIntersections(surf); // done in mris_make_surfaces
    MRISsaveVertexPositions(surf, PIAL_VERTICES) ;
    if(longitudinal) {
      // In the longitudinal case, reset the starting position to be
      // slightly inside the orig_pial between final white (insurfpath) and 
      // init pial. Couldn't this cause intersections?
      MRISblendXYZandTXYZ(surf, 0.75f, 0.25f);
    }
    // This will be used to keep track of which vertices found pial
    // surface in previous cycle.  Should already be clear, just
    // including from mris_make_surfaces
    MRISclearMark2s(surf) ;  
    MRIScomputeMetricProperties(surf);
    MRISfaceMetric(surf,0);
    MRISedgeMetric(surf,0);
    MRIScornerMetric(surf,0);
  }

  timer.reset() ;
  printf("n_averages %d\n",n_averages);
  for (i = 0 ;  n_averages >= n_min_averages ; n_averages /= 2, current_sigma /= 2, i++) {

    printf("Iteration %d =========================================\n",i);
    printf("n_averages=%d, current_sigma=%g\n",n_averages,current_sigma); fflush(stdout);

    if(seg && surftype == GRAY_WHITE){
      printf("Freezing midline and others\n");  fflush(stdout);
      // This rips the midline vertices so that they do not move (thus
      // "fix"). It may also rip some vertices near the putamen and
      // maybe near lesions. It may set v->marked2 which will
      // influence the cortex.label. At onset, it will set all marked
      // and marked2 to 0. At then end it will set all marked=0.
      // It does not unrip any vertices, so, unless they are unripped
      // at some other point, the number of ripped vertices will 
      // increase.
      MRISripMidline(surf, seg, invol, hemi, surftype, 0) ;
    }

    parms.sigma = current_sigma ;
    parms.n_averages = n_averages ;
    
    if(surftype == GRAY_CSF){
      printf("Masking bright non-wm for pial surface mid_gray = %g\n",adgws.MID_GRAY);
      // Replace brightborder voxels with MID_GRAY (can be done once?)
      // Why would you want to do this? The brightborder voxels are >= 100,
      // but this would make them look like GM and force the pial outside of them.
      // This would happen, eg, for a bright vessel in cortex
      MRImask(invol, mri_labeled, invol, BRIGHT_BORDER_LABEL, adgws.MID_GRAY);
      // Replace bright voxels with 255 (this gets changed below)
      // Not sure why this is needed except that it gets set to 0 below which
      // could look like a strong negative gradient
      MRImask(invol, mri_labeled, invol, BRIGHT_LABEL, 255) ;
      /* From mris_make_surfaces: Replace bright stuff such as eye
        sockets with 255.  Simply zeroing it out would make the border
        always go through the sockets, and ignore subtle local minima
        in intensity at the border of the sockets.  Will set to 0
        after border values have been computed so that it doesn't mess
        up gradients.  */
    }

    printf("Computing border values \n");
    // The outputs are set in each vertex structure:
    //   v->val2 = current_sigma; // smoothing level along gradient used to find the target
    //   v->val  = max_mag_val; // intensity at target location
    //   v->d = max_mag_dist;   // dist to target along normal
    //   v->mean = max_mag;     // derivative at target intensity
    //   v->marked = 1;         // vertex has good data
    //   v->targx = v->x + v->nx * v->d; // same for y and z
    // Note "max_thickness" here is really a limit on the distance that CBV will
    // search along the normal (inside and out)
    MRIScomputeBorderValues(surf, invol, NULL, inside_hi,border_hi,border_low,outside_low,outside_hi,
			    current_sigma, 2*max_thickness, parms.fp, surftype, NULL, 0, parms.flags,seg,-1,-1) ;
    // Note: 3rd input (NULL) was "mri_smooth" in mris_make_surfaces, but
    // this was always a copy of the input (mri_T1 or invol); it is not used in CBV

    if(surftype == GRAY_CSF){
      // Replace bright voxels with 0 (mask them out)
      // This undoes some of the masking above (or vice versa)
      MRImask(invol, mri_labeled, invol, BRIGHT_LABEL, 0) ;
    }

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

    /*There are frequently regions of gray whose intensity is fairly
      flat. We want to make sure the surface settles at the innermost
      edge of this region, so on the first pass, set the target
      intensities artificially high so that the surface will move
      all the way to white matter before moving outwards to seek the
      border (I know it's a hack, but it improves the surface in
      a few areas. The alternative is to explicitly put a gradient-seeking
      term in the cost functional instead of just using one to find
      the target intensities).
    */

    INTEGRATION_PARMS_copy(&old_parms, &parms) ;

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
    MRISpositionSurface(surf, invol, invol, &parms);

    //sprintf(tmpstr,"lh.place.postpos%02d",i);
    //MRISwrite(surf, tmpstr);

    old_parms.start_t = parms.start_t ;
    INTEGRATION_PARMS_copy(&parms, &old_parms) ;

    if(!n_averages)
      break ; 

  } // end major loop placing the white surface using

  MRISremoveIntersections(surf); //matches mris_make_surface

  printf("\n\n");
  printf("Writing output to %s\n",outsurfpath);
  err = MRISwrite(surf,outsurfpath);
  if(err){
    printf("ERROR: writing to %s\n",outsurfpath);
    exit(1);
  }

  msec = timer.milliseconds() ;
  printf("#ET# mris_place_surface %5.2f minutes\n", (float)msec/(60*1000.0f));
  printf("#VMPC# mris_make_surfaces VmPeak  %d\n",GetVmPeak());
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
    else if(!strcmp(option, "--long")) longitudinal = 1;
    else if(!strcmp(option, "--white")) surftype = GRAY_WHITE;
    else if(!strcmp(option, "--pial"))  surftype = GRAY_CSF;
    else if(!strcmp(option, "--use-aparc")) UseAParc = 1;
    else if(!strcmp(option, "--i")){
      if(nargc < 1) CMDargNErr(option,1);
      insurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--o")){
      if(nargc < 1) CMDargNErr(option,1);
      outsurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--init")){
      if(nargc < 1) CMDargNErr(option,1);
      initsurfpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--nsmooth")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nsmoothsurf);
      nargsused = 1;
    }
    else if(!strcmp(option, "--sd")){
      if(nargc < 1) CMDargNErr(option,1);
      printf("using %s as SUBJECTS_DIR...\n", pargv[0]) ;
      setenv("SUBJECTS_DIR",pargv[0],1);
      nargsused = 1;
    }
    else if(!strcmp(option, "--s")){
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
    else if(!strcasecmp(option, "--max-thickness")){
      // Note "max_thickness" here is really a limit on the distance that CBV will
      // search along the normal (inside and out)
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%f",&max_thickness);
      nargsused = 1;
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
  if(surftype == -1){
    printf("ERROR: must specify surface type --white or --pial\n");
    exit(1);
  }

  return;
}



/* --------------------------------------------- */
static void print_usage(void) 
{
printf("\n");
printf("USAGE: ./mris_place_surface\n");
printf(" --s subject hemi insurfname outsurfname\n");
printf("\n");
}


/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n");
  printf("\n");
  exit(1) ;
}

/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  return;
}

