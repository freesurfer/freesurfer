/**
 * @file  mris_register_label_map.c
 * @brief register the resting state patterns from one label map to another based on priors
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2013/04/09 12:29:19 $
 *    $Revision: 1.1 $
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



/*!
\file mris_register_label_map.c
\brief program to register an individual subject's resting state map to a group average
\author Bruce Fischl

*/


// $Id: mris_register_label_map.c,v 1.1 2013/04/09 12:29:19 fischl Exp $

/*
  BEGINHELP

  ENDHELP
*/

/*
  BEGINUSAGE

  ENDUSAGE
*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>

#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mrisurf.h"

static MRI *sample_fixed(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats) ;
static double compute_error(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_mov_avg, MRI *mri_stats) ;
static double sample_stats(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats, double x, double y, double z, double *pvar) ;
static double sample_vertex_error(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_lh_fixed, MRI_SURFACE *mris_rh_fixed, MHT *mht, MRI *mri_stats, MRI *mri_corrmat_mov, MRI *mri_label_avg_mov, int vno, double x, double y, double z);
static int compute_warp_gradient(MRI_SURFACE *mris_lh_mov, MRI_SURFACE *mris_rh_mov, MRI_SURFACE *mris_lh_fixed,  MRI_SURFACE *mris_rh_fixed,  MHT *mht_lh, MHT *mht_rh, MRI *mri_cmat_permuted,  MRI *mri_stats, MRI *mri_label_avg, MRI *mri_corrmat_mov, double dt) ;
static int  warp_surface(MRI_SURFACE *mris_lh_mov, MRI_SURFACE *mris_rh_mov, MRI_SURFACE *mris_lh_fixed, MRI_SURFACE *mris_rh_fixed, MRI *mri_cmat_mov, MRI *mri_stats, double tol, LABEL *area, int offset, double dt) ;
static int compute_vertex_permutation(MRI_SURFACE *mris_mov_mov, MRI_SURFACE *mris_fixed, MHT *mht,  int *vertices) ;
static MRI *average_within_label(MRI *mri_cmat, LABEL *area, int offset, MRI *mri_avg) ;
static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;
static MRI *compute_mean_and_variance(MRI **mri_label_avg, int nsubjects)  ;

static char vcid[] = "$Id: mris_register_label_map.c,v 1.1 2013/04/09 12:29:19 fischl Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

static char *cmat_name = NULL ;
static char *trgsubject = NULL ;
static char  *label_name = NULL ;
static char *prior_name = NULL ;
static char *TempVolFile=NULL;
static char *SUBJECTS_DIR = NULL;
static char **subjects = NULL ;
static int nsubjects = 0 ;
static char *hemi = NULL ;
static int offset = 0 ;
static char *output_name = NULL ;
static int write_diags = 1 ;
static double tol = 0.01 ;   // terminate when error func % change is smaller than this
static double dt = .1 ;

#define MAX_SUBJECTS 1000

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) 
{
  int         nargs, n, ic;
  char        *subject, fname[STRLEN] ;
  MRI         *mri_cmat, *mri_prior, *mri_stats, *mri_label_avg[MAX_SUBJECTS], *mri_cmat_mov ;
  LABEL       *labels[MAX_SUBJECTS], *target_area ;
  MRI_SURFACE *mris_lh_mov, *mris_rh_mov, *mris_lh_fixed, *mris_rh_fixed ;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);
  dump_options(stdout);

  if (SUBJECTS_DIR == NULL)
  {
    SUBJECTS_DIR = getenv("SUBJECTS_DIR");
    if (SUBJECTS_DIR == NULL) {
      printf("ERROR: SUBJECTS_DIR not defined in environment\n");
      exit(1);
    }
  }

  for (n = 0 ; n < nsubjects ; n++)
  {
    subject = subjects[n] ;
    printf("processing subject %s\n", subject) ;
    sprintf(fname, "%s/%s/surf/%s", SUBJECTS_DIR, subject, cmat_name) ;
    mri_cmat = MRIread(fname) ;
    if (!mri_cmat)
      ErrorExit(ERROR_NOFILE, "%s: could not read cmat from %s", Progname, fname) ;

    sprintf(fname, "%s/%s/label/%s.%s", SUBJECTS_DIR, subject, hemi, label_name) ;
    labels[n] = LabelRead("", fname) ;
    if (!labels[n])
      ErrorExit(ERROR_NOFILE, "%s: could not read label from %s", Progname, fname) ;
    if (stricmp(hemi, "rh") == 0)
      offset = mri_cmat->width/2 ;  // lower right block is for rh
    mri_label_avg[n] = average_within_label(mri_cmat, labels[n], offset, NULL) ;
    MRIfree(&mri_cmat) ;

    if (write_diags)
    {
      sprintf(fname, "%s.%s.%s.label_avg.mgz", hemi, subject, output_name) ;
      mri_label_avg[n]->width /= 2 ;
      MRIwrite(mri_label_avg[n], fname) ;
      mri_label_avg[n]->width *= 2 ;
    }
  }
  mri_stats = compute_mean_and_variance(mri_label_avg, nsubjects) ;
  if (write_diags)
  {
    sprintf(fname, "%s.%s.label_avg.mgz", hemi, output_name) ;
    mri_stats->width /= 2 ;
    MRIwriteFrame(mri_stats, fname, 0) ;
    mri_stats->width *= 2 ;
  }
  ic  = 5 ;
  switch (MRInvox(mri_label_avg[0])/2)
  {
  case 10242: ic = 5 ; break ;
  default: ErrorExit(ERROR_UNSUPPORTED, "only IC 5 supported (nvertices = %d detected)\n", MRInvox(mri_label_avg[0])/2);
  }

  sprintf(fname, "%s/fsaverage%d/surf/lh.sphere", SUBJECTS_DIR, ic) ;
  mris_lh_mov = MRISread(fname) ;
  if (!mris_lh_mov)
    ErrorExit(ERROR_NOFILE, "%s: could not read *surface from %s", Progname, fname) ;
  sprintf(fname, "%s/fsaverage%d/surf/rh.sphere", SUBJECTS_DIR, ic) ;
  mris_rh_mov = MRISread(fname) ;
  if (!mris_rh_mov)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s", Progname, fname) ;

  mris_lh_fixed = MRISread(fname) ;
  if (!mris_lh_fixed)
    ErrorExit(ERROR_NOFILE, "%s: could not read *surface from %s", Progname, fname) ;
  sprintf(fname, "%s/fsaverage%d/surf/rh.sphere", SUBJECTS_DIR, ic) ;
  mris_rh_fixed = MRISread(fname) ;
  if (!mris_rh_fixed)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s", Progname, fname) ;

  MRISsetNeighborhoodSize(mris_lh_mov, 2) ;
  MRISsetNeighborhoodSize(mris_rh_mov, 2) ;
  MRISsetNeighborhoodSize(mris_lh_fixed, 2) ;
  MRISsetNeighborhoodSize(mris_rh_fixed, 2) ;
  MRIScomputeSecondFundamentalForm(mris_lh_mov) ;
  MRIScomputeSecondFundamentalForm(mris_rh_mov) ;
  MRIScomputeSecondFundamentalForm(mris_lh_fixed) ;
  MRIScomputeSecondFundamentalForm(mris_rh_fixed) ;

  sprintf(fname, "%s/fsaverage%d/label/%s.%s", SUBJECTS_DIR, ic, hemi, prior_name) ;
  mri_prior = MRIread(fname) ;
  if (!mri_prior)
    ErrorExit(ERROR_NOFILE, "%s: could not prior overlay from %s", Progname, fname) ;

  sprintf(fname, "%s/%s/surf/%s", SUBJECTS_DIR, trgsubject, cmat_name) ;
  mri_cmat_mov = MRIread(fname) ;
  if (!mri_cmat_mov)
      ErrorExit(ERROR_NOFILE, "%s: could not read movable cmat from %s", Progname, fname) ;

  sprintf(fname, "%s/fsaverage%d/label/%s.MT.thresh.label", SUBJECTS_DIR, ic, hemi) ;
  target_area = LabelRead("", fname) ;
  if (!target_area)
      ErrorExit(ERROR_NOFILE, "%s: could not read movable cmat from %s", Progname, fname) ;


  warp_surface(mris_lh_mov, mris_rh_mov, mris_lh_fixed, mris_rh_fixed, mri_cmat_mov, mri_stats, tol, target_area, offset, dt) ;
  sprintf(fname, "%s/fsaverage%d/surf/lh.sphere.%s", SUBJECTS_DIR, ic, output_name) ;
  printf("writing output surface to %s\n", fname) ;
  MRISwrite(mris_lh_mov, fname) ;
  sprintf(fname, "%s/fsaverage%d/surf/rh.sphere.%s", SUBJECTS_DIR, ic, output_name) ;
  printf("writing output surface to %s\n", fname) ;
  MRISwrite(mris_rh_mov, fname) ;
  return 0;
}
/* ------ Doxygen markup starts on the line below ---- */
/*!
\fn int parse_commandline(int argc, char **argv)
\brief Parses the command-line arguments
\param argc - number of command line arguments
\param argv - pointer to a character pointer
*/
/* ------ Doxygen markup ends on the line above ---- */
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
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--sdir"))
    {
      SUBJECTS_DIR = pargv[0] ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--prior"))
    {
      prior_name = pargv[0] ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--v"))
    {
      Gdiag_no = atoi(pargv[0]) ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--tol"))
    {
      tol = atof(pargv[0]) ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--trgsubject"))
    {
      trgsubject = pargv[0] ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--output"))
    {
      output_name = pargv[0] ;
      printf("setting output name = %s\n", output_name) ;
      fflush(stdout) ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--hemi"))
    {
      hemi = pargv[0] ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--subjects"))
    {
      int nextarg ;

      for (nextarg = 0 ; nextarg < nargc ; nextarg++)
	if (pargv[nextarg][0] == '-' && pargv[nextarg][1] == '-')
	  break ;
      printf("parsing --subjects with argc = %d\n", nextarg) ;
      fflush(stdout) ;
      nsubjects = nargsused = nextarg ;
      subjects = pargv ;
    }
    else if (!strcasecmp(option, "--label"))
    {
      label_name = pargv[0] ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--cmat"))
    {
      cmat_name = pargv[0] ;
      nargsused = 1 ;
    }
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;

    else if (!strcasecmp(option, "--temp-vol")) {
      if (nargc < 1) CMDargNErr(option,1);
      TempVolFile = pargv[0];
      nargsused = 1;
    } else {
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
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void usage_exit(void)
\brief Prints usage and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_usage(void)
\brief Prints usage and returns (does not exit)
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --subjects  SUBJECT_LIST: list of training subjects \n");
  printf("   --trgsubject SUBJECT: name of target subject \n");
  printf("   --prior      PRIOR: name of prior surface overlay \n");
  printf("   --label      LABEL: name of label for each subject \n");
  printf("   --temp-vol volfile : template volume \n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --sdir      SUBJECTS_DIR\n");
  printf("   --version   print out version and exit\n");
  printf("   --v         VNO  debug this vertex\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_help(void)
\brief Prints help and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_help(void) {
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void print_version(void)
\brief Prints version and exits
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void check_options(void)
\brief Checks command-line options
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
static void check_options(void) {
  return;
}

/* -- Doxygen markup starts on the line below (this line not needed for Doxygen) -- */
/*!
\fn static void dump_options(FILE *fp)
\brief Prints command-line options to the given file pointer
\param FILE *fp - file pointer
*/
/* ------ Doxygen markup ends on the line above  (this line not needed for Doxygen) -- */
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

static MRI *
compute_mean_and_variance(MRI **mri_label_avg, int nsubjects) 
{
  MRI    *mri_stats, *mri ;
  int    n, x, y, z, f ;
  float  val, mean, var ;

  mri = mri_label_avg[0] ;
  mri_stats = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, 2) ;
  MRIcopyHeader(mri, mri_stats) ;
  for (n = 0 ; n < nsubjects ; n++)
  {
    mri = mri_label_avg[n] ;
    for (x = 0 ; x < mri->width ; x++)
      for (y = 0 ; y < mri->height ; y++)
	for (z = 0 ; z < mri->depth ; z++)
	  for (f = 0 ; f < mri->nframes ; f++)
	  {
	    val = MRIgetVoxVal(mri, x, y, z, f) ;
	    mean = MRIgetVoxVal(mri_stats, x, y, z, 0) ;
	    var = MRIgetVoxVal(mri_stats, x, y, z, 1) ;
	    mean += val ; var += (val*val) ;
	    MRIsetVoxVal(mri_stats, x, y, z, 0, mean) ;
	    MRIsetVoxVal(mri_stats, x, y, z, 1, var) ;
	  }
  }

  for (x = 0 ; x < mri->width ; x++)
    for (y = 0 ; y < mri->height ; y++)
      for (z = 0 ; z < mri->depth ; z++)
      {
	mean = MRIgetVoxVal(mri_stats, x, y, z, 0) ;
	var = MRIgetVoxVal(mri_stats, x, y, z, 1) ;
	mean /= nsubjects ;
	var = var / nsubjects - mean*mean ;
	if (var < .0001)
	  var = .01 ; 
	MRIsetVoxVal(mri_stats, x, y, z, 0, mean) ;
	MRIsetVoxVal(mri_stats, x, y, z, 1, var) ;
      }
  return(mri_stats) ;
}

static MRI *
average_within_label(MRI *mri_cmat, LABEL *area, int offet, MRI *mri_avg) 
{
  int  n, vno, x, y, z ;
  float avg, val ;

  if (mri_avg == NULL)
  {
    mri_avg = MRIallocSequence(mri_cmat->width, mri_cmat->height, mri_cmat->depth, MRI_FLOAT, 1) ;
    MRIcopyHeader(mri_cmat, mri_avg) ;
  }

  for (x = 0 ; x < mri_cmat->width ; x++)
    for (y = 0 ; y < mri_cmat->height ; y++)
      for (z = 0 ; z < mri_cmat->depth ; z++)
      {
	for (avg = 0.0, n = 0 ; n < area->n_points ; n++)
	{
	  vno = area->lv[n].vno + offset ;
	  if (vno == x)
	    continue ;  // don't include diagonal in estimate
	  
	  val = MRIgetVoxVal(mri_cmat, x, y, z, vno) ;
	  if (x == Gdiag_no)
	  {
	    printf("vno %d, map %d, adding %2.2f\n", x, vno, val) ;
	    DiagBreak() ;
	  }
	  avg += val ;
	}
	avg /= (area->n_points-1) ;  // didn't  include diagonal
	if (fabs(avg) > 1000)
	  DiagBreak() ;
	MRIsetVoxVal(mri_avg, x, y, z, 0, avg) ;
      }

  return(mri_avg) ;
}


static int
compute_vertex_permutation(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht,  int *vertices)
{
  int    vno ;
  VERTEX *v, *vfixed ;

  for (vno = 0 ; vno < mris_mov->nvertices ; vno++)
  {
    v = &mris_mov->vertices[vno] ;
    vfixed = MHTfindClosestVertex(mht, mris_fixed, v) ;
    vertices[vno] = vfixed - mris_fixed->vertices ;
    
  }    
  
  return(NO_ERROR) ;
}

static MRI *
permute_corrmat(MRI *mri_corrmat_src, int *lh_vertices, int *rh_vertices, MRI *mri_corrmat_dst)
{
  int vno_src, vno_dst, r, offset ;
  float val ;

  if (mri_corrmat_dst == NULL)
    mri_corrmat_dst = MRIclone(mri_corrmat_src, NULL) ;

  offset = mri_corrmat_src->width/2 ;
  for (vno_src = 0 ; vno_src < mri_corrmat_src->width ; vno_src++)
  {
    if (vno_src < offset)   // lh
    {
      vno_dst = lh_vertices[vno_src] ;
      if (vno_src != vno_dst)
	DiagBreak() ;
    }
    else
    {
      vno_dst = rh_vertices[vno_src-offset] ;
      if (vno_src-offset != vno_dst)
	DiagBreak() ;
    }

    // swap rows
    for (r = 0 ; r <  mri_corrmat_src->width ; r++)
    {
      val = MRIgetVoxVal(mri_corrmat_src, vno_src, 0,0, r) ;
      MRIsetVoxVal(mri_corrmat_dst, vno_dst,0, 0, r, val) ;
    }
    // swap cols
    for (r = 0 ; r <  mri_corrmat_src->width ; r++)
    {
      val = MRIgetVoxVal(mri_corrmat_src, r, 0,0, vno_src) ;
      MRIsetVoxVal(mri_corrmat_dst, r,0, 0, vno_dst, val) ;
    }
  }

  return(mri_corrmat_dst) ;
}


#define MAX_ITERS 100
    
static double
compute_error(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_mov_avg, MRI *mri_stats)
{
  int     vno ;
  double  sse ;
  double  val, mean, var ;
  VERTEX  *v ;

  for (sse = 0.0, vno = 0 ; vno < mris_mov->nvertices ; vno++)
  {
    v = &mris_mov->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    mean = sample_stats(mris_mov, mris_fixed, mht, mri_stats, v->x, v->y, v->z, &var) ;
    val = MRIgetVoxVal(mri_mov_avg, vno, 0,0, 0) ;
    if (FZERO(var))
      var = 1.0 ;
    sse += SQR(mean-val)/var ;
  }
  return(sqrt(sse/mri_mov_avg->width)) ;
}

static double
sample_vertex_error(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_lh_fixed, MRI_SURFACE *mris_rh_fixed, MHT *mht, MRI *mri_stats, MRI *mri_corrmat_mov, MRI *mri_label_avg_mov, int vno, double x, double y, double z)
{
  double mean, var, val, error ;

  mean = sample_stats(mris_mov, mris_lh_fixed, mht, mri_stats, x, y, z, &var) ;
  val = MRIgetVoxVal(mri_label_avg_mov, vno, 0, 0, 0) ;
  error = SQR(val-mean)/var ;
  
  return(error) ;
}
static MRI *
sample_fixed(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats) 
{
  MRI    *mri ;
  double mean, var ;
  int    vno ;
  VERTEX *v ;

  mri = MRIalloc(mris_mov->nvertices, 1, 1, MRI_FLOAT) ;
  MRIcopyHeader(mri_stats, mri) ;
  for (vno = 0 ; vno < mris_mov->nvertices ; vno++)
  {
    v = &mris_mov->vertices[vno] ;
    mean = sample_stats(mris_mov, mris_fixed, mht, mri_stats, v->x, v->y, v->z, &var) ;
    MRIsetVoxVal(mri, vno, 0, 0, 0, mean) ;
  }
  return(mri) ;
}
static double
sample_stats(MRI_SURFACE *mris_mov, MRI_SURFACE *mris_fixed, MHT *mht, MRI *mri_stats, double x, double y, double z, double *pvar)
{
  int   fno ;
  FACE  *face ;
  double fdist, val0, val1, val2, mean, var ;

  MHTfindClosestFaceGeneric(mht, mris_fixed, x, y, z, 8000, -1, 0, &face, &fno, &fdist);

  // interpolate mean
  val0 = MRIgetVoxVal(mri_stats, face->v[0], 0, 0,0) ;
  val1 = MRIgetVoxVal(mri_stats, face->v[1], 0, 0,0) ;
  val2 = MRIgetVoxVal(mri_stats, face->v[2], 0, 0,0) ;
  mean = MRISsampleFace(mris_fixed, fno, CURRENT_VERTICES, x, y, z, val0, val1, val2) ;
  if (fabs(val0) > 100 || fabs(val1) > 100 || fabs(val2) > 100 || fabs(mean)>100)
    DiagBreak() ;

  // interpolate variance
  val0 = MRIgetVoxVal(mri_stats, face->v[0], 0, 0,1) ;
  val1 = MRIgetVoxVal(mri_stats, face->v[1], 0, 0,1) ;
  val2 = MRIgetVoxVal(mri_stats, face->v[2], 0, 0,1) ;
  var = MRISsampleFace(mris_fixed, fno, CURRENT_VERTICES, x, y, z, val0, val1, val2) ;
  if (fabs(val0) > 100 || fabs(val1) > 100 || fabs(val2) > 100 || fabs(var)>100)
    DiagBreak() ;

  *pvar = var ;
  return(mean) ;
}


static int
compute_warp_gradient(MRI_SURFACE *mris_lh_mov, MRI_SURFACE *mris_rh_mov, MRI_SURFACE *mris_lh_fixed,  MRI_SURFACE *mris_rh_fixed,  MHT *mht_lh, MHT *mht_rh, MRI *mri_cmat_permuted,  MRI *mri_stats, MRI *mri_label_avg, MRI *mri_corrmat_mov, double dt)
{
  int    vno ;
  VERTEX *v ;
  double  e1p, e1m, e2p, e2m, x, y, z, dx, dy, dz ;

  for (vno = 0 ; vno < mris_lh_mov->nvertices ; vno++)
  {
    v = &mris_lh_mov->vertices[vno] ;

    // sample vertex energy in plus and minus two principal directions
    x = v->x + v->e1x*dt ; y = v->y + v->e1y*dt ; z = v->z + v->e1z*dt ; 
    e1p  = sample_vertex_error(mris_lh_mov, mris_lh_fixed, mris_rh_fixed, mht_lh, mri_stats, mri_corrmat_mov, mri_label_avg, vno, x, y, z) ;
    x = v->x - v->e1x*dt ; y = v->y - v->e1y*dt ; z = v->z - v->e1z*dt ; 
    e1m  = sample_vertex_error(mris_lh_mov, mris_lh_fixed, mris_rh_fixed, mht_lh, mri_stats, mri_corrmat_mov, mri_label_avg, vno, x, y, z) ;
    x = v->x + v->e2x*dt ; y = v->y + v->e2y*dt ; z = v->z + v->e2z*dt ; 
    e2p  = sample_vertex_error(mris_lh_mov, mris_lh_fixed, mris_rh_fixed, mht_lh, mri_stats, mri_corrmat_mov, mri_label_avg, vno, x, y, z) ;
    x = v->x - v->e2x*dt ; y = v->y - v->e2y*dt ; z = v->z - v->e2z*dt ; 
    e2m  = sample_vertex_error(mris_lh_mov, mris_lh_fixed, mris_rh_fixed, mht_lh, mri_stats, mri_corrmat_mov, mri_label_avg, vno, x, y, z) ;
    if (e1p < e1m && e1p < e2p && e1p < e2m)   // move in e1p dir
    {
      dx = v->e1x*dt ; dy = v->e1y*dt ; dz = v->e1z*dt ;
    }
    else if (e1m < e2p && e1m < e2m)
    {
      dx = -v->e1x*dt ; dy = -v->e1y*dt ; dz = -v->e1z*dt ;
    }
    else if (e2m < e2p)
    {
      dx = -v->e2x*dt ; dy = -v->e2y*dt ; dz = -v->e2z*dt ;
    }
    else
    {
      dx = v->e2x*dt ; dy = v->e2y*dt ; dz = v->e2z*dt ;
    }
    v->dx += dx ; v->dy += dy ; v->dz += dz ;
  }

  return(NO_ERROR) ;
}

static int
warp_surface(MRI_SURFACE *mris_lh_mov, MRI_SURFACE *mris_rh_mov, MRI_SURFACE *mris_lh_fixed, MRI_SURFACE *mris_rh_fixed, MRI *mri_cmat_mov, MRI *mri_stats, double tol, LABEL *area, int offset, double dt) 
{
  int  *lh_vertices, *rh_vertices, iter = 0 ;
  MHT  *mht_lh, *mht_rh, *mht_lh_faces, *mht_rh_faces ;
  MRI  *mri_cmat_permuted = NULL, *mri_mov_avg ;
  double last_error, error, pct_change ;

  lh_vertices = (int *)calloc(mris_lh_fixed->nvertices, sizeof(int)) ;
  rh_vertices = (int *)calloc(mris_rh_fixed->nvertices, sizeof(int)) ;
  mht_lh = MHTfillVertexTableRes(mris_lh_fixed, NULL, CURRENT_VERTICES, ceil(mris_lh_fixed->avg_vertex_dist));
  mht_rh = MHTfillVertexTableRes(mris_rh_fixed, NULL, CURRENT_VERTICES, ceil(mris_rh_fixed->avg_vertex_dist));
  mht_lh_faces = MHTfillTableAtResolution(mris_lh_fixed, NULL, CURRENT_VERTICES, ceil(mris_lh_fixed->avg_vertex_dist)) ;
  mht_rh_faces = MHTfillTableAtResolution(mris_rh_fixed, NULL, CURRENT_VERTICES, ceil(mris_rh_fixed->avg_vertex_dist)) ;
  compute_vertex_permutation(mris_lh_mov, mris_lh_fixed, mht_lh, lh_vertices) ;
  compute_vertex_permutation(mris_rh_mov, mris_rh_fixed, mht_rh, rh_vertices) ;
  mri_cmat_permuted = permute_corrmat(mri_cmat_mov, lh_vertices, rh_vertices, mri_cmat_permuted) ;
  mri_mov_avg = average_within_label(mri_cmat_permuted, area, offset, NULL) ;
  last_error = compute_error(mris_lh_mov, mris_lh_fixed, mht_lh_faces, mri_mov_avg, mri_stats) ;
  if (write_diags)
  {
    char fname[STRLEN] ;
    MRI  *mri_fixed ; 

    mri_fixed = sample_fixed(mris_lh_mov, mris_lh_fixed, mht_lh_faces, mri_stats) ;
    sprintf(fname, "lh.label_avg.%3.3d.mgz",  iter) ;
    printf("writing snapshot to %s\n", fname) ;
    MRIwrite(mri_fixed, fname) ;
    MRIfree(&mri_fixed) ;
  }

  printf("%3.3d: error = %2.3f\n", iter, last_error) ;
  do
  {
    compute_warp_gradient(mris_lh_mov, mris_rh_mov, mris_lh_fixed, mris_rh_fixed, mht_lh_faces, mht_rh_faces, mri_cmat_permuted, mri_stats, mri_mov_avg, mri_cmat_permuted, dt) ;
    MRISapplyGradient(mris_lh_mov, 1) ; MRISprojectOntoSphere(mris_lh_mov, mris_lh_mov, mris_lh_mov->radius) ;
    MRISapplyGradient(mris_rh_mov, 1) ; MRISprojectOntoSphere(mris_rh_mov, mris_rh_mov, mris_rh_mov->radius) ;

    mri_cmat_permuted = permute_corrmat(mri_cmat_mov, lh_vertices, rh_vertices, mri_cmat_permuted) ;
    average_within_label(mri_cmat_permuted, area, offset, mri_mov_avg) ;
    if (write_diags)
    {
      char fname[STRLEN] ;
      MRI  *mri_fixed ; 

      sprintf(fname, "lh.label_avg.%3.3d.mgz",  iter+1) ;
      printf("writing snapshot to %s\n", fname) ;
      mri_fixed = sample_fixed(mris_lh_mov, mris_lh_fixed, mht_lh_faces, mri_stats) ;
      MRIwrite(mri_fixed, fname) ;
      MRIfree(&mri_fixed) ;
    }
    compute_vertex_permutation(mris_lh_mov, mris_lh_fixed, mht_lh, lh_vertices) ;
    compute_vertex_permutation(mris_rh_mov, mris_rh_fixed, mht_rh, rh_vertices) ;
    error = compute_error(mris_lh_mov, mris_lh_fixed, mht_lh_faces, mri_mov_avg, mri_stats) ;
    pct_change = 100*(last_error - error) / last_error ;
    printf("%3.3d: error = %2.3f (%2.3f%%)\n", iter+1, error, pct_change) ;
    last_error = error ;
  } while (pct_change > tol && iter++ < MAX_ITERS) ;
    

  free(lh_vertices) ; free(rh_vertices) ; MHTfree(&mht_lh) ; MHTfree(&mht_rh) ;
  return(NO_ERROR) ;
}
