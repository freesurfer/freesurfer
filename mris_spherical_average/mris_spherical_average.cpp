/**
 * @brief spherical averaging of labels, curvs and vals.
 *
 * apply spherical averaging to various scalars (see Fischl et al, HBM, 1999)
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"

#include "mri.h"
#include "mrisurf.h"
#include "mrisurf_project.h"
#include "mrisurf_vals.h"
#include "mrishash.h"
#include "icosahedron.h"

#include "error.h"
#include "diag.h"
#include "proto.h"
#include "tags.h"
#include "label.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static const char *surf_dir = "surf" ;
static int erode = 0 ;
static int dilate = 0 ;
static float threshold = 0 ;
static int reassign = 0 ;
static int normalize_flag = 0 ;
static int condition_no = 0 ;
static int stat_flag = 0 ;
static char *output_subject_name = NULL ;
static int navgs = 0 ;
static char *ohemi = NULL ;
static char *osurf = NULL ;
static const char *orig_name = "white" ;
static int segment = 0 ;  // not implemented yet
static char *mask_name = NULL ;

static int compute_average_label_area = 0 ;
static int which_ic = 7 ;
static char *sdir = NULL ;
static char *osdir = NULL ;
static double logodds_slope = 0.1 ;

static int spatial_prior_avgs = 0 ;
static char *spatial_prior_fname = NULL ;
static char dir[STRLEN] = "";

int
main(int argc, char *argv[])
{
  char            **av, *out_fname, *surf_name, fname[STRLEN],
                  *hemi, *cp, *data_fname ;
  int             ac, nargs, i, which, nsubjects ;
  double          max_len, mean, sigma ;
  MRI_SURFACE     *mris, *mris_avg ;
  MRIS_HASH_TABLE *mht = NULL ;
  LABEL           *area, *area_avg = NULL ;
  float           average_label_area = 0 ;

  std::string cmdline = getAllInfo(argc, argv, "mris_spherical_average");

  nargs = handleVersionOption(argc, argv, "mris_spherical_average");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (!sdir)
  {
    sdir = getenv("SUBJECTS_DIR") ;
    if (!sdir)
      ErrorExit(ERROR_BADPARM, "%s: no SUBJECTS_DIR in envoronment.\n",
                Progname);
  }
  if (!osdir)
  {
    osdir = sdir ;
  }
  /*
     command line: <which> <fname> <hemi> <spherical surf> <subject 1> ...
                   <output>
  */
  if (argc < 7)
  {
    usage_exit() ;
  }

  which = -1 ;
  if (!stricmp(argv[1], "coords"))
  {
    which = VERTEX_COORDS ;
  }
  else if (!stricmp(argv[1], "vals"))
  {
    which = VERTEX_VALS ;
    if (strlen(dir) == 0)
    {
      strcpy(dir, "label") ;
    }
  }
  else if (!stricmp(argv[1], "area"))
  {
    which = VERTEX_AREA ;
    if (strlen(dir) == 0)
    {
      strcpy(dir, "surf") ;
    }
  }
  else if (!stricmp(argv[1], "curv"))
  {
    which = VERTEX_CURV ;
    if (strlen(dir) == 0)
    {
      strcpy(dir, "surf") ;
    }
  }
  else if (!stricmp(argv[1], "label"))
  {
    which = VERTEX_LABEL ;
    if (strlen(dir) == 0)
    {
      strcpy(dir, "label") ;
    }
  }
  else if (!stricmp(argv[1], "logodds"))
  {
    which = VERTEX_LOGODDS ;
    if (strlen(dir) == 0)
    {
      strcpy(dir, "label") ;
    }
  }
  else
  {
    usage_exit() ;
  }

  data_fname = argv[2] ;
  hemi = argv[3] ;
  if (!ohemi)
  {
    ohemi = hemi ;
  }
  surf_name = argv[4] ;
  if (!osurf)
  {
    osurf = surf_name ;
  }
  out_fname = argv[argc-1] ;

  cp = getenv("FREESURFER_HOME") ;
  if (!cp)
  {
    ErrorExit(ERROR_BADPARM,
              "%s: FREESURFER_HOME not defined in environment",
              Progname);
  }

  sprintf(fname, "%s/lib/bem/ic%d.tri", cp, which_ic) ;
  mris_avg = ICOread(fname) ;
  if (!mris_avg)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not read ico from %s",Progname,fname) ;
  }

  if (which != VERTEX_LABEL)
    MRISclearWhichAndVal2(mris_avg, which) ;

#define FIRST_SUBJECT 5
  for (nsubjects = 0, i = FIRST_SUBJECT ; i < argc-1 ; i++, nsubjects++)
  {
    fprintf(stderr, "processing subject %s...\n", argv[i]) ;
    sprintf(fname, "%s/%s/%s/%s.%s", sdir, argv[i], surf_dir, hemi, surf_name) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    MRISaddCommandLine(mris, cmdline) ;
    if (which == VERTEX_COORDS)
    {
      printf("reading surface coords from %s...\n", orig_name) ;
      if (MRISreadVertexPositions(mris, orig_name) != NO_ERROR)
      {
        ErrorExit(ERROR_BADPARM,
                  "could not read surface positions from %s",
                  orig_name) ;
      }
      MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    }
    // read orig coords in case we need to assign vertices
    if (which == VERTEX_LABEL || which == VERTEX_LOGODDS)
    {
      if (MRISreadOriginalProperties(mris, orig_name) != NO_ERROR)
      {
        ErrorExit(ERROR_BADPARM,
                  "could not read surface positions from %s",
                  orig_name) ;
      }
    }
    MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    if (i == FIRST_SUBJECT)  /* scale the icosahedron up */
    {
      MRISprojectOntoSphere(mris_avg, mris_avg, DEFAULT_RADIUS) ;
      mean = MRIScomputeVertexSpacingStats(mris_avg,
                                           &sigma,
                                           NULL,
                                           &max_len,
                                           NULL,
                                           NULL,
                                           CURRENT_VERTICES);
      if (max_len > mean+3*sigma)
      {
        max_len = mean+3*sigma ;
      }
      mht = MHTcreateVertexTable_Resolution(mris_avg, CURRENT_VERTICES,2*max_len);
    }

    switch (which)
    {
    case VERTEX_VALS:
    {
      char fname[STRLEN] ;
      sprintf(fname,"%s/%s/%s/%s.%s", sdir, argv[i], dir, hemi, data_fname) ;
      if (MRISreadValues(mris, fname) != NO_ERROR)
      {
        ErrorExit(ERROR_BADFILE,
                  "%s: could not read val file %s.\n",
                  Progname, fname);
      }
      if (mask_name)
      {
	LABEL *area ;
	sprintf(fname,"%s/%s/label/%s.%s", sdir, argv[i], hemi, mask_name) ;
	area = LabelRead(NULL, fname) ;
	if (!area)
	  ErrorPrintf(ERROR_BADFILE,"%s: could not read label file %s for %s (%s).\n",
		      Progname, mask_name, argv[i], fname);
	else
	{
	  LabelDilate(area, mris, dilate, CURRENT_VERTICES) ;
	  LabelSetVals(mris, area, 0) ;
	  LabelFree(&area) ;
	}
      }

      MRIScopyValuesToCurvature(mris) ;
      if (threshold > 0)
      {
        MRISbinarizeCurvature(mris, threshold, 0, 1, 0) ;
      }
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize_flag)
      {
        MRISnormalizeCurvature(mris, NORM_MEAN) ;
      }
      MRIScopyCurvatureToValues(mris) ; // MRIScombine will use v->val
      break ;
    }
    case VERTEX_LOGODDS:
      if (i == FIRST_SUBJECT)
      {
        area_avg = LabelAlloc(mris_avg->nvertices, NULL, data_fname) ;
      }
      if (strchr(data_fname, '/') != NULL)
      {
        strcpy(fname, data_fname) ;
      }
      else
      {
        sprintf(fname, "%s/%s/%s/%s", sdir, argv[i], dir, data_fname) ;
      }
      area = LabelRead(NULL, fname) ;
      if (!area)
        ErrorExit(ERROR_BADFILE,"%s: could not read label file %s for %s.\n",
                  Progname, data_fname, argv[i]);
      if (reassign)
      {
        LabelUnassign(area) ;
      }
      LabelFillUnassignedVertices(mris, area, ORIG_VERTICES) ;
      if (argc-1-FIRST_SUBJECT > 1)
      {
        LabelSetStat(area, 1) ;
      }
      else
      {
        printf("only %d subject - copying statistics...\n",
               argc-1-FIRST_SUBJECT );
      }
      if (navgs > 0)
      {
        int i ;
        LabelMarkStats(area, mris) ;
        LabelFree(&area) ;
        MRISaverageMarkedStats(mris, navgs) ;
        area = LabelFromMarkedSurface(mris) ;
        for (i = 0 ; i < area->n_points ; i++)
          if (FZERO(area->lv[i].stat))
          {
            DiagBreak() ;
          }
      }
      MRISlogOdds(mris, area, logodds_slope) ;
      break ;
    case VERTEX_LABEL:
      if (i == FIRST_SUBJECT)
      {
        area_avg = LabelAlloc(mris_avg->nvertices, NULL, data_fname) ;
      }
      if (strchr(data_fname, '/') != NULL)
      {
        strcpy(fname, data_fname) ;
      }
      else
      {
        sprintf(fname, "%s/%s/%s/%s", sdir, argv[i], dir, data_fname) ;
      }
      area = LabelRead(NULL, fname) ;
      if (!area)
        ErrorExit(ERROR_BADFILE,"%s: could not read label file %s for %s.\n",
                  Progname, data_fname, argv[i]);
      if (reassign)
      {
        LabelUnassign(area) ;
      }
      LabelFillUnassignedVertices(mris, area, ORIG_VERTICES) ;
      if (argc-1-FIRST_SUBJECT > 1)
      {
        LabelSetStat(area, 1) ;
      }
      else
      {
        printf("only %d subject - copying statistics...\n",
               argc-1-FIRST_SUBJECT );
      }
      if (navgs > 0)
      {
        int i ;
        LabelMarkStats(area, mris) ;
        LabelFree(&area) ;
        MRISaverageMarkedStats(mris, navgs) ;
        area = LabelFromMarkedSurface(mris) ;
        for (i = 0 ; i < area->n_points ; i++)
          if (FZERO(area->lv[i].stat))
          {
            DiagBreak() ;
          }
      }
      average_label_area += LabelArea(area, mris) ;
      area_avg = LabelSphericalCombine(mris, area, mht, mris_avg, area_avg) ;
      break ;
    case VERTEX_CURVATURE:
      if (MRISreadCurvatureFile(mris, data_fname) != NO_ERROR)
        ErrorExit(ERROR_BADFILE,"%s: could not read curvature file %s.\n",
                  Progname, data_fname);
      MRISaverageCurvatures(mris, navgs) ;
      if (normalize_flag)
      {
        MRISnormalizeCurvature(mris, NORM_MEAN) ;
      }
      break ;
    case VERTEX_AREA:
      if (MRISreadOriginalProperties(mris, data_fname) != NO_ERROR)
        ErrorExit(ERROR_BADFILE,"%s: could not read surface file %s.\n",
                  Progname, data_fname);
#if 0
      fprintf(stderr, "total orig area = %2.1f 10^3 mm\n",
              mris->orig_area/1000) ;
#endif
      break ;
    default:

      break ;
    }
    if (which != VERTEX_LABEL)
    {
      MRIScombine(mris, mris_avg, mht, which) ;
    }
    if (i < argc-2)
    {
      MRISfree(&mris) ;
    }
  }
  if (which != VERTEX_LABEL)
  {
    MRISnormalize(mris_avg, nsubjects, which) ;
  }
  else
  {
    average_label_area /= nsubjects ;
    printf("average label area =  %2.2fmm^2\n", average_label_area) ;
  }
  if (mht)
  {
    MHTfree(&mht) ;
  }

  if (output_subject_name)
  {
    sprintf(fname, "%s/%s/%s/%s.%s", osdir,output_subject_name,surf_dir,ohemi,osurf);
    fprintf(stderr, "reading output surface %s...\n", fname) ;
    MRISfree(&mris) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    if (MRISreadOriginalProperties(mris, orig_name) != NO_ERROR)
    {
      ErrorExit(ERROR_BADPARM,
                "could not read surface positions from %s",
                orig_name) ;
    }
    MRISprojectOntoSphere(mris, mris, DEFAULT_RADIUS) ;
  }

  if (which != VERTEX_LABEL)
    MRISclearWhichAndVal2(mris, which) ;
  mean = MRIScomputeVertexSpacingStats(mris, &sigma, NULL, &max_len, NULL,NULL,
                                       CURRENT_VERTICES);
  if (max_len > mean+3*sigma)
  {
    max_len = mean+3*sigma ;
  }
  mht = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 2*max_len);
  if (which != VERTEX_LABEL)
  {
    MRISsphericalCopy(mris_avg, mris, mht, which) ;
  }
  else
  {
    LabelFree(&area) ;
    area = LabelSphericalCombine(mris_avg, area_avg, mht, mris, NULL) ;
    LabelRemoveDuplicates(area) ;
  }
  MHTfree(&mht) ;
  if (which == VERTEX_AREA)
  {
    MRISorigAreaToCurv(mris) ;
  }
  if (stat_flag)    /* write out summary statistics files */
  {
    int    vno ;
    VERTEX *v ;
    float  dof ;
    FILE   *fp ;

    sprintf(fname, "%s/sigavg%d-%s.w", out_fname, condition_no, ohemi);
    fprintf(stderr, "writing output means to %s\n", fname) ;
    MRISwriteCurvatureToWFile(mris, fname) ;

    /* change variances to squared standard errors */
    dof = nsubjects ;
    if (!FZERO(dof)) for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag)
        {
          continue ;
        }
        v->curv = v->val2 / dof ;   /* turn it into a standard error */
      }

    sprintf(fname, "%s/sigvar%d-%s.w", out_fname, condition_no, ohemi);
    fprintf(stderr, "writing output variances to %s\n", fname) ;
    MRISwriteCurvatureToWFile(mris, fname) ;

    /* write out dof file */
    sprintf(fname, "%s/sigavg%d.dof", out_fname, condition_no) ;
    fp = fopen(fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open dof file %s\n",
                Progname,fname);
    fprintf(stderr, "writing dof file %s\n", fname) ;
    fprintf(fp, "%d\n", (int)dof) ;
    fclose(fp) ;
  }
  else
  {
    if (Gdiag & DIAG_SHOW)
    {
      fprintf(stderr,"writing blurred pattern to surface to %s\n",out_fname);
    }
    switch (which)
    {
    case VERTEX_LABEL:
      LabelToOriginal(area, mris) ;  // use orig coords in output surface
      printf("writing label with %d points to %s...\n", area->n_points,
             out_fname) ;
      if (normalize_flag)
      {
        LabelNormalizeStats(area, (float)nsubjects) ;
      }
      if (spatial_prior_fname)
      {
        MRISclearMarks(mris) ;
        LabelMarkStats(area, mris) ;
        MRIScopyStatsToValues(mris) ;
        MRISaverageVals(mris, spatial_prior_avgs) ;
        MRISwriteValues(mris, spatial_prior_fname) ;
      }
       if (threshold > 0)
      {
        LabelThreshold(area, threshold) ;
      }
      if (dilate > 0 || erode > 0)
      {
        MRISsaveVertexPositions(mris, TMP_VERTICES) ;
        MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
      }
      if (dilate > 0)
      {
        LabelDilate(area, mris, dilate, CURRENT_VERTICES) ;
      }
      if (erode > 0)
      {
        LabelErode(area, mris, erode) ;
      }
      if (dilate > 0 || erode > 0)
      {
        MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
      }
      if (compute_average_label_area)
      {
	LABEL *area_saved ;
	double thresh, best_thresh, surface_area, best_area ;

	area_saved = LabelCopy(area, NULL) ;
	best_area = LabelArea(area_saved, mris) ; best_thresh = 0.0 ;
	for (thresh = 0 ; thresh <= 1 ; thresh += .01)
	{
	  LabelThreshold(area_saved, thresh) ;
	  surface_area = LabelArea(area_saved, mris) ;
	  if (fabs(surface_area - average_label_area) < fabs(best_area - average_label_area))
	  {
	    best_area = surface_area ;
	    best_thresh = thresh ;
	  }
	}
	LabelThreshold(area, best_thresh) ;
	printf("threshold that best approximates average area of %2.2f mm^2 is %2.2f (%2.2f mm^2)\n", average_label_area, best_thresh, LabelArea(area,mris));
      }

      LabelWrite(area, out_fname) ;
      break ;
    case VERTEX_LOGODDS:
    case VERTEX_VALS:
      //      MRIScopyCurvatureToValues(mris) ;
      MRISwriteValues(mris, out_fname) ;
      break ;
    case VERTEX_AREA:
    case VERTEX_CURV:
      MRISwriteCurvature(mris, out_fname) ;
      break ;
    case VERTEX_COORDS:
      MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
      MRISwrite(mris, out_fname) ;
      break ;
    default:
      break ;
    }
  }

  MRISfree(&mris) ;
  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "ohemi"))
  {
    ohemi = argv[2] ;
    fprintf(stderr, "output hemisphere = %s\n", ohemi) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "surf_dir"))
  {
    surf_dir = argv[2] ;
    fprintf(stderr, "using %s instead surf directory\n", surf_dir) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "lslope"))
  {
    logodds_slope = atof(argv[2]) ;
    printf("using logodds slope = %2.2f\n", logodds_slope) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "ic"))
  {
    which_ic = atoi(argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "average_area"))
  {
    compute_average_label_area = 1 ;
    printf("computing threshold that yields surface area of average label closest to average of individual surface areas\n") ;
  }
  else if (!stricmp(option, "sdir"))
  {
    sdir = argv[2] ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "dir"))
  {
    strcpy(dir, argv[2]) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "mask"))
  {
    mask_name = argv[2] ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "prior"))
  {
    spatial_prior_avgs = atoi(argv[2]) ;
    spatial_prior_fname = argv[3] ;
    nargs = 2 ;
    printf("blurring label priors and writing output to %s\n",
           spatial_prior_fname) ;
  }
  else if (!stricmp(option, "segment"))
  {
    segment = 1 ;
    printf("segmenting input labels and only using largest connected component\n") ;
  }
  else if (!stricmp(option, "osdir"))
  {
    osdir = argv[2] ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "orig"))
  {
    orig_name = argv[2] ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "erode"))
  {
    erode = atoi(argv[2]) ;
    nargs = 1 ;
    printf("eroding label %d times before writing\n", erode) ;
  }
  else if (!stricmp(option, "dilate"))
  {
    dilate = atoi(argv[2]) ;
    nargs = 1 ;
    printf("dilating label %d times before writing\n", dilate) ;
  }
  else if (!stricmp(option, "close"))
  {
    dilate = atoi(argv[2]) ;
    erode = dilate ;
    nargs = 1 ;
    printf("closing label %d times before writing\n", dilate) ;
  }
  else if (!stricmp(option, "reassign"))
  {
    reassign = 1 ;
    printf("recomputing label vertex assignments\n") ;
  }
  else if (!stricmp(option, "osurf"))
  {
    osurf = argv[2] ;
    fprintf(stderr, "output surface = %s\n", osurf) ;
    nargs = 1 ;
  }
  else switch (toupper(*option))
    {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      printf("debugging vertex %d\n", Gdiag_no) ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      fprintf(stderr, "blurring thickness for %d iterations\n",navgs);
      nargs = 1 ;
      break ;
    case 'T':
      threshold = atof(argv[2]) ;
      printf("thresholding label stat at %2.3f before writing\n", threshold) ;
      nargs = 1 ;
      break ;
    case 'O':
      output_subject_name = argv[2] ;
      nargs = 1 ;
      fprintf(stderr, "painting output onto subject %s.\n",
              output_subject_name);
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'S':   /* write out stats */
      stat_flag = 1 ;
      condition_no = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "writing out summary statistics as condition %d\n",
              condition_no) ;
      break ;
    case 'N':
      normalize_flag = 1 ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr,
          "usage: %s [option] <which> <fname> <hemi> <spherical surf> "
          "<subject 1> ... <output>\n", Progname) ;
  fprintf(stderr, "where which is one of\n"
          "\tcoords\n"
          "\tlabel\n"
          "\tvals\n"
          "\tcurv\n"
          "\tarea\n") ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr,
          "\nThis program will add a template into an average surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-segment        only use largest connected component of label\n");
  fprintf(stderr, "-n              normalize output so it can be interpreted as a probability\n") ;
  fprintf(stderr, "-orig  <name>   use <name> as original surface position default=orig\n");
  fprintf(stderr, "-o  <output subject name>   use <output subject> as the space to write the results in instead of the last subject given\n");
  fprintf(stderr, "-osdir  <output subject dir>   use <output subject dir> as the subjects dir for the output subject'\n");
  fprintf(stderr, "-sdir  <subjects dir>   use <subject dir> as the subjects dir'\n");
  fprintf(stderr, "-average_area   compute threshold for label that will give the average label apporximately the average surface area\n") ;
  fprintf(stderr, "-s <cond #>     generate summary statistics and write\n"
          "                them into sigavg<cond #>-<hemi>.w and\n"
          "                sigvar<cond #>-<hemi>.w.\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

