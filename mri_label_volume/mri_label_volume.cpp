/**
 * @brief program for computing the volume of a set of labeled voxels
 *
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
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "timer.h"
#include "cma.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static const char *log_fname = NULL ;
static void usage_exit(int code) ;

static int quiet = 0 ;
static int spread_sheet = 0 ;
static int partial_volume = 0 ;
static MRI *mri_vals ;  /* for use in partial volume calculation */
static int in_label = -1 ;
static int out_label = -1 ;
static char  *subject_name = NULL ;
static int all_flag = 0 ;
static int compute_pct = 0 ;
static char *brain_fname = NULL ;
static char *icv_fname = NULL ;
#define MAX_COLS 10000
static char *col_strings[MAX_COLS] ;
static int ncols  = 0 ;
static double atlas_icv = -1 ;


int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs, msec, minutes, label, volume, seconds, i ;
  Timer start ;
  MRI    *mri ;
  FILE   *log_fp ;
  double  vox_volume, brain_volume ;

  nargs = handleVersionOption(argc, argv, "mri_label_volume");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if ((all_flag && argc < 2) || (all_flag == 0 && argc < 3))
    usage_exit(1) ;

  // Load the segmentation
  mri = MRIread(argv[1]) ;
  if (!mri)
    ErrorExit
    (ERROR_NOFILE, "%s: could not read volume from %s", Progname,argv[1]) ;

  // Volume of a single voxel
  vox_volume = mri->xsize * mri->ysize * mri->zsize ;

  // Not sure why you would want to do this, but here it is
  if (in_label >= 0) MRIreplaceValues(mri, mri, in_label, out_label) ;

  // This looks like it is for debugging
  if (all_flag) {
    int   nvox ;
    float volume ;
    nvox = MRItotalVoxelsOn(mri, WM_MIN_VAL) ;
    volume = nvox * mri->xsize * mri->ysize * mri->zsize ;
    printf("total volume = %d voxels, %2.1f mm^3\n", nvox, volume) ;
    exit(0) ;
  }

  // Compute the volume of the brain in one of
  // four ways (or don't use brain volume)
  if (brain_fname) {
    // (1) Load the brain volume and count voxels over a threshold
    MRI *mri_brain = MRIread(brain_fname) ;
    if (mri_brain == NULL)
      ErrorExit
      (ERROR_BADPARM,
       "%s: could not read brain volume from %s\n", Progname,brain_fname) ;

    brain_volume = (double)MRItotalVoxelsOn(mri_brain, WM_MIN_VAL) ;
    MRIfree(&mri_brain) ;
    brain_volume *= (mri->xsize * mri->ysize * mri->zsize) ;
  } else if (atlas_icv > 0) {
    // (2) Use ICV from atlas. See below when option is processed.
    brain_volume = atlas_icv;
  } else if (compute_pct) {
    // (3) Count voxels in labels that are brain
    for (brain_volume = 0.0, label = 0 ; label <= MAX_CMA_LABEL ; label++) {
      if (!IS_BRAIN(label)) continue ;
      brain_volume += (double)MRIvoxelsInLabel(mri, label) ;
    }
    brain_volume *= (mri->xsize * mri->ysize * mri->zsize) ;
  } else if (icv_fname) {
    // (4) Use ICV supplied in a file
    FILE *fp ;
    fp = fopen(icv_fname, "r") ;
    if (fp == NULL)
      ErrorExit
      (ERROR_NOFILE,
       "%s: could not open ICV file %s\n", Progname, icv_fname) ;
    if (fscanf(fp, "%lf", &brain_volume) != 1)
      ErrorExit
      (ERROR_NOFILE,
       "%s: could not read ICV from %s\n", Progname, icv_fname) ;
    fclose(fp) ;
    printf("using intra-cranial volume = %2.1f\n", brain_volume) ;
  } else {
    // (5) Just use brain_volume=1 (ie, don't try to take it into account)
    brain_volume = 1.0 ;
  }

  // For spread sheet, print first col as subj name.
  // The next col is the brain volume.
  // Not sure what col_strings are.
  if (spread_sheet) {
    int i ;

    if (log_fname == NULL) log_fname = "area_volumes.log" ;
    log_fp = fopen(log_fname, "a+") ;

    fprintf(log_fp, "%s  ", subject_name) ;
    if (icv_fname ||
        compute_pct ||
        atlas_icv > 0) fprintf(log_fp, "%f ", brain_volume)  ;

    for (i = 0 ;i < ncols ; i++) fprintf(log_fp, "%s ", col_strings[i])  ;
    fclose(log_fp) ;
  }

  // The rest of the args are label indices
  // (or "brain" means to use all labels in brain)
  for (i = 2 ; i < argc ; i++) {

    // choke and die if this argument is an option...
    // this prog doesnt support args at the end of the command
    if (ISOPTION(*argv[i])) {
      fprintf(stderr,"ERROR: options not allowed at end of command:\n");
      usage_exit(1) ;
    }

    // Compute the label volume
    if (stricmp(argv[i], "brain") == 0) {
      // All labels in brain
      volume = 0 ;
      for (label = 0 ; label <= MAX_CMA_LABEL ; label++) {
        if (!IS_BRAIN(label)) continue ;
        if (partial_volume)
          volume += MRIvoxelsInLabelWithPartialVolumeEffects
            (mri, mri_vals, label, NULL, NULL) ;
        else
          volume += MRIvoxelsInLabel(mri, label) ;
      }
      label = -1 ;
    } else {
      // Label-by-label
      label = atoi(argv[i]) ;
      if (quiet == 0)
        printf("processing label %d...\n", label) ;

      if (partial_volume)
        volume = MRIvoxelsInLabelWithPartialVolumeEffects
          (mri, mri_vals, label, NULL, NULL) ;
      else
        volume = MRIvoxelsInLabel(mri, label) ;
    }

    // Open the logfile for appending
    if (log_fname) {
      char fname[STRLEN] ;

      sprintf(fname, log_fname, label) ;
      printf("logging to %s...\n", fname) ;
      log_fp = fopen(fname, "a+") ;
      if (!log_fp)ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                              Progname, fname) ;
    } else log_fp = NULL ;

    // print percent volume, or ...
    if (compute_pct || icv_fname || atlas_icv > 0) {
      printf("%d voxels (%2.1f mm^3) in label %d, "
             "%%%2.6f of %s volume (%2.0f)\n",
             volume, volume*vox_volume,label,
             100.0*(float)volume/(float)brain_volume,
             atlas_icv > 0 ? "eTIV" : "brain",
             brain_volume) ;
      if (log_fp) {
        if (spread_sheet)
          fprintf(log_fp,"%2.6f ",
                  100.0*(float)volume*vox_volume/(float)brain_volume) ;
        else
          fprintf(log_fp,"%2.6f\n",
                  100.0*(float)volume*vox_volume/(float)brain_volume) ;

        fclose(log_fp) ;
      }
    }
    // Print actual volume (instead of percent)
    else {
      if (quiet == 0)
        printf("%d (%2.1f mm^3) voxels in label %d\n",
               volume,volume*vox_volume, label) ;
      else
        printf("%2.1f\n", volume*vox_volume) ;
      if (log_fp) {
        if (spread_sheet)
          fprintf(log_fp,"%2.6f ", (float)volume*vox_volume) ;
        else
          fprintf(log_fp,"%2.1f\n", vox_volume*(float)volume) ;
        fclose(log_fp) ;
      }
    }

  } // End loop over labels

  // For spread sheet, add a newline
  if (spread_sheet) {
    log_fp = fopen(log_fname, "a+") ;
    fprintf(log_fp,  "\n") ;
    fclose(log_fp) ;
  }

  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;

  if (DIAG_VERBOSE_ON)
    fprintf(stderr, "overlap calculation took %d minutes and %d seconds.\n",
            minutes, seconds) ;

  exit(0) ;
  return(0) ;
}


/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "ICV")) {
    icv_fname = argv[2] ;
    printf("reading ICV from %s\n", icv_fname) ;
    nargs = 1 ;
  } else if (!stricmp(option, "PV")) {
    partial_volume = 1 ;
    nargs = 1 ;
    mri_vals = MRIread(argv[2]) ;
    if (mri_vals == NULL)
      ErrorExit
        (ERROR_NOFILE,
         "%s: could not read intensity volume %s", Progname, argv[3]) ;
    printf("including partial volume effects in calculations\n") ;
  } else if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
  } else if (!stricmp(option, "atlas_icv") ||
             !stricmp(option, "eTIV") ||
             !stricmp(option, "eTIV_matdat")) {
    double atlas_det; // filled-in by MRIestimateTIV
    double eTIV_scale_factor = (double)atoi(argv[3]);
    atlas_icv = MRIestimateTIV(argv[2],eTIV_scale_factor,&atlas_det);
    printf("using eTIV from atlas transform of %2.0f cm^3\n",
           atlas_icv/(10*10*10)) ;
    nargs = 2 ;

    if (!stricmp(option, "eTIV_matdat")) {
      /* create matlab-readable data file (append to
         existing file named "det_eTIV_matdat.m") */
      FILE *dat_fp = fopen("det_eTIV_matdat.m", "a+") ;
      if (dat_fp) {
        fprintf(dat_fp,
                "eTIV_scale_factor = %f;\n", eTIV_scale_factor);
        fprintf(dat_fp,
                "det_eTIV(length(det_eTIV)+1).det = %f;\n", atlas_det);
        fprintf(dat_fp,
                "det_eTIV(length(det_eTIV)).eTIV = %f;\n", atlas_icv);
        fprintf(dat_fp,
                "det_eTIV(length(det_eTIV)).id = '%s';\n\n\n", argv[4]);
        fclose(dat_fp);
      }
      nargs = 3 ;
    }
  } else switch (toupper(*option)) {
  case 'C':
    if (ncols >=  MAX_COLS)
      ErrorExit
        (ERROR_NOMEMORY,
         "%s: too many columns specified (max=%d)\n", Progname, ncols) ;
    col_strings[ncols++] = argv[2] ;
    nargs = 1 ;
    break ;
  case 'Q':
    quiet = 1 ;
    break;
  case 'S':
    spread_sheet = 1  ;
    subject_name = argv[2] ;
    nargs = 1 ;
    break  ;
  case 'A':
    all_flag = 1 ;
    printf("computing volume of all non-zero voxels\n") ;
    break ;
  case 'T':
    in_label = atoi(argv[2]) ;
    out_label = atoi(argv[3]) ;
    nargs = 2 ;
    printf("translating label %d to label %d\n", in_label, out_label) ;
    break ;
  case 'B':
    brain_fname = argv[2] ;
    compute_pct = 1 ;
    nargs = 1 ;
    printf("reading brain volume from %s...\n", brain_fname) ;
    break ;
  case 'L':
    log_fname = argv[2] ;
    nargs = 1 ;
    /*    fprintf(stderr, "logging results to %s\n", log_fname) ;*/
    break ;
  case 'P':
    compute_pct = 1 ;
    break ;
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}


/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s [options] <volume> <label 1> <label 2> ...\n", Progname) ;
  printf("valid options are:\n") ;
  printf("  -pv <fname>   - compute partial volume effects "
         "using intensity volume <fname>\n") ;
  printf("  -icv <fname>  - normalize by the intracranial "
         "volume in <fname>\n") ;
  printf("  -s <subject>  - output in spreadsheet mode, "
         "including <subject> name in file\n") ;
  printf("  -a            - compute volume of all non-zero "
         "voxels (e.g. for computing brain volume)\n") ;
  printf("  -t <in> <out> - replace label <in> with label <out>. "
         "Useful for compute e.g. whole hippo vol\n") ;
  printf("  -b <brain vol>- compute the brain volume from "
         "<brain vol> and normalize by it\n") ;
  printf("  -p            - compute volume as a %% of all non-zero labels\n") ;
  printf("  -l <fname>    - log results to <fname>\n") ;
  printf("  -atlas_icv <fname> <scalefactor> - specify LTA or XFM atlas "
         "transform file\n  and scale factor to "
         "use for ICV correction (c.f. Buckner et al. 2004)\n");
  printf("  -eTIV <fname> <scalefactor> - same as -atlas_icv\n");
  printf("  -eTIV_matdat <fname> <scalefactor> <subject> - same as -eTIV,\n"
         "  and generate matlab data appending <subject> to structure\n");
  exit(code) ;
}
