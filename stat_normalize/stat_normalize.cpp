/*
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "stats.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static float resolution = 8.0f ;
static float fov = 256.0f ;
static int   coordinate_system = TALAIRACH_COORDS ;
static char  *hemi ;
static char  *surf_name ;   /* used if in surface-based coordinates */

int
main(int argc, char *argv[]) {
  char        *in_prefix, *out_prefix, out_fname[100], name[100],
    path[100], fname[100], *cp, subjects_dir[100] ;
  const char *coord_name;
  int         n, nargs, ino, event ;
  SV          *sv, *sv_avg = NULL ;
  MRI_SURFACE *mris ;

  nargs = handleVersionOption(argc, argv, "stat_normalize");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  /* print out command-line */
  for (n=0; n < argc; n++) printf("%s ",argv[n]);
  printf("\n");

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    print_help() ;
  cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM,
              "%s: SUBJECTS_DIR not specified in the environment", Progname) ;
  strcpy(subjects_dir, cp) ;

  switch (coordinate_system) {
  default:
  case TALAIRACH_COORDS:
    coord_name = "Talairach" ;
    break ;
  case SPHERICAL_COORDS:
    coord_name = "Spherical" ;
    break ;
  case ELLIPSOID_COORDS:
    coord_name = "Ellipsoid" ;
    break ;
  }

  out_prefix = argv[argc-1] ;

  for (ino = 1 ; ino < argc-1 ; ino++) {
    /* for each path/prefix specified, go through all slices */
    in_prefix = argv[ino] ;
    fprintf(stderr, "reading stat volume %s.\n", in_prefix) ;
    FileNamePath(in_prefix, path) ;
    FileNameOnly(in_prefix, name) ;

    sv = StatReadVolume2(in_prefix) ;
    if (!sv)
      ErrorExit(ERROR_NOFILE, "%s: could not read stat files %s",
                Progname, in_prefix) ;

    if (!sv_avg)
      sv_avg = StatAllocStructuralVolume(sv, fov, resolution, coord_name) ;
    switch (coordinate_system) {
    default:
    case TALAIRACH_COORDS:
      StatAccumulateTalairachVolume(sv_avg, sv) ;
      break ;
    case SPHERICAL_COORDS:
    case ELLIPSOID_COORDS:
      int req = snprintf(fname, 100, "%s/%s/surf/%s.orig",
			 subjects_dir, sv->reg->name, hemi) ;   
      if( req >= 100 ) {
	std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      fprintf(stderr, "reading surface %s\n", fname) ;
      mris = MRISread(fname) ;
      if (!mris)
        ErrorExit(ERROR_NOFILE,
                  "%s: could not read surface %s",Progname,fname);

      /* load the coordinates in the canonical
         surface space for this surf. */
      if (MRISreadCanonicalCoordinates(mris, surf_name) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read canonical surface %s.",
                  Progname, surf_name) ;
      fprintf(stderr, "adding %s to the average\n", sv->reg->name) ;
      StatAccumulateSurfaceVolume(sv_avg, sv, mris) ;
      MRISfree(&mris) ;
      break ;
    }

    StatFree(&sv) ;
  }


  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    for (event = 0 ; event < sv_avg->nevents ; event++) {

      sprintf(out_fname, "avg%d.mgh", event) ;
      MRIwrite(sv_avg->mri_avgs[event], out_fname) ;
      sprintf(out_fname, "std%d.mgh", event) ;
      MRIwrite(sv_avg->mri_stds[event], out_fname) ;
    }

  fprintf(stderr, "writing average volume to %s.\n", out_prefix) ;
  StatWriteVolume(sv_avg, out_prefix) ;
  StatFree(&sv_avg) ;

  exit(0) ;
  return(0) ;  /* for ansi */
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
  if (!stricmp(option, "-help"))         print_help() ;
  else if (!stricmp(option, "-version")) print_version() ;
  else switch (toupper(*option)) {
    case 'E':
      coordinate_system = ELLIPSOID_COORDS ;
      hemi = argv[2] ;
      surf_name = argv[3] ;
      nargs = 2 ;
      break ;
    case 'S':
      coordinate_system = SPHERICAL_COORDS ;
      hemi = argv[2] ;
      surf_name = argv[3] ;
      nargs = 2 ;
      break ;
    case '?':
    case 'U':
      usage_exit() ;
      exit(1) ;
      break ;
    case 'R':
      sscanf(argv[2], "%f", &resolution) ;
      printf("INFO: settting resolution to %f\n",resolution);
      nargs = 1 ;
      break ;
    case 'X':
      stats_talxfm = argv[2];
      printf("INFO: using %s\n",stats_talxfm);
      nargs = 1 ;
      break ;
    case 'I':
      stats_fixxfm = 1;
      printf("INFO: flag giving to devolve tal xfm\n");
      break ;
    case 'C':
      statnorm_float2int = float2int_code(argv[2]);
      if (statnorm_float2int < 0) {
        printf("ERROR: float2int code %s unrecognized\n",argv[2]);
        exit(1);
      }
      printf("INFO: using %s float2int\n",argv[2]);
      nargs = 1 ;
      break ;
    case 'f':
      sscanf(argv[2], "%f", &fov) ;
      nargs = 1 ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf
  (stderr,
   "usage: %s [options] <input sv prefix> <output sv prefix>\n",
   Progname) ;
  fprintf(stderr, "options are:\n") ;
  fprintf
  (stderr,
   "\t-r <resolution>            - set output resolution (def=8mm)\n") ;
  fprintf
  (stderr,
   "\t-f <field of view>         - set output field of view (def=256)\n");
  fprintf
  (stderr,
   "\t-S <hemisphere> <surface>  - average in spherical coordinates\n");
  fprintf
  (stderr,
   "\t-x xfmfile - use subjid/mri/transforms/xfmfile instead of\n");
  fprintf
  (stderr,
   "\t-i - fix xfm for non-zero center of orig volume\n");
  fprintf
  (stderr,
   "\t-c float2int - <tkregister>, round\n");
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will convert average a sequence of\n") ;
  fprintf(stderr,"volume-based statistics in Talairach space:\n\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

