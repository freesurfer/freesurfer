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
#include "macros.h"
#include "utils.h"

static char vcid[]="$Id: mris_errors.c,v 1.3 1997/12/10 23:23:11 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
int MRISareaErrors(MRI_SURFACE *mris) ;

char *Progname ;
static MRI_SURFACE  *mris ;

static int patch_flag = 0 ;
static int write_flag = 0 ;

int
main(int argc, char *argv[])
{
  char         *cp, **av, *in_fname, *out_fname, fname[100], path[100], 
               name[100], hemi[100] ;
  int          ac, nargs ;

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

  if (argc < 2)
    usage_exit() ;

  in_fname = argv[1] ;
#if 0
  out_fname = argv[2] ;
  cp = strrchr(out_fname, '.') ;
#endif

  if (patch_flag)   /* read in orig surface before reading in patch */
  {
    FileNamePath(in_fname, path) ;
    cp = strrchr(in_fname, '.') ;
    if (cp)
    {
      strncpy(hemi, cp-2, 2) ;
      hemi[2] = 0 ;
    }
    else
      strcpy(hemi, "lh") ;
    sprintf(fname, "%s/%s.smoothwm", path, hemi) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    MRIScomputeMetricProperties(mris) ;
    MRISstoreMetricProperties(mris) ;
    FileNameOnly(in_fname, name) ;
    if (MRISreadPatch(mris, name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read patch file %s",
                Progname, name) ;
  }
  else
  {
    mris = MRISread(in_fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, in_fname) ;

    MRISstoreCurrentPositions(mris) ;
    MRISreadVertexPositions(mris, "smoothwm") ;
    MRIScomputeMetricProperties(mris) ;
    MRISstoreMetricProperties(mris) ;
    MRISrestoreOldPositions(mris) ;
  }

  MRIScomputeMetricProperties(mris) ;  /* recompute areas and normals */
  if (write_flag)
  {
    MRISareaErrors(mris) ;
    MRISangleErrors(mris) ;
  }

  if (1)
  {
    sprintf(fname, "%s.area_error", in_fname) ;
    printf("writing area errors to %s\n", fname) ;
    MRISwriteAreaError(mris, fname) ;
    sprintf(fname, "%s.angle_error", in_fname) ;
    printf("writing angle errors to %s\n", fname) ;
    MRISwriteAngleError(mris, fname) ;
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
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option))
  {
  case 'W':
    write_flag = 1 ;
    break ;
  case 'P':
    patch_flag = 1 ;
    nargs = 0 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
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
          "usage: %s [options] <input image file>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program will unfold an MRI on the surface of an ellipsoid.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

int
MRISareaErrors(MRI_SURFACE *mris)
{
  int      fno, tno, max_f = -1 ;
  FACE     *face ;
  float    ferror, max_ferror, total_error, total_sq_error,
           error, mean_error, std_error, pct_error, n ;

  MRISupdateEllipsoidSurface(mris) ;
  total_error = total_sq_error = max_ferror = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    ferror = 0.0f ;
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      error = face->area[tno] - face->orig_area[tno] ;
      pct_error = error / face->orig_area[tno] * 100.0f ;
      printf("%d %2.3f %2.3f %2.3f %2.1f\n", fno, face->orig_area[tno], 
             face->area[tno], error, pct_error) ;
      total_sq_error += (error * error) ;
      ferror += fabs(error) ;
      total_error += error ;
    }
    if (ferror >= max_ferror)
    {
      max_ferror = ferror ;
      max_f = fno ;
    }
  }

  n = (float)(2*mris->nfaces) ;
  mean_error = total_error / n ;
  std_error = sqrt(total_sq_error / (float)n - mean_error*mean_error) ;
  fprintf(stderr, "max error occurs at %d, error = %2.3f\n",max_f, max_ferror);
  fprintf(stderr, "mean error = %2.3f, std = %2.3f\n", mean_error, std_error);
  return(NO_ERROR) ;
}


int
MRISangleErrors(MRI_SURFACE *mris)
{
  int      fno, tno, max_f = -1, ano ;
  FACE     *face ;
  FILE     *fp ;
  float    ferror, max_ferror, total_error, total_sq_error,
           error, mean_error, std_error, pct_error, n ;

  fp = fopen("angle.err", "w") ;
  MRISupdateEllipsoidSurface(mris) ;
  total_error = total_sq_error = max_ferror = 0.0f ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    ferror = 0.0f ;
    for (tno = 0 ; tno < TRIANGLES_PER_FACE ; tno++)
    {
      for (ano = 0 ; ano < ANGLES_PER_TRIANGLE ; ano++)
      {
        error = deltaAngle(face->angle[tno][ano], face->orig_angle[tno][ano]);
        pct_error = error / face->orig_angle[tno][ano] * 100.0f ;
        fprintf(fp, "%d %2.3f %2.3f %2.3f %2.1f\n", fno, 
                DEGREES(face->orig_angle[tno][ano]), 
                DEGREES(face->angle[tno][ano]), DEGREES(error), pct_error) ;
        total_sq_error += (error * error) ;
        ferror += fabs(error) ;
        total_error += error ;
      }
    }
    if (ferror >= max_ferror)
    {
      max_ferror = ferror ;
      max_f = fno ;
    }
  }

  n = (float)(2*mris->nfaces) ;
  mean_error = total_error / n ;
  std_error = sqrt(total_sq_error / (float)n - mean_error*mean_error) ;
  fprintf(stderr, "max angle error occurs at %d, error = %2.3f\n",
          max_f, DEGREES(max_ferror));
  fprintf(stderr, "mean angle error = %2.3f, std = %2.3f\n", 
          DEGREES(mean_error), DEGREES(std_error));
  fclose(fp) ;
  return(NO_ERROR) ;
}

