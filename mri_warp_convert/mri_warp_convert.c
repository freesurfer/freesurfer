/**
 * @file  mri_warp_convert.c
 * @brief convert an FSL volume warp into the M3Z format.
 *
 * simple wrapper around GCAMremoveSingularitiesAndReadWarpFromMRI().
 *
 */
/*
 * Original Author: jonathan polimeni
 * CVS Revision Info:
 *    $Author: jonp $
 *    $Date: 2012/02/13 23:05:52 $
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "gcamorph.h"
#include "mri.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
static void usage_exit(int code) ;

char *basename( char *path );

char*
basename (char* path)
{
  char *ptr = strrchr (path, '/');
  return ptr ? ptr + 1 : (char*)path;
}


int
main(int argc, char *argv[])
{
  char         **av, *out_name ;
  int          ac, nargs ;
  int          msec, minutes, seconds ;
  struct timeb start ;
  GCA_MORPH    *gcam ;
  MRI          *mri = NULL ;
  MATRIX       *m;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_warp_convert.c,v 1.1 2012/02/13 23:05:52 jonp Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = basename(argv[0]) ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
  {
    usage_exit(1) ;
  }

  // mri_convert expects a relative warp
  fprintf(stdout, "[%s]:  reading warp file '%s'\n", Progname, argv[1]);
  fprintf(stdout, "assuming RELATIVE warp convention\n");

  // TODO: add support for absolute warps as well, add option to specify convention


  mri = MRIread(argv[1]) ;
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read warp volume %s\n", Progname,argv[1]) ;

  m = MRIgetVoxelToRasXform(mri) ;

  // NOTE: this assumes a standard siemens image orientation in which
  // case a neurological orientation means that the first frame is
  // flipped

  if ( MatrixDeterminant(m) > 0 )
    {
      fprintf(stdout, "non-negative Jacobian determinant -- converting to radiological ordering\n");
    }
  {
    // 2012/feb/08: tested with anisotropic voxel sizes

    MRI *mri2 = NULL ;
    int c=0,r=0,s=0;
    float v;

    mri2 = MRIcopy(mri,NULL);
    for(c=0; c < mri->width; c++)
      {
        for(r=0; r < mri->height; r++)
          {
            for(s=0; s < mri->depth; s++)
              {
                // only flip first frame (by negating relative shifts)
                v = MRIgetVoxVal(mri, c,r,s,0) / mri->xsize;
                if ( MatrixDeterminant(m) > 0 )
                  MRIsetVoxVal(    mri2,c,r,s,0,-v);
                else
                  MRIsetVoxVal(    mri2,c,r,s,0, v);

                v = MRIgetVoxVal(mri, c,r,s,1) / mri->ysize;
                MRIsetVoxVal(    mri2,c,r,s,1, v);

                v = MRIgetVoxVal(mri, c,r,s,2) / mri->zsize;
                MRIsetVoxVal(    mri2,c,r,s,2, v);

              }
          }
      }
    MRIfree(&mri);
    mri = mri2;

  }
  MatrixFree(&m) ;


  // this does all the work! (gcamorph.c)
  gcam = GCAMalloc(mri->width, mri->height, mri->depth) ;
  GCAMinitVolGeom(gcam, mri, mri) ;

  // not sure if removing singularities is ever a bad thing
#if 1
  GCAMremoveSingularitiesAndReadWarpFromMRI(gcam, mri) ;
#else
  GCAMreadWarpFromMRI(gcam, mri) ;
#endif

  fprintf(stdout, "[%s]:  writing warp file '%s'\n", Progname, argv[2]);

  out_name = argv[2] ;
  GCAMwrite(gcam, out_name) ;

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "conversion took %d minutes and %d seconds.\n", minutes, seconds) ;
  exit(0) ;
  return(0) ;
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
  switch (toupper(*option))
  {
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
usage_exit(int code)
{
  printf("usage: %s [options] <input FSL warp volume> <output warp volume>.m3z\n", Progname) ;
  exit(code) ;
}




