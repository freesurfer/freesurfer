/**
 * @file  mri_cc_ma_fill.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:14 $
 *    $Revision: 1.3 $
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


///////////////////////////////////////////
// mri_cc_ma_fill.c
//
// written by Peng Yu
// date: 11/11/04
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2011/03/02 00:04:14 $
// Revision       : $Revision: 1.3 $
////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include "mri.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "proto.h"
#include "timer.h"
#include "cma.h"
#include "version.h"
#include "error.h"

//static char vcid[] = "$Id: mri_cc_ma_fill.c,v 1.3 2011/03/02 00:04:14 nicks Exp $";


typedef struct medial_axis_type_ {
  float x,y;            /* curr position */
  float nx,ny;         /* curr normal */
  float dx, dy ;     /* current change in position */
  float odx, ody ;  /* last change of position (for momentum) */
  float radius ;
}
atom_type, MEDATOM ;


int             main(int argc, char *argv[]) ;
static int      get_option(int argc, char *argv[]) ;
const char            *Progname ;

int
main(int argc, char *argv[]) {
  int         nargs, msec, length = 100, n=0, i=0, j=0, k=0;
  int         cc_x=0, width, height, depth;
  Timer then ;
  MRI         *mri_cc, *mri_ma, *mri_out;
  FILE        *fp;
  float       x, y, val;
  MEDATOM     *medial_axis;

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    ErrorExit(ERROR_BADPARM,
              "usage: %s <cc_volume> <P_value> <output_volume> ", Progname);

  then.reset() ;

  fprintf(stdout, "reading corpus callosum volume from %s\n", argv[1]);
  mri_cc = MRIread(argv[1]) ;
  if (!mri_cc)
    ErrorExit(ERROR_NOFILE, "%s: could not read cc volume %s", Progname, argv[1]) ;

  if ((fp = fopen(argv[2], "r")) == NULL) {
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "p value input: file %s does not exist!", argv[2]));
  }
  fprintf(stdout, "reading p value from %s\n",argv[2]);

  width = mri_cc->width;
  height = mri_cc->height;
  depth = mri_cc->depth;

  /*find mid-sagittal plane*/
  for (k = floor(depth/2)-2; k < floor(depth/2)+2; k++) {
    for (j = floor(height/2)-2; j < floor(height/2)+2 ; j++) {
      for (i = 0 ; i < mri_cc->width ; i++) {
        if MRIvox(mri_cc,i,j,k)>0
          {
            cc_x=i;
            break;
          }
        }
    }
  }

  medial_axis = (MEDATOM *)calloc(length, sizeof(MEDATOM));
  for (i=0; i<length; i++) {
    fscanf(fp, "%d %f %f %f %f\n", &n, &medial_axis[i].x, &medial_axis[i].y, &medial_axis[i].radius, &val);
    fprintf(stdout, "%d   %f   %f   %f   %f\n", n, medial_axis[i].x, medial_axis[i].y, medial_axis[i].radius, val);
  }

  mri_out = MRIcopy(mri_cc, NULL) ;
  MRIvalueFill(mri_out, 0) ;

  for (k = 0; k < depth ; k++)
    for (j = 0; j < height ; j++)
      if ( MRIvox(mri_cc, cc_x, j, k)>0 ) {
        nearest = 1000;
        for (i=0; i<length; i++) {
          x = medial_axis[i].x - 0.5 ;
          y = medial_axis[i].y - 0.5 ;
          dist = sqrt((k-x)*(k-x) + (j-y)*(j-y));
          if ( dist<nearest )
            MRIvox(mri_out,cc_x,j,k) = nint(100*medial_axis[i].radius);
          else if ( dist==nearest )
            fprintf(stdout, "vertex %d %d %d has more than one nearest neighbor!", cc_x, j, k);
        }
      }

  MRIwrite(mri_out,argv[3]);
  fprintf(stdout, "writing output volue to %s", argv[3]);

  MRIfree(&mri_cc);
  free(medial_axis);
  msec = then.milliseconds() ;
  fclose(fp);
  exit(0) ;
  return(0) ;
}



static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */

  switch (toupper(*option)) {
  case '?':
  case 'U':
    fprintf(stdout,
            "usage: %s <input medial axis text file> <output text file> <PD volume> <T1 volume>\n",
            Progname) ;
    exit(1) ;
    break ;
  default:
    fprintf(stdout, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }
  return(nargs) ;
}






