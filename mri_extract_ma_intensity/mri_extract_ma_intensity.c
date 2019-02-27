/**
 * @file  mri_extract_ma_intensity.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:15 $
 *    $Revision: 1.4 $
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
// mri_extract_ma_intensity.c
//
// written by Peng Yu
// date: 11/05/04
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2011/03/02 00:04:15 $
// Revision       : $Revision: 1.4 $
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
#include "mrimorph.h"
#include "timer.h"
#include "mrinorm.h"
#include "cma.h"
#include "version.h"
#include "error.h"
#include "volume_io/geom_structs.h"
#include "transform.h"
#include "talairachex.h"
#include "matrix.h"
#include "mriTransform.h"

//static char vcid[] = "$Id: mri_extract_ma_intensity.c,v 1.4 2011/03/02 00:04:15 nicks Exp $";

int             main(int argc, char *argv[]) ;
static int      get_option(int argc, char *argv[]) ;
const char            *Progname ;

static double cc_tal_x = 126 ;
static double cc_tal_y = 107 ;
static double cc_tal_z = 96 ;

int
main(int argc, char *argv[]) {
  int         nargs, msec, n_sample = 100, n=0, i=0;
  Timer then ;
  MRI         *mri_pd, *mri_t1;
  FILE        *fp_in, *fp_out;
  float       x, y, radius, val=0;
  double      val_pd=0, val_t1=0;

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    ErrorExit(ERROR_BADPARM,
              "usage: %s <medial_axis_file> <output_file> <input_volume_1> <input_volume_2> ", Progname);

  then.reset() ;

  if ((fp_in = fopen(argv[1], "r")) == NULL) {
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "medial axis measurement: file %s does not exist!", argv[1]));
  }
  fprintf(stdout, "reading medial axis atoms from %s\n",argv[1]);


  if ((fp_out = fopen(argv[1], "w")) == NULL) {
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "intensity output: file %s does not exist!", argv[2]));
  }
  fprintf(stdout, "reading medial axis atoms from %s\n",argv[2]);

  if ( argc > 2) {
    fprintf(stdout, "reading proton density map from %s\n", argv[3]);
    mri_pd = MRIread(argv[3]) ;
    if (!mri_pd)
      ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname, argv[3]) ;
  }

  if ( argc >3 ) {
    fprintf(stdout, "reading T1 map from %s\n", argv[4]);
    mri_t1 = MRIread(argv[4]) ;
    if (!mri_t1)
      ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname, argv[4]) ;
  }

  for (i=0; i<n_sample; i++) {
    fscanf(fp_in, "%d   %f   %f   %f   %f\n", &n, &x, &y, &radius, &val);
    MRIsampleVolumeFrameType(mri_pd, cc_tal_x, y, x, 0, SAMPLE_TRILINEAR, &val_pd) ;
    MRIsampleVolumeFrameType(mri_t1, cc_tal_x, y, x, 0, SAMPLE_TRILINEAR, &val_t1) ;
    fprintf(fp_out, "%d   %f   %f   %f   %f   %f   %f\n", n, x, y, radius, val, val_pd, val_t1);
  }

  MRIfree(&mri_t1) ;
  MRIfree(&mri_pd);
  msec = then.milliseconds() ;
  fclose(fp_in);
  fclose(fp_out);
  exit(0) ;
  return(0) ;
}



static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */

  if (!stricmp(option, "seed")) {
    cc_tal_x = atof(argv[2]) ;
    cc_tal_y = atof(argv[3]) ;
    cc_tal_z = atof(argv[4]) ;
    nargs = 3 ;
    fprintf(stderr, "cc seed at (%f %f %f)\n",
            cc_tal_x, cc_tal_y, cc_tal_z) ;
  } else switch (toupper(*option)) {
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






