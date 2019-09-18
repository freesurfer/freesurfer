/**
 * @file  mris_pval_fill.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:33 $
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
// mris_pval_fill.c
//
// written by Peng Yu
// date: 11/11/04
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2011/03/02 00:04:33 $
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
#include "mrisurf.h"
#include "icosahedron.h"


//static char vcid[] = "$Id: mris_pval_fill.c,v 1.4 2011/03/02 00:04:33 nicks Exp $";

int             main(int argc, char *argv[]) ;
static int      get_option(int argc, char *argv[]) ;
const char            *Progname ;

int
main(int argc, char *argv[]) {
  int         nargs, msec, i=0, order=7;
  Timer then ;
  MRIS        *mris_in, *mris_out;
  MRI_SP      *mrisp ;
  FILE        *fp_p;
  float       val, temp;

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
              "usage: %s <cc_volume> <medial axis file> <P_value> <output_volume> ", Progname);

  then.reset() ;

  fprintf(stdout, "reading surface from %s\n", argv[1]);
  mris_in = MRISread(argv[1]) ;
  if (!mris_in)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, argv[1]) ;

  MRISreadOriginalProperties(mris_in, argv[2]) ;
  fprintf(stdout, "Reading original surface from %s orig area is %f\n", argv[2],mris_in->orig_area);

  mris_out = ReadIcoByOrder(order, 100);
  //for (m = 0; m<mris_out->nvertices; m++)
  // mris_out->vertices[m].nsize=1;
  mrisp = MRISPalloc(1, 3);
  MRIScoordsToParameterization(mris_in, mrisp, 1, ORIGINAL_VERTICES) ;
  //MRISPblur(mrisp, mrisp, 1, 0);
  //MRISPblur(mrisp, mrisp, 1, 1);
  //MRISPblur(mrisp, mrisp, 1, 2);
  MRIScoordsFromParameterization(mrisp, mris_out, ORIGINAL_VERTICES) ;

  if ((fp_p = fopen(argv[3], "r")) == NULL) {
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "p value input: file %s does not exist!", argv[3]));
  }
  fprintf(stdout, "reading p value from %s\n",argv[3]);


  for (i=0; i<IcoNVtxsFromOrder(order-1)*3; i++) {
    fscanf(fp_p,"%f %f\n", &val, &temp);
    if (val<0.05)
      mris_out->vertices[i].curv=1-val;
    fprintf(stdout, "%d   %f \n", i, val);
  }


  MRISwriteCurvature(mris_out, argv[4]) ;
  fprintf(stdout, "writing output surface to %s\n", argv[4]);

  MRISfree(&mris_in);
  MRISfree(&mris_out);
  fclose(fp_p);
  msec = then.milliseconds() ;

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






