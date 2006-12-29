/**
 * @file  mri_add_xform_to_header.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:04 $
 *    $Revision: 1.7 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "volume_io.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "version.h"
#include "fio.h"
#include "cmdargs.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void usage_exit(void);

char *Progname ;

int verbose = 0 ;
int CopyNameOnly = 0;

int
main(int argc, char *argv[]) {
  char   **av ;
  int    ac, nargs ;
  MRI    *mri ;
  char   *xform_fname, *in_fname, *out_fname ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_add_xform_to_header.c,v 1.7 2006/12/29 02:09:04 nicks Exp $", "$Name:  $");

  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 1) usage_exit();

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
    ErrorExit(ERROR_BADPARM, "%s: no transform name specified", Progname) ;
  xform_fname = argv[1] ;

  if (argc < 3)
    ErrorExit(ERROR_BADPARM, "%s: no input name specified", Progname) ;
  in_fname = argv[2] ;

  if (argc < 4)  out_fname = in_fname ;
  else           out_fname = argv[3] ;

  if (verbose)  fprintf(stderr, "reading from %s...", in_fname) ;

  // we have two cases, in_fname is just a directory name or .mgz
  if (fio_IsDirectory(in_fname))mri = MRIreadInfo(in_fname) ; // must be old COR volume
  else if (fio_FileExistsReadable(in_fname)) {
    char *ext = fio_extension(in_fname);
    if (ext==0)
      ErrorExit(ERROR_BADPARM, "%s: no extension found", Progname) ;
    printf("INFO: extension is %s\n", ext);
    if (strcmp(ext, "mgz")==0 || strcmp(ext, "mgh")==0)
      mri = MRIread(in_fname);      // mgh or mgz
    else {
      ErrorExit(ERROR_BADPARM,
                "%s: currently only .mgz or .mgh saves transform name", Progname) ;
    }
  }
  if (!mri)
    ErrorExit(ERROR_NO_FILE, "%s: could not open source file %s",
              Progname, xform_fname) ;

  if (! CopyNameOnly) {
    // why do we need to load the transform at this time
    // mri is removed anyway???? -- good point, added -s for noload
    if (input_transform_file(xform_fname, &mri->transform) != OK)
      ErrorPrintf(ERROR_NO_MEMORY,
                  "%s: could not read xform file '%s'\n", Progname, xform_fname);
    // my guess is just to verify the validity of the transform?
    mri->linear_transform = get_linear_transform_ptr(&mri->transform) ;
    mri->inverse_linear_transform =
      get_inverse_linear_transform_ptr(&mri->transform) ;
    mri->free_transform = 1 ;
  }
  strcpy(mri->transform_fname, xform_fname) ;

  if (verbose)
    fprintf(stderr, "done.\nwriting to %s...", out_fname) ;
  // this writes COR-.info only
  if (fio_IsDirectory(out_fname))  MRIwriteInfo(mri, out_fname) ;
  else     MRIwrite(mri, out_fname);  // currently only mgh format write xform info

  if (verbose)  fprintf(stderr, "done.\n") ;

  MRIfree(&mri) ;
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
  switch (toupper(*option)) {
  case 'C':
    CopyNameOnly = 1;
    break ;
  case 'V':
    verbose = !verbose ;
    break ;
  case 'N':
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    printf("usage: %s [xform file name] [input directory] [output directory]\n",
           argv[0]) ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s <options> xfmfile invol outvol \n",Progname) ;
  printf("\n");
  printf("   -v : verbose \n");
  printf("   -c : do not try to load the xfmfile, just copy name\n");
  printf("\n");
}
