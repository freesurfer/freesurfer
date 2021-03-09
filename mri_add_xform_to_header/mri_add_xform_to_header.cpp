/**
 * @brief Just adds specified xform to the volume header
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

#include "utils.h"
#include "mri.h"
#include "macros.h"
#include "minc.h"
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

const char *Progname ;

int verbose = 0 ;
int CopyNameOnly = 0;

int
main(int argc, char *argv[])
{
  char   **av ;
  int    ac, nargs ;
  MRI    *mri=NULL ;
  char   *xform_fname=NULL, *in_fname=NULL, *out_fname=NULL ;

  nargs = handleVersionOption(argc, argv, "mri_add_xform_to_header");

  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 1)
  {
    usage_exit();
  }

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
  {
    ErrorExit(ERROR_BADPARM, "%s: no transform name specified", Progname) ;
  }
  xform_fname = argv[1] ;

  if (argc < 3)
  {
    ErrorExit(ERROR_BADPARM, "%s: no input name specified", Progname) ;
  }
  in_fname = argv[2] ;

  if (argc < 4)
  {
    out_fname = in_fname ;
  }
  else
  {
    out_fname = argv[3] ;
  }

  if (verbose)
  {
    fprintf(stderr, "reading from %s...", in_fname) ;
  }

  // we have two cases, in_fname is just a directory name or .mgz
  if (fio_IsDirectory(in_fname))
  {
    mri = MRIreadInfo(in_fname) ;  // must be old COR volume
  }
  else if (fio_FileExistsReadable(in_fname))
  {
    char *ext = fio_extension(in_fname);
    if (ext==0)
    {
      ErrorExit(ERROR_BADPARM, "%s: no extension found", Progname) ;
    }
    printf("INFO: extension is %s\n", ext);
    if (strcmp(ext, "mgz")==0 || strcmp(ext, "mgh")==0)
    {
      mri = MRIread(in_fname);  // mgh or mgz
    }
    else
    {
      ErrorExit(ERROR_BADPARM,
                "%s: currently only .mgz or .mgh saves transform name",
                Progname) ;
    }
  }
  if (!mri)
    ErrorExit(ERROR_NO_FILE, "%s: could not open source file %s",
              Progname, in_fname) ;

  if (! CopyNameOnly)
  {
    // why do we need to load the transform at this time
    // mri is removed anyway???? -- good point, added -s for noload
    if (input_transform_file(xform_fname, &mri->transform) != OK)
      ErrorPrintf(ERROR_NO_MEMORY,
                  "%s: could not read xform file '%s'\n",
                  Progname, xform_fname);
    // my guess is just to verify the validity of the transform?
    mri->linear_transform = get_linear_transform_ptr(&mri->transform) ;
    mri->inverse_linear_transform =
      get_inverse_linear_transform_ptr(&mri->transform) ;
    mri->free_transform = 1 ;
  }
  strcpy(mri->transform_fname, xform_fname) ;

  if (verbose)
  {
    fprintf(stderr, "done.\nwriting to %s...", out_fname) ;
  }
  // this writes COR-.info only
  if (fio_IsDirectory(out_fname))
  {
    MRIwriteInfo(mri, out_fname) ;
  }
  else
  {
    MRIwrite(mri, out_fname);  // currently only mgh format write xform info
  }

  if (verbose)
  {
    fprintf(stderr, "done.\n") ;
  }

  MRIfree(&mri) ;
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
  if (!stricmp(option,"-help")||!stricmp(option,"-usage"))
  {
    usage_exit();
  }
  else switch (toupper(*option))
    {
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
    case 'H':
    case 'U':
      print_usage();
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
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

/* --------------------------------------------- */
#include "mri_add_xform_to_header.help.xml.h"
static void print_usage(void)
{
  outputHelpXml(mri_add_xform_to_header_help_xml,
                mri_add_xform_to_header_help_xml_len);
}
