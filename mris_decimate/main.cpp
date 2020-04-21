/*
 * Original Author: Dan Ginsburg (@ Children's Hospital Boston)
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

/// \brief Brief description
/// Reduce the number of vertices and faces in a surface.
///
/// \b NAME
///
/// mris_decimate
///
/// \b SYNPOSIS
///
/// mris_decimate [options] <input surface> <output surface>
///
/// \b DESCRIPTION
///
///  This tool reduces the number of triangles in a surface using the
///  the GNU Triangulated Surface Library documented at:
///
///           http://gts.sourceforge.net/reference/book1.html
///
///  Please see the GTS documentation for details on the decimation algorithm.
///  mris_decimate will read in an existing surface and write out a new one
///  with less triangles.  The decimation level and other options can be provided
///  on the command-line.
///



#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <sstream>
#include <iostream>

#include "mris_decimate.h"



#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"




///
//  Function Prototypes
//
static int  get_option(int argc, char **argv);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

///
//  Global Variables
//
const char *Progname = NULL;
char *cmdline;
int debug=0;
int checkoptsonly=0;
struct utsname uts;
DECIMATION_OPTIONS gDecimationOptions;

/////////////////////////////////////////////////////////////////////
//
//  Public Functions
//
//

//////////////////////////////////////////////////////////////////
//
//  Private Functions
//
//

///
/// Callback to print out status messages during decimation
///
void DecimateProgressCallback(float percent, const char *msg, void *userData)
{
  std::cout << std::string(msg) << std::endl;
}


///
/// \fn Main entrypoint for mris_decimate
/// \return 0 on succesful run, 1 on error
///
int main(int argc, char *argv[])
{
  // Initialize Decimation options
  memset(&gDecimationOptions, 0, sizeof(DECIMATION_OPTIONS));
  gDecimationOptions.desiredNumFaces = -1;
  gDecimationOptions.desiredFaceArea = -1;
  gDecimationOptions.decimationLevel = 0.5; // Default decimation level if not specified

  // This flag is to run code to sort output vertices. This was needed
  // when compiling GTS with hashes instead of btrees to make the
  // output deterministic. Even then it did not always give the same
  // output when the input was ?h.orig.nofix, because the underlying
  // algorithm did not give the same vertices. The algorithm appears
  // to give deterministic output when using btrees, so the sorting
  // feature is off by default.
  gDecimationOptions.SortVertices = 0;

  char *in_fname, out_fpath[STRLEN] ;
  int nargs;
  MRI_SURFACE *mris;
  double avgfacearea ;

  nargs = handleVersionOption(argc, argv, "main");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  Progname = argv[0] ;
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);

  if (argc < 3)
  {
    usage_exit();
  }

  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
  {
    usage_exit();
  }

  if(gDecimationOptions.desiredNumFaces > 0 && gDecimationOptions.desiredFaceArea > 0){
    printf("ERROR: cannot set -n and -a\n");
    exit(1);
  }

  in_fname = argv[1] ;
  FileNameAbsolute(argv[2], out_fpath);

  mris = MRISread(in_fname) ;
  if (!mris)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;
  }

  std::cout << "Original Surface Number of vertices: " << mris->nvertices << std::endl;
  std::cout << "Original Surface Number of faces: " << mris->nfaces << std::endl;

  if(gDecimationOptions.desiredNumFaces > 0){
    gDecimationOptions.decimationLevel = (float)gDecimationOptions.desiredNumFaces/mris->nfaces;
    printf("Setting decimation level to %f based on %d/%d\n",gDecimationOptions.decimationLevel,
	   gDecimationOptions.desiredNumFaces,mris->nfaces);
  }

  MRIScomputeMetricProperties(mris);
  avgfacearea = mris->total_area/mris->nfaces;
  printf("Average Face Area of input is %8.6f\n",avgfacearea);
  if(gDecimationOptions.desiredFaceArea > 0){
    gDecimationOptions.decimationLevel = avgfacearea/gDecimationOptions.desiredFaceArea;
    printf("Setting decimation level to %f based on %g/%g\n",gDecimationOptions.decimationLevel,
	   avgfacearea,gDecimationOptions.desiredFaceArea);
    if(gDecimationOptions.decimationLevel > 1.0){
      printf("  INFO: decimation level > 1, so setting to 1. There will be no change to the output\n");
    }
  }

  dump_options(stdout);

  // Decimate the surface
  decimateSurface(&mris, gDecimationOptions, DecimateProgressCallback);

  MRIScomputeMetricProperties(mris);
  avgfacearea = mris->total_area/mris->nfaces;
  printf("Average Face Area of output is %8.6f\n",avgfacearea);

  // Write out the results
  MRISwrite(mris, out_fpath);
  MRISfree(&mris);

  return 0;
}

///
/// \fn int get_option(int argc, char **argv)
/// \brief Parses a command-line argument
/// \param argc - number of command line arguments
/// \param argv - pointer to a character pointer
///
static int get_option(int argc, char *argv[])
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
  else switch (toupper(*option))
    {
    case 'D':
      gDecimationOptions.decimationLevel = atof(argv[2]) ;
      printf("using decimation = %2.2f\n", gDecimationOptions.decimationLevel) ;
      nargs = 1 ;
      break ;
    case 'A':
      gDecimationOptions.desiredFaceArea = atof(argv[2]) ;
      printf("desired face area = %7.6f\n", gDecimationOptions.desiredFaceArea);
      nargs = 1 ;
      break ;
    case 'N':
      gDecimationOptions.desiredNumFaces = atoi(argv[2]) ;
      printf("desired number of vertices = %d\n", gDecimationOptions.desiredNumFaces);
      nargs = 1 ;
      break ;
    case 'Q':
      gDecimationOptions.SortVertices = 1;
      printf("turning on vertex sorting\n");
      break ;
    case 'M':
      gDecimationOptions.setMinimumAngle = true;
      gDecimationOptions.minimumAngle = atof(argv[2]);
      printf("using minimumAngle = %f\n", gDecimationOptions.minimumAngle) ;
      nargs = 1;
      break;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}


///
/// \fn static void usage_exit(void)
/// \brief Prints usage and exits
///
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

///
/// \fn static void print_usage(void)
/// \brief Prints usage and returns (does not exit)
///
#include "mris_decimate.help.xml.h"
static void print_usage(void)
{
  outputHelpXml(mris_decimate_help_xml,
                mris_decimate_help_xml_len);
}

///
/// \fn static void print_help(void)
/// \brief Prints help and exits
///
static void print_help(void)
{
  print_usage() ;
  exit(1) ;
}

///
/// \fn static void print_version(void)
/// \brief Prints version and exits
///
static void print_version(void)
{
  std::cout << getVersion() << std::endl;
  exit(1) ;
}


///
/// \fn static void dump_options(FILE *fp)
/// \brief Prints command-line options to the given file pointer
/// \param FILE *fp - file pointer
///
static void dump_options(FILE *fp)
{
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());

  fprintf(fp,"\ndecimationLevel            %f\n", gDecimationOptions.decimationLevel);
  if (gDecimationOptions.setMinimumAngle)
  {
    fprintf(fp,"minimumAngle               %f\n", gDecimationOptions.minimumAngle) ;
  }
  else
  {
    fprintf(fp,"minimumAngle               %f\n", 1.0f) ;
  }

  return;
}
