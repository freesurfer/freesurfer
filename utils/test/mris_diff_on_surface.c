/*
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


/* Compute the difference of two data defined on the same surface mesh */
/* The write function somehow will screen out vextices with value exactly
 * equal to zero! No wonder my output file gets smaller
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
#include "mri.h"
#include "macros.h"
#include "mrishash.h"
#include "mri_identify.h"
#include "icosahedron.h"
#include "version.h"

#define MAX_DATA_NUMBERS 200


int main(int argc, char *argv[]) ;

int framesave = 0;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *srctypestring = NULL;
int srctype = MRI_VOLUME_TYPE_UNKNOWN;
char *trgtypestring = NULL;
int trgtype = MRI_VOLUME_TYPE_UNKNOWN;

int negflag = 0;
int debugflag = 0;
int debugvtx = 0;
int pathflag = 0;

const char *Progname ;


int main(int argc, char *argv[])
{
  char **av, *surf_name, *out_prefix, *fname;

  int nargs, ac, i, nsubjects, total, index;

  double scalar, std, tmp, maxV, minV, meanV;

  MRI *SrcVals[2], *AvgVals;

  MRI_SURFACE *BaseSurf;

  nargs = handleVersionOption(argc, argv, "mris_diff_on_surface");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

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

  /* command line: <surf> <datafile 1> <datafile 2> <output prefix> */

  if (argc != 5)
    usage_exit();

  surf_name = argv[1];
  out_prefix = argv[argc - 1];

  if (srctypestring == NULL || trgtypestring == NULL)
  {
    printf("Please specify both input and output data type!\n");
    usage_exit();
  }

  printf("Reading underlying surface file\n");
  BaseSurf = MRISread(surf_name);
  if (!BaseSurf)
    ErrorExit(ERROR_NOFILE, "%s:could not read surface %s", Progname, surf_name);

  printf("Base surface has %d vertices\n", BaseSurf->nvertices);


  /* Read in the first data file */
  fname = argv[2];
  /* only two data types are supported */
  if (!strcmp(srctypestring,"curv"))
  { /* curvature file */
    if (MRISreadCurvatureFile(BaseSurf, fname) != 0)
    {
      printf("ERROR: reading curvature file\n");
      exit(1);
    }
    SrcVals[0] = MRIcopyMRIS(NULL, BaseSurf, 0, "curv");
  }
  else if (!strcmp(srctypestring,"paint") || !strcmp(srctypestring,"w"))
  {
    MRISreadValues(BaseSurf,fname);
    SrcVals[0] = MRIcopyMRIS(NULL, BaseSurf, 0, "val");
  }
  else
  {
    printf("ERROR: unknown data file format\n");
    exit(1);
  }

  if (SrcVals[0] == NULL)
  {
    fprintf(stderr, "ERROR loading data values from %s\n", fname);
  }

  /* Read in the second data file */
  fname = argv[3];
  /* only two data types are supported */
  if (!strcmp(srctypestring,"curv"))
  { /* curvature file */
    if (MRISreadCurvatureFile(BaseSurf, fname) != 0)
    {
      printf("ERROR: reading curvature file\n");
      exit(1);
    }
    SrcVals[1] = MRIcopyMRIS(NULL, BaseSurf, 0, "curv");
  }
  else if (!strcmp(srctypestring,"paint") || !strcmp(srctypestring,"w"))
  {
    MRISreadValues(BaseSurf,fname);
    SrcVals[1] = MRIcopyMRIS(NULL, BaseSurf, 0, "val");
  }
  else
  {
    printf("ERROR: unknown data file format\n");
    exit(1);
  }

  if (SrcVals[1] == NULL)
  {
    fprintf(stderr, "ERROR loading data values from %s\n", fname);
  }

  if (debugflag)
  {
    for (i=0; i < 2; i++)
    {
      printf("Data%d at vertex %d has value %g\n",i, debugvtx,  MRIFseq_vox(SrcVals[i], debugvtx, 0, 0, 0));
    }
  }

#if 0
  AvgVals = MRIclone(SrcVals[0], NULL);

  if (negflag) /* Add the two data sets */
    AvgVals = MRIadd(SrcVals[0], SrcVals[1], AvgVals);
  else /* Data1 - Data2 */
    AvgVals = MRIsubtract(SrcVals[0], SrcVals[1], AvgVals);
#endif

  AvgVals = MRIcopy(SrcVals[0], NULL);

  if (negflag)
  {
    for (index=0; index < BaseSurf->nvertices; index++)
    {
      MRIFseq_vox(AvgVals, index, 0, 0, 0) =  MRIFseq_vox(SrcVals[0], index, 0, 0, 0) +
                                              MRIFseq_vox(SrcVals[1], index, 0, 0, 0);
    }
  }
  else
  {
    for (index=0; index < BaseSurf->nvertices; index++)
    {
      MRIFseq_vox(AvgVals, index, 0, 0, 0) =  MRIFseq_vox(SrcVals[0], index, 0, 0, 0) -
                                              MRIFseq_vox(SrcVals[1], index, 0, 0, 0);
    }
  }

  maxV = -1000.0;
  minV = 1000.0;
  meanV=0.0;

  for (index=0; index < BaseSurf->nvertices; index++)
  {
    scalar = MRIFseq_vox(AvgVals, index, 0, 0, 0);
    if (maxV < scalar) maxV = scalar;
    if (minV > scalar) minV = scalar;
    meanV += scalar;
  }

  meanV /= BaseSurf->nvertices;

  printf("Output max = %g, min = %g, mean = %g\n", maxV, minV, meanV);

  if (debugflag)
  {
    printf("Output at vertex %d has value %g\n", debugvtx,  MRIFseq_vox(AvgVals, debugvtx, 0, 0, 0));
  }

  if (pathflag)
    sprintf(fname, "%s", out_prefix);
  else
  {
    if (negflag)
      sprintf(fname, "%s.sum.w", out_prefix) ;
    else
      sprintf(fname, "%s.diff.w", out_prefix) ;
  }

  if (!strcmp(trgtypestring,"paint") || !strcmp(trgtypestring,"w"))
  {

    /* This function will remove a zero-valued vertices */
    /* Make sense, since default value is considered as zero */
    /* But it will confuse the processing with matlab! */
    /* So I copy the data to the curv field to force every value is
     *  written out
     */
    /* MRIScopyMRI(BaseSurf, AvgVals, framesave, "val");*/
    /* MRISwriteValues(BaseSurf,fname); */
    MRIScopyMRI(BaseSurf, AvgVals, framesave, "curv");
    MRISwriteCurvatureToWFile(BaseSurf,fname);

  }
  else
  {
    fprintf(stderr, "ERROR unknown output file format.\n");
  }

  /* Free memories */
  MRISfree(&BaseSurf);
  MRIfree(&AvgVals);
  for (i=0; i < 2; i++)
  {
    MRIfree(&SrcVals[i]);
  }

  return 0;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

/* --------------------------------------------- */
static void print_usage(void)
{
  fprintf(stdout, "USAGE: %s [options] surface data1 data2 output_prefix\n",Progname) ;
  fprintf(stdout, "\n");
  fprintf(stdout, "Options:\n");
  fprintf(stdout, "   -src_type %%s input surface data format\n");
  fprintf(stdout, "   -trg_type  %%s output format\n");
  fprintf(stdout, "   -neg  take negative of data2, thus compute sum!\n");
  fprintf(stdout, "\n");
  std::cout << getVersion() << std::endl;
  printf("\n");

}

/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf(
    "This program computes the difference of two data sets defined on the \n"
    "the same surface mesh. Result = data1 - data2. \n"
    "\n"
    "OPTIONS\n"
    "\n"
    "  -src_type typestring\n"
    "\n"
    "    Format type string. Can be either curv (for FreeSurfer curvature file), \n"
    "    paint or w (for FreeSurfer paint files)."
    "\n"
    "  -trg_type typestring\n"
    "\n"
    "    Format type string. Can be paint or w (for FreeSurfer paint files)."
    "\n"
    "  -neg\n"
    "    Take negative of data2, thus results in the sum of the two data sets."
    "\n");

  exit(1) ;
}


/* --------------------------------------------- */
static void print_version(void)
{
  fprintf(stdout, "%s\n", getVersion().c_str()) ;
  exit(1) ;
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
  else if (!stricmp(option, "src_type"))
  {
    srctypestring = argv[2];
    srctype = string_to_type(srctypestring);
    nargs = 1 ;
  }
  else if (!stricmp(option, "trg_type"))
  {
    trgtypestring = argv[2];
    trgtype = string_to_type(srctypestring);
    nargs = 1 ;
  }
  else if (!stricmp(option, "abspath"))
  {
    pathflag = 1;
  }
  else if (!stricmp(option, "neg"))
  {
    negflag = 1 ;
  }
  else if (!stricmp(option, "debug"))
  {
    debugflag = 1;
    debugvtx = atoi(argv[2]);
    nargs = 1;
  }
  else
  {
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    print_help() ;
    exit(1) ;
  }

  return(nargs) ;
}
