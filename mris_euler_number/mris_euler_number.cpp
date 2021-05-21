/**
 * @brief calculates the Euler number of a surface (=2 if perfect)
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
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static float curv_thresh = 2.0f ;
static int patch_flag = 0 ;
char *outfile=NULL;

int
main(int argc, char *argv[])
{
  char         **av, *in_fname, fname[STRLEN] ;
  int          ac, nargs, nvertices, nfaces, nedges, eno, dno ;
  MRI_SURFACE  *mris ;

  nargs = handleVersionOption(argc, argv, "mris_euler_number");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
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

  if (argc < 2)
  {
    usage_exit() ;
  }

  in_fname = argv[1] ;

  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;
  eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
  fprintf(stdout, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
          nvertices, nedges, nfaces, eno, 1-eno/2) ;

  fprintf(stdout, "      F =2V-4:          %d %s= %d-4 (%d)\n",
          nfaces, nfaces == 2*nvertices-4 ? "" : "!", 2*nvertices,
          2*nvertices-4-nfaces) ;
  fprintf(stdout, "      2E=3F:            %d %s= %d (%d)\n",
          2*nedges, 2*nedges == 3*nfaces ? "" : "!", 3*nfaces,
          2*nedges-3*nfaces) ;

  dno = MRIStopologicalDefectIndex(mris) ;
  fprintf(stdout, "\ntotal defect index = %d\n", dno) ;
  if(outfile && !patch_flag){
    // write out number of holes
    FILE *fp;
    fp = fopen(outfile,"w");
    fprintf(fp,"%5d\n",1-eno/2);
    fclose(fp);
  }

  if (patch_flag)
  {
    MRISremoveTopologicalDefects(mris, curv_thresh) ;
    fprintf(stdout, "\nafter editing:\n") ;

    eno = MRIScomputeEulerNumber(mris, &nvertices, &nfaces, &nedges) ;
    fprintf(stdout, "euler # = v-e+f = 2g-2: %d - %d + %d = %d --> %d holes\n",
            nvertices, nedges, nfaces, eno, 2-eno) ;

    fprintf(stdout, "      F =2V-4:          %d %s= %d-4 (%d)\n",
            nfaces, nfaces == 2*nvertices-4 ? "" : "!", 2*nvertices,
            2*nvertices-4-nfaces) ;
    fprintf(stdout, "      2E=3F:            %d %s= %d (%d)\n",
            2*nedges, 2*nedges == 3*nfaces ? "" : "!", 3*nfaces,
            2*nedges-3*nfaces) ;

    dno = MRIStopologicalDefectIndex(mris) ;
    fprintf(stdout, "total defect index = %d\n", dno) ;

    sprintf(fname, "%s.edit", in_fname) ;
    fprintf(stdout, "writing out patched surface to %s\n", fname) ;
    MRISwritePatch(mris, fname) ;
    if(outfile){
      // write out number of holes
      FILE *fp;
      fp = fopen(outfile,"w");
      fprintf(fp,"%5d\n",2-eno);
      fclose(fp);
    }
  }


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
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else switch (toupper(*option))
    {
    case 'P':
      patch_flag = 1 ;
      break ;
    case '?':
    case 'H':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'T':
      curv_thresh = (float)atof(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'O':
      outfile = argv[2];
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
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

#include "mris_euler_number.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mris_euler_number_help_xml,
                mris_euler_number_help_xml_len);
}

static void
print_help(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

