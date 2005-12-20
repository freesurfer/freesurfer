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
#include "fio.h"
#include "version.h"

static char vcid[] = 
"$Id: mris_convert.c,v 1.17 2005/12/20 22:54:13 nicks Exp $";


/*-------------------------------- CONSTANTS -----------------------------*/

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int convertToWFile(char *in_fname, char *out_fname) ;
static int convertFromWFile(char *in_fname, char *out_fname) ;
static int writeAsciiCurvFile(MRI_SURFACE *mris, char *out_fname) ;

/*-------------------------------- DATA ----------------------------*/

char *Progname ;

static int talairach_flag = 0 ;
static int patch_flag = 0 ;
static int read_orig_positions = 0 ;
static int w_file_dst_flag = 0 ;
static int w_file_src_flag = 0 ;
static int curv_file_flag = 0 ;
static char *curv_fname ;
static char *orig_surf_name = NULL ;

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[])
{
  MRI_SURFACE  *mris ;
  char         **av, *in_fname, *out_fname, fname[STRLEN], hemi[10],
    *cp, path[STRLEN], *dot, ext[STRLEN] ;
  int          ac, nargs ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: mris_convert.c,v 1.17 2005/12/20 22:54:13 nicks Exp $", "$Name:  $");
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

  if (argc < 3)
    usage_exit() ;

  in_fname = argv[1] ;
  out_fname = argv[2] ;

  dot = strrchr(out_fname, '.') ;
  if (dot)
    {
      strcpy(ext, dot+1) ;
      if (!stricmp(ext, "W")) 
        w_file_dst_flag = 1 ;
    }

  if (w_file_dst_flag)
    {
      convertToWFile(in_fname, out_fname) ;
      exit(0) ;
    }

  dot = strrchr(in_fname, '.') ;
  if (dot)
    {
      strcpy(ext, dot+1) ;
      if (!stricmp(ext, "W")) 
        w_file_src_flag = 1 ;
    }

  if (w_file_src_flag)
    {
      convertFromWFile(in_fname, out_fname) ;
      exit(0) ;
    }

  if (patch_flag)   /* read in orig surface before reading in patch */
    {
      char name[100] ;

      FileNamePath(in_fname, path) ;
      FileNameOnly(in_fname, name) ;
      cp = strchr(name, '.') ;
      if (cp)
        {
          strncpy(hemi, cp-2, 2) ;
          hemi[2] = 0 ;
        }
      else
        strcpy(hemi, "lh") ;

      sprintf(fname, "%s/%s.orig", path, hemi) ;
      mris = MRISread(fname) ;
      if (!mris)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, fname) ;
      if (MRISreadPatch(mris, in_fname) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read patch file %s",
                  Progname, in_fname) ;
      if (read_orig_positions)
        {
          if (MRISreadVertexPositions(mris, orig_surf_name) != NO_ERROR)
            ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                      Progname, orig_surf_name) ;
        }
    }
  else
    {
      mris = MRISread(in_fname) ;
      if (!mris)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, in_fname) ;
    }

  if (curv_file_flag)
    {
      int type ;
                
      MRISreadCurvatureFile(mris, curv_fname) ;
      type = MRISfileNameType(out_fname) ;
      if (type == MRIS_ASCII_FILE)
        writeAsciiCurvFile(mris, out_fname) ;
      else
        MRISwriteCurvature(mris, out_fname) ;
    }
  else if (mris->patch)
    MRISwritePatch(mris, out_fname) ;
  else
    MRISwrite(mris, out_fname) ;

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
    case 'C':
      curv_file_flag = 1 ;
      curv_fname = argv[2] ;
      nargs = 1 ;
      break ;
    case 'O':
      read_orig_positions = 1 ;
      orig_surf_name = argv[2] ;
      nargs = 1 ;
      break ;
    case 'P':
      patch_flag = 1 ;
      nargs = 0 ;
      break ;
    case 'T':
      talairach_flag = 1 ;
      break ;
    case '?':
    case 'U':
    case 'H':
      print_help() ;
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
  print_help() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <input surface file> <output surface file>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
          "\nThis program will convert an MRI surface to ascii, and "
	  "vice-versa.\n") ;
  fprintf(stderr, 
          "\nvalid options are:\n") ;
  fprintf(stderr, 
          "  -p                input is a patch not a full surface\n") ;
  fprintf(stderr, 
          "  -c <curv file>    input is curvature file (must still "
          "specify surface)\n\n") ;
  fprintf(stderr, 
          "Surface and curvature files can be ascii or binary.\n") ;
  fprintf(stderr, 
          "Ascii file is assumed if filename ends with .asc\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}
static int
convertToWFile(char *in_fname, char *out_fname)
{
  FILE   *infp, *outfp ;
  char   line[300], *cp ;
  int    vno, l = 0, num, ilat ;
  float  val ;

  fprintf(stderr, "writing w file %s...\n", out_fname) ;
  outfp = fopen(out_fname,"wb");
  if (outfp==NULL) 
    ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",Progname,out_fname) ;

  infp = fopen(in_fname,"rb");
  if (infp==NULL) 
    ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",Progname,in_fname) ;


  cp = fgetl(line, 299, infp) ;
  ilat = atoi(cp) ; /* not used at the moment */
  cp = fgetl(line, 299, infp) ;
  num = atoi(cp) ;
  fwrite2(0,outfp);
  fwrite3(num,outfp);

  while ((cp = fgetl(line, 299, infp)) != NULL)
    {
      l++ ;
      if (sscanf(cp, "%d %f", &vno, &val) != 2)
        {
          ErrorPrintf(ERROR_BADFILE,
                      "%s: could not scan parms from line %d: %s.\n",
                      Progname, l, cp) ;
          val = 0.0f ;   /* don't know what it is... */
        }
      fwrite3(vno,outfp);
      fwriteFloat(val, outfp) ;
    }
  fclose(infp) ;
  fclose(outfp) ;
  return(NO_ERROR) ;
}


static int
convertFromWFile(char *in_fname, char *out_fname)
{
  FILE   *infp, *outfp ;
  int    vno, num, ilat, i ;
  float  val, lat ;

  fprintf(stderr, "writing ascii w file %s...\n", out_fname) ;
  outfp = fopen(out_fname,"wb");
  if (outfp==NULL) 
    ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",Progname,out_fname) ;

  infp = fopen(in_fname,"rb");
  if (infp==NULL) 
    ErrorExit(ERROR_NOFILE, "%s: Can't create file %s\n",Progname,in_fname) ;

  fread2(&ilat,infp);
  lat = (float)ilat/10.0;
  fprintf(outfp, "%2.3f\n", lat) ;
  fread3(&num,infp);
  fprintf(outfp, "%d\n", num) ;
  for (i=0;i<num;i++)
    {
      fread3(&vno,infp);
      val = freadFloat(infp) ;
      fprintf(outfp, "%d  %f\n", vno, val) ;
    }
  fclose(outfp); fclose(infp) ;

  return(NO_ERROR) ;
}


static int
writeAsciiCurvFile(MRI_SURFACE *mris, char *out_fname)
{
  FILE   *fp ;
  int    vno ;
  VERTEX *v ;

  fp = fopen(out_fname, "w") ;
  if (!fp)
    ErrorExit(ERROR_BADFILE, "%s could not open output file %s.\n",
              Progname, out_fname) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      fprintf(fp, "%3.3d %2.5f %2.5f %2.5f %2.5f\n",
              vno, v->x, v->y, v->z, v->curv) ;
    }

  fclose(fp) ;
  return(NO_ERROR) ;
}

