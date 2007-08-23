/**
 * @file  mris_convert.c
 * @brief Format conversions of surface files and scalar overlay files
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2007/08/23 17:44:16 $
 *    $Revision: 1.24 $
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
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "macros.h"
#include "fio.h"
#include "version.h"
#include "matrix.h"
#include "transform.h"

int MRISmatrixMultiply(MRIS *mris, MATRIX *M);

//------------------------------------------------------------------------
static char vcid[] =
  "$Id: mris_convert.c,v 1.24 2007/08/23 17:44:16 greve Exp $";

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
static char *talxfmsubject = NULL;
static int patch_flag = 0 ;
static int read_orig_positions = 0 ;
static int w_file_dst_flag = 0 ;
static int w_file_src_flag = 0 ;
static int curv_file_flag = 0 ;
static char *curv_fname ;
static char *orig_surf_name = NULL ;
static double scale=0;
static int rescale=0;  // for rescaling group average surfaces
static int output_normals=0;
static MATRIX *XFM=NULL;
static int write_vertex_neighbors = 0;
int MRISwriteVertexNeighborsAscii(MRIS *mris, char *out_fname);

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[]) {
  MRI_SURFACE  *mris ;
  char         **av, *in_fname, *out_fname, fname[STRLEN], hemi[10],
  *cp, path[STRLEN], *dot, ext[STRLEN] ;
  int          ac, nargs ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
          (argc, argv,
           "$Id: mris_convert.c,v 1.24 2007/08/23 17:44:16 greve Exp $",
           "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3) usage_exit() ;

  in_fname = argv[1] ;
  out_fname = argv[2] ;

  if (talxfmsubject && curv_file_flag) {
    printf("ERROR: cannot specify -t and -c\n");
    exit(1);
  }

  // check whether outputis a .w file
  dot = strrchr(out_fname, '.') ;
  if (dot) {
    strcpy(ext, dot+1) ;
    if (!stricmp(ext, "W")) w_file_dst_flag = 1 ;
  }

  if (w_file_dst_flag) {
    convertToWFile(in_fname, out_fname) ; //???
    exit(0) ;
  }

  // check whether input is a .w file
  dot = strrchr(in_fname, '.') ;
  if (dot) {
    strcpy(ext, dot+1) ;
    if (!stricmp(ext, "W")) w_file_src_flag = 1 ;
  }

  if (w_file_src_flag) {
    convertFromWFile(in_fname, out_fname) ; //???
    exit(0) ;
  }

  if (patch_flag) {   /* read in orig surface before reading in patch */
    char name[100] ;

    FileNamePath(in_fname, path) ;
    FileNameOnly(in_fname, name) ;
    cp = strchr(name, '.') ;
    if (cp) {
      strncpy(hemi, cp-2, 2) ;
      hemi[2] = 0 ;
    } else
      strcpy(hemi, "lh") ;

    sprintf(fname, "%s/%s.orig", path, hemi) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    if (MRISreadPatch(mris, in_fname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read patch file %s",
                Progname, in_fname) ;
    if (read_orig_positions) {
      if (MRISreadVertexPositions(mris, orig_surf_name) != NO_ERROR)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, orig_surf_name) ;
    }
  } else {
    mris = MRISread(in_fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, in_fname) ;
  }

  if (talxfmsubject) {
    XFM = DevolveXFM(talxfmsubject, NULL, NULL);
    if (XFM == NULL) exit(1);
    printf("Applying talairach transform\n");
    MatrixPrint(stdout,XFM);
    MRISmatrixMultiply(mris,XFM);
  }

  if (rescale) {
    if (mris->group_avg_surface_area == 0) {
      printf("ERROR: cannot rescale a non-group surface\n");
      exit(1);
    }
    scale = sqrt((double)mris->group_avg_surface_area/mris->total_area);
  }

  if (scale > 0) {
    printf("scale = %lf\n",scale);
    MRISscale(mris,scale);
    MRIScomputeMetricProperties(mris);
  }

  if (curv_file_flag) {
    int type ;

    MRISreadCurvatureFile(mris, curv_fname) ;
    type = MRISfileNameType(out_fname) ;
    if (type == MRIS_ASCII_FILE)
      writeAsciiCurvFile(mris, out_fname) ;
    else
      MRISwriteCurvature(mris, out_fname) ;
  } else if (mris->patch)
    MRISwritePatch(mris, out_fname) ;
  else if (output_normals)
    MRISwriteNormalsAscii(mris, out_fname) ;
  else if (write_vertex_neighbors)
    MRISwriteVertexNeighborsAscii(mris, out_fname) ;
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
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option)) {
    case 'C':
      curv_file_flag = 1 ;
      curv_fname = argv[2] ;
      nargs = 1 ;
      break ;
    case 'N':
      output_normals = 1;
      break ;
    case 'O':
      read_orig_positions = 1 ;
      orig_surf_name = argv[2] ;
      nargs = 1 ;
      break ;
    case 'S':
      sscanf(argv[2],"%lf",&scale);
      nargs = 1 ;
      break ;
    case 'R':
      rescale = 1;
      break ;
    case 'P':
      patch_flag = 1 ;
      nargs = 0 ;
      break ;
    case 'T':
      talairach_flag = 1 ;
      talxfmsubject = argv[2] ;
      nargs = 1 ;
      break ;
    case 'V':
      write_vertex_neighbors = 1;
      nargs = 0 ;
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
usage_exit(void) {
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "Usage: %s [options] <input surface file> <output surface file>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  printf(
    "\nThis program will convert an MRI surface to ascii, and "
    "vice-versa.\n") ;
  printf( "\nValid options are:\n") ;
  printf( "  -p                input is a patch, not a full surface\n") ;
  printf( "  -c <scalar file>  input is scalar overlay file (must still\n"
          "                    specify surface)\n") ;
  printf( "  -o origname       read orig positions\n") ;
  printf( "  -s scale          scale vertex xyz by scale\n") ;
  printf( "  -r                rescale vertex xyz so total area is\n"
          "                    same as group average\n") ;
  printf( "  -t subject        apply talairach xfm of subject to\n"
          "                    vertex xyz\n");
  printf( "  -n                output is an ascii file where vertex data\n") ;
  printf( "                    is the surface normal vector\n") ;
  printf( "  -v Writes out neighbors of a vertex in each row. The first\n");
  printf( "     column is the vertex number, the 2nd col is the number of neighbors,\n");
  printf( "     the remaining cols are the vertex numbers of the neighbors.  \n");
  printf( "     Note: there can be a different number of neighbors for each vertex.\n") ;
  printf( "\n") ;
  printf( "Surface and scalar files can be ascii or binary.\n") ;
  printf( "Ascii file is assumed if filename ends with .asc\n") ;
  printf( "\n") ;
  printf( "EXAMPLES:\n") ;
  printf( "\n");
  printf( "Convert a surface file to ascii:\n");
  printf( "  mris_convert lh.white lh.white.asc\n") ;
  printf( "\n");
  printf( "Write vertex neighbors  to ascii:\n");
  printf( "  mris_convert -v lh.white lh.white.neighbors.asc\n") ;
  printf( "\n");
  printf( "Convert a surface file to ascii (vertices are surface normals):\n");
  printf( "  mris_convert -n lh.white lh.white.normals.asc\n") ;
  printf( "\n");
  printf( "Apply talairach xfm to white surface, save as binary:\n");
  printf( "  mris_convert -t bert lh.white lh.white.tal\n") ;
  printf( "\n");
  printf( "Convert a scalar overlay file to ascii:\n");
  printf( "  mris_convert -c lh.thickness lh.white lh.thickness.asc\n") ;
  printf( "\n") ;
  printf( "See also mri_surf2surf\n") ;
  exit(1) ;
}

static void
print_version(void) {
  printf( "%s\n", vcid) ;
  exit(1) ;
}

static int
convertToWFile(char *in_fname, char *out_fname) {
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

  while ((cp = fgetl(line, 299, infp)) != NULL) {
    l++ ;
    if (sscanf(cp, "%d %f", &vno, &val) != 2) {
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
convertFromWFile(char *in_fname, char *out_fname) {
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
  for (i=0;i<num;i++) {
    fread3(&vno,infp);
    val = freadFloat(infp) ;
    fprintf(outfp, "%d  %f\n", vno, val) ;
  }
  fclose(outfp);
  fclose(infp) ;

  return(NO_ERROR) ;
}


static int
writeAsciiCurvFile(MRI_SURFACE *mris, char *out_fname) {
  FILE   *fp ;
  int    vno ;
  VERTEX *v ;

  fp = fopen(out_fname, "w") ;
  if (!fp)
    ErrorExit(ERROR_BADFILE, "%s could not open output file %s.\n",
              Progname, out_fname) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    fprintf(fp, "%3.3d %2.5f %2.5f %2.5f %2.5f\n",
            vno, v->x, v->y, v->z, v->curv) ;
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}

/*
  \fn int MRISwriteVertexNeighborsAscii(MRIS *mris, char *out_fname)
  \brief Writes out neighbors of a vertex in each row. The first
  column is the vertex number, the 2nd col is the number of neighbors,
  the remaining cols are the vertex numbers of the neighbors.
*/

int MRISwriteVertexNeighborsAscii(MRIS *mris, char *out_fname)
{
  int vno, nnbrs, nbrvno;
  FILE *fp;

  fp = fopen(out_fname,"w");
  if(fp == NULL){
    printf("ERROR: opening %s\n",out_fname);
    exit(1);
  }

  for(vno=0; vno < mris->nvertices; vno++){
    nnbrs = mris->vertices[vno].vnum;
    fprintf(fp,"%6d %2d   ",vno,nnbrs);
    for (nbrvno = 0; nbrvno < nnbrs; nbrvno++)
      fprintf(fp,"%6d ",mris->vertices[vno].v[nbrvno]);
    fprintf(fp,"\n");
  }
  fclose(fp);

  return(0);
}

