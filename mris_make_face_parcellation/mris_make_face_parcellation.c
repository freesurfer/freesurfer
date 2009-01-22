/**
 * @file  mris_make_face_parcellation.c
 * @brief make a parcellation where each unit is assigned based on the face it maps to in a specified spherical surface
 *        (usually an icosahedral one)
 *
 * 
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2009/01/22 14:01:39 $
 *    $Revision: 1.2 $
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
#include "tags.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "mrishash.h"

static char vcid[] =
  "$Id: mris_make_face_parcellation.c,v 1.2 2009/01/22 14:01:39 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

int
main(int argc, char *argv[]) {
  char               **av, *in_fname, *ico_fname, *out_fname, path[STRLEN], ico_name[STRLEN] ;
  int                ac, nargs ;
  float              scale ;
  MRI_SURFACE        *mris, *mris_ico ;
  float              radius ;
  int                fno, r, g, b, vno, annot ;
  double             fdist ;
  char               cmdline[CMD_LINE_LEN] ;
  FACE               *face ;
  MHT                *mht ;
  VERTEX             *v ;

  make_cmd_version_string
  (argc, argv,
   "$Id: mris_make_face_parcellation.c,v 1.2 2009/01/22 14:01:39 fischl Exp $",
   "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
          (argc, argv,
           "$Id: mris_make_face_parcellation.c,v 1.2 2009/01/22 14:01:39 fischl Exp $",
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

  if (argc < 4)
    print_help() ;

  in_fname = argv[1] ;
  ico_fname = argv[2] ;
  out_fname = argv[3] ;
  FileNamePath(out_fname, path) ;
  FileNameOnly(ico_fname, ico_name) ;

  mris = MRISread(in_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_fname) ;
  mris_ico = MRISread(ico_fname) ;
  if (!mris_ico)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, ico_fname) ;
  if (mris_ico->nfaces < 256)
    scale = 256 / mris_ico->nfaces ;
  else
    scale = 1 ;

  MRISaddCommandLine(mris, cmdline) ;

  mris->ct = CTABalloc(mris_ico->nfaces) ;
  strcpy (mris->ct->fname, ico_fname);

  printf("parcellating hemisphere into %d units\n", mris_ico->nfaces) ;
  for (fno = 0 ; fno < mris_ico->nfaces ; fno++)
  {
    int f = nint(scale * fno), m ;

    // don't let the color be too close to 0 by scaling them
    r = (nint(scale*f) / (256*256)) ;
    g = (nint(scale*f) / 256) ;
    b = (nint(scale*f) % 256) ;
    m = fno % 10 ;

    // try to avoid having adjacent triangles with similar colors
    switch (m)
    {
    default:
    case 0: r += 128 ; break ;
    case 1: g += 128 ; break ;
    case 2: b += 128 ; break ;
    case 3: r = 255-r ; break ;
    case 4: g = 255-g ; break ;
    case 5: b = 255-b ; break ;
    case 6: b += 128 ; r = 255-r ; break ;
    case 7: g += 128 ; b = 255-b ;break ;
    case 8: g += 128 ; b = 255-b ;break ; r += 128 ;
    case 9: g += 128 ; b = 255-b ;break ; r = 255-r ;
    }

    if (r < 0)
      r = 0 ;
    if (g < 0)
      g = 0 ;
    if (b < 0)
      b = 0 ;
    r = r % 256 ; g = g % 256 ; b = b %256 ;
    if (r > 255 || g > 255 || b > 255)
      DiagBreak() ;
    sprintf (mris->ct->entries[fno]->name, "%s face %d", ico_name, fno);
    mris->ct->entries[fno]->ri = r ;
    mris->ct->entries[fno]->gi = g ;
    mris->ct->entries[fno]->bi = b ;
    mris->ct->entries[fno]->ai = 255;
    
    /* Now calculate the float versions. */
    mris->ct->entries[fno]->rf =
      (float)mris->ct->entries[fno]->ri / 255.0;
    mris->ct->entries[fno]->gf =
      (float)mris->ct->entries[fno]->gi / 255.0;
    mris->ct->entries[fno]->bf =
      (float)mris->ct->entries[fno]->bi / 255.0;
    mris->ct->entries[fno]->af =
      (float)mris->ct->entries[fno]->ai / 255.0;
  }

  radius = MRISaverageRadius(mris) ;
  MRISscaleBrain(mris_ico, mris_ico, radius / mris_ico->radius) ;
  MRIScomputeMetricProperties(mris_ico) ;
  mht = MHTfillTableAtResolution(mris_ico, NULL, CURRENT_VERTICES, 1.0);
  //  mht = MHTfillTableAtResolution(mris_ico, NULL, CURRENT_VERTICES, mris_ico->avg_vertex_dist/10);

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (((vno % (mris->nvertices/10))) == 0)
      printf("%2.1f%% done\n", 100.0*(float)vno / mris->nvertices) ;
    v = &mris->vertices[vno] ;
    MHTfindClosestFaceGeneric(mht, mris, v->x, v->y, v->z, mris_ico->avg_vertex_dist/10, 
                              mris_ico->avg_vertex_dist/10, &face, &fno, &fdist) ;
    if (fno < 0)
    {
      fno = mhtBruteForceClosestFace(mris_ico, v->x, v->y, v->z, 
                                      CURRENT_VERTICES, NULL);    
      if (fno  < 0)
      {
        printf("warning: v %d not found in MHT\n", vno) ;
        continue ;
      }
    }
    CTABannotationAtIndex(mris->ct, fno, &annot);
    v->annotation = annot ;
  }

  if (Gdiag & DIAG_SHOW)
    fprintf(stderr, "writing annotation to %s\n", out_fname) ;
  MRISwriteAnnotation(mris, out_fname) ;
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
  else if (!stricmp(option, "-version")){
    print_version() ;
  } else switch (toupper(*option)) {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      printf("debugging vertex %d\n", Gdiag_no) ;
      nargs = 1 ;
      break ;
    case 'W':
      Gdiag |= DIAG_WRITE ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      print_usage() ;
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
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <input surface> <output surface>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
    "\nThis generates a parcellation based on which icosahedral face each vertex maps to\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

