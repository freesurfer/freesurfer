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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "mrisurf.h"
#include "mrishash.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;
static const char *orig_surf_name = "sphere.reg" ;
static const char *inf_surf_name = "inflated" ;
int MRISmapCuts(MRI_SURFACE *mris_in, MRI_SURFACE *mris_out) ;

static int dilate = 0 ;

int
main(int argc, char *argv[]) {
  char         **av, in_surf_fname[STRLEN], *in_patch_fname, *out_patch_fname, hemi[STRLEN] ;
  int          ac, nargs;
  char         path[STRLEN], out_surf_fname[STRLEN], *cp ;
  int          msec, minutes, seconds ;
  Timer start ;
  MRI_SURFACE  *mris_in, *mris_out ;

  nargs = handleVersionOption(argc, argv, "mris_map_cuts");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;


  in_patch_fname = argv[1] ;
  out_patch_fname = argv[2] ;
  FileNamePath(in_patch_fname, path) ;
  cp = strrchr(in_patch_fname, '/') ;
  if (!cp)
    cp = in_patch_fname ;
  cp = strchr(cp, '.') ;
  if (cp)
  {
    strncpy(hemi, cp-2, 2) ;
    hemi[2] = 0 ;
  }
  else
    strcpy(hemi, "lh") ;
  int req = snprintf(in_surf_fname, STRLEN, "%s/%s.%s", path, hemi, orig_surf_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  FileNamePath(out_patch_fname, path) ;
  cp = strrchr(out_patch_fname, '/') ;
  if (!cp)
    cp = out_patch_fname ;
  cp = strchr(cp, '.') ;
  if (cp)
  {
    strncpy(hemi, cp-2, 2) ;
    hemi[2] = 0 ;
  }
  else
    strcpy(hemi, "lh") ;
  req = snprintf(out_surf_fname, STRLEN, "%s/%s.%s", path, hemi, orig_surf_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  mris_in = MRISread(in_surf_fname) ;
  mris_out = MRISread(out_surf_fname) ;
  MRISsaveVertexPositions(mris_in, CANONICAL_VERTICES) ;
  MRISsaveVertexPositions(mris_out, CANONICAL_VERTICES) ;
  if (MRISreadVertexPositions(mris_out, inf_surf_name)  != NO_ERROR)
    ErrorExit(ERROR_BADPARM, "%s: could not inflated surface %s",
              Progname, inf_surf_name) ;

  if (MRISreadPatch(mris_in, in_patch_fname) != NO_ERROR)
    ErrorExit(ERROR_BADPARM, "%s: could not read patch file %s",
              Progname, in_patch_fname) ;
  MRISmapCuts(mris_in, mris_out) ;
  if (dilate)
  {
    printf("dilating patch %d times\n", dilate) ;
    MRISdilateRipped(mris_out, dilate) ;
    printf("%d valid vertices (%2.1f %% of total)\n",
           MRISvalidVertices(mris_out), 100.0*MRISvalidVertices(mris_out)/mris_out->nvertices) ;
  }

  printf("writing output to %s\n", out_patch_fname) ;
  MRISwritePatch(mris_out, out_patch_fname) ;
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "cut mapping took %d minutes"
          " and %d seconds.\n", minutes, seconds) ;
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
  case 'D':
    dilate = atoi(argv[2]) ;
    nargs = 1 ;
    printf("dilating output cuts %d times\n", dilate) ;
    break ;
  case '?':
  case 'U':
    usage_exit(0) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s [options] <input patch> <output patch>",
         Progname) ;
#if 0
  printf(
    "\tf <f low> <f hi> - apply specified filter (not implemented yet)\n"
  );
  printf("\tn - noise-sensitivity normalize inverse (default=1)") ;
#endif
  exit(code) ;
}

/*
  expects the two surfaces to have the CANONICAL_VERTICES field set to the 
  sphere.reg positions.
*/
int
MRISmapCuts(MRI_SURFACE *mris_in, MRI_SURFACE *mris_out)
{
  MHT    *mht_in, *mht_out ;
  int    vno_out /*, vno_in, fno_in, fno_out, n, ripflag*/ ;
#if 0
  FACE   *f_in, *f_out ;
#endif

  MRISstoreRipFlags(mris_in) ;
  MRISunrip(mris_in) ;

  mht_in  = MHTcreateVertexTable(mris_in,  CANONICAL_VERTICES) ;
  mht_out = MHTcreateVertexTable(mris_out, CANONICAL_VERTICES) ;

#if 0
  for (vno_in = 0 ; vno_in < mris_in->nvertices ; vno_in++)
  {
    if (vno_in == Gdiag_no)
      DiagBreak() ;
    v_in = &mris_in->vertices[vno_in] ;
    if (v_in->oripflag == 0)
      continue ;
    v_out = MHTfindClosestVertexInTable(mht_out, mris_out, v_in->cx, v_in->cy, v_in->cz, 1) ;
    if (v_out == NULL)
      DiagBreak() ;
    else
      v_out->ripflag = 1 ;
  }

  for (vno_out = 0 ; vno_out < mris_out->nvertices ; vno_out++)
  {
    if (vno_out == Gdiag_no)
      DiagBreak() ;
    v_out = &mris_out->vertices[vno_out] ;
    v_in = MHTfindClosestVertexInTable(mht_in, mris_in, v_out->cx, v_out->cy, v_out->cz, 1) ;
    if (v_in == NULL)
      DiagBreak() ;
    else if (v_in->oripflag)
      v_out->ripflag = 1 ;
  }

  for (fno_in = 0 ; fno_in < mris_in->nfaces ; fno_in++)
  {
    f_in = &mris_in->faces[fno_in] ;
    if (f_in->ripflag == 0)
      continue ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      vno_in = f_in->v[n] ;
      if (vno_in == Gdiag_no)
        DiagBreak() ;
      v_in = &mris_in->vertices[vno_in] ;
      v_out = MHTfindClosestVertexInTable(mht_out, mris_out, v_in->cx, v_in->cy, v_in->cz, 1) ;
      if (v_out == NULL)
        DiagBreak() ;
      else
        v_out->ripflag = 1 ;
    }
  }

  for (fno_out = 0 ; fno_out < mris_out->nfaces ; fno_out++)
  {
    f_out = &mris_out->faces[fno_out] ;
    ripflag = 0 ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      vno_out = f_out->v[n] ;
      if (vno_out == Gdiag_no)
        DiagBreak() ;
      v_out = &mris_out->vertices[vno_out] ;
      v_in = MHTfindClosestVertexInTable(mht_in, mris_in, v_out->cx, v_out->cy, v_out->cz, 1) ;
      if (v_in == NULL)
        DiagBreak() ;
      else
        if (v_in->ripflag)
          ripflag = 1 ;
    }
    f_out->ripflag = ripflag ;
  }

  for (fno_out = 0 ; fno_out < mris_out->nfaces ; fno_out++)
  {
    double cx, cy, cz ;
    int    n ;

    f_out = &mris_out->faces[fno_out] ;
    ripflag = 0 ;
    
    cx = cy = cz = 0.0 ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
    {
      vno_out = f_out->v[n] ;
      if (vno_out == Gdiag_no)
        DiagBreak() ;
      v_out = &mris_out->vertices[vno_out] ;
      cx += v_out->cx ; cy += v_out->cy ; cz += v_out->cz ;

    }
    cx /= VERTICES_PER_FACE ;
    cy /= VERTICES_PER_FACE ;
    cz /= VERTICES_PER_FACE ;
    v_in = MHTfindClosestVertexInTable(mht_in, mris_in, cx, cy, cz, 1) ;
    if (v_in == NULL)
      DiagBreak() ;
    else
      if (v_in->ripflag)
      {
        f_out->ripflag = 1 ;
        for (n = 0 ; n < VERTICES_PER_FACE ; n++)
          mris_out->vertices[f_out->v[n]].ripflag = 1 ;
      }
  }
#endif

#define STEP_SIZE 0.05
  for (vno_out = 0 ; vno_out < mris_out->nvertices ; vno_out++)
  {
    double d, dx, dy, dz ;
    int    n ;

    VERTEX_TOPOLOGY const * const v_outt = &mris_out->vertices_topology[vno_out] ;
    VERTEX                * const v_out  = &mris_out->vertices         [vno_out] ;
    
    for (n = 0 ; n < v_outt->vnum ; n++)
    {
      VERTEX const * const vn = &mris_out->vertices[v_outt->v[n]] ;
      dx = vn->cx - v_out->cx ; 
      dy = vn->cy - v_out->cy ; 
      dz = vn->cz - v_out->cz ; 
      for (d = 0.0 ; d <= 1.0 ; d+= STEP_SIZE)
      {
        VERTEX const * const v_in = MHTfindClosestVertexInTable(mht_in, mris_in, 
                                           v_out->cx+dx*d, 
                                           v_out->cy+dy*d, 
                                           v_out->cz+dz*d, 
                                           1) ;
        if (v_in == NULL)
          DiagBreak() ;
        else if (v_in->oripflag)
        {
          v_out->ripflag = 1 ;
          break ;
        }
      }
    }
  }

  MRISsetRipInFacesWithRippedVertices(mris_out) ;
  MRISrestoreRipFlags(mris_in) ;

  MHTfree(&mht_in) ; MHTfree(&mht_out) ;
  return(NO_ERROR) ;
}

