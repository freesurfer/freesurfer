/**
 * @brief program for merging two parcellations into one, taking some divisions from one and some from the other
 *
 */
/*
 * Original Author: Bruce Fischl
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "timer.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "label.h"
#include "MARS_DT_Boundary.h"


int main(int argc, char *argv[]) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int merge_annotations(COLOR_TABLE *ct, MRI_SURFACE *mris1, MRI_SURFACE *mris2, MRI_SURFACE *mris) ;

const char *Progname ;
static char          fsdir[STRLEN] = "" ;

int
main(int argc, char *argv[])
{
  char          **av, *parc1, *parc2, *oname, surf_name[STRLEN], path[STRLEN],hemi[STRLEN],*cp,fname[STRLEN];
  int           ac, nargs, msec ;
  MRI_SURFACE   *mris1, *mris2 ;
  Timer then ;
  COLOR_TABLE   *ct ;

  nargs = handleVersionOption(argc, argv, "mris_merge_parcellations");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
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

  if (argc < 4)
    usage_exit() ;
  then.reset() ;

  if (strlen(fsdir) == 0)
  {
    cp = getenv("FREESURFER_HOME") ;
    if (cp == NULL)
      ErrorExit(ERROR_BADPARM, "FRESURFER_HOME must be defined in the environment") ;
    strcpy(fsdir, cp) ;
  }
  int req = snprintf(fname, STRLEN, "%s/FreeSurferColorLUT.txt", fsdir) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  ct = CTABreadASCII(fname) ;
  if (ct == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read color table from %s", Progname, fname) ;

  parc1 = argv[1] ;
  parc2 = argv[2] ;
  oname = argv[3] ;
  FileNamePath(parc1, path);
  FileNameOnly(parc1, fname) ;
  cp = strstr(fname, "h.") ;
  if (cp == NULL)
    ErrorExit(ERROR_UNSUPPORTED, "%s: could not scan hemisphere from fname %s", Progname, fname) ;
  strncpy(hemi, cp-1, 2) ; hemi[2] = 0 ;
  req = snprintf(surf_name, STRLEN, "%s/../surf/%s.orig", path, hemi) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  mris1 = MRISread(surf_name) ;
  if (mris1 == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, surf_name) ;
  if (MRISreadAnnotation(mris1, parc1) != NO_ERROR)
    ErrorExit(ERROR_BADFILE, "%s:could not open annotation %s", Progname, parc1);
  if (mris1->ct == NULL)
    ErrorExit(ERROR_BADFILE, "Annotation %s does not contain a color table", parc1);

  mris2 = MRISread(surf_name) ;
  if (mris2 == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, surf_name) ;
  if (MRISreadAnnotation(mris2, parc2) != NO_ERROR)
    ErrorExit(ERROR_BADFILE, "%s:could not open annotation %s", Progname, parc2);
  if (mris2->ct == NULL)
    ErrorExit(ERROR_BADFILE, "Annotation %s does not contain a color table", parc2);

  merge_annotations(ct, mris1, mris2, mris2) ;

  printf("writing output parcellation to %s\n", oname) ;
  MRISwriteAnnotation(mris2, oname) ;
  msec = then.milliseconds() ;
  printf("parcellation merging took %2.1f minutes\n", (float)msec/(60*1000.0f));

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
  else if (!stricmp(option, "fsdir"))
  {
    strcpy(fsdir, argv[2]) ;
    printf("using %s as FREESURFER_HOME", fsdir) ;
  }
  else switch (toupper(*option))
    {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      print_help() ;
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

static void
print_usage(void)
{
  printf("usage: %s [options] <surface> <label1> <label2>\n",
         Progname) ;
  printf("\t-a <annot name>    compute pairwise HD between all annotations\n");
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr,
          "This program computes the Hausdorff distance between two labels on "
          "the surface\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr,
          "(see the source code!)\n") ;
  exit(1) ;
}


static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

/*
  these are actually the same surface, the 1st one has Rahul's parcellation, the
  2nd one Christophe's. Grab the cingulate subdivisions from Rahul's and put them in Christophe's.
*/
static int
merge_annotations(COLOR_TABLE *ct, MRI_SURFACE *mris1, MRI_SURFACE *mris2, MRI_SURFACE *mris)
{
  int       vno, annot, s_cingulate, g_cingulate, s_pericallosal ;
  int       caudal_acc, posterior_cingulate, rostral_acc ;
  int       s_caudal_acc, s_posterior_cingulate, s_rostral_acc, s_caudal_peri, s_rostral_peri, s_posterior_peri;
  int       g_caudal_acc, g_posterior_cingulate, g_rostral_acc, filled, n ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    CTABprintASCII(ct, stdout) ;
  if (mris2->hemisphere == RIGHT_HEMISPHERE)
  {
    caudal_acc = CTABentryNameToAnnotation("ctx-rh-caudalanteriorcingulate", ct) ;
    if (caudal_acc < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-caudalanteriorcingulate") ;
    posterior_cingulate = CTABentryNameToAnnotation("ctx-rh-posteriorcingulate", ct) ;
    if (posterior_cingulate < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-posteriorcingulate") ;
    rostral_acc = CTABentryNameToAnnotation("ctx-rh-rostralanteriorcingulate", ct) ;
    if (rostral_acc < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-rostralanteriorcingulate") ;

    s_caudal_acc = CTABentryNameToAnnotation("ctx-rh-S_cingulate-caudal_ACC", ct) ;
    if (s_caudal_acc < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-S_cingulate-caudal_ACC") ;
    s_rostral_acc = CTABentryNameToAnnotation("ctx-rh-S_cingulate-rostral_ACC", ct);
    if (s_rostral_acc < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-S_cingulate-rostral_ACC") ;
    s_posterior_cingulate = CTABentryNameToAnnotation("ctx-rh-S_cingulate-posterior", ct) ;
    if (s_posterior_cingulate < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-S_cingulate-posterior") ;

    g_caudal_acc = CTABentryNameToAnnotation("ctx-rh-G_cingulate-caudal_ACC", ct) ;
    if (g_caudal_acc < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-G_cingulate-caudal_ACC") ;
    g_rostral_acc = CTABentryNameToAnnotation("ctx-rh-G_cingulate-rostral_ACC", ct);
    if (g_rostral_acc < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-G_cingulate-rostral_ACC") ;
    g_posterior_cingulate = CTABentryNameToAnnotation("ctx-rh-G_cingulate-posterior", ct) ;
    if (g_posterior_cingulate < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-G_cingulate-posterior") ;

    s_caudal_peri = CTABentryNameToAnnotation("ctx-rh-S_pericallosal-caudal", ct) ;
    if (s_caudal_peri < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-S_pericallosal-caudal") ;
    s_rostral_peri = CTABentryNameToAnnotation("ctx-rh-S_pericallosal-rostral", ct);
    if (s_rostral_peri < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-S_pericallosal-rostral") ;
    s_posterior_peri = CTABentryNameToAnnotation("ctx-rh-S_pericallosal-posterior", ct) ;
    if (s_posterior_peri < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-S_pericallosal-posterior") ;

    s_cingulate = CTABentryNameToAnnotation("ctx-rh-S_cingulate-Main_part_and_Intracingulate", ct) ;
    if (s_cingulate < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-S_cingulate-Main_part_and_Intracingulate") ;
    g_cingulate = CTABentryNameToAnnotation("ctx-rh-G_cingulate-Main_part", ct) ;
    if (g_cingulate < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-G_cingulate-Main_part") ;
    s_pericallosal = CTABentryNameToAnnotation("ctx-rh-S_pericallosal", ct) ;
    if (s_pericallosal < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-rh-S_pericallosal") ;
  }
  else // left hemi
  {
    caudal_acc = CTABentryNameToAnnotation("ctx-lh-caudalanteriorcingulate", ct) ;
    if (caudal_acc < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-caudalanteriorcingulate") ;
    posterior_cingulate = CTABentryNameToAnnotation("ctx-lh-posteriorcingulate", ct) ;
    if (posterior_cingulate < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-posteriorcingulate") ;
    rostral_acc = CTABentryNameToAnnotation("ctx-lh-rostralanteriorcingulate", ct) ;
    if (rostral_acc < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-rostralanteriorcingulate") ;

    s_caudal_acc = CTABentryNameToAnnotation("ctx-lh-S_cingulate-caudal_ACC", ct) ;
    if (s_caudal_acc < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-S_cingulate-caudal_ACC") ;
    s_rostral_acc = CTABentryNameToAnnotation("ctx-lh-S_cingulate-rostral_ACC", ct);
    if (s_rostral_acc < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-S_cingulate-rostral_ACC") ;
    s_posterior_cingulate = CTABentryNameToAnnotation("ctx-lh-S_cingulate-posterior", ct) ;
    if (s_posterior_cingulate < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-S_cingulate-posterior") ;

    g_caudal_acc = CTABentryNameToAnnotation("ctx-lh-G_cingulate-caudal_ACC", ct) ;
    if (g_caudal_acc < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-G_cingulate-caudal_ACC") ;
    g_rostral_acc = CTABentryNameToAnnotation("ctx-lh-G_cingulate-rostral_ACC", ct) ;
    if (g_rostral_acc < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-G_cingulate-rostral_ACC") ;
    g_posterior_cingulate = CTABentryNameToAnnotation("ctx-lh-G_cingulate-posterior", ct);
    if (g_posterior_cingulate < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-G_cingulate-posterior") ;

    s_caudal_peri = CTABentryNameToAnnotation("ctx-lh-S_pericallosal-caudal", ct) ;
    if (s_caudal_peri < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-S_pericallosal-caudal") ;
    s_rostral_peri = CTABentryNameToAnnotation("ctx-lh-S_pericallosal-rostral", ct);
    if (s_rostral_peri < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-S_pericallosal-rostral") ;
    s_posterior_peri = CTABentryNameToAnnotation("ctx-lh-S_pericallosal-posterior", ct) ;
    if (s_posterior_peri < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-S_pericallosal-posterior") ;

    s_cingulate = CTABentryNameToAnnotation("ctx-lh-S_cingulate-Main_part_and_Intracingulate", ct) ;
    if (s_cingulate < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-S_cingulate-Main_part_and_Intracingulate") ;
    g_cingulate = CTABentryNameToAnnotation("ctx-lh-G_cingulate-Main_part", ct) ;
    if (g_cingulate < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-G_cingulate-Main_part") ;
    s_pericallosal = CTABentryNameToAnnotation("ctx-lh-S_pericallosal", ct) ;
    if (s_pericallosal < 0)
      ErrorExit(ERROR_UNSUPPORTED, "Cannot find annotation %s in LUT", "ctx-lh-S_pericallosal") ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX * const v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    annot = mris2->vertices[vno].annotation ;
    if  (mris1->vertices[vno].annotation == caudal_acc)
    {
      if (annot == s_cingulate /* || annot == s_pericallosal*/)
        annot = s_caudal_acc ; 
      else if (annot == g_cingulate)
        annot = g_caudal_acc ; 
    }
    else if  (mris1->vertices[vno].annotation == rostral_acc)
    {
      if (annot == s_cingulate /* || annot == s_pericallosal*/)
        annot = s_rostral_acc ;
      else if (annot == g_cingulate)
        annot = g_rostral_acc ; 
    }
    else if (mris1->vertices[vno].annotation == posterior_cingulate)
    {
      if (annot == s_cingulate /* || annot == s_pericallosal*/)
        annot = s_posterior_cingulate ; 
      else if (annot == g_cingulate)
        annot = g_posterior_cingulate ; 
    }
    if (annot < 0)
      DiagBreak() ;
    v->annotation = annot ;
  }

  // now diffuse cing sulcal labels into the rest of the cingulate sulcus in Christophe's labels
  if (Gdiag_no >= 0)
  {
    int index ;
    CTABfindAnnotation(mris1->ct, mris1->vertices[Gdiag_no].annotation, &index);
    printf("mris1->annot[%d] = %s (%d)\n", Gdiag_no, mris1->ct->entries[index]->name, mris1->vertices[Gdiag_no].annotation) ;
    CTABfindAnnotation(mris2->ct, mris2->vertices[Gdiag_no].annotation, &index);
    printf("mris2->annot[%d] = %s (%d)\n", Gdiag_no, mris2->ct->entries[index]->name, mris2->vertices[Gdiag_no].annotation) ;
  }

  MRISclearMarks(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    VERTEX * const v = &mris->vertices[vno] ;
    if (v->annotation == s_cingulate)
      v->marked = 1 ;
  }

  do
  {
    filled = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->marked == 1) 
      {
        for (n = 0 ; n < vt->vnum ; n++)
        {
          VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
          if (vn->marked < 0)
            continue ;
          
          if (vn->annotation == s_posterior_cingulate)
            v->annotation = s_posterior_cingulate ;
          else if (vn->annotation == s_caudal_acc)
            v->annotation = s_caudal_acc ;
          else if (vn->annotation == s_rostral_acc)
            v->annotation = s_rostral_acc ;
          else
            continue ;

          v->marked = -1 ;     // marked in this cycle - don't use it until done with this iter
          filled++ ;
          break ;
        }
      }
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      VERTEX * const v = &mris->vertices[vno] ;
      if (v->marked == -1)
        v->marked = 0 ;
    }
    printf("%d vertices filled\n", filled) ;
  } while (filled > 0) ;
  mris->ct = ct ;

  // split up the pericallosal sulcus based on cingulate labels
  MRISclearMarks(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (vno == Gdiag_no)
      DiagBreak() ;
    VERTEX * const v = &mris->vertices[vno] ;
    if (v->annotation == s_pericallosal)
      v->marked = 1 ;
  }
  
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->marked == 1) 
    {
      for (n = 0 ; n < vt->vnum ; n++)
      {
        VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
        if (vn->marked < 0)
          continue ;
        
        if (vn->annotation == g_posterior_cingulate)
          v->annotation = s_posterior_peri ;
        else if (vn->annotation == g_caudal_acc)
          v->annotation = s_caudal_peri ;
        else if (vn->annotation == g_rostral_acc)
          v->annotation = s_rostral_peri ;
        else if (vn->annotation == g_cingulate)
          v->annotation = s_pericallosal ;
        else
          continue ;
        
        v->marked = -1 ;     // marked in this cycle - don't use it until done with this iter
        break ;
      }
    }
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX * const v = &mris->vertices[vno] ;
    if (v->marked == -1)
      v->marked = 0 ;
  }
  do
  {
    filled = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->marked == 1) 
      {
        for (n = 0 ; n < vt->vnum ; n++)
        {
          VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
          if (vn->marked < 0 || vn->marked == 1)
            continue ;
          
          if (vn->annotation == s_posterior_peri)
            v->annotation = s_posterior_peri ;
          else if (vn->annotation == s_caudal_peri)
            v->annotation = s_caudal_peri ;
          else if (vn->annotation == s_rostral_peri)
            v->annotation = s_rostral_peri ;
          else if (vn->annotation == s_pericallosal)
            v->annotation = s_pericallosal ;
          else
            continue ;
          
          v->marked = -1 ;     // marked in this cycle - don't use it until done with this iter
          filled++ ;
          break ;
        }
      }
    }
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      VERTEX * const v = &mris->vertices[vno] ;
      if (v->marked == -1)
        v->marked = 0 ;
    }
    printf("%d vertices filled\n", filled) ;
  } while (filled > 0) ;

  // now remove the ctx-?h- from the parcellation unit names
  {
    int  i, max_i = 0 ;
    CTE  *cte ;
    char buf[STRLEN] ;

    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      int index ;
      CTABfindAnnotation(mris->ct, mris->vertices[vno].annotation, &index);
      if (index > max_i)
        max_i = index ;
    }
      

    for (i = 0 ; i < mris->ct->nentries ; i++)
    {
      cte = mris->ct->entries[i] ;
      if (cte == NULL)
        continue ;
      if ((strncmp(cte->name, "ctx-lh-", 7) == 0) ||
          (strncmp(cte->name, "ctx-rh-", 7) == 0))
      {
        strcpy(buf, cte->name) ;
        strcpy(cte->name, buf+7) ;
      }
    }
    mris->ct->nentries = max_i ;
  }
  return(NO_ERROR) ;
}

