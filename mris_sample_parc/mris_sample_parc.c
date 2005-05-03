
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
#include "annotation.h"
#include "version.h"

static char vcid[] = "$Id: mris_sample_parc.c,v 1.12 2005/05/03 20:24:49 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  translate_indices_to_annotations(MRI_SURFACE *mris, char *translation_fname) ;


char *Progname ;
static char *color_table_fname = NULL ;
static int mode_filter = 0 ;
static char *surf_name = WHITE_MATTER_NAME ;
static char *thickness_name = "thickness" ;
static char sdir[STRLEN] ;
static char *translation_fname = "cma_parcellation_colors.txt" ;
static int wsize = 7 ;
static int unknown_label = -1 ;
static int fix_topology = 1 ;
static int fix_label_topology(MRI_SURFACE *mris) ;
static int resegment_label(MRI_SURFACE *mris, LABEL *segment) ;
static float proj_mm = 0.0 ;

#define MAX_LABEL 1000

int
main(int argc, char *argv[])
{
  char          **av, *hemi, *subject_name, *cp, fname[STRLEN], *parc_name, *annot_name ;
  int           ac, nargs, vno ;
  MRI_SURFACE   *mris ;
  MRI           *mri_parc ;
  VERTEX        *v ;
  double        d ;
  Real          x, y, z, xw, yw, zw ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_sample_parc.c,v 1.12 2005/05/03 20:24:49 fischl Exp $", "$Name:  $");
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

  if (argc < 4)
    usage_exit() ;

  subject_name = argv[1] ;
  hemi = argv[2] ; 
  parc_name = argv[3] ;
  annot_name = argv[4] ;

  if (strlen(sdir) == 0)  /* if not specified explicitly as option */
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, 
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  sprintf(fname, "%s/%s/mri/%s", sdir, subject_name, parc_name) ;
  printf("reading parcellation volume from %s...\n", fname) ;
  mri_parc = MRIread(fname) ;
  if (!mri_parc)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject_name, hemi, surf_name) ;
  printf("reading input surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;

	if (FZERO(proj_mm))
	{
		if (MRISreadCurvatureFile(mris, thickness_name) != NO_ERROR)
			ErrorExit(ERROR_NOFILE, "%s: could not read thickness file %s",
								Progname, thickness_name) ;
	}

	if (color_table_fname)
		mris->ct = CTABread(color_table_fname) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

		if (!FZERO(proj_mm))
			d = proj_mm ;
		else
			d = v->curv*.5 ;  /* halfway out */
    x = v->x+d*v->nx ; y = v->y+d*v->ny ; z = v->z+d*v->nz ;
		MRIworldToVoxel(mri_parc, x, y, z, &xw, &yw, &zw) ;
    MRIsurfaceRASToVoxel(mri_parc, x, y, z, &xw, &yw, &zw) ;
		v->annotation = v->val = 
			MRIfindNearestNonzero(mri_parc, wsize, xw, yw, zw) ;
  }
  if (unknown_label >= 0)
  {
    LABEL **labels, *label ;
    int   nlabels, i, biggest_label, most_vertices, nzero ;

#define TMP_LABEL 1000
    for (nzero = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->annotation == 0)
      {
        v->annotation = TMP_LABEL;
        nzero++ ;
      }
    }
    printf("%d unknown vertices found\n", nzero) ;
    MRISsegmentAnnotated(mris, &labels, &nlabels, 10) ;
    most_vertices = 0 ; biggest_label = -1 ;
    for (i = 0 ; i < nlabels ; i++)
    {
      label = labels[i] ;
      if (mris->vertices[label->lv[0].vno].annotation == TMP_LABEL)
      {
        if (label->n_points > most_vertices)
        {
          biggest_label = i ; most_vertices = label->n_points ;
        }
      }
    }
    if (biggest_label >= 0)
    {
      label = labels[biggest_label] ;
      printf("replacing label # %d with %d vertices (vno=%d) with label %d\n",
             biggest_label, label->n_points, label->lv[0].vno, unknown_label) ;
      for (i = 0 ; i < label->n_points ; i++)
      {
        v = &mris->vertices[label->lv[i].vno] ;
        v->annotation = v->val = unknown_label ;
      }
    }
    for (nzero = vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->annotation == TMP_LABEL)
      {
        v->annotation = 0;
        nzero++ ;
      }
    }
    printf("after replacement, %d unknown vertices found\n", nzero) ;
    MRISmodeFilterZeroVals(mris) ;  /* get rid of the rest of the unknowns by mode filtering */
    for (i = 0 ; i < nlabels ; i++)
      LabelFree(&labels[i]) ;
    free(labels) ;
  }

  if (fix_topology)
    fix_label_topology(mris) ;

  if (mode_filter)
  {
    printf("mode filtering sample labels...\n") ;
#if 0
    MRISmodeFilterZeroVals(mris) ;
#else
    MRISmodeFilterVals(mris, mode_filter) ;
#endif
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->annotation = v->val ;
    }
  }

  /* this will fill in the v->annotation field from the v->val ones */
  translate_indices_to_annotations(mris, translation_fname) ;

  printf("writing annotation to %s...\n", annot_name) ;
  MRISwriteAnnotation(mris, annot_name) ; 
	/*  MRISreadAnnotation(mris, fname) ;*/
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
  else if (!stricmp(option, "sdir"))
  {
    strcpy(sdir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as SUBJECTS_DIR\n", sdir) ;
  }
  else if (!stricmp(option, "surf"))
  {
		surf_name = argv[2] ;
    nargs = 1 ;
    printf("using %s as surface name\n", surf_name) ;
  }
  else if (!stricmp(option, "proj"))
  {
		proj_mm = atof(argv[2]) ;
    nargs = 1 ;
    printf("projecting %2.2f mm along surface normal\n", proj_mm) ;
  }
  else if (!stricmp(option, "file"))
  {
    translation_fname = argv[2] ;
    nargs = 1 ;
    printf("using %s as translation fname\n", sdir) ;
  }
  else if (!stricmp(option, "ct"))
  {
		color_table_fname = argv[2] ;
    nargs = 1 ;
    printf("embedding color table %s into output annot file\n", color_table_fname) ;
		translation_fname = color_table_fname ;
  }
  else switch (toupper(*option))
  {
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'F':
    mode_filter = atoi(argv[2]) ;
    nargs = 1 ;
    printf("applying mode filter %d times to parcellation\n", mode_filter) ; 
    break ;
  case 'W':
    wsize = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using window size=%d for sampling\n", wsize) ;
    break ;
  case 'T':
    thickness_name = argv[2] ;
    nargs = 1 ;
    printf("using thickness file %s\n", thickness_name) ;
    break ;
  case 'U':
    unknown_label = atoi(argv[2]) ;
    printf("changing largest connected unknown region to label %d\n", unknown_label) ;
    nargs = 1 ;
    break ;
  case '?':
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
usage_exit(void)
{
  print_help() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <subject name> <hemi> <parc name> <output annot>\n", 
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
          "\nThis program samples a volumetric parcellation onto a surface\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr,
          "-t <thickness file>          - use specified file for computing "
          "thickness statistics\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

#if 0
static int
translate_indices_to_annotations(MRI_SURFACE *mris, char *translation_fname)
{
  char   fname[STRLEN], *cp, line[STRLEN], name[STRLEN], **names ;
  FILE   *fp ;
  int    nlines, i, vno, *r, *g, *b, *indices ;
  VERTEX *v ;

  cp = getenv("FREESURFER_HOME") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "%s: FREESURFER_HOME not defined in environment", Progname) ;
  sprintf(fname, "%s/%s", cp, translation_fname) ;
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "%s: could not read translation file %s", Progname, fname) ;

  nlines = 0 ;
  cp = fgetl(line, STRLEN-1, fp) ;
  while (cp != NULL)
  {
    nlines++ ;
    cp = fgetl(line, STRLEN-1, fp) ;
  }
  printf("%d lines found in file\n", nlines) ;
  rewind(fp) ;

  r = (int *)calloc(nlines, sizeof(int)) ;
  g = (int *)calloc(nlines, sizeof(int)) ;
  b = (int *)calloc(nlines, sizeof(int)) ;
  indices = (int *)calloc(nlines, sizeof(int)) ;
  names = (char **)calloc(nlines, sizeof(char *)) ;
  if (!r || !g || !b || !indices || !names)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d-len internal buffers",
              Progname, nlines) ;

  for (i = 0 ; i < nlines ; i++)
  {
    cp = fgetl(line, STRLEN-1, fp) ;
    sscanf(cp, "%d %s %d %d %d %*d",
           &indices[i], name, &r[i], &g[i], &b[i]) ;
    names[i] = (char *)calloc(strlen(name)+1, sizeof(char)) ;
    strcpy(names[i], name) ;
    printf("parsing parcellation unit %s...\n", name) ;
  }
  fclose(fp) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
  }
  return(NO_ERROR) ;
}
#else
static int
translate_indices_to_annotations(MRI_SURFACE *mris, char *translation_fname)
{
  int    vno ;
  VERTEX *v ;

  read_named_annotation_table(translation_fname) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->annotation = index_to_annotation((int)v->val) ;
  }
  return(NO_ERROR) ;
}
#endif
static int
fix_label_topology(MRI_SURFACE *mris)
{
  int    i, vno, nsegments, most_vertices, max_label, label, j, iter, nchanged=0,
         max_index ;
  LABEL **segments ;
  VERTEX *v ;

  for (max_label = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->annotation > max_label)
      max_label = v->annotation ;
  }

  iter = 0 ;
  do
  {
    nchanged = 0 ;
    MRISsegmentAnnotated(mris, &segments, &nsegments, 0) ;
    
    for (i = 0 ; i < nsegments ; i++)
    {
      if (segments[i] == NULL)
        continue ;
      most_vertices = segments[i]->n_points ; max_index = i ;
      label = mris->vertices[segments[i]->lv[0].vno].annotation ;
      
      /* find label with most vertices */
      for (j = 0 ; j < nsegments ; j++)
      {
        if (j == i || segments[j] == NULL)
          continue ;
        if (mris->vertices[segments[j]->lv[0].vno].annotation != label)
          continue ;
        if (segments[j]->n_points > most_vertices)
        {
          most_vertices = segments[j]->n_points ;
          max_index = j ;
        }
      }
      
      /* resegment all others connected-components with this label */
      for (j = 0 ; j < nsegments ; j++)
      {
        if (j == max_index || segments[j] == NULL)
          continue ;
        if (mris->vertices[segments[j]->lv[0].vno].annotation != label)
          continue ;
        resegment_label(mris, segments[j]) ;
        nchanged++ ;
        LabelFree(&segments[j]) ;
      }
    }
    for (i = 0 ; i < nsegments ; i++)
      if (segments[i])
        LabelFree(&segments[i]) ;

    free(segments) ;
    printf("pass %d: %d segments changed\n", iter+1, nchanged) ;
  } while (nchanged > 0 && iter++ < 10) ;

  MRISclearMarks(mris) ;
  return(NO_ERROR) ;
}

static int
resegment_label(MRI_SURFACE *mris, LABEL *segment)
{
  int    histo[MAX_LABEL], i, n, vno, ino, index, max_histo, max_index, nchanged, lno, label ;
  VERTEX *v, *vn ;

  label = mris->vertices[segment->lv[0].vno].annotation ;
  for (ino  = 0 ; ino < 100 ; ino++)
  {
    nchanged = 0 ;
    for (lno = 0 ; lno < segment->n_points ; lno++)
    {
      vno = segment->lv[lno].vno ;
      v = &mris->vertices[vno] ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->val != label || v->ripflag)
        continue ;   /* already changed */

      memset(histo, 0, sizeof(histo)) ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        index = (int)nint(vn->val) ;
        if (index < 0 || index > MAX_LABEL)
          continue ;
        if (vn->val != label)  /* don't care about same label */
          histo[index]++ ;
      }
      max_histo = histo[0] ; max_index = -1 ;
      for (i = 0 ; i < MAX_LABEL ; i++)
      {
        if (histo[i] > max_histo)
        {
          max_histo = histo[i] ;
          max_index = i ;
        }
      }
      if (max_index >= 0)
        v->valbak = max_index ;
    }
    for (lno = 0 ; lno < segment->n_points ; lno++)
    {
      vno = segment->lv[lno].vno ;
      v = &mris->vertices[vno] ;
      if (v->ripflag || v->val != label)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->val != v->valbak)
        nchanged++ ;
      v->val = v->annotation = v->valbak ;
    }

    /*    printf("iter %d: %d changed\n", ino, nchanged) ;*/
    if (!nchanged)
      break ;
  }
  return(NO_ERROR) ;
}
