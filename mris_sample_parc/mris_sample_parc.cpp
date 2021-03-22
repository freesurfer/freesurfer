/**
 * @brief program for sampling a volumetric parcellation onto a surface model
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
#include "annotation.h"
#include "version.h"
#include "cma.h"
#include "gca.h"
#include "mrishash.h"


int main(int argc, char *argv[]) ;

static int  replace_vertices_with_label(MRI_SURFACE *mris,
                                        MRI *mri,
                                        int label,
                                        double proj_mm);
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  translate_indices_to_annotations(MRI_SURFACE *mris,
    const char *translation_fname) ;
static int  fix_label_topology(MRI_SURFACE *mris, int nvertices) ;
static int  resegment_label(MRI_SURFACE *mris, LABEL *segment) ;
int MRIsampleParcellationToSurface(MRI_SURFACE *mris, MRI *mri_parc) ;

const char *Progname ;
static int   avgs = 0 ;
static int nclose = 0 ;
static char *color_table_fname = NULL ;
static LABEL *cortex_label = NULL ; // limit surface area calc to cortex.label
static int   mode_filter = 0 ;
static const char *surf_name = WHITE_MATTER_NAME ;
static const char *thickness_name = "thickness" ;
static char  sdir[STRLEN] ;
static const char *translation_fname = "cma_parcellation_colors.txt" ;
static int   wsize = 7 ;
static int   unknown_label = -1 ;
static int   fix_topology = -1 ;  // < 0 means do all
static float proj_mm = 0.0 ;
static float proj_frac = 0.5 ;
static int   replace_label = 0;
#define MAX_TRANS 100
static int   ntrans = 0 ;
static int   trans_in[MAX_TRANS] ;
static int   trans_out[MAX_TRANS] ;

static int sample_from_vol_to_surf = 0 ;
static char  *mask_fname = NULL ;
static int   mask_val ;
static int label_index = -1 ;

int
main(int argc, char *argv[]) {
  char          **av, *hemi, *subject_name, *cp, fname[STRLEN];
  char          *parc_name, *annot_name ;
  int           ac, nargs, vno, i ;
  MRI_SURFACE   *mris ;
  MRI           *mri_parc ;
  VERTEX        *v ;
  double        d ;
  double        x, y, z, xw, yw, zw ;

  nargs = handleVersionOption(argc, argv, "mris_sample_parc");
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

  if (parc_name[0] == '/')  // full path specified
    strcpy(fname, parc_name) ;
  else
    sprintf(fname, "%s/%s/mri/%s", sdir, subject_name, parc_name) ;
  printf("reading parcellation volume from %s...\n", fname) ;
  mri_parc = MRIread(fname) ;
  if (!mri_parc)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

  if (mask_fname) {
    MRI *mri_mask, *mri_tmp ;

    mri_tmp = MRIread(mask_fname) ;
    if (mri_tmp == NULL)
      ErrorExit(ERROR_BADPARM, "%s: could not load mask volume %s", Progname, mask_fname) ;
    mri_mask = MRIclone(mri_tmp, NULL) ;
    MRIcopyLabel(mri_tmp, mri_mask, mask_val) ;
    MRIdilate(mri_mask, mri_mask) ;
    MRIdilate(mri_mask, mri_mask) ;
    MRIdilate(mri_mask, mri_mask) ;
    MRIdilate(mri_mask, mri_mask) ;
    MRIfree(&mri_tmp) ;
    mri_tmp = MRIclone(mri_parc, NULL) ;
    MRIcopyLabeledVoxels(mri_parc, mri_mask, mri_tmp, mask_val) ;
    MRIfree(&mri_parc) ;
    mri_parc = mri_tmp ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_parc, "p.mgz") ;
    MRIfree(&mri_mask) ;
  }

  for (i = 0 ; i < ntrans ; i++) {
    MRIreplaceValues(mri_parc, mri_parc, trans_in[i], trans_out[i]) ;
  }
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject_name, hemi, surf_name) ;
  printf("reading input surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  if (avgs > 0)
    MRISaverageVertexPositions(mris, avgs) ;

  if (FZERO(proj_mm)) {
    if (MRISreadCurvatureFile(mris, thickness_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read thickness file %s",
                Progname, thickness_name) ;
  }

  if (color_table_fname) {
    mris->ct = CTABreadASCII(color_table_fname) ;
    if (mris->ct == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read color file %s",
                Progname, color_table_fname) ;
  }

  if (sample_from_vol_to_surf) // sample from volume to surface */
  {
    MRIsampleParcellationToSurface(mris, mri_parc) ;
    if (cortex_label)
      LabelMaskSurface(cortex_label,mris);
  } else  /* sample from surface to volume */
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;

      if (!FZERO(proj_mm))
        d = proj_mm ;
      else
        d = v->curv*proj_frac ;  /* halfway out */
      x = v->x+d*v->nx ;
      y = v->y+d*v->ny ;
      z = v->z+d*v->nz ;
      MRIsurfaceRASToVoxel(mri_parc, x, y, z, &xw, &yw, &zw) ;
      v->annotation = v->val =
                        MRIfindNearestNonzero(mri_parc, wsize, xw, yw, zw, ((float)wsize-1)/2) ;
      if (v->val == 0xffffffff)
        DiagBreak() ;
    }
  }
  if (replace_label)
    replace_vertices_with_label(mris, mri_parc, replace_label, proj_mm);
  if (unknown_label >= 0) {
    LABEL **labels, *label ;
    int   nlabels, i, biggest_label, most_vertices, nzero ;

#define TMP_LABEL 1000
    for (nzero = vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->annotation == 0) {
        v->annotation = TMP_LABEL;
        nzero++ ;
      }
    }
    printf("%d unknown vertices found\n", nzero) ;
    MRISsegmentAnnotated(mris, &labels, &nlabels, 10) ;
    most_vertices = 0 ;
    biggest_label = -1 ;
    for (i = 0 ; i < nlabels ; i++) {
      label = labels[i] ;
      if (mris->vertices[label->lv[0].vno].annotation == TMP_LABEL) {
        if (label->n_points > most_vertices) {
          biggest_label = i ;
          most_vertices = label->n_points ;
        }
      }
    }
    if (biggest_label >= 0) {
      label = labels[biggest_label] ;
      printf("replacing label # %d with %d vertices "
             "(vno=%d) with label %d\n",
             biggest_label,
             label->n_points,
             label->lv[0].vno,
             unknown_label) ;
      for (i = 0 ; i < label->n_points ; i++) {
        v = &mris->vertices[label->lv[i].vno] ;
        v->annotation = v->val = unknown_label ;
      }
    }
    for (nzero = vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->annotation == TMP_LABEL) {
        v->annotation = 0;
        nzero++ ;
      }
    }
    printf("after replacement, %d unknown vertices found\n", nzero) ;
    MRISmodeFilterZeroVals(mris) ;  /* get rid of the rest
                                    of the unknowns by mode filtering */
    for (i = 0 ; i < nlabels ; i++)
      LabelFree(&labels[i]) ;
    free(labels) ;
  }

  MRIScopyValsToAnnotations(mris) ;
  if (fix_topology != 0)
    fix_label_topology(mris, fix_topology) ;

  if (mode_filter) {
    printf("mode filtering sample labels...\n") ;
#if 0
    MRISmodeFilterZeroVals(mris) ;
#else
    MRISmodeFilterVals(mris, mode_filter) ;
#endif
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->annotation = v->val ;
    }
  }

  /* this will fill in the v->annotation field from the v->val ones */
  translate_indices_to_annotations(mris, translation_fname) ;

  if (label_index >= 0)
  {
    int index ;
    LABEL *area ;

    printf("writing label to %s...\n", annot_name) ;
    MRISclearMarks(mris) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      if (vno == Gdiag_no)
        DiagBreak() ;
      v = &mris->vertices[vno] ;
      if (v->annotation > 0)
        DiagBreak() ;
      CTABfindAnnotation(mris->ct, v->annotation, &index);
      if (index == label_index)
        v->marked = 1 ;
    }
    area = LabelFromMarkedSurface(mris) ;
    if (nclose > 0)
    {
      LabelDilate(area, mris, nclose, CURRENT_VERTICES) ;
      LabelErode(area, mris, nclose) ;
    }
    LabelWrite(area, annot_name) ;
  }
  else
  {
    printf("writing annotation to %s...\n", annot_name) ;
    MRISwriteAnnotation(mris, annot_name) ;
  }
  /*  MRISreadAnnotation(mris, fname) ;*/
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
  else if (!stricmp(option, "sdir")) {
    strcpy(sdir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as SUBJECTS_DIR\n", sdir) ;
  } else if (!stricmp(option, "label")) {
    label_index = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using label %s (%d) and writing output in label format\n",
           cma_label_to_name(label_index), label_index) ;
  } else if (!stricmp(option, "surf")) {
    surf_name = argv[2] ;
    nargs = 1 ;
    printf("using %s as surface name\n", surf_name) ;
  }
  else if (!stricmp(option, "cortex"))
  {
    cortex_label = LabelRead(NULL, argv[2]) ;
    if (cortex_label == NULL)
    {
      ErrorExit(ERROR_NOFILE, "") ;
    }
    nargs = 1 ;
    printf("INFO: using %s as mask to calc cortex "
           "NumVert, SurfArea and MeanThickness.\n",
           argv[2]);
  } else if (!stricmp(option, "mask")) {
    mask_fname = argv[2] ;
    mask_val = atoi(argv[3]) ;
    nargs = 2 ;
    printf("using %d as a mask in %s\n", mask_val, mask_fname) ;
  } else if (!stricmp(option, "vol2surf")) {
    sample_from_vol_to_surf = 1 ;
    printf("sampling from volume to surface...\n") ;
  } else if (!stricmp(option, "close")) {
    nclose = atoi(argv[2]) ;
    printf("applying %dth order morphological close to label\n", nclose) ;
    nargs = 1 ;
  } else if (!stricmp(option, "fix")) {
    fix_topology = atoi(argv[2]) ;
    printf("fixing topology of all labels smaller "
           "than %d vertices\n",fix_topology);
    nargs = 1 ;
  } else if (!stricmp(option, "replace")) {
    replace_label = atoi(argv[2]) ;
    nargs = 1 ;
    printf("replacing label %s (%d) with deeper ones\n",
           cma_label_to_name(replace_label), replace_label) ;
  } else if (!stricmp(option, "trans")) {
    if (ntrans >= MAX_TRANS)
      ErrorExit(ERROR_NOMEMORY, "%s: too many translations (%d)\n",
                Progname,ntrans);
    trans_in[ntrans] = atof(argv[2]) ;
    trans_out[ntrans] = atof(argv[3]) ;
    nargs = 2 ;
    printf("translating %s (%d) to %s (%d)\n",
           cma_label_to_name(trans_in[ntrans]), trans_in[ntrans],
           cma_label_to_name(trans_out[ntrans]), trans_out[ntrans]) ;
    ntrans++ ;
  } else if (!stricmp(option, "projmm") || !stricmp(option, "proj")) {
    proj_mm = atof(argv[2]) ;
    nargs = 1 ;
    printf("projecting %2.2f mm along surface normal\n", proj_mm) ;
  } else if (!stricmp(option, "projfrac")) {
    proj_frac = atof(argv[2]) ;
    nargs = 1 ;
    printf("projecting %2.2f %% along surface normal\n", proj_frac) ;
  } else if (!stricmp(option, "file")) {
    translation_fname = argv[2] ;
    nargs = 1 ;
    printf("using %s as translation fname\n", sdir) ;
  } else if (!stricmp(option, "ct")) {
    color_table_fname = argv[2] ;
    nargs = 1 ;
    printf("embedding color table %s into output annot file\n",
           color_table_fname) ;
    translation_fname = color_table_fname ;
  } else switch (toupper(*option)) {
    case 'V':
      if (argc==2) {
        print_help();
        exit(1);
      }
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'F':
      mode_filter = atoi(argv[2]) ;
      nargs = 1 ;
      printf("applying mode filter %d times to parcellation\n", mode_filter) ;
      break ;
    case 'A':
      avgs = atoi(argv[2]) ;
      nargs = 1 ;
      printf("smoothing surface %d times\n", avgs) ;
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
      if (argc==2) {
        print_help();
        exit(1);
      }
      unknown_label = atoi(argv[2]) ;
      printf("changing largest connected unknown region to label %d\n",
             unknown_label) ;
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
usage_exit(void) {
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "Usage:\n"
          "------\n"
          "\n%s [options] <subject name> <hemi> <parc name> "
          "<output annot>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program samples a volumetric parcellation "
          "onto a surface. \nManual labeling can be carried out ");
  fprintf(stderr,
          "directly on surface models \nusing drawing tools in tksurfer, "
          "or volumetrically in tkmedit, \nthen sampled onto the ");
  fprintf(stderr,
          "surfaces using mris_sample_parc.\n") ;
  fprintf(stderr,
          "mris_ca_train is used to create an atlas from a set of \n"
          "annotated subjects. The information output by mris_ca_train\n");
  fprintf(stderr,
          "is then used by mris_ca_label to automatically assign a\n"
          "neuroanatomical label to each location on a cortical "
          "surface model.\n\n");
  fprintf(stderr,
          "Required args:\n"
          "--------------\n\n") ;
  fprintf(stderr,
          "  <subject name>       the subject id\n\n");
  fprintf(stderr,
          "  <hemi>               hemisphere: rh or lh\n\n");
  fprintf(stderr,
          "  <parc name>          parcellation filename\n\n");
  fprintf(stderr,
          "  <output annot>       annotation filename\n\n");
  fprintf(stderr,
          "Optional args:\n"
          "--------------\n\n");
  fprintf(stderr,
          "  -sdir <directory>    use <directory> as subjects directory \n"
          "                       (default: $SUBJECTS_DIR)\n\n");
  fprintf(stderr,
          "  -surf <filename>     use <filename> as surface "
          "(default: 'white')\n\n");
  fprintf(stderr,
          "  -fix <number>        fix topology of all labels smaller \n"
          "                       than <number> vertices (default=-1, "
          "do all)\n\n");
  fprintf(stderr,
          "  -replace <number>    replace label <number> with deeper "
          "ones\n\n");
  fprintf(stderr,
          "  -trans <number_in> <number_out>      translate <number_in> to \n"
          "                                       <number_out>\n\n");
  fprintf(stderr,
          "  -cortex <cortex label file>          mask regions outside of the \n"
          "                                       specified cortex label\n\n");
  fprintf(stderr,
          "  -projmm <number>     project <number> millimeters along \n"
          "                       surface normal (default=0.0)\n\n");
  fprintf(stderr,
          "  -proj <number>       same as -projmm\n\n");
  fprintf(stderr,
          "  -projfrac <number>   project <number> percent along surface \n"
          "                       normal (default=0.5)\n\n");
  fprintf(stderr,
          "  -file <filename>     use <filename> as translation \n"
          "                       (default: 'cma_parcellation_colors.txt')"
          "\n\n");
  fprintf(stderr,
          "  -ct <filename>       embed color table <filename> into output \n"
          "                       annotation file\n\n");
  fprintf(stderr,
          "  -v <number>          diagnostic level (default=0)\n\n");
  fprintf(stderr,
          "  -f <number>          apply mode filter <number> times to \n"
          "                       parcellation (default=0)\n\n");
  fprintf(stderr,
          "  -a <number>          smooth surface <number> times "
          "(default=0)\n\n");
  fprintf(stderr,
          "  -w <number>          use window size <number> for sampling\n"
          "                       (default=7)\n\n");
  fprintf(stderr,
          "  -t <filename>        use thickness file <filename> \n"
          "                       (default: 'thickness')\n\n");
  fprintf(stderr,
          "  -u <number>          change largest connected unknown region to\n"
          "                       label <number> (default: don't change)\n\n");
  fprintf(stderr,
          "  --help               print help info\n\n");
  fprintf(stderr,
          "  --version            print version info\n");
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

#if 0
static int
translate_indices_to_annotations(MRI_SURFACE *mris, const char *translation_fname) {
  char   fname[STRLEN], *cp, line[STRLEN], name[STRLEN], **names ;
  FILE   *fp ;
  int    nlines, i, vno, *r, *g, *b, *indices ;
  VERTEX *v ;

  cp = getenv("FREESURFER_HOME") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM,
              "%s: FREESURFER_HOME not defined in environment", Progname) ;
  sprintf(fname, "%s/%s", cp, translation_fname) ;
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE,
              "%s: could not read translation file %s", Progname, fname) ;

  nlines = 0 ;
  cp = fgetl(line, STRLEN-1, fp) ;
  while (cp != NULL) {
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

  for (i = 0 ; i < nlines ; i++) {
    cp = fgetl(line, STRLEN-1, fp) ;
    sscanf(cp, "%d %s %d %d %d %*d",
           &indices[i], name, &r[i], &g[i], &b[i]) ;
    names[i] = (char *)calloc(strlen(name)+1, sizeof(char)) ;
    strcpy(names[i], name) ;
    printf("parsing parcellation unit %s...\n", name) ;
  }
  fclose(fp) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
  }
  return(NO_ERROR) ;
}
#else
static int
translate_indices_to_annotations(MRI_SURFACE *mris, const char *translation_fname) {
  int    vno ;
  VERTEX *v ;

  read_named_annotation_table(translation_fname) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->annotation = index_to_annotation((int)v->val) ;
    if (v->annotation == -1)
      v->annotation = (unknown_label > 0) ? unknown_label : 0 ;
  }
  return(NO_ERROR) ;
}
#endif
static int
fix_label_topology(MRI_SURFACE *mris, int nvertices) {
  int    i, vno, nsegments, most_vertices, max_label;
  int    label, j, iter, nchanged=0, max_index, gdiag_seg ;
  LABEL **segments ;
  VERTEX *v ;

  if (nvertices < 0)
    nvertices = mris->nvertices+1 ;
  for (max_label = vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->annotation > max_label)
      max_label = v->annotation ;
  }

  iter = 0 ;
  do {
    nchanged = 0 ;
    MRISsegmentAnnotated(mris, &segments, &nsegments, 0) ;
    gdiag_seg = -1 ;
    if (Gdiag_no >= 0)  // find segment with Gdiag_no in it
    {
      for (i = 0 ; gdiag_seg < 0 && i < nsegments ; i++) 
      {
	for (j = 0 ; j < segments[i]->n_points ; j++)
	  if (segments[i]->lv[j].vno == Gdiag_no)
	  {
	    gdiag_seg = i ;
	    break ;
	  }
      }
    }

    
    for (i = 0 ; i < nsegments ; i++) {
      if (segments[i] == NULL)
        continue ;
      most_vertices = segments[i]->n_points ;
      max_index = i ;
      label = mris->vertices[segments[i]->lv[0].vno].annotation ;

      /* find label with most vertices */
      for (j = 0 ; j < nsegments ; j++) {
        if (j == i || segments[j] == NULL)
          continue ;
        if (mris->vertices[segments[j]->lv[0].vno].annotation != label)
          continue ;
        if (segments[j]->n_points > most_vertices) {
          most_vertices = segments[j]->n_points ;
          max_index = j ;
        }
      }

      /* resegment all others connected-components with this label */
      for (j = 0 ; j < nsegments ; j++) {
        if (j == max_index || segments[j] == NULL ||
            (segments[j]->n_points > nvertices))
          continue ;
        if (mris->vertices[segments[j]->lv[0].vno].annotation != label)
          continue ;
	if (j == gdiag_seg)
	  DiagBreak() ;
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
  } while (nchanged > 0 && iter++ < 10) ;  // can oscillate forever, so terminate after 10

  MRISclearMarks(mris) ;
  return(NO_ERROR) ;
}

static int
resegment_label(MRI_SURFACE *mris, LABEL *segment) {
  int    *histo, max_label, i, n, vno, ino, index;
  int    max_histo, max_index, nchanged, lno, label ;

  max_label = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX const * const v = &mris->vertices[vno] ;
    if (v->val > max_label)
      max_label = v->val ;
  }
  histo = (int *)calloc(max_label+1, sizeof(*histo)) ;
  if (histo == NULL)
    ErrorExit(ERROR_NOMEMORY, "resegment_label: could not allocate %d element histogram", max_label+1) ;

  label = mris->vertices[segment->lv[0].vno].annotation ;
  for (ino  = 0 ; ino < 100 ; ino++) {
    nchanged = 0 ;
    for (lno = 0 ; lno < segment->n_points ; lno++) {
      vno = segment->lv[lno].vno ;
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
      VERTEX                * const v  = &mris->vertices         [vno];
      if (vno == Gdiag_no)
        DiagBreak() ;
      if (v->val != label || v->ripflag)
        continue ;   /* already changed */

      memset(histo, 0, (max_label+1)*sizeof(*histo)) ;
      for (n = 0 ; n < vt->vnum ; n++) {
        VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
        index = (int)nint(vn->val) ;
        if (index < 0)
          continue ;
        if (vn->val != label)  /* don't care about same label */
          histo[index]++ ;
      }
      max_histo = histo[0] ;
      max_index = -1 ;
      for (i = 0 ; i <= max_label ; i++) {
        if (histo[i] > max_histo) {
          max_histo = histo[i] ;
          max_index = i ;
        }
      }
      if (max_index >= 0)
        v->valbak = max_index ;
    }
    for (lno = 0 ; lno < segment->n_points ; lno++) {
      vno = segment->lv[lno].vno ;
      VERTEX * const v = &mris->vertices[vno];
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
  free(histo) ;
  return(NO_ERROR) ;
}
static int
replace_vertices_with_label(MRI_SURFACE *mris,
                            MRI *mri,
                            int label,
                            double proj_mm) {
  int      vno, new_label, num=0 ;
  VERTEX   *v ;
  double   d, x, y, z, xw, yw, zw ;
  MRI      *mri_tmp ;

  mri_tmp = MRIreplaceValues(mri, NULL, label, 0) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag || v->annotation != label)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

#define MAX_SEARCH_LEN 6 // mm

    if (proj_mm > 0) {
      for (d = proj_mm ; d <= MAX_SEARCH_LEN ; d += proj_mm) {
        x = v->x+d*v->nx ;
        y = v->y+d*v->ny ;
        z = v->z+d*v->nz ;
        MRIsurfaceRASToVoxel(mri, x, y, z, &xw, &yw, &zw) ;
//        new_label = (int)MRIfindNearestNonzero(mri, wsize, xw, yw, zw, ((float)wsize-1)/2) ;
        new_label = (int)MRIfindNearestNonzero(mri_tmp, wsize, xw, yw, zw, -1) ;
        if (new_label != label) {
          v->annotation = v->val = new_label ;
          num++ ;
          break ;
        }
      }
    } else {
      for (d = proj_mm ; d >= -MAX_SEARCH_LEN ; d += proj_mm) {
        x = v->x+d*v->nx ;
        y = v->y+d*v->ny ;
        z = v->z+d*v->nz ;
        MRIsurfaceRASToVoxel(mri, x, y, z, &xw, &yw, &zw) ;
        new_label = (int)MRIfindNearestNonzero(mri_tmp, wsize, xw, yw, zw, ((float)wsize-1)/2) ;
        if (new_label != label && new_label > 0) {
          v->annotation = v->val = new_label ;
          num++ ;
          break ;
        }
      }
    }
    if (v->val == label)   // couldn't find a new label for it - replace it with 0 to mark it for later reprocessing
    {
      v->val = v->annotation = 0 ;
      if (vno == Gdiag_no)
	printf("replacing vertex %d label %d with 0 for later processing\n", vno, Gdiag_no) ;
    }
  }

  MRIfree(&mri_tmp) ;
  printf("%d vertex labels replaced\n", num) ;
  return(NO_ERROR) ;
}

int
MRIsampleParcellationToSurface(MRI_SURFACE *mris, MRI *mri_parc) {
  int             min_label, max_label, **label_histo, l, vno, nlabels, x, y, z, max_l ;
  float           fmin, fmax, max_count, d ;
  MRIS_HASH_TABLE *mht ;
  VERTEX          *v ;
  double          xs, ys, zs, xv, yv, zv, val ;
  MRI             *mri_parc_unused ;

  mri_parc_unused = MRIcopy(mri_parc, NULL) ;
  MRIvalRange(mri_parc, &fmin, &fmax) ;
  min_label = (int)floor(fmin) ;
  max_label = (int)ceil(fmax) ;
  nlabels = max_label - min_label + 1 ;

  label_histo = (int **)calloc(mris->nvertices, sizeof(int *)) ;
  if (label_histo == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not create label frequency histo", Progname) ;
  for (vno = 0 ; vno < mris->nvertices ;  vno++) {
    label_histo[vno] = (int *)calloc(nlabels, sizeof(int)) ;
    if (label_histo[vno] == NULL)
      ErrorExit(ERROR_NOMEMORY, "%s: could not create label frequency histo[%d] with %d bins",
                Progname, vno, nlabels) ;
  }

  mht = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 8.0) ;

  MRISclearMarks(mris) ;

  // build histograms at each vertex
  for (x = 0 ; x < mri_parc->width ; x++) {
    for (y = 0 ; y < mri_parc->height ; y++) {
      for (z = 0 ; z < mri_parc->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        l = (int)MRIgetVoxVal(mri_parc, x, y, z, 0) ;
        if (l == 0)
          continue ;
        MRIvoxelToSurfaceRAS(mri_parc, x, y, z, &xs, &ys, &zs) ;
        v = MHTfindClosestVertexInTable(mht, mris, xs, ys, zs, 0) ;
        if (v == NULL)
          continue ;
        if (sqrt(SQR(v->x-xs) + SQR(v->y-ys) + SQR(v->z-zs)) > 3)
          continue ;
        MRIsetVoxVal(mri_parc_unused, x, y, z, 0, 0) ;
        vno = v-mris->vertices ;
        if (vno == Gdiag_no)
        {
          printf("v %d: sampling from (%d, %d, %d) - %d\n",
                 vno, x, y, z, l);
          DiagBreak() ;
        }
        label_histo[vno][l-min_label]++ ;
      }
    }
  }

  MRIwrite(mri_parc_unused, "pu.mgz") ;
  for (vno = 0 ; vno < mris->nvertices ;  vno++) {
    if (vno == Gdiag_no)
      DiagBreak() ;
    max_l = 0 ;
    max_count = 0 ;
    for (l = 0 ; l < nlabels ; l++) {
      if (label_histo[vno][l] > max_count) {
        max_count = label_histo[vno][l] ;
        max_l = l+min_label ;
      }
    }
    v = &mris->vertices[vno] ;
    v->val = v->annotation = max_l ;
    if (max_count > 0)
      v->marked = 1 ;
  }
  for (vno = 0 ; vno < mris->nvertices ;  vno++) {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->marked)
      continue ;  // found something here

    for (d = 0 ; d <= 2 ; d += 0.25) {
      xs = v->x + d*v->nx;
      ys = v->y + d*v->ny;
      zs = v->z + d*v->nz;
      MRIsurfaceRASToVoxel(mri_parc, xs, ys, zs, &xv, &yv, &zv);
      MRIsampleVolumeType(mri_parc, xv, yv, zv, &val, SAMPLE_NEAREST) ;
      l = (int)nint(val) ;
      if (l > 0) {
        v->val = v->annotation = l ;
        break ;
      }
    }
  }

  MHTfree(&mht) ;
  for (vno = 0 ; vno < mris->nvertices ;  vno++)
    free(label_histo[vno]) ;
  free(label_histo) ;
  MRIfree(&mri_parc_unused) ;
  return(NO_ERROR) ;
}

