
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/utsname.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "fio.h"
#include "version.h"
#include "colortab.h"

static char vcid[] =
"$Id: mris_anatomical_stats.c,v 1.40 2006/11/03 18:29:31 fischl Exp $";

int main(int argc, char *argv[]) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
double MRISmeasureTotalWhiteMatterVolume(MRI *mri) ;
double MRISmeasureCorticalGrayMatterVolume(MRI_SURFACE *mris) ;
int MRIScomputeCurvatureStats(MRI_SURFACE *mris, double *pavg, double *pvar,
                              float ignore_below, float ignore_above) ;
double MRIScomputeAbsoluteCurvature(MRI_SURFACE *mris) ;
int    MRISrestoreSurface(MRI_SURFACE *mris) ;
int    MRIScountVertices(MRI_SURFACE *mris);
#if 0
int    MRISreadAnnotFile(MRI_SURFACE *mris, char *fname) ;
int    MRISripVerticesWithMark(MRI_SURFACE *mris, int mark) ;
int    MRISripVerticesWithoutMark(MRI_SURFACE *mris, int mark) ;
int    MRISreplaceMarks(MRI_SURFACE *mris, int in_mark, int out_mark) ;
#endif
int    MRISripVerticesWithAnnotation(MRI_SURFACE *mris, int annotation) ;
int    MRISripVerticesWithoutAnnotation(MRI_SURFACE *mris, int annotation) ;
int    MRISreplaceAnnotations(MRI_SURFACE *mris,
                              int in_annotation,
                              int out_annotation) ;
char *Progname ;
static double sigma = 0.0f ;
static float ignore_below = 0 ;
static float ignore_above = 20 ;
static char *label_name = NULL ;
static char *annotation_name = NULL ;
static char *thickness_name = "thickness" ;
static int histo_flag = 0 ;
static char *gray_histo_name ;
static char *mri_name = "T1" ;
static int noheader = 0 ;
static char *log_file_name = NULL ;
static int tabular_output_flag = 0;
static char sdir[STRLEN] = "" ;
static int MGZ = 1; // for use with MGZ format
static char *tablefile=NULL;
static char *annotctabfile=NULL; // for outputing the color table
static FILE *fp=NULL;
static int nsmooth = 0;
static char *white_name = "white" ;
static char *pial_name = "pial" ;

int
main(int argc, char *argv[])
{
  char          **av, *hemi, *sname, *cp, fname[STRLEN], *surf_name ;
  int           ac, nargs, vno ;
  MRI_SURFACE   *mris ;
  MRI           *mri_wm, *mri_kernel = NULL, *mri_orig ;
  double        gray_volume, wm_volume, thickness_mean, thickness_var,
    mean_abs_mean_curvature, mean_abs_gaussian_curvature, ici, fi ;
  int           annotation = 0 ;
  FILE          *log_fp = NULL ;
  VERTEX        *v ;
  HISTOGRAM     *histo_gray ;
  int           ct_index;
  int           n_vertices = -1;
  MRI           *ThicknessMap = NULL;
  struct utsname uts;
  char *cmdline;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option
    (argc, argv,
     "$Id: mris_anatomical_stats.c,v 1.40 2006/11/03 18:29:31 fischl Exp $",
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  cmdline = argv2cmdline(argc,argv);
  uname(&uts);

  mean_abs_mean_curvature = mean_abs_gaussian_curvature = gray_volume = 0.0 ;
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

  sname = argv[1] ;
  if (strlen(sdir) == 0)
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM,
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  hemi = argv[2] ;
  if (argc > 3)
    surf_name = argv[3] ;
  else
    surf_name = WHITE_MATTER_NAME ;

  if (sigma > 0.0)
    mri_kernel = MRIgaussian1d(sigma, 100) ;
  sprintf(fname, "%s/%s/mri/wm", sdir, sname) ;
  if(MGZ) strcat(fname, ".mgz");
  fprintf(stderr, "reading volume %s...\n", fname) ;
  mri_wm = MRIread(fname) ;
  if (!mri_wm)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

  if (mri_kernel)
  {
    fprintf(stderr, "smoothing brain volume with sigma = %2.3f\n", sigma) ;
    MRIconvolveGaussian(mri_wm, mri_wm, mri_kernel) ;
#if 0
    fprintf(stderr, "smoothing wm volume with sigma = %2.3f\n", sigma) ;
    MRIconvolveGaussian(mri_wm, mri_wm, mri_kernel) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri_wm, "/tmp/wm_smooth.mnc") ;
#endif
    MRIfree(&mri_kernel);
  }

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, surf_name) ;
  fprintf(stderr, "reading input surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  // read in white and pial surfaces
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, pial_name) ;
  fprintf(stderr, "reading input pial surface %s...\n", fname) ;
  if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, white_name) ;
  fprintf(stderr, "reading input white surface %s...\n", fname) ;
  if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;

  MRIScomputeMetricProperties(mris) ;
  wm_volume = MRISmeasureTotalWhiteMatterVolume(mri_wm) ;
#if 0
  fprintf(stderr, "measuring gray matter thickness of surface...\n") ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, GRAY_MATTER_NAME) ;
  fprintf(stderr, "reading pial surface position from %s...\n", fname) ;
  MRISreadVertexPositions(mris, fname) ;
#else
  MRISreadCurvatureFile(mris, thickness_name) ;

  if(nsmooth > 0){
    printf("Smooth thickness map with %d iterations on surface\n", nsmooth);
    ThicknessMap = MRIcopyMRIS(NULL, mris, 0, "curv");
    if(ThicknessMap == NULL){
      printf("Unable to copy thickness data to a MRI volume \n");
    }else{
      MRISsmoothMRI(mris, ThicknessMap, nsmooth, NULL,ThicknessMap);
      MRIScopyMRI(mris, ThicknessMap, 0, "curv");
      MRIfree(&ThicknessMap);
    }
  }

#endif

#if 0
  fprintf(stderr, "measuring cortical thickness...") ;
  MRISmeasureCorticalThickness(mris, mri_wm) ;
#endif
  MRIScopyCurvatureToImagValues(mris) ;   /* save thickness measures */

  fprintf(stderr, "done.\ncomputing second fundamental form...") ;
  MRISsetNeighborhoodSize(mris, 2) ;
  MRIScomputeSecondFundamentalForm(mris) ;

  if(annotation_name){
    if(MRISreadAnnotation(mris, annotation_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s:  could  not read annotation file %s",
                Progname, annotation_name) ;
    if(annotctabfile != NULL && mris->ct != NULL){
      printf("Saving annotation colortable %s\n",annotctabfile);
      CTABwriteFileASCII(mris->ct,annotctabfile);
    }
  }
  printf(" ... done.\n") ;

  if (label_name)
  {
    LABEL  *area ;
    char   fname[STRLEN] ;

    sprintf(fname, "%s/%s/label/%s", sdir, sname, label_name) ;

    area = LabelRead(NULL, fname) ;
    if (!area)
      ErrorExit(ERROR_NOFILE, "%s: could not read label file %s\n", sname) ;
    LabelRipRestOfSurface(area, mris) ;
    LabelFree(&area) ;
    MRIScomputeMetricProperties(mris) ;
  }

  if (histo_flag)
  {
    sprintf(fname, "%s/%s/mri/%s", sdir, sname, mri_name) ;
    if(MGZ) strcat(fname, ".mgz");
    fprintf(stderr, "reading volume %s...\n", fname) ;
    mri_orig = MRIread(fname) ;
    if (!mri_orig)
      ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
                Progname, fname) ;
    histo_gray = HISTOalloc(256) ;
  }
  else
  {
    histo_gray = NULL ; mri_orig = NULL ;
  }

  if (log_file_name)
  {
    log_fp = fopen(log_file_name, "a") ;
    if (!log_fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open log file %s", Progname,
                log_file_name) ;
  }

#define SHOW_WHITE_MATTER_VOLUME 0
#if SHOW_WHITE_MATTER_VOLUME
  fprintf(stdout, "total white matter volume               = %2.0f mm^3\n",wm_volume) ;
#endif

  if (annotation_name && tabular_output_flag)
  {
    fprintf(stdout, "\n");
    fprintf(stdout, "table columns are:\n");
    fprintf(stdout, "    number of vertices\n");
    fprintf(stdout, "    total surface area (mm^2)\n");
    fprintf(stdout, "    total gray matter volume (mm^3)\n");
    fprintf(stdout, "    average cortical thickness "
            "+- standard deviation (mm)\n");
    fprintf(stdout, "    integrated rectified mean curvature\n");
    fprintf(stdout, "    integrated rectified Gaussian curvature\n");
    fprintf(stdout, "    folding index\n");
    fprintf(stdout, "    intrinsic curvature index\n");
    fprintf(stdout, "    structure name\n");
    fprintf(stdout, "\n");
  }

  if(annotation_name && tablefile != NULL){
    fp = fopen(tablefile,"w");
    fprintf(fp,"# Table of FreeSurfer cortical "
            "parcellation anatomical statistics \n");
    fprintf(fp,"# \n");
    fprintf(fp,"# CreationTime %s\n",VERcurTimeStamp());
    fprintf(fp,"# generating_program %s\n",Progname);
    fprintf(fp,"# cvs_version %s\n",vcid);
    fprintf(fp,"# mrisurf.c-cvs_version %s\n",MRISurfSrcVersion());
    fprintf(fp,"# cmdline %s\n",cmdline);
    fprintf(fp,"# sysname  %s\n",uts.sysname);
    fprintf(fp,"# hostname %s\n",uts.nodename);
    fprintf(fp,"# machine  %s\n",uts.machine);
    fprintf(fp,"# user     %s\n",VERuser());
    fprintf(fp,"# \n");
    fprintf(fp,"# SUBJECTS_DIR %s\n",sdir);
    fprintf(fp,"# anatomy_type surface\n");
    fprintf(fp,"# subjectname %s\n",sname);
    fprintf(fp,"# hemi %s\n",hemi);
    fprintf(fp,"# AnnotationFile %s\n",annotation_name);
    fprintf(fp,"# AnnotationFileTimeStamp %s\n",
            VERfileTimeStamp(annotation_name));
#if SHOW_WHITE_MATTER_VOLUME
    fprintf(fp,"# TotalWhiteMatterVolume  %2.0f mm^3\n",wm_volume) ;
#endif

    fprintf(fp,"# Measure Cortex, NumVert, Number of Vertices, %d, unitless\n",
            mris->nvertices);
    fprintf(fp,"# Measure Cortex, SurfArea, Surface Area,  %g, mm^2\n",
            mris->total_area);

    fprintf(fp,"# NTableCols 10\n");

    fprintf(fp,"# TableCol  1 ColHeader StructName\n");
    fprintf(fp,"# TableCol  1 FieldName Structure Name\n");
    fprintf(fp,"# TableCol  1 Units     NA\n");

    fprintf(fp,"# TableCol  2 ColHeader NumVert\n");
    fprintf(fp,"# TableCol  2 FieldName Number of Vertices\n");
    fprintf(fp,"# TableCol  2 Units     unitless\n");

    fprintf(fp,"# TableCol  3 ColHeader SurfArea\n");
    fprintf(fp,"# TableCol  3 FieldName Surface Area\n");
    fprintf(fp,"# TableCol  3 Units     mm^2\n");

    fprintf(fp,"# TableCol  4 ColHeader GrayVol\n");
    fprintf(fp,"# TableCol  4 FieldName Gray Matter Volume\n");
    fprintf(fp,"# TableCol  4 Units     mm\n");

    fprintf(fp,"# TableCol  5 ColHeader ThickAvg \n");
    fprintf(fp,"# TableCol  5 FieldName Average Thickness\n");
    fprintf(fp,"# TableCol  5 Units     mm\n");

    fprintf(fp,"# TableCol  6 ColHeader ThickStd\n");
    fprintf(fp,"# TableCol  6 FieldName Thickness StdDev\n");
    fprintf(fp,"# TableCol  6 Units     mm \n");

    fprintf(fp,"# TableCol  7 ColHeader MeanCurv\n");
    fprintf(fp,"# TableCol  7 FieldName Integrated Rectified"
            " Mean Curvature\n");
    fprintf(fp,"# TableCol  7 Units     mm^-1\n");

    fprintf(fp,"# TableCol  8 ColHeader GausCurv \n");
    fprintf(fp,"# TableCol  8 FieldName Integrated Rectified"
            " Gaussian Curvature\n");
    fprintf(fp,"# TableCol  8 Units     mm^-2\n");

    fprintf(fp,"# TableCol  9 ColHeader  FoldInd\n");
    fprintf(fp,"# TableCol  9 FieldName  Folding Index \n");
    fprintf(fp,"# TableCol  9 Units      unitless \n");

    fprintf(fp,"# TableCol 10 ColHeader CurvInd\n");
    fprintf(fp,"# TableCol 10 FieldName Intrinsic Curvature Index\n");
    fprintf(fp,"# TableCol 10 Units     unitless\n");

    fprintf(fp,"# ColHeaders StructName NumVert SurfArea GrayVol "
            "ThickAvg ThickStd MeanCurv GausCurv FoldInd CurvInd\n");
    fclose(fp);
  }

  if (annotation_name) // new stuff
  {
#define MAX_INDICES 1000
    double  areas[MAX_INDICES],
      volumes[MAX_INDICES], thicknesses[MAX_INDICES],
      avg_thick, volume, thickness_vars[MAX_INDICES], std ;
    int     v0_index, v1_index, v2_index, fno, m, i, dofs[MAX_INDICES] ;
    VERTEX  *v0, *v1, *v2 ;
    FACE    *f ;

    memset(areas, 0, sizeof(areas)) ;
    memset(volumes, 0, sizeof(volumes)) ;
    memset(thicknesses, 0, sizeof(thicknesses)) ;
    memset(dofs, 0, sizeof(dofs)) ;
    memset(thickness_vars, 0, sizeof(thickness_vars)) ;

    MRIScomputeMetricProperties(mris) ;

    // first do white surface
    MRISsaveVertexPositions(mris, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      f = &mris->faces[fno] ;
      if (f->ripflag)
        continue ;
      v0 = &mris->vertices[f->v[0]] ;
      CTABfindAnnotation(mris->ct, v0->annotation,&v0_index);
      v1 = &mris->vertices[f->v[1]] ;
      CTABfindAnnotation(mris->ct, v1->annotation,&v1_index);
      v2 = &mris->vertices[f->v[2]] ;
      CTABfindAnnotation(mris->ct, v2->annotation,&v2_index);
      for (avg_thick = 0.0, m = 0 ; m < VERTICES_PER_FACE ; m++)
      {
        vno = f->v[m] ;
        v = &mris->vertices[vno] ;
        avg_thick += v->imag_val ;
      }
      avg_thick /= VERTICES_PER_FACE ;
      volume = (avg_thick * f->area) ;
      volumes[v0_index] += volume/VERTICES_PER_FACE ;
      volumes[v1_index] += volume/VERTICES_PER_FACE ;
      volumes[v2_index] += volume/VERTICES_PER_FACE ;
    }

    // now do pial surface
    MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      f = &mris->faces[fno] ;
      if (f->ripflag)
        continue ;
      v0 = &mris->vertices[f->v[0]] ;
      CTABfindAnnotation(mris->ct, v0->annotation,&v0_index);
      v1 = &mris->vertices[f->v[1]] ;
      CTABfindAnnotation(mris->ct, v1->annotation,&v1_index);
      v2 = &mris->vertices[f->v[2]] ;
      CTABfindAnnotation(mris->ct, v2->annotation,&v2_index);
      for (avg_thick = 0.0, m = 0 ; m < VERTICES_PER_FACE ; m++)
      {
        vno = f->v[m] ;
        v = &mris->vertices[vno] ;
        avg_thick += v->imag_val ;
      }
      avg_thick /= VERTICES_PER_FACE ;
      volume = (avg_thick * f->area) ;
      volumes[v0_index] += volume/VERTICES_PER_FACE ;
      volumes[v1_index] += volume/VERTICES_PER_FACE ;
      volumes[v2_index] += volume/VERTICES_PER_FACE ;
    }

    // area should just be surface area os
    // specified surface (in ORIG_VERTICES)
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      f = &mris->faces[fno] ;
      if (f->ripflag)
        continue ;
      v0 = &mris->vertices[f->v[0]] ;
      CTABfindAnnotation(mris->ct, v0->annotation,&v0_index );
      v1 = &mris->vertices[f->v[1]] ;
      CTABfindAnnotation(mris->ct, v1->annotation,&v1_index );
      v2 = &mris->vertices[f->v[2]] ;
      CTABfindAnnotation(mris->ct, v2->annotation,&v2_index );
      areas[v0_index] += f->area/VERTICES_PER_FACE ;
      areas[v1_index] += f->area/VERTICES_PER_FACE ;
      areas[v2_index] += f->area/VERTICES_PER_FACE ;
    }

    // compute thickness for each annotation
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v0 = &mris->vertices[vno] ;
      CTABfindAnnotation(mris->ct, v0->annotation,&v0_index );
      thicknesses[v0_index] += v0->imag_val ;
      dofs[v0_index]++ ;
    }

    for (i = 0 ; i < MAX_INDICES ; i++)
    {
      if (dofs[i] == 0)
        continue ;
      thicknesses[i] /= dofs[i] ;
    }
    // compute thickness variance for each annotation
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v0 = &mris->vertices[vno] ;
      CTABfindAnnotation(mris->ct, v0->annotation,&v0_index );
      std = (v0->imag_val-thicknesses[v0_index]);
      thickness_vars[v0_index] += std*std ;
    }

    for (i = 0 ; i < MAX_INDICES ; i++)
    {
      if (dofs[i] == 0)
        continue ;
      if (CTABannotationAtIndex(mris->ct, i, &annotation) != NO_ERROR)
        continue ;
      MRISripVerticesWithoutAnnotation(mris, annotation) ;
      MRISuseMeanCurvature(mris) ;
      mean_abs_mean_curvature = MRIScomputeAbsoluteCurvature(mris) ;

      MRISuseGaussianCurvature(mris) ;
      mean_abs_gaussian_curvature = MRIScomputeAbsoluteCurvature(mris) ;
      MRIScomputeCurvatureIndices(mris, &ici, &fi) ;
      MRISrestoreSurface(mris) ;

      volumes[i] /= 2 ;
      thickness_vars[i] /= dofs[i] ;
      /* output */

      if(annotation_name && tablefile != NULL){
        fp = fopen(tablefile,"a");
        fprintf(fp, "%-40s", mris->ct->entries[i]->name);
        fprintf(fp, "%5d", dofs[i]);
        fprintf(fp, "  %5.0f", areas[i]) ;
        fprintf(fp, "  %5.0f", volumes[i]) ;
        fprintf(fp, "  %5.3f %5.3f",
                thicknesses[i], sqrt(thickness_vars[i])) ;
        fprintf(fp, "  %8.3f", mean_abs_mean_curvature) ;
        fprintf(fp, "  %8.3f", mean_abs_gaussian_curvature) ;
        fprintf(fp, "  %7.3f", fi);
        fprintf(fp, "  %6.3f",ici);
        fprintf(fp, "\n");
        fclose(fp);
      }

      if(tabular_output_flag)
      {
        fprintf(stdout, "%5d", dofs[i]);
        fprintf(stdout, "  %5.0f", areas[i]) ;
        fprintf(stdout, "  %5.0f", volumes[i]) ;
        fprintf(stdout, "  %5.3f %5.3f",
                thicknesses[i], sqrt(thickness_vars[i])) ;
        fprintf(stdout, "  %8.3f", mean_abs_mean_curvature) ;
        fprintf(stdout, "  %8.3f", mean_abs_gaussian_curvature) ;
        fprintf(stdout, "  %7.3f", fi);
        fprintf(stdout, "  %6.3f",ici);
        fprintf(stdout, "  %s", mris->ct->entries[i]->name);
        fprintf(stdout, "\n");
      }
      else
      {
        if (mris->ct == NULL)
          ErrorExit
            (ERROR_BADFILE,
             "%s: no color table loaded - cannot translate annot  file",
             Progname);
        fprintf(stdout,
                "structure is \"%s\"\n", mris->ct->entries[i]->name);
        fprintf(stdout,
                "number of vertices                      = %d\n",
                dofs[i]);
        fprintf(stdout,
                "total surface area                      = %2.0f mm^2\n",
                areas[i]) ;
        fprintf(stdout,
                "total gray matter volume                = %2.0f mm^3\n",
                volumes[i]) ;
        fprintf(stdout,
                "average cortical thickness              = %2.3f mm "
                "+- %2.3f mm\n",
                thicknesses[i], sqrt(thickness_vars[i])) ;
        fprintf(stdout,
                "average integrated rectified mean curvature     = "
                "%2.3f\n",
                mean_abs_mean_curvature) ;
        fprintf(stdout,
                "average integrated rectified Gaussian curvature = "
                "%2.3f\n",
                mean_abs_gaussian_curvature) ;
        fprintf(stdout,
                "folding index                           = %2.3f\n", fi);
        fprintf(stdout,
                "intrinsic curvature index               = %2.3f\n",ici);
      }
    }
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    exit(0) ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (!histo_flag && annotation_name == NULL)
      break ;
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (histo_flag)
    {
      double d, xw, yw, zw, x, y, z, thickness, val ;

      thickness = v->imag_val ;

      for (d = 0.5 ; d < (thickness-0.5) ; d += 0.1)
      {
        x = v->x+d*v->nx ; y = v->y+d*v->ny ; z = v->z+d*v->nz ;
        // MRIworldToVoxel(mri_orig, x, y, z, &xw, &yw, &zw) ;
        MRIsurfaceRASToVoxel(mri_orig, x, y, z, &xw, &yw, &zw) ;
        MRIsampleVolume(mri_orig, xw, yw, zw, &val) ;
        HISTOaddSample(histo_gray, val, 0, 255) ;
      }
    }
    if (annotation_name)
    {

      annotation = v->annotation ;
      if (mris->ct && Gdiag_no >= 0){
        CTABfindAnnotation(mris->ct, annotation,&ct_index );
        if (ct_index == Gdiag_no) /* 6 is ectorhinal */
          DiagBreak() ;
      }

      MRISripVerticesWithoutAnnotation(mris, annotation) ;

      n_vertices = MRIScountVertices(mris);

      MRIScomputeMetricProperties(mris) ;

      MRIScopyCurvatureFromImagValues(mris) ;  /* restore
                                                  thickness measures */
      MRIScomputeCurvatureStats(mris, &thickness_mean, &thickness_var,
                                ignore_below, ignore_above) ;

      gray_volume = MRISmeasureCorticalGrayMatterVolume(mris) ;

      MRISuseMeanCurvature(mris) ;
      mean_abs_mean_curvature = MRIScomputeAbsoluteCurvature(mris) ;

      MRISuseGaussianCurvature(mris) ;
      mean_abs_gaussian_curvature = MRIScomputeAbsoluteCurvature(mris) ;

      MRIScomputeCurvatureIndices(mris, &ici, &fi) ;

      /* output */

      if(annotation_name && tablefile != NULL){
        fp = fopen(tablefile,"a");
        CTABfindAnnotation(mris->ct, annotation,&ct_index );
        if(ct_index < 0)
          fprintf(fp, "  ** annotation %08x", annotation);
        else
          fprintf(fp, "%-40s", mris->ct->entries[ct_index]->name);
        fprintf(fp, "%5d", n_vertices);
        fprintf(fp, "  %5.0f", mris->total_area) ;
        fprintf(fp, "  %5.0f", gray_volume) ;
        fprintf(fp, "  %5.3f %5.3f", thickness_mean, sqrt(thickness_var)) ;
        fprintf(fp, "  %8.3f", mean_abs_mean_curvature) ;
        fprintf(fp, "  %8.3f", mean_abs_gaussian_curvature) ;
        fprintf(fp, "  %7.3f", fi);
        fprintf(fp, "  %6.3f",ici);
        fprintf(fp, "\n");
        fclose(fp);
      }

      if(tabular_output_flag)
      {

        fprintf(stdout, "%5d", n_vertices);
        fprintf(stdout, "  %5.0f", mris->total_area) ;
        fprintf(stdout, "  %5.0f", gray_volume) ;
        fprintf(stdout, "  %5.3f +- %5.3f",
                thickness_mean, sqrt(thickness_var)) ;
        fprintf(stdout, "  %8.3f", mean_abs_mean_curvature) ;
        fprintf(stdout, "  %8.3f", mean_abs_gaussian_curvature) ;
        fprintf(stdout, "  %7.3f", fi);
        fprintf(stdout, "  %6.3f",ici);

        CTABfindAnnotation(mris->ct, annotation,&ct_index );
        if(ct_index < 0)
          fprintf(stdout, "  ** annotation %08x", annotation);
        else
          fprintf(stdout, "  %s", mris->ct->entries[ct_index]->name);

        fprintf(stdout, "\n");

      }
      else
      {

        if (mris->ct == NULL)
          ErrorExit
            (ERROR_BADFILE,
             "%s: no color table loaded - cannot translate annot  file",
             Progname);
        CTABfindAnnotation(mris->ct, annotation,&ct_index );

        if(ct_index < 0)
          fprintf
            (stdout,
             "statistics for unknown label: annotation %d (%08x) "
             "(%d %d %d)\n",
             annotation,
             annotation,
             annotation&0xff,
             (annotation>>8)&0xff,
             (annotation>>16)&0xff) ;
        else
          fprintf(stdout,
                  "structure is \"%s\"\n",
                  mris->ct->entries[ct_index]->name);
        fprintf(stdout,
                "number of vertices                      = %d\n",
                n_vertices);
        fprintf(stdout,
                "total surface area                      = %2.0f mm^2\n",
                mris->total_area) ;
        fprintf(stdout,
                "total gray matter volume                = %2.0f mm^3\n",
                gray_volume) ;
        fprintf(stdout,
                "average cortical thickness              = %2.3f mm "
                "+- %2.3f mm\n",
                thickness_mean, sqrt(thickness_var)) ;
        fprintf(stdout,
                "average integrated rectified mean curvature     = "
                "%2.3f\n",
                mean_abs_mean_curvature) ;
        fprintf(stdout,
                "average integrated rectified Gaussian curvature = "
                "%2.3f\n",
                mean_abs_gaussian_curvature) ;
        fprintf(stdout,
                "folding index                           = %2.3f\n", fi);
        fprintf(stdout,
                "intrinsic curvature index               = %2.3f\n",ici);
      }

      MRISrestoreSurface(mris) ;
      MRISreplaceAnnotations(mris, annotation, -1) ;
      MRISripVerticesWithAnnotation(mris, -1) ;

    }
  }

  if (annotation_name == NULL)
  {
    MRIScopyCurvatureFromImagValues(mris) ;  /* restore thickness measures */
    MRIScomputeCurvatureStats(mris, &thickness_mean, &thickness_var,
                              ignore_below, ignore_above) ;

    fprintf(stdout,
            "total surface area                      = %2.0f mm^2\n",
            mris->total_area) ;

    gray_volume = MRISmeasureCorticalGrayMatterVolume(mris) ;
    fprintf(stdout,
            "total gray matter volume                = %2.0f mm^3\n",
            gray_volume) ;

    fprintf
      (stdout,
       "average cortical thickness              = %2.3f mm +- %2.3f mm\n",
       thickness_mean, sqrt(thickness_var)) ;

    MRISuseMeanCurvature(mris) ;
    mean_abs_mean_curvature = MRIScomputeAbsoluteCurvature(mris) ;
    fprintf
      (stdout,
       "average integrated rectified mean curvature     = %2.3f\n",
       mean_abs_mean_curvature) ;
    MRISuseGaussianCurvature(mris) ;
    mean_abs_gaussian_curvature = MRIScomputeAbsoluteCurvature(mris) ;
    fprintf
      (stdout,
       "average integrated rectified Gaussian curvature = %2.3f\n",
       mean_abs_gaussian_curvature) ;
    MRIScomputeCurvatureIndices(mris, &ici, &fi) ;
    fprintf(stdout, "folding index                           = %2.3f\n", fi);
    fprintf
      (stdout, "intrinsic curvature index               = %2.3f\n", ici);
  }
  if (log_fp)
  {
#if SHOW_WHITE_MATTER_VOLUME
    if (!noheader)
      fprintf(log_fp, "%% %s: <wm vol> <surf area> <gray vol> "
              "<thick mean> <thick var> <integ rect. mean curv> "
              "<integ rect. Gauss curv> <fold index> <intr curv ind>\n",
              sname) ;
    fprintf
      (log_fp,
       "%2.0f\t%2.0f\t%2.0f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\n",
       wm_volume,
       mris->total_area,
       gray_volume,
       thickness_mean,
       sqrt(thickness_var),
       mean_abs_mean_curvature,
       mean_abs_gaussian_curvature,
       fi,
       ici) ;
#else
    if (!noheader)
      fprintf(log_fp, "%% %s: <surf area> <gray vol> "
              "<thick mean> <thick var> <integ rect. mean curv> "
              "<integ rect. Gauss curv> <fold index> <intr curv ind>\n",
              sname) ;
    fprintf
      (log_fp,
       "%2.0f\t%2.0f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\t%2.3f\n",
       mris->total_area,
       gray_volume,
       thickness_mean,
       sqrt(thickness_var),
       mean_abs_mean_curvature,
       mean_abs_gaussian_curvature,
       fi,
       ici) ;
#endif
    fclose(log_fp) ;
  }
  if (histo_flag)
  {
    fprintf(stderr, "plotting gray matter histogram to file %s...\n",
            gray_histo_name) ;
    HISTOplot(histo_gray, gray_histo_name) ;
    MRIfree(&mri_orig) ;
    HISTOfree(&histo_gray) ;
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "log"))
  {
    log_file_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "outputting results to %s...\n", log_file_name) ;
  }
  else if (!stricmp(option, "nsmooth"))
  {
    nsmooth = atoi(argv[2]);
    nargs = 1;
    printf("Smooth thickness by %d steps before using it \n", nsmooth);
  }
  else if (!stricmp(option, "noheader"))
  {
    noheader = 1 ;
    printf("supressing printing of headers to log file\n") ;
  }
  else if (!stricmp(option, "white"))
  {
    white_name = argv[2] ;
    nargs = 1 ;
    printf("using %s as white matter surface name\n", white_name) ;
  }
  else if (!stricmp(option, "pial"))
  {
    pial_name = argv[2] ;
    nargs = 1 ;
    printf("using %s as pial matter surface name\n", pial_name) ;
  }
  else if (!stricmp(option, "sdir"))
  {
    char str[STRLEN] ;
    strcpy(sdir, argv[2]) ;
    printf("using  %s as  SUBJECTS_DIR...\n", sdir)  ;
    nargs = 1 ;
    sprintf(str, "SUBJECTS_DIR=%s", sdir) ;
    putenv(str) ;
  }
  else if (!stricmp(option, "mgz"))
  {
    MGZ = 1;
    printf("INFO: assuming MGZ format for volumes.\n");
  }
  else if (!stricmp(option, "COR"))
  {
    MGZ = 0;
    printf("INFO: assuming COR format for volumes.\n");
  }
  else switch (toupper(*option))
  {
  case 'T':
    thickness_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "using thickness file %s.\n", thickness_name) ;
    break ;
  case 'L':
    label_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "limiting computations to label %s.\n", label_name) ;
    break ;
  case 'M':
    mri_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "computing histograms on intensity values from %s...\n",
            mri_name) ;
    break ;
  case 'H':
    histo_flag = 1 ;
    gray_histo_name = argv[2] ;
    nargs = 1 ;
    fprintf
      (stderr,
       "writing histograms of intensity distributions to %s...\n",
       gray_histo_name);
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'A':
    annotation_name = argv[2] ;
    nargs = 1 ;
    fprintf(stderr, "computing statistics for each annotation in %s.\n",
            annotation_name) ;
    break ;
  case 'C':
    annotctabfile = argv[2] ;
    nargs = 1 ;
    break ;
  case 'I':
    ignore_below = atof(argv[2]) ;
    ignore_above = atof(argv[3]) ;
    fprintf(stderr,
            "only considering thicknesses in the range [%2.1f,%2.1f].\n",
            ignore_below, ignore_above) ;
    nargs = 2 ;
    break ;
  case 'B':
    tabular_output_flag = 1;
    nargs = 0;
    break;
  case 'F':
    tablefile = argv[2] ;
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
usage_exit(void)
{
  print_help() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr,
          "usage: %s [options] <subject name> <hemi> [<surface name>]\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf
    (stderr,
     "\nThis program measures a variety of anatomical properties\n") ;
  fprintf
    (stderr, "\nvalid options are:\n\n") ;
  fprintf
    (stderr,
     "-i  <low thresh> <hi thresh> - only consider thicknesses in\n"
     "                               the specified range.\n") ;
  fprintf
    (stderr,
     "-l <label file>              - limit calculations to specified "
     "label\n") ;
  fprintf
    (stderr,
     "-t <thickness file>          - use specified file for computing "
     "thickness statistics\n") ;
  fprintf
    (stderr,
     "-a <annotation file>         - compute properties for each label\n"
     "                               in the annotation file separately"
     "\n") ;
  fprintf
    (stderr,
     "-b                           - tabular output\n");
  fprintf
    (stderr,
     "-f tablefile  - table output to a file (different format than -b) \n");
  fprintf
    (stderr,
     "-log <log>    - will write the stats into a file named <log>\n");
  fprintf
    (stderr,
     "-nsmooth <#>   - will smooth the thicknessmap # of "
     "iterations before using it \n");
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

double
MRISmeasureTotalWhiteMatterVolume(MRI *mri)
{
  double  total_volume, voxel_volume ;
  int     x, y, z, width, height, depth ;
  BUFTYPE *psrc ;

  width = mri->width ; height = mri->height ; depth = mri->depth ;
  voxel_volume = mri->xsize * mri->ysize * mri->zsize ;
  for (total_volume = 0.0, y = 0 ; y < height ; y++)
	{
		for (z = 0 ; z < depth ; z++)
		{
			psrc = &MRIvox(mri, 0, y, z) ;
			for (x = 0 ; x < width ; x++)
			{
				if (*psrc++ > 10)
					total_volume += voxel_volume ;
			}
		}
	}
  return(total_volume) ;
}

int
MRIScomputeCurvatureStats(MRI_SURFACE *mris, double *pavg, double *pvar,
                          float ignore_below, float ignore_above)
{
  VERTEX    *v ;
  int       vno ;
  double    mean, var, n ;

  for (n = mean = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
	{
		v = &mris->vertices[vno] ;
		if (v->ripflag)
			continue ;
		if (v->curv < ignore_below || v->curv > ignore_above)
			continue ;
		mean += v->curv ;
		n += 1.0 ;
	}

  mean /= n ;
  for (var = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
	{
		v = &mris->vertices[vno] ;
		if (v->ripflag)
			continue ;
		if (v->curv < ignore_below || v->curv > ignore_above)
			continue ;
		var += (v->curv - mean) * (v->curv - mean) ;
	}

  var /= n ;
  *pavg = mean ; *pvar = var ;
  return(NO_ERROR) ;
}
double
MRISmeasureCorticalGrayMatterVolume(MRI_SURFACE *mris)
{
  FACE      *f ;
  VERTEX    *v ;
  int       fno, m, vno ;
  double    mean, total, volume, n, avg_thick ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  for (n = total = 0.0, fno = 0 ; fno < mris->nfaces ; fno++)
    {
      f = &mris->faces[fno] ;
      if (f->ripflag)
        continue ;
      for (avg_thick = 0.0, m = 0 ; m < VERTICES_PER_FACE ; m++)
        {
          vno = f->v[m] ;
          v = &mris->vertices[vno] ;
          avg_thick += v->curv ;
        }
      avg_thick /= VERTICES_PER_FACE ;
      volume = (avg_thick * f->area) ;
      total += volume ;
      n += 1.0 ;
    }

  MRISrestoreVertexPositions(mris, PIAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      f = &mris->faces[fno] ;
      if (f->ripflag)
        continue ;
      for (avg_thick = 0.0, m = 0 ; m < VERTICES_PER_FACE ; m++)
        {
          vno = f->v[m] ;
          v = &mris->vertices[vno] ;
          avg_thick += v->curv ;
        }
      avg_thick /= VERTICES_PER_FACE ;
      volume = (avg_thick * f->area) ;
      total += volume ;
      n += 1.0 ;
    }

  total /= 2 ;  // average of white and pial surface areas
  if (n > 0)
    mean = total / n ;
  else
    mean = 0 ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  return(total) ;
}


double
MRIScomputeAbsoluteCurvature(MRI_SURFACE *mris)
{
  VERTEX    *v ;
  int       vno ;
  double    total, n ;

  for (n = total = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      total += (fabs(v->curv) * v->area) ;
      n += 1.0 ;
    }

  return(total/n) ;
}

#if 0
int
MRISreadAnnotFile(MRI_SURFACE *mris, char *name)
{
  int    i,j,vno,num;
  FILE   *fp;
  char   fname[100] ;

  MRISbuildFileName(mris, name, fname) ;
  fp = fopen(fname,"r");
  if (fp==NULL)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE, "Annotion file %s not found",fname));
  for (vno=0;vno<mris->nvertices;vno++)
    mris->vertices[vno].marked = -1 ;
  num = freadInt(fp);
  /*  printf("num=%d\n",num);*/
  for (j=0;j<num;j++)
    {
      vno = freadInt(fp);
      i = freadInt(fp);
      if (vno>=mris->nvertices||vno<0)
        printf("vertex index out of range: %d i=%d\n",vno,i);
      else
        mris->vertices[vno].marked = i;
    }
  fclose(fp);
  return(NO_ERROR) ;
}

int
MRISripVerticesWithMark(MRI_SURFACE *mris, int mark)
{
  int      vno ;
  VERTEX   *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (v->marked == mark)
        v->ripflag = 1 ;
    }
  MRISripFaces(mris) ;
  return(NO_ERROR) ;
}
int
MRISripVerticesWithoutMark(MRI_SURFACE *mris, int mark)
{
  int      vno ;
  VERTEX   *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (v->marked != mark)
        v->ripflag = 1 ;
    }
  MRISripFaces(mris) ;
  return(NO_ERROR) ;
}

int
MRISreplaceMarks(MRI_SURFACE *mris, int in_mark, int out_mark)
{
  int      vno ;
  VERTEX   *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (v->marked == in_mark)
        v->marked = out_mark ;
    }
  return(NO_ERROR) ;
}

#endif

int
MRISripVerticesWithAnnotation(MRI_SURFACE *mris, int annotation)
{
  int      vno ;
  VERTEX   *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (v->annotation == annotation)
        v->ripflag = 1 ;
    }
  MRISripFaces(mris) ;
  return(NO_ERROR) ;
}
int
MRISripVerticesWithoutAnnotation(MRI_SURFACE *mris, int annotation)
{
  int      vno ;
  VERTEX   *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
#if 0
      if (v->ripflag)
        continue ;
#endif
      if (v->annotation != annotation)
        v->ripflag = 1 ;
      else
        v->ripflag = 0 ;
    }
  MRISripFaces(mris) ;
  return(NO_ERROR) ;
}

int
MRISreplaceAnnotations(MRI_SURFACE *mris,
                       int in_annotation,
                       int out_annotation)
{
  int      vno ;
  VERTEX   *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (v->annotation == in_annotation)
        v->annotation = out_annotation ;
    }
  return(NO_ERROR) ;
}

int
MRISrestoreSurface(MRI_SURFACE *mris)
{
  int      vno, fno ;
  VERTEX   *v ;
  FACE     *f ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      v->ripflag = 0 ;
    }
  for (fno = 0 ; fno < mris->nfaces ; fno++)
    {
      f = &mris->faces[fno] ;
      f->ripflag = 0 ;
    }
  return(NO_ERROR) ;
}

int
MRIScountVertices(MRI_SURFACE *mris)
{

  int vno, n;

  n = 0;
  for(vno = 0;vno < mris->nvertices;vno++)
    if(!mris->vertices[vno].ripflag)
      n++;

  return(n);

}

