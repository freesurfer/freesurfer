/**
 * @brief extracts an array ("a variable") from surface-registration template
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
#include "label.h"
#include "mri_identify.h"
#include "fsinit.h"
#include "mri2.h"


int main(int argc, char *argv[]) ;

static MRI_SP *mrispComputeCorrelations(MRI_SP *mrisp, MRI_SP *mrisp_contra);
MRI *mrisComputeLabelCorrelations(MRI_SURFACE *mris, LABEL *area, MRI *mri_overlay, MRI *mri_corr) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static MRI_SURFACE *mris_contra = NULL ;
static MRI *mri_contra_overlay = NULL ;

static int coords = -1 ;
static int compute_corr = 0 ;
static int spherical_corr = 0 ;
static char *clabel_fname = NULL;  // contra hemi label
static char *label_fname = NULL;
static char *seed_label_fname = NULL;
static int normalize = 0 ;
static int navgs = 0 ;
static float sigma=0;
static int barycentric = 0 ;

static int frame_to_read = -1 ;

static char subjects_dir[STRLEN] ;

static float scale = 1 ;
static int nframes = 1 ;

int
main(int argc, char *argv[])
{
  char         **av, *out_fname;
  int          ac, nargs, file_type ;
  char         *in_surf, *in_overlay;
  MRI_SURFACE  *mris;
  MRI_SP       *mrisp, *mrisp_contra = NULL ;
  MRI          *mri_overlay ;

  nargs = handleVersionOption(argc, argv, "mrisp_write");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  FSinit() ;
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
  {
    usage_exit() ;
  }

  in_surf = argv[1] ;
  in_overlay = argv[2] ;
  out_fname = argv[3] ;

  fprintf(stderr, "reading surface from %s...\n", in_surf) ;
  mris = MRISread(in_surf) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, in_surf) ;
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;

  if (coords >= 0)
  {
    MRISsaveVertexPositions(mris, TMP_VERTICES) ;
    if (MRISreadVertexPositions(mris, in_overlay) != NO_ERROR)
      ErrorExit(Gerror, NULL) ;
    MRISsaveVertexPositions(mris, coords) ;
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    mrisp = MRIScoordsToParameterization(mris, NULL, scale, coords) ;
    printf("writing coordinate parameterization to %s\n", out_fname);
    MRISPwrite(mrisp,out_fname) ;
    exit(0) ;
  }
  file_type = mri_identify(in_overlay);
  if (file_type == MRI_MGH_FILE || file_type == NIFTI1_FILE || file_type == NII_FILE)
  {
    int   frame ;
    LABEL *area, *carea ;

    if (frame_to_read >= 0)
      mri_overlay = MRIreadEx(in_overlay, frame_to_read) ;
    else
      mri_overlay = MRIread(in_overlay) ;
    if (mri_overlay == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface-encoded volume file from %s", Progname, in_overlay) ;
    {
      MRI *mri_tmp ;
      int reshapefactor = mri_overlay->height * mri_overlay->depth;
      
      mri_tmp =mri_reshape(mri_overlay, reshapefactor * mri_overlay->width, 1, 1, mri_overlay->nframes);
      MRIfree(&mri_overlay) ;
      mri_overlay = mri_tmp ;
    }

    if (compute_corr) // store a frame of correlations for each vertex in the label in the mrisp
    {
      MRI *mri_corr ;

      area = LabelRead(NULL, seed_label_fname) ;
      if (area == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read label file from %s", seed_label_fname) ;
      printf("computing %d surface correlations\n", area->n_points) ;
      mri_corr = mrisComputeLabelCorrelations(mris, area, mri_overlay, NULL) ;
      MRIfree(&mri_overlay) ;
      mri_overlay = mri_corr ;
    }

    printf("processing surface-encoded volume file with %d frames\n", mri_overlay->nframes) ;
    if (mris_contra)  // compute cross-hemi correlations also
      mrisp_contra = MRISPalloc(scale, mri_contra_overlay->nframes) ;
    mrisp = MRISPalloc(scale, mri_overlay->nframes) ;
    if (label_fname)
    {
      area = LabelRead(NULL, label_fname) ;
      if (area == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read label file from %s", label_fname) ;
    }
    else
      area = NULL ;

    if (clabel_fname)
    {
      carea = LabelRead(NULL, clabel_fname) ;
      if (carea == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read label file from %s", clabel_fname) ;
    }
    else
      carea = NULL ;

    if (mri_overlay->width != mris->nvertices)
      ErrorExit(ERROR_UNSUPPORTED, "%s: overlay width (%d) does not match number of vertices in surface (%d)", Progname,mri_overlay->width, mris->nvertices);

    for (frame = 0 ; frame < mri_overlay->nframes ; frame++)
    {
      printf("\rframe %3.3d of %3.3d", frame, mri_overlay->nframes) ;
      MRISimportValFromMRI(mris, mri_overlay, frame);
      MRIScopyValuesToCurvature(mris) ;
      if (normalize)
	MRISnormalizeCurvature(mris, NORM_MEAN);

      if (label_fname)  // if compute_corr==TRUE then label is a set of seeds and don't mask
	LabelMaskSurfaceCurvature(area, mris);

      MRISaverageCurvatures(mris, navgs) ;
      if (label_fname)  // if compute_corr==TRUE then label is a set of seeds and don't mask
	LabelMaskSurfaceCurvature(area, mris);

      if (frame == Gz)
      {
	extern int DEBUG_U, DEBUG_V ;
	DEBUG_U = Gx ; DEBUG_V = Gy ;
      }
      if (barycentric)
	MRIStoParameterizationBarycentric(mris, mrisp, scale, frame) ;
      else
	MRIStoParameterization(mris, mrisp, scale, frame) ;
      if (mrisp_contra)  // doing cross-hemi correlations
      {
	MRISimportValFromMRI(mris_contra, mri_contra_overlay, frame);
	MRIScopyValuesToCurvature(mris_contra) ;
	if (normalize)
	  MRISnormalizeCurvature(mris_contra, NORM_MEAN);

	if (clabel_fname)
	  LabelMaskSurfaceCurvature(carea, mris_contra);

	MRISaverageCurvatures(mris_contra, navgs) ;

	if (clabel_fname)
	  LabelMaskSurfaceCurvature(carea, mris_contra);

	if (barycentric)
	  MRIStoParameterizationBarycentric(mris_contra, mrisp_contra, scale, frame) ;
	else
	  MRIStoParameterization(mris_contra, mrisp_contra, scale, frame) ;
      }
    }
    if (spherical_corr)
    {
      MRI_SP *mrisp_sphere ;
      mrisp_sphere = mrispComputeCorrelations(mrisp, mrisp_contra);
      MRISPfree(&mrisp) ;
      mrisp = mrisp_sphere ;
    }
    printf("\n") ;
  }
  else // process a 'curvature' file like thickness with a single frame
  {
    LABEL *area, *carea ;
    mrisp = MRISPalloc(scale, 1) ;
    file_type = mri_identify(in_overlay);
    if (file_type == MGH_LABEL_FILE)  // read in a label and create an overlay from it
    {
      int n ;

      printf("reading input label from %s\n", in_overlay) ;
      area = LabelRead(NULL, in_overlay) ;
      if (area == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not load label file %s to generate overlay",
		  in_overlay) ;
      if (LabelMaxStat(area) > 0)
      {
	printf("max stat %2.2f - transferring stats\n", LabelMaxStat(area)) ;
	for (n = 0 ; n < area->n_points ; n++)
	  mris->vertices[area->lv[n].vno].curv = area->lv[n].stat ;
      }
      else
      {
	printf("max stat %2.2f - setting to 1\n", LabelMaxStat(area)) ;
	for (n = 0 ; n < area->n_points ; n++)
	  mris->vertices[area->lv[n].vno].curv = 1 ;
      }

      LabelFree(&area) ;
    }
    else
    {
      printf("reading overlay from %s\n", in_overlay) ;
      if (MRISreadCurvatureFile(mris, in_overlay) != NO_ERROR)
	ErrorExit(ERROR_NOFILE, "%s: could not read input overlay %s", Progname, in_overlay) ;
    }

    if (label_fname)
    {
      area = LabelRead(NULL, label_fname) ;
      if (area == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read label file from %s", label_fname) ;
      LabelMaskSurfaceCurvature(area, mris);
    }
    if (clabel_fname)
    {
      if (mris_contra == NULL)
	ErrorExit(ERROR_NOFILE, "no contra surface specified (use -contra or remove -c <label>)");

      carea = LabelRead(NULL, clabel_fname) ;
      if (carea == NULL)
	ErrorExit(ERROR_NOFILE, "%s: could not read label file from %s", label_fname) ;
      LabelMaskSurfaceCurvature(carea, mris_contra);
    }
    
    if (normalize)
      MRISnormalizeCurvature(mris, NORM_MEAN);
    
    MRISaverageCurvatures(mris, navgs) ;
    if (label_fname)
    {
      LabelMaskSurfaceCurvature(area, mris);
      LabelFree(&area) ;
    }
    if (clabel_fname)
    {
      LabelMaskSurfaceCurvature(carea, mris_contra);
      LabelFree(&carea) ;
    }
    
    MRIStoParameterization(mris, mrisp, scale, 0) ;
  }
  if (sigma> 0)
  {
    int f ;
    MRI_SP *mrisp_dst ;

    printf("applying spherical convolution with sigma = %2.1f\n", sigma) ;
    mrisp_dst = MRISPclone(mrisp) ;
    for (f = 0 ; f < mrisp->Ip->num_frame ; f++)
      MRISPblur(mrisp, mrisp_dst, sigma, f) ;
    MRISPfree(&mrisp) ;
    mrisp = mrisp_dst ;
  }
  printf("writing output file to %s\n", out_fname) ;

  MRISPwrite(mrisp, out_fname) ;
  
  MRISPfree(&mrisp) ;
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
  int    nargs = 0 ;
  char   *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "SDIR"))
  {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  }
  else if (!stricmp(option, "COORDS"))
  {
    if (!stricmp(argv[2], "white"))
      coords = WHITE_VERTICES ;
    else if (!stricmp(argv[2], "pial"))
      coords = PIAL_VERTICES ;
    else
      ErrorExit(ERROR_UNSUPPORTED, "Unknown coords value %s", argv[2]) ;
    nargs = 1 ;
    printf("writing coords %s (%d) into parameterizaton\n", argv[2], coords) ;
  }
  else if (!stricmp(option, "DEBUG_VOXEL"))
  {
    extern int DEBUG_U, DEBUG_V ;
    DEBUG_U = Gx = atoi(argv[2]) ;
    DEBUG_V = Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "CONTRA"))
  {
    mris_contra = MRISread(argv[2]) ;
    if (mris_contra == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read contra sphere from %s",
		Progname, argv[3]) ;
    mri_contra_overlay = MRIread(argv[3]) ;
    if (mri_contra_overlay == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read contra overlay volume from %s",
		Progname, argv[3]) ;
    nargs = 2 ;
    printf("reading contra hemi overlay and sphere for cross-hemi calculations\n");
  }
  else if (!stricmp(option, "CORR"))
  {
    seed_label_fname = argv[2] ;
    printf("computing vertex correlations inside label %s\n", seed_label_fname) ;
    compute_corr = 1 ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "SPCORR"))
  {
    spherical_corr = 1 ;
    printf("computing correlations in spherical map\n") ;
  }
  else if (!stricmp(option, "FRAME"))
  {
    frame_to_read = atoi(argv[2]) ;
    nargs = 1 ;
    printf("extracting frame %d from input volume\n", frame_to_read) ;
  }
  else if (!stricmp(option, "BARYCENTRIC") || !stricmp(option, "BARY"))
  {
    printf("computing spherical mapping using barycentric interpolation\n") ;
    barycentric = 1 ;
    nargs = 0 ;
  }
  else if (!stricmp(option, "sigma"))
  {
    sigma = atof(argv[2]) ;
    nargs = 1 ;
    printf("blurring maps with Cartesian kernel with sigma=%2.2F\n", sigma) ;
  }
  else if (!stricmp(option, "NFRAMES")) // not implemented yet
  {
    nframes = atoi(argv[2]) ;
    nargs = 1 ;
    printf("writing out %d frames - NOT IMPLEMENTED YET\n", nframes) ;
    exit(1) ;
  }
  else if (!stricmp(option, "scale"))
  {
    scale = atof(argv[2]);
    printf("scaling width/height MRISP by %2.2f\n", scale) ;
    nargs=1;
  }
  else switch (toupper(*option))
    {
    case 'L':
      label_fname = argv[2] ;
      printf("masking label %s\n", label_fname) ;
      nargs = 1; 
      break ;
    case 'C':
      clabel_fname = argv[2] ;
      printf("masking contra label %s\n", clabel_fname) ;
      nargs = 1; 
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "averaging overlay %d times...\n", navgs) ;
      break ;
    case 'N':
      normalize = 1 ;
      fprintf(stderr, "normalizing curvature by variance.\n") ;
      break ;
    case 'W':
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'H':
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
  print_usage() ;
  exit(1) ;
}

#include "mrisp_write.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mrisp_write_help_xml,
                mrisp_write_help_xml_len);
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

MRI *
mrisComputeLabelCorrelations(MRI_SURFACE *mris, LABEL *area, MRI *mri_overlay, MRI *mri_corr)
{
  int    lvno  ;

  if (mri_corr == NULL)
    mri_corr = MRIallocSequence(mris->nvertices, 1, 1, MRI_FLOAT, area->n_points);
  if (mri_corr == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d x %d correlation matrix", mris->nvertices, area->n_points);

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (lvno = 0 ; lvno < area->n_points ; lvno++)
  {
    int    vno1, vno2, t  ;
    double norm1, norm2, dot, val1, val2 ;

    vno1 = area->lv[lvno].vno ;
    for (vno2 = 0 ; vno2 < mris->nvertices ; vno2++)
    {
      dot = norm1 = norm2 =  0 ;
      for (t = 0 ; t < mri_overlay->nframes ; t++)
      {
	val1 = MRIgetVoxVal(mri_overlay, vno1, 0, 0, t) ;
	val2 = MRIgetVoxVal(mri_overlay, vno2, 0, 0, t) ;
	dot += val1*val2 ;
	norm1 += val1*val1 ;
	norm2 += val2*val2 ;
      }
      if (FZERO(norm1) || FZERO(norm2))
	dot = 0 ;
      else
	dot = dot / ((sqrt(norm1) * sqrt(norm2))) ;
      MRIsetVoxVal(mri_corr, vno2, 0, 0, lvno, dot) ;
    }
  }

  return(mri_corr) ;
}

static MRI_SP *
mrispComputeCorrelations(MRI_SP *mrisp, MRI_SP *mrisp_contra)
{
  MRI_SP *mrisp_sphere, *mrisp_debug ;
  int    nframes, x, y, width, height, t ;
  double **norms, mean, val, **cnorms, ****corrs = nullptr;

  mrisp = MRISPclone(mrisp) ;  // we will modify this one and free it later
  width = mrisp->Ip->cols ; height = mrisp->Ip->rows ;
  nframes = width*height ;
  if (mrisp_contra)
  {
    mrisp_contra = MRISPclone(mrisp_contra) ;  // we will modify this one and free it later
    nframes *= 2 ;
  }
  mrisp_sphere = MRISPalloc(mrisp->scale, nframes);
  
  mrisp_debug = NULL ;
  if (getenv("MRISP_DEBUG"))
  {
    char *cp = getenv("MRISP_DEBUG") ;
    printf("using file %s for debugging\n", cp) ;
    mrisp_debug = MRISPread(cp) ;
    if (mrisp_debug == NULL)
      ErrorExit(ERROR_NOFILE, "could not open %s", cp) ;
  }

  // first compute means and make timecourses zero mean
  norms = (double **)calloc(width, sizeof(double *));
  if (!norms)
    ErrorExit(ERROR_NOFILE, "mrispComputeCorrelations: could not allocate norm buffer", Progname);
//(48,51) and (48,115)
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      mean = 0.0 ;
      for (t = 0 ; t < mrisp->Ip->num_frame ; t++)
      {
	val = *IMAGEFseq_pix(mrisp->Ip, x, y, t) ;
	mean += val ;
      }
      mean /= mrisp->Ip->num_frame ;
      for (t = 0 ; t < mrisp->Ip->num_frame ; t++)
	*IMAGEFseq_pix(mrisp->Ip, x, y, t) -= mean ;
    }
  }

  // now compute norms of zeromean timecourses
  for (x = 0 ; x < width ; x++)
  {
    norms[x] = (double *)calloc(height, sizeof(double));
    if (!norms[x])
      ErrorExit(ERROR_NOFILE, "mrispComputeCorrelations: could not allocate norm buffer", Progname);
    for (y = 0 ; y < height ; y++)
    {
      for (t = 0 ; t < mrisp->Ip->num_frame ; t++)
      {
	val = *IMAGEFseq_pix(mrisp->Ip, x, y, t) ;
	norms[x][y] += val*val ;
      }
      norms[x][y] = sqrt(norms[x][y]) ;
    }
  }

  // now compute correlations
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (x = 0 ; x < width ; x++)
  {
    int x1, y1, y, frame, t ;
    double norm1, norm2, dot ;

    for (y = 0 ; y < height ; y++)
    {
      if (x == 48 && (y == 51 || y == 115))
	  DiagBreak() ;
      if (y == 48 && (x == 51 || x == 115))
	  DiagBreak() ;
      norm1 = norms[x][y];
      if (FZERO(norm1))
	continue ;
      for (x1 = 0 ; x1 < width ; x1++)
	for (y1 = 0 ; y1 < height ; y1++)
	{
	  double v1, v2 ;
	  if (x == 0 && y == 0 && x1 == 10 && y1 == 42) // vertices 0 and 1000 in fsaverage5
	    DiagBreak() ;
	  frame = x1*height + y1 ;  // images are stored in column-major format
	  frame = y1*width + x1 ;  // images are stored in row-major format
	  norm2 = norms[x1][y1];
	  if (FZERO(norm2))
	    continue ;
	  if (mrisp_debug && *IMAGEFseq_pix(mrisp_debug->Ip, x, y, frame)>0)
	  {
	    *IMAGEFseq_pix(mrisp_sphere->Ip, x, y, frame) = 1 ;
	    continue ;
	  }
	  for (dot = 0.0, t = 0 ; t < mrisp->Ip->num_frame ; t++)
	  {
	    v1 = *IMAGEFseq_pix(mrisp->Ip, x, y, t) ;
	    v2 = *IMAGEFseq_pix(mrisp->Ip, x1, y1, t);
	    dot += v1 * v2;
	  }
	  *IMAGEFseq_pix(mrisp_sphere->Ip, x, y, frame) = dot / (norm1*norm2);
	}
    }
  }

  if (mrisp_contra)
  {
    printf("\ncomputing cross-hemi correlations\n") ;
    // first compute means and make timecourses zero mean
    cnorms = (double **)calloc(width, sizeof(double *));
    if (!cnorms)
      ErrorExit(ERROR_NOFILE, "mrispComputeCorrelations: could not allocate norm buffer", Progname);
    for (x = 0 ; x < width ; x++)
    {
      cnorms[x] = (double *)calloc(height, sizeof(double));
      if (!cnorms[x])
	ErrorExit(ERROR_NOFILE, "mrispComputeCorrelations: could not allocate norm buffer", Progname);
      for (y = 0 ; y < height ; y++)
      {
	mean = 0.0 ;
	for (t = 0 ; t < mrisp_contra->Ip->num_frame ; t++)
	{
	  val = *IMAGEFseq_pix(mrisp_contra->Ip, x, y, t) ;
	  mean += val ;
	  cnorms[x][y] += val*val ;
	}
	cnorms[x][y] = sqrt(cnorms[x][y]) ;
	mean /= mrisp_contra->Ip->num_frame ;
	for (t = 0 ; t < mrisp_contra->Ip->num_frame ; t++)
	  *IMAGEFseq_pix(mrisp_contra->Ip, x, y, t) -= mean ;
      }
    }

    // now compute norms from zero mean timecourses
    cnorms = (double **)calloc(width, sizeof(double *));
    if (!cnorms)
      ErrorExit(ERROR_NOFILE, "mrispComputeCorrelations: could not allocate norm buffer", Progname);
    for (x = 0 ; x < width ; x++)
    {
      cnorms[x] = (double *)calloc(height, sizeof(double));
      if (!cnorms[x])
	ErrorExit(ERROR_NOFILE, "mrispComputeCorrelations: could not allocate norm buffer", Progname);
      for (y = 0 ; y < height ; y++)
      {
	for (t = 0 ; t < mrisp_contra->Ip->num_frame ; t++)
	{
	  val = *IMAGEFseq_pix(mrisp_contra->Ip, x, y, t) ;
	  cnorms[x][y] += val*val ;
	}
	cnorms[x][y] = sqrt(cnorms[x][y]) ;
      }
    }
    corrs = (double ****)calloc(width, sizeof(double)) ;
    for (x = 0 ; x < width ; x++)
    {
      corrs[x] = (double ***)calloc(height, sizeof(double));
      for (y = 0 ; y < height ; y++)
      {
	int x1, y1, frame ;
	corrs[x][y] = (double **)calloc(2*width, sizeof(double));
	for (frame = x1 = 0 ; x1 < 2*width ; x1++)
	{
	  corrs[x][y][x1] = (double *)calloc(height, sizeof(double));
	  for (y1 = 0 ; y1 < height ; y1++, frame++)
	    corrs[x][y][x1][y1] = *IMAGEFseq_pix(mrisp_sphere->Ip, x, y, frame);
	}
      }
    }
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (x = 0 ; x < width ; x++)
    {
      int x1, y1, y, frame, t ;
      double norm1, norm2, dot ;

      for (y = 0 ; y < height ; y++)
      {
	norm1 = norms[x][y];
	if (FZERO(norm1))
	  continue ;
	for (x1 = 0 ; x1 < width ; x1++)
	  for (y1 = 0 ; y1 < height ; y1++)
	  {
	    frame = (width*height)+x1*height + y1 ; // contra data starts after ipsi data
	    frame = (width*height)+y1*width + x1 ; // contra data starts after ipsi data
	    norm2 = cnorms[x1][y1];
	    if (FZERO(norm2))
	      continue ;
	    for (dot = 0.0, t = 0 ; t < mrisp->Ip->num_frame ; t++)
	      dot += *IMAGEFseq_pix(mrisp->Ip, x, y, t) * *IMAGEFseq_pix(mrisp_contra->Ip, x1, y1, t);
	    *IMAGEFseq_pix(mrisp_sphere->Ip, x, y, frame) = dot / (norm1*norm2);
	    corrs[x][y][x1+width][y1] = *IMAGEFseq_pix(mrisp_sphere->Ip, x, y, frame) ;
	    if (x == Gx && y == Gy)
	      DiagBreak() ;
	  }
      }
    }
    for (x = 0 ; x < width ; x++)
      free(cnorms[x]);
    free(cnorms) ;
    for (x = 0 ; x < width ; x++)
    {
      for (y = 0 ; y < height ; y++)
      {
	int x1 ;
	for (x1 = 0 ; x1 < 2*width ; x1++)
	  free(corrs[x][y][x1]) ;
	free(corrs[x][y]) ;
      }
    }
    free(corrs[x]) ;
  }
  free(corrs) ;

  for (x = 0 ; x < width ; x++)
    free(norms[x]);
  free(norms) ;
  MRISPfree(&mrisp) ; // not passed version - the one we allocated
  if (mrisp_contra)
    MRISPfree(&mrisp_contra) ; // not passed version - the one we allocated
  return(mrisp_sphere) ;
}

