

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "utils.h"
#include "gca.h"
#include "cma.h"
#include "mrinorm.h"

char         *Progname ;

static char *normalized_transformed_sample_fname = NULL ;
static char *mask_fname = NULL ;
static char *sample_fname = NULL ;
static char *ctl_point_fname = NULL ;
static int novar = 0 ;

static float min_prior = 0.6 ;
static FILE *diag_fp = NULL ;


static int get_option(int argc, char *argv[]) ;

static char *renormalization_fname = NULL ;
static double TR = 0.0, TE = 0.0, alpha = 0.0 ;
static char *tissue_parms_fname = NULL ;
static char *example_T1 = NULL ;
static char *example_segmentation = NULL ;
static double min_region_prior(GCA *gca, int xp, int yp, int zp, int wsize, int label) ;
static GCA_SAMPLE *find_control_points(GCA *gca, GCA_SAMPLE *gcas, int total_nsamples, 
                                       int *pnorm_samples, int nregions, int label,
                                       MRI *mri_in, TRANSFORM *transform, double min_prior,
                                       double ctrl_point_pct) ;

static GCA_SAMPLE *gcas_concatenate(GCA_SAMPLE *gcas1, GCA_SAMPLE *gcas2, int n1, int n2);
static int  gcas_bounding_box(GCA_SAMPLE *gcas, int nsamples, int *pxmin, int *pymin, int *pzmin, 
                              int *pxmax, int *pymax, int *pzmax, int label) ;
static int  uniform_region(MRI *mri, int x, int y, int z, int wsize, GCA_SAMPLE *gcas) ;

/* 
   command line consists of three inputs:

   argv[1]  - input volume
   argv[2]  - atlas (gca)
   argv[3]  - transform (lta/xfm/m3d)
   argv[4]  - output volume
*/

#define DEFAULT_CTL_POINT_PCT   .25
static double ctl_point_pct = DEFAULT_CTL_POINT_PCT ;

static int normalization_structures[] =
{
  Left_Cerebral_White_Matter,
  Right_Cerebral_White_Matter,
  Left_Cerebellum_White_Matter,
  Right_Cerebellum_White_Matter,
  Brain_Stem
} ;

#define NSTRUCTURES (sizeof(normalization_structures) / sizeof(normalization_structures[0]))

static int nregions = 3 ;  /* divide each struct into 3x3x3 regions */

int
main(int argc, char *argv[])
{
  char         *gca_fname, *in_fname, *out_fname, **av, *xform_fname ;
  MRI          *mri_in, *mri_norm = NULL ;
  GCA          *gca ;
  int          ac, nargs, nsamples, msec, minutes, seconds, i, struct_samples, norm_samples, n ;
  struct timeb start ;
  GCA_SAMPLE   *gcas, *gcas_norm = NULL, *gcas_struct ;
  TRANSFORM    *transform = NULL ;


  Progname = argv[0] ;

  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    ErrorExit(ERROR_BADPARM, 
              "usage: %s <in brain> <template> <output file name>\n",
              Progname) ;

  in_fname = argv[1] ;
  gca_fname = argv[2] ;
  xform_fname = argv[3] ;
  out_fname = argv[4] ;

  TimerStart(&start) ;
  printf("reading '%s'...\n", gca_fname) ;
  fflush(stdout) ;
  gca = GCAread(gca_fname) ;
  if (gca == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open GCA %s.\n",
              Progname, gca_fname) ;

  printf("reading transform from '%s'...\n", xform_fname) ;
  fflush(stdout) ;
  transform = TransformRead(xform_fname) ;
  if (!transform)
    ErrorExit(ERROR_BADPARM, "%s: could not open xform file %s", Progname,xform_fname) ;

  if (novar)
    GCAunifyVariance(gca) ;

  if (renormalization_fname)
  {
    FILE   *fp ;
    int    *labels, nlines, i ;
    float  *intensities, f1, f2 ;
    char   *cp, line[STRLEN] ;

    fp = fopen(renormalization_fname, "r") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not read %s",
                Progname, renormalization_fname) ;

    cp = fgetl(line, 199, fp) ;
    nlines = 0 ;
    while (cp)
    {
      nlines++ ;
      cp = fgetl(line, 199, fp) ;
    }
    rewind(fp) ;
    printf("reading %d labels from %s...\n", nlines,renormalization_fname) ;
    labels = (int *)calloc(nlines, sizeof(int)) ;
    intensities = (float *)calloc(nlines, sizeof(float)) ;
    cp = fgetl(line, 199, fp) ;
    for (i = 0 ; i < nlines ; i++)
    {
      sscanf(cp, "%e  %e", &f1, &f2) ;
      labels[i] = (int)f1 ; intensities[i] = f2 ;
      if (labels[i] == Left_Cerebral_White_Matter)
        DiagBreak() ;
      cp = fgetl(line, 199, fp) ;
    }
    GCArenormalizeIntensities(gca, labels, intensities, nlines) ;
    free(labels) ; free(intensities) ;
  }




  printf("reading '%s'...\n", in_fname) ;
  fflush(stdout) ;
  mri_in = MRIread(in_fname) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not open input volume %s.\n",
              Progname, in_fname) ;

  if (mask_fname)
  {
    MRI *mri_mask ;

    mri_mask = MRIread(mask_fname) ;
    if (!mri_mask)
      ErrorExit(ERROR_NOFILE, "%s: could not open mask volume %s.\n",
                Progname, mask_fname) ;

    MRImask(mri_in, mri_mask, mri_in, 0, 0) ;
    MRIfree(&mri_mask) ;
  }
  if (alpha > 0)
    mri_in->flip_angle = alpha ;
  if (TR > 0)
    mri_in->tr = TR ;
  if (TE > 0)
    mri_in->te = TE ;

  if (example_T1)
  {
    MRI *mri_T1, *mri_seg ;

    mri_seg = MRIread(example_segmentation) ;
    if (!mri_seg)
      ErrorExit(ERROR_NOFILE,"%s: could not read example segmentation from %s",
                Progname, example_segmentation) ;
    mri_T1 = MRIread(example_T1) ;
    if (!mri_T1)
      ErrorExit(ERROR_NOFILE,"%s: could not read example T1 from %s",
                Progname, example_T1) ;
    printf("scaling atlas intensities using specified examples...\n") ;
    MRIeraseBorderPlanes(mri_seg) ;
    GCArenormalizeToExample(gca, mri_seg, mri_T1) ;
    MRIfree(&mri_seg) ; MRIfree(&mri_T1) ;
  }

  if (tissue_parms_fname)   /* use FLASH forward model */
    GCArenormalizeToFlash(gca, tissue_parms_fname, mri_in) ;

#if 0
  GCAhistoScaleImageIntensities(gca, mri_in) ;
#endif
  gcas = GCAfindAllSamples(gca, &nsamples) ;
  printf("using %d sample points...\n", nsamples) ;
  GCAcomputeSampleCoords(gca, mri_in, gcas, nsamples, transform) ;
  if (sample_fname)
    GCAtransformAndWriteSamples(gca, mri_in, gcas, nsamples, sample_fname, transform) ;
  

  for (n = 3 ; n <= nregions ; n++)
  {
    for (norm_samples = i = 0 ; i < NSTRUCTURES ; i++)
    {
      printf("finding control points in %s....\n", cma_label_to_name(normalization_structures[i])) ;
      gcas_struct = find_control_points(gca, gcas, nsamples, &struct_samples, n,
                                        normalization_structures[i], mri_in, transform, min_prior,
                                        ctl_point_pct) ;
      if (i)
      {
        GCA_SAMPLE *gcas_tmp ;
        gcas_tmp = gcas_concatenate(gcas_norm, gcas_struct, norm_samples, struct_samples) ;
        free(gcas_norm) ;
        norm_samples += struct_samples ;
        gcas_norm = gcas_tmp ;
      }
      else
      {
        gcas_norm = gcas_struct ; norm_samples = struct_samples ;
      }
    }
    
    printf("using %d total control points for intensity normalization...\n", norm_samples) ;
    if (normalized_transformed_sample_fname)
      GCAtransformAndWriteSamples(gca, mri_in, gcas_norm, norm_samples, 
                                  normalized_transformed_sample_fname, 
                                  transform) ;

    if (mri_norm)
      MRIfree(&mri_norm) ;
    mri_norm = GCAnormalizeSamples(mri_in, gca, gcas_norm, norm_samples,
                                   transform, ctl_point_fname) ;
    if (Gdiag & DIAG_WRITE)
    {
      char fname[STRLEN] ;
      sprintf(fname, "norm%d.mgh", n) ;
      printf("writing normalized volume to %s...\n", fname) ;
      MRIwrite(mri_norm, fname) ;
      sprintf(fname, "norm_samples%d.mgh", n) ;
      GCAtransformAndWriteSamples(gca, mri_in, gcas_norm, norm_samples, 
                                  fname, transform) ;

    }
    MRIcopy(mri_norm, mri_in) ;
  }

  printf("writing normalized volume to %s...\n", out_fname) ;
  if (MRIwrite(mri_norm, out_fname)  != NO_ERROR)
    ErrorExit(ERROR_BADFILE, "%s: could not write normalized volume to %s",
              Progname, out_fname);

  MRIfree(&mri_norm) ;


#if 0
  if (gca)
    GCAfree(&gca) ;
#endif
  if (mri_in)
    MRIfree(&mri_in) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("normalization took %d minutes and %d seconds.\n", 
          minutes, seconds) ;
  if (diag_fp)
    fclose(diag_fp) ;
  exit(0) ;
  return(0) ;
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
  StrUpper(option) ;
  if (!strcmp(option, "FSAMPLES"))
  {
    sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing control points to %s...\n", sample_fname) ;
  }
  else if (!strcmp(option, "MASK"))
  {
    mask_fname = argv[2] ;
    nargs = 1 ;
    printf("using MR volume %s to mask input volume...\n", mask_fname) ;
  }
  else if (!strcmp(option, "DIAG"))
  {
    diag_fp = fopen(argv[2], "w") ;
    if (!diag_fp)
      ErrorExit(ERROR_NOFILE, "%s: could not open diag file %s for writing",
                Progname, argv[2]) ;
    printf("opening diag file %s for writing\n", argv[2]) ;
    nargs = 1 ;
  }
  else if (!strcmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    nargs = 3 ;
  }
  else if (!strcmp(option, "TR"))
  {
    TR = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TR=%2.1f msec\n", TR) ;
  }
  else if (!strcmp(option, "EXAMPLE"))
  {
    example_T1 = argv[2] ;
    example_segmentation = argv[3] ;
    printf("using %s and %s as example T1 and segmentations respectively.\n",
           example_T1, example_segmentation) ;
    nargs = 2 ;
  }
  else if (!strcmp(option, "TE"))
  {
    TE = atof(argv[2]) ;
    nargs = 1 ;
    printf("using TE=%2.1f msec\n", TE) ;
  }
  else if (!strcmp(option, "ALPHA"))
  {
    nargs = 1 ;
    alpha = RADIANS(atof(argv[2])) ;
    printf("using alpha=%2.0f degrees\n", DEGREES(alpha)) ;
  }
  else if (!strcmp(option, "NSAMPLES"))
  {
    normalized_transformed_sample_fname = argv[2] ;
    nargs = 1 ;
    printf("writing  transformed normalization control points to %s...\n", 
            normalized_transformed_sample_fname) ;
  }
  else if (!strcmp(option, "RENORM"))
  {
    renormalization_fname = argv[2] ;
    nargs = 1 ;
    printf("renormalizing using predicted intensity values in %s...\n",
           renormalization_fname) ;
  }
  else if (!strcmp(option, "FLASH"))
  {
    tissue_parms_fname = argv[2] ;
    nargs = 1 ;
    printf("using FLASH forward model and tissue parms in %s to predict"
           " intensity values...\n", tissue_parms_fname) ;
  }
  else if (!strcmp(option, "PRIOR"))
  {
    min_prior = atof(argv[2]) ;
    nargs = 1 ;
    printf("using prior threshold %2.2f\n", min_prior) ;
  }
  else if (!stricmp(option, "NOVAR"))
  {
    novar = 1 ;
    printf("not using variance estimates\n") ;
  }
  else switch (*option)
  {
  case 'W':
    Gdiag |= DIAG_WRITE ;
    break ;
  case 'N':
    nregions = atoi(argv[2]) ;
    printf("using %d regions/struct for normalization\n", nregions) ;
    nargs = 1 ;
    break ;
  case 'F':
    ctl_point_fname = argv[2] ;
    nargs = 1 ;
    printf("reading manually defined control points from %s\n", ctl_point_fname) ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    printf("usage: %s <in volume> <atlas> <transform> <normalized volume>\n", 
           argv[0]) ;
    exit(1) ;
    break ;
  case 'P':
    ctl_point_pct = atof(argv[2]) ;
    nargs = 1 ;
    printf("using top %2.1f%% wm points as control points....\n",
           100.0*ctl_point_pct) ;
    break ;
  default:
    printf("unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static GCA_SAMPLE *
find_control_points(GCA *gca, GCA_SAMPLE *gcas_total,
                    int total_samples, int *pnorm_samples, int nregions, int label,
                    MRI *mri_in, TRANSFORM *transform, double min_prior, double ctl_point_pct)
{
  int        i, j, *ordered_indices, nsamples, xmin, ymin, zmin, xmax, ymax, zmax, xv,yv,zv,
             x, y, z, xi, yi, zi, region_samples, used_in_region, wsize=3, histo_peak ;
  GCA_SAMPLE *gcas, *gcas_region, *gcas_norm ;
  double     mean, var, val ;
  HISTOGRAM  *histo, *hsmooth ;
  GC1D       *gc ;

  histo = HISTOalloc(256) ; hsmooth = HISTOalloc(256) ;
  for (nsamples = i = 0 ; i < total_samples ; i++)
  {
    if (gcas_total[i].label != label)
      continue ;
    nsamples++ ;
  }

  *pnorm_samples = 0 ;
  printf("found %d control points for structure...\n", nsamples) ;
  gcas = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE)) ;
  gcas_region = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE)) ;
  gcas_norm = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE)) ;
  if (!gcas || !gcas_region || !gcas_norm)
    ErrorExit(ERROR_NOMEMORY, "find_control_points: could not allocate %d samples\n",nsamples);

  for (j = i = 0 ; i < total_samples ; i++)
  {
    if (gcas_total[i].label != label)
      continue ;
    memmove(&gcas[j], &gcas_total[i], sizeof(GCA_SAMPLE)) ;
    j++ ;
  }
  ordered_indices = (int *)calloc(nsamples, sizeof(int)) ;

  gcas_bounding_box(gcas, nsamples, &xmin, &ymin, &zmin, &xmax, &ymax, &zmax, label) ;
  printf("bounding box (%d, %d, %d) --> (%d, %d, %d)\n",
         xmin, ymin, zmin, xmax, ymax, zmax) ;
  for (x = 0 ; x < nregions ; x++)
  {
    for (y = 0 ; y < nregions ; y++)
    {
      for (z = 0 ; z < nregions ; z++)
      {
        /* only process samples in this region */
        for (region_samples = i = 0 ; i < nsamples ; i++)
        {
          xi = (int)(nregions*(gcas[i].x - xmin) / (xmax-xmin+1)) ;
          yi = (int)(nregions*(gcas[i].y - ymin) / (ymax-ymin+1)) ;
          zi = (int)(nregions*(gcas[i].z - zmin) / (zmax-zmin+1)) ;
          if ((xi < 0 || xi >= nregions) ||
              (yi < 0 || yi >= nregions) ||
              (zi < 0 || zi >= nregions))
            DiagBreak() ;
          xv = gcas[i].x ; yv = gcas[i].y ; zv = gcas[i].z ;
          if (xv == Gx && yv == Gy && zv == Gz)
            DiagBreak() ;
          if (sqrt(SQR(xv-Gx)+SQR(yv-Gy)+SQR(zv-Gz)) < 4)
            DiagBreak() ;
          if (xi != x || yi != y || zi != z || gcas[i].prior < min_prior)
            continue ;

          if (min_region_prior(gca, gcas[i].xp, gcas[i].yp, gcas[i].zp, wsize, label) < min_prior)
            continue ;
          if (uniform_region(mri_in, xv, yv, zv, wsize, &gcas[i]) == 0)
            continue ;
          memmove(&gcas_region[region_samples], &gcas[i], sizeof(GCA_SAMPLE)) ;
          region_samples++ ;
          if (gcas[i].x == Gx && gcas[i].y == Gy && gcas[i].z == Gz)
            DiagBreak() ;
        }

        printf("\t%d total samples found in region (%d, %d, %d)\n", 
               region_samples,x, y,z) ;
        if (region_samples < 8) /* can't reliably estimate statistics */
          continue ;
        HISTOclear(histo, histo) ;
        /* compute mean and variance of label within this region */
        for (mean = var = 0.0, i = 0 ; i < region_samples ; i++)
        {
          val = MRIvox(mri_in, gcas_region[i].x,gcas_region[i].y,gcas_region[i].z) ;
          histo->counts[(int)val]++ ;
          mean += val ;
          var += (val*val) ;
        }

        HISTOsmooth(histo, hsmooth, 2) ;
        histo_peak = HISTOfindLastPeakRelative(hsmooth, 3, .25) ;

          
        if (histo_peak < 0)
          continue ;

        for (mean = var = 0.0, i = 0 ; i < region_samples ; i++)
        {
          val = MRIvox(mri_in, gcas_region[i].x,gcas_region[i].y,gcas_region[i].z) ;
          mean += val ;
          var += (val*val) ;
        }

        mean /= (double)region_samples ;
        var = var / (double)region_samples - mean*mean ;
        mean = histo_peak ;
        printf("\tlabel %s: %2.1f +- %2.1f\n", cma_label_to_name(label), mean, sqrt(var)) ;

        /* ignore GCA mean and variance - use image instead (otherwise bias field will mess us up) */
        for (i = 0 ; i < region_samples ; i++)
        {
          gcas_region[i].mean = mean ;
          /*          gcas_region[i].var = var ;*/
        }

        GCAcomputeLogSampleProbability(gca, gcas_region, mri_in, transform, region_samples) ;
        GCArankSamples(gca, gcas_region, region_samples, ordered_indices) ;
#if 0
        /* use detected peak as normalization value for whole region */
        used_in_region = 1 ; j = ordered_indices[0] ;
        MRIvox(mri_in, gcas_region[j].x, gcas_region[j].y, gcas_region[j].z) = histo_peak ;
        memmove(&gcas_norm[*pnorm_samples], &gcas_region[j], sizeof(GCA_SAMPLE)) ;
        (*pnorm_samples)++ ;
#else
#if 1
        GCAremoveOutlyingSamples(gca, gcas_region, mri_in, transform, region_samples, 2.0) ;
#endif
        for (used_in_region = i = 0 ; i < region_samples ; i++)
        {
          j = ordered_indices[i] ;
          if (gcas_region[j].label != label)  /* it was an outlier */
            continue ;
          memmove(&gcas_norm[*pnorm_samples], &gcas_region[j], sizeof(GCA_SAMPLE)) ;
          (*pnorm_samples)++ ; used_in_region++ ;
        }
        if ((used_in_region <= 0) && region_samples>0)
        {
          j = ordered_indices[0] ;
          /*          gcas_region[j].label = label ;*/
          printf("forcing use of sample %d @ (%d, %d, %d)\n", j,
                 gcas_region[j].x, gcas_region[j].y, gcas_region[j].z) ;
          memmove(&gcas_norm[*pnorm_samples], &gcas_region[j], sizeof(GCA_SAMPLE)) ;
          (*pnorm_samples)++ ; used_in_region++ ;
        }
#endif
        printf("\t%d samples used in region\n", used_in_region) ;
      }
    }
  }

  /* put gca means back into samples */
  for (i = 0 ; i < *pnorm_samples ; i++)
  {
    gc = GCAfindPriorGC(gca, gcas_norm[i].xp, gcas_norm[i].yp, gcas_norm[i].zp,
                        gcas_norm[i].label) ; 
    if (gc)
    {
      gcas_norm[i].mean = gc->mean ;
      gcas_norm[i].var = gc->var ;
    }
  }
  HISTOfree(&histo) ; HISTOfree(&hsmooth) ;
  free(gcas_region) ;
  free(gcas) ;
  return(gcas_norm) ;
}

static GCA_SAMPLE *
gcas_concatenate(GCA_SAMPLE *gcas1, GCA_SAMPLE *gcas2, int n1, int n2)
{
  GCA_SAMPLE *gcas ;
  int        i ;

  gcas = (GCA_SAMPLE *)calloc(n1+n2, sizeof(GCA_SAMPLE)) ;
  if (!gcas)
    ErrorExit(ERROR_NOMEMORY, "gcas_concatenate: could not allocate %d samples",n1+n2) ;

  for (i = 0 ; i < n1 ; i++)
    memmove(&gcas[i], &gcas1[i], sizeof(GCA_SAMPLE)) ;
  for (i = 0 ; i < n2 ; i++)
    memmove(&gcas[i+n1], &gcas2[i], sizeof(GCA_SAMPLE)) ;

  return(gcas) ;
}

static int
gcas_bounding_box(GCA_SAMPLE *gcas, int nsamples, int *pxmin, int *pymin, int *pzmin, 
                  int *pxmax, int *pymax, int *pzmax, int label)
{
  int   i, xmin, ymin, zmin, xmax, ymax, zmax ;


  xmax = ymax = zmax = -1 ;
  xmin = ymin = zmin = 1000000 ;
  for (i = 0 ; i < nsamples ; i++)
  {
    if (gcas[i].x < xmin)
      xmin = gcas[i].x ;
    if (gcas[i].y < ymin)
      ymin = gcas[i].y ;
    if (gcas[i].z < zmin)
      zmin = gcas[i].z ;

    if (gcas[i].x > xmax)
      xmax = gcas[i].x ;
    if (gcas[i].y > ymax)
      ymax = gcas[i].y ;
    if (gcas[i].z > zmax)
      zmax = gcas[i].z ;
  }

  *pxmin = xmin ; *pymin = ymin ; *pzmin = zmin ;
  *pxmax = xmax ; *pymax = ymax ; *pzmax = zmax ;
  return(NO_ERROR) ;
}

static double
min_region_prior(GCA *gca, int xp, int yp, int zp, int wsize, int label)
{
  int       whalf, xi, yi, zi, xk, yk, zk ;
  double    min_prior, prior ;
  GCA_PRIOR *gcap ;

  min_prior = 1.0 ; whalf = (wsize-1)/2 ;
  for (xi = -whalf ; xi <= whalf ; xi++)
  {
    xk = xp+xi ;
    if (xk < 0 || xk >= gca->prior_width)
      continue ;
    for (yi = -whalf ; yi <= whalf ; yi++)
    {
      yk = yp+yi ;
      if (yk < 0 || yk >= gca->prior_height)
        continue ;
      for (zi = -whalf ; zi <= whalf ; zi++)
      {
        zk = zp+zi ;
        if (zk < 0 || zk >= gca->prior_depth)
          continue ;
        gcap = &gca->priors[xk][yk][zk] ;
        prior = getPrior(gcap, label) ;
        if (prior < min_prior)
          min_prior = prior ;
      }
    }
  }

  return(min_prior) ;
}

static int
uniform_region(MRI *mri, int x, int y, int z, int wsize, GCA_SAMPLE *gcas)
{
  int   xk, yk, zk, whalf, xi, yi, zi ;
  float val0, val, sigma ;

  whalf = (wsize-1)/2 ;
  sigma = sqrt(gcas->var) ;
  val0 = (float)MRIvox(mri, x, y, z) ;
  if (sigma < 0.05*val0)   /* don't let it be too small */
    sigma = 0.05*val0 ;
  if (sigma > 0.1*val0)    /* don't let it be too big */
    sigma = 0.1*val0 ;

  for (xk = -whalf ; xk <= whalf ; xk++)
  {
    xi = mri->xi[x+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->yi[y+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++)
      {
        zi = mri->zi[z+zk] ;
        val = MRIvox(mri, xi, yi,  zi) ;
        if (fabs(val-val0) > 1.5*sigma)
          return(0) ;
      }
    }
  }
  
  return(1) ;
}

