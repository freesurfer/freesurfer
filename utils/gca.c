#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mri.h"
#include "error.h"
#include "diag.h"
#include "utils.h"
#include "gca.h"
#include "proto.h"
#include "fio.h"
#include "transform.h"
#include "macros.h"
#include "utils.h"
#include "cma.h"

int Ggca_label = -1 ;
int Ggca_x = -1 ;
int Ggca_y = -1 ;
int Ggca_z = -1 ;

/* this is the hack section */
double PRIOR_FACTOR = 0.1 ;
#define LABEL_UNDETERMINED   255

static int total_pruned = 0 ;


#define MIN_VAR  (2*2)   /* should make this configurable */
#define BIG_AND_NEGATIVE            -10000000.0
#define UNKNOWN_DIST                4  /* within 4 mm of some known label */
#define GCA_VERSION                  2.0
#define DEFAULT_MAX_LABELS_PER_GCAN  4

static double FLASHforwardModel(double flip_angle, double TR, double PD, 
                                double T1) ;
static int    MRIorderIndices(MRI *mri, short *x_indices, short *y_indices, 
                              short *z_indices) ;
static double gcaVoxelGibbsLogLikelihood(GCA *gca, MRI *mri_labels, 
                                         MRI *mri_inputs, int x, int y, int z, 
                                         LTA *lta, double gibbs_coef) ;
static double gcaNbhdGibbsLogLikelihood(GCA *gca, MRI *mri_labels, 
                                        MRI *mri_inputs, int x, int y, int z, 
                                        LTA *lta, double gibbs_coef) ;
static GCA_SAMPLE *gcaExtractLabelAsSamples(GCA *gca, MRI *mri_labeled, 
                                            LTA *lta, 
                                            int *pnsamples, int label) ;
int MRIcomputeVoxelPermutation(MRI *mri, short *x_indices, short *y_indices,
                               short *z_indices);
GCA *gcaAllocMax(int ninputs, float spacing, int width, int height, int depth, 
                 int max_labels) ;
static int GCAupdateNode(GCA *gca, MRI *mri, int xn, int yn, int zn,
                         float val,int label, GCA *gca_prune);
static int GCAupdateNodeGibbsPriors(GCA *gca,MRI*mri,int xn,int yn,int zn,
                              int x, int y, int z, int label);
static int different_nbr_labels(GCA *gca, int x, int y, int z, int wsize,
                                 int label) ;
static int gcaRegionStats(GCA *gca, int x0, int y0, int z0, int wsize,
                   float *priors, float *vars, float *means) ;
static int gcaFindBestSample(GCA *gca, int x, int y, int z, int best_label,
                             int wsize, GCA_SAMPLE *gcas);
static int gcaFindClosestMeanSample(GCA *gca, float mean, float min_prior,
                                    int x, int y, int z, 
                                    int label, int wsize, GCA_SAMPLE *gcas);
static GC1D *gcaFindHighestPriorGC(GCA *gca, int x, int y, int z,int label,
                                 int wsize) ;
static double gcaGibbsImageLogLikelihood(GCA *gca, MRI *mri_labels, 
                                         MRI *mri_inputs, LTA *lta) ;

static GC1D *gcaFindGC(GCA *gca, int x, int y, int z,int label) ;
static int mriFillRegion(MRI *mri, int x,int y,int z,int fill_val,int whalf);
static int gcaFindMaxPriors(GCA *gca, float *max_priors) ;
static int gcaFindIntensityBounds(GCA *gca, float *pmin, float *pmax) ;
static int dump_gcan(GCA_NODE *gcan, FILE *fp, int verbose) ;
static GCA_NODE *findSourceGCAN(GCA *gca, MRI *mri_src,LTA *lta,
                                int x,int y,int z);
static float gcaComputeLabelStats(GCA *gca, int target_label, float *pvar);
static float getLabelMean(GCA_NODE *gcan, int label, float *pvar) ;
static int   borderVoxel(MRI *mri, int x, int y, int z) ;
static int   GCAmaxLikelihoodBorderLabel(GCA *gca, MRI *mri_inputs, 
                                         MRI *mri_labels, LTA *lta,
                                         int x, int y, int z, float min_ratio);


/* arrays for indexing 6-connected neighbors */
static int xnbr_offset[] = { 1, -1, 0, 0,  0, 0, 0} ;
static int ynbr_offset[] = { 0, 0,  1, -1, 0, 0, 0} ;
static int znbr_offset[] = { 0, 0,  0, 0,  1, -1,0} ;

static int
dump_gcan(GCA_NODE *gcan, FILE *fp, int verbose)
{
  int       n, i, j ;
  GC1D      *gc ;

  for (n = 0 ; n < gcan->nlabels ; n++)
  {
    gc = &gcan->gcs[n] ;
    fprintf(fp, "%d: label %s, mean %2.1f, std %2.1f, prior %2.3f\n",
            n, cma_label_to_name(gcan->labels[n]), gc->mean, sqrt(gc->var), 
            gc->prior) ;
    if (verbose) for (i = 0 ; i < GIBBS_NEIGHBORS ; i++)
    {
      fprintf(fp, "\tnbr %d (%d,%d,%d): %d labels\n",
              i, xnbr_offset[i], ynbr_offset[i], znbr_offset[i], 
              gc->nlabels[i]) ;
      for (j = 0 ; j < gc->nlabels[i] ; j++)
      {
        fprintf(fp, "\t\tlabel %s, prior %2.3f\n", 
                cma_label_to_name(gc->labels[i][j]),
                gc->label_priors[i][j]) ;
      }
    }
  }
  return(NO_ERROR) ;
}


int
GCAvoxelToNode(GCA *gca, MRI *mri, int xv, int yv, int zv, int *pxn, 
               int *pyn, int *pzn)
{
  float   xscale, yscale, zscale ;

#if 1
  xscale = mri->xsize / gca->spacing ;
  yscale = mri->ysize / gca->spacing ;
  zscale = mri->zsize / gca->spacing ;
#else
  xscale = 1.0f / gca->spacing ;
  yscale = 1.0f / gca->spacing ;
  zscale = 1.0f / gca->spacing ;
#endif
  *pxn = (int)(xv * xscale) ;
  if (*pxn < 0)
    *pxn = 0 ;
  else if (*pxn >= gca->width)
    *pxn = gca->width-1 ;

  *pyn = (int)(yv * yscale) ;
  if (*pyn < 0)
    *pyn = 0 ;
  else if (*pyn >= gca->height)
    *pyn = gca->height-1 ;

  *pzn = (int)(zv * zscale) ;
  if (*pzn < 0)
    *pzn = 0 ;
  else if (*pzn >= gca->depth)
    *pzn = gca->depth-1 ;

  return(NO_ERROR) ;
}

int
GCAsourceVoxelToNode(GCA *gca, MRI *mri, LTA *lta,int xv, int yv, int zv, 
                     int *pxn, int *pyn, int *pzn)
{
  static VECTOR *v_input, *v_canon = NULL ;
  float   xt, yt, zt ;

  if (!v_canon)
  {
    v_input = VectorAlloc(4, MATRIX_REAL) ;
    v_canon = VectorAlloc(4, MATRIX_REAL) ;
    *MATRIX_RELT(v_input, 4, 1) = 1.0 ;
    *MATRIX_RELT(v_canon, 4, 1) = 1.0 ;
  }

  V3_X(v_input) = (float)xv*mri->xsize ; 
  V3_Y(v_input) = (float)yv*mri->ysize ; 
  V3_Z(v_input) = (float)zv*mri->zsize ; 
  MatrixMultiply(lta->xforms[0].m_L, v_input, v_canon) ;
  /* should change this to remove cast, but backwards comp for now */
  xt = nint(V3_X(v_canon)/mri->xsize) ; 
  yt = nint(V3_Y(v_canon)/mri->ysize) ; 
  zt = nint(V3_Z(v_canon)/mri->zsize) ; 
  GCAvoxelToNode(gca, mri, xt, yt, zt, pxn, pyn, pzn) ;

  return(NO_ERROR) ;
}

int
GCAnodeToVoxel(GCA *gca, MRI *mri, int xn, int yn, int zn, int *pxv, 
               int *pyv, int *pzv)
{
  float   xscale, yscale, zscale, offset ;

#if 1
  xscale = mri->xsize / gca->spacing ;
  yscale = mri->ysize / gca->spacing ;
  zscale = mri->zsize / gca->spacing ;
#else
  xscale = 1.0f / gca->spacing ;
  yscale = 1.0f / gca->spacing ;
  zscale = 1.0f / gca->spacing ;
#endif
  offset = (gca->spacing-1.0)/2.0 ;
  *pxv = (int)(xn / xscale + offset) ;
  if (*pxv < 0)
    *pxv = 0 ;
  else if (*pxv >= mri->width)
    *pxv = mri->width-1 ;

  *pyv = (int)(yn / yscale + offset) ;
  if (*pyv < 0)
    *pyv = 0 ;
  else if (*pyv >= mri->height)
    *pyv = mri->height-1 ;

  *pzv = (int)(zn / zscale + offset) ;
  if (*pzv < 0)
    *pzv = 0 ;
  else if (*pzv >= mri->depth)
    *pzv = mri->depth-1 ;

  return(NO_ERROR) ;
}


GCA  *
GCAalloc(int ninputs, float spacing, int width, int height, int depth)
{
  return(gcaAllocMax(ninputs, spacing, width,height,depth,
                     DEFAULT_MAX_LABELS_PER_GCAN));
}

GCA *
gcaAllocMax(int ninputs, float spacing, int width, int height, int depth, 
            int max_labels)
{
  GCA      *gca ;
  GCA_NODE *gcan ;
  int      x, y, z ;

  gca = calloc(1, sizeof(GCA)) ;
  if (!gca)
    ErrorExit(ERROR_NOMEMORY, "GCAalloc: could not allocate struct") ;

  gca->ninputs = ninputs ;
  gca->spacing = spacing ;

  gca->width = (int)ceil((float)width/spacing) ;
  gca->height = (int)ceil((float)height/spacing) ;
  gca->depth = (int)ceil((float)depth/spacing) ;

  gca->nodes = (GCA_NODE ***)calloc(gca->width, sizeof(GCA_NODE **)) ;
  if (!gca->nodes)
    ErrorExit(ERROR_NOMEMORY, "GCAalloc: could not allocate nodes") ;
  for (x = 0 ; x < gca->width ; x++)
  {
    gca->nodes[x] = (GCA_NODE **)calloc(gca->height, sizeof(GCA_NODE *)) ;
    if (!gca->nodes[x])
      ErrorExit(ERROR_NOMEMORY, "GCAalloc: could not allocate %dth **",x) ;

    for (y = 0 ; y < gca->height ; y++)
    {
      gca->nodes[x][y] = (GCA_NODE *)calloc(gca->depth, sizeof(GCA_NODE)) ;
      if (!gca->nodes[x][y])
        ErrorExit(ERROR_NOMEMORY,"GCAalloc: could not allocate %d,%dth *",x,y);
      for (z = 0 ; z < gca->depth ; z++)
      {
        gcan = &gca->nodes[x][y][z] ;
        gcan->max_labels = max_labels ;

        if (max_labels > 0)
        {
          /* allocate new ones */
          gcan->gcs = (GC1D *)calloc(gcan->max_labels, sizeof(GC1D)) ;
          if (!gcan->gcs)
            ErrorExit(ERROR_NOMEMORY, "GCANalloc: couldn't allocate gcs to %d",
                      gcan->max_labels) ;
          gcan->labels = (char *)calloc(gcan->max_labels, sizeof(char)) ;
          if (!gcan->labels)
            ErrorExit(ERROR_NOMEMORY,
                      "GCANalloc: couldn't allocate labels to %d",
                      gcan->max_labels) ;
        }
      }
    }
  }
  
  return(gca) ;
}

int
GCAfree(GCA **pgca)
{
  GCA  *gca ;
  int  x, y, z ;

  gca = *pgca ;
  *pgca = NULL ;

  for (x = 0 ; x < gca->width ; x++)
  {
    for (y = 0 ; y < gca->height ; y++)
    {
      for (z = 0 ; z < gca->depth ; z++)
      {
        GCANfree(&gca->nodes[x][y][z]) ;
      }
      free(gca->nodes[x][y]) ;
    }
    free(gca->nodes[x]) ;
  }

  free(gca->nodes) ;
  free(gca) ;
  return(NO_ERROR) ;
}

int
GCANfree(GCA_NODE *gcan)
{
  if (gcan->nlabels)
  {
    free(gcan->labels) ;
    free(gcan->gcs) ;
  }
  return(NO_ERROR) ;
}
int
GCAtrain(GCA *gca, MRI *mri_inputs, MRI *mri_labels, LTA *lta, GCA *gca_prune)
{
  int    x, y, z, width, height, depth, label, val, xn, yn, zn, i,
         xnbr, ynbr, znbr, xn_nbr, yn_nbr, zn_nbr ;


  /* convert transform to voxel coordinates */

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_labels->width ; height = mri_labels->height; 
  depth = mri_labels->depth ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        label = MRIvox(mri_labels, x, y, z) ;
#if 0
        if (!label)
          continue ;
#endif

        val = MRIvox(mri_inputs, x, y, z) ;

        GCAsourceVoxelToNode(gca, mri_inputs, lta, x, y, z, &xn, &yn, &zn) ;

        if (GCAupdateNode(gca, mri_inputs, xn, yn, zn, 
                          (float)val,label,gca_prune) == NO_ERROR)
          GCAupdateNodeGibbsPriors(gca, mri_labels, xn, yn, zn, x, y,z, label);

        if (xn == Ggca_x && yn == Ggca_y && zn == Ggca_z && 
            (label == Ggca_label || Ggca_label < 0))
        {
          GC1D *gc ;
          gc = gcaFindGC(gca, xn, yn, zn, label) ;
          if (gc)
            printf("voxel(%d,%d,%d) = %d --> node(%d,%d,%d), "
                   "label %s (%d), mean %2.1f\n",
                   x, y, z, val, xn, yn, zn, cma_label_to_name(label),label,
                   gc->mean / gc->prior) ;
        }
        for (i = 0 ; i < GIBBS_NEIGHBORS ; i++)
        {
          xnbr = x+xnbr_offset[i] ; 
          ynbr = y+ynbr_offset[i] ; 
          znbr = z+znbr_offset[i] ; 
          GCAsourceVoxelToNode(gca, mri_inputs, lta, 
                               xnbr, ynbr, znbr,&xn_nbr,&yn_nbr,&zn_nbr) ;
          if (xn_nbr == xn && yn_nbr == yn && zn_nbr == zn)
            continue ;  /* only update if it is a different node */

          if (GCAupdateNode(gca, mri_inputs, xn_nbr, yn_nbr, zn_nbr, 
                            (float)val,label,gca_prune) == NO_ERROR)
            GCAupdateNodeGibbsPriors(gca, mri_labels, xn_nbr, yn_nbr, zn_nbr, 
                                     x, y,z, label);
        }
      }
    }
  }

  return(NO_ERROR) ;
}
int
GCAwrite(GCA *gca, char *fname)
{
  FILE     *fp ;
  int      x, y, z, n, i, j ; 
  GCA_NODE *gcan ;
  GC1D     *gc ;

  fp  = fopen(fname, "wb") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "GCAwrite: could not open GCA %s for writing",fname)) ;

  fwriteFloat(GCA_VERSION, fp) ;
  fwriteFloat(gca->spacing, fp) ;
  fwriteInt(gca->width,fp);fwriteInt(gca->height,fp); fwriteInt(gca->depth,fp);
  fwriteInt(gca->ninputs,fp) ;

  for (x = 0 ; x < gca->width ; x++)
  {
    for (y = 0 ; y < gca->height ; y++)
    {
      for (z = 0 ; z < gca->depth ; z++)
      {
        if (x == 139 && y == 103 && z == 139)  /* wm should be pallidum */
          DiagBreak() ;
        gcan = &gca->nodes[x][y][z] ;
        fwriteInt(gcan->nlabels, fp) ;
        fwriteInt(gcan->total_training, fp) ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          fputc((int)gcan->labels[n],fp) ;
          fwriteFloat(gc->mean, fp) ;
          fwriteFloat(gc->var, fp) ;
          fwriteFloat(gc->prior, fp) ;
          for (i = 0 ; i < GIBBS_NEIGHBORS ; i++)
          {
            fwriteInt(gc->nlabels[i], fp) ;
            for (j = 0 ; j < gc->nlabels[i] ; j++)
            {
              fwriteInt((int)gc->labels[i][j], fp) ;
              fwriteFloat(gc->label_priors[i][j], fp) ;
            }
          }
        }
      }
    }
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}

GCA *
GCAread(char *fname)
{
  FILE     *fp ;
  int      x, y, z, n, i, j ; 
  GCA      *gca ;
  GCA_NODE *gcan ;
  GC1D     *gc ;
  float    version, spacing ;
  int      width, height, depth, ninputs ;

  fp  = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL,
                (ERROR_NOFILE,
                 "GCAread: could not open GCA %s for reading",fname)) ;

  version = freadFloat(fp) ;
  spacing = freadFloat(fp) ;
  width = freadInt(fp); height = freadInt(fp); depth = freadInt(fp);
  ninputs = freadInt(fp) ;

  gca = gcaAllocMax(ninputs, spacing, spacing*width, spacing*height, 
                    spacing*depth, 0) ;
  if (!gca)
    ErrorReturn(NULL, (Gdiag, NULL)) ;

  for (x = 0 ; x < gca->width ; x++)
  {
    for (y = 0 ; y < gca->height ; y++)
    {
      for (z = 0 ; z < gca->depth ; z++)
      {
        if (x == 28 && y == 39 && z == 39)
          DiagBreak() ;
        gcan = &gca->nodes[x][y][z] ;
        gcan->nlabels = freadInt(fp) ;
        gcan->total_training = freadInt(fp) ;
        gcan->labels = (char *)calloc(gcan->nlabels, sizeof(char)) ;
        if (!gcan->labels)
          ErrorExit(ERROR_NOMEMORY, "GCAread(%s): could not allocate %d "
                    "labels @ (%d,%d,%d)", fname, gcan->nlabels, x, y, z) ;
        gcan->gcs = (GC1D *)calloc(gcan->nlabels, sizeof(GC1D)) ;
        if (!gcan->gcs)
          ErrorExit(ERROR_NOMEMORY, "GCAread(%s); could not allocated %d gcs "
                    "@ (%d,%d,%d)", fname, gcan->nlabels, x, y, z) ;

        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          gcan->labels[n] = (char)fgetc(fp) ;
          gc->mean = freadFloat(fp) ;
          gc->var = freadFloat(fp) ;
          gc->prior = freadFloat(fp) ;
          for (i = 0 ; i < GIBBS_NEIGHBORS ; i++)
          {
            gc->nlabels[i] = freadInt(fp) ;

#if 1
            /* allocate new ones */
            gc->label_priors[i] = 
              (float *)calloc(gc->nlabels[i],sizeof(float));
            if (!gc->label_priors[i])
              ErrorExit(ERROR_NOMEMORY, "GCAread(%s): "
                        "couldn't expand gcs to %d", fname,gc->nlabels) ;
            gc->labels[i] = (char *)calloc(gc->nlabels[i], sizeof(char)) ;
            if (!gc->labels)
              ErrorExit(ERROR_NOMEMORY, 
                        "GCAread(%s): couldn't expand labels to %d",
                        fname, gc->nlabels[i]) ;
#endif
            for (j = 0 ; j < gc->nlabels[i] ; j++)
            {
              gc->labels[i][j] = (char)freadInt(fp) ;
              gc->label_priors[i][j] = freadFloat(fp) ;
            }
          }
        }
      }
    }
  }

  fclose(fp) ;
  return(gca) ;
}

static int
GCAupdateNode(GCA *gca, MRI *mri, int xn, int yn, int zn, float val, int label,
              GCA *gca_prune)
{
  int      n ;
  GCA_NODE *gcan ;
  GC1D     *gc ;


  if (label >= MAX_CMA_LABEL)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, 
                 "GCAupdateNode(%d, %d, %d, %2.1f, %d): label out of range",
                 xn, yn, zn, val, label)) ;

  if (xn == 23 && yn == 30 && zn == 32)
    DiagBreak() ;
  if (xn == 23 && yn == 30 && zn == 32 && label == 41)
    DiagBreak() ;
  if (xn == 23 && yn == 30 && zn == 32 && label == 41 && val < 60)
    DiagBreak() ;


  if (label > 0 && gca_prune != NULL)
  {
    GCA_NODE *gcan_prune ;
    GC1D     *gc_prune ;
    float    mean, std ;
    int      nprune ;

     gcan_prune = &gca_prune->nodes[xn][yn][zn] ;

     for (nprune = 0 ; nprune < gcan_prune->nlabels ; nprune++)
     {
       if (gcan_prune->labels[nprune] == label)
         break ;
     }
     if (nprune >= gcan_prune->nlabels)
       ErrorPrintf(ERROR_BADPARM, "WARNING: pruning GCA at (%d,%d,%d) doesn't "
                   "contain label %d", xn, yn, zn, label) ;
     gc_prune = &gcan_prune->gcs[nprune] ;

     mean = gc_prune->mean ;
     std = sqrt(gc_prune->var) ;
     if ((fabs(val - mean) / std) > 2) /* more than 2 stds from mean */
     {
       if (xn == 23 && yn == 30 && zn == 32)
       {
         printf(
                "pruning val %2.0f, label %d @ (%d,%d,%d),u=%2.1f, std=%2.1f\n"
                ,val, label, xn,yn,zn,mean, std) ;
         DiagBreak() ;
       }
       total_pruned++ ;
       return(ERROR_BAD_PARM) ;
     }
  }
    

  gcan = &gca->nodes[xn][yn][zn] ;

  for (n = 0 ; n < gcan->nlabels ; n++)
  {
    if (gcan->labels[n] == label)
      break ;
  }
  if (n >= gcan->nlabels)  /* have to allocate a new classifier */
  {

    if (n >= gcan->max_labels)
    {
      int  old_max_labels ;
      char *old_labels ;
      GC1D *old_gcs ;

      old_max_labels = gcan->max_labels ; 
      gcan->max_labels += 2 ;
      old_labels = gcan->labels ;
      old_gcs = gcan->gcs ;

      /* allocate new ones */
      gcan->gcs = (GC1D *)calloc(gcan->max_labels, sizeof(GC1D)) ;
      if (!gcan->gcs)
        ErrorExit(ERROR_NOMEMORY, "GCANupdateNode: couldn't expand gcs to %d",
                  gcan->max_labels) ;
      gcan->labels = (char *)calloc(gcan->max_labels, sizeof(char)) ;
      if (!gcan->labels)
        ErrorExit(ERROR_NOMEMORY, "GCANupdateNode: couldn't expand labels to %d",
                  gcan->max_labels) ;

      /* copy the old ones over */
      memmove(gcan->gcs, old_gcs, old_max_labels*sizeof(GC1D)) ;
      memmove(gcan->labels, old_labels, old_max_labels*sizeof(char)) ;

      /* free the old ones */
      free(old_gcs) ; free(old_labels) ;
    }
    gcan->nlabels++ ;
  }

  gc = &gcan->gcs[n] ;

  /* these will be updated when training is complete */
  gc->mean += val ; gc->var += val*val ; gc->prior += 1.0f ;
  gcan->total_training++ ;
  
  gcan->labels[n] = label ;

  return(NO_ERROR) ;
}
int
GCAcompleteTraining(GCA *gca)
{
  int      x, y, z, n, total_nodes, total_gcs, i, j ;
  float    nsamples ;
  GCA_NODE *gcan ;
  GC1D     *gc ;

  total_nodes = gca->width*gca->height*gca->depth ; total_gcs = 0 ;
  for (x = 0 ; x < gca->width ; x++)
  {
    for (y = 0 ; y < gca->height ; y++)
    {
      for (z = 0 ; z < gca->depth ; z++)
      {
        gcan = &gca->nodes[x][y][z] ;
        total_gcs += gcan->nlabels ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          float var ;

          gc = &gcan->gcs[n] ;
          nsamples = gc->prior ;  /* prior is count of # of occurences now */

          for (i = 0 ; i < GIBBS_NEIGHBORS ; i++)
          {
            for (j = 0 ; j < gc->nlabels[i] ; j++)
              gc->label_priors[i][j] /= (float)nsamples ;
          }

          gc->mean /= nsamples ;
          var = gc->var / nsamples - gc->mean*gc->mean ;
          if (var < -0.1)
            DiagBreak() ;
          if (var < MIN_VAR)
            var = MIN_VAR ;
          gc->var = var ;
          gc->prior /= (float)gcan->total_training ;
        }
        if (x == 63 && y == 107 && z == 120 && DIAG_VERBOSE_ON)
          dump_gcan(gcan, stdout, 0) ;
      }
    }
  }
  printf("%d classifiers: %2.1f per node\n", 
         total_gcs, (float)total_gcs/(float)total_nodes) ;
  if (total_pruned > 0)
  {
    printf("%d samples pruned during training\n", total_pruned) ;
    total_pruned = 0 ;
  }
  return(NO_ERROR) ;
}

MRI  *
GCAlabel(MRI *mri_inputs, GCA *gca, MRI *mri_dst, LTA *lta)
{
  int      x, y, z, width, height, depth, label, val,
           xn, yn, zn, n ;
  GCA_NODE *gcan ;
  GC1D     *gc ;
  float    dist, max_p, p ;

  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_inputs, NULL) ;
    if (!mri_dst)
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst") ;
    MRIcopyHeader(mri_inputs, mri_dst) ;
  }


  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_inputs->width ; height = mri_inputs->height; 
  depth = mri_inputs->depth ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == 95 && y == 86 && z == 99)  /* wm should be ventricle (1337)*/
          DiagBreak() ; 
        if (x == 151 && y == 117 && z == 121)  /* wm should be hippo(1428)*/
          DiagBreak() ; 
        if (x == 92 && y == 115 && z == 117)  /* was 79 */
          DiagBreak() ; 

        if (x == width/2 && y == height/2 && z == depth/2)
          DiagBreak() ;

        GCAsourceVoxelToNode(gca, mri_inputs, lta, x, y, z, &xn, &yn, &zn) ;
        val = MRIvox(mri_inputs, x, y, z) ;
        if (x == 153 && y == 119 && z == 117)  /* wm should be hippo (1484) */
        {
          Gx = xn ; Gy = yn ; Gz = zn ;
          DiagBreak() ;
        }

        gcan = &gca->nodes[xn][yn][zn] ;
        label = 0 ; max_p = BIG_AND_NEGATIVE ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          
          /* compute 1-d Mahalanobis distance */
          dist = (val-gc->mean) ;
          if (FZERO(gc->var))  /* make it a delta function */
          {
            if (FZERO(dist))
              p = 1.0 ;
            else
              p = 0.0 ;
          }
          else
            p = 1 / sqrt(gc->var * 2 * M_PI) * exp(-dist*dist/gc->var) ;

          p *= gc->prior ;
          if (p > max_p)
          {
            max_p = p ;
            label = gcan->labels[n] ;
          }
        }
        if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
        {
          printf("(%d, %d, %d): T1=%d, label %s (%d), p=%2.2e\n",
                 x, y, z, MRIvox(mri_inputs,x,y,z), 
                 cma_label_to_name(label), label, max_p) ;
        }
        MRIvox(mri_dst, x, y, z) = label ;
      }
    }
  }

  return(mri_dst) ;
}

MRI  *
GCAlabelProbabilities(MRI *mri_inputs, GCA *gca, MRI *mri_dst, LTA *lta)
{
  int      x, y, z, width, height, depth, label, val,
           xn, yn, zn, n ;
  GCA_NODE *gcan ;
  GC1D     *gc ;
  double   dist, max_p, p, total_p ;

  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_inputs, NULL) ;
    if (!mri_dst)
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst") ;
    MRIcopyHeader(mri_inputs, mri_dst) ;
  }

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_inputs->width ; height = mri_inputs->height; 
  depth = mri_inputs->depth ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == 152 && y == 126 && z == 127)
          DiagBreak() ;
        if (x == 63 && y == 107 && z == 120)
          DiagBreak() ; 

        val = MRIvox(mri_inputs, x, y, z) ;
        GCAsourceVoxelToNode(gca, mri_inputs, lta, x, y, z, &xn, &yn, &zn) ;
        
        gcan = &gca->nodes[xn][yn][zn] ;
        label = 0 ; max_p = BIG_AND_NEGATIVE ;
        for (total_p = 0.0, n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          
          /* compute 1-d Mahalanobis distance */
          dist = (val-gc->mean) ;
          if (FZERO(gc->var))  /* make it a delta function */
          {
            if (FZERO(dist))
              p = 1.0 ;
            else
              p = 0.0 ;
          }
          else
            p = 1 / sqrt(gc->var * 2 * M_PI) * exp(-dist*dist/gc->var) ;

          p *= gc->prior ;
          if (p > max_p)
          {
            max_p = p ;
            label = gcan->labels[n] ;
          }
          total_p += p ;
        }
        max_p = 255.0* max_p / total_p ; if (max_p > 255) max_p = 255 ;
        MRIvox(mri_dst, x, y, z) = (BUFTYPE)max_p ;
      }
    }
  }

  return(mri_dst) ;
}
MRI  *
GCAcomputeProbabilities(MRI *mri_inputs, GCA *gca, MRI *mri_labels, 
                        MRI *mri_dst, LTA *lta)
{
  int      x, y, z, width, height, depth, label, val,
           xn, yn, zn, n ;
  GCA_NODE *gcan ;
  GC1D     *gc ;
  double   dist, label_p, p, total_p ;

  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_inputs, NULL) ;
    if (!mri_dst)
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst") ;
    MRIcopyHeader(mri_inputs, mri_dst) ;
  }


  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_inputs->width ; height = mri_inputs->height; 
  depth = mri_inputs->depth ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == 152 && y == 126 && z == 127)
          DiagBreak() ;
        if (x == 63 && y == 107 && z == 120)
          DiagBreak() ; 

        GCAsourceVoxelToNode(gca, mri_inputs, lta, x, y, z, &xn, &yn, &zn) ;
        val = MRIvox(mri_inputs, x, y, z) ;
        label = MRIvox(mri_labels, x, y, z) ;

        gcan = &gca->nodes[xn][yn][zn] ;
        for (label_p = total_p = 0.0, n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          
          /* compute 1-d Mahalanobis distance */
          dist = (val-gc->mean) ;
          if (FZERO(gc->var))  /* make it a delta function */
          {
            if (FZERO(dist))
              p = 1.0 ;
            else
              p = 0.0 ;
          }
          else
            p = 1 / sqrt(gc->var * 2 * M_PI) * exp(-dist*dist/gc->var) ;

          p *= gc->prior ;
          if (label == gcan->labels[n])
          {
            label_p = p ;
          }
          total_p += p ;
        }
        label_p = 255.0* label_p / total_p ; if (label_p > 255) label_p = 255 ;
        MRIvox(mri_dst, x, y, z) = (BUFTYPE)label_p ;
      }
    }
  }

  return(mri_dst) ;
}

#define STARTING_T   500

MRI  *
GCAannealUnlikelyVoxels(MRI *mri_inputs, GCA *gca, MRI *mri_dst,LTA *lta, 
          int max_iter, MRI  *mri_fixed)
{
  int    x, y, z, width, depth, height, *x_indices, *y_indices, *z_indices,
         nindices, index, iter, nchanged, xn, yn, zn, n, nbad,
         old_label  ;
  double log_likelihood, T, delta_E, p, rn, new_ll, old_ll,
         total_likelihood ;
  GCA_NODE *gcan ;
  MRI      *mri_bad ;

  mri_bad = MRIclone(mri_inputs, NULL) ;

  width = mri_dst->width ; height = mri_dst->height ; 
  depth = mri_dst->depth ;
  

  for (nindices = x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (mri_fixed && MRIvox(mri_fixed, x, y, z) > 0)
          continue ;
        log_likelihood = 
          gcaVoxelGibbsLogLikelihood(gca, mri_dst, mri_inputs, x, y, z,lta,
                                     PRIOR_FACTOR);
        if (log_likelihood < BIG_AND_NEGATIVE/2)
        {
          MRIvox(mri_bad, x, y, z) = 1 ;
          nindices++ ;
        }
      }
    }
  }

  printf("annealing %d voxels...\n", nindices) ;
  x_indices = (int *)calloc(nindices, sizeof(int)) ;
  y_indices = (int *)calloc(nindices, sizeof(int)) ;
  z_indices = (int *)calloc(nindices, sizeof(int)) ;
  for (index = x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (MRIvox(mri_bad, x, y, z) > 0)
        {
          x_indices[index] = x ; y_indices[index] = y ; z_indices[index] = z ;
          index++ ;
        }
      }
    }
  }

  MRIfree(&mri_bad) ;
  T = STARTING_T ; iter = 0 ;
  do
  {
    total_likelihood = 0.0 ;
    for (nbad = nchanged = index = 0 ; index < nindices ; index++)
    {
      x = x_indices[index] ; y = y_indices[index] ; z = z_indices[index] ;
      if (x == 155 && y == 126 && z == 128)
        DiagBreak() ;
      
      /* find the node associated with this coordinate and classify */
      GCAsourceVoxelToNode(gca, mri_inputs, lta, x, y, z, &xn, &yn, &zn) ;
      gcan = &gca->nodes[xn][yn][zn] ;
      
      if (gcan->nlabels == 1)
        continue ;
      n = (int)randomNumber(0.0, (double)gcan->nlabels-0.0001) ;
      if (gcan->labels[n] == MRIvox(mri_dst, x, y, z))
        continue ;
      old_ll = 
        gcaNbhdGibbsLogLikelihood(gca, mri_dst, mri_inputs, x, y, z, lta,
                                  PRIOR_FACTOR) ;
      old_label = MRIvox(mri_dst, x, y, z) ;
      MRIvox(mri_dst, x, y, z) = gcan->labels[n] ;
      new_ll = 
        gcaNbhdGibbsLogLikelihood(gca, mri_dst, mri_inputs, x, y, z, lta,
                                  PRIOR_FACTOR) ;
      delta_E = new_ll - old_ll ;
      p = exp(delta_E / T) ;
      rn = randomNumber(0.0, 1.0) ;
      
      if (p > rn)
      {
        if (new_ll < BIG_AND_NEGATIVE/2)
          nbad++ ;
        nchanged++ ;
        total_likelihood += new_ll ;
      }
      else
      {
        total_likelihood += old_ll ;
        if (old_ll < BIG_AND_NEGATIVE/2)
          nbad++ ;
        MRIvox(mri_dst, x, y, z) = old_label ;
      }
    }
    T = T * 0.99 ;
    fprintf(stderr, "%03d: T = %2.2f, nchanged %d, nbad = %d, ll=%2.2f\n",
            iter, T, nchanged, nbad, total_likelihood/(double)nindices) ;
    if (!nchanged)
      break ;
  } while (iter++ < max_iter) ;

  free(x_indices) ; free(y_indices) ; free(z_indices) ;
  return(mri_dst) ;
}

GCA  *
GCAreduce(GCA *gca_src)
{
  GCA       *gca_dst ;
  int       xs, ys, zs, xd, yd, zd, swidth, sheight, sdepth,
            ns, nd, dwidth, dheight, ddepth, dof ;
  float     spacing = gca_src->spacing ;
  GCA_NODE  *gcan_src, *gcan_dst ;
  GC1D      *gc_src, *gc_dst ;

  swidth = gca_src->width ; sheight = gca_src->height; sdepth = gca_src->depth;
  
  gca_dst = 
    GCAalloc(gca_src->ninputs, spacing*2, spacing*swidth,
             spacing*gca_src->height,spacing*sdepth) ;


  dwidth = gca_dst->width ; dheight = gca_dst->height; ddepth = gca_dst->depth;

  for (zs = 0 ; zs < sdepth ; zs++)
  {
    for (ys = 0 ; ys < sheight ; ys++)
    {
      for (xs = 0 ; xs < swidth ; xs++)
      {
        xd = xs/2 ; yd = ys/2 ; zd = zs/2 ;
        if (xd == 15 && yd == 22 && zd == 35)
          DiagBreak() ;
        gcan_src = &gca_src->nodes[xs][ys][zs] ;
        gcan_dst = &gca_dst->nodes[xd][yd][zd] ;
        
        for (ns = 0 ; ns < gcan_src->nlabels ; ns++)
        {
          gc_src = &gcan_src->gcs[ns] ;
          dof = gc_src->prior * gcan_src->total_training ;
          
          /* find label in destination */
          for (nd = 0 ; nd < gcan_dst->nlabels ; nd++)
          {
            if (gcan_dst->labels[nd] == gcan_src->labels[ns])
              break ;
          }
          if (nd >= gcan_dst->nlabels)  /* couldn't find it */
          {
            if (nd >= gcan_dst->max_labels)
            {
              int  old_max_labels ;
              char *old_labels ;
              GC1D *old_gcs ;
              
              old_max_labels = gcan_dst->max_labels ; 
              gcan_dst->max_labels += 2 ;
              old_labels = gcan_dst->labels ;
              old_gcs = gcan_dst->gcs ;
              
              /* allocate new ones */
              gcan_dst->gcs = 
                (GC1D *)calloc(gcan_dst->max_labels, sizeof(GC1D)) ;
              if (!gcan_dst->gcs)
                ErrorExit(ERROR_NOMEMORY, 
                          "GCANreduce: couldn't expand gcs to %d",
                          gcan_dst->max_labels) ;
              gcan_dst->labels = 
                (char *)calloc(gcan_dst->max_labels, sizeof(char)) ;
              if (!gcan_dst->labels)
                ErrorExit(ERROR_NOMEMORY, 
                          "GCANupdateNode: couldn't expand labels to %d",
                          gcan_dst->max_labels) ;
              
              /* copy the old ones over */
              memmove(gcan_dst->gcs, old_gcs, old_max_labels*sizeof(GC1D));
              memmove(gcan_dst->labels, old_labels, 
                      old_max_labels*sizeof(char)) ;
              
              /* free the old ones */
              free(old_gcs) ; free(old_labels) ;
            }
            gcan_dst->nlabels++ ;
          }
          gc_dst = &gcan_dst->gcs[nd] ;
          gc_dst->mean += dof*gc_src->mean ;
          gc_dst->var += dof*(gc_src->var) ;
          gc_dst->prior += dof ;
          gcan_dst->total_training += dof ;
          gcan_dst->labels[nd] = gcan_src->labels[ns] ;
        }
      }
    }
  }

  for (xd = 0 ; xd < gca_dst->width ; xd++)
  {
    for (yd = 0 ; yd < gca_dst->height ; yd++)
    {
      for (zd = 0 ; zd < gca_dst->depth ; zd++)
      {
        gcan_dst = &gca_dst->nodes[xd][yd][zd] ;
        for (nd = 0 ; nd < gcan_dst->nlabels ; nd++)
        {
          float var ;

          gc_dst = &gcan_dst->gcs[nd] ;
          dof = gc_dst->prior ;/* prior is count of # of occurences now */
          gc_dst->mean /= dof ;
          var = gc_dst->var / dof ;
          if (var < -0.1)
            DiagBreak() ;
          if (var < MIN_VAR)
            var = MIN_VAR ;
          gc_dst->var = var ;
          gc_dst->prior /= (float)gcan_dst->total_training ;
        }
      }
    }
  }
  return(gca_dst) ;
}

#define MAX_LABELS_PER_GCAN              5
#define MAX_DIFFERENT_LABELS    500

typedef struct
{
  int   label ;
  float prob ;
} LABEL_PROB ;

static int compare_sort_probabilities(const void *plp1, const void *plp2);
MRI *
GCAclassify(MRI *mri_inputs,GCA *gca,MRI *mri_dst,LTA *lta,int max_labels)
{
  int        x, y, z, width, height, depth, val,
             xn, yn, zn, n ;
  GCA_NODE   *gcan ;
  GC1D       *gc ;
  float      dist, max_p, p, total_p ;
  LABEL_PROB label_probs[1000] ;

  if (max_labels > MAX_LABELS_PER_GCAN || max_labels <= 0)
    max_labels = MAX_LABELS_PER_GCAN ;

  if (!mri_dst)
  {
    int width = mri_inputs->width, height = mri_inputs->height, 
        depth = mri_inputs->depth ;

    mri_dst = MRIallocSequence(width, height, depth, MRI_UCHAR, 2*max_labels) ;
    if (!mri_dst)
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst") ;
    MRIcopyHeader(mri_inputs, mri_dst) ;
  }

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_inputs->width ; height = mri_inputs->height; 
  depth = mri_inputs->depth ;
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == 67 && y == 87 && z == 114)
          DiagBreak() ; 

        val = MRIvox(mri_inputs, x, y, z) ;

        /* find the node associated with this coordinate and classify */
        GCAsourceVoxelToNode(gca, mri_inputs, lta, x, y, z, &xn, &yn, &zn) ;
        gcan = &gca->nodes[xn][yn][zn] ;
        max_p = BIG_AND_NEGATIVE ;
        for (total_p = 0.0, n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          
          /* compute 1-d Mahalanobis distance */
          dist = (val-gc->mean) ;
          if (FZERO(gc->var))  /* make it a delta function */
          {
            if (FZERO(dist))
              p = 1.0 ;
            else
              p = 0.0 ;
          }
          else
            p = 1 / sqrt(gc->var * 2 * M_PI) * exp(-dist*dist/gc->var) ;

          p *= gc->prior ;
          total_p += p ;
          label_probs[n].prob = p ; label_probs[n].label = gcan->labels[n] ;
        }
        /* now sort the labels and probabilities */
        qsort(label_probs, gcan->nlabels, sizeof(LABEL_PROB), 
              compare_sort_probabilities) ;

        for (n = 0 ; n < max_labels ; n++)
        {
          if (n < gcan->nlabels)
          {
            MRIseq_vox(mri_dst, x, y, z, n*2) = label_probs[n].label ;
            MRIseq_vox(mri_dst, x, y, z, n*2+1) = 
              (BUFTYPE)nint(255.0*label_probs[n].prob/total_p) ;
          }
          else
          {
            MRIseq_vox(mri_dst, x, y, z, n*2) = 255 ;
            MRIseq_vox(mri_dst, x, y, z, n*2+1) = 0 ;
          }
        }
      }
    }
  }

  return(mri_dst) ;
}
 
static int
compare_sort_probabilities(const void *plp1, const void *plp2)
{
  LABEL_PROB  *lp1, *lp2 ;

  lp1 = (LABEL_PROB *)plp1 ;
  lp2 = (LABEL_PROB *)plp2 ;

  if (lp1->prob > lp2->prob)
    return(1) ;
  else if (lp1->prob < lp2->prob)
    return(-1) ;

  return(0) ;
}

float
GCAcomputeLogSampleProbability(GCA *gca, GCA_SAMPLE *gcas, 
                               MRI *mri_inputs, MATRIX *m_L, int nsamples)
{
  int        x, y, z, xt, yt, zt, width, height, depth, val,
             xn, yn, zn, i ;
  MATRIX     *m_L_inv ;
  VECTOR     *v_input, *v_canon ;
  GCA_NODE   *gcan ;
  GC1D       *gc ;
  float      dist ;
  double     total_log_p, log_p ;

  m_L_inv = MatrixInverse(m_L, NULL) ;
  if (!m_L_inv)
  {
    MatrixPrint(stderr, m_L) ;
    ErrorExit(ERROR_BADPARM, 
              "GCAcomputeLogSampleProbability: matrix is not invertible") ;
  }
  v_input = VectorAlloc(4, MATRIX_REAL) ;
  v_canon = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v_input, 4, 1) = 1.0 ; *MATRIX_RELT(v_canon, 4, 1) = 1.0 ;

  /* go through each GC in the sample and compute the probability of
     the image at that point.
  */
  width = mri_inputs->width ; height = mri_inputs->height; 
  depth = mri_inputs->depth ;
  for (total_log_p = 0.0, i = 0 ; i < nsamples ; i++)
  {
    if (i == Gdiag_no)
      DiagBreak() ;
    if (Gdiag_no == gcas[i].label)
      DiagBreak() ;

    xn = gcas[i].xn ; yn = gcas[i].yn ; zn = gcas[i].zn ; 
    GCAnodeToVoxel(gca, mri_inputs, xn, yn, zn, &xt, &yt, &zt) ;
    V3_X(v_canon) = (float)xt ;
    V3_Y(v_canon) = (float)yt ;
    V3_Z(v_canon) = (float)zt ;
    MatrixMultiply(m_L_inv, v_canon, v_input) ;
    x = nint(V3_X(v_input)); y = nint(V3_Y(v_input)); z = nint(V3_Z(v_input));
    if (x < 0) x = 0 ;
    if (y < 0) y = 0 ;
    if (z < 0) z = 0 ;
    if (x >= width)  x = width-1 ;
    if (y >= height) y = height-1 ;
    if (z >= depth)  z = depth-1 ;
    val = MRIvox(mri_inputs, x, y, z) ;

    gcan = &gca->nodes[xn][yn][zn] ;
    gc = &gcan->gcs[gcas[i].n] ;
#define TRIM_DISTANCES 0
#if TRIM_DISTANCES
         
    dist = (val-gc->mean) ;
#define TRIM_DIST 20
    if (abs(dist) > TRIM_DIST)
      dist = TRIM_DIST ;
    log_p =
      -log(sqrt(gc->var)) - 
      0.5 * (dist*dist/gc->var) +
      log(gc->prior) ;
#else
    dist = (val-gcas[i].mean) ;
    log_p =
      -log(sqrt(gcas[i].var)) - 
      0.5 * (dist*dist/gcas[i].var) +
      log(gcas[i].prior) ;
#endif
    total_log_p += log_p ;
    gcas[i].log_p = log_p ;

    if (!finite(total_log_p))
    {
      DiagBreak() ;
      fprintf(stderr, 
              "total log p not finite at (%d, %d, %d) n = %d, "
              "var=%2.2f\n", x, y, z, gcas[i].n, gc->var) ;
    }
  }

  MatrixFree(&m_L_inv) ; MatrixFree(&v_canon) ; MatrixFree(&v_input) ;
  return((float)total_log_p) ;
}
/*
  compute the probability of the image given the transform and the class
  stats.
*/
float
GCAcomputeLogImageProbability(GCA *gca, MRI *mri_inputs, MRI *mri_labels,
                              LTA *lta)
{
  int        x, y, z, width, height, depth, val,
             xn, yn, zn, n, label ;
  GCA_NODE   *gcan ;
  GC1D       *gc ;
  float      dist ;
  double     total_log_p ;


  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_inputs->width ; height = mri_inputs->height; 
  depth = mri_inputs->depth ;
  for (total_log_p = 0.0, x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == 85 && y == 89 && z == 135)
          DiagBreak() ; 

        val = MRIvox(mri_inputs, x, y, z) ;
        label = MRIvox(mri_labels, x, y, z) ;

        /* find the node associated with this coordinate and classify */
        GCAsourceVoxelToNode(gca, mri_inputs, lta, x, y, z, &xn, &yn, &zn) ;
        gcan = &gca->nodes[xn][yn][zn] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          if (gcan->labels[n] == label)
            break ;
        }
        if (n < gcan->nlabels)
        {
          gc = &gcan->gcs[n] ;
         
          dist = (val-gc->mean) ;
          total_log_p +=
            -log(sqrt(gc->var)) - 
            0.5 * (dist*dist/gc->var) +
            log(gc->prior) ;
          if (!finite(total_log_p))
          {
            DiagBreak() ;
            fprintf(stderr, 
                    "total log p not finite at (%d, %d, %d) n = %d, "
                    "var=%2.2f\n", x, y, z, n, gc->var) ;
          }
        }
      }
    }
  }

  return((float)total_log_p) ;
}

#if 0
int
GCAsampleStats(GCA *gca, MRI *mri, int class, 
               Real x, Real y, Real z,
               Real *pmean, Real *pvar, Real *pprior)
{
  int    xm, xp, ym, yp, zm, zp, width, height, depth, n ;
  Real   val, xmd, ymd, zmd, xpd, ypd, zpd ;  /* d's are distances */
  Real   prior, mean, var, wt ;
  GCAN   *gcan ;
  GC1D   *gc ;

  width = gca->width ; height = gca->height ; depth = gca->depth ;

  xm = MAX((int)x, 0) ;
  xp = MIN(width-1, xm+1) ;
  ym = MAX((int)y, 0) ;
  yp = MIN(height-1, ym+1) ;
  zm = MAX((int)z, 0) ;
  zp = MIN(depth-1, zm+1) ;

  xmd = x - (float)xm ;
  ymd = y - (float)ym ;
  zmd = z - (float)zm ;
  xpd = (1.0f - xmd) ;
  ypd = (1.0f - ymd) ;
  zpd = (1.0f - zmd) ;

  prior = mean = var = 0.0 ;


  gcan = &gca->nodes[xm][ym][zm] ;
  wt = xpd * ypd * zpd 
  for (n = 0 ; n < gcan->nlabels ; n++)
  {
    if (gcan->labels[n] == class)
      break ;
  }
  if (n < gcan->nlabels)   /* found it */
  {
    gc = &gcan->gcs[n] ;
    prior += wt * gc->prior ;
  }

  gcan = &gca->nodes[xp][ym][zm] ;
  wt = xmd * ypd * zpd 

  gcan = &gca->nodes[xm][yp][zm] ;
  wt = xpd * ymd * zpd ;

  gcan = &gca->nodes[xp][yp][zm] ;
  wt = xmd * ymd * zpd ;

  gcan = &gca->nodes[xm][ym][zp] ;
  wt = xpd * ypd * zmd ;

  gcan = &gca->nodes[xp][ym][zp] ;
  wt = xmd * ypd * zmd ;

  gcan = &gca->nodes[xm][yp][zp] ;
  wt = xpd * ymd * zmd ;

  gcan = &gca->nodes[xp][yp][zp] ;
  wt = xmd * ymd * zmd ;
  return(NO_ERROR) ;
}
#endif  

#define MIN_SPACING   8

int
GCAtransformSamples(GCA *gca_src, GCA *gca_dst, GCA_SAMPLE *gcas, int nsamples)
{
  int       scale, i, label, xd, yd, zd, xs, ys, zs, n, xk, yk, zk, 
            xd0, yd0, zd0, min_y_i ;
  float     max_p, min_v, vscale, min_y ;
  GCA_NODE  *gcan ;
  GC1D      *gc ;
  MRI       *mri_found ;

  vscale = 1 ;
  mri_found = MRIalloc(gca_dst->width*vscale, gca_dst->height*vscale, 
                       gca_dst->depth*vscale, MRI_UCHAR) ;

  scale = gca_src->spacing / gca_dst->spacing ;
  min_y = 10000 ; min_y_i = -1 ;
  for (i = 0 ; i < nsamples ; i++)
  {
    label = gcas[i].label ; 
    xs = gcas[i].xn ; ys = gcas[i].yn ; zs = gcas[i].zn ; 
    xd0 = (int)(xs * scale) ; yd0 = (int)(ys * scale) ; 
    zd0 = (int)(zs * scale) ;
    max_p = -1.0 ;  /* find appropriate label with highest prior in dst */
    min_v =  10000000.0f ;
    for (xk = -scale/2 ; xk <= scale/2 ; xk++)
    {
      xd = MIN(MAX(0, xd0+xk), gca_dst->width-1) ;
      for (yk = -scale/2 ; yk <= scale/2 ; yk++)
      {
        yd = MIN(MAX(0, yd0+yk), gca_dst->height-1) ;
        for (zk = -scale/2 ; zk <= scale/2 ; zk++)
        {
          zd = MIN(MAX(0, zd0+zk), gca_dst->height-1) ;
          if (MRIvox(mri_found, (int)(xd*vscale), (int)(yd*vscale), 
                     (int)(zd*vscale)))
            continue ;
          gcan = &gca_dst->nodes[xd][yd][zd] ;
          for (n = 0 ; n < gcan->nlabels ; n++)
          {
            gc = &gcan->gcs[n] ;
            if (gcan->labels[n] == label && 
                (gc->prior > max_p || 
                 (FEQUAL(gc->prior,max_p) && gc->var < min_v)))
            {
              min_v = gc->var ;
              max_p = gc->prior ;
              gcas[i].xn = xd ; gcas[i].yn = yd ; gcas[i].zn = zd ; 
              gcas[i].n = n ;
            }
          }
        }
      }
    }
    if (max_p < 0)
    {
      fprintf(stderr, "WARNING: label %d not found at (%d,%d,%d)\n",
              label, gcas[i].xn, gcas[i].yn, gcas[i].zn) ;
      DiagBreak() ;
    }
    MRIvox(mri_found, (int)(gcas[i].xn*vscale),
           (int)(gcas[i].yn*vscale), (int)(gcas[i].zn*vscale))= 1;
    if (gcas[i].yn < min_y)
    {
      min_y = gcas[i].yn ;
      min_y_i = i ;
    }
  }

  i = min_y_i ;
  fprintf(stderr, "min_y = (%d, %d, %d) at i=%d, label=%d, n=%d\n",
          gcas[i].xn, gcas[i].yn, gcas[i].zn, i, gcas[i].label, gcas[i].n) ;
  MRIfree(&mri_found) ;
  return(NO_ERROR) ;
}

/* don't use a label with prior less than this */
#define MIN_MAX_PRIORS 0.5


static int exclude_classes[] = 
{
  0 /*1, 6, 21, 22, 23, 24, 25, 30, 57, 61, 62, 63*/
} ;


#define TILES     3
#define MAX_PCT   .1

static int compare_gca_samples(const void *pc1, const void *pc2) ;

static int
compare_gca_samples(const void *pgcas1, const void *pgcas2)
{
  register GCA_SAMPLE *gcas1, *gcas2 ;

  gcas1 = (GCA_SAMPLE *)pgcas1 ;
  gcas2 = (GCA_SAMPLE *)pgcas2 ;

/*  return(c1 > c2 ? 1 : c1 == c2 ? 0 : -1) ;*/
  if (FEQUAL(gcas1->prior, gcas2->prior))
  {
#if 0
    if (FEQUAL(gcas1->var, gcas2->var))
    {
      int  zv1, zv2 ;

      zv1 = gcas1->zn%10 ; zv2 = gcas2->zn%10 ;
      if (zv1 == zv2)
      {
        int  yv1, yv2 ;

        yv1 = gcas1->yn%10 ; zv2 = gcas2->yn%10 ;      
        return(yv1-yv2) ;
      }

      return(zv1-zv2) ;
    }
#endif

    if (gcas1->var > gcas2->var)
      return(1) ;
    else
      return(-1) ;
  }

  if (gcas1->prior > gcas2->prior)
    return(-1) ;
  else if (gcas1->prior < gcas2->prior)
    return(1) ;

  return(0) ;
}
#define MAX_SPACING   16  /* mm */

GCA_SAMPLE *
GCAfindStableSamplesByLabel(GCA *gca, int nsamples, float min_prior)
{
  GCA_SAMPLE *gcas, *ordered_labels[MAX_DIFFERENT_LABELS], *gcas2 ;
  GCA_NODE   *gcan ;
  GC1D       *gc ;
  int        found[MAX_DIFFERENT_LABELS], x, y, z, width, height, depth, n,
             label, nfound, samples_added[MAX_DIFFERENT_LABELS] ;
  float      histo[MAX_DIFFERENT_LABELS], spacing, scale ;
  int        volume[TILES][TILES], max_class, total, extra, total_found,
             label_counts[MAX_DIFFERENT_LABELS], 
             current_index[MAX_DIFFERENT_LABELS],
             ordered_label_counts[MAX_DIFFERENT_LABELS], i, index,
             *x_indices, *y_indices, *z_indices, nindices ;

  memset(histo, 0, sizeof(histo)) ;
  memset(label_counts, 0, sizeof(label_counts)) ;
  memset(samples_added, 0, sizeof(samples_added)) ;
  memset(current_index, 0, sizeof(current_index)) ;
  memset(ordered_label_counts, 0, sizeof(ordered_label_counts)) ;
  memset(found, 0, sizeof(found)) ;
  memset(volume, 0, sizeof(volume)) ;
  gcas = calloc(nsamples, sizeof(GCA_SAMPLE )) ;
  width = gca->width ; height = gca->height ; depth = gca->depth ;

  /* compute the max priors and min variances for each class */
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        gcan = &gca->nodes[x][y][z] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          label = gcan->labels[n] ;
          if (label == Gdiag_no)
            DiagBreak() ;
          histo[label] += (float)nint(gcan->total_training*gc->prior) ;
          label_counts[label]++ ;
        }
      }
    }
  }

  for (max_class = MAX_DIFFERENT_LABELS ; max_class >= 0 ; max_class--)
    if (histo[max_class] > 0)
      break ;
  fprintf(stderr, "max class = %d\n", max_class) ;

  /* count total # of samples */
  for (total = 0, n = 1 ; n <= max_class ; n++)
    total += histo[n] ;

  fprintf(stderr, "%d total training samples found\n", total) ;

  /* turn histogram into samples/class */
  for (n = 1 ; n <= max_class ; n++)
    histo[n] = (float)nint(0.25+histo[n]*(float)nsamples/total) ;

  for (n = 0 ; n < sizeof(exclude_classes)/sizeof(exclude_classes[0]) ; n++)
    histo[exclude_classes[n]] = 0 ;

  /* crop max # per class */
  for (n = 1 ; n <= max_class ; n++)
    if (histo[n] > nint(MAX_PCT*nsamples))
      histo[n] = nint(MAX_PCT*nsamples) ;
      
  for (extra = 0, n = 1 ; n <= max_class ; n++)
    extra += histo[n] ;
  extra = nsamples-extra ;  /* from rounding */
  printf("%d extra samples for redistribution...\n", extra) ;

  /* first add to classes with only one sample */
  for (n = 1 ; extra > 0 && n <= max_class ; n++)
  {
    if (histo[n] == 1)
    {
      histo[n]++ ;
      extra-- ;
    }
  }  

  while (extra > 0)  /* add 1 to each class */
  {
    for (n = 1 ; extra > 0 && n <= max_class ; n++)
    {
      if (histo[n] >= 1)
      {
        histo[n]++ ;
        extra-- ;
      }
    }
  }

  {
    FILE *fp ;
    fp = fopen("classes.dat", "w") ;
    for (n = 1 ; n <= max_class ; n++)
      if (!FZERO(histo[n]))
        fprintf(fp, "%d  %2.1f\n", n, histo[n]) ;
    fclose(fp) ;
  }

  /* allocate arrays that will be used for sorting */
  for (n = 1 ; n <= max_class ; n++)
    if (histo[n] > 0)
    {
      ordered_labels[n] = 
        (GCA_SAMPLE *)calloc(label_counts[n], sizeof(GCA_SAMPLE)) ;
      if (!ordered_labels[n])
        ErrorExit(ERROR_NO_MEMORY, 
                  "GCAfindStableSamplesByLabel: could not allocate %d samples "
                  "for label %d", label_counts[n], n) ;
    }


  nindices = width*height*depth ;
  x_indices = (int *)calloc(nindices, sizeof(int)) ;
  y_indices = (int *)calloc(nindices, sizeof(int)) ;
  z_indices = (int *)calloc(nindices, sizeof(int)) ;

  for (i = 0 ; i < nindices ; i++)
  {
    x_indices[i] = i % width ;
    y_indices[i] = (i/width) % height ;
    z_indices[i] = (i / (width*height)) % depth ;
  }
  for (i = 0 ; i < nindices ; i++)
  {
    int tmp ;
    index = (int)randomNumber(0.0, (double)(nindices-0.0001)) ;

    tmp = x_indices[index] ; 
    x_indices[index] = x_indices[i] ; x_indices[i] = tmp ;

    tmp = y_indices[index] ; 
    y_indices[index] = y_indices[i] ; y_indices[i] = tmp ;

    tmp = z_indices[index] ; 
    z_indices[index] = z_indices[i] ; z_indices[i] = tmp ;
  }

  for (index = 0 ; index < nindices ; index++)
  {
    x = x_indices[index] ; y = y_indices[index] ; z = z_indices[index] ;
    gcan = &gca->nodes[x][y][z] ;
    for (n = 0 ; n < gcan->nlabels ; n++)
    {
      gc = &gcan->gcs[n] ;
      label = gcan->labels[n] ;
      if (histo[label] > 0)
      {
        i = ordered_label_counts[label] ;
        ordered_label_counts[label]++ ;
        ordered_labels[label][i].xn = x ;
        ordered_labels[label][i].yn = y ;
        ordered_labels[label][i].zn = z ;
        ordered_labels[label][i].n = n ;
        ordered_labels[label][i].label = label ;
        ordered_labels[label][i].prior = gc->prior ;
        ordered_labels[label][i].var = gc->var ;
        ordered_labels[label][i].mean = gc->mean ;
      }
    }
  }


  free(x_indices) ; free(y_indices) ; free(z_indices) ;

  /* sort each list of labels by prior, then by variance */
  for (n = 1 ; n <= max_class ; n++)
    if (histo[n] > 0)
    {
      qsort(ordered_labels[n], ordered_label_counts[n], sizeof(GCA_SAMPLE), 
            compare_gca_samples) ;
    }

  total_found = nfound = 0 ;

  for (spacing = MAX_SPACING ; spacing >= gca->spacing ; spacing /= 2)
  {
    MRI  *mri_found ;
    int  xv, yv, zv ;

    for (n = 1 ; n <= max_class ; n++) 
      current_index[n] = 0 ;

    scale = gca->spacing / spacing ;
    mri_found = MRIalloc(width*scale, height*scale, depth*scale, MRI_UCHAR);
    for (i = 0 ; i < total_found ; i++)
    {
      xv = gcas[i].xn*scale ; yv = gcas[i].yn*scale ; zv = gcas[i].zn*scale ;
      MRIvox(mri_found, xv, yv, zv) = 1 ;
    }

    do
    {

      nfound = 0 ;
      for (n = 1 ; n <= max_class ; n++) 
      {
        if (samples_added[n] < histo[n])/* add another sample from this class*/
        {
          while ((ordered_labels[n][current_index[n]].label < 0) &&
                 current_index[n] < ordered_label_counts[n])
            current_index[n]++ ;
          if (current_index[n] < ordered_label_counts[n])
          {
            gcas2 = &ordered_labels[n][current_index[n]] ;
            current_index[n]++ ; 
            xv = gcas2->xn*scale ; yv = gcas2->yn*scale ; zv = gcas2->zn*scale;
            
            if (n == Gdiag_no)
              DiagBreak() ;
            
            if (!MRIvox(mri_found, xv, yv, zv)) /* none at this location yet */
            {
              if (gcas2->label == Gdiag_no)
                DiagBreak() ;
              samples_added[n]++ ;
              gcas[total_found+nfound++] = *gcas2 ;
              MRIvox(mri_found, xv, yv, zv) = gcas2->label ;
              gcas2->label = -1 ;
              if (nfound+total_found >= nsamples)
                break ;
            }
          }
        }
      }

      total_found += nfound ;
    } while ((nfound > 0) && total_found < nsamples) ;
    MRIfree(&mri_found) ;
  } 

  if (total_found != nsamples)
  {
    ErrorPrintf(ERROR_BADPARM, "could only find %d samples!\n", total_found) ;
  }
  return(gcas) ;
}


GCA_SAMPLE *
GCAfindContrastSamples(GCA *gca, int *pnsamples, int min_spacing,
                       float min_prior)
{
  GCA_SAMPLE *gcas ;
  GC1D       *gc ;
  int        x, y, z, width, height, depth, label, nfound, node_stride,
             label_counts[MAX_DIFFERENT_LABELS], best_label, nzeros,
             xi, yi, zi, xk, yk, zk, i, j, k, labels[3][3][3], found ;
  float      max_prior, total_mean,
             priors[3][3][3][MAX_DIFFERENT_LABELS], best_mean, best_sigma,
             means[3][3][3][MAX_DIFFERENT_LABELS], mean, sigma,
             vars[3][3][3][MAX_DIFFERENT_LABELS] ;

  memset(label_counts, 0, sizeof(label_counts)) ;
  width = gca->width ; height = gca->height ; depth = gca->depth ;
  gcas = calloc(width*height*depth, sizeof(GCA_SAMPLE)) ;
  
  node_stride = min_spacing / gca->spacing ;

  total_mean = 0.0 ;

  /* compute the max priors and min variances for each class */
  for (nzeros = nfound = x = 0 ; x < width ; x += node_stride)
  {
    for (y = 0 ; y < height ; y += node_stride)
    {
      for (z = 0 ; z < depth ; z += node_stride)
      {
        if (abs(x-31)<=node_stride &&
            abs(y-22)<=node_stride &&
            abs(z-36)<=node_stride)
          DiagBreak() ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xi = x+xk ; i = xk+1 ;
          if (xi < 0 || xi >= width)
            continue ;
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yi = y+yk ; j = yk + 1 ;
            if (yi < 0 || yi >= height)
              continue ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zi = z+zk ; k = zk+1 ;
              if (zi < 0 || zi >= depth)
                continue ;
              gcaRegionStats(gca, xi, yi, zi, 
                             node_stride/3, 
                             priors[i][j][k], 
                             vars[i][j][k], 
                             means[i][j][k]) ;
              if (priors[i][j][k][4] > .5*min_prior ||  /* left lat ven */
                  priors[i][j][k][5] > .5*min_prior ||  /* inf left lat ven */
                  priors[i][j][k][14] > .5*min_prior ||  /* 3rd ven */
                  priors[i][j][k][15] > .5*min_prior ||  /* 4th ven */
                  priors[i][j][k][43] > .5*min_prior ||
                  priors[i][j][k][44] > .5*min_prior)
                DiagBreak() ;
              max_prior = 0 ;

              /* find highest prior label */
              for (label = 0 ; label < MAX_DIFFERENT_LABELS ; label++)
              {
                if (priors[i][j][k][label] > max_prior)
                {
                  max_prior = priors[i][j][k][label] ;
                  labels[i][j][k] = label ;
                }
              }
            }
          }
        }

        /* search nbrhd for high-contrast stable label */
        best_label = labels[1][1][1] ; found = 0 ;
        best_mean = means[1][1][1][best_label] ; 
        best_sigma = sqrt(vars[1][1][1][best_label]) ;

        gc = gcaFindHighestPriorGC(gca, x, y, z, best_label,node_stride/3);
        if (!gc || gc->prior < min_prior)
          continue ;

        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xi = x+xk ; i = xk+1 ;
          if (xi < 0 || xi >= width)
            continue ;
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yi = y+yk ; j = yk + 1 ;
            if (yi < 0 || yi >= height)
              continue ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zi = z+zk ; k = zk+1 ;
              if (zi < 0 || zi >= depth)
                continue ;
              if (i == 1 && j == 1 && k == 1)
                continue ;
              label = labels[i][j][k] ;
              if (priors[i][j][k][label] < min_prior || label == best_label)
                continue ;
              
              mean = means[i][j][k][label] ; 
              sigma = sqrt(vars[i][j][k][label]) ;
              gc = gcaFindHighestPriorGC(gca, xi, yi, zi, label,node_stride/3);
              if (!gc || gc->prior < min_prior)
                continue ;
              if (abs(best_mean - mean) > 4*abs(best_sigma+sigma))
              {
#if 0
                if (gcaFindBestSample(gca, xi, yi, zi, label, node_stride/3, 
                         &gcas[nfound]) == NO_ERROR)
#else
                if (gcaFindClosestMeanSample(gca, mean, min_prior,
                                             xi, yi, zi, label, 
                                             node_stride/3, &gcas[nfound]) == 
                    NO_ERROR)
#endif
                {
                  if (label == Gdiag_no)
                    DiagBreak() ;
                  found = 1 ;
                  if (label > 0)
                    total_mean += means[i][j][k][label] ;
                  else
                    nzeros++ ;
                  label_counts[label]++ ;
                  nfound++ ;
                }
              }
            }
          }
        }

        if (found)   /* add central label */
        {
          if (best_label == Gdiag_no)
            DiagBreak() ;
#if 0
          if (gcaFindBestSample(gca, x, y, z, best_label, node_stride/3, 
                                &gcas[nfound]) == NO_ERROR)
#else
          if (gcaFindClosestMeanSample(gca, best_mean, min_prior, x, y, z, 
                                       best_label, node_stride/3, 
                                       &gcas[nfound]) == 
              NO_ERROR)
#endif
          {
            if (best_label > 0)
              total_mean += means[1][1][1][best_label] ;
            else
              nzeros++ ;
            label_counts[best_label]++ ;
            nfound++ ;
          }
        }
      }
    }
  }

  fprintf(stderr, "total sample mean = %2.1f (%d zeros)\n", 
          total_mean/((float)nfound-nzeros), nzeros) ;

  {
    int  n ;
    FILE *fp ;

    fp = fopen("classes.dat", "w") ;
    for (n = 0 ; n < MAX_DIFFERENT_LABELS ; n++)
      if (label_counts[n] > 0)
        fprintf(fp, "%d  %d\n", n, label_counts[n]) ;
    fclose(fp) ;
  }

  *pnsamples = nfound ;

  return(gcas) ;
}
#if 1
GCA_SAMPLE *
GCAfindStableSamples(GCA *gca, int *pnsamples, int min_spacing,float min_prior)
{
  GCA_SAMPLE *gcas ;
  int        xi, yi, zi, width, height, depth, label, nfound,
             label_counts[MAX_DIFFERENT_LABELS], best_label, nzeros ;
  float      max_prior, total_mean, /*mean_dist,*/ best_mean_dist,
             priors[MAX_DIFFERENT_LABELS], means[MAX_DIFFERENT_LABELS], 
             vars[MAX_DIFFERENT_LABELS], max_priors[MAX_DIFFERENT_LABELS],
             node_stride, x, y, z, min_unknown, max_unknown ;
  MRI        *mri_filled ;

#define MIN_UNKNOWN_DIST  8

  gcaFindMaxPriors(gca, max_priors) ;
  gcaFindIntensityBounds(gca, &min_unknown, &max_unknown) ;
  printf("bounding unknown intensity as < %2.1f or > %2.1f\n", 
         min_unknown, max_unknown) ;

  memset(label_counts, 0, sizeof(label_counts)) ;
  width = gca->width ; height = gca->height ; depth = gca->depth ;
  gcas = calloc(width*height*depth, sizeof(GCA_SAMPLE)) ;

  mri_filled = MRIalloc(width*gca->spacing+1,height*gca->spacing+1, 
                        depth*gca->spacing+1,MRI_UCHAR);
  
  node_stride = (float)min_spacing / (float)gca->spacing ;

  total_mean = 0.0 ;

  /* compute the max priors and min variances for each class */
  for (nzeros = nfound = x = 0 ; x < width ; x += node_stride)
  {
    xi = nint(x) ;
    for (y = 0 ; y < height ; y += node_stride)
    {
      yi = nint(y) ;
      for (z = 0 ; z < depth ; z += node_stride)
      {
        zi = nint(z) ;
        if (abs(x-31)<=node_stride &&
            abs(y-22)<=node_stride &&
            abs(z-36)<=node_stride)
          DiagBreak() ;
        gcaRegionStats(gca, x, y, z, node_stride/2, priors, vars, means) ;

        if (priors[4] > .5*min_prior ||  /* left lat ven */
            priors[5] > .5*min_prior ||  /* inf left lat ven */
            priors[14] > .5*min_prior ||  /* 3rd ven */
            priors[15] > .5*min_prior ||  /* 4th ven */
            priors[43] > .5*min_prior ||
            priors[44] > .5*min_prior)
          DiagBreak() ;
            
        best_mean_dist = 0.0 ;
        max_prior = -1 ;

        if ((different_nbr_labels(gca, x, y, z, 
                                  ceil(UNKNOWN_DIST/gca->spacing),0) > 0) &&
            (priors[0] >= min_prior) &&
            (priors[0] >= 0.9*max_priors[0]))
          best_label = 0 ;
        else 
          best_label = -1 ;

        for (label = 1 ; label < MAX_DIFFERENT_LABELS ; label++)
        {
          if ((priors[label] < min_prior) || 
              (priors[label] < .9*max_priors[label]))
            continue ;

          if ((best_label == 0) ||
              (priors[label] > max_prior) ||
              (FEQUAL(priors[label], max_prior) && 
               label_counts[best_label] > label_counts[label]))
          {
            best_label = label ;
            max_prior = priors[label] ;
          }
        }
#if 1
        if (nfound > 0)
        {
          double p = randomNumber(0, 1.0) ;
          if (p < ((double)label_counts[best_label] / nfound))
            continue ;
        }
#endif
        if (best_label >= 0)
        {
          if (best_label == 0)
          {
#if 1
            if (gcaFindBestSample(gca, x, y, z, best_label, node_stride/2, 
                                  &gcas[nfound]) == NO_ERROR)
#else
            gcaFindClosestMeanSample(gca, 255, min_prior, x, y, z, 0,
                                     node_stride/2, &gcas[nfound]);
            if (means[0] > 100)
#endif
            {
              int  xv, yv, zv ;
              
              if (((means[0] <= min_unknown) ||
                   (means[0] >= max_unknown))) /* disabled check */
              {
                xv = (int)((float)x*gca->spacing) ;
                yv = (int)((float)y*gca->spacing) ;
                zv = (int)((float)z*gca->spacing) ;
                if (MRIvox(mri_filled, xv, yv, zv) == 0)
                {
                  mriFillRegion(mri_filled, xv, yv, zv, 1, MIN_UNKNOWN_DIST) ;
                  /*                MRIvox(mri_filled, xv, yv, zv) = 1 ;*/
                  nzeros++ ;
                  label_counts[best_label]++ ;
                  nfound++ ;
                }
              }
              else
                DiagBreak() ;
            }
          }
          else
          {
            /* find best gc with this label */
            if (gcaFindBestSample(gca, x, y, z, best_label, node_stride/2, 
                                  &gcas[nfound]) == NO_ERROR)
            {
              total_mean += means[best_label] ;
              label_counts[best_label]++ ;
              nfound++ ;
            }
          }
        }
      }
    }
  }

  fprintf(stderr, "total sample mean = %2.1f (%d zeros)\n", 
          total_mean/((float)nfound-nzeros), nzeros) ;

  {
    int  n ;
    FILE *fp ;

    fp = fopen("classes.dat", "w") ;
    for (n = 0 ; n < MAX_DIFFERENT_LABELS ; n++)
      if (label_counts[n] > 0)
        fprintf(fp, "%d  %d\n", n, label_counts[n]) ;
    fclose(fp) ;
  }

  *pnsamples = nfound ;

  MRIfree(&mri_filled) ;
  return(gcas) ;
}
#else
GCA_SAMPLE *
GCAfindStableSamples(GCA *gca, int *pnsamples, int min_spacing,float min_prior)
{
  GCA_SAMPLE *gcas ;
  GCA_NODE   *gcan ;
  GC1D       *gc, *best_gc ;
  int        x, y, z, width, height, depth, n, label, nfound, node_stride,
             label_counts[MAX_DIFFERENT_LABELS], best_n, best_label,
             xk, yk, zk, xi, yi, zi, best_x, best_y, best_z, nzeros ;
  float      max_prior, total_mean, mean_dist, best_mean_dist ;

  memset(label_counts, 0, sizeof(label_counts)) ;
  width = gca->width ; height = gca->height ; depth = gca->depth ;
  gcas = calloc(width*height*depth, sizeof(GCA_SAMPLE)) ;
  
  node_stride = min_spacing / gca->spacing ;

  total_mean = 0.0 ;

  /* compute the max priors and min variances for each class */
  for (nzeros = nfound = x = 0 ; x < width ; x += node_stride)
  {
    for (y = 0 ; y < height ; y += node_stride)
    {
      for (z = 0 ; z < depth ; z += node_stride)
      {
        if (abs(x-31)<=node_stride &&
            abs(y-22)<=node_stride &&
            abs(z-36)<=node_stride)
          DiagBreak() ;
        best_n = best_label = -1 ; best_gc = NULL ; best_mean_dist = 0.0 ;
        max_prior = -1 ; best_x = best_y = best_z = -1 ;
        for (xk = 0 ; xk < node_stride ; xk++)
        {
          xi = MIN(x+xk, width-1) ;
          for (yk = 0 ; yk < node_stride ; yk++)
          {
            yi = MIN(y+yk, height-1) ;
            for (zk = 0 ; zk < node_stride ; zk++)
            {
              zi = MIN(z+zk, depth-1) ;
              if (xi == 31 && yi == 22 && zi == 36)
                DiagBreak() ;
          
              gcan = &gca->nodes[xi][yi][zi] ;
              for (n = 0 ; n < gcan->nlabels ; n++)
              {
                gc = &gcan->gcs[n] ;
                label = gcan->labels[n] ;
                if (label == Gdiag_no)
                  DiagBreak() ;
                if ((label <= 0 && best_label >= 0) || gc->prior < min_prior)
                  continue ;

                if (label == 0 && 
                    !different_nbr_labels(gca, xi, yi, zi, 
                                        (int)ceil(UNKNOWN_DIST/gca->spacing),
                                          0))
                  continue ;
                
                if (nfound > nzeros)
                {
                  mean_dist = (gc->mean-total_mean/((float)nfound-nzeros)) ;
                  mean_dist = sqrt(mean_dist * mean_dist) ;
                  if (best_label >= 0 && 2*mean_dist < best_mean_dist)
                    continue ;
                }
                else
                  mean_dist = 0 ;

                if ((best_label == 0) ||
                    (gc->prior > max_prior) ||
                    (FEQUAL(gc->prior, max_prior) && 
                     label_counts[best_label] > label_counts[label]) ||
                    (mean_dist > 2*best_mean_dist)
#if 0
                    || label_counts[best_label] > 5*label_counts[label]
#endif
                    )
                {
                  if (label == Gdiag_no)
                    DiagBreak() ;
                  best_label = label ; best_gc = gc ;
                  max_prior = gc->prior ; best_n = n ;
                  if (nfound > nzeros)
                  {
                    best_mean_dist = 
                      (gc->mean-total_mean/(float)(nfound-nzeros)) ;
                    best_mean_dist = sqrt(best_mean_dist * best_mean_dist) ;
                  }
                  best_x = xi ; best_y = yi ; best_z = zi ;
                }
              }
            }
          }
        }

        if (best_label == 0 && 
            !different_nbr_labels(gca, best_x, best_y,best_z,
                                  (int)ceil(UNKNOWN_DIST/gca->spacing),0))
          continue ;
                
        if (x == 4 && y == 3 && z == 6)
          DiagBreak() ;
        if (best_n >= 0)
        {
          gcas[nfound].xn = best_x ; gcas[nfound].yn = best_y ;
          gcas[nfound].zn = best_z ;
          gcas[nfound].label = best_label ;
          gcas[nfound].n = best_n ;
          gcas[nfound].var = best_gc->var ;
          gcas[nfound].mean = best_gc->mean ;
          gcas[nfound].prior = best_gc->prior ;
          if (best_label > 0)
            total_mean += best_gc->mean ;
          else
            nzeros++ ;
          label_counts[best_label]++ ;
          nfound++ ;
        }
      }
    }
  }

  fprintf(stderr, "total sample mean = %2.1f (%d zeros)\n", 
          total_mean/((float)nfound-nzeros), nzeros) ;

  {
    int  n ;
    FILE *fp ;

    fp = fopen("classes.dat", "w") ;
    for (n = 0 ; n < MAX_DIFFERENT_LABELS ; n++)
      if (label_counts[n] > 0)
        fprintf(fp, "%d  %d\n", n, label_counts[n]) ;
    fclose(fp) ;
  }

  *pnsamples = nfound ;

  return(gcas) ;
}
#endif
int
GCAwriteSamples(GCA *gca, MRI *mri, GCA_SAMPLE *gcas, int nsamples,
                char *fname)
{
  int    n, xv, yv, zv ;
  MRI    *mri_dst ;

  mri_dst = MRIclone(mri, NULL) ;

  for (n = 0 ; n < nsamples ; n++)
  {
    if (gcas[n].label == Gdiag_no)
      DiagBreak() ;
    GCAnodeToVoxel(gca, mri_dst, gcas[n].xn, gcas[n].yn, gcas[n].zn, 
                   &xv, &yv, &zv) ;
    if (gcas[n].label > 0)
      MRIvox(mri_dst, xv, yv, zv) = gcas[n].label ;
    else
      MRIvox(mri_dst, xv, yv, zv) = 29 ;  /* Left undetermined - visible */
    if (DIAG_VERBOSE_ON && Gdiag & DIAG_SHOW)
      printf("label %d: (%d, %d, %d) <-- (%d, %d, %d)\n",
             gcas[n].label,gcas[n].xn,gcas[n].yn,gcas[n].zn, xv, yv, zv) ;
  }
  MRIwrite(mri_dst, fname) ;
  MRIfree(&mri_dst) ;
  return(NO_ERROR) ;
}
MRI *
GCAmri(GCA *gca, MRI *mri)
{
  int      width, height, depth, x, y, z, xn, yn, zn, n ;
  float    val ;
  GC1D     *gc ;
  GCA_NODE *gcan ;

  if (!mri)
    mri = MRIalloc(gca->width, gca->height, gca->depth, MRI_UCHAR) ;
  width = mri->width ; height = mri->height ; depth = mri->depth ;

  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        GCAvoxelToNode(gca, mri, x, y, z, &xn, &yn, &zn) ;
        gcan = &gca->nodes[xn][yn][zn] ;
        for (val = 0.0, n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          val += gc->mean*gc->prior ;
        }
        switch (mri->type)
        {
        default:
          ErrorReturn(NULL,
                      (ERROR_UNSUPPORTED, 
                       "GCAmri: unsupported image type %d", mri->type)) ;
        case MRI_UCHAR:
          MRIvox(mri, x, y, z) = (unsigned char)val ;
          break ;
        case MRI_SHORT:
          MRISvox(mri, x, y, z) = (short)val ;
          break ;
        }
      }
    }
  }
  return(mri) ;
}
MRI *
GCAlabelMri(GCA *gca, MRI *mri, int label, LTA *lta)
{
  int      width, height, depth, x, y, z, xn, yn, zn, n ;
  float    val ;
  GC1D     *gc ;
  GCA_NODE *gcan ;
  MRI      *mri_norm ;

  if (!mri)
    mri = MRIalloc(gca->width, gca->height, gca->depth, MRI_UCHAR) ;
  width = mri->width ; height = mri->height ; depth = mri->depth ;
  mri_norm = MRIalloc(width, height, depth, MRI_SHORT) ;

  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (x == width/2 && y == height/2 && z == depth/2)
          DiagBreak() ;
        if (x == 19 && y == 14 && z == 15)
          DiagBreak() ;

        GCAsourceVoxelToNode(gca, mri, lta, x, y, z, &xn, &yn, &zn) ;
        gcan = &gca->nodes[xn][yn][zn] ;
        if (gcan->nlabels < 1)
          continue ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          if (gcan->labels[n] == label)
            break ;
        }
        if (n >= gcan->nlabels)
          continue ;
        val = gc->mean ;
        switch (mri->type)
        {
        default:
          ErrorReturn(NULL,
                      (ERROR_UNSUPPORTED, 
                       "GCAlabelMri: unsupported image type %d", mri->type)) ;
        case MRI_UCHAR:
          MRIvox(mri, x, y, z) = (unsigned char)val ;
          break ;
        case MRI_SHORT:
          MRISvox(mri, x, y, z) = (short)val ;
          break ;
        case MRI_FLOAT:
          MRIFvox(mri, x, y, z) = (float)val ;
          break ;
        }
        MRISvox(mri_norm, x, y, z) += 1 ;
      }
    }
  }
  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        if (MRISvox(mri_norm, x, y, z) > 0)
          MRISvox(mri, x, y, z) = 
            nint((float)MRISvox(mri,x,y,z)/(float)MRISvox(mri_norm,x,y,z)) ;
      }
    }
  }
  MRIfree(&mri_norm) ;
  return(mri) ;
}

static int 
different_nbr_labels(GCA *gca, int x, int y, int z, int wsize, int label)
{
  int      xk, yk, zk, xi, yi, zi, nbrs, n, width, height, depth ;
  GCA_NODE *gcan ;

  width = gca->width ; height = gca->height ; depth = gca->depth ;
  for (nbrs = 0, xk = -wsize ; xk <= wsize ; xk++)
  {
    xi = x+xk ;
    if (xi < 0 || xi >= width)
      continue ;

    for (yk = -wsize ; yk <= wsize ; yk++)
    {
      yi = y+yk ;
      if (yi < 0 || yi >= height)
        continue ;

      for (zk = -wsize ; zk <= wsize ; zk++)
      {
        zi = z+zk ;
        if (zi < 0 || zi >= depth)
          continue ;
        
        gcan = &gca->nodes[xi][yi][zi] ;

        for (n = 0 ; n < gcan->nlabels ; n++)
          if (gcan->labels[n] != label)
            nbrs++ ;
      }
    }
  }
  return(nbrs) ;
}
static int
gcaRegionStats(GCA *gca, int x, int y, int z, int wsize,
               float *priors, float *vars, float *means)
{
  int        xk, yk, zk, xi, yi, zi, width, height, depth, label, n,
             total_training ;
  GCA_NODE  *gcan  ;
  GC1D      *gc ;
  float     dof ;

  if (x == 28 && y == 24 && z == 36)
    DiagBreak() ;

  memset(priors, 0, MAX_DIFFERENT_LABELS*sizeof(float)) ;
  memset(vars, 0, MAX_DIFFERENT_LABELS*sizeof(float)) ;
  memset(means, 0, MAX_DIFFERENT_LABELS*sizeof(float)) ;

  width = gca->width ; height = gca->height ; depth = gca->depth ;
  total_training = 0 ;
  for (xk = -wsize ; xk <= wsize ; xk++)
  {
    xi = x+xk ;
    if (xi < 0 || xi >= width)
      continue ;

    for (yk = -wsize ; yk <= wsize ; yk++)
    {
      yi = y+yk ;
      if (yi < 0 || yi >= height)
        continue ;

      for (zk = -wsize ; zk <= wsize ; zk++)
      {
        zi = z+zk ;
        if (zi < 0 || zi >= depth)
          continue ;
        
        gcan = &gca->nodes[xi][yi][zi] ;
        total_training += gcan->total_training ;

        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          label = gcan->labels[n] ;
          if (label == Gdiag_no)
            DiagBreak() ;
          gc = &gcan->gcs[n] ;
          dof = (gc->prior * (float)gcan->total_training) ;
          
          vars[label] +=  dof * gc->var ;
          means[label] += dof * gc->mean ;
          priors[label] += dof ;
        }
      }
    }
  }

  for (label = 0 ; label < MAX_DIFFERENT_LABELS ; label++)
  {
    dof = priors[label] ;
    if (dof > 0)
    {
      vars[label] /= dof ;
      priors[label] /= total_training ;
      means[label] /= dof ;
    }
  }
  return(NO_ERROR) ;
}
static int
gcaFindBestSample(GCA *gca, int x, int y, int z,int label,int wsize,
                  GCA_SAMPLE *gcas)
{
  int        xk, yk, zk, xi, yi, zi, width, height, depth, n,
             best_n, best_x, best_y, best_z ;
  GCA_NODE   *gcan  ;
  GC1D       *gc, *best_gc ;
  float      max_prior, min_var ;

  width = gca->width ; height = gca->height ; depth = gca->depth ;
  max_prior = 0.0 ; min_var = 10000.0f ; best_gc = NULL ;
  best_x = best_y = best_z = -1 ; best_n = -1 ;
  for (xk = -wsize ; xk <= wsize ; xk++)
  {
    xi = x+xk ;
    if (xi < 0 || xi >= width)
      continue ;

    for (yk = -wsize ; yk <= wsize ; yk++)
    {
      yi = y+yk ;
      if (yi < 0 || yi >= height)
        continue ;

      for (zk = -wsize ; zk <= wsize ; zk++)
      {
        zi = z+zk ;
        if (zi < 0 || zi >= depth)
          continue ;
        
        gcan = &gca->nodes[xi][yi][zi] ;

        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          if (gcan->labels[n] != label)
            continue ;
          DiagBreak() ;
          gc = &gcan->gcs[n] ;
          if (gc->prior > max_prior ||
              (FEQUAL(gc->prior, max_prior) && (gc->var < min_var)) ||
              (label == 0 && gc->mean > best_gc->mean))
          {
            max_prior = gc->prior ; min_var = gc->var ;
            best_gc = gc ; best_x = xi ; best_y = yi ; best_z = zi ; 
            best_n = n ;
          }
        }
      }
    }
  }

  if (best_x < 0)
  {
    ErrorPrintf(ERROR_BADPARM, 
                "could not find GC1D for label %d at (%d,%d,%d)\n", 
                label, x, y, z) ;
    return(ERROR_BADPARM) ;
  }
  if (best_x == 145/4 && best_y == 89/4 && best_z == 125/4)
    DiagBreak() ;
  gcas->xn = best_x ; gcas->yn = best_y ; gcas->zn = best_z ;
  gcas->label = label ; gcas->n = best_n ; gcas->mean = best_gc->mean ;
  gcas->var = best_gc->var ; gcas->prior = best_gc->prior ;

  return(NO_ERROR) ;
}
static GC1D *
gcaFindHighestPriorGC(GCA *gca, int x, int y, int z,int label,int wsize)
{
  int        xk, yk, zk, xi, yi, zi, width, height, depth, n ;
  GCA_NODE   *gcan  ;
  GC1D       *gc, *best_gc ;
  float      max_prior ;

  width = gca->width ; height = gca->height ; depth = gca->depth ;
  max_prior = 0.0 ; best_gc = NULL ;
  for (xk = -wsize ; xk <= wsize ; xk++)
  {
    xi = x+xk ;
    if (xi < 0 || xi >= width)
      continue ;

    for (yk = -wsize ; yk <= wsize ; yk++)
    {
      yi = y+yk ;
      if (yi < 0 || yi >= height)
        continue ;

      for (zk = -wsize ; zk <= wsize ; zk++)
      {
        zi = z+zk ;
        if (zi < 0 || zi >= depth)
          continue ;
        
        gcan = &gca->nodes[xi][yi][zi] ;

        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          if (gcan->labels[n] != label)
            continue ;
          gc = &gcan->gcs[n] ;
          if (gc->prior > max_prior)
          {
            max_prior = gc->prior ; 
            best_gc = gc ; 
          }
        }
      }
    }
  }

  if (best_gc == NULL)
  {
    ErrorPrintf(ERROR_BADPARM, 
                "could not find GC for label %d at (%d,%d,%d)\n", 
                label, x, y, z) ;
    return(NULL) ;
  }

  return(best_gc) ;
}

static int
gcaFindClosestMeanSample(GCA *gca, float mean, float min_prior, int x, int y, 
                         int z, int label, int wsize, GCA_SAMPLE *gcas)
{
  int        xk, yk, zk, xi, yi, zi, width, height, depth, n,
             best_n, best_x, best_y, best_z ;
  GCA_NODE   *gcan  ;
  GC1D       *gc, *best_gc ;
  float      min_dist ;

  width = gca->width ; height = gca->height ; depth = gca->depth ;
  min_dist = 1000000.0f ; best_gc = NULL ;
  best_x = best_y = best_z = -1 ; best_n = -1 ;
  for (xk = -wsize ; xk <= wsize ; xk++)
  {
    xi = x+xk ;
    if (xi < 0 || xi >= width)
      continue ;

    for (yk = -wsize ; yk <= wsize ; yk++)
    {
      yi = y+yk ;
      if (yi < 0 || yi >= height)
        continue ;

      for (zk = -wsize ; zk <= wsize ; zk++)
      {
        zi = z+zk ;
        if (zi < 0 || zi >= depth)
          continue ;
        
        gcan = &gca->nodes[xi][yi][zi] ;

        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          if (gcan->labels[n] != label || gc->prior < min_prior)
            continue ;
          if (abs(gc->mean - mean) < min_dist)
          {
            min_dist = abs(gc->mean-mean) ;
            best_gc = gc ; best_x = xi ; best_y = yi ; best_z = zi ; 
            best_n = n ;
          }
        }
      }
    }
  }

  if (best_x < 0)
  {
    ErrorPrintf(ERROR_BADPARM, 
                "could not find GC for label %d at (%d,%d,%d)\n", 
                label, x, y, z) ;
    return(ERROR_BADPARM) ;
  }
  if (best_x == 141/4 && best_y == 37*4 && best_z == 129*4)
    DiagBreak() ;
  gcas->xn = best_x ; gcas->yn = best_y ; gcas->zn = best_z ;
  gcas->label = label ; gcas->n = best_n ; gcas->mean = best_gc->mean ;
  gcas->var = best_gc->var ; gcas->prior = best_gc->prior ;

  return(NO_ERROR) ;
}
int
GCAtransformAndWriteSamples(GCA *gca, MRI *mri, GCA_SAMPLE *gcas, 
                            int nsamples,char *fname,LTA *lta)
{
  int    n, xv, yv, zv, label ;
  MRI    *mri_dst ;
  MATRIX *m_L_inv ;
  VECTOR *v_input, *v_canon ;

  m_L_inv = MatrixInverse(lta->xforms[0].m_L, NULL) ;
  if (!m_L_inv)
  {
    MatrixPrint(stderr, lta->xforms[0].m_L) ;
    ErrorExit(ERROR_BADPARM, 
              "GCAcomputeLogSampleProbability: matrix is not invertible") ;
  }
  mri_dst = MRIclone(mri, NULL) ;
  
  v_input = VectorAlloc(4, MATRIX_REAL) ;
  v_canon = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v_input, 4, 1) = 1.0 ; *MATRIX_RELT(v_canon, 4, 1) = 1.0 ;

  for (n = 0 ; n < nsamples ; n++)
  {
    if (gcas[n].label == Gdiag_no)
      DiagBreak() ;
    GCAnodeToVoxel(gca, mri_dst, gcas[n].xn, gcas[n].yn, gcas[n].zn, 
                   &xv, &yv, &zv) ;
    V3_X(v_canon) = (float)xv ; V3_Y(v_canon) = (float)yv ;
    V3_Z(v_canon) = (float)zv ;
    MatrixMultiply(m_L_inv, v_canon, v_input) ;
    xv = nint(V3_X(v_input)); yv = nint(V3_Y(v_input)); zv=nint(V3_Z(v_input));
    if (xv < 0) xv = 0 ;
    if (xv >= mri->width) xv = mri->width-1 ;
    if (yv < 0) yv = 0 ;
    if (yv >= mri->height) yv = mri->height-1 ;
    if (zv < 0) zv = 0 ;
    if (zv >= mri->depth) zv = mri->depth-1 ;
#if 0
    if (gcas[n].label > 0)
      MRIvox(mri_dst, xv, yv, zv) = gcas[n].label ;
    else
      MRIvox(mri_dst, xv, yv, zv) = 29 ;  /* Left undetermined - visible */
#else
    if (gcas[n].label > 0)
      label = gcas[n].label ;
    else
      label = 29 ;  /* Left undetermined - visible */
    mriFillRegion(mri_dst, xv, yv, zv, label, 0) ;
#endif
    gcas[n].x = xv ;
    gcas[n].y = yv ;
    gcas[n].z = zv ;
    if (DIAG_VERBOSE_ON && Gdiag & DIAG_SHOW)
      printf("label %d: (%d, %d, %d) <-- (%d, %d, %d)\n",
             gcas[n].label,gcas[n].xn,gcas[n].yn,gcas[n].zn, xv, yv, zv) ;
  }
  fprintf(stderr, "writing samples to %s...\n", fname) ;
  MRIwrite(mri_dst, fname) ;
  MRIfree(&mri_dst) ;

  MatrixFree(&m_L_inv) ; MatrixFree(&v_canon) ; MatrixFree(&v_input) ;

  return(NO_ERROR) ;
}

static int
mriFillRegion(MRI *mri, int x, int y, int z, int fill_val, int whalf)
{
  int   xi, xk, yi, yk, zi, zk ;

  for (xk = -whalf ; xk <= whalf ; xk++)
  {
    xi = mri->xi[x+xk] ;
    for (yk = -whalf ; yk <= whalf ; yk++)
    {
      yi = mri->xi[y+yk] ;
      for (zk = -whalf ; zk <= whalf ; zk++)
      {
        zi = mri->zi[z+zk] ;
        MRIvox(mri, xi, yi, zi) = fill_val ;
      }
    }
  }
  return(NO_ERROR) ;
}

static int
gcaFindIntensityBounds(GCA *gca, float *pmin, float *pmax)
{
  int      width, depth, height, x, y, z, n, label ;
  GCA_NODE *gcan ;
  GC1D     *gc ;
  float    mn, mx, offset ;

  mn = 100000.0f ; mx = -mn ;
  width = gca->width ; height = gca->height ; depth = gca->depth ;

  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        gcan = &gca->nodes[x][y][z] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          label = gcan->labels[n] ;
          if (label == 0)  /* don't include unknowns */
            continue ;
          if (gc->prior < 0.1)
            continue ;  /* exclude unlikely stuff (errors in labeling) */
          offset = 0.5*sqrt(gc->var) ;
          if (gc->mean + offset > mx)
          {
            mx = gc->mean + offset ;
#if 0
            if (mx > 130)
            {
              printf("label %d, mean %2.1f, std %2.1f, prior %2.3f\n",
                      gcan->labels[n], gc->mean, sqrt(gc->var), gc->prior) ;
              
            }
#endif
          }
          if (gc->mean  < mn)
          {
            mn = gc->mean ;
#if 0
            if (mn < 10)
            {
              printf("label %d, mean %2.1f, std %2.1f, prior %2.3f\n",
                      gcan->labels[n], gc->mean, sqrt(gc->var), gc->prior) ;
              
            }
#endif
          }
          if (mn < 5)
            mn = 5 ;
          if (mx > 225)
            mx = 225 ;
        }
      }
    }
  }
  *pmin = mn ; *pmax = mx ;
  return(NO_ERROR) ;
}
static int
gcaFindMaxPriors(GCA *gca, float *max_priors)
{
  int      width, depth, height, x, y, z, n, label ;
  GCA_NODE *gcan ;
  GC1D     *gc ;

  memset(max_priors, 0, MAX_DIFFERENT_LABELS*sizeof(float)) ;
  width = gca->width ; height = gca->height ; depth = gca->depth ;

  for (x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        gcan = &gca->nodes[x][y][z] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          label = gcan->labels[n] ;
          if (gc->prior > max_priors[label])
            max_priors[label] = gc->prior ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}

static int xn_bad = 15 ;
static int yn_bad = 24 ;
static int zn_bad = 27 ;
static int xl_bad = 56 ;
static int yl_bad = 99 ;
static int zl_bad = 106 ;
static int label_bad = 42 ;
static int nbr_label_bad = 42 ;
static int bad_i = 0 ;

static int
GCAupdateNodeGibbsPriors(GCA *gca, MRI*mri, int xn, int yn, int zn,
                   int xl, int yl, int zl, int label)
{
  int       n, i, xnbr, ynbr, znbr, nbr_label ;
  GCA_NODE *gcan ;
  GC1D     *gc ;


  if (xl == xl_bad && yl == yl_bad && zl == zl_bad)
    DiagBreak() ;
  if ((xn == 16 && yn == 26 && zn == 27 && label == 0))
    DiagBreak() ;
  if (((xn == xn_bad && yn == yn_bad && zn == zn_bad && label == label_bad)))
    DiagBreak() ;

  gcan = &gca->nodes[xn][yn][zn] ;

  for (n = 0 ; n < gcan->nlabels ; n++)
  {
    if (gcan->labels[n] == label)
      break ;
  }
  if (n >= gcan->nlabels)  /* have to allocate a new classifier */
    ErrorExit(ERROR_BADPARM, "gca(%d, %d, %d): could not find label %d",
              xn, yn, zn, label) ;

  gc = &gcan->gcs[n] ;
  for (i = 0 ; i < GIBBS_NEIGHBORHOOD ; i++)
  {
    /* coordinates of neighboring point */
    xnbr = mri->xi[xl+xnbr_offset[i]] ;
    ynbr = mri->yi[yl+ynbr_offset[i]] ;
    znbr = mri->zi[zl+znbr_offset[i]] ;

    nbr_label = MRIvox(mri, xnbr, ynbr, znbr) ;
    if (xn == xn_bad && yn == yn_bad && zn == zn_bad && label == label_bad &&
        nbr_label == nbr_label_bad)
      DiagBreak() ;
    if (xn == xn_bad && yn == yn_bad && zn == zn_bad && label == label_bad &&
        nbr_label == nbr_label_bad && i == bad_i)
      DiagBreak() ;

    /* now see if this label exists already as a nbr */
    for (n = 0 ; n < gc->nlabels[i] ; n++)
      if (gc->labels[i][n] == nbr_label)
        break ;

    if (xn == 16 && yn == 26 && zn == 27 && n > 0)
      DiagBreak() ;
    if (((xn == 16 && yn == 25 && zn == 27 && label == 42)) && (n > 0))
      DiagBreak() ;

    if (n >= gc->nlabels[i])   /* not there - reallocate stuff */
    {
#if 1
      char *old_labels ;
      float *old_label_priors ;

      if (xl == 54 && yl == 105 && zl == 112)
        DiagBreak() ;

      old_labels = gc->labels[i] ;
      old_label_priors = gc->label_priors[i] ;

      /* allocate new ones */
      gc->label_priors[i] = (float *)calloc(gc->nlabels[i]+1, sizeof(float)) ;
      if (!gc->label_priors[i])
        ErrorExit(ERROR_NOMEMORY, "GCAupdateNodeGibbsPriors: "
                  "couldn't expand gcs to %d" ,gc->nlabels[i]+1) ;
      gc->labels[i] = (char *)calloc(gc->nlabels[i]+1, sizeof(char)) ;
      if (!gc->labels[i])
        ErrorExit(ERROR_NOMEMORY, 
                  "GCANupdateNode: couldn't expand labels to %d",
                  gc->nlabels[i]+1) ;

      if (gc->nlabels[i] > 0)  /* copy the old ones over */
      {
        memmove(gc->label_priors[i], old_label_priors, 
                gc->nlabels[i]*sizeof(float)) ;
        memmove(gc->labels[i], old_labels, gc->nlabels[i]*sizeof(char)) ;

        /* free the old ones */
        free(old_label_priors) ; free(old_labels) ;
      }
      gc->labels[i][gc->nlabels[i]++] = nbr_label ;
#else
      if (n >= MAX_NBR_LABELS)   /* not there - reallocate stuff */
      {
        ErrorPrintf(ERROR_NOMEMORY, 
                    "GCAupdateNodeGibbsPriors: gca(%d,%d,%d) has more than %d "
                    "different neighbors", xn, yn, zn, MAX_NBR_LABELS);
        continue ;
      }
      else
        gc->labels[i][gc->nlabels[i]++] = nbr_label ;
#endif
    }
    gc->label_priors[i][n] += 1.0f ;
  }

  return(NO_ERROR) ;
}

#if 0
MRI  *
GCAreclassifyUsingGibbsPriors(MRI *mri_inputs, GCA *gca, MRI *mri_dst,LTA *lta,
                            int max_iter)
{
  int      x, y, z, width, height, depth, label, val, iter,
           xn, yn, zn, n, i, j, nchanged, xnbr, ynbr, znbr, nbr_label,
           index, nindices ;
  short    *x_indices, *y_indices, *z_indices ;
  GCA_NODE *gcan ;
  GC1D     *gc ;
  double   dist, max_p, p, prior, ll, lcma = 0.0, new_ll, old_ll ;
  MRI      *mri_changed, *mri_probs ;

  nindices = mri_dst->width * mri_dst->height * mri_dst->depth ;
  x_indices = (short *)calloc(nindices, sizeof(short)) ;
  y_indices = (short *)calloc(nindices, sizeof(short)) ;
  z_indices = (short *)calloc(nindices, sizeof(short)) ;
  if (!x_indices || !y_indices || !z_indices)
    ErrorExit(ERROR_NOMEMORY, "GCAreclassifyUsingGibbsPriors: "
              "could not allocate index set") ;


  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_inputs, NULL) ;
    if (!mri_dst)
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst") ;
    MRIcopyHeader(mri_inputs, mri_dst) ;
  }

  mri_changed = MRIclone(mri_dst, NULL) ;


  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_inputs->width ; height = mri_inputs->height; 
  depth = mri_inputs->depth ; iter = 0 ;
  for (x = 0 ; x < width ; x++)
    for (y = 0 ; y < height ; y++)
        for (z = 0 ; z < depth ; z++)
          MRIvox(mri_changed,x,y,z) = 1 ;

#if 0
  {
    MRI   *mri_cma ;
    char  fname[STRLEN], *cp ;

    strcpy(fname, mri_inputs->fname) ;
    cp = strrchr(fname, '/') ;
    if (cp)
    {
      strcpy(cp+1, "parc") ;
      mri_cma = MRIread(fname) ;
      if (mri_cma)
      {

        ll = gcaGibbsImageLogLikelihood(gca, mri_dst, mri_inputs, lta) ;
        lcma = gcaGibbsImageLogLikelihood(gca, mri_cma, mri_inputs, lta) ;
        lcma /= (double)(width*depth*height) ;
        ll /= (double)(width*depth*height) ;
        fprintf(stderr, "image ll: %2.3f (CMA=%2.3f)\n", ll, lcma) ;
        MRIfree(&mri_cma) ;
      }
    }
  }
#endif

  do
  {
    if (iter == 0)
    {
      mri_probs = GCAlabelProbabilities(mri_inputs, gca, NULL, lta) ;
      MRIorderIndices(mri_probs, x_indices, y_indices, z_indices) ;
      MRIfree(&mri_probs) ;
    }
    else
      MRIcomputeVoxelPermutation(mri_inputs, x_indices, y_indices,
                                 z_indices) ;
      
    nchanged = 0 ;
    for (index = 0 ; index < nindices ; index++)
    {
      x = x_indices[index] ; y = y_indices[index] ; z = z_indices[index] ;
      GCAsourceVoxelToNode(gca, mri_inputs, lta, x, y, z, &xn, &yn, &zn) ;
      if (x == 100 && y == 104 && z == 130)
        DiagBreak() ;

#if 1
      if (MRIvox(mri_changed, x, y, z) == 0)
        continue ;
#endif

      if (x == 63 && y == 107 && z == 120)
        DiagBreak() ; 
          
      val = MRIvox(mri_inputs, x, y, z) ;
          
          
      /* find the node associated with this coordinate and classify */
      GCAsourceVoxelToNode(gca, mri_inputs, lta, x, y, z, &xn, &yn, &zn) ;
      gcan = &gca->nodes[xn][yn][zn] ;
      label = 0 ; max_p = BIG_AND_NEGATIVE ;
#if 0
      if (x == 141 && y == 62 && z == 126 && iter == 0)
      {
        printf("val = %d\n", val) ;
        dump_gcan(gcan, stdout, 0) ;
      }
#endif
      for (n = 0 ; n < gcan->nlabels ; n++)
      {
        gc = &gcan->gcs[n] ;
        
        /* compute 1-d Mahalanobis distance */
        dist = (val-gc->mean) ;
#define USE_LOG_PROBS  1
        p = -log(sqrt(gc->var)) - .5*(dist*dist/gc->var) + log(gc->prior);
        
        for (prior = 0.0f, i = 0 ; i < GIBBS_NEIGHBORS ; i++)
        {
          xnbr = mri_dst->xi[x+xnbr_offset[i]] ;
          ynbr = mri_dst->yi[y+ynbr_offset[i]] ;
          znbr = mri_dst->zi[z+znbr_offset[i]] ;
          nbr_label = MRIvox(mri_dst, xnbr, ynbr, znbr) ;
          for (j = 0 ; j < gc->nlabels[i] ; j++)
          {
            if (nbr_label == gc->labels[i][j])
              break ;
          }
          if (j < gc->nlabels[i])
          {
            if (!FZERO(gc->label_priors[i][j]))
              prior += log(gc->label_priors[i][j]) ;
            else
              prior += BIG_AND_NEGATIVE ;
          }
          else
          {
            if (x == 141 && y == 62 && z == 126)
              DiagBreak() ;
            prior += BIG_AND_NEGATIVE ; /* never occurred - make it unlikely */
          }
        }
        p += prior ;
        if (p > max_p)
        {
          max_p = p ;
          label = gcan->labels[n] ;
        }
      }
      
      if (FZERO(max_p))
      {
        label = 0 ; max_p = BIG_AND_NEGATIVE ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          
          /* compute 1-d Mahalanobis distance */
          dist = (val-gc->mean) ;
          if (FZERO(gc->var))  /* make it a delta function */
          {
            if (FZERO(dist))
              p = 1.0 ;
            else
              p = 0.0 ;
          }
          else
            p = 1 / sqrt(gc->var * 2 * M_PI) * exp(-dist*dist/gc->var) ;
          
          p *= gc->prior ;
          if (p > max_p)
          {
            max_p = p ;
            label = gcan->labels[n] ;
          }
        }
      }
      
      if (x == Ggca_x && y == Ggca_y && z == Ggca_z)
      {
        printf(
             "(%d, %d, %d): old label %s (%d), new label %s (%d) (p=%2.3f)\n",
             x, y, z, cma_label_to_name(MRIvox(mri_dst,x,y,z)),
             MRIvox(mri_dst,x,y,z), cma_label_to_name(label), label, max_p) ;
        if (label == Ggca_label)
        {
          if (DIAG_VERBOSE_ON)
            dump_gcan(gcan, stderr, 1) ;
          DiagBreak() ;
        }
      }
      
      if (MRIvox(mri_dst, x, y, z) != label)
      {
        int old_label = MRIvox(mri_dst, x, y, z) ;
        if (x == 100 && y == 104 && z == 130)
          DiagBreak() ;
        old_ll = 
          gcaNbhdGibbsLogLikelihood(gca, mri_dst, mri_inputs, x, y, z, lta,
                                    PRIOR_FACTOR) ;
        MRIvox(mri_dst, x, y, z) = label ;
        new_ll = 
          gcaNbhdGibbsLogLikelihood(gca, mri_dst, mri_inputs, x, y, z, lta,
                                    PRIOR_FACTOR) ;
        if (new_ll > old_ll)
        {
          MRIvox(mri_changed, x, y, z) = 1 ;
          nchanged++ ;
        }
        else
        {
          MRIvox(mri_dst, x, y, z) = old_label ;
          MRIvox(mri_changed, x, y, z) = 0 ;
        }
      }
      else
        MRIvox(mri_changed, x, y, z) = 0 ;
    }
    ll = gcaGibbsImageLogLikelihood(gca, mri_dst, mri_inputs, lta) ;
    ll /= (double)(width*depth*height) ;
    if (!FZERO(lcma))
      printf("pass %d: %d changed. image ll: %2.3f (CMA=%2.3f)\n", 
             iter+1, nchanged, ll, lcma) ;
    else
      printf("pass %d: %d changed. image ll: %2.3f\n", 
             iter+1, nchanged, ll) ;
    MRIdilate(mri_changed, mri_changed) ;

  } while (nchanged > 0 && iter++ < max_iter) ;

  free(x_indices) ; free(y_indices) ; free(z_indices) ;
  MRIfree(&mri_changed) ;


  return(mri_dst) ;
}
#endif

MRI  *
GCAanneal(MRI *mri_inputs, GCA *gca, MRI *mri_dst,LTA *lta, 
          int max_iter)
{
  int      x, y, z, width, height, depth, label, val, iter,
           xn, yn, zn, n, nchanged,
           index, nindices, old_label ;
  short    *x_indices, *y_indices, *z_indices ;
  GCA_NODE *gcan ;
  double   ll, lcma = 0.0, old_ll, new_ll, min_ll ;
  MRI      *mri_changed, *mri_probs ;

  printf("performing simulated annealing...\n") ;

  nindices = mri_dst->width * mri_dst->height * mri_dst->depth ;
  x_indices = (short *)calloc(nindices, sizeof(short)) ;
  y_indices = (short *)calloc(nindices, sizeof(short)) ;
  z_indices = (short *)calloc(nindices, sizeof(short)) ;
  if (!x_indices || !y_indices || !z_indices)
    ErrorExit(ERROR_NOMEMORY, "GCAanneal: "
              "could not allocate index set") ;


  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_inputs, NULL) ;
    if (!mri_dst)
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst") ;
    MRIcopyHeader(mri_inputs, mri_dst) ;
  }

  mri_changed = MRIclone(mri_dst, NULL) ;

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_inputs->width ; height = mri_inputs->height; 
  depth = mri_inputs->depth ; iter = 0 ;
  for (x = 0 ; x < width ; x++)
    for (y = 0 ; y < height ; y++)
        for (z = 0 ; z < depth ; z++)
          MRIvox(mri_changed,x,y,z) = 1 ;

  old_ll = gcaGibbsImageLogLikelihood(gca, mri_dst, mri_inputs, lta) ;
  old_ll /= (double)(width*depth*height) ;
#if 0
  {
    MRI   *mri_cma ;
    char  fname[STRLEN], *cp ;

    strcpy(fname, mri_inputs->fname) ;
    cp = strrchr(fname, '/') ;
    if (cp)
    {
      strcpy(cp+1, "parc") ;
      mri_cma = MRIread(fname) ;
      if (mri_cma)
      {

        lcma = gcaGibbsImageLogLikelihood(gca, mri_cma, mri_inputs, lta) ;
        lcma /= (double)(width*depth*height) ;
        fprintf(stderr, "image ll: %2.3f (CMA=%2.3f)\n", old_ll, lcma) ;
        MRIfree(&mri_cma) ;
      }
    }
  }
#endif

  do
  {
    if (iter == 0)
    {
      mri_probs = GCAlabelProbabilities(mri_inputs, gca, mri_dst, lta) ;
      MRIorderIndices(mri_probs, x_indices, y_indices, z_indices) ;
      MRIfree(&mri_probs) ;
    }
    else
      MRIcomputeVoxelPermutation(mri_inputs, x_indices, y_indices,
                                 z_indices) ;
      
    nchanged = 0 ;
    for (index = 0 ; index < nindices ; index++)
    {
      x = x_indices[index] ; y = y_indices[index] ; z = z_indices[index] ;
      if (x == 155 && y == 126 && z == 128)
        DiagBreak() ;

#if 1
      if (MRIvox(mri_changed, x, y, z) == 0)
        continue ;
#endif

      if (x == 63 && y == 107 && z == 120)
        DiagBreak() ; 
          
      val = MRIvox(mri_inputs, x, y, z) ;
          
      /* find the node associated with this coordinate and classify */
      GCAsourceVoxelToNode(gca, mri_inputs, lta, x, y, z, &xn, &yn, &zn) ;
      gcan = &gca->nodes[xn][yn][zn] ;


      label = old_label = MRIvox(mri_dst, x, y, z) ;
      min_ll = gcaNbhdGibbsLogLikelihood(gca, mri_dst, mri_inputs, x, y,z,lta,
                                         PRIOR_FACTOR);

      for (n = 0 ; n < gcan->nlabels ; n++)
      {
        if (gcan->labels[n] == old_label)
          continue ;
        MRIvox(mri_dst, x, y, z) = gcan->labels[n] ;
        new_ll = 
          gcaNbhdGibbsLogLikelihood(gca, mri_dst, mri_inputs, x, y,z,lta,
                                    PRIOR_FACTOR);
        if (new_ll > min_ll)
        {
          min_ll = new_ll ; label = gcan->labels[n] ;
        }
      }
      if (label != old_label)
      {
        nchanged++ ;
        MRIvox(mri_changed, x, y, z) = 1 ;
      }
      else
        MRIvox(mri_changed, x, y, z) = 0 ;
      MRIvox(mri_dst, x, y, z) = label ;
      
#if 0
      if (x == 155 && y == 126 && z == 128)
      {
        printf("(%d, %d, %d): old label %d, new label %d (p=%2.3f)\n",
               x, y, z, old_label, label, min_ll) ;
        if (label == 53)
        {
          if (DIAG_VERBOSE_ON && old_label != label)
            dump_gcan(gcan, stderr, 0) ;
          DiagBreak() ;
        }
      }
#endif
      
    }
    if (nchanged > 10000)
    {
      ll = gcaGibbsImageLogLikelihood(gca, mri_dst, mri_inputs, lta) ;
      ll /= (double)(width*depth*height) ;
      if (!FZERO(lcma))
        printf("pass %d: %d changed. image ll: %2.3f (CMA=%2.3f)\n", 
               iter+1, nchanged, ll, lcma) ;
      else
        printf("pass %d: %d changed. image ll: %2.3f\n", iter+1, nchanged, ll);
    }
    else
      printf("pass %d: %d changed.\n", iter+1, nchanged) ;
    MRIdilate(mri_changed, mri_changed) ;

  } while (nchanged > 0 && iter++ < max_iter) ;

  free(x_indices) ; free(y_indices) ; free(z_indices) ;
  MRIfree(&mri_changed) ;

  return(mri_dst) ;
}

char *gca_write_fname = NULL ;
int gca_write_iterations = 0 ;

#define MAX_PRIOR_FACTOR 0.5
#define MIN_PRIOR_FACTOR 0.2

MRI  *
GCAreclassifyUsingGibbsPriors(MRI *mri_inputs, GCA *gca, MRI *mri_dst,LTA *lta,
                            int max_iter, MRI *mri_fixed)
{
  int      x, y, z, width, height, depth, label, val, iter,
           xn, yn, zn, n, nchanged,
           index, nindices, old_label, fixed ;
  short    *x_indices, *y_indices, *z_indices ;
  GCA_NODE *gcan ;
  double   ll, lcma = 0.0, old_ll, new_ll, min_ll ;
  MRI      *mri_changed, *mri_probs /*, *mri_zero */ ;

  fixed = (mri_fixed != NULL) ;

#if 0
  GCAannealUnlikelyVoxels(mri_inputs, gca, mri_dst, lta, max_iter*100,
                          mri_fixed) ;
  {
    char  fname[STRLEN], *cp ;

    strcpy(fname, mri_inputs->fname) ;
    cp = strrchr(fname, '/') ;
    if (cp)
    {
      strcpy(cp+1, "anneal") ;
      fprintf(stderr, "writing results of annealing to %s...\n", fname) ;
      MRIwrite(mri_dst, fname) ;
    }
  }
#endif

  nindices = mri_dst->width * mri_dst->height * mri_dst->depth ;
  x_indices = (short *)calloc(nindices, sizeof(short)) ;
  y_indices = (short *)calloc(nindices, sizeof(short)) ;
  z_indices = (short *)calloc(nindices, sizeof(short)) ;
  if (!x_indices || !y_indices || !z_indices)
    ErrorExit(ERROR_NOMEMORY, "GCAreclassifyUsingGibbsPriors: "
              "could not allocate index set") ;


#if 0
  mri_zero = MRIclone(mri_inputs, NULL) ;
#endif
  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_inputs, NULL) ;
    if (!mri_dst)
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst") ;
    MRIcopyHeader(mri_inputs, mri_dst) ;
  }

  mri_changed = MRIclone(mri_dst, NULL) ;


  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_inputs->width ; height = mri_inputs->height; 
  depth = mri_inputs->depth ; iter = 0 ;
  for (x = 0 ; x < width ; x++)
    for (y = 0 ; y < height ; y++)
        for (z = 0 ; z < depth ; z++)
          MRIvox(mri_changed,x,y,z) = 1 ;

  old_ll = gcaGibbsImageLogLikelihood(gca, mri_dst, mri_inputs, lta) ;
  old_ll /= (double)(width*depth*height) ;
#if 0
  {
    MRI   *mri_cma ;
    char  fname[STRLEN], *cp ;

    strcpy(fname, mri_inputs->fname) ;
    cp = strrchr(fname, '/') ;
    if (cp)
    {
      strcpy(cp+1, "parc") ;
      mri_cma = MRIread(fname) ;
      if (mri_cma)
      {
        lcma = gcaGibbsImageLogLikelihood(gca, mri_cma, mri_inputs, lta) ;
        lcma /= (double)(width*depth*height) ;
        fprintf(stderr, "image ll: %2.3f (CMA=%2.3f)\n", old_ll, lcma) ;
        MRIfree(&mri_cma) ;
      }
    }
  }
#endif

  PRIOR_FACTOR = MAX_PRIOR_FACTOR ;
  do
  {
    if (iter == 0)
    {
      /*      char  fname[STRLEN], *cp ;*/
      /*      int   nfixed ;*/

      if (gca_write_iterations)
      {
        char fname[STRLEN] ;
        sprintf(fname, "%s%03d.mgh", gca_write_fname, iter) ;
        printf("writing snapshot to %s...\n", fname) ;
        MRIwrite(mri_dst, fname) ;
      }
      mri_probs = GCAlabelProbabilities(mri_inputs, gca, NULL, lta) ;
#if 0
      for (nfixed = x = 0 ; x < width ; x++)
        for (y = 0 ; y < height ; y++)
          for (z = 0 ; z < depth ; z++)
          {
            if (MRIvox(mri_probs, x, y, z) >= nint(.9*255))
            {
              nfixed++ ;
              MRIvox(mri_fixed, x, y, z) = 1 ;
            }
          }
      
      fprintf(stderr, "%d new fixed points added\n", nfixed) ;
#endif
#if 0
      strcpy(fname, mri_inputs->fname) ;
      cp = strrchr(fname, '/') ;
      strcpy(cp+1, "probs") ;
      fprintf(stderr, "writing label probabilities to %s...\n", fname) ;
      MRIwrite(mri_probs, fname) ;
#endif
      MRIorderIndices(mri_probs, x_indices, y_indices, z_indices) ;
      MRIfree(&mri_probs) ;
    }
    else
      MRIcomputeVoxelPermutation(mri_inputs, x_indices, y_indices,
                                 z_indices) ;
      
    nchanged = 0 ;
    for (index = 0 ; index < nindices ; index++)
    {
      x = x_indices[index] ; y = y_indices[index] ; z = z_indices[index] ;
      if (x == 139 && y == 103 && z == 139)  /* wm should be pallidum */
        DiagBreak() ;
      if (x == 92 && y == 115 && z == 117)
        DiagBreak() ;
      if (mri_fixed && MRIvox(mri_fixed, x, y, z))
        continue ;

      if (MRIvox(mri_changed, x, y, z) == 0)
        continue ;

      if (x == 156 && y == 124 && z == 135)  /* wm should be amygdala */
        DiagBreak() ; 
          
      val = MRIvox(mri_inputs, x, y, z) ;
          
      /* find the node associated with this coordinate and classify */
      GCAsourceVoxelToNode(gca, mri_inputs, lta, x, y, z, &xn, &yn, &zn) ;
      gcan = &gca->nodes[xn][yn][zn] ;
      if (gcan->nlabels == 1)
        continue ;

      label = old_label = MRIvox(mri_dst, x, y, z) ;
      min_ll = gcaNbhdGibbsLogLikelihood(gca, mri_dst, mri_inputs, x, y,z,lta,
                                         PRIOR_FACTOR);

#if 0
      if (min_ll < BIG_AND_NEGATIVE/2 && mri_zero)
        MRIvox(mri_zero, x, y, z) = 255 ;
#endif
      for (n = 0 ; n < gcan->nlabels ; n++)
      {
        if (gcan->labels[n] == old_label)
          continue ;
        MRIvox(mri_dst, x, y, z) = gcan->labels[n] ;
        new_ll = 
          gcaNbhdGibbsLogLikelihood(gca, mri_dst, mri_inputs, x, y,z,lta,
                                    PRIOR_FACTOR);
        if (new_ll > min_ll)
        {
          min_ll = new_ll ; label = gcan->labels[n] ;
        }
      }
      if (x == Ggca_x && y == Ggca_y && z == Ggca_z && 
          (label == Ggca_label || old_label == Ggca_label || Ggca_label < 0))
      {
        printf(
         "(%d, %d, %d): old label %s (%d), new label %s (%d) (log(p)=%2.3f)\n",
         x, y, z, cma_label_to_name(old_label), old_label, 
         cma_label_to_name(label), label, min_ll) ;
        if (label == Right_Caudate)
          DiagBreak() ;
      }
#if 1
      if (((old_label == Left_Cerebral_White_Matter) && 
          (label == Left_Cerebral_Cortex)) ||
          ((old_label == Right_Cerebral_White_Matter) && 
           (label == Right_Cerebral_Cortex)))  /* don't let wm become gm */
        label = old_label ;
#endif
        if (label != old_label)

      {
        nchanged++ ;
        MRIvox(mri_changed, x, y, z) = 1 ;
      }
      else
        MRIvox(mri_changed, x, y, z) = 0 ;
      MRIvox(mri_dst, x, y, z) = label ;
      
      
    }
    if (nchanged > 10000 && iter < 2)
    {
      ll = gcaGibbsImageLogLikelihood(gca, mri_dst, mri_inputs, lta) ;
      ll /= (double)(width*depth*height) ;
      if (!FZERO(lcma))
        printf("pass %d: %d changed. image ll: %2.3f (CMA=%2.3f), PF=%2.3f\n", 
               iter+1, nchanged, ll, lcma, PRIOR_FACTOR) ;
      else
        printf("pass %d: %d changed. image ll: %2.3f, PF=%2.3f\n", 
               iter+1, nchanged, ll, PRIOR_FACTOR) ;
    }
    else
      printf("pass %d: %d changed.\n", iter+1, nchanged) ;
    MRIdilate(mri_changed, mri_changed) ;

#if 0
    if (!iter && DIAG_VERBOSE_ON)
    {
      char  fname[STRLEN], *cp ;
      /*      int   nvox ;*/
      
      strcpy(fname, mri_inputs->fname) ;
      cp = strrchr(fname, '/') ;
      if (cp)
      {
        strcpy(cp+1, "zero") ;
        nvox = MRIvoxelsInLabel(mri_zero, 255) ;
        fprintf(stderr, "writing %d low probability points to %s...\n", 
                nvox, fname) ;
        MRIwrite(mri_zero, fname) ;
        MRIfree(&mri_zero) ;
      }
    }
#endif
#define MIN_CHANGED 2000
    if (nchanged <= MIN_CHANGED)
    {
      for (x = 0 ; x < width ; x++)
        for (y = 0 ; y < height ; y++)
          for (z = 0 ; z < depth ; z++)
            MRIvox(mri_changed,x,y,z) = 1 ;
      if (fixed)
      {
        printf("removing fixed flag...\n") ;
        if (mri_fixed)
          MRIclear(mri_fixed) ;
        fixed = 0 ;
      }
      else
      {
        PRIOR_FACTOR /= 2 ;
        fprintf(stderr, "setting PRIOR_FACTOR to %2.4f\n", PRIOR_FACTOR) ;
      }
      if (gca_write_iterations < 0)
      {
        char fname[STRLEN] ;
        static int fno = 0 ;

        sprintf(fname, "%s%03d.mgh", gca_write_fname, fno+1) ; fno++ ;
        printf("writing snapshot to %s...\n", fname) ;
        MRIwrite(mri_dst, fname) ;
      }
    }
    if ((gca_write_iterations > 0) && !(iter % gca_write_iterations))
    {
      char fname[STRLEN] ;
      sprintf(fname, "%s%03d.mgh", gca_write_fname, iter+1) ;
      printf("writing snapshot to %s...\n", fname) ;
      MRIwrite(mri_dst, fname) ;
    }
  } while ((nchanged > MIN_CHANGED || PRIOR_FACTOR > MIN_PRIOR_FACTOR) && 
           (iter++ < max_iter)) ;

#if 0
  {
    char  fname[STRLEN], *cp ;
    int   nzero, n_nonzero ;

    strcpy(fname, mri_inputs->fname) ;
    cp = strrchr(fname, '/') ;
    if (cp)
    {
      strcpy(cp+1, "indices") ;
      mri_probs = GCAcomputeProbabilities(mri_inputs, gca, mri_dst,NULL, lta) ;
      MRIorderIndices(mri_probs, x_indices, y_indices, z_indices) ;
      for (nzero = index = 0 ; index < nindices ; index++, nzero++)
      {
        x = x_indices[index] ; y = y_indices[index] ; z = z_indices[index] ;
        if (MRIvox(mri_probs, x, y, z) != 255)
          break ;
        MRIvox(mri_probs, x, y, z) = 0 ;
      }
      n_nonzero = nindices - nzero ;
      for ( ; index < nindices ; index++)
      {
        x = x_indices[index] ; y = y_indices[index] ; z = z_indices[index] ;
        if (MRIvox(mri_probs, x, y, z) == 255)
          MRIvox(mri_probs, x, y, z) = 0 ;
        else
        {
          MRIvox(mri_probs, x, y, z) = 
            100 * (float)(n_nonzero-index)/n_nonzero ;
        }
      }
      MRIwrite(mri_probs, fname) ;
      MRIfree(&mri_probs) ;
    }
  }
#endif


  free(x_indices) ; free(y_indices) ; free(z_indices) ;
  MRIfree(&mri_changed) ;

  return(mri_dst) ;
}
int
MRIcomputeVoxelPermutation(MRI *mri, short *x_indices, short *y_indices,
                           short *z_indices)
{
  int width, height, depth, tmp, nindices, i, index ;

  width = mri->width, height = mri->height ; depth = mri->depth ;
  nindices = width*height*depth ;

  for (i = 0 ; i < nindices ; i++)
  {
    x_indices[i] = i % width ;
    y_indices[i] = (i/width) % height ;
    z_indices[i] = (i / (width*height)) % depth ;
  }
  for (i = 0 ; i < nindices ; i++)
  {
    index = (int)randomNumber(0.0, (double)(nindices-0.0001)) ;

    tmp = x_indices[index] ; 
    x_indices[index] = x_indices[i] ; x_indices[i] = tmp ;

    tmp = y_indices[index] ; 
    y_indices[index] = y_indices[i] ; y_indices[i] = tmp ;

    tmp = z_indices[index] ; 
    z_indices[index] = z_indices[i] ; z_indices[i] = tmp ;
  }
  return(NO_ERROR) ;
}

#if 0
static int
gcaGibbsSort(GCA *gca, MRI *mri_labels, MRI *mri_inputs,
                           MATRIX *m_L)
{
  int    x, y, z, width, depth, height ;
  double total_log_likelihood, log_likelihood ;
  MRI    *mri_probs ;

  width = mri_labels->width ; height = mri_labels->height ; 
  depth = mri_labels->depth ;
  mri_probs = MRIclone(mri_labels, NULL) ;
  
  for (total_log_likelihood = 0.0, x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        log_likelihood = 
          gcaVoxelGibbsLogLikelihood(gca, mri_labels, mri_inputs, x, y, z,lta,
                                     PRIOR_FACTOR);
        log_likelihood *= 20 ;
        if (log_likelihood > 255)
          log_likelihood = 255 ;
        MRIvox(mri_probs,x,y,z) = log_likelihood ;
        total_log_likelihood += log_likelihood ;
      }
    }
  }
  MRIorderIndices(mri_probs, x_indices, y_indices, z_indices) ;
  MRIfree(&mri_probs) ;
  return(NO_ERROR) ;
}
#endif

static double
gcaGibbsImageLogLikelihood(GCA *gca, MRI *mri_labels, MRI *mri_inputs,
                           LTA *lta)
{
  int    x, y, z, width, depth, height ;
  double total_log_likelihood, log_likelihood ;

  width = mri_labels->width ; height = mri_labels->height ; 
  depth = mri_labels->depth ;
  
  for (total_log_likelihood = 0.0, x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++)
      {
        log_likelihood = 
          gcaVoxelGibbsLogLikelihood(gca, mri_labels, mri_inputs, x, y, z,lta,
                                     PRIOR_FACTOR);
        total_log_likelihood += log_likelihood ;
      }
    }
  }
  return(total_log_likelihood) ;
}
static double
gcaNbhdGibbsLogLikelihood(GCA *gca, MRI *mri_labels, MRI *mri_inputs, int x, 
                      int y, int z, LTA *lta, double gibbs_coef)
{
  double total_log_likelihood/*, log_likelihood*/ ;
  /*  int    i, xnbr, ynbr, znbr ;*/


  total_log_likelihood = 
    gcaVoxelGibbsLogLikelihood(gca, mri_labels, mri_inputs, x, y, z, lta,
                               gibbs_coef) ;

#if 0
  for (i = 0 ; i < GIBBS_NEIGHBORS ; i++)
  {
    xnbr = mri_inputs->xi[x+xnbr_offset[i]] ;
    ynbr = mri_inputs->yi[y+ynbr_offset[i]] ;
    znbr = mri_inputs->zi[z+znbr_offset[i]] ;
    log_likelihood = 
      gcaVoxelGibbsLogLikelihood(gca, mri_labels, mri_inputs, xnbr, ynbr, 
                                 znbr, lta, gibbs_coef) ;
    total_log_likelihood += log_likelihood ;
  }
#endif

  return(total_log_likelihood) ;
}

static double
gcaVoxelGibbsLogLikelihood(GCA *gca, MRI *mri_labels, MRI *mri_inputs, int x, 
                      int y, int z, LTA *lta, double gibbs_coef)
{
  double    log_likelihood, dist, nbr_prior ;
  int       xn, yn, zn, xnbr, ynbr, znbr, nbr_label, label, val,
            i,j, n;
  GCA_NODE *gcan ;
  GC1D     *gc ;


  val = MRIvox(mri_inputs, x, y, z) ;
  label = MRIvox(mri_labels, x, y, z) ;

  /* find the node associated with this coordinate and classify */
  GCAsourceVoxelToNode(gca, mri_inputs, lta, x, y, z, &xn, &yn, &zn) ;
  gcan = &gca->nodes[xn][yn][zn] ;

  for (n = 0 ; n < gcan->nlabels ; n++)
  {
    if (gcan->labels[n] == label)
      break ;
  }
  if (n >= gcan->nlabels)
    return(BIG_AND_NEGATIVE) ;

  gc = &gcan->gcs[n] ;
  
  /* compute 1-d Mahalanobis distance */
  dist = (val-gc->mean) ;
  log_likelihood = 
    -log(sqrt(gc->var)) - .5*(dist*dist/gc->var) ;

  nbr_prior = 0.0 ;
  for (i = 0 ; i < GIBBS_NEIGHBORS ; i++)
  {
    xnbr = mri_labels->xi[x+xnbr_offset[i]] ;
    ynbr = mri_labels->yi[y+ynbr_offset[i]] ;
    znbr = mri_labels->zi[z+znbr_offset[i]] ;
    nbr_label = MRIvox(mri_labels, xnbr, ynbr, znbr) ;
    for (j = 0 ; j < gc->nlabels[i] ; j++)
    {
      if (nbr_label == gc->labels[i][j])
        break ;
    }
    if (j < gc->nlabels[i])
    {
      if (!FZERO(gc->label_priors[i][j]))
        nbr_prior += log(gc->label_priors[i][j]) ;
      else
        nbr_prior += BIG_AND_NEGATIVE ;
    }
    else   /* never occurred - make it unlikely */
    {
      if (x == 92 && y == 115 && z == 117)
        DiagBreak() ;
      nbr_prior += BIG_AND_NEGATIVE ;
    }
  }
  log_likelihood += (gibbs_coef * nbr_prior + log(gc->prior)) ;

  return(log_likelihood) ;
}
static int compare_sort_mri(const void *plp1, const void *plp2);
typedef struct
{
  unsigned char x, y, z, val ;
} SORT_VOXEL ;

static int
MRIorderIndices(MRI *mri, short *x_indices, short *y_indices, short *z_indices)
{
  int         width, height, depth, nindices, index, x, y, z ;
  SORT_VOXEL  *sort_voxels ;

  width = mri->width, height = mri->height ; depth = mri->depth ;
  nindices = width*height*depth ;

  sort_voxels = (SORT_VOXEL *)calloc(nindices, sizeof(SORT_VOXEL)) ;
  if (!sort_voxels)
    ErrorExit(ERROR_NOMEMORY,"MRIorderIndices: could not allocate sort table");

  for (index = x = 0 ; x < width ; x++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (z = 0 ; z < depth ; z++, index++)
      {
        sort_voxels[index].x = x ;
        sort_voxels[index].y = y ;
        sort_voxels[index].z = z ;
        sort_voxels[index].val = MRIvox(mri, x, y, z) ;
      }
    }
  }
  qsort(sort_voxels, nindices, sizeof(SORT_VOXEL), compare_sort_mri) ;

  for (index = 0 ; index < nindices ; index++)
  {
    x_indices[index] = sort_voxels[index].x ;
    y_indices[index] = sort_voxels[index].y ;
    z_indices[index] = sort_voxels[index].z ;
  }
  
  free(sort_voxels) ;
  return(NO_ERROR) ;
}

static int
compare_sort_mri(const void *psv1, const void *psv2)
{
  SORT_VOXEL  *sv1, *sv2 ;

  sv1 = (SORT_VOXEL *)psv1 ;
  sv2 = (SORT_VOXEL *)psv2 ;

  if (sv1->val > sv2->val)
    return(1) ;
  else if (sv1->val < sv2->val)
    return(-1) ;

  return(0) ;
}
MRI *
GCAbuildMostLikelyVolume(GCA *gca, MRI *mri)
{
  int       x,  y, z, xn, yn, zn, width, depth, height, n ;
  GCA_NODE *gcan ;
  double   max_prior, max_mean ;
  GC1D     *gc ;

  width = mri->width ; depth = mri->depth ; height = mri->height ;

  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        GCAvoxelToNode(gca, mri, x, y, z, &xn, &yn, &zn) ;
        gcan = &gca->nodes[xn][yn][zn] ;
        max_prior = gcan->gcs[0].prior ; max_mean = gcan->gcs[0].mean ;
        for (n = 1 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          if (gc->prior > max_prior)
          {
            max_prior = gc->prior ; max_mean = gc->mean ;
          }
        }
        MRIvox(mri, x, y, z) = (int)max_mean ;
      }
    }
  }

  return(mri) ;
}

#if 1
static GC1D *
gcaFindGC(GCA *gca, int x, int y, int z,int label)
{
  int        n ;
  GCA_NODE   *gcan  ;

  gcan = &gca->nodes[x][y][z] ;

  for (n = 0 ; n < gcan->nlabels ; n++)
  {
    if (gcan->labels[n] == label)
      return(&gcan->gcs[n]) ;
  }

  return(NULL) ;
}
#endif

#include "mrisegment.h"
/* each segment must be at least this much of total to be retained */
#define MIN_SEG_PCT  0.05   

static int gcaReclassifySegment(GCA *gca, MRI *mri_inputs, MRI *mri_labels,
                                MRI_SEGMENT *mseg, int old_label, LTA *lta);
static int gcaReclassifyVoxel(GCA *gca, MRI *mri_inputs, MRI *mri_labels,
                              int x, int y, int z, int old_label, LTA *lta);
MRI *
GCAconstrainLabelTopology(GCA *gca, MRI *mri_inputs,MRI *mri_src, MRI *mri_dst,
                          LTA *lta)
{
  int              i, j, nvox /*, x, y, z, width, height, depth*/ ;
  MRI_SEGMENTATION *mriseg ;

  mri_dst = MRIcopy(mri_src, mri_dst) ;

  for (i = 1 ; i < 255 ; i++)
  {
    nvox = MRIvoxelsInLabel(mri_dst, i) ;
    if (!nvox)
      continue ;
    /*    printf("label %03d: %d voxels\n", i, nvox) ;*/
    mriseg = MRIsegment(mri_src, (float)i, (float)i) ;
    /*    printf("\t%d segments:\n", mriseg->nsegments) ;*/
    for (j = 0 ; j < mriseg->nsegments ; j++)
    {
      /*      printf("\t\t%02d: %d voxels", j, mriseg->segments[j].nvoxels) ;*/
      if ((float)mriseg->segments[j].nvoxels / (float)nvox < MIN_SEG_PCT)
      {
        /*        printf(" - reclassifying...") ;*/
        gcaReclassifySegment(gca,mri_inputs,mri_dst, &mriseg->segments[j], i,
                             lta);
      }
      /*      printf("\n") ;*/
    }
    MRIsegmentFree(&mriseg) ;
  }

#if 0
  width = mri_dst->width ; height = mri_dst->height ; 
  depth  = mri_dst->depth ; 
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (x == 144 && y == 118 && z == 127)
          DiagBreak() ;
        if (MRIvox(mri_dst, x, y, z) == LABEL_UNDETERMINED)
          gcaReclassifyVoxel(gca, mri_inputs, mri_dst, 
                             x, y, z, LABEL_UNDETERMINED, lta) ;
      }
    }
  }
#endif

  return(mri_dst) ;
}

static int
gcaReclassifySegment(GCA *gca, MRI *mri_inputs, MRI *mri_labels,
                     MRI_SEGMENT *mseg, int old_label, LTA *lta)
{
  int   i ;

  for (i = 0 ; i < mseg->nvoxels ; i++)
  {
#if 1
    gcaReclassifyVoxel(gca, mri_inputs, mri_labels, 
                       mseg->voxels[i].x, mseg->voxels[i].y, mseg->voxels[i].z,
                       old_label, lta) ;
#else
    MRIvox(mri_labels,mseg->voxels[i].x,mseg->voxels[i].y,mseg->voxels[i].z) =
      LABEL_UNDETERMINED ;
#endif
  }

  return(NO_ERROR) ;
}

static int
gcaReclassifyVoxel(GCA *gca, MRI *mri_inputs, MRI *mri_labels,
                     int x, int y, int z, int old_label, LTA *lta)
{
  int     nbr_labels[255], xi, yi, zi, xk, yk, zk, i, new_label ;
  double  max_p, p ;

  if (x == 144 && y == 118 && z == 127)
    DiagBreak() ;
  memset(nbr_labels, 0, sizeof(nbr_labels)) ;
  for (zk = -1 ; zk <= 1 ; zk++)
  {
    zi = mri_labels->zi[z+zk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      yi = mri_labels->yi[y+yk] ;
      for (xk = -1 ; xk <= 1 ; xk++)
      {
        xi = mri_labels->xi[x+xk] ;
        nbr_labels[MRIvox(mri_labels, xi, yi, zi)]++ ;
      }
    }
  }
  new_label = 0 ; max_p = 10*BIG_AND_NEGATIVE ;
  nbr_labels[old_label] = 0 ;
  for (i = 0 ; i <= 255 ; i++)
  {
    if (nbr_labels[i] > 0)
    {
      MRIvox(mri_labels, x, y, z) = i ;
      p = gcaVoxelGibbsLogLikelihood(gca, mri_labels, mri_inputs, x, y, z,lta,
                                     PRIOR_FACTOR);
      if (p > max_p)
      {
        max_p = p ;
        new_label = i ;
      }
    }
  }

  MRIvox(mri_labels, x, y, z) = new_label ;
  return(NO_ERROR) ;
}
#define MAX_VENTRICLE_ITERATIONS  10
MRI *
GCAexpandVentricle(GCA *gca, MRI *mri_inputs, MRI *mri_src,
                   MRI *mri_dst, LTA *lta, int target_label)
{
  int      nchanged, x, y, z, width, height, depth, xn, yn, zn, xi, yi, zi,
           xk, yk, zk, label, val, total_changed, i ;
  GCA_NODE *gcan ;
  float    v_mean, v_var, label_mean, label_var, pv, plabel, dist ;
  MRI      *mri_tmp ;

  /* compute label mean and variance */

  v_mean = gcaComputeLabelStats(gca, target_label, &v_var) ;
  printf("ventricle intensity = %2.1f +- %2.1f\n", v_mean, sqrt(v_var)) ;

  if (mri_src != mri_dst)
    mri_dst = MRIcopy(mri_src, mri_dst) ;

  mri_tmp = MRIcopy(mri_dst, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;

  i = total_changed = 0 ;
  do
  {
    nchanged = 0 ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == 140 && y == 79 && z == 136)
            DiagBreak() ; /* wm should be ventricle (1482)*/
          if (x == 140 && y == 111 && z == 139)
            DiagBreak() ;  /* should be wm */
          if (x == 138 && y == 103 && z == 139)
            DiagBreak() ;  /* should be pallidum */
          
          label = MRIvox(mri_dst, x, y, z) ;
          
          if (label == target_label)
          {
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              for (yk = -1 ; yk <= 1 ; yk++)
              {
                for (xk = -1 ; xk <= 1 ; xk++)
                {
                  if (fabs(xk) + fabs(yk) + fabs(zk) > 1)
                    continue ;
                  xi = mri_src->xi[x+xk] ;
                  yi = mri_src->yi[y+yk] ;
                  zi = mri_src->zi[z+zk] ;
                  label = MRIvox(mri_dst, xi, yi, zi) ;
                  if (label != target_label) /* should it be changed? */
                  {
                    GCAsourceVoxelToNode(gca, mri_dst, lta, xi, yi, zi, 
                                         &xn, &yn, &zn) ;
                    gcan = &gca->nodes[xn][yn][zn] ;
                    label_mean = getLabelMean(gcan, label, &label_var) ;
            
                    if (xi == 140 && yi == 79 && zi == 136)
                      DiagBreak() ;   /* v should be wm */
                    if (xi == 145 && yi == 130 && zi == 87)
                      DiagBreak() ;   
                    val = MRIvox(mri_inputs, xi, yi, zi) ;
                    dist = val - label_mean ;
                    if (!FZERO(label_var))
                      plabel = 1/sqrt(label_var * 2 * M_PI)*
                        exp(-dist*dist/label_var);
                    else
                    {
                      if (!FZERO(dist))
                        plabel = 0.0f ;
                      else
                        plabel = 1.0f ;
                    }
                    dist = val - v_mean ;
                    if (!FZERO(v_var))
                      pv = 1/sqrt(v_var * 2 * M_PI) * exp(-dist*dist/v_var);
                    else
                    {
                      if (!FZERO(dist))
                        pv = 0.0f ;
                      else
                        pv = 1.0f ;
                    }
                    if (pv > 5*plabel)
                    {
                      if (xi == 140 && yi == 79 && zi == 136)
                        DiagBreak() ;   /* v should be wm */
                      nchanged++ ;
                      MRIvox(mri_tmp, xi, yi, zi) = target_label ;
                    }
                  }
                }
              }
            }
          } 
        }
      }
    }
    MRIcopy(mri_tmp, mri_dst) ;
    total_changed += nchanged ;
    if (++i > MAX_VENTRICLE_ITERATIONS)
      break ;
  } while (nchanged > 0) ;

  MRIfree(&mri_tmp) ;
  printf("%d labels changed to %d...\n", total_changed, target_label) ;
  return(mri_dst) ;
}
#define MAX_CORTICAL_ITERATIONS  10
MRI *
GCAexpandCortex(GCA *gca, MRI *mri_inputs, MRI *mri_src,
                MRI *mri_dst, LTA *lta)
{
  int      nchanged, x, y, z, width, height, depth, xn, yn, zn, xi, yi, zi,
           xk, yk, zk, label, val, total_changed, i, wm_nbr, gray_nbr ;
  GCA_NODE *gcan ;
  float    wm_mean, wm_var, gray_mean, gray_var, 
           label_mean, label_var, ldist, wdist, gdist ;
  MRI      *mri_tmp ;

  /* compute label mean and variance */

  wm_mean = gcaComputeLabelStats(gca, Left_Cerebral_White_Matter, &wm_var) ;
  gray_mean = gcaComputeLabelStats(gca, Left_Cerebral_Cortex, &gray_var) ;
  printf("cortex mean - gray %2.1f +- %2.1f, white %2.0f +- %2.0f\n", 
         gray_mean, sqrt(gray_var), wm_mean, sqrt(wm_var)) ;

  if (mri_src != mri_dst)
    mri_dst = MRIcopy(mri_src, mri_dst) ;

  mri_tmp = MRIcopy(mri_dst, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;

  i = total_changed = 0 ;
  do
  {
    nchanged = 0 ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == 140 && y == 79 && z == 136)
            DiagBreak() ; /* wm should be ventricle (1482)*/
          if (x == 140 && y == 111 && z == 139)
            DiagBreak() ;  /* should be wm */
          if (x == 78 && y == 79 && z == 93) 
            DiagBreak() ;  /* unknown, should be wm */
          
          label = MRIvox(mri_dst, x, y, z) ;
          
          if (label == Unknown ||
              label == Left_Cerebral_Cortex ||
              label == Right_Cerebral_Cortex ||
              label == Left_Cerebral_White_Matter ||
              label == Right_Cerebral_White_Matter
              )
          {
            wm_nbr = gray_nbr = 0 ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              for (yk = -1 ; yk <= 1 ; yk++)
              {
                for (xk = -1 ; xk <= 1 ; xk++)
                {
                  if (fabs(xk) + fabs(yk) + fabs(zk) > 1)
                    continue ;
                  xi = mri_src->xi[x+xk] ;
                  yi = mri_src->yi[y+yk] ;
                  zi = mri_src->zi[z+zk] ;
                  label = MRIvox(mri_dst, xi, yi, zi) ;
                  if (label ==  Left_Cerebral_Cortex ||
                      label ==  Right_Cerebral_Cortex)
                    gray_nbr = label ;
                  else
                  if (label ==  Left_Cerebral_White_Matter ||
                      label ==  Right_Cerebral_White_Matter)
                    wm_nbr = label ;
                }
              }
            }
            if (!wm_nbr && !gray_nbr)
              continue ;

            if (x == 78 && y == 79 && z == 93) 
              DiagBreak() ;  /* unknown, should be wm */

            val = MRIvox(mri_inputs, x, y, z) ;
            if (wm_nbr)
              wdist = fabs(val - wm_mean) ;
            else
              wdist = 100000 ;
            if (gray_nbr)
              gdist = fabs(val - gray_mean) ;
            else
              gdist = 100000 ;

            if (wdist > sqrt(wm_var))
              wdist = 10000 ; /* hack - don't label unlikely white */
            if (gdist > sqrt(gray_var))
              gdist = 10000 ;  /* hack - don't label unlikely gray */
            
            GCAsourceVoxelToNode(gca, mri_dst, lta, x, y, z,  &xn, &yn, &zn) ;
            gcan = &gca->nodes[xn][yn][zn] ;
            label_mean = getLabelMean(gcan, label, &label_var) ;
            
            if (x == 140 && y == 79 && z == 136)
              DiagBreak() ;   /* v should be wm */
            ldist = fabs(val - label_mean) ;

            ldist *= .75 ;  /* bias towards retaining label */

            if (wdist < gdist)   /* might change to wm */
            {
              if (wdist < ldist)
              {
                if (x == 140 && y == 79 && z == 136)
                  DiagBreak() ;   /* v should be wm */
                nchanged++ ;
                MRIvox(mri_tmp, x, y, z) = wm_nbr ;
              }
            }
            else                 /* might change to gm */
            {
              if (gdist < ldist)
              {
                if (x == 140 && y == 79 && z == 136)
                  DiagBreak() ;   /* v should be wm */
                nchanged++ ;
                MRIvox(mri_tmp, x, y, z) = gray_nbr ;
              }
            }
          } 
        }
      }
    }
    MRIcopy(mri_tmp, mri_dst) ;
    total_changed += nchanged ;
    if (++i > MAX_CORTICAL_ITERATIONS)
      break ;
  } while (nchanged > 0) ;

  MRIfree(&mri_tmp) ;
  printf("%d labels changed to cortex...\n", total_changed) ;
  return(mri_dst) ;
}
MRI   *
GCAexpandLabelIntoWM(GCA *gca, MRI *mri_inputs, MRI *mri_src,
               MRI *mri_dst, LTA *lta, MRI *mri_fixed, int target_label)
{
  int      nchanged, x, y, z, width, height, depth, xn, yn, zn, xi, yi, zi,
           xk, yk, zk, nbr_label, n, label, val, total_changed, i ;
  GCA_NODE *gcan, *gcan_nbr ;
  float    label_mean, wm_mean, prior ;
  MRI      *mri_tmp ;

  if (mri_src != mri_dst)
    mri_dst = MRIcopy(mri_src, mri_dst) ;

  mri_tmp = MRIcopy(mri_dst, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;

  i = total_changed = 0 ;
  do
  {
    nchanged = 0 ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == 140 && y == 111 && z == 139)
            DiagBreak() ;  /* should be wm */
          if (x == 138 && y == 103 && z == 139)
            DiagBreak() ;  /* should be pallidum */
          
          label = MRIvox(mri_dst, x, y, z) ;
          GCAsourceVoxelToNode(gca, mri_dst, lta, x, y, z, &xn, &yn, &zn) ;
          gcan = &gca->nodes[xn][yn][zn] ;
          
          if (label == target_label)
          {
            label_mean = wm_mean = 0.0 ;
            for (n = 0 ; n < gcan->nlabels ; n++)
            {
              if (gcan->labels[n] == target_label)
                label_mean = gcan->gcs[n].mean ;
              else if 
                ((gcan->labels[n] == Left_Cerebral_White_Matter) ||
                 (gcan->labels[n] == Right_Cerebral_White_Matter))
                wm_mean = gcan->gcs[n].mean ;
            }
            
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              for (yk = -1 ; yk <= 1 ; yk++)
              {
                for (xk = -1 ; xk <= 1 ; xk++)
                {
                  if (fabs(xk) + fabs(yk) + fabs(zk) > 1)
                    continue ;
                  xi = mri_src->xi[x+xk] ;
                  yi = mri_src->yi[y+yk] ;
                  zi = mri_src->zi[z+zk] ;
                  nbr_label = MRIvox(mri_dst, xi, yi, zi) ;
                  if ((nbr_label == Right_Cerebral_White_Matter) ||
                      (nbr_label == Left_Cerebral_White_Matter))
                  {
                    if (xi == 140 && yi == 111 && zi == 139)
                      DiagBreak() ;  /* should be wm */
                    if (xi == 138 && yi == 103 && zi == 139)
                      DiagBreak() ;  /* should be pallidum */
                    val = MRIvox(mri_inputs, xi, yi, zi) ;
                    if (fabs(val-wm_mean) > fabs(val-label_mean))
                    {
          
                      GCAsourceVoxelToNode(gca, mri_dst, lta, xi, yi, zi, 
                                           &xn, &yn, &zn) ;
                      gcan_nbr = &gca->nodes[xn][yn][zn] ;
                      for (prior = 0.0f, n = 0 ; n < gcan_nbr->nlabels ; n++)
                      {
                        if (gcan_nbr->labels[n] == target_label)
                        {
                          prior = gcan_nbr->gcs[n].prior ;
                          break ;
                        }
                      }
#define PRIOR_THRESH 0.01
                      if (prior >= PRIOR_THRESH)  /* target is possible */
                      {
                        nchanged++ ;
                        MRIvox(mri_tmp, xi, yi, zi) = target_label ;
                        /*                        MRIvox(mri_fixed, xi, yi, zi) = 0 ;*/
                      }
                    }
                  }
                }
              }
            }
          } 
        }
      }
    }
    MRIcopy(mri_tmp, mri_dst) ;
    total_changed += nchanged ;
    if (++i >= 1)
      break ;
  } while (nchanged > 0) ;

  MRIfree(&mri_tmp) ;
  printf("%d labels changed to %d...\n", total_changed, target_label) ;
  return(mri_dst) ;
}


int
GCArankSamples(GCA *gca, GCA_SAMPLE *gcas, int nsamples, int *ordered_indices)
{
  LABEL_PROB  *label_probs ;
  int         i ;

  label_probs = (LABEL_PROB *)calloc(nsamples, sizeof(LABEL_PROB)) ;
  for (i = 0 ; i < nsamples ; i++)
  {
    label_probs[i].label = i ; label_probs[i].prob =  gcas[i].log_p ;
  }

  /* now sort the samples by probability */
  qsort(label_probs, nsamples, sizeof(LABEL_PROB), 
              compare_sort_probabilities) ;

  for (i = 0 ; i < nsamples ; i++)
  {
    ordered_indices[i] = label_probs[nsamples-(i+1)].label ;
  }

  free(label_probs) ;
  return(NO_ERROR) ;
}

#include "mrinorm.h"
MRI *
GCAnormalizeSamples(MRI *mri_in, GCA *gca, GCA_SAMPLE *gcas, int nsamples, 
                    LTA *lta)
{
  MRI    *mri_dst, *mri_ctrl, *mri_bias ;
  int    xv, yv, zv, n, x, y, z, width, height, depth, out_val ;
  MATRIX *m_L_inv ;
  VECTOR *v_input, *v_canon ;
  float   bias ;

  m_L_inv = MatrixInverse(lta->xforms[0].m_L, NULL) ;
  if (!m_L_inv)
  {
    MatrixPrint(stderr, lta->xforms[0].m_L) ;
    ErrorExit(ERROR_BADPARM, 
              "GCAcomputeLogSampleProbability: matrix is not invertible") ;
  }
  mri_dst = MRIclone(mri_in, NULL) ;
  mri_ctrl = MRIclone(mri_in, NULL) ;
  v_input = VectorAlloc(4, MATRIX_REAL) ;
  v_canon = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v_input, 4, 1) = 1.0 ; *MATRIX_RELT(v_canon, 4, 1) = 1.0 ;
  mri_bias = MRIalloc(mri_in->width,mri_in->height,mri_in->depth,MRI_SHORT);
  if (!mri_bias)
    ErrorExit(ERROR_NOMEMORY, 
              "GCAnormalize: could not allocate (%d,%d,%d,2) bias image",
              mri_in->width,mri_in->height,mri_in->depth) ;
              

  width = mri_in->width ; height = mri_in->height ; depth = mri_in->depth ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        MRISvox(mri_bias, x,y,z) = 1000 ;
      }
    }
  }
  for (n = 0 ; n < nsamples ; n++)
  {
    GCAnodeToVoxel(gca, mri_dst, gcas[n].xn, gcas[n].yn, gcas[n].zn, 
                   &xv, &yv, &zv) ;
    V3_X(v_canon) = (float)xv ; V3_Y(v_canon) = (float)yv ;
    V3_Z(v_canon) = (float)zv ;
    MatrixMultiply(m_L_inv, v_canon, v_input) ;
    xv = nint(V3_X(v_input)); yv = nint(V3_Y(v_input)); zv=nint(V3_Z(v_input));
    if (xv < 0) xv = 0 ;
    if (xv >= mri_in->width) xv = mri_in->width-1 ;
    if (yv < 0) yv = 0 ;
    if (yv >= mri_in->height) yv = mri_in->height-1 ;
    if (zv < 0) zv = 0 ;
    if (zv >= mri_in->depth) zv = mri_in->depth-1 ;

    if (xv == 181 && yv == 146 && zv == 128)
      DiagBreak() ;
    if (xv == 181 && yv == 135 && zv == 121)
      DiagBreak() ;
    if (gcas[n].label == 29 || gcas[n].label == 61)
    {
      gcas[n].label = 0 ;
      DiagBreak() ;
    }
    if (gcas[n].label > 0)
    {
      MRIvox(mri_ctrl, xv, yv, zv) = CONTROL_MARKED ;
      bias = 1000.0f*((float)gcas[n].mean/MRIvox(mri_in, xv, yv, zv)) ;
      if (bias < 900)
        bias = 900 ;
      if (bias > 1100)
        bias = 1100 ;
      MRISvox(mri_bias, xv, yv, zv) = (short)nint(bias) ;
    }
    else
      MRIvox(mri_ctrl, xv, yv, zv) = CONTROL_NONE ;
  }

  MRIbuildVoronoiDiagram(mri_bias, mri_ctrl, mri_bias) ;
  /*  MRIwrite(mri_bias, "bias.mgh") ;*/
#if 0
  {
    MRI *mri_kernel, *mri_smooth ;

    mri_kernel = MRIgaussian1d(2.0f, 100) ;
    mri_smooth = MRIconvolveGaussian(mri_bias, NULL, mri_kernel) ;
    MRIfree(&mri_bias) ; mri_bias = mri_smooth ; MRIfree(&mri_kernel) ;
  }
#else
  MRIsoapBubble(mri_bias, mri_ctrl, mri_bias, 10) ;
#endif
  /*  MRIwrite(mri_bias, "smooth_bias.mgh") ;*/


  width = mri_in->width ; height = mri_in->height ; depth = mri_in->depth ;
  for (z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        bias = MRISvox(mri_bias, x, y, z)/1000.0f ;
        out_val = nint((float)MRIvox(mri_in, x, y, z)*bias) ;
        if (out_val < 0)
          out_val = 0 ;
        else if (out_val > 255)
          out_val = 255 ;
        MRIvox(mri_dst, x, y, z) = (BUFTYPE)out_val ;
      }
    }
  }

  MatrixFree(&m_L_inv) ; MatrixFree(&v_canon) ; MatrixFree(&v_input) ;
  MRIfree(&mri_bias) ; MRIfree(&mri_ctrl) ;
  return(mri_dst) ;
}
int
GCAnodeToSourceVoxel(GCA *gca, MRI *mri, MATRIX *m_L, int xn, int yn, int zn, 
                    int *pxv, int *pyv, int *pzv)
{
  static VECTOR  *v_input, *v_canon = NULL ;
  int   width, height, depth, xt, yt, zt, xv, yv, zv ;
  static MATRIX  *m_L_inv = NULL ;

  m_L_inv = MatrixInverse(m_L, m_L_inv) ;
  if (!m_L_inv)
  {
    MatrixPrint(stderr, m_L) ;
    ErrorExit(ERROR_BADPARM, 
              "GCAcomputeLogSampleProbability: matrix is not invertible") ;
  }

  width = mri->width ; height = mri->height ;  depth = mri->depth ;
  if (!v_canon)
  {
    v_input = VectorAlloc(4, MATRIX_REAL) ;
    v_canon = VectorAlloc(4, MATRIX_REAL) ;
    *MATRIX_RELT(v_input, 4, 1) = 1.0 ; *MATRIX_RELT(v_canon, 4, 1) = 1.0 ;
  }
  GCAnodeToVoxel(gca, mri, xn, yn, zn, &xt, &yt, &zt) ;
  V3_X(v_canon) = (float)xt*mri->xsize ; 
  V3_Y(v_canon) = (float)yt*mri->ysize ;
  V3_Z(v_canon) = (float)zt*mri->zsize ;
  MatrixMultiply(m_L_inv, v_canon, v_input) ;
  xv = nint(V3_X(v_input)/mri->xsize); 
  yv = nint(V3_Y(v_input)/mri->ysize); 
  zv = nint(V3_Z(v_input)/mri->zsize);
  if (xv < 0) xv = 0 ;
  if (yv < 0) yv = 0 ;
  if (zv < 0) zv = 0 ;
  if (xv >= width)  xv = width-1 ;
  if (yv >= height) yv = height-1 ;
  if (zv >= depth)  zv = depth-1 ;

  *pxv = xv ; *pyv = yv ; *pzv = zv ;
  return(NO_ERROR) ;
}

float 
GCAlabelProbability(MRI *mri_src, GCA *gca, LTA *lta,
                          int x, int y, int z, int label)
{
  int      xn, yn, zn, n, val ;
  GCA_NODE *gcan ;
  float    label_mean, dist, label_var, plabel;

  if (x == 91 && y == 85 && z == 99) /* wm should be ventricle (1337)*/
    DiagBreak() ; 
  if (x == 140 && y == 111 && z == 139)
    DiagBreak() ;  /* should be wm */
  if (x == 138 && y == 103 && z == 139)
    DiagBreak() ;  /* should be pallidum */
  
  GCAsourceVoxelToNode(gca, mri_src, lta, x, y, z, &xn, &yn, &zn) ;
  gcan = &gca->nodes[xn][yn][zn] ;

  label_var = label_mean = 0.0 ;
  for (n = 0 ; n < gcan->nlabels ; n++)
  {
    if (gcan->labels[n] == label)
    {
      label_mean = gcan->gcs[n].mean ;
      label_var = gcan->gcs[n].var ;
    }
  }
    
  val = MRIvox(mri_src, x, y, z) ;
  dist = val - label_mean ;
  if (!FZERO(label_var))
    plabel = 
      1/sqrt(label_var * 2 * M_PI) * exp(-dist*dist/label_var);
  else
  {
    if (!FZERO(dist))
      plabel = 0.0f ;
    else
      plabel = 1.0f ;
  }

  return(plabel) ;
}

static GCA_NODE *
findSourceGCAN(GCA *gca, MRI *mri_src, LTA *lta, int x, int y, int z)
{
  int      xn, yn, zn ;
  GCA_NODE *gcan ;

  GCAsourceVoxelToNode(gca, mri_src, lta, x, y, z, &xn, &yn, &zn) ;
  gcan = &gca->nodes[xn][yn][zn] ;
  return(gcan) ;
}
static float
getLabelMean(GCA_NODE *gcan, int label, float *pvar)
{
  int   n ;
  float mean = -1.0f ;

  if (*pvar)
    *pvar = 0.0f ;
  for (n = 0 ; n < gcan->nlabels ; n++)
  {
    if (gcan->labels[n] == label)
    {
      mean = gcan->gcs[n].mean ;
      if (pvar)
        *pvar = gcan->gcs[n].var ;
      break ;
    }
  }
  return(mean) ;
}

MRI   *
GCAmaxLikelihoodBorders(GCA *gca, MRI *mri_inputs, MRI *mri_src,
                        MRI *mri_dst, LTA *lta, int max_iter, float min_ratio)
{
  int      nchanged, x, y, z, width, height, depth, label, total_changed, i ;
  MRI      *mri_tmp ;

  if (mri_src != mri_dst)
    mri_dst = MRIcopy(mri_src, mri_dst) ;

  mri_tmp = MRIcopy(mri_dst, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;

  for (total_changed = i = 0 ; i < max_iter ; i++)
  {
    nchanged = 0 ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == 99 && y == 129 && z == 127)
            DiagBreak() ;  /* gray should be wm */
          if (x == 98 && y == 124 && z == 127)
            DiagBreak() ;  /* wm should be hippo */

          if (borderVoxel(mri_dst,x,y,z))
          {
            label = GCAmaxLikelihoodBorderLabel(gca, mri_inputs, mri_dst,lta,
                                                x, y, z, min_ratio) ;
            if (x == Ggca_x && y == Ggca_y && z == Ggca_z && 
                (label == Ggca_label || MRIvox(mri_tmp,x,y,z) == Ggca_label || 
                 Ggca_label < 0))
              printf(
                   "MLE (%d, %d, %d): old label %s (%d), new label %s (%d)\n",
                   x, y, z, cma_label_to_name(MRIvox(mri_tmp,x,y,z)),
                   MRIvox(mri_tmp,x,y,z), cma_label_to_name(label),
                   label) ;
            if (label != MRIvox(mri_dst, x, y, z))
            {
              nchanged++ ;
              MRIvox(mri_tmp, x, y, z) = label ;
            }
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_dst) ;
    total_changed += nchanged ;
    if (!nchanged)
      break ;
  } 

  MRIfree(&mri_tmp) ;
  printf("%d border labels changed to MLE ...\n", total_changed) ;
  return(mri_dst) ;
}
static int
borderVoxel(MRI *mri, int x, int y, int z)
{
  int   xi, yi, zi, xk, yk, zk, label ;

  label = MRIvox(mri, x, y, z) ;

  for (xk = -1 ; xk <= 1 ; xk++)
  {
    xi = mri->xi[x+xk] ;
    for (yk = -1 ; yk <= 1 ; yk++)
    {
      for (zk = -1 ; zk <= 1 ; zk++)
      {
        if (abs(xk)+abs(yk)+abs(zk) != 1)
          continue ;
        yi = mri->yi[y+yk] ;
        zi = mri->zi[z+zk] ;
        if (MRIvox(mri, xi, yi, zi) != label)
          return(1) ;
      }
    }
  }
  return(0) ;
}

static int
GCAmaxLikelihoodBorderLabel(GCA *gca, MRI *mri_inputs, MRI *mri_labels, 
                       LTA *lta, int x, int y, int z, float min_ratio)
{
  float    mean, var, p, max_p, dist, log_likelihood ;
  int      label, i, xi, yi, zi, best_label, val, orig_label ;
  GCA_NODE *gcan ;

  if (x == 99 && y == 129 && z == 127)
    DiagBreak() ;  /* gray should be wm */
  if (x == 98 && y == 124 && z == 127)
    DiagBreak() ;  /* wm should be hippo */
  
  val = MRIvox(mri_inputs, x, y, z) ;
  gcan = findSourceGCAN(gca, mri_inputs, lta, x, y, z) ;
  orig_label = best_label = MRIvox(mri_labels, x, y, z) ;
  mean = getLabelMean(gcan, best_label, &var) ;
  dist = val - mean ;
  max_p = 1 / sqrt(var * 2 * M_PI) * exp(-dist*dist/var) ;

  for (i = 0 ; i < GIBBS_NEIGHBORS ; i++)
  {
    xi = mri_inputs->xi[x+xnbr_offset[i]] ;
    yi = mri_inputs->yi[y+ynbr_offset[i]] ;
    zi = mri_inputs->zi[z+znbr_offset[i]] ;
    gcan = findSourceGCAN(gca, mri_inputs, lta, xi, yi, zi) ;
    label = MRIvox(mri_labels, xi, yi, zi) ;
    mean = getLabelMean(gcan, label, &var) ;
    dist = val - mean ;
    p = 1 / sqrt(var * 2 * M_PI) * exp(-dist*dist/var) ;
    if ((best_label == orig_label && p > min_ratio*max_p) ||
        (best_label != orig_label && p > max_p))
    {
      max_p = p ;
      best_label = label ;
    }
  }

  /* test to make sure that it is not an impossible Gibbs configuration */
  if (best_label != MRIvox(mri_labels, x, y, z))
  {
    label = MRIvox(mri_labels, x, y, z) ;
    MRIvox(mri_labels, x, y, z) = best_label ;  /* test potential new label */
    log_likelihood = 
      gcaNbhdGibbsLogLikelihood(gca, mri_labels, mri_inputs, x, y,z,
                                lta, 1) ;
    if (log_likelihood < BIG_AND_NEGATIVE/2) 
      best_label = label ;  /* impossible - revert to old */
    MRIvox(mri_labels, x, y, z) = label ;  /* revert back to old label */
  }
  return(best_label) ;
}


static float
gcaComputeLabelStats(GCA *gca, int target_label, float *pvar)
{
  int      x, y, z, n ;
  double   mean, var, dof, total_dof ;
  GC1D     *gc ;
  GCA_NODE *gcan ;

  mean = var = total_dof = 0.0 ;
  for (x = 0 ; x < gca->width ; x++)
  {
    for (y = 0 ; y < gca->height ; y++)
    {
      for (z = 0 ; z < gca->depth ; z++)
      {
        gcan = &gca->nodes[x][y][z] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          if (gcan->labels[n] == target_label)
          {
            gc = &gcan->gcs[n] ;
            dof = gc->prior * gcan->total_training ;
            mean += dof*gc->mean ;
            var += dof*gc->var ;
            total_dof += dof ;
          }
        }
      }
    }
  }

  if (total_dof > 0.0)
  {
    mean /= total_dof ;
    var /= total_dof ;
  }
  if (pvar)
    *pvar = var ;
  return(mean) ;
}
int
GCAaccumulateTissueStatistics(GCA *gca, MRI *mri_T1,MRI *mri_PD,
                              MRI *mri_labeled, LTA *lta)
{
  int              x, y, z, n, label, biggest_label, T1, PD, xp, yp, zp ;
  GCA_NODE         *gcan ;
  GCA_TISSUE_PARMS *gca_tp ;
  VECTOR           *v_parc, *v_T1 ;
  static VECTOR *v_ras_cor = NULL, *v_ras_flash ;
  MATRIX           *m_L ;
  
  v_parc = VectorAlloc(4, MATRIX_REAL) ;
  v_T1 = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v_parc, 4, 1) = 1.0 ;
  *MATRIX_RELT(v_T1, 4, 1) = 1.0 ;

  /* first build a list of all labels that exist */
  for (biggest_label = x = 0 ; x < gca->width ; x++)
  {
    for (y = 0 ; y < gca->height ; y++)
    {
      for (z = 0 ; z < gca->height ; z++)
      {
        gcan = &gca->nodes[x][y][z] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gca->tissue_parms[(int)gcan->labels[n]].label = gcan->labels[n] ;
          if (gcan->labels[n] > biggest_label)
            biggest_label = gcan->labels[n] ;
        }
      }
    }
  }

  if (lta)
    m_L = lta->xforms[0].m_L ;
  else
    m_L = NULL ;

  for (label = 0 ; label <= biggest_label ; label++)
  {
    if (gca->tissue_parms[label].label <= 0)
      continue ;
    gca_tp = &gca->tissue_parms[label] ;
    
    for (z = 0 ; z < mri_T1->depth ; z++)
    {
      V3_Z(v_parc) = z ;
      for (y = 0 ; y < mri_T1->height ; y++)
      {
        V3_Y(v_parc) = y ;
        for (x = 0 ; x < mri_T1->width ; x++)
        {
          if (MRIvox(mri_labeled, x, y, z) != label)
            continue ;
          if (borderVoxel(mri_labeled, x, y, z))
            continue ;
          if (lta)
          {
            MATRIX *m_tmp ;
            V3_X(v_parc) = x ;
            MatrixMultiply(m_L, v_parc, v_T1) ;
            xp = nint(V3_X(v_T1)) ; yp = nint(V3_Y(v_T1)) ; 
            zp = nint(V3_Z(v_T1)) ;
            
            m_tmp = MRIgetVoxelToRasXform(mri_labeled) ;
            v_ras_cor = MatrixMultiply(m_tmp, v_parc, v_ras_cor);
            MatrixFree(&m_tmp) ;
            m_tmp = MRIgetVoxelToRasXform(mri_T1) ;
            v_ras_flash = MatrixMultiply(m_tmp, v_T1, v_ras_flash);
            MatrixFree(&m_tmp) ;
            if (!x && !y && !z && 0)
            {
              MatrixPrint(stdout, v_ras_cor) ;
              MatrixPrint(stdout, v_ras_flash) ;
            }

            if ((xp < 0 || xp >= mri_T1->width) ||
                (yp < 0 || yp >= mri_T1->height) ||
                (zp < 0 || zp >= mri_T1->depth))
              continue ;
          }
          else
          {
            xp = x ; yp = y ; zp = z ;
          }

          T1 = MRISvox(mri_T1, xp, yp, zp) ;
          PD = MRISvox(mri_PD, xp, yp, zp) ;
          gca_tp->total_training++ ;
          gca_tp->T1_mean += T1 ;
          gca_tp->T1_var += T1*T1 ;
          gca_tp->PD_mean += PD ;
          gca_tp->PD_var += PD*PD ;
        }
      }
    }
  }
  return(NO_ERROR) ;
}

int
GCAnormalizeTissueStatistics(GCA *gca)
{
  int              n ;
  double           nsamples ;
  GCA_TISSUE_PARMS *gca_tp ;

  for (n = 0 ; n < MAX_GCA_LABELS ; n++)
  {
    gca_tp = &gca->tissue_parms[n] ;
    if (gca_tp->total_training <= 0)
      continue ;
    nsamples = gca_tp->total_training ;
    gca_tp->T1_mean /= nsamples ;
    gca_tp->PD_mean /= nsamples ;
    gca_tp->T2_mean /= nsamples ;
    gca_tp->T1_var = 
      gca_tp->T1_var / nsamples - gca_tp->T1_mean*gca_tp->T1_mean;
    gca_tp->PD_var = 
      gca_tp->PD_var / nsamples - gca_tp->PD_mean*gca_tp->PD_mean;
    printf("%30.30s: T1=%4d +- %4d, PD=%4d +- %4d \n", 
           cma_label_to_name(n), 
           nint(gca_tp->T1_mean), nint(sqrt(gca_tp->T1_var)),
           nint(gca_tp->PD_mean), nint(sqrt(gca_tp->PD_var))) ;
  }

  return(NO_ERROR) ;
}


char *
cma_label_to_name(int label)
{
  static char name[100] ;

  if (label == Unknown)
    return("Unknown") ;
  if (label == Left_Cerebral_Exterior)
    return("Left_Cerebral_Exterior") ;
  if (label == Left_Cerebral_White_Matter)
    return("Left_Cerebral_White_Matter") ;
  if (label == Left_Cerebral_Cortex)
    return("Left_Cerebral_Cortex") ;
  if (label == Left_Lateral_Ventricle)
    return("Left_Lateral_Ventricle") ;
  if (label == Left_Inf_Lat_Vent)
    return("Left_Inf_Lat_Vent") ;
  if (label == Left_Cerebellum_Exterior)
    return("Left_Cerebellum_Exterior") ;
  if (label == Left_Cerebellum_White_Matter)
    return("Left_Cerebellum_White_Matter") ;
  if (label == Left_Cerebellum_Cortex)
    return("Left_Cerebellum_Cortex") ;
  if (label == Left_Thalamus)
    return("Left_Thalamus") ;
  if (label == Left_Thalamus_Proper)
    return("Left_Thalamus_Proper") ;
  if (label == Left_Caudate)
    return("Left_Caudate") ;
  if (label == Left_Putamen)
    return("Left_Putamen") ;
  if (label == Left_Pallidum)
    return("Left_Pallidum") ;
  if (label == Third_Ventricle)
    return("Third_Ventricle") ;
  if (label == Fourth_Ventricle)
    return("Fourth_Ventricle") ;
  if (label == Brain_Stem)
    return("Brain_Stem") ;
  if (label == Left_Hippocampus)
    return("Left_Hippocampus") ;
  if (label == Left_Amygdala)
    return("Left_Amygdala") ;
  if (label == Left_Insula)
    return("Left_Insula") ;
  if (label == Left_Operculum)
    return("Left_Operculum") ;
  if (label == Line_1)
    return("Line_1") ;
  if (label == Line_2)
    return("Line_2") ;
  if (label == Line_3)
    return("Line_3") ;
  if (label == CSF)
    return("CSF") ;
  if (label == Left_Lesion)
    return("Left_Lesion") ;
  if (label == Left_Accumbens_area)
    return("Left_Accumbens_area") ;
  if (label == Left_Substancia_Nigra)
    return("Left_Substancia_Nigra") ;
  if (label == Left_VentralDC)
    return("Left_VentralDC") ;
  if (label == Left_undetermined)
    return("Left_undetermined") ;
  if (label == Left_vessel)
    return("Left_vessel") ;
  if (label == Left_choroid_plexus)
    return("Left_choroid_plexus") ;
  if (label == Left_F3orb)
    return("Left_F3orb") ;
  if (label == Left_lOg)
    return("Left_lOg") ;
  if (label == Left_aOg)
    return("Left_aOg") ;
  if (label == Left_mOg)
    return("Left_mOg") ;
  if (label == Left_pOg)
    return("Left_pOg") ;
  if (label == Left_Stellate)
    return("Left_Stellate") ;
  if (label == Left_Porg)
    return("Left_Porg") ;
  if (label == Left_Aorg)
    return("Left_Aorg") ;
  if (label == Right_Cerebral_Exterior)
    return("Right_Cerebral_Exterior") ;
  if (label == Right_Cerebral_White_Matter)
    return("Right_Cerebral_White_Matter") ;
  if (label == Right_Cerebral_Cortex)
    return("Right_Cerebral_Cortex") ;
  if (label == Right_Lateral_Ventricle)
    return("Right_Lateral_Ventricle") ;
  if (label == Right_Inf_Lat_Vent)
    return("Right_Inf_Lat_Vent") ;
  if (label == Right_Cerebellum_Exterior)
    return("Right_Cerebellum_Exterior") ;
  if (label == Right_Cerebellum_White_Matter)
    return("Right_Cerebellum_White_Matter") ;
  if (label == Right_Cerebellum_Cortex)
    return("Right_Cerebellum_Cortex") ;
  if (label == Right_Thalamus)
    return("Right_Thalamus") ;
  if (label == Right_Thalamus_Proper)
    return("Right_Thalamus_Proper") ;
  if (label == Right_Caudate)
    return("Right_Caudate") ;
  if (label == Right_Putamen)
    return("Right_Putamen") ;
  if (label == Right_Pallidum)
    return("Right_Pallidum") ;
  if (label == Right_Hippocampus)
    return("Right_Hippocampus") ;
  if (label == Right_Amygdala)
    return("Right_Amygdala") ;
  if (label == Right_Insula)
    return("Right_Insula") ;
  if (label == Right_Operculum)
    return("Right_Operculum") ;
  if (label == Right_Lesion)
    return("Right_Lesion") ;
  if (label == Right_Accumbens_area)
    return("Right_Accumbens_area") ;
  if (label == Right_Substancia_Nigra)
    return("Right_Substancia_Nigra") ;
  if (label == Right_VentralDC)
    return("Right_VentralDC") ;
  if (label == Right_undetermined)
    return("Right_undetermined") ;
  if (label == Right_vessel)
    return("Right_vessel") ;
  if (label == Right_choroid_plexus)
    return("Right_choroid_plexus") ;
  if (label == Right_F3orb)
    return("Right_F3orb") ;
  if (label == Right_lOg)
    return("Right_lOg") ;
  if (label == Right_aOg)
    return("Right_aOg") ;
  if (label == Right_mOg)
    return("Right_mOg") ;
  if (label == Right_pOg)
    return("Right_pOg") ;
  if (label == Right_Stellate)
    return("Right_Stellate") ;
  if (label == Right_Porg)
    return("Right_Porg") ;
  if (label == Right_Aorg)
    return("Right_Aorg") ;
  sprintf(name, "UNKNOWN (%d)", label) ;
  return(name) ;
}
MRI *
GCArelabel_cortical_gray_and_white(GCA *gca, MRI *mri_inputs, 
                                   MRI *mri_src, MRI *mri_dst,LTA *lta)
{
  int      nchanged, x, y, z, width, height, depth, total_changed,
           label, xn, yn, zn, left, new_wm, new_gray ;
  MRI      *mri_tmp ;
  GCA_NODE *gcan ;
  float    wm_mean, wm_var, gray_mean, gray_var, gray_dist, wm_dist, val ;

  if (mri_src != mri_dst)
    mri_dst = MRIcopy(mri_src, mri_dst) ;

  mri_tmp = MRIcopy(mri_dst, NULL) ;

  width = mri_src->width ; height = mri_src->height ; depth = mri_src->depth ;

  total_changed = new_wm = new_gray = 0 ;
  do
  {
    nchanged = 0 ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == 99 && y == 129 && z == 127)
            DiagBreak() ;  /* gray should be wm */
          if (x == 153 && y == 124 && z == 121)
            DiagBreak() ;  /* gm should be wm */

          label = MRIvox(mri_src, x, y, z) ;
          if (label != Left_Cerebral_Cortex && 
              label != Left_Cerebral_White_Matter &&
              label != Right_Cerebral_Cortex && 
              label != Right_Cerebral_White_Matter)
            continue ;

          val = (float)MRIvox(mri_inputs, x, y, z) ;
          GCAsourceVoxelToNode(gca, mri_dst, lta, x, y, z, &xn, &yn, &zn) ;
          gcan = &gca->nodes[xn][yn][zn] ;
          if (label == Left_Cerebral_Cortex || 
              label == Left_Cerebral_White_Matter)
          {
            left = 1 ;
            gray_mean = getLabelMean(gcan, Left_Cerebral_Cortex, &gray_var) ;
            wm_mean = getLabelMean(gcan,Left_Cerebral_White_Matter,&wm_var);
          }
          else
          {
            left = 0 ;
            gray_mean = getLabelMean(gcan, Right_Cerebral_Cortex, &gray_var) ;
            wm_mean =getLabelMean(gcan,Right_Cerebral_White_Matter,&wm_var);
          }
#if 0
          gray_dist = fabs(gray_mean - val) / sqrt(gray_var) ;
          wm_dist = fabs(wm_mean - val) / sqrt(wm_var) ;
#else
          gray_dist = fabs(gray_mean - val) ;
          wm_dist = fabs(wm_mean - val) ;
#endif
          if (gray_dist < wm_dist)
            label = left ? Left_Cerebral_Cortex : Right_Cerebral_Cortex ;
          else
            label = 
              left ? Left_Cerebral_White_Matter : Right_Cerebral_White_Matter ;
          if (label != MRIvox(mri_dst, x, y, z))
          {
            if (label == Left_Cerebral_Cortex ||label == Right_Cerebral_Cortex)
              new_gray++ ;
            else
              new_wm++ ;
            nchanged++ ;
            MRIvox(mri_tmp, x, y, z) = label ;
          }
        }
      }
    }
    MRIcopy(mri_tmp, mri_dst) ;
    total_changed += nchanged ;
  } while (nchanged > 0) ;

  MRIfree(&mri_tmp) ;
  printf("%d gm and wm labels changed (%%%2.0f to gray, %%%2.0f to white)\n", 
         total_changed, 100.0f*(float)new_gray/total_changed,
         100.0f*(float)new_wm/total_changed) ;
  return(mri_dst) ;
}
int
GCAdump(GCA *gca,MRI *mri,int x, int y, int z, LTA *lta, FILE *fp, int verbose)
{
  int      xn, yn, zn ;
  GCA_NODE *gcan ;

  GCAsourceVoxelToNode(gca, mri, lta, x, y, z, &xn, &yn, &zn) ;
  printf("\nGCA node at voxel (%d, %d, %d) --> node (%d, %d, %d)\n",
         x, y, z, xn, yn, zn) ;
  gcan = &gca->nodes[xn][yn][zn] ;
  dump_gcan(gcan, fp, verbose) ;
  return(NO_ERROR) ;
}

static GCA_SAMPLE *
gcaExtractLabelAsSamples(GCA *gca, MRI *mri_labeled, LTA *lta,
                         int *pnsamples, int label)
{
  int         i, nsamples, width, height, depth, x, y, z, xn, yn, zn, n ;
  GCA_SAMPLE  *gcas ;
  GCA_NODE    *gcan ;
  GC1D        *gc ;

  width = mri_labeled->width ; 
  height = mri_labeled->height ; 
  depth = mri_labeled->depth ; 

  for (nsamples = z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (MRIvox(mri_labeled, x, y, z) != label)
          continue ;
        nsamples++ ;
      }
    }
  }

  gcas = (GCA_SAMPLE *)calloc(nsamples, sizeof(GCA_SAMPLE)) ;
  if (!gcas)
    ErrorExit(ERROR_NOMEMORY, 
              "gcaExtractLabelAsSamples(%d): could not allocate %d samples\n",
              label, nsamples) ;


  for (i = z = 0 ; z < depth ; z++)
  {
    for (y = 0 ; y < height ; y++)
    {
      for (x = 0 ; x < width ; x++)
      {
        if (MRIvox(mri_labeled, x, y, z) != label)
          continue ;

        GCAsourceVoxelToNode(gca, mri_labeled, lta, x, y, z, &xn, &yn, &zn) ;
        gcan = &gca->nodes[xn][yn][zn] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          if (gcan->labels[n] == label)
            break ;
        }

        gc = &gcan->gcs[n] ;
        gcas[i].label = label ;
        gcas[i].xn = xn ; gcas[i].yn = yn ; gcas[i].zn = zn ; 
        gcas[i].x = x ;   gcas[i].y = y ;   gcas[i].z = z ; 
        gcas[i].var = gc->var ;
        gcas[i].mean = gc->mean ;
        gcas[i].prior = gc->prior ;
        i++ ;
      }
    }
  }

  *pnsamples = nsamples ;
  return(gcas) ;
}


#define SAMPLE_PCT       0.20
#define MIN_SAMPLES      10
#define MEAN_RESOLUTION  8

int
GCArenormalizeLabels(MRI *mri_in, MRI *mri_labeled, GCA *gca, LTA *lta)
{
  int              x, y, z, n, label, biggest_label, nsamples, xv, yv, zv,
                   *ordered_indices, i, index, width, height, depth ;
  float            mean, var, val ;
  GCA_NODE         *gcan ;
  MATRIX           *m_L ;
  GCA_SAMPLE       *gcas ;
  GC1D             *gc ;
  MRI              *mri_means, *mri_control, *mri_tmp ;
  char             fname[STRLEN] ;
  
  /* first build a list of all labels that exist */
  for (biggest_label = x = 0 ; x < gca->width ; x++)
  {
    for (y = 0 ; y < gca->height ; y++)
    {
      for (z = 0 ; z < gca->height ; z++)
      {
        gcan = &gca->nodes[x][y][z] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          if (gcan->labels[n] > biggest_label)
            biggest_label = gcan->labels[n] ;
        }
      }
    }
  }

  if (lta)
    m_L = lta->xforms[0].m_L ;
  else
    m_L = NULL ;

  for (label = 1 ; label <= biggest_label ; label++)
  {
    gcas = gcaExtractLabelAsSamples(gca, mri_labeled,lta,&nsamples,label);
    if (!nsamples)
      continue ;
    if (nsamples < MIN_SAMPLES/SAMPLE_PCT)
    {
      free(gcas) ;
      continue ;
    }
    ordered_indices = (int *)calloc(nsamples, sizeof(int)) ;
    GCAcomputeLogSampleProbability(gca, gcas, mri_in, m_L, nsamples) ;
    GCArankSamples(gca, gcas, nsamples, ordered_indices) ;

    if (nint(nsamples*SAMPLE_PCT) < MIN_SAMPLES)
      nsamples = MIN_SAMPLES ;
    else
      nsamples = nint(SAMPLE_PCT * (float)nsamples) ;
    

    for (var = mean = 0.0f, i = 0 ; i < nsamples ; i++)
    {
      index = ordered_indices[i] ;
      val = (float)MRIvox(mri_in, gcas[index].x, gcas[index].y, gcas[index].z);
      mean += val ;
      var += val*val ;
    }
    mean /= (float)nsamples ; var = var / nsamples - mean*mean ; 
    var = sqrt(var);
    printf("label %s: using %d samples to estimate mean = %2.1f +- %2.1f\n",
           cma_label_to_name(label), nsamples, mean, var) ;
    
    width = nint(mri_in->width / MEAN_RESOLUTION) ;
    height = nint(mri_in->height / MEAN_RESOLUTION) ;
    depth = nint(mri_in->depth / MEAN_RESOLUTION) ;
    mri_means = MRIalloc(width, height, depth, MRI_FLOAT) ;
    mri_control = MRIalloc(width, height, depth, MRI_SHORT) ;
    MRIsetResolution(mri_means, 
                     MEAN_RESOLUTION,MEAN_RESOLUTION,MEAN_RESOLUTION) ;
    MRIsetResolution(mri_control, 
                     MEAN_RESOLUTION,MEAN_RESOLUTION,MEAN_RESOLUTION) ;

    GCAlabelMri(gca, mri_means, label, lta) ;
    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
    {
      sprintf(fname, "%s_label.mgh", cma_label_to_name(label)) ;
      MRIwrite(mri_means, fname) ;
    }

    for (mean = 0.0f, n = z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (MRIFvox(mri_means, x, y, z) > 1)
          {
            n++ ;
            mean += MRIFvox(mri_means, x, y, z) ;
          }
        }
      }
    }
    mean /= (float)n ;
    printf("mean GCA value %2.1f\n", mean) ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (MRIFvox(mri_means, x, y, z) < 1)
          {
            MRIFvox(mri_means, x, y, z) = mean ;
          }
        }
      }
    }
          
    for (i = 0 ; i < nsamples ; i++)
    {
      index = ordered_indices[i] ;
      GCAnodeToSourceVoxel(gca, mri_means, m_L, 
                           gcas[index].xn, gcas[index].yn, gcas[index].zn, 
                           &x, &y, &z) ;
      val = (float)MRIvox(mri_in, gcas[index].x, gcas[index].y, gcas[index].z);
      if (x == 19 && y == 14 && z == 15)
        DiagBreak() ;
      if (MRISvox(mri_control, x, y, z) == 0)
        MRIFvox(mri_means, x, y, z) = val ;
      else
        MRIFvox(mri_means, x, y, z) += val ;
      MRISvox(mri_control, x, y, z)++ ;
    }
        
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == 19 && y == 14 && z == 15)
            DiagBreak() ;
          if (MRISvox(mri_control, x, y, z) > 0)
          {
            MRIFvox(mri_means,x,y,z) = 
              nint((float)MRIFvox(mri_means,x,y,z)/
                   (float)MRISvox(mri_control,x,y,z)) ;
            MRISvox(mri_control, x, y, z) = 1 ;
          }
        }
      }
    }

    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
    {
      sprintf(fname, "%s_means.mgh", cma_label_to_name(label)) ;
      MRIwrite(mri_means, fname) ;
    }

    mri_tmp = MRIalloc(width, height, depth, MRI_SHORT) ;
    MRIcopy(mri_means, mri_tmp) ;
    MRIfree(&mri_means) ; mri_means = mri_tmp ;

    mri_tmp = MRIalloc(width, height, depth, MRI_UCHAR) ;
    MRIcopy(mri_control, mri_tmp) ;
    MRIfree(&mri_control) ; mri_control = mri_tmp ;

    MRIsoapBubble(mri_means, mri_control, mri_means, 10) ;
    MRIclear(mri_control) ;
    MRIsoapBubble(mri_means, mri_control, mri_means, 1) ;
    
    if (DIAG_VERBOSE_ON && Gdiag & DIAG_WRITE)
    {
      sprintf(fname, "%s_soap.mgh", cma_label_to_name(label)) ;
      MRIwrite(mri_means, fname) ;

      sprintf(fname, "%s_control.mgh", cma_label_to_name(label)) ;
      MRIwrite(mri_control, fname) ;
    }


    width = gca->width ; height = gca->height ; depth = gca->depth ;
    for (z = 0 ; z < depth ; z++)
    {
      for (y = 0 ; y < height ; y++)
      {
        for (x = 0 ; x < width ; x++)
        {
          if (x == Gx && y == Gy && z == Gz)
            DiagBreak() ;
          gcan = &gca->nodes[x][y][z] ;
          for (n = 0 ; n < gcan->nlabels ; n++)
          {
            if (gcan->labels[n] == label)
            {
              if (x == Gx && y == Gy && z == Gz)
                DiagBreak() ;
              gc = &gcan->gcs[n] ;
              GCAnodeToVoxel(gca, mri_means, x, y, z, &xv, &yv, &zv) ;
              gc->mean = (float)MRISvox(mri_means, xv, yv, zv);
              break ;
            }
          }
        }
      }
    }

    MRIfree(&mri_means) ; MRIfree(&mri_control) ;
    free(gcas) ; free(ordered_indices) ;
  }

  return(NO_ERROR) ;
}

int
GCArenormalizeIntensities(GCA *gca, int *labels, float *intensities, int num)
{
  float     wm_mean, scale_factor, gray_mean ;
  int       xn, yn, zn, n, i, label ;
  GCA_NODE  *gcan ;
  GC1D      *gc ;
  double    *means, *wts ;

  means = (double *)calloc(num, sizeof(double)) ;
  wts = (double *)calloc(num, sizeof(double)) ;
  if (!means || !wts)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d wt and mean vectors",
              Progname, num) ;

  /* compute overall white matter mean to use as anchor for rescaling */
  for (zn = 0 ; zn < gca->depth ; zn++)
  {
    for (yn = 0 ; yn < gca->height ; yn++)
    {
      for (xn = 0 ; xn < gca->width ; xn++)
      {
        gcan = &gca->nodes[xn][yn][zn] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          /* find index in lookup table for this label */
          label = gcan->labels[n] ;
          for (i = 0 ; i < num ; i++)
            if (label == labels[i])
              break ;

          if (i >= num)
            continue ;  /* shouldn't happen */

          gc = &gcan->gcs[n] ;
          wts[i] += gc->prior ;
          means[i] += gc->mean*gc->prior ;
                           
        }
      }
    }
  }
  gray_mean = scale_factor = wm_mean = -1 ;
  for (i = 0 ; i < num ; i++)
  {
    if (!FZERO(wts[i]))
      means[i] /= wts[i] ;
    if (labels[i] == Left_Cerebral_White_Matter)
    {
      wm_mean = means[i] ;
      scale_factor = wm_mean / intensities[i] ;
    }
    if (labels[i] == Left_Cerebral_Cortex)
      gray_mean = means[i] ;
  }
  if (wm_mean < 0)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "GCArenormalizeIntensities: "
                 "could not find white matter in intensity table")) ;
  printf("mean wm intensity = %2.1f, scaling gray to %2.2f (from %2.2f)\n", 
         wm_mean, scale_factor*gray_mean, gray_mean) ;

  /* now go through each label and scale it by gca wm to norm wm ratio */
  /* compute overall white matter mean to use as anchor for rescaling */
  for (zn = 0 ; zn < gca->depth ; zn++)
  {
    for (yn = 0 ; yn < gca->height ; yn++)
    {
      for (xn = 0 ; xn < gca->width ; xn++)
      {
        gcan = &gca->nodes[xn][yn][zn] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          label = gcan->labels[n] ;
          if (label == Gdiag_no)
            DiagBreak() ;

          /* find index in lookup table for this label */
          for (i = 0 ; i < num ; i++)
            if (label == labels[i])
              break ;
          if (i >= num)
            continue ;
          
          gc = &gcan->gcs[n] ;
          gc->mean = scale_factor*intensities[i]*gc->mean/means[i] ;
        }
      }
    }
  }

  free(wts) ;
  free(means) ;
  return(NO_ERROR) ;
}

int
GCAunifyVariance(GCA *gca)
{
  int        xn, yn, zn, n ;
  GCA_NODE   *gcan ;
  GC1D       *gc ;
 
  for (zn = 0 ; zn < gca->depth ; zn++)
  {
    for (yn = 0 ; yn < gca->height ; yn++)
    {
      for (xn = 0 ; xn < gca->width ; xn++)
      {
        gcan = &gca->nodes[xn][yn][zn] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          gc->var = 25.0f ;
        }
      }
    }
  }

  return(NO_ERROR) ;
}


double
GCAlabelMean(GCA *gca, int label)
{
  int       xn, yn, zn, n ;
  GCA_NODE  *gcan ;
  GC1D      *gc ;
  double    wt, mean ;

  /* compute overall white matter mean to use as anchor for rescaling */
  for (wt = mean = 0.0, zn = 0 ; zn < gca->depth ; zn++)
  {
    for (yn = 0 ; yn < gca->height ; yn++)
    {
      for (xn = 0 ; xn < gca->width ; xn++)
      {
        gcan = &gca->nodes[xn][yn][zn] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          /* find index in lookup table for this label */
          if (gcan->labels[n] != label)
            continue ;
          gc = &gcan->gcs[n] ;
          wt += gc->prior ;
          mean += gc->mean*gc->prior ;
                           
        }
      }
    }
  }
  mean /= wt ;
  return(mean) ;
}

int
GCAregularizeConditionalDensities(GCA *gca, float smooth)
{
  int       xn, yn, zn, n, i, label, max_label ;
  GCA_NODE  *gcan ;
  GC1D      *gc ;
  double    *means, *wts ;

  means = (double *)calloc(MAX_GCA_LABELS, sizeof(double)) ;
  wts = (double *)calloc(MAX_GCA_LABELS, sizeof(double)) ;
  if (!means || !wts)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d wt and mean vectors",
              Progname, MAX_GCA_LABELS) ;

  /* compute overall mean for each class */
  for (max_label = 1, zn = 0 ; zn < gca->depth ; zn++)
  {
    for (yn = 0 ; yn < gca->height ; yn++)
    {
      for (xn = 0 ; xn < gca->width ; xn++)
      {
        gcan = &gca->nodes[xn][yn][zn] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          /* find index in lookup table for this label */
          label = gcan->labels[n] ;

          gc = &gcan->gcs[n] ;
          wts[label] += gc->prior ;
          means[label] += gc->mean*gc->prior ;
          if (label > max_label)
            max_label = label ;
        }
      }
    }
  }

  for (i = 0 ; i <= max_label ; i++)
  {
    if (!FZERO(wts[i]))
      means[i] /= wts[i] ;
  }

  /* now impose regularization */
  for (zn = 0 ; zn < gca->depth ; zn++)
  {
    for (yn = 0 ; yn < gca->height ; yn++)
    {
      for (xn = 0 ; xn < gca->width ; xn++)
      {
        gcan = &gca->nodes[xn][yn][zn] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          /* find index in lookup table for this label */
          label = gcan->labels[n] ;
          if (label <= 0)
            continue ;

          gc = &gcan->gcs[n] ;
          gc->mean = means[label]*smooth + gc->mean*(1.0f-smooth) ;
        }
      }
    }
  }

  free(wts) ; free(means) ;
  return(NO_ERROR) ;
}

int
GCAhistogramTissueStatistics(GCA *gca, MRI *mri_T1,MRI *mri_PD,
                              MRI *mri_labeled, LTA *lta, char *fname)
{
  int              x, y, z, n, label, biggest_label, T1, PD, xp, yp, zp ;
  GCA_NODE         *gcan ;
  GCA_TISSUE_PARMS *gca_tp ;
  VECTOR           *v_parc, *v_T1 ;
  static VECTOR *v_ras_cor = NULL, *v_ras_flash ;
  MATRIX           *m_L ;
  FILE             *fp ;
  
  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "GCAhistogramTissueStatistics: could not open %s",
              fname) ;
  
  v_parc = VectorAlloc(4, MATRIX_REAL) ;
  v_T1 = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v_parc, 4, 1) = 1.0 ;
  *MATRIX_RELT(v_T1, 4, 1) = 1.0 ;

  /* first build a list of all labels that exist */
  for (biggest_label = x = 0 ; x < gca->width ; x++)
  {
    for (y = 0 ; y < gca->height ; y++)
    {
      for (z = 0 ; z < gca->height ; z++)
      {
        gcan = &gca->nodes[x][y][z] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gca->tissue_parms[(int)gcan->labels[n]].label = gcan->labels[n] ;
          if (gcan->labels[n] > biggest_label)
            biggest_label = gcan->labels[n] ;
        }
      }
    }
  }

  if (lta)
    m_L = lta->xforms[0].m_L ;
  else
    m_L = NULL ;

  for (label = 0 ; label <= biggest_label ; label++)
  {
    if (gca->tissue_parms[label].label <= 0)
      continue ;
    gca_tp = &gca->tissue_parms[label] ;
    
    for (z = 0 ; z < mri_T1->depth ; z++)
    {
      V3_Z(v_parc) = z ;
      for (y = 0 ; y < mri_T1->height ; y++)
      {
        V3_Y(v_parc) = y ;
        for (x = 0 ; x < mri_T1->width ; x++)
        {
          if (MRIvox(mri_labeled, x, y, z) != label)
            continue ;
          if (borderVoxel(mri_labeled, x, y, z))
            continue ;
          if (lta)
          {
            MATRIX *m_tmp ;
            V3_X(v_parc) = x ;
            MatrixMultiply(m_L, v_parc, v_T1) ;
            xp = nint(V3_X(v_T1)) ; yp = nint(V3_Y(v_T1)) ; 
            zp = nint(V3_Z(v_T1)) ;
            
            m_tmp = MRIgetVoxelToRasXform(mri_labeled) ;
            v_ras_cor = MatrixMultiply(m_tmp, v_parc, v_ras_cor);
            MatrixFree(&m_tmp) ;
            m_tmp = MRIgetVoxelToRasXform(mri_T1) ;
            v_ras_flash = MatrixMultiply(m_tmp, v_T1, v_ras_flash);
            MatrixFree(&m_tmp) ;
            if (!x && !y && !z && 0)
            {
              MatrixPrint(stdout, v_ras_cor) ;
              MatrixPrint(stdout, v_ras_flash) ;
            }

            if ((xp < 0 || xp >= mri_T1->width) ||
                (yp < 0 || yp >= mri_T1->height) ||
                (zp < 0 || zp >= mri_T1->depth))
              continue ;
          }
          else
          {
            xp = x ; yp = y ; zp = z ;
          }

          T1 = MRISvox(mri_T1, xp, yp, zp) ;
          PD = MRISvox(mri_PD, xp, yp, zp) ;
          fprintf(fp, "%d %d %d\n", label, T1, PD) ;
          gca_tp->total_training++ ;
          gca_tp->T1_mean += T1 ;
          gca_tp->T1_var += T1*T1 ;
          gca_tp->PD_mean += PD ;
          gca_tp->PD_var += PD*PD ;
        }
      }
    }
  }

  fclose(fp) ;
  return(NO_ERROR) ;
}
int
GCArenormalizeToFlash(GCA *gca, char *tissue_parms_fname, MRI *mri)
{
  FILE     *fp ;
  char     *cp, line[STRLEN] ;
  int      labels[MAX_GCA_LABELS], nlabels ;
  float   intensities[MAX_GCA_LABELS], TR, alpha, T1, PD ;
  
  TR = mri->tr ; alpha = mri->flip_angle ;

  fp = fopen(tissue_parms_fname, "r") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE, 
                 "GCArenormalizeToFlash: could not open tissue parms file %s",
                 tissue_parms_fname)) ;

  cp = fgetl(line, STRLEN-1, fp) ; nlabels = 0 ;
  while (cp)
  {
    if (sscanf(cp, "%d %f %f", &labels[nlabels], &T1, &PD)
        != 3)
      ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE, 
                 "GCArenormalizeToFlash: could not parse %dth line %s in %s",
                 nlabels+1, cp, tissue_parms_fname)) ;
                
    intensities[nlabels] = FLASHforwardModel(alpha, TR, PD, T1) ;

    nlabels++ ;
    cp = fgetl(line, STRLEN-1, fp) ;
  }
  fclose(fp) ;
  GCArenormalizeIntensities(gca, labels, intensities, nlabels) ;
  return(NO_ERROR) ;
}

static double
FLASHforwardModel(double flip_angle, double TR, double PD, double T1)
{
  double FLASH, E1 ;
  double  CFA, SFA ;


  CFA = cos(flip_angle) ; SFA = sin(flip_angle) ;
  E1 = exp(-TR/T1) ;
      
  FLASH = PD * SFA ;
  if (!DZERO(T1))
    FLASH *= (1-E1)/(1-CFA*E1);
  return(FLASH) ;
}
int
GCAmeanFilterConditionalDensities(GCA *gca, float navgs)
{
  int       xn, yn, zn, xn1, yn1, zn1, n, label, max_label, niter ;
  GCA_NODE  *gcan ;
  GC1D      *gc ;
  double    mean, wt, var ;
  MRI       *mri_means, *mri_vars ;

  mri_means = MRIalloc(gca->width, gca->height, gca->depth, MRI_FLOAT) ;
  mri_vars = MRIalloc(gca->width, gca->height, gca->depth, MRI_FLOAT) ;

  /* compute overall mean for each class */
  for (max_label = 1, zn = 0 ; zn < gca->depth ; zn++)
  {
    for (yn = 0 ; yn < gca->height ; yn++)
    {
      for (xn = 0 ; xn < gca->width ; xn++)
      {
        gcan = &gca->nodes[xn][yn][zn] ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          /* find index in lookup table for this label */
          label = gcan->labels[n] ;

          gc = &gcan->gcs[n] ;
          if (label > max_label)
            max_label = label ;
        }
      }
    }
  }

  /* now impose regularization */
  for (niter = 0 ; niter < navgs ; niter++)
  {
    for (label = 0 ; label <= max_label ; label++)
    {
      if (IS_UNKNOWN(label))
        continue ;
      for (zn = 0 ; zn < gca->depth ; zn++)
      {
        for (yn = 0 ; yn < gca->height ; yn++)
        {
          for (xn = 0 ; xn < gca->width ; xn++)
          {
            mean = var = wt = 0.0 ;
            for (xn1 = xn-1 ; xn1 <= xn+1 ; xn1++)
            {
              if (xn1 < 0 || xn1 >= gca->width)
                continue ;
              for (yn1 = yn-1 ; yn1 <= yn+1 ; yn1++)
              {
                if (yn1 < 0 || yn1 >= gca->height)
                  continue ;
                for (zn1 = zn-1 ; zn1 <= zn+1 ; zn1++)
                {
                  if (zn1 < 0 || zn1 >= gca->depth)
                    continue ;
                  if (xn1 == xn && yn1 == yn && zn1 == zn)
                    DiagBreak() ;
                  gcan = &gca->nodes[xn1][yn1][zn1] ;
                  for (n = 0 ; n < gcan->nlabels ; n++)
                  {
                    /* find index in lookup table for this label */
                    if (gcan->labels[n] != label)
                      continue ;
                  
                    gc = &gcan->gcs[n] ;
                    mean += gc->prior*gc->mean ;
                    var += gc->prior*gc->var ;
                    wt += gc->prior ;
                    if (!finite(mean / wt))
                      DiagBreak() ;
                    if (!finite(var / wt))
                      DiagBreak() ;
                    break ;
                  }
                }
              }
            }
            if (FZERO(wt))
              continue ;   /* label didn't occur here */
            if (!finite(mean / wt))
              DiagBreak() ;
            if (!finite(var / wt))
              DiagBreak() ;
            MRIFvox(mri_means, xn, yn, zn) = mean / wt ;
            MRIFvox(mri_vars, xn, yn, zn) = var / wt ;
            if (!finite(MRIFvox(mri_means,xn,yn,zn)))
              DiagBreak() ;
            if (!finite(MRIFvox(mri_vars, xn, yn, zn)))
              DiagBreak() ;
          }
        }
      }
      
      for (zn = 0 ; zn < gca->depth ; zn++)
      {
        for (yn = 0 ; yn < gca->height ; yn++)
        {
          for (xn = 0 ; xn < gca->width ; xn++)
          {
            gcan = &gca->nodes[xn][yn][zn] ;
            for (n = 0 ; n < gcan->nlabels ; n++)
            {
              if (gcan->labels[n] != label)
                continue ;
              gc = &gcan->gcs[n] ;
              gc->mean = MRIFvox(mri_means, xn, yn, zn) ;
              gc->var = MRIFvox(mri_vars, xn, yn, zn) ;
              if (!finite(gc->mean))
                DiagBreak() ;
              if (!finite(gc->var))
                DiagBreak() ;
              break ;
            }
          }
        }
      }
    }
  }

  MRIfree(&mri_vars) ; MRIfree(&mri_means) ;
  return(NO_ERROR) ;
}

#define MIN_BIN  50
#define OFFSET_SIZE  25
#define BOX_SIZE   60   /* mm */
#define HALF_BOX   (BOX_SIZE/2)
int
GCAhistoScaleImageIntensities(GCA *gca, MRI *mri)
{
  float      x0, y0, z0 ;
  int        mri_peak ;
  double     wm_mean ;
  HISTOGRAM *h_mri, *h_smooth ;
  MRI_REGION box ;

  MRIfindApproximateSkullBoundingBox(mri, 50, &box) ;
  x0 = box.x+box.dx/3 ; y0 = box.y+box.dy/3 ; z0 = box.z+box.dz/2 ;
  printf("using (%.0f, %.0f, %.0f) as brain centroid...\n",x0, y0, z0) ;
  box.x = x0 - HALF_BOX*mri->xsize ; box.dx = BOX_SIZE*mri->xsize ;
  box.y = y0 - HALF_BOX*mri->ysize ; box.dy = BOX_SIZE*mri->ysize ;
  box.z = z0 - HALF_BOX*mri->zsize ; box.dz = BOX_SIZE*mri->zsize ;
  wm_mean = GCAlabelMean(gca, Left_Cerebral_White_Matter) ;
  wm_mean = (wm_mean + GCAlabelMean(gca, Right_Cerebral_White_Matter))/2.0 ;
  printf("mean wm in atlas = %2.0f, using box (%d,%d,%d) --> (%d, %d,%d) "
         "to find MRI wm\n", wm_mean, box.x, box.y, box.z, 
         box.x+box.dx-1,box.y+box.dy-1, box.z+box.dz-1) ;
  
  h_mri = MRIhistogramRegion(mri, 0, NULL, &box) ; 
  HISTOclearBins(h_mri, h_mri, 0, MIN_BIN) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    HISTOplot(h_mri, "mri.histo") ;
  mri_peak = HISTOfindLastPeak(h_mri, HISTO_WINDOW_SIZE,MIN_HISTO_PCT);
  mri_peak = h_mri->bins[mri_peak] ;
  printf("before smoothing, mri peak at %d\n", mri_peak) ;
  h_smooth = HISTOsmooth(h_mri, NULL, 2) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    HISTOplot(h_smooth, "mri_smooth.histo") ;
  mri_peak = HISTOfindLastPeak(h_smooth, HISTO_WINDOW_SIZE,MIN_HISTO_PCT);
  mri_peak = h_mri->bins[mri_peak] ;
  printf("after smoothing, mri peak at %d, scaling input intensities "
         "by %2.3f\n", mri_peak, wm_mean/mri_peak) ;
  HISTOfree(&h_smooth) ; HISTOfree(&h_mri) ;
  MRIscalarMul(mri, mri, wm_mean/mri_peak) ;

  return(NO_ERROR) ;
}

int
GCAhisto(GCA *gca, int nbins, int **pcounts)
{
  int   *counts, x, y, z ;
  GCA_NODE *gcan ;

  *pcounts = counts = (int *)calloc(nbins+1, sizeof(int)) ;

  for (x = 0 ; x < gca->width ; x++)
  {
    for (y = 0 ; y < gca->height ; y++)
    {
      for (z = 0 ; z < gca->depth ; z++)
      {
        gcan = &gca->nodes[x][y][z] ;
        if (gcan->nlabels == 0 || (gcan->nlabels == 1 && 
                                   gcaFindGC(gca,x,y,z,Unknown) != NULL))
          continue ;
        counts[gcan->nlabels]++ ;
      }
    }
  }

  return(NO_ERROR) ;
}

