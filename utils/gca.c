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

#define GCA_VERSION         1.0
#define DEFAULT_MAX_LABELS  4

GCA *gcaAllocMax(int ninputs, float spacing, int length, int max_labels) ;
static int GCAupdateNode(GCA *gca, int xt, int yt, int zt,float val,int label);

int
GCAvoxelToNode(GCA *gca, int xv, int yv, int zv, int *pxn, int *pyn, int *pzn)
{
  int spacing = gca->spacing, spm1 = gca->spacing-1 ;

  *pxn = (xv - spm1) / spacing ;
  if (*pxn < 0)
    *pxn = 0 ;
  else if (*pxn >= gca->width)
    *pxn = gca->width-1 ;

  *pyn = (yv - spm1) / spacing ;
  if (*pyn < 0)
    *pyn = 0 ;
  else if (*pyn >= gca->height)
    *pyn = gca->height-1 ;

  *pzn = (zv - spm1) / spacing ;
  if (*pzn < 0)
    *pzn = 0 ;
  else if (*pzn >= gca->depth)
    *pzn = gca->depth-1 ;

  return(NO_ERROR) ;
}


GCA *
GCAalloc(int ninputs, float spacing, int length)
{
  return(gcaAllocMax(ninputs, spacing, length, DEFAULT_MAX_LABELS)) ;
}

GCA *
gcaAllocMax(int ninputs, float spacing, int length, int max_labels)
{
  GCA      *gca ;
  GCA_NODE *gcan ;
  int      x, y, z ;

  gca = calloc(1, sizeof(GCA)) ;
  if (!gca)
    ErrorExit(ERROR_NOMEMORY, "GCAalloc: could not allocate struct") ;

  gca->ninputs = ninputs ;
  gca->spacing = spacing ;

  length = (length - (spacing-1)) / spacing + 1 ;
  gca->width = gca->height = gca->depth = length ;/*nint(length / spacing)*/ ;

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
          gcan->gcs = (GC *)calloc(gcan->max_labels, sizeof(GC)) ;
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
GCAtrain(GCA *gca, MRI *mri_inputs, MRI *mri_labels, LTA *lta)
{
  int    x, y, z, xt, yt, zt, width, height, depth, label, val,xi,yi,zi ;
  MATRIX *m_L, *m_label_to_input, *m_voxel_to_ras, *m_ras_to_voxel ;
  VECTOR *v_input, *v_canon, *v_label ;

  v_input = VectorAlloc(4, MATRIX_REAL) ;
  v_canon = VectorAlloc(4, MATRIX_REAL) ;
  v_label = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v_input, 4, 1) = 1.0 ;
  *MATRIX_RELT(v_canon, 4, 1) = 1.0 ;
  *MATRIX_RELT(v_label, 4, 1) = 1.0 ;

  m_voxel_to_ras = MRIgetVoxelToRasXform(mri_labels) ;
  m_ras_to_voxel = MRIgetRasToVoxelXform(mri_inputs) ;
  m_label_to_input = MatrixMultiply(m_ras_to_voxel, m_voxel_to_ras, NULL) ;
  MatrixFree(&m_voxel_to_ras) ; MatrixFree(&m_ras_to_voxel) ;

  /* convert transform to voxel coordinates */
#if 0
  m_L = MRIrasXformToVoxelXform(mri_inputs,mri_inputs,lta->xforms[0].m_L,NULL);
#else
  m_L = lta->xforms[0].m_L ;
#endif

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_labels->width ; height = mri_labels->height; 
  depth = mri_labels->depth ;
  for (x = 0 ; x < width ; x++)
  {
    V3_X(v_label) = (float)x ;
    for (y = 0 ; y < height ; y++)
    {
      V3_Y(v_label) = (float)y ;
      for (z = 0 ; z < depth ; z++)
      {
        label = MRIvox(mri_labels, x, y, z) ;
#if 0
        if (!label)
          continue ;
#endif
        V3_Z(v_label) = (float)z ;
        MatrixMultiply(m_label_to_input, v_label, v_input) ;
        xi = nint(V3_X(v_input)) ;
        yi = nint(V3_Y(v_input)) ;
        zi = nint(V3_Z(v_input)) ;
        val = MRIvox(mri_inputs, xi, yi, zi) ;

        if (xi == 110 && yi == 97 && zi == 128)
          DiagBreak() ;  /* wm in 1332 */

        MatrixMultiply(m_L, v_input, v_canon) ;
        xt = nint(V3_X(v_canon)) ;
        yt = nint(V3_Z(v_canon)) ;
        zt = nint(V3_Y(v_canon)) ;
        if (xt < 0) xt = 0 ;
        if (yt < 0) yt = 0 ;
        if (zt < 0) zt = 0 ;
        if (xt >= width)  xt = width-1 ;
        if (yt >= height) yt = height-1 ;
        if (zt >= depth)  zt = depth-1 ;

        GCAupdateNode(gca, xt, yt, zt, (float)val, label) ;
      }
    }
  }

  VectorFree(&v_input) ; VectorFree(&v_canon) ; VectorFree(&v_label) ;
#if 0
  MatrixFree(&m_L) ;
#endif
  return(NO_ERROR) ;
}
int
GCAwrite(GCA *gca, char *fname)
{
  FILE     *fp ;
  int      x, y, z, n ; 
  GCA_NODE *gcan ;
  GC       *gc ;

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
        if (x == 25 && y == 59 && z == 35)
          DiagBreak() ;
        gcan = &gca->nodes[x][y][z] ;
        fwriteInt(gcan->nlabels, fp) ;
        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          fputc((int)gcan->labels[n],fp) ;
          fwriteFloat(gc->mean, fp) ;
          fwriteFloat(gc->var, fp) ;
          fwriteFloat(gc->prior, fp) ;
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
  int      x, y, z, n ; 
  GCA      *gca ;
  GCA_NODE *gcan ;
  GC       *gc ;
  float    version, spacing ;
  int      width, height, depth, ninputs ;

  fp  = fopen(fname, "rb") ;
  if (!fp)
    ErrorReturn(NULL,
                (ERROR_NOFILE,
                 "GCAread: could not open GCA %s for writing",fname)) ;

  version = freadFloat(fp) ;
  spacing = freadFloat(fp) ;
  width = freadInt(fp); height = freadInt(fp); depth = freadInt(fp);
  ninputs = freadInt(fp) ;

  gca = gcaAllocMax(ninputs, spacing, spacing*width, 0) ;
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
        gcan->labels = (char *)calloc(gcan->nlabels, sizeof(char)) ;
        if (!gcan->labels)
          ErrorExit(ERROR_NOMEMORY, "GCAread(%s); could not allocate %d "
                    "labels @ (%d,%d,%d)", fname, x, y, z) ;
        gcan->gcs = (GC *)calloc(gcan->nlabels, sizeof(GC)) ;
        if (!gcan->gcs)
          ErrorExit(ERROR_NOMEMORY, "GCAread(%s); could not allocated %d gcs "
                    "@ (%d,%d,%d)", fname, x, y, z) ;

        for (n = 0 ; n < gcan->nlabels ; n++)
        {
          gc = &gcan->gcs[n] ;
          gcan->labels[n] = (char)fgetc(fp) ;
          gc->mean = freadFloat(fp) ;
          gc->var = freadFloat(fp) ;
          gc->prior = freadFloat(fp) ;
        }
      }
    }
  }

  fclose(fp) ;
  return(gca) ;
}

static int
GCAupdateNode(GCA *gca, int xt, int yt, int zt, float val, int label)
{
  int      xn, yn, zn, n ;
  GCA_NODE *gcan ;
  GC       *gc ;

  GCAvoxelToNode(gca, xt, yt, zt, &xn, &yn, &zn) ;

  if (xn == 13 && yn == 23 && zn == 33)
    DiagBreak() ;

  gcan = &gca->nodes[xn][yn][zn] ;

  gcan->total_training++ ;
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
      GC   *old_gcs ;

      old_max_labels = gcan->max_labels ; 
      gcan->max_labels += 2 ;
      old_labels = gcan->labels ;
      old_gcs = gcan->gcs ;

      /* allocate new ones */
      gcan->gcs = (GC *)calloc(gcan->max_labels, sizeof(GC)) ;
      if (!gcan->gcs)
        ErrorExit(ERROR_NOMEMORY, "GCANupdateNode: couldn't expand gcs to %d",
                  gcan->max_labels) ;
      gcan->labels = (char *)calloc(gcan->max_labels, sizeof(char)) ;
      if (!gcan->labels)
        ErrorExit(ERROR_NOMEMORY, "GCANupdateNode: couldn't expand labels to %d",
                  gcan->max_labels) ;

      /* copy the old ones over */
      memmove(gcan->gcs, old_gcs, old_max_labels*sizeof(GC)) ;
      memmove(gcan->labels, old_labels, old_max_labels*sizeof(char)) ;

      /* free the old ones */
      free(old_gcs) ; free(old_labels) ;
    }
    gcan->nlabels++ ;
  }

  gc = &gcan->gcs[n] ;

  /* these will be updated when training is complete */
  gc->mean += val ; gc->var += val*val ; gc->prior += 1.0f ;
  
  gcan->labels[n] = label ;

  return(NO_ERROR) ;
}
#define MIN_VAR  (5*5)   /* should make this configurable */
int
GCAcompleteTraining(GCA *gca)
{
  int      x, y, z, n, total_nodes, total_gcs ;
  float    nsamples ;
  GCA_NODE *gcan ;
  GC       *gc ;

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
          gc->mean /= nsamples ;
          var = gc->var / nsamples - gc->mean*gc->mean ;
          if (var < -0.1)
            DiagBreak() ;
          if (var < MIN_VAR)
            var = MIN_VAR ;
          gc->var = var ;
          gc->prior /= (float)gcan->total_training ;
        }
      }
    }
  }
  printf("%d classifiers: %2.1f per node\n", 
         total_gcs, (float)total_gcs/(float)total_nodes) ;
  return(NO_ERROR) ;
}

MRI  *
GCAlabel(MRI *mri_inputs, GCA *gca, MRI *mri_dst, LTA *lta)
{
  int      x, y, z, xt, yt, zt, width, height, depth, label, val,
           xn, yn, zn, n ;
  MATRIX   *m_L ;
  VECTOR   *v_input, *v_canon ;
  GCA_NODE *gcan ;
  GC       *gc ;
  float    dist, max_p, p ;

  if (!mri_dst)
  {
    mri_dst = MRIclone(mri_inputs, NULL) ;
    if (!mri_dst)
      ErrorExit(ERROR_NOMEMORY, "GCAlabel: could not allocate dst") ;
    MRIcopyHeader(mri_inputs, mri_dst) ;
  }

  v_input = VectorAlloc(4, MATRIX_REAL) ;
  v_canon = VectorAlloc(4, MATRIX_REAL) ;
  *MATRIX_RELT(v_input, 4, 1) = 1.0 ;
  *MATRIX_RELT(v_canon, 4, 1) = 1.0 ;

  /* convert transform to voxel coordinates */
#if 0
  m_L = MRIrasXformToVoxelXform(mri_inputs,mri_inputs,lta->xforms[0].m_L,NULL);
#else
  m_L = lta->xforms[0].m_L ;
#endif

  /* go through each voxel in the input volume and find the canonical
     voxel (and hence the classifier) to which it maps. Then update the
     classifiers statistics based on this voxel's intensity and label.
  */
  width = mri_inputs->width ; height = mri_inputs->height; 
  depth = mri_inputs->depth ;
  for (x = 0 ; x < width ; x++)
  {
    V3_X(v_input) = (float)x ;
    for (y = 0 ; y < height ; y++)
    {
      V3_Y(v_input) = (float)y ;
      for (z = 0 ; z < depth ; z++)
      {
        if (x == 97 && y == 97 && z == 128)
          DiagBreak() ; 

        V3_Z(v_input) = (float)z ;
        val = MRIvox(mri_inputs, x, y, z) ;

        /* compute coordinates in canonical space */
        MatrixMultiply(m_L, v_input, v_canon) ;
        xt = nint(V3_X(v_canon)) ;
        yt = nint(V3_Z(v_canon)) ;
        zt = nint(V3_Y(v_canon)) ;
        if (xt < 0) xt = 0 ;
        if (yt < 0) yt = 0 ;
        if (zt < 0) zt = 0 ;
        if (xt >= width)  xt = width-1 ;
        if (yt >= height) yt = height-1 ;
        if (zt >= depth)  zt = depth-1 ;

        /* find the node associated with this coordinate and classify */
        GCAvoxelToNode(gca, xt, yt, zt, &xn, &yn, &zn) ;
        gcan = &gca->nodes[xn][yn][zn] ;
        label = 0 ; max_p = -10000000.0 ;
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
        MRIvox(mri_dst, x, y, z) = label ;
      }
    }
  }

#if 0
  MatrixFree(&m_L) ;
#endif

  return(mri_dst) ;
}

