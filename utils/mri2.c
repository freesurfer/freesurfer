/**
 * @file  mri2.c
 * @brief more routines for loading, saving, and operating on MRI structures
 *
 */
/*
 * Original Author: Douglas N. Greve
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2012/06/11 17:46:18 $
 *    $Revision: 1.80 $
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
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <float.h>
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mrisurf.h"
#include "mriBSpline.h"
#include "fio.h"
#include "stats.h"
#include "corio.h"
#include "bfileio.h"
#include "proto.h"
#include "mri2.h"
#include "sig.h"
#include "cma.h"
#include "chronometer.h"

//#define MRI2_TIMERS


#include "affine.h"

#ifdef FS_CUDA
#include "mrivol2vol_cuda.h"
#endif

/* overwrite generic nint to speed up execution 
  Make sure NOT use calls like nint(f++) in this file!
*/
#define nint(f)  (f<0?((int)(f-0.5)):((int)(f+0.5)))

#define VERBOSE_MODE
#undef VERBOSE_MODE

/*-------------------------------------------------------------
  mri_load_bvolume() -- same as bf_ldvolume() but returns an
  MRI structure instead of a BF_DATA. See bfileio.h.
  -------------------------------------------------------------*/
MRI * mri_load_bvolume(char *bfstem)
{
  BF_DATA *bfvol;
  MRI *vol;
  int r,c,s,f;
  float val;

  /* first load as a BF_DATA sturcture */
  bfvol = bf_ldvolume(bfstem);
  if (bfvol == NULL) return(NULL);

  /* allocate the MRI */
  vol = MRIallocSequence(bfvol->ncols,bfvol->nrows,
                         bfvol->nslcs,MRI_FLOAT,bfvol->nfrms);
  if (vol==NULL)
  {
    bf_freebfd(&bfvol);
    fprintf(stderr,"mri_load_bvolume(): could not alloc vol\n");
    return(NULL);
  }

  /* copy data from the BF_DATA struct to the ARRAY4D*/
  for (r=0;r<bfvol->nrows;r++)
  {
    for (c=0;c<bfvol->ncols;c++)
    {
      for (s=0;s<bfvol->nslcs;s++)
      {
        for (f=0;f<bfvol->nfrms;f++)
        {
          val = BF_GETVAL(bfvol,r,c,s,f);
          MRIFseq_vox(vol,c,r,s,f) = val;
        }
      }
    }
  }

  bf_freebfd(&bfvol);
  return(vol);
}

/*-------------------------------------------------------------
  mri_save_as_bvolume() - same as bf_svvolume() but takes an MRI
  sturucture as input. See also bfileio.h.
  -------------------------------------------------------------*/
int mri_save_as_bvolume(MRI *vol, char *stem, int svendian, int svtype)
{
  BF_DATA *bfvol;
  int r,c,s,f;
  float val;

  /* allocate a temporary BF_DATA struct */
  bfvol = bf_allocbfd(vol->height, vol->width,  vol->depth, vol->nframes);
  if (bfvol == NULL) return(1);

  /* copy data from ARRAY4D to BF_DATA */
  for (r=0;r<bfvol->nrows;r++)
  {
    for (c=0;c<bfvol->ncols;c++)
    {
      for (s=0;s<bfvol->nslcs;s++)
      {
        for (f=0;f<bfvol->nfrms;f++)
        {
          val = MRIFseq_vox(vol,c,r,s,f);
          BF_SETVAL(val,bfvol,r,c,s,f);
        }
      }
    }
  }

  /* save the BF_DATA volume */
  bf_svvolume(bfvol,stem,svendian,svtype);

  bf_freebfd(&bfvol);
  return(0);
}
/*-------------------------------------------------------------
  mri_load_bvolume_frame() -- loads a single frame as an MRI.
  -------------------------------------------------------------*/
MRI * mri_load_bvolume_frame(char *bfstem, int frameno)
{
  BF_DATA *bfvol;
  MRI *vol;
  int r,c,s;
  float val;

  /* first load as a BF_DATA sturcture */
  bfvol = bf_ldvolume(bfstem);
  if (bfvol == NULL) return(NULL);

  if (frameno >= bfvol->nfrms)
  {
    fprintf(stderr,"ERROR: mri_load_bvolume_frame(): frameno = %d, exceeds "
            "number of frames in %s = %d\n",frameno,bfstem,bfvol->nfrms);
    bf_freebfd(&bfvol);
    return(NULL);
  }

  /* allocate the MRI */
  vol = MRIallocSequence(bfvol->ncols,bfvol->nrows,
                         bfvol->nslcs,MRI_FLOAT,1);
  if (vol==NULL)
  {
    bf_freebfd(&bfvol);
    fprintf(stderr,"mri_load_bvolume_frame(): could not alloc vol\n");
    return(NULL);
  }

  /* copy data from the BF_DATA struct to the ARRAY4D*/
  for (r=0;r<bfvol->nrows;r++)
  {
    for (c=0;c<bfvol->ncols;c++)
    {
      for (s=0;s<bfvol->nslcs;s++)
      {
        val = BF_GETVAL(bfvol,r,c,s,frameno);
        MRIFseq_vox(vol,c,r,s,0) = val;
      }
    }
  }

  bf_freebfd(&bfvol);
  return(vol);
}
/*-------------------------------------------------------------
  mri_save_as_cor() - basically the same as sv_cor() (corio.c) but
  takes an MRI structure as input.  The MRi has much more
  flexibility than the COR structure which is limited to 256^3, one frame,
  and type uchar.  If any MRI dimension exceeds 256, it is truncated;
  if it is less than 256, the extra is filled with zeros. If the MRI
  has more than one frame, the frame can be set by the frame argument.
  If the values in the MRI are beyond the range of 0-255 (ie, that
  of unsigned char), the MRI can be rescaled to 0-255 by setting
  the rescale flag to non-zero. See also corio.h.
  -------------------------------------------------------------*/
int mri_save_as_cor(MRI *vol,  char *cordir, int frame, int rescale)
{
  unsigned char **COR;
  int r,c,s;
  int rmax, cmax, smax;
  float val;

  if (frame >= vol->nframes)
  {
    fprintf(stderr,"mri_save_as_cor(): frame = %d, must be <= %d\n",
            frame,vol->nframes);
    return(1);
  }

  /* make sure the output directory is writable */
  if (! cordir_iswritable(cordir)) return(1);

  /* allocate a temporary COR volume */
  COR = alloc_cor();
  if (COR==NULL) return(1);

  /* rescale to 0-255 (range of uchar) */
  if (rescale) mri_rescale(vol,0,255,vol);

  /* make sure maximum subscript does not
     exceed either range */
  if (vol->height < 256) rmax = vol->height;
  else                  rmax = 256;
  if (vol->width < 256)  cmax = vol->width;
  else                  cmax = 256;
  if (vol->depth < 256)  smax = vol->depth;
  else                  smax = 256;

  /* copy data from ARRAY4D to COR */
  for (r=0; r < rmax ; r++)
  {
    for (c=0; c < cmax ; c++)
    {
      for (s=0; s < smax ; s++)
      {
        val = MRIFseq_vox(vol,c,r,s,0);
        CORVAL(COR,r,c,s) = (unsigned char) val;
      }
    }
  }

  /* save the COR volume */
  sv_cor(COR, cordir);

  free_cor(&COR);
  return(0);
}
/*------------------------------------------------------------
  mri_rescale() -- rescales array to be min <= val <= max. Uses
  outvol if non-null, otherwise allocates an output volume. Can
  be done in-place. Rescales across all frames.
  ------------------------------------------------------------*/
MRI *mri_rescale(MRI *vol, float min, float max, MRI *outvol)
{
  int r,c,s,f;
  float val, volmin, volmax, volrange, range;
  MRI *tmpvol;

  if (outvol != NULL) tmpvol = outvol;
  else
  {
    tmpvol = MRIallocSequence(vol->width,vol->height,vol->depth,
                              MRI_FLOAT,vol->nframes);
    if (tmpvol == NULL) return(NULL);
  }

  /* find the minimum and maximum */
  volmin = MRIFseq_vox(vol,0,0,0,0);
  volmax = MRIFseq_vox(vol,0,0,0,0);
  for (r=0;r<vol->height;r++)
  {
    for (c=0;c<vol->width;c++)
    {
      for (s=0;s<vol->depth;s++)
      {
        for (f=0;f<vol->nframes;f++)
        {
          val = MRIFseq_vox(vol,c,r,s,f);
          if (volmin > val) volmin = val;
          if (volmax < val) volmax = val;
        }
      }
    }
  }
  volrange = volmax - volmin;
  range    = max - min;

  /* rescale to the new range */
  for (r=0;r<vol->height;r++)
  {
    for (c=0;c<vol->width;c++)
    {
      for (s=0;s<vol->depth;s++)
      {
        for (f=0;f<vol->nframes;f++)
        {
          val = MRIFseq_vox(vol,c,r,s,f);
          val = range*val + min;
          MRIFseq_vox(tmpvol,c,r,s,f) = val;
        }
      }
    }
  }

  printf("volmin = %g, volmax = %g\n",volmin,volmax);

  return(tmpvol);
}
/*------------------------------------------------------------
  mri_minmax() -- gets min and max values of volume
  ------------------------------------------------------------*/
int mri_minmax(MRI *vol, float *min, float *max)
{
  int r,c,s,f;
  float val;

  *min = MRIFseq_vox(vol,0,0,0,0);
  *max = MRIFseq_vox(vol,0,0,0,0);
  for (r=0;r<vol->height;r++)
  {
    for (c=0;c<vol->width;c++)
    {
      for (s=0;s<vol->depth;s++)
      {
        for (f=0;f<vol->nframes;f++)
        {
          val = MRIFseq_vox(vol,c,r,s,f);
          if (*min > val) *min = val;
          if (*max < val) *max = val;
        }
      }
    }
  }
  printf("volmin = %g, volmax = %g\n",*min,*max);
  return(0);
}
/*--------------------------------------------------------
  mri_framepower() -- raises the value in each frame to that
  indicated by the framepower vector.  This is a pre-processing
  step for averaging across voxels when the data contain both
  averages and standard deviations. The stddevs need to be
  squared before spatial averaging and then square-rooted before
  saving again.
  -----------------------------------------------------------*/
int mri_framepower(MRI *vol, float *framepower)
{
  int r,c,s,f;
  float val;

  for (f=0;f<vol->nframes;f++)
  {

    /* power = 1 -- don't do anything */
    if ( fabs(framepower[f]-1.0) < .00001  )
      continue;

    /* power = 0.5 -- use sqrt() */
    if ( fabs(framepower[f]-0.5) < .00001  )
    {
      for (r=0;r<vol->height;r++)
      {
        for (c=0;c<vol->width;c++)
        {
          for (s=0;s<vol->depth;s++)
          {
            val = MRIFseq_vox(vol,c,r,s,f);
            MRIFseq_vox(vol,c,r,s,f) = sqrt(val);
          }
        }
      }
      continue;
    }

    /* power = 2 -- use val*val */
    if ( fabs(framepower[f]-0.5) < .00001  )
    {
      for (r=0;r<vol->height;r++)
      {
        for (c=0;c<vol->width;c++)
        {
          for (s=0;s<vol->depth;s++)
          {
            val = MRIFseq_vox(vol,c,r,s,f);
            MRIFseq_vox(vol,c,r,s,f) = val*val;
          }
        }
      }
      continue;
    }

    /* generic: use pow() -- least efficient */
    for (r=0;r<vol->height;r++)
    {
      for (c=0;c<vol->width;c++)
      {
        for (s=0;s<vol->depth;s++)
        {
          val = MRIFseq_vox(vol,c,r,s,f);
          MRIFseq_vox(vol,c,r,s,f) = pow(val,framepower[f]);
        }
      }
    }

  } /* end loop over frames */

  return(0);
}
/*-------------------------------------------------------------------
  mri_binarize() - converts each element to either 0 or 1 depending
  upon whether the value of the element is above or below the
  threshold. If the invert flag is set to 1, then the binarization is
  inverted. Can be done in-place (ie, vol=volbin). If volbin is NULL,
  then a new MRI vol is allocated and returned. tail is either
  positive, negative, or absolute (NULL=positive). nover is the
  number of voxels in vol that meet the threshold criteria
  (ie, the number of 1's in volbin).
  -----------------------------------------------------------------*/
MRI *mri_binarize(MRI *vol, float thresh, char *tail, int invert,
                  MRI *volbin, int *nover)
{
  int r,c,s,f, tailcode;
  float val;
  int b, onval, offval;
  MRI *voltmp;

  if (tail == NULL) tail = "positive";

  /* check the first 3 letters of tail */
  if (!strncasecmp(tail,"positive",3)) tailcode = 1;
  else if (!strncasecmp(tail,"negative",3)) tailcode = 2;
  else if (!strncasecmp(tail,"absolute",3)) tailcode = 3;
  else
  {
    fprintf(stderr,"mri_binarize: tail = %s unrecoginzed\n",tail);
    return(NULL);
  }

  if (volbin == NULL)
  {
    voltmp = MRIallocSequence(vol->width,vol->height,vol->depth,
                              MRI_FLOAT,vol->nframes);
    if (voltmp == NULL) return(NULL);
  }
  else voltmp = volbin;

  if (!invert)
  {
    onval  = 1;
    offval = 0;
    printf("NOT INVERTING\n");
  }
  else
  {
    onval  = 0;
    offval = 1;
    printf("INVERTING\n");
  }

  *nover = 0;
  for (r=0;r<vol->height;r++)
  {
    for (c=0;c<vol->width;c++)
    {
      for (s=0;s<vol->depth;s++)
      {
        for (f=0;f<vol->nframes;f++)
        {
          //val = MRIFseq_vox(vol,c,r,s,f);
          val = MRIgetVoxVal(vol,c,r,s,f);
          switch (tailcode)
          {
          case 2:
            val = -val;
            break;
          case 3:
            val = fabs(val);
            break;
          }
          if (val > thresh) b = onval;
          else             b = offval;
          if (b) (*nover) ++;
          //MRIFseq_vox(voltmp,c,r,s,f) = b;
          MRIsetVoxVal(voltmp,c,r,s,f,b);
        }
      }
    }
  }

  return(voltmp);
}
/*--------------------------------------------------------
  mri_load_cor_as_float()
  --------------------------------------------------------*/
MRI *mri_load_cor_as_float(char *cordir)
{
  MRI *ucvol;
  MRI *vol;
  int r,c,s;
  float val;

  /* read in the cor as unsigned char */
  ucvol = MRIread(cordir);
  if (ucvol == NULL) return(NULL);

  /* allocate a float volume */
  vol = MRIallocSequence(256,256,256,MRI_FLOAT,1);
  if (vol == NULL) return(NULL);

  for (r=0;r<vol->height;r++)
  {
    for (c=0;c<vol->width;c++)
    {
      for (s=0;s<vol->depth;s++)
      {
        val = (float)(MRIseq_vox(ucvol,c,r,s,0));
        MRIFseq_vox(vol,c,r,s,0) = val;
      }
    }
  }

  MRIfree(&ucvol);

  return(vol);
}
/* ---------------------------------------- */
/* not tested */
MRI *mri_load_wfile(char *wfile)
{
  FILE *fp;
  int i,ilat, num, vtx, nvertices;
  int *vtxnum;
  float *wval;
  MRI *w;

  fp = fopen(wfile,"r");
  if (fp==NULL)
  {
    fprintf(stderr,"ERROR: Progname: mri_load_wfile():\n");
    fprintf(stderr,"Could not open %s\n",wfile);
    fprintf(stderr,"(%s,%d,%s)\n",__FILE__, __LINE__,__DATE__);
    return(NULL);
  }

  fread2(&ilat,fp);
  fread3(&num,fp);

  vtxnum = (int *)   calloc(sizeof(int),   num);
  wval   = (float *) calloc(sizeof(float), num);

  for (i=0;i<num;i++)
  {
    fread3(&vtxnum[i],fp);
    wval[i] = freadFloat(fp) ;
  }
  fclose(fp);

  nvertices = vtxnum[num-1] + 1;

  w = MRIallocSequence(nvertices,1,1,MRI_FLOAT,1);
  for (i=0;i<num;i++)
  {
    vtx = vtxnum[i];
    MRIFseq_vox(w,vtx,0,0,0) = wval[i];
  }

  free(vtxnum);
  free(wval);
  return(w);
}
/*------------------------------------------------------------
  mri_sizeof() - returns the size of the data type of the MRI
  volume (in number of bytes).
  ------------------------------------------------------------*/
size_t mri_sizeof(MRI *vol)
{
  size_t bytes = 0;

  switch (vol->type)
  {
  case MRI_UCHAR:
    bytes = sizeof(BUFTYPE) ;
    break ;
  case MRI_SHORT:
    bytes = sizeof(short);
    break;
  case MRI_FLOAT:
    bytes = sizeof(float) ;
    break ;
  case MRI_INT:
    bytes = sizeof(int) ;
    break ;
  case MRI_LONG:
    bytes = sizeof(long) ;
    break ;
  }

  return(bytes);
}
/*------------------------------------------------------------
  mri_reshape() --
  ------------------------------------------------------------*/
MRI *mri_reshape(MRI *vol, int ncols, int nrows, int nslices, int nframes)
{
  MRI *outvol;
  int r,c,s,f, nv1, nv2;
  int r2,c2,s2,f2;

  if (vol->nframes == 0) vol->nframes = 1;

  nv1 = vol->width * vol->height * vol->depth * vol->nframes;
  nv2 = ncols*nrows*nslices*nframes;

  if (nv1 != nv2)
  {
    printf("ERROR: mri_reshape: number of elements cannot change\n");
    printf("  nv1 = %d, nv1 = %d\n",nv1,nv2);
    return(NULL);
  }

  outvol = MRIallocSequence(ncols,nrows,nslices,vol->type,nframes);
  if (outvol == NULL) return(NULL);

  MRIcopyHeader(vol, outvol); /* does not change dimensions */

  //printf("vol1: %d %d %d %d %d\n",vol->width,vol->height,
  // vol->depth,vol->nframes,mri_sizeof(vol));
  //printf("vol2: %d %d %d %d %d\n",outvol->width,outvol->height,
  // outvol->depth,outvol->nframes,mri_sizeof(outvol));

  c2 = 0;
  r2 = 0;
  s2 = 0;
  f2 = 0;
  for (f = 0; f < vol->nframes; f++)
  {
    for (s = 0; s < vol->depth;   s++)
    {
      for (r = 0; r < vol->height;  r++)
      {
        for (c = 0; c < vol->width;   c++)
        {

          switch (vol->type)
          {
          case MRI_UCHAR:
            MRIseq_vox(outvol,c2,r2,s2,f2) = MRIseq_vox(vol,c,r,s,f);
            break ;
          case MRI_SHORT:
            MRISseq_vox(outvol,c2,r2,s2,f2) = MRISseq_vox(vol,c,r,s,f);
            break;
          case MRI_FLOAT:
            MRIFseq_vox(outvol,c2,r2,s2,f2) = MRIFseq_vox(vol,c,r,s,f);
            break ;
          case MRI_INT:
            MRIIseq_vox(outvol,c2,r2,s2,f2) = MRIIseq_vox(vol,c,r,s,f);
            break ;
          case MRI_LONG:
            MRILseq_vox(outvol,c2,r2,s2,f2) = MRILseq_vox(vol,c,r,s,f);
            break ;
          }

          c2++;
          if (c2 == ncols)
          {
            c2 = 0;
            r2++;
            if (r2 == nrows)
            {
              r2 = 0;
              s2++;
              if (s2 == nslices)
              {
                s2 = 0;
                f2++;
              }
            }
          }
        }
      }
    }
  }

  return(outvol);
}

/*---------------------------------------------------------------
  MRIvol2Vol() - samples the values of one volume into that of
  another. Handles multiple frames. Can do nearest-neighbor,
  trilinear, and sinc interpolation (sinc may be flaky).
  Use this function instead of vol2vol_linear() in resample.c.

  Vt2s is the 4x4 matrix which converts CRS in the target to
  CRS in the source (ie, it is a vox2vox). If it is NULL, then
  the vox2vox is computed from the src and targ vox2ras matrices,
  assuming the src and targ share the same RAS.

  InterpCode is either: SAMPLE_NEAREST, SAMPLE_TRILINEAR, or
  SAMPLE_SINC.

  param is a generic parameter. For sinc, param is the hw parameter,
  otherwise, it currently has no meaning.
  ---------------------------------------------------------------*/
int MRIvol2Vol(MRI *src, MRI *targ, MATRIX *Vt2s,
               int InterpCode, float param)
{
#ifdef FS_CUDA
  int cudaReturn;
#else
  int   ct,  rt,  st,  f;
  int   ics, irs, iss;
  float fcs, frs, fss;
  double rval;
  float *valvect;
#endif
  int sinchw;
  MATRIX *V2Rsrc=NULL, *invV2Rsrc=NULL, *V2Rtarg=NULL;
  int FreeMats=0;
  MRI_BSPLINE * bspline = NULL;

#ifdef VERBOSE_MODE
  Chronometer tTotal, tSample;

  printf( "%s: Begin\n", __FUNCTION__ );

  printf( "Sources sizes are w=%i h=%i d=%i f=%i\n",
	  src->width, src->height, src->depth, src->nframes );
  printf( "src type is %i\n", src->type );

  printf( "Target sizes are w=%i h=%i d=%i f=%i\n",
	  src->width, src->height, src->depth, src->nframes );
  printf( "targ type is %i\n", targ->type );

  InitChronometer( &tTotal );
  InitChronometer( &tSample );
  StartChronometer( &tTotal );
#endif

  if (src->nframes != targ->nframes)
  {
    printf("ERROR: MRIvol2vol: source and target have different number "
           "of frames\n");
    return(1);
  }

  // Compute vox2vox matrix based on vox2ras of src and target.
  // Assumes that src and targ have same RAS space.
  if (Vt2s == NULL)
  {
    V2Rsrc = MRIxfmCRS2XYZ(src,0);
    invV2Rsrc = MatrixInverse(V2Rsrc,NULL);
    V2Rtarg = MRIxfmCRS2XYZ(targ,0);
    Vt2s = MatrixMultiply(invV2Rsrc,V2Rtarg,NULL);
    FreeMats = 1;
  }
  if (Gdiag_no > 0)
  {
    printf("MRIvol2Vol: Vt2s Matrix (%d)\n",FreeMats);
    MatrixPrint(stdout,Vt2s);
  }

  sinchw = nint(param);
#ifndef FS_CUDA
  valvect = (float *) calloc(sizeof(float),src->nframes);
#endif

#ifdef VERBOSE_MODE
  StartChronometer( &tSample );
#endif

#ifdef FS_CUDA
  cudaReturn = MRIvol2vol_cuda( src, targ,
				Vt2s,
				InterpCode,
				param );
  if( cudaReturn != 0 ) {
    fprintf( stderr, "%s: CUDA call failed!\n", __FUNCTION__ );
    exit( EXIT_FAILURE );
  }
#else

  if (InterpCode == SAMPLE_CUBIC_BSPLINE)
    bspline = MRItoBSpline(src,NULL,3);

  for (ct=0; ct < targ->width; ct++)
  {
    for (rt=0; rt < targ->height; rt++)
    {
      for (st=0; st < targ->depth; st++)
      {

        /* Column in source corresponding to CRS in Target */
        fcs = Vt2s->rptr[1][1] * ct + Vt2s->rptr[1][2] * rt +
              Vt2s->rptr[1][3] * st + Vt2s->rptr[1][4] ;
        ics = nint(fcs);
        if (ics < 0 || ics >= src->width) continue;

        /* Row in source corresponding to CRS in Target */
        frs = Vt2s->rptr[2][1] * ct + Vt2s->rptr[2][2] * rt +
              Vt2s->rptr[2][3] * st + Vt2s->rptr[2][4] ;
        irs = nint(frs);
        if (irs < 0 || irs >= src->height) continue;

        /* Slice in source corresponding to CRS in Target */
        fss = Vt2s->rptr[3][1] * ct + Vt2s->rptr[3][2] * rt +
              Vt2s->rptr[3][3] * st + Vt2s->rptr[3][4] ;
        iss = nint(fss);
        if (iss < 0 || iss >= src->depth) continue;

        /* Assign output volume values */
        if (InterpCode == SAMPLE_TRILINEAR)
          MRIsampleSeqVolume(src, fcs, frs, fss, valvect,
                             0, src->nframes-1) ;
        else
        {
          for (f=0; f < src->nframes ; f++)
          {
            switch (InterpCode)
            {
            case SAMPLE_NEAREST:
              valvect[f] = MRIgetVoxVal(src,ics,irs,iss,f);
              break ;
            case SAMPLE_CUBIC_BSPLINE:
              MRIsampleBSpline(bspline, fcs, frs, fss, f, &rval);
              valvect[f] = rval;
              break ;
            case SAMPLE_SINC:      /* no multi-frame */
              MRIsincSampleVolume(src, fcs, frs, fss, sinchw, &rval) ;
              valvect[f] = rval;
              break ;
            default:
              printf("ERROR: MRIvol2vol: interpolation method %i unknown\n",InterpCode);
              exit(1);

            }
          }
        }

        for (f=0; f < src->nframes; f++)
          MRIsetVoxVal(targ,ct,rt,st,f,valvect[f]);

      } /* target col */
    } /* target row */
    exec_progress_callback(ct, targ->width, 0, 1);
  } /* target slice */
#endif

#ifdef VERBOSE_MODE
  StopChronometer( &tSample );
#endif

#ifndef FS_CUDA
  free(valvect);
#endif
  if (FreeMats)
  {
    MatrixFree(&V2Rsrc);
    MatrixFree(&invV2Rsrc);
    MatrixFree(&V2Rtarg);
    MatrixFree(&Vt2s);
  }
  
  if (bspline) MRIfreeBSpline(&bspline);

#ifdef VERBOSE_MODE
  StopChronometer( &tTotal );


  printf( "Timings ------------\n" );
  printf( "  tSample : %9.3f ms\n", GetChronometerValue( &tSample ) );
  printf( "Total     : %9.3f ms\n", GetChronometerValue( &tTotal ) );
  printf( "%s: Done\n", __FUNCTION__ );
#endif

  return(0);
}
/*-------------------------------------------------------*/
/*!
  \fn int MRIvol2VolTkReg(MRI *mov, MRI *targ, MATRIX *Rtkreg,
      int InterpCode, float param)
  \brief Applies a tkregister matrix to a volume. Computes
    the vox2vox matrix then calls MRIvol2Vol().
  \param mov - source volume
  \param targ - target volume. Must be fully alloced.
  \param Rtkreg - tkreg ras2ras transform from target to source
    (assumes identity if NULL)
  \param InterCode - SAMPLE_NEAREST, SAMPLE_TRILINEAR, SAMPLE_SINC.
    Sinc probably does not work.
 */
int MRIvol2VolTkReg(MRI *mov, MRI *targ, MATRIX *Rtkreg,
                    int InterpCode, float param)
{
  MATRIX *vox2vox = NULL;
  MATRIX *Tmov, *invTmov, *Ttarg;
  int err;

  if (Rtkreg != NULL)
  {
    // TkReg Vox2RAS matrices
    Tmov      = MRIxfmCRS2XYZtkreg(mov);
    invTmov   = MatrixInverse(Tmov,NULL);
    Ttarg    = MRIxfmCRS2XYZtkreg(targ);
    // vox2vox = invTmov*R*Ttarg
    vox2vox = MatrixMultiply(invTmov,Rtkreg,vox2vox);
    MatrixMultiply(vox2vox,Ttarg,vox2vox);
  }
  else vox2vox = NULL;

  // resample
  err = MRIvol2Vol(mov,targ,vox2vox,InterpCode,param);

  if (vox2vox)
  {
    MatrixFree(&vox2vox);
    MatrixFree(&Tmov);
    MatrixFree(&invTmov);
    MatrixFree(&Ttarg);
  }

  return(err);
}

/*-----------------------------------------------------------*/
/*!
  \fn MRI *MRIvol2VolTLKernel(MRI *src, MRI *targ, MATRIX *Vt2s)
  \brief Computes the trilinear interpolation kernel at each voxel.
  \param src - source volume
  \param targ - target volume
  \param Vt2s - vox2vox transform from target to source (can be NULL)
 */
MRI *MRIvol2VolTLKernel(MRI *src, MRI *targ, MATRIX *Vt2s)
{
  int   ct,  rt,  st,  f;
  int   ics, irs, iss;
  float fcs, frs, fss;
  double *kvect=NULL;
  MATRIX *V2Rsrc=NULL, *invV2Rsrc=NULL, *V2Rtarg=NULL;
  int FreeMats=0;
  MRI *kernel;

  kernel = MRIallocSequence(targ->width, targ->height, targ->depth,
                            MRI_FLOAT, 8);
  if (kernel == NULL) return(NULL);
  MRIcopyHeader(targ,kernel);

  // Compute vox2vox matrix based on vox2ras of src and target.
  // Assumes that src and targ have same RAS space.
  if (Vt2s == NULL)
  {
    V2Rsrc = MRIxfmCRS2XYZ(src,0);
    invV2Rsrc = MatrixInverse(V2Rsrc,NULL);
    V2Rtarg = MRIxfmCRS2XYZ(targ,0);
    Vt2s = MatrixMultiply(invV2Rsrc,V2Rtarg,NULL);
    FreeMats = 1;
  }

  for (ct=0; ct < targ->width; ct++)
  {
    for (rt=0; rt < targ->height; rt++)
    {
      for (st=0; st < targ->depth; st++)
      {

        /* Column in source corresponding to CRS in Target */
        fcs = Vt2s->rptr[1][1] * ct + Vt2s->rptr[1][2] * rt +
              Vt2s->rptr[1][3] * st + Vt2s->rptr[1][4] ;
        ics = nint(fcs);
        if (ics < 0 || ics >= src->width) continue;

        /* Row in source corresponding to CRS in Target */
        frs = Vt2s->rptr[2][1] * ct + Vt2s->rptr[2][2] * rt +
              Vt2s->rptr[2][3] * st + Vt2s->rptr[2][4] ;
        irs = nint(frs);
        if (irs < 0 || irs >= src->height) continue;

        /* Slice in source corresponding to CRS in Target */
        fss = Vt2s->rptr[3][1] * ct + Vt2s->rptr[3][2] * rt +
              Vt2s->rptr[3][3] * st + Vt2s->rptr[3][4] ;
        iss = nint(fss);
        if (iss < 0 || iss >= src->depth) continue;

        kvect = MRItrilinKernel(src, fcs, frs, fss, kvect);
        for (f=0; f < 8 ; f++) MRIsetVoxVal(targ,ct,rt,st,f,kvect[f]);

      } /* target col */
    } /* target row */
  } /* target slice */

  free(kvect);
  if (FreeMats)
  {
    MatrixFree(&V2Rsrc);
    MatrixFree(&invV2Rsrc);
    MatrixFree(&V2Rtarg);
    MatrixFree(&Vt2s);
  }

  return(0);
}
/*---------------------------------------------------------------
  MRIdimMismatch() - checks whether the dimensions on two volumes
  are inconsistent. Returns:
    0 no mismatch
    1 columns mismatch
    2 rows mismatch
    3 slices mismatch
    4 frames mismatch - note: frameflag must = 1 to check frames
  ---------------------------------------------------------------*/
int MRIdimMismatch( const MRI *v1, const MRI *v2, int frameflag )
{
  if (v1->width  != v2->width)  return(1);
  if (v1->height != v2->height) return(2);
  if (v1->depth  != v2->depth)  return(3);
  if (frameflag && v1->nframes != v2->nframes)  return(4);
  return(0);
}

/*---------------------------------------------------------------
  MRIfdr2vwth() - computes the voxel-wise threshold needed to realize
  the given False Discovery Rate (FDR) based on the values in the
  given frame. Optionally thresholds the the MRI volume.

  frame - 0-based frame number of input volume to use
  fdr - false dicovery rate, between 0 and 1, eg: .05
  signid -
      0 = use all values regardless of sign
     +1 = use only positive values
     -1 = use only negative values
     If a vertex does not meet the sign criteria, its val2 is 0
  log10flag - interpret vol values as -log10(p)
  mask - binary mask volume. If the mask value at a voxel is 1,
     then its the input voxel value will be used to compute the
     threshold (if it also meets the sign criteria). If the mask is
     0, then ovol will be set to 0 at that location. Pass NULL if
     not using a mask.
  vwth - voxel-wise threshold between 0 and 1. If log10flag is set,
     then vwth = -log10(vwth). ovol voxels with p values
     GREATER than vwth will be 0. Note that this is the same
     as requiring -log10(p) > -log10(vwth).
  ovol - output volume. Pass NULL if no output volume is needed.

  So, for the ovol to be set to something non-zero, the vol must
  meet the sign, mask, and threshold criteria. If vol meets all
  the criteria, then ovol=vol (ie, no log10 transforms). vol
  itself is not altered, unless it is also passed as the ovol.
  Note: only frame=0 in the ovol is used.

  Return Values:
    0 - everything is OK
    1 - something went wrong

  Ref: http://www.sph.umich.edu/~nichols/FDR/FDR.m
  Thresholding of Statistical Maps in Functional Neuroimaging Using
  the False Discovery Rate.  Christopher R. Genovese, Nicole A. Lazar,
  Thomas E. Nichols (2002).  NeuroImage 15:870-878.

  *----------------------------------------------------*/
int MRIfdr2vwth(MRI *vol, int frame, double fdr, int signid,
                int log10flag, MRI *mask, double *vwth, MRI *ovol)
{
  double *p=NULL, val=0.0, valnull=0.0, maskval;
  int Nv, np, c, r ,s, maskflag=0;

  if (vol->nframes <= frame)
  {
    printf("ERROR: MRIfdr2vwth: frame = %d, must be < nframes = %d\n",
           frame,vol->nframes);
    return(1);
  }
  if (vol->type != MRI_FLOAT)
  {
    printf("ERROR: MRIfdr2vwth: input volume is not of type MRI_FLOAT\n");
    return(1);
  }
  if (ovol != NULL)
  {
    if (ovol->type != MRI_FLOAT)
    {
      printf("ERROR: MRIfdr2vwth: output volume is not of type MRI_FLOAT\n");
      return(1);
    }
    if (MRIdimMismatch(vol, ovol, 0))
    {
      printf("ERROR: MRIfdr2vwth: output/input dimension mismatch\n");
      return(1);
    }
  }
  if (mask != NULL)
  {
    if (MRIdimMismatch(vol, mask, 0))
    {
      printf("ERROR: MRIfdr2vwth: mask/input dimension mismatch\n");
      return(1);
    }
    maskflag = 1;
  }

  Nv = vol->width * vol->height * vol->depth;
  if (log10flag) valnull = 0;
  else          valnull = 1;

  p = (double *) calloc(Nv,sizeof(double));
  np = 0;
  for (c=0; c < vol->width; c++)
  {
    for (r=0; r < vol->height; r++)
    {
      for (s=0; s < vol->depth; s++)
      {

        if (maskflag)
        {
          // Must be in the mask if using a mask
          maskval = MRIgetVoxVal(mask,c,r,s,0);
          if (maskval < 0.5)  continue;
        }
        val = MRIFseq_vox(vol,c,r,s,frame);

        // Check the sign
        if (signid == -1 && val > 0) continue;
        if (signid == +1 && val < 0) continue;

        // Get value, convert from log10 if needed
        val = fabs(val);
        if (log10flag) val = pow(10,-val);

        p[np] = val;
        np++;
      }
    }
  }
  printf("MRIfdr2vwth: np = %d, nv = %d\n",np,Nv);

  // Check that something met the match criteria,
  // otherwise return an error
  if (np==0)
  {
    printf("ERROR: MRIfdr2vwth: no matching voxels found\n");
    free(p);
    return(1);
  }

  *vwth = fdr2vwth(p,np,fdr);
  free(p);

  if (ovol == NULL)
  {
    if (log10flag) *vwth = -log10(*vwth);
    return(0);
  }

  // Perform the thresholding
  for (c=0; c < vol->width; c++)
  {
    for (r=0; r < vol->height; r++)
    {
      for (s=0; s < vol->depth; s++)
      {

        if (maskflag)
        {
          maskval = MRIgetVoxVal(mask,c,r,s,0);
          if (maskval < 0.5)
          {
            // Set to null if out of mask
            MRIFseq_vox(ovol,c,r,s,0) = valnull;
            continue;
          }
        }

        val = MRIFseq_vox(vol,c,r,s,frame);

        if (signid == -1 && val > 0)
        {
          // Set to null if wrong sign
          MRIFseq_vox(ovol,c,r,s,0) = valnull;
          continue;
        }
        if (signid == +1 && val < 0)
        {
          // Set to null if wrong sign
          MRIFseq_vox(ovol,c,r,s,0) = valnull;
          continue;
        }

        val = fabs(val);
        if (log10flag) val = pow(10,-val);

        if (val > *vwth)
        {
          // Set to null if greather than thresh
          MRIFseq_vox(ovol,c,r,s,0) = valnull;
          continue;
        }

        // Otherwise, it meets all criteria, so
        // pass the original value through
        MRIFseq_vox(ovol,c,r,s,0) =
          MRIFseq_vox(vol,c,r,s,frame);
      }
    }
  }
  if (log10flag) *vwth = -log10(*vwth);

  return(0);
}


/*-------------------------------------------------------------------
  MRIcovarianceMatrix() - computes the cross-frame (temporal)
  covariance matrix. Returns M=D*D/Nv' where D is an Nt-by-Nv data
  set, where Nt is the number of frames/timepoints and Nv is the
  number of voxels/vertices which equals ncols*nrows*nslices in the
  mri strucutre. If mask is non-null, then voxels/vertices whose value
  in the mask are less than 0.5 are excluded from the computation of
  the covariance matrix, and Nv becomes the number points in the mask.
  ------------------------------------------------------------------*/
MATRIX *MRIcovarianceMatrix(MRI *mri, MRI *mask)
{
  int UseMask = 0, nmask,f1,f2;
  int r,c,s;
  double sum,v1,v2;
  MATRIX *M;

  // Handle masking
  if (mask != NULL)
  {
    // count number of points in the mask
    nmask = 0;
    for (c=0;c<mri->width;c++)
    {
      for (r=0;r<mri->height;r++)
      {
        for (s=0;s<mri->depth;s++)
        {
          if (MRIgetVoxVal(mask,c,r,s,0) > 0.5) nmask++;
        }
      }
    }
    //printf("Number of voxels in the mask %d\n",nmask);
    if (nmask == 0)
    {
      printf("ERROR: no voxels in mask\n");
      return(NULL);
    }
    UseMask = 1;
  }
  else
  {
    // Otherwise use all voxels/vertices
    nmask = mri->width * mri->height * mri->depth;
    UseMask = 0;
  }

  // Allocate the covariance matrix
  M = MatrixAlloc(mri->nframes,mri->nframes,MATRIX_REAL);

  // Compute the covariance matrix
  for (f1=0;f1<mri->nframes;f1++)
  {
    //printf("f1 = %d\n",f1);
    for (f2=f1;f2<mri->nframes;f2++)
    {
      sum = 0;
      for (c=0;c<mri->width;c++)
      {
        for (r=0;r<mri->height;r++)
        {
          for (s=0;s<mri->depth;s++)
          {
            if (UseMask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
            v1 = MRIgetVoxVal(mri,c,r,s,f1);
            v2 = MRIgetVoxVal(mri,c,r,s,f2);
            sum += (v1*v2);
            //printf("%g %g %g\n",v1,v2,sum);
          } //s
        } //r
      } //s
      M->rptr[f1+1][f2+1] = sum/nmask;
      M->rptr[f2+1][f1+1] = sum/nmask;
    } // f1
  } // f2

  //MatrixPrint(stdout,M);
  return(M);
}

/*--------------------------------------------------------------------
  MRIpca() - computes the SVD/PCA of the input volume D.  pU is the
  matrix of temporal EVs, pS is the vector singular values, and pV arg
  the spatial EVs (reshaped back into a volume).  Note: this uses the
  SVD code from matrix.c, which is not very accurate. If mask is
  non-null, uses only the voxels for which mask>0.5. When complete,
  D = U*S*V';
  -------------------------------------------------------------------*/
int MRIpca(MRI *D, MATRIX **pU, VECTOR **pS, MRI **pV, MRI *mask)
{
  int dim, dim_real, nvoxels, c,r,s,f, UseMask, nmask=0;
  MATRIX *M, *VV, *UinvS, *Fd, *Fv;
  VECTOR *S2;
  double *sum2, v;

  nvoxels = D->width * D->height * D->depth;
  dim = MIN(nvoxels,D->nframes);

  // Count the number of voxels in the mask
  if (mask != NULL)
  {
    UseMask = 1;
    for (c=0;c<D->width;c++)
      for (r=0;r<D->height;r++)
        for (s=0;s<D->depth;s++)
          if (UseMask && MRIgetVoxVal(mask,c,r,s,0) > 0.5) nmask++;
  }
  else
  {
    UseMask = 0;
    nmask = nvoxels;
  }

  M = MRIcovarianceMatrix(D, mask);
  if (M==NULL) return(1);
  //MatrixWriteTxt("cvm.dat",M);

  // Compute the SVD of the Temporal Cov Matrix
  S2 = RVectorAlloc(D->nframes, MATRIX_REAL) ;
  *pU = MatrixCopy(M,NULL); // It's done in-place so make a copy
  VV = MatrixSVD(*pU, S2, NULL); // M = U*S2*VV, VV = U';
  //MatrixWriteTxt("s2.dat",S2);

  *pS = RVectorAlloc(D->nframes, MATRIX_REAL) ;
  for (c=1; c <= dim; c++)
  {
    S2->rptr[1][c] *= nmask; // Needs this
    (*pS)->rptr[1][c] = sqrt(S2->rptr[1][c]);
    // Exclude dims for which the singular value is very small
    if (S2->rptr[1][c] < EPSILON) break;
  }
  dim_real = c-1;

  // Compute U*inv(S) (note square root to go from S2 to S)
  UinvS = MatrixAlloc(D->nframes,dim_real,MATRIX_REAL);
  for (c=1; c <= dim_real; c++)
  {
    for (f=1; f <= D->nframes; f++)
      UinvS->rptr[f][c] = (double)(*pU)->rptr[f][c]/sqrt(S2->rptr[1][c]);
  }
  //printf("UinvS %d -------------------------\n",dim_real);
  //MatrixPrint(stdout,UinvS);

  // Allocate the spatial EVs
  *pV = MRIallocSequence(D->width, D->height,D->depth,MRI_FLOAT, dim_real);
  MRIcopyHeader(D,*pV);

  // Compute V = D'*U*inv(S)
  sum2 = (double *) calloc(dim_real,sizeof(double));
  Fd = MatrixAlloc(1,D->nframes,MATRIX_REAL);
  Fv = MatrixAlloc(1,dim_real,MATRIX_REAL);
  for (c=0;c<D->width;c++)
  {
    for (r=0;r<D->height;r++)
    {
      for (s=0;s<D->depth;s++)
      {
        if (UseMask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
        for (f=0; f < D->nframes; f++)
          Fd->rptr[1][f+1] = MRIgetVoxVal(D,c,r,s,f);
        MatrixMultiply(Fd,UinvS,Fv);
        for (f=0; f < dim_real; f++)
        {
          MRIsetVoxVal(*pV,c,r,s,f,Fv->rptr[1][f+1]);
          sum2[f] += (Fv->rptr[1][f+1])*(Fv->rptr[1][f+1]);
        }
      }
    }
  }

  /* Normalize to a unit vector */
  for (c=0;c<D->width;c++)
  {
    for (r=0;r<D->height;r++)
    {
      for (s=0;s<D->depth;s++)
      {
        if (UseMask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
        for (f=0; f < dim_real; f++)
        {
          v = MRIgetVoxVal(*pV,c,r,s,f);
          MRIsetVoxVal(*pV,c,r,s,f,v/sqrt(sum2[f]));
        }
      }
    }
  }

  free(sum2);
  MatrixFree(&M);
  MatrixFree(&VV);
  MatrixFree(&UinvS);
  MatrixFree(&Fd);
  MatrixFree(&Fv);
  MatrixFree(&S2);

  return(0);
}

/*--------------------------------------------------------------
  PrintPCAStats() - prints pca summary statistics to a stream.
  The summary stats are: (1) nth EV (2) var spanned by nth,
  (3) var spanned by 1-nth EVs, (3) percent var spanned by nth,
  (4) perent var spanned by 1-nth EVs,
  --------------------------------------------------------------*/
int PrintPCAStats(FILE *fp, MATRIX *Spca)
{
  int n;
  double totvar,v,vsum;

  totvar=0.0;
  for (n=1; n <= Spca->cols; n++)
    totvar += (Spca->rptr[1][n] * Spca->rptr[1][n]);

  vsum = 0.0;
  for (n=1; n <= Spca->cols; n++)
  {
    v = (Spca->rptr[1][n] * Spca->rptr[1][n]);
    vsum += v;
    fprintf(fp,"%3d   %8.2f %9.2f   %6.2f  %6.2f \n",
            n,v,vsum,100*v/totvar,100*vsum/totvar);
  }
  return(0);
}

/*--------------------------------------------------------------
  WritePCAStats() - writes pca summary statistics to a file.
  see PrintPCAStats() for more info.
  --------------------------------------------------------------*/
int WritePCAStats(char *fname, MATRIX *Spca)
{
  FILE *fp;
  fp = fopen(fname,"w");
  if (fp == NULL)
  {
    printf("ERROR: opening %s\n",fname);
    return(1);
  }
  PrintPCAStats(fp, Spca);
  return(0);
}
/*---------------------------------------------------------------
  MRIsqrt() - computes sqrt(fabs(v)). Calls MRIsquraRoot().
  ---------------------------------------------------------------*/
MRI *MRIsqrt(MRI *invol, MRI *outvol)
{
  outvol = MRIsquareRoot(invol,NULL,outvol);
  return(outvol);
}
/*------------------------------------------------------*/
/*!
 \fn MRI *MRImax(MRI *mri1, MRI *mri2, MRI *out)
 \brief Computes the voxel-by-voxel max between the input MRIs.
 \param MRI *mri1 - first input
 \param MRI *mri2 - second input
 \param MRI *out  - output (can be NULL)
 \return MRI * - pointer to output MRI
 Both inputs must be of the same size, but they can be of
 different data types. If the output is not specified, then
 it gets the data type of mri1. If it is specified, it can
 be any data type.
*/
MRI *MRImax(MRI *mri1, MRI *mri2, MRI *out)
{
  int c, r, s, f, n, ncols, nrows, nslices,nframes;
  void   *pmri1=NULL,*pmri2=NULL, *pout=NULL;
  double  v1=0,v2=0,v;
  int sz1, sz2,szout;

  ncols   = mri1->width;
  nrows   = mri1->height;
  nslices = mri1->depth;
  nframes = mri1->nframes;

  if (out==NULL)
  {
    out = MRIallocSequence(ncols, nrows, nslices, mri1->type, nframes);
    if (out==NULL)
    {
      printf("ERROR: MRImax: could not alloc output\n");
      return(NULL);
    }
    MRIcopyHeader(mri1,out); // ordinarily would need to change nframes
  }
  if (out->width != ncols   || out->height != nrows ||
      out->depth != nslices || out->nframes != nframes)
  {
    printf("ERROR: MRImax: dimension mismatch\n");
    return(NULL);
  }

  // Number of bytes in the mri data types
  sz1   = MRIsizeof(mri1->type);
  sz2   = MRIsizeof(mri2->type);
  szout = MRIsizeof(out->type);

  n = 0;
  for (f=0; f<nframes; f++)
  {
    for (s=0; s<nslices; s++)
    {
      for (r=0; r<nrows; r++)
      {
        // Pointers to the start of the column
        pmri1  = (void *) mri1->slices[n][r];
        pmri2  = (void *) mri2->slices[n][r];
        pout   = (void *) out->slices[n][r];
        for (c=0; c<ncols; c++)
        {

          v1 = MRIptr2dbl(pmri1, mri1->type);
          v2 = MRIptr2dbl(pmri2, mri2->type);
          if (v1 > v2) v = v1;
          else        v = v2;
          MRIdbl2ptr(v, pout, out->type);

          pmri1 += sz1;
          pmri2 += sz2;
          pout  += szout;
        } // cols
      } // rows
      n++;
    } // slices
  } // frames

  return(out);
}


/*---------------------------------------------------------------
  MRImaxAbsDiff() - finds the voxel where the two volumes differ
  the most.
  ---------------------------------------------------------------*/
double MRImaxAbsDiff(MRI *vol1, MRI *vol2,
                     int *cmax, int *rmax, int *smax, int *fmax)
{
  int c,r,s,f;
  double v1,v2,maxdiff;

  maxdiff = 0.0;
  for (c=0;c<vol1->width;c++)
  {
    for (r=0;r<vol1->height;r++)
    {
      for (s=0;s<vol1->depth;s++)
      {
        for (f=0; f < vol1->nframes; f++)
        {

          v1 = MRIgetVoxVal(vol1,c,r,s,f);
          v2 = MRIgetVoxVal(vol2,c,r,s,f);
          if (maxdiff < fabs(v1-v2))
          {
            maxdiff = fabs(v1-v2);
            *cmax = c;
            *rmax = r;
            *smax = s;
            *fmax = f;
          }

        }
      }
    }
  }
  return(maxdiff);
}
/* --------------------------------------------------------------- */
MRI *MRImultiplyConst(MRI *src, double vconst, MRI *dst)
{
  int r,c,s,f;
  double v;

  if (dst==NULL)
  {
    dst = MRIallocSequence(src->width,src->height,src->depth,
                           MRI_FLOAT,src->nframes) ;
    MRIcopyHeader(src, dst);
  }

  for (c=0; c < src->width; c++)
  {
    for (r=0; r < src->height; r++)
    {
      for (s=0; s < src->depth; s++)
      {
        for (f=0; f < src->nframes; f++)
        {
          v = MRIgetVoxVal(src,c,r,s,f);
          MRIsetVoxVal(dst,c,r,s,f,v*vconst);
        }
      }
    }
  }

  return(dst);
}
/* --------------------------------------------------------------- */
MRI *MRIaddConst(MRI *src, double vconst, MRI *dst)
{
  int r,c,s,f;
  double v;

  if (dst==NULL)
  {
    dst = MRIallocSequence(src->width,src->height,src->depth,
                           MRI_FLOAT,src->nframes) ;
    MRIcopyHeader(src, dst);
  }

  for (c=0; c < src->width; c++)
  {
    for (r=0; r < src->height; r++)
    {
      for (s=0; s < src->depth; s++)
      {
        for (f=0; f < src->nframes; f++)
        {
          v = MRIgetVoxVal(src,c,r,s,f);
          MRIsetVoxVal(dst,c,r,s,f,v+vconst);
        }
      }
    }
  }

  return(dst);
}

/*--------------------------------------------------------------------
  MRIframeBinarize() - creates a binary mask of voxels for which the
  abs(mri->val) of all frames are > thresh. If input mask is not null,
  then that mask is pruned, ie, if a voxel was not in the input mask,
  it will not be in the return mask. If it was in the input mask but
  it does not meet the frame criteria in mri, then it will be set to
  0. Note: if mask is not null, it's values will be changed.
  --------------------------------------------------------------------*/
MRI *MRIframeBinarize(MRI *mri, double thresh, MRI *mask)
{
  int c,r,s,f,n,premask;
  double val,m;

  premask = 1;
  if (!mask)
  {
    mask = MRIcloneBySpace(mri,MRI_FLOAT,1);
    MRIclear(mask);
    premask = 0;
  }

  for (c=0; c < mri->width; c++)
  {
    for (r=0; r < mri->height; r++)
    {
      for (s=0; s < mri->depth; s++)
      {
        if (premask)
        {
          m = MRIgetVoxVal(mask,c,r,s,0);
          if (m < 0.5) continue;
        }
        n = 0;
        for (f=0; f < mri->nframes; f++)
        {
          val = MRIgetVoxVal(mri,c,r,s,f);
          if (fabs(val) > thresh) n++;
        }
        if (n == mri->nframes) MRIsetVoxVal(mask,c,r,s,0,1);
        else                  MRIsetVoxVal(mask,c,r,s,0,0);
      }
    }
  }
  return(mask);
}
/*!
  \fn MRI *MRIexp(MRI *mri, double a, double b, MRI *mask, MRI *out)
  \brief Computes a*exp(b*mri). If a mask is supplied, then values
         outside the mask (ie, mask < 0.5) are set to 0 (if out=NULL)
         or the previous value of out.
 */
MRI *MRIexp(MRI *mri, double a, double b, MRI *mask, MRI *out)
{
  int c, r, s, f;
  double val, valout, m;
  int err;

  if (out==NULL)
  {
    out = MRIcloneBySpace(mri,MRI_FLOAT, -1);
    if (out==NULL)
    {
      printf("ERROR: MRIexp: could not alloc\n");
      return(NULL);
    }
  }
  else
  {
    err = MRIdimMismatch(mri, out, 1);
    if (err)
    {
      printf("ERROR: MRIexp(): output dimension mismatch (%d)\n",err);
      return(NULL);
    }
    if (out->type != MRI_FLOAT)
    {
      printf("ERROR: MRIexp(): structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }

  for (c=0; c < mri->width; c++)
  {
    for (r=0; r < mri->height; r++)
    {
      for (s=0; s < mri->depth; s++)
      {
        if (mask)
        {
          m = MRIgetVoxVal(mask,c,r,s,0);
          if (m < 0.5) continue;
        }
        for (f=0; f < mri->nframes; f++)
        {
          val = MRIgetVoxVal(mri,c,r,s,f);
          valout = a*exp(b*val);
          MRIsetVoxVal(out,c,r,s,f,valout);
        }
      }
    }
  }
  return(out);
}

/*!
  \fn MRI *MRIsum(MRI *mri1,MRI *mri2, double a,double b, MRI *mask,MRI *out)
  \brief Computes a*mri1 + b*mri2. If a mask is supplied, then values
         outside the mask (ie, mask < 0.5) are set to 0 (if out=NULL)
         or the previous value of out.
 */
MRI *MRIsum(MRI *mri1, MRI *mri2, double a, double b, MRI *mask, MRI *out)
{
  int c, r, s, f;
  double val1, val2, valout, m;
  int err;

  err = MRIdimMismatch(mri1, mri2, 1);
  if (err)
  {
    printf("ERROR: MRIsum(): input dimension mismatch (%d)\n",err);
    return(NULL);
  }

  if (out==NULL)
  {
    out = MRIcloneBySpace(mri1,MRI_FLOAT, -1);
    if (out==NULL)
    {
      printf("ERROR: MRIsum: could not alloc\n");
      return(NULL);
    }
  }
  else
  {
    err = MRIdimMismatch(mri1, out, 1);
    if (err)
    {
      printf("ERROR: MRIsum(): output dimension mismatch (%d)\n",err);
      return(NULL);
    }
  }

  for (c=0; c < mri1->width; c++)
  {
    for (r=0; r < mri1->height; r++)
    {
      for (s=0; s < mri1->depth; s++)
      {
        if (mask)
        {
          m = MRIgetVoxVal(mask,c,r,s,0);
          if (m < 0.5) continue;
        }
        for (f=0; f < mri1->nframes; f++)
        {
          val1 = MRIgetVoxVal(mri1,c,r,s,f);
          val2 = MRIgetVoxVal(mri2,c,r,s,f);
          valout = a*val1 + b*val2;
          MRIsetVoxVal(out,c,r,s,f,valout);
        }
      }
    }
  }
  return(out);
}
/*!
  \fn MRI *MRIvote(MRI *in, MRI *mask, MRI *vote)
  \brief Select the most frequently occuring value measured
     across frames in each voxel.
  \param vote - has 2 frames:
     (1) Most freqently occuring value
     (2) Fraction of occurances (ie noccurances/nframes)
  Note: input will be sorted in asc order
 */
MRI *MRIvote(MRI *in, MRI *mask, MRI *vote)
{
  int c, r, s, f, f0, ncols, nrows, nslices,nframes;
  float m;
  double vmax,v,v0;
  int runlen, runlenmax;
  MRI *sorted;

  printf("MRIvote: sorting\n");
  sorted = MRIsort(in,mask,in); // this sorts the input
  if (sorted == NULL) return(NULL);
  printf("MRIvote: done sorting\n");

  ncols   = in->width;
  nrows   = in->height;
  nslices = in->depth;
  nframes = in->nframes;

  if (vote==NULL)
  {
    vote = MRIallocSequence(ncols, nrows, nslices, in->type, 2);
    if (vote==NULL)
    {
      printf("ERROR: MRIvote: could not alloc\n");
      return(NULL);
    }
    MRIcopyHeader(in,vote);
    vote->nframes = 2;
  }
  if (in->type != vote->type)
  {
    printf("ERROR: MRIvote: type mismatch\n");
    return(NULL);
  }
  if (vote->width != ncols   || vote->height  != nrows ||
      vote->depth != nslices || vote->nframes != 2)
  {
    printf("ERROR: MRIvote: dimension mismatch\n");
    return(NULL);
  }

  for (c=0; c<ncols; c++)
  {
    for (r=0; r<nrows; r++)
    {
      for (s=0; s<nslices; s++)
      {
        if (mask)
        {
          m = MRIgetVoxVal(mask,c,r,s,0);
          if (m < 0.5) continue;
        }
        vmax = 0;
        runlenmax = 0;
        v0 = MRIgetVoxVal(sorted,c,r,s,0); // value at start of run
        f0 = 0;                            // frame at start of run
        f = 1;
        while (f < nframes)
        {
          v = MRIgetVoxVal(sorted,c,r,s,f);
          if (v0 != v)
          {
            // new value is different than that of run start
            runlen = f - f0; // runlength for v0
            if (runlenmax < runlen)
            {
              runlenmax = runlen;
              vmax = v0;
            }
            v0 = v;
            f0 = f;
          }
          f++;
        }
        // Need to do this one more time in case last value
        // has the longest run
        runlen = f - f0;
        if (runlenmax < runlen)
        {
          runlenmax = runlen;
          vmax = v0;
          v0 = v;
          f0 = f;
        }
        MRIsetVoxVal(vote,c,r,s,0,vmax);
        MRIsetVoxVal(vote,c,r,s,1,(double)runlenmax/nframes);

      } // slices
    } // rows
  } // cols

  return(vote);
}



/* MRImakeVox2VoxReg() - takes a target volume, a movable volume, a
   registration method, and an optional registration filename, and
   will make an mriTransform in which A is the target non-comformed
   index space, ARAS is the target tkRegRAS space, B is the movable
   tkRegRAS space (registered), and B is the movable non-conformed
   index space. If *transform is NULL, this will make a new
   mriTransform and return it, otherwise it will set the values in the
   *transform passed in.

   In the use case of a functional volume "being registered onto" an
   anatomical volume, the anatomical volume should be the targ
   parameter, and the functional should be the mov parameter.
 */
int MRImakeVox2VoxReg(MRI* targ, MRI* mov,
                      int regtype, char* regname,
                      mriTransformRef* transform)
{
  int retcode=0;
  char *cur_char, *base_end;
  int err;
  MATRIX *targ_idx_to_tkregras=NULL;
  MATRIX *mov_idx_to_tkregras=NULL;
  struct stat file_info;
  fMRI_REG *reg_info=NULL;
  char regpath[1000];
  char fullregname[1000];
  MATRIX *targ_tkregras_to_mov_tkregras=NULL;
  Trns_tErr trnscode;

  if (NULL==targ)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,"MRImakeVox2VoxReg: targ was NULL"));

  if (NULL==mov)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,"MRImakeVox2VoxReg: mov was NULL"));

  if (regtype < VOX2VOXREGTYPE_FILE || regtype > VOX2VOXREGTYPE_IDENTITY)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,"MRImakeVox2VoxReg: invalid reg type"));

  if (VOX2VOXREGTYPE_FILE==regtype &&
      NULL==regname)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,"MRImakeVox2VoxReg: reg type was FILE but "
                 "regname was NULL"));

  /* We need to build three matrices to fill in the three transforms
     in the mriTransform object:

     targ (A) ----------> tkRegRAS unregistered (ARAS)
                             |
                             |
                             V
      mov (B) ----------> tkRegRAS registered (BRAS)

     A->ARAS and B->BRAS take each volume from native index space to
     tkRegRAS space, an RAS space in which the center of the world is
     the center of the volume. We use MRIxfmXRS2XYZtkreg to get these.

     For the ARAS->BRAS matrix, we're either going to use a matrix
     from a registration file or one generated by MRItkRegMtx.
  */

  /* Create the targ A->RAS matrix. */
  targ_idx_to_tkregras = MRIxfmCRS2XYZtkreg(targ);
  if (NULL==targ_idx_to_tkregras)
  {
    print("ERROR: MRImakeVox2VoxReg: Couldn't create targ_idx_to_tkregras\n");
    goto error;
  }

  /* Create the mov B->RAS matrix. */
  mov_idx_to_tkregras = MRIxfmCRS2XYZtkreg(mov);
  if (NULL==mov_idx_to_tkregras)
  {
    print("ERROR: MRImakeVox2VoxReg: Couldn't create mov_idx_to_tkregras\n");
    goto error;
  }

  /* Now we build the ARAS->BRAS matrix, which is the
     registration. Switch on our registration type. */
  switch (regtype)
  {
  case VOX2VOXREGTYPE_FILE:
  case VOX2VOXREGTYPE_FIND:

    /* If we're reading a file, copy the file from the input or
    generate one from our data file location. */
    if (VOX2VOXREGTYPE_FILE==regtype)
    {
      strncpy(fullregname, regname, sizeof(fullregname));
    }
    else if (VOX2VOXREGTYPE_FIND==regtype)
    {
      /* Copy the movable volume name and find the last / in the
         file name. From there, copy in "register.dat" for our file
         name. */
      strncpy(regpath, mov->fname, sizeof(regpath));
      cur_char = regpath;
      base_end = regpath;
      while (NULL!=cur_char && '\0' != *cur_char)
      {
        if ('/' == *cur_char)
          base_end = cur_char;
        cur_char++;
      }
      *base_end = '\0';
      snprintf(fullregname,sizeof(fullregname),
               "%s/%s",regpath,"register.dat");
    }

    /* Check that the file exists. */
    err = stat(fullregname,&file_info);
    if (0!=err)
    {
      printf("ERROR: MRImakeVox2VoxReg: Couldn't find registration %s\n",
             fullregname);
      goto error;
    }

    /* Check if it's a regular file. */
    if (!S_ISREG(file_info.st_mode))
    {
      printf("ERROR: MRImakeVox2VoxReg: Couldn't open registration %s\n",
             fullregname);
      goto error;
    }

    /* Read the registration */
    reg_info = StatReadRegistration(fullregname);
    if (NULL==reg_info)
    {
      printf("ERROR: MRImakeVox2VoxReg: Problem reading registration %s\n",
             fullregname);
      goto error;
    }

    /* Copy the registration matrix. */
    targ_tkregras_to_mov_tkregras = MatrixCopy(reg_info->mri2fmri,NULL);

    break;

  case VOX2VOXREGTYPE_IDENTITY:

    /* Use MRItkRegMtx to generate an identity registration between
    the two volumes. */
    targ_tkregras_to_mov_tkregras = MRItkRegMtx(targ,mov,NULL);

    break;
  }

  /* Now look at *transform and create a new one if it doesn't
     exist. */
  if (*transform==NULL)
  {
    trnscode=Trns_New(transform);
    if (Trns_tErr_NoErr!=trnscode)
    {
      printf("ERROR: MRImakeVox2VoxReg: Error creating mriTransform\n");
      goto error;
    }
  }

#if 0
  printf("targ_idx_to_tkregras:\n");
  MatrixPrint(stdout,targ_idx_to_tkregras);
  printf("mov_idx_to_tkregras:\n");
  MatrixPrint(stdout,mov_idx_to_tkregras);
  printf("targ_tkregras_to_mov_tkregras:\n");
  MatrixPrint(stdout,targ_tkregras_to_mov_tkregras);
#endif

  /* Copy the matrices in. */
  trnscode=Trns_CopyAtoRAS(*transform,targ_idx_to_tkregras);
  if (Trns_tErr_NoErr!=trnscode)
  {
    printf("ERROR: MRImakeVox2VoxReg: Error copying A to RAS\n");
    goto error;
  }
  trnscode=Trns_CopyBtoRAS(*transform,mov_idx_to_tkregras);
  if (Trns_tErr_NoErr!=trnscode)
  {
    printf("ERROR: MRImakeVox2VoxReg: Error copying B to RAS\n");
    goto error;
  }
  trnscode=Trns_CopyARAStoBRAS(*transform,targ_tkregras_to_mov_tkregras);
  if (Trns_tErr_NoErr!=trnscode)
  {
    printf("ERROR: MRImakeVox2VoxReg: Error copying ARAS to BRAS\n");
    goto error;
  }

  goto cleanup;
error:

  if (0==retcode)
    retcode = -1;

cleanup:

  if (NULL!=targ_idx_to_tkregras)
    MatrixFree(&targ_idx_to_tkregras);

  if (NULL!=mov_idx_to_tkregras)
    MatrixFree(&mov_idx_to_tkregras);

  if (NULL!=reg_info)
    StatFreeRegistration(&reg_info);

  if (NULL!=targ_tkregras_to_mov_tkregras)
    MatrixFree(&targ_tkregras_to_mov_tkregras);

  return(retcode);
}
/*!
  \fn double MRIsum2All(MRI *mri)
  \brief squares and then sums all the voxels.
*/
double MRIsum2All(MRI *mri)
{
  int c,r,s,f;
  double sum2all,val;

  sum2all = 0;
  for (c=0; c < mri->width; c++)
  {
    for (r=0; r < mri->height; r++)
    {
      for (s=0; s < mri->depth; s++)
      {
        for (f=0; f < mri->nframes; f++)
        {
          val = MRIgetVoxVal(mri,c,r,s,f);
          sum2all += (val*val);
        }
      }
    }
  }
  return(sum2all);
}
/*!
  \fn MRI *MRIsquare(MRI *in, MRI *mask, MRI *out)
  \brief Squares the value at each voxel. Values outside
  of the mask are set to 0. mask can be NULL.
*/
MRI *MRIsquare(MRI *in, MRI *mask, MRI *out)
{
  int c,r,s,f;
  double val, mval;

  if (out == NULL)  out = MRIclone(in, NULL);

  mval = 1;
  for (c=0; c < in->width; c++)
  {
    for (r=0; r < in->height; r++)
    {
      for (s=0; s < in->depth; s++)
      {
        if (mask) mval = MRIgetVoxVal(mask,c,r,s,0);
        if (mask)
        {
          val = MRIgetVoxVal(mask,c,r,s,0);
          if (val < 0.5) continue;
        }
        for (f=0; f < in->nframes; f++)
        {
          if (mval > 0.5) val = MRIgetVoxVal(in,c,r,s,f);
          else val = 0.0;
          MRIsetVoxVal(out,c,r,s,f,val*val);
        }
      }
    }
  }
  return(out);
}
/*!
  \fn MRI *MRIsquareRoot(MRI *in, MRI *mask, MRI *out)
  \brief Square root of the fabs(value) at each voxel.
  Values outside of the mask are set to 0.  mask can be
  NULL.
*/
MRI *MRIsquareRoot(MRI *in, MRI *mask, MRI *out)
{
  int c,r,s,f;
  double val, mval;

  if (out == NULL)
  {
    out = MRIallocSequence(in->width, in->height,
                           in->depth,MRI_FLOAT, in->nframes);
    MRIcopyHeader(in,out);
    out->type = MRI_FLOAT;
  }

  mval = 1;
  for (c=0; c < in->width; c++)
  {
    for (r=0; r < in->height; r++)
    {
      for (s=0; s < in->depth; s++)
      {
        if (mask) mval = MRIgetVoxVal(mask,c,r,s,0);
        for (f=0; f < in->nframes; f++)
        {
          if (mval > 0.5)
          {
            val = MRIgetVoxVal(in,c,r,s,f);
          }
          else val = 0.0;
          MRIsetVoxVal(out,c,r,s,f,sqrt(fabs(val)));
        }
      }
    }
  }
  return(out);
}
/*---------------------------------------------------------------*/
/*!
  \fn MRI *MRIchecker(MRI *mri, MRI *checker)
  \brief Creates a checkerboard pattern. Adjacent columns differ
  by 1. Adjacent rows differ by 2. Adjacent slices differ by a
  sign. Adjacent frames differ by a factor of 10. This is for
  visualizing the voxel sampling. Set fthresh = .1, fmax=3.5,
  and use linear blending.
*/
MRI *MRIchecker(MRI *mri, MRI *checker)
{
  int c,r,s,f;
  double cval=0, rval=0, sval=0, fval=0;

  if (checker == NULL)
  {
    checker = MRIallocSequence(mri->width,mri->height,mri->depth,
                               MRI_FLOAT,mri->nframes);
    if (checker == NULL) return(NULL);
    MRIcopyHeader(mri,checker);
    checker->type = MRI_FLOAT;
  }

  // Fill even cols with
  for (c=0; c < mri->width; c++)
  {
    if (c%2 == 0) cval = 1;
    else         cval = 2;
    for (r=0; r < mri->height; r++)
    {
      if (r%2 == 0) rval = cval;
      else         rval = cval+2;
      for (s=0; s < mri->depth; s++)
      {
        if (s%2 == 0) sval = +rval;
        else         sval = -rval;
        for (f=0; f < mri->nframes; f++)
        {
          if (f%2 == 0) fval = sval;
          else         fval = 10*sval;
          //printf("%d %d %d %d  %g\n",c,r,s,f,fval);
          MRIsetVoxVal(checker,c,r,s,f,fval);
        }
      }
    }
  }
  return(checker);
}

/*-----------------------------------------------------------*/
/*!
  \fn MRI *MRIvol2VolDelta(MRI *mov, MRI *targ, MATRIX *Rt2s)
  \brief Computes the amount by which each voxel moves
  \param mov - source volume
  \param targ - target volume
  \param Rt2s - ras2ras transform from target to source (can be NULL, but why?)
 */
MRI *MRIvol2VolDelta(MRI *mov, MRI *targ, MATRIX *Rt2s)
{
  int   ct,  rt,  st;
  double dx, dy, dz;
  MATRIX *targCRS,*targRAS=NULL, *movRAS=NULL, *targVox2RAS, *targVox2movRAS;
  int FreeMats=0;
  MRI *delta;

  delta = MRIallocSequence(targ->width,
                           targ->height,
                           targ->depth,
                           MRI_FLOAT, 3);
  if (delta == NULL) return(NULL);
  MRIcopyHeader(targ,delta);

  // Compute ras2ras matrix based on vox2ras of mov and target.
  // Assumes that mov and targ have same RAS space.
  if (Rt2s == NULL)
  {
    Rt2s = MRItkRegMtx(targ,mov,NULL);
    FreeMats = 1;
  }
  targVox2RAS = MRIxfmCRS2XYZtkreg(targ);
  targVox2movRAS = MatrixMultiply(Rt2s,targVox2RAS,NULL);

  targCRS = MatrixAlloc(4,1,MATRIX_REAL);
  targCRS->rptr[4][1] = 1;

  for (ct=0; ct < targ->width; ct++)
  {
    for (rt=0; rt < targ->height; rt++)
    {
      for (st=0; st < targ->depth; st++)
      {
        targCRS->rptr[1][1] = ct;
        targCRS->rptr[2][1] = rt;
        targCRS->rptr[3][1] = st;
        targRAS = MatrixMultiply(targVox2RAS,   targCRS,targRAS);
        movRAS  = MatrixMultiply(targVox2movRAS,targCRS,movRAS);

        dx = targRAS->rptr[1][1] - movRAS->rptr[1][1];
        dy = targRAS->rptr[2][1] - movRAS->rptr[2][1];
        dz = targRAS->rptr[3][1] - movRAS->rptr[3][1];

        MRIsetVoxVal(delta,ct,rt,st,0,dx);
        MRIsetVoxVal(delta,ct,rt,st,1,dy);
        MRIsetVoxVal(delta,ct,rt,st,2,dz);

      } /* target col */
    } /* target row */
  } /* target slice */

  if (FreeMats) MatrixFree(&Rt2s);
  MatrixFree(&targCRS);
  MatrixFree(&targRAS);
  MatrixFree(&movRAS);
  MatrixFree(&targVox2RAS);
  MatrixFree(&targVox2movRAS);

  return(delta);
}

/*-----------------------------------------------------------*/
/*!
  \fn MRI *MRIcrop(MRI *mri, int c1, int r1, int s1, int c2, int r2, int s2)
  \brief Crop volume, keep correct goemetry.
  \param mri - volume to be cropped
  \param crs1 - col, row, slice to start cropping (inclusive)
  \param crs2 - col, row, slice to end cropping (inclusive)
 */
MRI *MRIcrop(MRI *mri, int c1, int r1, int s1, int c2, int r2, int s2)
{
  int c, r, s, f, Nc, Nr, Ns;
  MRI *crop;
  MATRIX *Vox2RAS, *crs, *P0;
  double v;

  if (c1 < 0 || c1 >= mri->width ||
      r1 < 0 || r1 >= mri->height ||
      s1 < 0 || s1 >= mri->depth)
  {
    printf("MRIcrop(): start point %d %d %d out of range\n",c1,r1,s1);
    return(NULL);
  }

  if (c2 < 0 || c2 >= mri->width ||
      r2 < 0 || r2 >= mri->height ||
      s2 < 0 || s2 >= mri->depth)
  {
    printf("MRIcrop(): end point %d %d %d out of range\n",c2,r2,s2);
    return(NULL);
  }

  // Size of cropped volume. +1 to make inclusive.
  Nc = c2-c1+1;
  Nr = r2-r1+1;
  Ns = s2-s1+1;

  crop = MRIallocSequence(Nc,Nr,Ns,mri->type,mri->nframes);
  MRIcopyHeader(mri, crop);

  // Compute location of 1st vox in cropped volume
  Vox2RAS = MRIxfmCRS2XYZ(mri,0);
  crs = MatrixAlloc(4,1,MATRIX_REAL);
  crs->rptr[1][1] = c1;
  crs->rptr[2][1] = r1;
  crs->rptr[3][1] = s1;
  crs->rptr[4][1] =  1;
  P0 = MatrixMultiply(Vox2RAS,crs,NULL);

  // Update header geometry for cropped
  MRIp0ToCRAS(crop, P0->rptr[1][1], P0->rptr[2][1], P0->rptr[3][1]);

  // Fill the value
  for (c=0; c < crop->width; c++)
  {
    for (r=0; r < crop->height; r++)
    {
      for (s=0; s < crop->depth; s++)
      {
        for (f=0; f < crop->nframes; f++)
        {
          v = MRIgetVoxVal(mri,c+c1,r+r1,s+s1,f);
          MRIsetVoxVal(crop,c,r,s,f,v);
        }
      }
    }
  }

  MatrixFree(&Vox2RAS);
  MatrixFree(&crs);
  MatrixFree(&P0);

  return(crop);
}

/*-----------------------------------------------------------*/
/*!
  \fn MRI *MRIuncrop(MRI *mri,MRI *crop,
  int c1,int r1,int s1, int c2, int r2, int s2)
  \brief Uncrops volume, keeps correct goemetry.
  \param mri - template for full, uncropped volume
  \param crs1 - col, row, slice in template at
  which cropping started (inclusive)
  \param crs2 - col, row, slice in template at
  which cropping ended   (inclusive)
 */
MRI *MRIuncrop(MRI *mri, MRI *crop,
               int c1, int r1, int s1,
               int c2, int r2, int s2)
{
  int c, r, s, f;
  MRI *uncrop;
  double v;

  if (c1 < 0 || c1 >= mri->width || r1 < 0 || r1 >= mri->height ||
      s1 < 0 || s1 >= mri->depth)
  {
    printf("MRIuncrop(): start point %d %d %d out of range\n",c1,r1,s1);
    return(NULL);
  }

  if (c2 < 0 || c2 >= mri->width || r2 < 0 || r2 >= mri->height ||
      s2 < 0 || r2 >= mri->depth)
  {
    printf("MRIuncrop(): end point %d %d %d out of range\n",c2,r2,s2);
    return(NULL);
  }

  uncrop = MRIcloneBySpace(mri, crop->type, crop->nframes);

  // Fill the values
  for (c=0; c < crop->width; c++)
  {
    for (r=0; r < crop->height; r++)
    {
      for (s=0; s < crop->depth; s++)
      {
        for (f=0; f < crop->nframes; f++)
        {
          v = MRIgetVoxVal(crop,c,r,s,f);
          MRIsetVoxVal(uncrop,c+c1,r+r1,s+s1,f,v);
        }
      }
    }
  }

  return(uncrop);
}
/* ----------------------------------------------------------*/
/*!
  \fn MRI *MRIreverseSlices(MRI *in, MRI *out)
  \brief Reverses the order of the slices and updates the
    vox2ras matrix.
*/
MRI *MRIreverseSlices(MRI *in, MRI *out)
{
  int c,r,s,f;
  MATRIX *M, *invM, *Sin, *Sout;
  double v;

  if (in == out)
  {
    printf("ERROR: MRIreverseSlices(): cannot be done in-place\n");
    return(NULL);
  }

  out = MRIcopy(in,out);
  if (out == NULL) return(NULL);

  // vox2ras for the input
  Sin = MRIxfmCRS2XYZ(in,0);

  // M converts inCRS to outCRS
  M = MatrixAlloc(4,4,MATRIX_REAL);
  M->rptr[1][1] = 1.0;
  M->rptr[2][2] = 1.0;
  // for reversal: sliceout = (Nslices-1) - slicein
  M->rptr[3][3] = -1.0;
  M->rptr[3][4] = in->depth-1.0;
  M->rptr[4][4] = 1.0;
  invM = MatrixInverse(M,NULL);
  if (invM == NULL)
  {
    printf("ERROR: inverting M\n");
    return(NULL);
  }

  // vox2ras for the output
  Sout = MatrixMultiply(Sin,invM,NULL);
  MRIsetVoxelToRasXform(out, Sout);
  MatrixFree(&Sin);
  MatrixFree(&M);
  MatrixFree(&invM);
  MatrixFree(&Sout);

  for (s=0; s < out->depth; s++)
  {
    for (c=0; c < out->width; c++)
    {
      for (r=0; r < out->height; r++)
      {
        for (f=0; f < out->nframes; f++)
        {
          v = MRIgetVoxVal(in,c,r,s,f);
          MRIsetVoxVal(out,c,r,(in->depth-1.0)-s,f,v);
        }
      }
    }
  }

  return(out);
}
/*----------------------------------------------------------------*/
MRI *MRIcutEndSlices(MRI *mri, int ncut)
{
  MRI *out;
  int nslices;
  int c,r,s,f,scut;
  double v;

  nslices = mri->depth - 2*ncut;
  if (nslices <= 0)
  {
    printf("ERROR: MRIcutEndSlices(): ncut = %d, input only has %d \n",
           ncut,mri->depth);
    return(NULL);
  }

  out = MRIallocSequence(mri->width,
                         mri->height,
                         nslices,
                         mri->type,
                         mri->nframes);
  MRIcopyHeader(mri,out);

  for (c = 0; c < mri->width; c ++)
  {
    for (r = 0; r < mri->height; r ++)
    {
      scut = 0;
      for (s = ncut; s < mri->depth - ncut; s++)
      {
        for (f = 0; f < mri->nframes; f++)
        {
          v = MRIgetVoxVal(mri,c,r,s,f);
          MRIsetVoxVal(out,c,r,scut,f,v);
        }
        scut ++;
      }
    }
  }
  return(out);
}

/* ----------------------------------------------------------*/
/*!
  \fn int *MRIsegIdList(MRI *seg, int *nlist, int frame)
  \brief Returns a list of the unique segmentation ids in
   the volume. The number in the list is *nlist. The volume need not
   be an int or char, but it is probably what it will be.
*/
int *MRIsegIdList(MRI *seg, int *nlist, int frame)
{
  int nvoxels,r,c,s,nth;
  int *tmplist = NULL;
  int *segidlist = NULL;

  nvoxels = seg->width * seg->height * seg->depth;
  tmplist = (int *) calloc(sizeof(int),nvoxels);

  // First, load all voxels into a list
  nth = 0;
  for (c=0; c < seg->width; c++)
  {
    for (r=0; r < seg->height; r++)
    {
      for (s=0; s < seg->depth; s++)
      {
        tmplist[nth] = (int) MRIgetVoxVal(seg,c,r,s,frame);
        nth++;
      }
    }
  }

  segidlist = unqiue_int_list(tmplist, nvoxels, nlist);
  free(tmplist);
  //for(nth=0; nth < *nlist; nth++)
  //printf("%3d %3d\n",nth,segidlist[nth]);
  return(segidlist);
}

/* ----------------------------------------------------------*/
/*!
  \fn double *MRIsegDice(MRI *seg1, MRI *seg2, int *nsegs, int **segidlist)
  \brief Computes dice coefficient for each segmentation. seg1 and seg2
  should have the same number of segs. Note: to get the name of
  the seg, CTABcopyName(ctab,segidlist[k],tmpstr,sizeof(tmpstr));
*/
double *MRIsegDice(MRI *seg1, MRI *seg2, int *nsegs, int **segidlist)
{
  int k,c,r,s,id1,id2,k1=0,k2=0;
  int nsegid1, *segidlist1;
  int nsegid2, *segidlist2;
  int *n1, *n2, *n12;
  double *dice;
  *nsegs = -1;

  // Extract a unique, sorted list of the ids
  segidlist1 = MRIsegIdList(seg1, &nsegid1,0);
  segidlist2 = MRIsegIdList(seg1, &nsegid2,0);

  if (nsegid1 != nsegid2)
  {
    printf("ERROR: MRIsegDice(): nsegs do not match %d %d\n",
           nsegid1,nsegid2);
    return(NULL);
  }
  printf("MRIsegDice(): found %d segs\n",nsegid1);
  *nsegs = nsegid1;

  n1  = (int *) calloc(nsegid1,sizeof(int));
  n2  = (int *) calloc(nsegid1,sizeof(int));
  n12 = (int *) calloc(nsegid1,sizeof(int));

  for (c=0; c < seg1->width; c++)
  {
    for (r=0; r < seg1->height; r++)
    {
      for (s=0; s < seg1->depth; s++)
      {
        // segid for 1st seg vol
        id1 = MRIgetVoxVal(seg1,c,r,s,0);
        for (k=0; k < nsegid1; k++)
        {
          if (id1 == segidlist1[k])
          {
            k1 = k;
            break;
          }
        }
        // segid for 2nd seg vol
        id2 = MRIgetVoxVal(seg2,c,r,s,0);
        for (k=0; k < nsegid2; k++)
        {
          if (id2 == segidlist2[k])
          {
            k2 = k;
            break;
          }
        }
        n1[k1]++;
        n2[k2]++;
        if (id1 == id2) n12[k1]++;
      }
    }
  }

  dice  = (double *) calloc(nsegid1,sizeof(double));
  for (k=0; k < nsegid1; k++)
    dice[k] = (double)n12[k]/((n1[k]+n2[k])/2.0);

  free(n1);
  free(n2);
  free(n12);
  free(segidlist2);
  *segidlist = segidlist1;

  return(dice);
}

/* ----------------------------------------------------------*/
/*!
  \fn MRI *MRIsegDiff(MRI *old, MRI *new, int *DiffFlag)
  \brief Determines differences between old and new segmentation.
  Voxels that are different will take the value of the new
  segmentation. Voxels that are NOT different will take the
  value of VOXEL_UNCHANGED (defined as 256 in cma.h and LUT.txt).
  This allows the diff to be loaded as a segmentation in tkmedit.
  If a difference was detected DiffFlag is set to 1, otherwise 0.
  Note that if there is no difference, all voxels will have the
  value of VOXEL_UNCHANGED. See also MRIsegMergeDiff().
*/
MRI *MRIsegDiff(MRI *old, MRI *new, int *DiffFlag)
{
  MRI *diff;
  int c,r,s;
  int vold, vnew, vdiff;

  diff = MRIallocSequence(new->width, new->height, new->depth, MRI_INT, 1);
  MRIcopyHeader(new,diff);

  *DiffFlag = 0;
  for (c=0; c < new->width; c++)
  {
    for (r=0; r < new->height; r++)
    {
      for (s=0; s < new->depth; s++)
      {
        vold = MRIgetVoxVal(old,c,r,s,0);
        vnew = MRIgetVoxVal(new,c,r,s,0);
        if (vold == vnew) vdiff = VOXEL_UNCHANGED;
        else
        {
          vdiff = vnew;
          *DiffFlag = 1;
        }
        MRIsetVoxVal(diff,c,r,s,0,vdiff);
      }
    }
  }

  return(diff);
}

/* ----------------------------------------------------------*/
/*!
  \fn MRI *MRIsegMergeDiff(MRI *old, MRI *diff)
  \brief Merges a segmentation with a "diff" segmentation to create a
  new segmentation. Voxels that have a diff value of VOXEL_UNCHANGED
  (cma.h) take their value from the "old" segmentation. Voxels that
  have a diff value other than VOXEL_UNCHANGED take their value from
  the "diff" segmentation.  See also MRIsegDiff().
*/
MRI *MRIsegMergeDiff(MRI *old, MRI *diff)
{
  MRI *new;
  int c,r,s;
  int vold, vdiff, vnew;

  new = MRIallocSequence(old->width, old->height, old->depth, MRI_INT, 1);
  MRIcopyHeader(old,new);

  for (c=0; c < new->width; c++)
  {
    for (r=0; r < new->height; r++)
    {
      for (s=0; s < new->depth; s++)
      {
        vold  = MRIgetVoxVal(old,c,r,s,0);
        vdiff = MRIgetVoxVal(diff,c,r,s,0);
        if (vdiff == VOXEL_UNCHANGED) vnew = vold;
        else                         vnew = MRIgetVoxVal(diff,c,r,s,0);
        MRIsetVoxVal(new,c,r,s,0,vnew);
      }
    }
  }

  return(new);
}

/* ----------------------------------------------------------*/
/*!
  \fn MRI *MRIhalfCosBias(MRI *in, double alpha, MRI *out)
  \brief Applies intensity bias using a half cosine. Each
    slice is rescaled separately. The scale factor is given
    by f = (alpha*cos(pi*s/(Ns-1))+(2-alpha))/2. Alpha is
    a control parameter. When 0, there is no bias. 1 is full
    bias.
*/
MRI *MRIhalfCosBias(MRI *in, double alpha, MRI *out)
{
  int c,r,s,f;
  double v,w;

  out = MRIcopy(in,out);

  for (s=0; s < out->depth; s++)
  {
    w = (alpha*cos(M_PI*s/(out->depth-1.0))+(2.0-alpha))/2.0;
    //printf("%2d %g\n",s,w);
    for (c=0; c < out->width; c++)
    {
      for (r=0; r < out->height; r++)
      {
        for (f=0; f < out->nframes; f++)
        {
          v = MRIgetVoxVal(in,c,r,s,f);
          MRIsetVoxVal(out,c,r,s,f,w*v);
        }
      }
    }
  }

  return(out);
}


int MRIvol2VolVSM(MRI *src, MRI *targ, MATRIX *Vt2s,
                  int InterpCode, float param, MRI *vsm)
{
  int   ct,  rt,  st,  f;
  int   ics, irs, iss, cvsm,rvsm;
  float fcs, frs, fss;
  float *valvect,drvsm;
  int sinchw;
  double rval,v;
  MATRIX *V2Rsrc=NULL, *invV2Rsrc=NULL, *V2Rtarg=NULL;
  MATRIX *crsT=NULL,*crsS=NULL;
  int FreeMats=0;

  if(DIAG_VERBOSE_ON)
    printf("Using MRIvol2VolVSM\n");

  if(src->nframes != targ->nframes){
    printf("ERROR: MRIvol2vol: source and target have different number "
           "of frames\n");
    return(1);
  }

  // Compute vox2vox matrix based on vox2ras of src and target.
  // Assumes that src and targ have same RAS space.
  if (Vt2s == NULL)
  {
    V2Rsrc = MRIxfmCRS2XYZ(src,0);
    invV2Rsrc = MatrixInverse(V2Rsrc,NULL);
    V2Rtarg = MRIxfmCRS2XYZ(targ,0);
    Vt2s = MatrixMultiply(invV2Rsrc,V2Rtarg,NULL);
    FreeMats = 1;
  }
  if (Gdiag_no > 0)
  {
    printf("MRIvol2Vol: Vt2s Matrix (%d)\n",FreeMats);
    MatrixPrint(stdout,Vt2s);
  }

  sinchw = nint(param);
  valvect = (float *) calloc(sizeof(float),src->nframes);

  MRI_BSPLINE * bspline = NULL;
  if (InterpCode == SAMPLE_CUBIC_BSPLINE)
    bspline = MRItoBSpline(src,NULL,3);
    
  crsT = MatrixAlloc(4,1,MATRIX_REAL);
  crsT->rptr[4][1] = 1;
  crsS = MatrixAlloc(4,1,MATRIX_REAL);
  for (ct=0; ct < targ->width; ct++)  {
    for (rt=0; rt < targ->height; rt++) {
      for (st=0; st < targ->depth; st++) {

        // Compute CRS in VSM space
        crsT->rptr[1][1] = ct;
        crsT->rptr[2][1] = rt;
        crsT->rptr[3][1] = st;
        crsS = MatrixMultiply(Vt2s,crsT,crsS);

        fcs = crsS->rptr[1][1];
        frs = crsS->rptr[2][1];
        fss = crsS->rptr[3][1];
        ics = nint(fcs);
        irs = nint(frs);
        iss = nint(fss);

        if(ics < 0 || ics >= src->width)  continue;
        if(irs < 0 || irs >= src->height) continue;
        if(iss < 0 || iss >= src->depth)  continue;

        if(vsm){
          /* Compute the voxel shift (converts from vsm
             space to mov space). This does a 3d interp to
             get vsm, not sure if really want a 2d*/
	  // Dont sample outside the BO mask
	  cvsm = floor(fcs);
	  rvsm = floor(frs);
	  if(cvsm < 0 || cvsm+1 >= src->width)  continue;
	  if(rvsm < 0 || rvsm+1 >= src->height) continue;
	  v = MRIgetVoxVal(vsm,cvsm,rvsm,iss,0);
	  if(fabs(v) < FLT_MIN) continue;
	  v = MRIgetVoxVal(vsm,cvsm+1,rvsm,iss,0);
	  if(fabs(v) < FLT_MIN) continue;
	  v = MRIgetVoxVal(vsm,cvsm,rvsm+1,iss,0);
	  if(fabs(v) < FLT_MIN) continue;
	  v = MRIgetVoxVal(vsm,cvsm+1,rvsm+1,iss,0);
	  if(fabs(v) < FLT_MIN) continue;
          MRIsampleSeqVolume(vsm, fcs, frs, fss, &drvsm, 0, 0);
	  if(drvsm == 0) continue;
          frs += drvsm;
          irs = nint(frs);
          if (irs < 0 || irs >= src->height) continue;
        }

        /* Assign output volume values */
        if (InterpCode == SAMPLE_TRILINEAR)
          MRIsampleSeqVolume(src, fcs, frs, fss, valvect,
                             0, src->nframes-1) ;
        else{
          for (f=0; f < src->nframes ; f++){
            switch (InterpCode){
            case SAMPLE_NEAREST:
              valvect[f] = MRIgetVoxVal(src,ics,irs,iss,f);
              break ;
            case SAMPLE_CUBIC_BSPLINE:     
              MRIsampleBSpline(bspline, fcs, frs, fss, f, &rval);
              valvect[f] = rval;
              break ;
            case SAMPLE_SINC:      /* no multi-frame */
              MRIsincSampleVolume(src, fcs, frs, fss, sinchw, &rval) ;
              valvect[f] = rval;
              break ;
            default:
              printf("ERROR: MRIvol2volVSM: interpolation method %i unknown\n",InterpCode);
              exit(1);
            }
          }
        }

        for (f=0; f < src->nframes; f++)
          MRIsetVoxVal(targ,ct,rt,st,f,valvect[f]);

      } /* target col */
    } /* target row */
  } /* target slice */

  free(valvect);
  MatrixFree(&crsS);
  MatrixFree(&crsT);
  if (bspline) MRIfreeBSpline(&bspline);
  if (FreeMats)
  {
    MatrixFree(&V2Rsrc);
    MatrixFree(&invV2Rsrc);
    MatrixFree(&V2Rtarg);
    MatrixFree(&Vt2s);
  }

  return(0);
}

/*---------------------------------------------------------------*/


MRI *MRIvol2surfVSM( const MRI *SrcVol,
                     const MATRIX *Rtk,
                     const MRI_SURFACE *TrgSurf,
                     const MRI *vsm, int InterpMethod, MRI *SrcHitVol,
                     float ProjFrac, int ProjType, int nskip, 
		     MRI *TrgVol )
{
  MATRIX *ras2vox, *vox2ras;
  AffineVector Scrs, Txyz;
  AffineMatrix ras2voxAffine;
  int   irow, icol, islc; /* integer row, col, slc in source */
  int cvsm,rvsm;
  float frow, fcol, fslc; /* float row, col, slc in source */
  float srcval, *valvect, rshift;
  int frm, vtx,nhits, err;
  double rval,val;
  float Tx, Ty, Tz;
  const VERTEX *v ;

#ifdef MRI2_TIMERS
  Chronometer tLoop;
  InitChronometer( &tLoop );
#endif

  if (vsm)  {
    err = MRIdimMismatch(vsm,SrcVol,0);
    if (err)    {
      printf("ERROR: MRIvol2surfVSM: vsm dimension mismatch %d\n",err);
      exit(1);
    }
  }

  vox2ras = MRIxfmCRS2XYZtkreg(SrcVol);
  ras2vox = MatrixInverse(vox2ras,NULL);
  if (Rtk != NULL) MatrixMultiply(ras2vox,Rtk,ras2vox);
  MatrixFree(&vox2ras);
  // ras2vox now converts surfacs RAS to SrcVol vox


  /* allocate a "volume" to hold the output */
  if(TrgVol == NULL){
    TrgVol = MRIallocSequence(TrgSurf->nvertices,1,1,MRI_FLOAT,SrcVol->nframes);
    if(TrgVol == NULL) return(NULL);
    MRIcopyHeader(SrcVol,TrgVol);
  } else {
    if(TrgVol->width != TrgSurf->nvertices || TrgVol->nframes != SrcVol->nframes){
      printf("ERROR: MRIvol2surfVSM: dimension mismatch (%d,%d), or (%d,%d)\n",
	     TrgVol->width,TrgSurf->nvertices,TrgVol->nframes,SrcVol->nframes);
      return(NULL);
    }
    // make sure all values are zero
    MRIconst(TrgVol->width,TrgVol->height,TrgVol->depth,1,0,TrgVol); 
  }
  // Dims here are meaningless, but setting to 1 means "volume" will be 
  // number of vertices.
  TrgVol->xsize = 1; 
  TrgVol->ysize = 1;
  TrgVol->zsize = 1;

  /* Zero the source hit volume */
  if (SrcHitVol != NULL)
    MRIconst(SrcHitVol->width,SrcHitVol->height,SrcHitVol->depth,
             1,0,SrcHitVol);

  srcval = 0;
  valvect = (float *) calloc(sizeof(float),SrcVol->nframes);
  nhits = 0;

  SetAffineMatrix( &ras2voxAffine, ras2vox );

  MRI_BSPLINE * bspline = NULL;
  if (InterpMethod == SAMPLE_CUBIC_BSPLINE)
    bspline = MRItoBSpline(SrcVol,NULL,3);


  /*--- loop through each vertex ---*/
#ifdef MRI2_TIMERS
  StartChronometer( &tLoop );
  unsigned int skipped = 0;
#endif
  for (vtx = 0; vtx < TrgSurf->nvertices; vtx+=nskip)
  {
    v = &TrgSurf->vertices[vtx] ;
    if( v->ripflag ) {
#ifdef MRI2_TIMERS
      skipped++;
#endif
      continue;
    }

    if (ProjFrac != 0.0)
    {
      if (ProjType == 0)
        ProjNormDist(&Tx,&Ty,&Tz,TrgSurf,vtx,ProjFrac);
      else
        ProjNormFracThick(&Tx,&Ty,&Tz,TrgSurf,vtx,ProjFrac);
    }
    else
    {
      Tx = v->x;
      Ty = v->y;
      Tz = v->z;
    }

    /* Load the Target xyz vector */
    SetAffineVector( &Txyz, Tx, Ty, Tz );
    /* Compute the corresponding Source col-row-slc vector */
    AffineMV( &Scrs, &ras2voxAffine, &Txyz );
    GetAffineVector( &Scrs, &fcol, &frow, &fslc );

    icol = nint(fcol);
    irow = nint(frow);
    islc = nint(fslc);

    /* check that the point is in the bounds of the volume */
    if (irow < 0 || irow >= SrcVol->height ||
        icol < 0 || icol >= SrcVol->width  ||
        islc < 0 || islc >= SrcVol->depth ) continue;

    if (vsm){
      /* Compute the voxel shift (converts from vsm
	 space to mov space). This does a 3d interp to
	 get vsm, not sure if really want a 2d*/
      // Dont sample outside the BO mask
      cvsm = floor(fcol);
      rvsm = floor(frow);
      if(cvsm < 0 || cvsm+1 >= vsm->width)  continue;
      if(rvsm < 0 || rvsm+1 >= vsm->height) continue;
      val = MRIgetVoxVal(vsm,cvsm,rvsm,islc,0);
      if(fabs(val) < FLT_MIN) continue;
      val = MRIgetVoxVal(vsm,cvsm+1,rvsm,islc,0);
      if(fabs(val) < FLT_MIN) continue;
      val = MRIgetVoxVal(vsm,cvsm,rvsm+1,islc,0);
      if(fabs(val) < FLT_MIN) continue;
      val = MRIgetVoxVal(vsm,cvsm+1,rvsm+1,islc,0);
      if(fabs(val) < FLT_MIN) continue;
      MRIsampleSeqVolume(vsm, fcol, frow, fslc, &rshift, 0, 0);
      if(rshift == 0) continue;
      frow += rshift;
      irow = nint(frow);
      if (irow < 0 || irow >= SrcVol->height) continue;
    }

#if 0
    if (Gdiag_no == vtx)
    {
      printf("diag -----------------------------\n");
      printf("vtx = %d  %g %g %g\n",vtx,Txyz->rptr[1][1],
             Txyz->rptr[2][1],Txyz->rptr[3][1]);
      printf("fCRS  %g %g %g\n",Scrs->rptr[1][1],
             Scrs->rptr[2][1],Scrs->rptr[3][1]);
      printf("CRS  %d %d %d\n",icol,irow,islc);
    }
#endif

    /* only gets here if it is in bounds */
    nhits ++;

    /* Assign output volume values */
    if (InterpMethod == SAMPLE_TRILINEAR)
    {
      MRIsampleSeqVolume(SrcVol, fcol, frow, fslc,
                         valvect, 0, SrcVol->nframes-1) ;
      if (Gdiag_no == vtx) printf("val = %f\n", valvect[0]) ;
      for (frm = 0; frm < SrcVol->nframes; frm++)
        MRIFseq_vox(TrgVol,vtx,0,0,frm) = valvect[frm];
    }
    else
    {
      for (frm = 0; frm < SrcVol->nframes; frm++)
      {
        switch (InterpMethod)
        {
        case SAMPLE_NEAREST:
          srcval = MRIgetVoxVal(SrcVol,icol,irow,islc,frm);
          break ;
        case SAMPLE_CUBIC_BSPLINE:
          MRIsampleBSpline(bspline, fcol, frow, fslc, frm, &rval);
          srcval = rval;
          break ;
        case SAMPLE_SINC:      /* no multi-frame */
          MRIsincSampleVolume(SrcVol, fcol, frow, fslc, 5, &rval) ;
          srcval = rval;
          break ;
        default:
          printf("ERROR: MRIvol2surfVSM: interpolation method %i unknown\n",InterpMethod);
          exit(1);
        } //switch
        MRIFseq_vox(TrgVol,vtx,0,0,frm) = srcval;
        if (Gdiag_no == vtx) printf("val[%d] = %f\n", frm, srcval) ;
      } // for
    }// else
    if (SrcHitVol != NULL) MRIFseq_vox(SrcHitVol,icol,irow,islc,0)++;
  }
#ifdef MRI2_TIMERS
  StopChronometer( &tLoop );
  printf( "%s: Main Loop complete in %6.3f ms (%6u %6u)\n",
          __FUNCTION__, GetChronometerValue( &tLoop ), skipped, nhits );
#endif

  MatrixFree(&ras2vox);
  free(valvect);
  if (bspline) MRIfreeBSpline(&bspline);

  //printf("vol2surf_linear: nhits = %d/%d\n",nhits,TrgSurf->nvertices);

  return(TrgVol);
}

int MRIvol2VolTkRegVSM(MRI *mov, MRI *targ, MATRIX *Rtkreg,
                       int InterpCode, float param, MRI *vsm)
{
  MATRIX *vox2vox = NULL;
  MATRIX *Tmov, *invTmov, *Ttarg;
  int err;

  if (Rtkreg != NULL)
  {
    // TkReg Vox2RAS matrices
    Tmov      = MRIxfmCRS2XYZtkreg(mov);
    invTmov   = MatrixInverse(Tmov,NULL);
    Ttarg    = MRIxfmCRS2XYZtkreg(targ);
    // vox2vox = invTmov*R*Ttarg
    vox2vox = MatrixMultiply(invTmov,Rtkreg,vox2vox);
    MatrixMultiply(vox2vox,Ttarg,vox2vox);
  }
  else vox2vox = NULL;

  // resample
  err = MRIvol2VolVSM(mov,targ,vox2vox,InterpCode,param,vsm);

  if (vox2vox)
  {
    MatrixFree(&vox2vox);
    MatrixFree(&Tmov);
    MatrixFree(&invTmov);
    MatrixFree(&Ttarg);
  }

  return(err);
}

/* ----------------------------------------------------------*/
/*!
  \fn MRI *MRIsegBoundary(MRI *seg)
  \brief Creates a segmentation boundary volume in which voxels
  whose neighbors are the same seg are set to 0.
*/
MRI *MRIsegBoundary(MRI *seg)
{
  MRI *boundary;
  int c,r,s,dc,dr,ds;
  int cseg, nseg, b;

  boundary = MRIclone(seg,NULL);

  for(c=1; c < seg->width-1; c++){
    for(r=1; r < seg->height-1; r++){
      for(s=1; s < seg->depth-1; s++){
	cseg = (int) MRIgetVoxVal(seg,c,r,s,0);
	if(cseg == 0) {
	  MRIsetVoxVal(boundary,c,r,s,0, 0);
	  continue;
	}
	b = 0;
	for(dc = -1; dc < 2; dc++){
	  for(dr = -1; dr < 2; dr++){
	    for(ds = -1; ds < 2; ds++){
	      nseg = (int) MRIgetVoxVal(seg,c+dc,r+dr,s+ds,0);
	      if(cseg != nseg){
		b=1; // It is a boundary voxel
		dc=2;dr=2;ds=2; // break from all three loops
	      }
	    }
	  }
	}
	if(b==0) MRIsetVoxVal(boundary,c,r,s,0, 0);
	else     MRIsetVoxVal(boundary,c,r,s,0, cseg);
      }
    }
  }

  return(boundary);
}

/*!
  \fn MRI *MRIcrs(MRI *in, MRI *out)
  \brief Creates a 1-frame volume with the value
  at frame 0 equal to the slice number.
*/
MRI *MRIsliceNo(MRI *in, MRI *out)
{
  int c,r,s;
  if(out == NULL){
    out = MRIalloc(in->width,in->height,in->depth,MRI_FLOAT);
    MRIcopyHeader(in,out);
  }

  for(c=0; c < in->width; c++){
    for(r=0; r < in->height; r++){
      for(s=0; s < in->depth; s++){
	MRIsetVoxVal(out,c,r,s,0, s);
      }
    }
  }

  return(out);
}

/*!
  \fn MRI *MRIcrs(MRI *in, MRI *out)
  \brief Creates a 1-frame volume with the value
  at frame 0 equal to the voxel index, eg, 
  sliceno*nrows*ncols + rowno*ncols + col.
*/
MRI *MRIindexNo(MRI *in, MRI *out)
{
  int c,r,s, index;
  if(out == NULL){
    out = MRIalloc(in->width,in->height,in->depth,MRI_FLOAT);
    MRIcopyHeader(in,out);
  }

  index = 0;
  for(s=0; s < in->depth; s++){
    for(r=0; r < in->height; r++){
      for(c=0; c < in->width; c++){
	MRIsetVoxVal(out,c,r,s,0, index);
	index ++;
      }
    }
  }

  return(out);
}

/* ----------------------------------------------------------*/
/*!
  \fn MRI *MRIcrs(MRI *in, MRI *out)
  \brief Creates a 3-frame volume with the value
  at frame 0 equal to the column number, frame 1
  the row, and frame 2 the slice.
*/
MRI *MRIcrs(MRI *in, MRI *out)
{
  int c,r,s;

  if(out == NULL){
    out = MRIallocSequence(in->width,in->height,in->depth,MRI_FLOAT,3);
    MRIcopyHeader(in,out);
  }

  for(s=0; s < in->depth; s++){
    for(r=0; r < in->height; r++){
      for(c=0; c < in->width; c++){
	MRIsetVoxVal(out,c,r,s,0, c);
	MRIsetVoxVal(out,c,r,s,1, r);
	MRIsetVoxVal(out,c,r,s,2, s);
      }
    }
  }

  return(out);
}

/*---------------------------------------------------------
  MRIsegStats() - computes statistics within a given
  segmentation. Returns the number of voxels in the
  segmentation.
  ---------------------------------------------------------*/
int MRIsegStats(MRI *seg, int segid, MRI *mri,int frame,
                float *min, float *max, float *range,
                float *mean, float *std)
{
  int id,nvoxels,r,c,s;
  double val, sum, sum2;

  *min = 0;
  *max = 0;
  sum  = 0;
  sum2 = 0;
  nvoxels = 0;
  for (c=0; c < seg->width; c++)
  {
    for (r=0; r < seg->height; r++)
    {
      for (s=0; s < seg->depth; s++)
      {
        id = (int) MRIgetVoxVal(seg,c,r,s,0);
        if (id != segid)
        {
          continue;
        }
        val =  MRIgetVoxVal(mri,c,r,s,frame);
        nvoxels++;
        if ( nvoxels == 1 )
        {
          *min = val;
          *max = val;
        }
        if (*min > val)
        {
          *min = val;
        }
        if (*max < val)
        {
          *max = val;
        }
        sum  += val;
        sum2 += (val*val);
      }
    }
  }

  *range = *max - *min;

  if (nvoxels != 0)
  {
    *mean = sum/nvoxels;
  }
  else
  {
    *mean = 0.0;
  }

  if (nvoxels > 1)
    *std = sqrt(((nvoxels)*(*mean)*(*mean) - 2*(*mean)*sum + sum2)/
                (nvoxels-1));
  else
  {
    *std = 0.0;
  }

  return(nvoxels);
}
/*------------------------------------------------------------*/
/*!
  \fn int MRIsegStatsRobust(MRI *seg, int segid, MRI *mri,int frame,
		      float *min, float *max, float *range,
		      float *mean, float *std, float Pct)
  \brief Computes stats based on the the middle 100-2*Pct values, ie,
         it trims Pct off the ends.
*/
int MRIsegStatsRobust(MRI *seg, int segid, MRI *mri,int frame,
		      float *min, float *max, float *range,
		      float *mean, float *std, float Pct)
{
  int id,nvoxels,r,c,s,k,m;
  double val, sum, sum2;
  float *vlist;

  *min = 0;
  *max = 0;
  *range = 0;
  *mean = 0;
  *std = 0;

  // Count number of voxels
  nvoxels = 0;
  for (c=0; c < seg->width; c++) {
    for (r=0; r < seg->height; r++)  {
      for (s=0; s < seg->depth; s++)  {
        id = (int) MRIgetVoxVal(seg,c,r,s,0);
        if (id != segid) continue;
        nvoxels++;
      }
    }
  }
  if(nvoxels == 0) return(nvoxels);

  // Load voxels into an array
  vlist = (float *) calloc(sizeof(float),nvoxels);
  nvoxels = 0;
  for (c=0; c < seg->width; c++) {
    for (r=0; r < seg->height; r++)  {
      for (s=0; s < seg->depth; s++)  {
        id = (int) MRIgetVoxVal(seg,c,r,s,0);
        if (id != segid) continue;
	vlist[nvoxels] = MRIgetVoxVal(mri,c,r,s,frame);
        nvoxels++;
      }
    }
  }
  // Sort the array
  qsort((void *) vlist, nvoxels, sizeof(float), compare_floats);

  // Compute stats excluding Pct of the values from each end
  sum  = 0;
  sum2 = 0;
  m = 0;
  //printf("Robust Indices: %d %d\n",(int)nint(Pct*nvoxels/100.0),(int)nint((100-Pct)*nvoxels/100.0));
  for(k=0; k < nvoxels; k++){
    if(k < Pct*nvoxels/100.0)       continue;
    if(k > (100-Pct)*nvoxels/100.0) continue;
    val = vlist[k];
    if(m == 0){
      *min = val;
      *max = val;
    }
    if (*min > val) *min = val;
    if (*max < val) *max = val;
    sum  += val;
    sum2 += (val*val);
    m = m + 1;
  }

  *range = *max - *min;
  *mean = sum/m;
  if(m > 1)
    *std = sqrt(((m)*(*mean)*(*mean) - 2*(*mean)*sum + sum2)/
                (m-1));
  else *std = 0.0;

  free(vlist);
  vlist = NULL;
  return(m);
}
/*---------------------------------------------------------
  MRIsegFrameAvg() - computes the average time course withing the
  given segmentation. Returns the number of voxels in the
  segmentation. favg must be preallocated to number of
  frames. favg = (double *) calloc(sizeof(double),mri->nframes);
  ---------------------------------------------------------*/
int MRIsegFrameAvg(MRI *seg, int segid, MRI *mri, double *favg)
{
  int id,nvoxels,r,c,s,f;
  double val;

  /* zero it out */
  for (f=0; f<mri->nframes; f++)
  {
    favg[f] = 0;
  }

  nvoxels = 0;
  for (c=0; c < seg->width; c++)
  {
    for (r=0; r < seg->height; r++)
    {
      for (s=0; s < seg->depth; s++)
      {
        id = (int) MRIgetVoxVal(seg,c,r,s,0);
        if (id != segid)
        {
          continue;
        }
        for (f=0; f<mri->nframes; f++)
        {
          val =  MRIgetVoxVal(mri,c,r,s,f);
          favg[f] += val;
        }
        nvoxels++;
      }
    }
  }

  if (nvoxels != 0)
    for (f=0; f<mri->nframes; f++)
    {
      favg[f] /= nvoxels;
    }

  return(nvoxels);
}

MRI *
MRImask_with_T2_and_aparc_aseg(MRI *mri_src, MRI *mri_dst, MRI *mri_T2, MRI *mri_aparc_aseg, float T2_thresh, int mm_from_exterior)
{
  int    x, y, z, nremoved, i ;
  MRI    *mri_bright, *mri_mask , *mri_tmp = NULL;

  mri_mask = MRIbinarize(mri_T2, NULL, T2_thresh, 255, 0) ;
  mri_bright = MRIcopy(mri_mask, NULL) ;


  if (mri_aparc_aseg)   // use T2 and aparc+aseg to remove non-brain stuff
  {
    MRIbinarize(mri_aparc_aseg, mri_aparc_aseg, 1, 0, 255) ;
    MRIdilate(mri_aparc_aseg, mri_aparc_aseg) ;
    MRInot(mri_aparc_aseg, mri_aparc_aseg) ;  // background now on, foreground off
    GetLargestCC6(mri_aparc_aseg) ;                 // remove disconnected background components
    MRIand(mri_mask, mri_aparc_aseg, mri_mask, 1) ;
    MRIopenN(mri_mask, mri_mask, 3) ;   // third order open will remove thin chains of bright T2 that are in the interior
  }
  else  // just use T2
  {
    GetLargestCC6(mri_mask) ;
  }


  MRInot(mri_mask, mri_mask) ;   // 0 now means mask it out and 1 means retain it

  for (i = nremoved = 0 ; i < mm_from_exterior ; i++)
  {
    mri_tmp =  MRIcopy(mri_mask, mri_tmp) ;
    for (x = 0 ; x < mri_mask->width ;  x++)
      for (y = 0 ; y < mri_mask->height ;  y++)
	for (z = 0 ; z < mri_mask->depth ;  z++)
	{
	  if (x == Gx && y == Gy && z == Gz)
	    DiagBreak() ;
	  if (MRIgetVoxVal(mri_mask, x, y, z, 0) == 0)  // already in the mask
	    continue ;
	  if (MRIgetVoxVal(mri_aparc_aseg, x, y, z, 0) == 0)  // too close to brain
	    continue ;
	  if (MRIgetVoxVal(mri_T2, x, y, z, 0) >= T2_thresh)  // bright in the T2
	  {
	    if (MRIneighborsOff(mri_mask, x, y, z, 1) > 0)   // touching the existing mask
	    {
	      if (x == Gx && y == Gy && z == Gz)
		DiagBreak() ;
	      MRIsetVoxVal(mri_tmp, x, y, z, 0, 0) ;   // add this voxel to the mask
	      nremoved++ ;
	    }
	  }
	}

    printf("%d T2-bright exterior voxels removed\n", nremoved) ;
    MRIcopy(mri_tmp, mri_mask) ;
  }


  mri_dst = MRImask(mri_dst, mri_mask, mri_dst, 0, 0) ;  // if mask == 0, then set dst as 0
  MRIfree(&mri_bright) ;  MRIfree(&mri_mask) ; MRIfree(&mri_tmp) ;
  return(mri_dst) ;
}

/* ----------------------------------------------------------*/
/*!
  \fn int *MRIsegmentationList(MRI *seg, int *pListLength)
  \brief Extracts a list of unique segmentation IDs from the volume. 
    Includes 0.
*/
int *MRIsegmentationList(MRI *seg, int *pListLength)
{
  int c,r,s,n,nvox;
  int *list, *voxlist;

  nvox = seg->width*seg->height*seg->depth;
  voxlist = (int *) calloc(nvox,sizeof(int));
  n = 0;
  for(s=0; s < seg->depth; s++){
    for(c=0; c < seg->width; c++){
      for(r=0; r < seg->height; r++){
	voxlist[n] = MRIgetVoxVal(seg,c,r,s,0);
	n++;
      }
    }
  }
  list = unqiue_int_list(voxlist, nvox, pListLength);
  printf("MRIsegmentationList(): found %d unique segmentations\n", *pListLength);
  if(Gdiag_no > 0) for(n=0; n<*pListLength; n++)  printf("%2d %5d\n",n,list[n]);
  free(voxlist);
  return(list);
}
