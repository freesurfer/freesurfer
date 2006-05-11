/*-------------------------------------------------------------------
  Name: mri2.c
  Author: Douglas N. Greve
  $Id: mri2.c,v 1.23 2006/05/11 20:30:16 kteich Exp $
  Purpose: more routines for loading, saving, and operating on MRI 
  structures.
  -------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "fio.h"
#include "stats.h"

#include "corio.h"
#include "bfileio.h"

#include "mri2.h"
#include "sig.h"

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
  if(bfvol == NULL) return(NULL);

  /* allocate the MRI */
  vol = MRIallocSequence(bfvol->ncols,bfvol->nrows,
       bfvol->nslcs,MRI_FLOAT,bfvol->nfrms);
  if(vol==NULL){
    bf_freebfd(&bfvol);
    fprintf(stderr,"mri_load_bvolume(): could not alloc vol\n");
    return(NULL);
  }

  /* copy data from the BF_DATA struct to the ARRAY4D*/
  for(r=0;r<bfvol->nrows;r++){
    for(c=0;c<bfvol->ncols;c++){
      for(s=0;s<bfvol->nslcs;s++){
  for(f=0;f<bfvol->nfrms;f++){
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
  if(bfvol == NULL) return(1);

  /* copy data from ARRAY4D to BF_DATA */
  for(r=0;r<bfvol->nrows;r++){
    for(c=0;c<bfvol->ncols;c++){
      for(s=0;s<bfvol->nslcs;s++){
  for(f=0;f<bfvol->nfrms;f++){
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
  if(bfvol == NULL) return(NULL);

  if(frameno >= bfvol->nfrms){
    fprintf(stderr,"ERROR: mri_load_bvolume_frame(): frameno = %d, exceeds "
      "number of frames in %s = %d\n",frameno,bfstem,bfvol->nfrms);
    bf_freebfd(&bfvol);
    return(NULL);
  }

  /* allocate the MRI */
  vol = MRIallocSequence(bfvol->ncols,bfvol->nrows,
       bfvol->nslcs,MRI_FLOAT,1);
  if(vol==NULL){
    bf_freebfd(&bfvol);
    fprintf(stderr,"mri_load_bvolume_frame(): could not alloc vol\n");
    return(NULL);
  }

  /* copy data from the BF_DATA struct to the ARRAY4D*/
  for(r=0;r<bfvol->nrows;r++){
    for(c=0;c<bfvol->ncols;c++){
      for(s=0;s<bfvol->nslcs;s++){
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
  
  if(frame >= vol->nframes){
    fprintf(stderr,"mri_save_as_cor(): frame = %d, must be <= %d\n",
      frame,vol->nframes);
    return(1);
  }

  /* make sure the output directory is writable */
  if(! cordir_iswritable(cordir)) return(1);

  /* allocate a temporary COR volume */
  COR = alloc_cor();
  if(COR==NULL) return(1);

  /* rescale to 0-255 (range of uchar) */
  if(rescale) mri_rescale(vol,0,255,vol);

  /* make sure maximum subscript does not 
     exceed either range */
  if(vol->height < 256) rmax = vol->height;
  else                  rmax = 256;
  if(vol->width < 256)  cmax = vol->width;
  else                  cmax = 256;
  if(vol->depth < 256)  smax = vol->depth;
  else                  smax = 256;

  /* copy data from ARRAY4D to COR */
  for(r=0; r < rmax ; r++){
    for(c=0; c < cmax ; c++){
      for(s=0; s < smax ; s++){
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

  if(outvol != NULL) tmpvol = outvol;
  else{
    tmpvol = MRIallocSequence(vol->width,vol->height,vol->depth,
            MRI_FLOAT,vol->nframes);
    if(tmpvol == NULL) return(NULL);
  }

  /* find the minimum and maximum */
  volmin = MRIFseq_vox(vol,0,0,0,0);
  volmax = MRIFseq_vox(vol,0,0,0,0);
  for(r=0;r<vol->height;r++){
    for(c=0;c<vol->width;c++){
      for(s=0;s<vol->depth;s++){
  for(f=0;f<vol->nframes;f++){
    val = MRIFseq_vox(vol,c,r,s,f);
    if(volmin > val) volmin = val;
    if(volmax < val) volmax = val;
  }
      }
    }
  }
  volrange = volmax - volmin;
  range    = max - min;
  
  /* rescale to the new range */
  for(r=0;r<vol->height;r++){
    for(c=0;c<vol->width;c++){
      for(s=0;s<vol->depth;s++){
  for(f=0;f<vol->nframes;f++){
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
  for(r=0;r<vol->height;r++){
    for(c=0;c<vol->width;c++){
      for(s=0;s<vol->depth;s++){
  for(f=0;f<vol->nframes;f++){
    val = MRIFseq_vox(vol,c,r,s,f);
    if(*min > val) *min = val;
    if(*max < val) *max = val;
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

  for(f=0;f<vol->nframes;f++){

    /* power = 1 -- don't do anything */
    if( fabs(framepower[f]-1.0) < .00001  )
      continue;

    /* power = 0.5 -- use sqrt() */
    if( fabs(framepower[f]-0.5) < .00001  ){
      for(r=0;r<vol->height;r++){
  for(c=0;c<vol->width;c++){
    for(s=0;s<vol->depth;s++){
      val = MRIFseq_vox(vol,c,r,s,f);
      MRIFseq_vox(vol,c,r,s,f) = sqrt(val);
    }
  }
      }
      continue;
    }

    /* power = 2 -- use val*val */
    if( fabs(framepower[f]-0.5) < .00001  ){
      for(r=0;r<vol->height;r++){
  for(c=0;c<vol->width;c++){
    for(s=0;s<vol->depth;s++){
      val = MRIFseq_vox(vol,c,r,s,f);
      MRIFseq_vox(vol,c,r,s,f) = val*val;
    }
  }
      }
      continue;
    }

    /* generic: use pow() -- least efficient */
    for(r=0;r<vol->height;r++){
      for(c=0;c<vol->width;c++){
  for(s=0;s<vol->depth;s++){
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

  if(tail == NULL) tail = "positive";

  /* check the first 3 letters of tail */
  if(!strncasecmp(tail,"positive",3)) tailcode = 1;
  else if(!strncasecmp(tail,"negative",3)) tailcode = 2;
  else if(!strncasecmp(tail,"absolute",3)) tailcode = 3;
  else{
    fprintf(stderr,"mri_binarize: tail = %s unrecoginzed\n",tail);
    return(NULL);
  }

  if(volbin == NULL){
    voltmp = MRIallocSequence(vol->width,vol->height,vol->depth,
            MRI_FLOAT,vol->nframes); 
    if(voltmp == NULL) return(NULL);
  }
  else voltmp = volbin;

  if(!invert){
    onval  = 1;
    offval = 0;
    printf("NOT INVERTING\n");
  }
  else{
    onval  = 0;
    offval = 1;
    printf("INVERTING\n");
  }

  *nover = 0;
  for(r=0;r<vol->height;r++){
    for(c=0;c<vol->width;c++){
      for(s=0;s<vol->depth;s++){
  for(f=0;f<vol->nframes;f++){
    val = MRIFseq_vox(vol,c,r,s,f); 
    switch(tailcode){
    case 2: val = -val; break;
    case 3: val = fabs(val); break;
    }
    if(val > thresh) b = onval;
    else             b = offval;
    if(b) (*nover) ++;
    MRIFseq_vox(voltmp,c,r,s,f) = b;
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
  if(ucvol == NULL) return(NULL);

  /* allocate a float volume */
  vol = MRIallocSequence(256,256,256,MRI_FLOAT,1);
  if(vol == NULL) return(NULL);

  for(r=0;r<vol->height;r++){
    for(c=0;c<vol->width;c++){
      for(s=0;s<vol->depth;s++){
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
  if (fp==NULL) {
    fprintf(stderr,"ERROR: Progname: mri_load_wfile():\n");
    fprintf(stderr,"Could not open %s\n",wfile);
    fprintf(stderr,"(%s,%d,%s)\n",__FILE__, __LINE__,__DATE__);
    return(NULL);
  }
  
  fread2(&ilat,fp);
  fread3(&num,fp);

  vtxnum = (int *)   calloc(sizeof(int),   num);
  wval   = (float *) calloc(sizeof(float), num);

  for (i=0;i<num;i++){
    fread3(&vtxnum[i],fp);
    wval[i] = freadFloat(fp) ;
  }
  fclose(fp);

  nvertices = vtxnum[num-1] + 1;

  w = MRIallocSequence(nvertices,1,1,MRI_FLOAT,1);
  for (i=0;i<num;i++){
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

  switch (vol->type){
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

  if(vol->nframes == 0) vol->nframes = 1;

  nv1 = vol->width * vol->height * vol->depth * vol->nframes;
  nv2 = ncols*nrows*nslices*nframes;

  if(nv1 != nv2){
    printf("ERROR: mri_reshape: number of elements cannot change\n");
    printf("  nv1 = %d, nv1 = %d\n",nv1,nv2);
    return(NULL);
  }

  outvol = MRIallocSequence(ncols,nrows,nslices,vol->type,nframes);
  if(outvol == NULL) return(NULL);

  MRIcopyHeader(vol, outvol); /* does not change dimensions */

  //printf("vol1: %d %d %d %d %d\n",vol->width,vol->height,
  // vol->depth,vol->nframes,mri_sizeof(vol));
  //printf("vol2: %d %d %d %d %d\n",outvol->width,outvol->height,
  // outvol->depth,outvol->nframes,mri_sizeof(outvol));

  c2 = 0; r2 = 0; s2 = 0; f2 = 0;
  for(f = 0; f < vol->nframes; f++){
    for(s = 0; s < vol->depth;   s++){
      for(r = 0; r < vol->height;  r++){
  for(c = 0; c < vol->width;   c++){

    switch (vol->type){
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
    if(c2 == ncols){
      c2 = 0; r2++;
      if(r2 == nrows){
        r2 = 0;  s2++;
        if(s2 == nslices){
    s2 = 0; f2++;
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
  int   ct,  rt,  st,  f;
  int   ics, irs, iss;
  float fcs, frs, fss;
  float *valvect;
  int sinchw;
  Real rval;
  MATRIX *V2Rsrc=NULL, *invV2Rsrc=NULL, *V2Rtarg=NULL;
  int FreeMats=0;

  if(src->nframes != targ->nframes){
    printf("ERROR: MRIvol2vol: source and target have different number "
	   "of frames\n");
    return(1);
  }

  // Compute vox2vox matrix based on vox2ras of src and target. 
  // Assumes that src and targ have same RAS space.
  if(Vt2s == NULL){
    V2Rsrc = MRIxfmCRS2XYZ(src,0);
    invV2Rsrc = MatrixInverse(V2Rsrc,NULL);
    V2Rtarg = MRIxfmCRS2XYZ(targ,0);
    Vt2s = MatrixMultiply(invV2Rsrc,V2Rtarg,NULL);
    FreeMats = 1;
  }
  if(Gdiag_no > 0){
    printf("MRIvol2Vol: Vt2s Matrix (%d)\n",FreeMats);
    MatrixPrint(stdout,Vt2s);
  }

  sinchw = nint(param);
  valvect = (float *) calloc(sizeof(float),src->nframes);

  for(ct=0; ct < targ->width; ct++){
    for(rt=0; rt < targ->height; rt++){
      for(st=0; st < targ->depth; st++){
	
	/* Column in source corresponding to CRS in Target */
	fcs = Vt2s->rptr[1][1] * ct + Vt2s->rptr[1][2] * rt +
	      Vt2s->rptr[1][3] * st + Vt2s->rptr[1][4] ;
	ics = nint(fcs);
	if(ics < 0 || ics >= src->width) continue;

	/* Row in source corresponding to CRS in Target */
	frs = Vt2s->rptr[2][1] * ct + Vt2s->rptr[2][2] * rt +
	      Vt2s->rptr[2][3] * st + Vt2s->rptr[2][4] ;
	irs = nint(frs);
	if(irs < 0 || irs >= src->height) continue;

	/* Slice in source corresponding to CRS in Target */
	fss = Vt2s->rptr[3][1] * ct + Vt2s->rptr[3][2] * rt +
	      Vt2s->rptr[3][3] * st + Vt2s->rptr[3][4] ;
	iss = nint(fss);
	if(iss < 0 || iss >= src->depth) continue;

	/* Assign output volume values */
	if(InterpCode == SAMPLE_TRILINEAR)
	  MRIsampleSeqVolume(src, fcs, frs, fss, valvect, 
			     0, src->nframes-1) ;
	else{
	  for(f=0; f < src->nframes ; f++){
	    switch(InterpCode){
	    case SAMPLE_NEAREST:
	      valvect[f] = MRIgetVoxVal(src,ics,irs,iss,f);
	      break ;
	    case SAMPLE_SINC:      /* no multi-frame */
	      MRIsincSampleVolume(src, fcs, frs, fss, sinchw, &rval) ;
	      valvect[f] = rval;
	      break ;
	    }
	  }
	}

	for(f=0; f < src->nframes; f++)
	  MRIsetVoxVal(targ,ct,rt,st,f,valvect[f]);

      } /* target col */
    } /* target row */
  } /* target slice */

  free(valvect);
  if(FreeMats){
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
int MRIdimMismatch(MRI *v1, MRI *v2, int frameflag)
{
  if(v1->width  != v2->width)  return(1);
  if(v1->height != v2->height) return(2);
  if(v1->depth  != v2->depth)  return(3);
  if(frameflag && v1->nframes != v2->nframes)  return(4);
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

  if(vol->nframes <= frame){
    printf("ERROR: MRIfdr2vwth: frame = %d, must be < nframes = %d\n",
	   frame,vol->nframes);
    return(1);
  }
  if(vol->type != MRI_FLOAT){
    printf("ERROR: MRIfdr2vwth: input volume is not of type MRI_FLOAT\n");
    return(1);
  }
  if(ovol != NULL){
    if(ovol->type != MRI_FLOAT){
      printf("ERROR: MRIfdr2vwth: output volume is not of type MRI_FLOAT\n");
      return(1);
    }
    if(MRIdimMismatch(vol, ovol, 0)){
      printf("ERROR: MRIfdr2vwth: output/input dimension mismatch\n");
      return(1);
    }
  }
  if(mask != NULL){
    if(MRIdimMismatch(vol, mask, 0)){
      printf("ERROR: MRIfdr2vwth: mask/input dimension mismatch\n");
      return(1);
    }
    maskflag = 1;
  }

  Nv = vol->width * vol->height * vol->depth;
  if(log10flag) valnull = 0;
  else          valnull = 1;

  p = (double *) calloc(Nv,sizeof(double));
  np = 0;
  for(c=0; c < vol->width; c++){
    for(r=0; r < vol->height; r++){
      for(s=0; s < vol->depth; s++){

	if(maskflag){
	  // Must be in the mask if using a mask
	  maskval = MRIgetVoxVal(mask,c,r,s,0);
	  if(maskval < 0.5)  continue;
	}
	val = MRIFseq_vox(vol,c,r,s,frame);

	// Check the sign
	if(signid == -1 && val > 0) continue;
	if(signid == +1 && val < 0) continue;

	// Get value, convert from log10 if needed
	val = fabs(val);
	if(log10flag) val = pow(10,-val);

	p[np] = val;
	np++;
      }
    }
  }
  printf("MRIfdr2vwth: np = %d, nv = %d\n",np,Nv);

  // Check that something met the match criteria, 
  // otherwise return an error
  if(np==0){
    printf("ERROR: MRIfdr2vwth: no matching voxels found\n");
    free(p);
    return(1);
  }

  *vwth = fdr2vwth(p,np,fdr);
  free(p);

  if(ovol == NULL){
    if(log10flag) *vwth = -log10(*vwth);
    return(0);
  }

  // Perform the thresholding
  for(c=0; c < vol->width; c++){
    for(r=0; r < vol->height; r++){
      for(s=0; s < vol->depth; s++){

	if(maskflag){
	  maskval = MRIgetVoxVal(mask,c,r,s,0);
	  if(maskval < 0.5){
	    // Set to null if out of mask
	    MRIFseq_vox(ovol,c,r,s,0) = valnull;
	    continue;
	  }
	}

	val = MRIFseq_vox(vol,c,r,s,frame);

	if(signid == -1 && val > 0){
	  // Set to null if wrong sign
	  MRIFseq_vox(ovol,c,r,s,0) = valnull;
	  continue;
	}
	if(signid == +1 && val < 0){
	  // Set to null if wrong sign
	  MRIFseq_vox(ovol,c,r,s,0) = valnull;
	  continue;
	}
	
	val = fabs(val);
	if(log10flag) val = pow(10,-val);
	
	if(val > *vwth){
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
  if(log10flag) *vwth = -log10(*vwth);

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
  if(mask != NULL){
    // count number of points in the mask
    nmask = 0;
    for(c=0;c<mri->width;c++){
      for(r=0;r<mri->height;r++){
	for(s=0;s<mri->depth;s++){
	  if(MRIgetVoxVal(mask,c,r,s,0) > 0.5) nmask++;
	}
      }
    }
    //printf("Number of voxels in the mask %d\n",nmask);
    if(nmask == 0){
      printf("ERROR: no voxels in mask\n");
      return(NULL);
    }
    UseMask = 1;
  }
  else{
    // Otherwise use all voxels/vertices
    nmask = mri->width * mri->height * mri->depth;
    UseMask = 0;
  }

  // Allocate the covariance matrix
  M = MatrixAlloc(mri->nframes,mri->nframes,MATRIX_REAL);

  // Compute the covariance matrix
  for(f1=0;f1<mri->nframes;f1++){
    //printf("f1 = %d\n",f1);
    for(f2=f1;f2<mri->nframes;f2++){
      sum = 0;
      for(c=0;c<mri->width;c++){
	for(r=0;r<mri->height;r++){
	  for(s=0;s<mri->depth;s++){
	    if(UseMask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
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
  if(mask != NULL){
    UseMask = 1;
    for(c=0;c<D->width;c++)
      for(r=0;r<D->height;r++)
	for(s=0;s<D->depth;s++)
	  if(UseMask && MRIgetVoxVal(mask,c,r,s,0) > 0.5) nmask++;
  }
  else{
    UseMask = 0;
    nmask = nvoxels;
  }

  M = MRIcovarianceMatrix(D, mask);
  if(M==NULL) return(1);
  //MatrixWriteTxt("cvm.dat",M);

  // Compute the SVD of the Temporal Cov Matrix
  S2 = RVectorAlloc(D->nframes, MATRIX_REAL) ;
  *pU = MatrixCopy(M,NULL); // It's done in-place so make a copy
  VV = MatrixSVD(*pU, S2, NULL); // M = U*S2*VV, VV = U';
  //MatrixWriteTxt("s2.dat",S2);

  *pS = RVectorAlloc(D->nframes, MATRIX_REAL) ;
  for(c=1; c <= dim; c++){
    S2->rptr[1][c] *= nmask; // Needs this
    (*pS)->rptr[1][c] = sqrt(S2->rptr[1][c]);
    // Exclude dims for which the singular value is very small
    if(S2->rptr[1][c] < EPSILON) break;
  }
  dim_real = c-1;

  // Compute U*inv(S) (note square root to go from S2 to S)
  UinvS = MatrixAlloc(D->nframes,dim_real,MATRIX_REAL);
  for(c=1; c <= dim_real; c++) {
    for(f=1; f <= D->nframes; f++)
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
  for(c=0;c<D->width;c++){
    for(r=0;r<D->height;r++){
      for(s=0;s<D->depth;s++){
	if(UseMask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	for(f=0; f < D->nframes; f++)
	  Fd->rptr[1][f+1] = MRIgetVoxVal(D,c,r,s,f);
	MatrixMultiply(Fd,UinvS,Fv);
	for(f=0; f < dim_real; f++){
	  MRIsetVoxVal(*pV,c,r,s,f,Fv->rptr[1][f+1]);
	  sum2[f] += (Fv->rptr[1][f+1])*(Fv->rptr[1][f+1]);
	}
      }
    }
  }

  /* Normalize to a unit vector */
  for(c=0;c<D->width;c++){
    for(r=0;r<D->height;r++){
      for(s=0;s<D->depth;s++){
	if(UseMask && MRIgetVoxVal(mask,c,r,s,0) < 0.5) continue;
	for(f=0; f < dim_real; f++){
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
  for(n=1; n <= Spca->cols; n++)
    totvar += (Spca->rptr[1][n] * Spca->rptr[1][n]);

  vsum = 0.0;
  for(n=1; n <= Spca->cols; n++){
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
  if(fp == NULL){
    printf("ERROR: opening %s\n",fname);
    return(1);
  }
  PrintPCAStats(fp, Spca);
  return(0);
}
/*---------------------------------------------------------------
  MRIsqrt() - computes sqrt(fabs(v))
  ---------------------------------------------------------------*/
MRI *MRIsqrt(MRI *invol, MRI *outvol)
{
  int c,r,s,f;
  double v;

  if(outvol == NULL){
    outvol = MRIallocSequence(invol->width, invol->height,
	      invol->depth,MRI_FLOAT, invol->nframes);
    MRIcopyHeader(invol,outvol);
    outvol->type = MRI_FLOAT;
  }
  // Should check that the dims are the same

  for(c=0;c<invol->width;c++){
    for(r=0;r<invol->height;r++){
      for(s=0;s<invol->depth;s++){
	for(f=0; f < invol->nframes; f++){
	  v = MRIgetVoxVal(invol,c,r,s,f);
	  MRIsetVoxVal(outvol,c,r,s,f,sqrt(fabs(v)));
	}
      }
    }
  }

  return(outvol);
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
  for(c=0;c<vol1->width;c++){
    for(r=0;r<vol1->height;r++){
      for(s=0;s<vol1->depth;s++){
	for(f=0; f < vol1->nframes; f++){

	  v1 = MRIgetVoxVal(vol1,c,r,s,f);
	  v2 = MRIgetVoxVal(vol2,c,r,s,f);
	  if(maxdiff < fabs(v1-v2)){
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

  if(dst==NULL){
    dst = MRIallocSequence(src->width,src->height,src->depth, 
			   MRI_FLOAT,src->nframes) ;
    MRIcopyHeader(src, dst);
  }

  for(c=0; c < src->width; c++){
    for(r=0; r < src->height; r++){
      for(s=0; s < src->depth; s++){
	for(f=0; f < src->nframes; f++){
	  v = MRIgetVoxVal(src,c,r,s,f);
	  MRIsetVoxVal(dst,c,r,s,f,v*vconst);
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
  if(!mask){
    mask = MRIcloneBySpace(mri,1);
    MRIclear(mask);
    premask = 0;
  }

  for(c=0; c < mri->width; c++){
    for(r=0; r < mri->height; r++){
      for(s=0; s < mri->depth; s++){
	if(premask){
	  m = MRIgetVoxVal(mask,c,r,s,0);
	  if(m < 0.5) continue;
	}
	n = 0;
	for(f=0; f < mri->nframes; f++){
	  val = MRIgetVoxVal(mri,c,r,s,f);
	  if(fabs(val) > thresh) n++;
	}
	if(n == mri->nframes) MRIsetVoxVal(mask,c,r,s,0,1);
	else                  MRIsetVoxVal(mask,c,r,s,0,0);
      }
    }
  }
  return(mask);
}

/* MRImakeVox2VoxReg() - takes a target volume, a movable volume, a
   registration method, and an optional registration filename, and
   will make an mriTransform in which A is the target non-comformed
   index space, ARAS is the target tkRegRAS space, B is the movable
   tkRegRAS space (registered), and B is the movable non-conformed
   index space. If *transform is NULL, this will make a new
   mriTransform and return it, otherwise it will set the values in the
   *transform passed in. */
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

  /* We need to build three matrices to fill in the three transforms
     in the mriTransform object:

     targ (A) ----------> tkRegRAS unregistered (ARAS)
                             |
                             |
                             V
      mov (B) ----------> tkRegRAS registered (BRAS)	    

     A->ARAS and B->BRAS take each volume from native index space to
     tkRegRAS space, an RAS space in which the center of the world is
     the center of the volume. We can build these manually by looking
     at the header info in the MRIs. They look like this:

                  [ -xsize    0     0     (xsize*(width-1))/2  ]
                  [    0      0   zsize  -(zsize*(depth-1))/2  ]
                  [    0   -ysize   0     (ysize*(height-1))/2 ]
                  [    0      0     0              1           ]

     For the ARAS->BRAS matrix, we're either going to use a matrix
     from a registration file or one generated by MRItkRegMtx.
  */

  /* Create the targ A->RAS matrix. */
  targ_idx_to_tkregras = MRIxfmCRS2XYZtkreg(targ);
  if(NULL==targ_idx_to_tkregras){
    print("ERROR: MRImakeVox2VoxReg: Couldn't create targ_idx_to_tkregras\n");
    goto error;
  }

  /* Create the mov A->RAS matrix. */
  mov_idx_to_tkregras = MRIxfmCRS2XYZtkreg(mov);
  if(NULL==mov_idx_to_tkregras){
    print("ERROR: MRImakeVox2VoxReg: Couldn't create mov_idx_to_tkregras\n");
    goto error;
  }

  /* Now we build the ARAS->BRAS matrix, which is the
     registration. Switch on our registration type. */
  switch(regtype)
    {
    case VOX2VOXREGTYPE_FILE:
    case VOX2VOXREGTYPE_FIND:

      /* If we're reading a file, copy the file from the input or
	 generate one from our data file location. */
      if(VOX2VOXREGTYPE_FILE==regtype){
	strncpy(fullregname, regname, sizeof(fullregname));
      }
      else if(VOX2VOXREGTYPE_FIND==regtype){
	/* Copy the movable volume name and find the last / in the
	   file name. From there, copy in "register.dat" for our file
	   name. */
	strncpy(regpath, mov->fname, sizeof(regpath));
	cur_char = regpath;
	base_end = regpath;
	while(NULL!=cur_char && '\0' != *cur_char){
	  if('/' == *cur_char) 
	    base_end = cur_char;
	  cur_char++;
	}
	*base_end = '\0';
	snprintf(fullregname,sizeof(fullregname),
		 "%s/%s",regpath,"register.dat");
      }

      /* Check that the file exists. */
      err = stat(fullregname,&file_info);
      if(0!=err){
	printf("ERROR: MRImakeVox2VoxReg: Couldn't find registration %s\n",
	       fullregname);
	goto error;
      }
      
      /* Check if it's a regular file. */
      if(!S_ISREG(file_info.st_mode)){
	printf("ERROR: MRImakeVox2VoxReg: Couldn't open registration %s\n",
	       fullregname);
	goto error;
      }

      /* Read the registration */
      reg_info = StatReadRegistration(fullregname);
      if(NULL==reg_info){
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
  if(*transform==NULL){
    trnscode=Trns_New(transform);
    if(Trns_tErr_NoErr!=trnscode){
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
  if(Trns_tErr_NoErr!=trnscode){
    printf("ERROR: MRImakeVox2VoxReg: Error copying A to RAS\n");
    goto error;
  }
  trnscode=Trns_CopyBtoRAS(*transform,mov_idx_to_tkregras);
  if(Trns_tErr_NoErr!=trnscode){
    printf("ERROR: MRImakeVox2VoxReg: Error copying B to RAS\n");
      goto error;
  }
  trnscode=Trns_CopyARAStoBRAS(*transform,targ_tkregras_to_mov_tkregras);
  if(Trns_tErr_NoErr!=trnscode){
    printf("ERROR: MRImakeVox2VoxReg: Error copying ARAS to BRAS\n");
    goto error;
  }
  
  goto cleanup;
 error:

  if(0==retcode)
    retcode = -1;

 cleanup:

  if(NULL!=targ_idx_to_tkregras)
    MatrixFree(&targ_idx_to_tkregras);

  if(NULL!=mov_idx_to_tkregras)
    MatrixFree(&mov_idx_to_tkregras);

  if(NULL!=reg_info)
    StatFreeRegistration(&reg_info);
  
  if(NULL!=targ_tkregras_to_mov_tkregras)
    MatrixFree(&targ_tkregras_to_mov_tkregras);

  return(retcode);
}
