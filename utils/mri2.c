/*-------------------------------------------------------------------
  Name: mri2.c
  Author: Douglas N. Greve
  $Id: mri2.c,v 1.11 2004/12/20 02:02:40 kteich Exp $
  Purpose: more routines for loading, saving, and operating on MRI 
  structures.
  -------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
  CRS in the source (ie, it is a vox2vox).

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

  if(src->nframes != targ->nframes){
    printf("ERROR: MRIvol2vol: source and target have different number "
	   "of frames\n");
    return(1);
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


