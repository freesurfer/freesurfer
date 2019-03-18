#include "fastmarching.h"

static MRI_REGION*mriFindLabel(MRI *mri,int label,int offset) {
  MRI_REGION *region;
  int i,j,k,nlabels;
  double xmin,ymax,zmin,xmax,ymin,zmax;
  double xw,yw,zw;

  region=(MRI_REGION*)calloc(1,sizeof(MRI_REGION));

  xmin=ymin=zmin=100000;
  xmax=ymax=zmax=-100000;
  for (nlabels=k=0;k<mri->depth;k++)
    for (j=0;j<mri->height;j++)
      for (i=0;i<mri->width;i++)
        if (MRIvox(mri,i,j,k)==label) {
          nlabels++;
          MRIvoxelToSurfaceRAS(mri, i, j, k, &xw, &yw, &zw) ;
          if (xmin>xw) xmin=xw;
          if (ymin>yw) ymin=yw;
          if (zmin>zw) zmin=zw;
          if (xmax<xw) xmax=xw;
          if (ymax<yw) ymax=yw;
          if (zmax<zw) zmax=zw;
        }

  if (nlabels==0) ErrorExit(1,"No labels %d could be found in the volume\n",nlabels);

  region->x=(int)floor(xmin-offset);
  region->y=(int)floor(ymin-offset);
  region->z=(int)floor(zmin-offset);

  region->dx=(int)ceil(xmax-xmin+2*offset);
  region->dy=(int)ceil(ymax-ymin+2*offset);
  region->dz=(int)ceil(zmax-zmin+2*offset);

  return region;
}

// the source mri mri_src is a float distance map with zero values inside.
// the resolution is mri_src->xsize
void MRISextractOutsideDistanceMap(MRIS *mris, MRI *mri_src, int label , int offset, float resolution, float max_distance) {
  MRI *mri_distance, *mri_fastmarching;

  /* resolution*/
  if ( resolution < 1.0f ) resolution=1.0f;
  /* distance */
  if ( max_distance < 1.0f ) max_distance=1.0f;
  /* offset */
  if (offset < 0 ) offset = 0;

  /* find the region of interest */
  MRI_REGION *region=mriFindLabel(mri_src,label,offset);

  /* compute */
  mri_distance=MRISbinarizeVolume(mris, region, resolution, 5.0f);

  mri_fastmarching=MRIclone(mri_distance,NULL);

  /* set values to zero */
  mapMRI_XYZ(mri_fastmarching,x,y,z) MRIFvox(mri_fastmarching,x,y,z)=0.0f;

  FastMarching<+1> fastmarching_out(mri_fastmarching, NULL);
  fastmarching_out.SetLimit(2*resolution*max_distance);
  fastmarching_out.InitForOutsideMatch(mri_distance,mri_src,label);
  fastmarching_out.Run(2*resolution*max_distance);

  mapMRI_XYZ(mri_fastmarching,x,y,z)
  if (MRIFvox(mri_distance,x,y,z)<0)
    MRIFvox(mri_fastmarching,x,y,z)=-0.1;

  /* extract values */
  for (int p=0;p<mris->nvertices;p++) {
    int i,j,k;
    float val;
    VERTEX *v=&mris->vertices[p];
    v->curv=1.0f;
    for (int n=0;n<5;n++) { //first valid point in the direction white->pial
      float t=n*0.1;
      i=iVOL(mri_fastmarching,((1-t)*v->whitex+t*v->pialx));
      j=jVOL(mri_fastmarching,((1-t)*v->whitey+t*v->pialy));
      k=kVOL(mri_fastmarching,((1-t)*v->whitez+t*v->pialz));

      if (i<0||i>mri_fastmarching->width-1||j<0||j>mri_fastmarching->height-1||k<0||k>mri_fastmarching->depth-1) continue;
      val=MRIFvox(mri_fastmarching,i,j,k)/resolution;
      if (val>=0) {
        v->curv=MIN(val,max_distance)/max_distance;
        break;
      }
    }
  }

  MRIfree(&mri_fastmarching);
  MRIfree(&mri_distance);
  free(region);

  return;
}
