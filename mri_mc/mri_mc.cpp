#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
/**
 * @brief generates a marching cubes triangulation
 *
 * Generates a marching cubes triangulation that is topologically consistent
 * with a specific choice of connectivity. By default, this program will only 
 * keep the main connected component.
 *
 */
/*
 * Original Author: Florent Segonne
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


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>

#include "fio.h"
#include "const.h"
#include "diag.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "MRIio_old.h"
#include "mri.h"
#include "mrisurf.h"
#include "version.h"
#include "tags.h"
#include "gca.h"
#include "MC.h"

#define MAXFACES    3000000
#define MAXVERTICES 1500000
#define FACEINCREASE 1.2
#define VERTEXINCREASE 1.2

const char *Progname;

typedef struct quad_face_type_ {
  int imnr,i,j,f;
  int num;
  int v[4];
}
quad_face_type;

typedef struct quad_vertex_type_ {
  float imnr,i,j;
  int num;
  int f[9];
}
quad_vertex_type;


typedef struct tesselation_parms_ {
  /*labeled volume*/
  MRI *mri;

  //for the marching cube tesselation
  int MC;

  int imgsize,width,height;

  unsigned char*** im;

  /*label information*/
  int number_of_labels;
  int* label_values;
  int *imin,*imax,*jmin,*jmax,*kmin,*kmax;
  int current_label;
  int ind;
  int xmin,xmax,ymin,ymax,zmin,zmax;

  /*final surfaces*/
  /*The surfaces are saved into a table of surfaces*/
  MRIS **mris_table;

  /*all labels are used for the tesselation if(all_flag) */
  int all_flag;

  /*to compute the tesselation consistent with a n topology*/
  int connectivity;

  /*face information*/
  int face_index;
  quad_face_type *face;
  int maxfaces;
  int *face_index_table0;
  int *face_index_table1;

  /*vertex information*/
  int vertex_index;
  quad_vertex_type *vertex;
  int maxvertices;
  //  int *vertex_index_table;

  /*used for the debugging*/
  char fname[STRLEN];

}
tesselation_parms;


static int downsample = 0 ;

/*initialization of tesselation_parms*/
/*note that not all the fields are allocated*/
void initTesselationParms(tesselation_parms *parms) {
  int i,j,k,n;
  int depth,width,height;
  int *xmin,*ymin,*zmin,*xmax,*ymax,*zmax;
  int val,label_nbr,*labels;
  unsigned char ***im;
  MRI *mri=parms->mri;

  depth=mri->depth;
  width=mri->width;
  height=mri->height;

  parms->width=2*MAX(width,height)+1;
  parms->height=parms->width;
  parms->imgsize=parms->width*parms->height;

  /* allocate and initialize the image*/
  im=(unsigned char***)malloc(sizeof(unsigned char**)*depth);
  if (!im)
    ErrorExit(ERROR_NO_MEMORY,"MRIStesselate: could not allocate the image");
  for (k=0;k<depth;k++) {
    im[k] = (unsigned char **)malloc(height*sizeof(char *));
    if (!im[k])
      ErrorExit(ERROR_NO_MEMORY,"MRIStesselate: could not allocate the image");
    for (j=0;j<height;j++) {
      im[k][j] = (unsigned char *)calloc(width,sizeof(char));
      if (!im[k][j])
        ErrorExit(ERROR_NO_MEMORY,
                  "MRIStesselate: could not allocate the image");
    }
  }
  parms->im=im;

  label_nbr=parms->number_of_labels;
  labels=parms->label_values;

  xmin=(int*)malloc(label_nbr*sizeof(int));
  ymin=(int*)malloc(label_nbr*sizeof(int));
  zmin=(int*)malloc(label_nbr*sizeof(int));
  xmax=(int*)malloc(label_nbr*sizeof(int));
  ymax=(int*)malloc(label_nbr*sizeof(int));
  zmax=(int*)malloc(label_nbr*sizeof(int));

  if ((!xmin) || (!ymin) || (!zmin) || (!xmax) || (!ymax) || (!zmax))
    ErrorExit(ERROR_NO_MEMORY,"MRIStesselate: table of labels");

  for (n=0;n<label_nbr;n++) {
    xmin[n]=ymin[n]=zmin[n]=100000;
    xmax[n]=ymax[n]=zmax[n]=0;
  }

  for (k=1;k<depth-1;k++)
    for (j=1;j<height-1;j++)
      for (i=1;i<width-1;i++)
        im[k][j][i]=MRIvox(mri,i,j,k);

  for (k=1;k<mri->depth-1;k++)
    for (j=1;j<mri->height-1;j++)
      for (i=1;i<mri->width-1;i++) {
        val=MRIvox(mri,i,j,k);
        for (n=0;n<label_nbr;n++)
          if (val==labels[n] || (parms->all_flag && val)) {
            if (i<xmin[n])  xmin[n]=i;
            if (j<ymin[n])  ymin[n]=j;
            if (k<zmin[n])  zmin[n]=k;

            if (i>xmax[n])  xmax[n]=i;
            if (j>ymax[n])  ymax[n]=j;
            if (k>zmax[n])  zmax[n]=k;
          }
      }

  /*be careful: the region is expanded by one voxel in the MRI referential*/
  for (n=0;n<label_nbr;n++) {
    xmin[n]--;
    ymin[n]--;
    zmin[n]--;
    xmax[n]++;
    ymax[n]++;
    zmax[n]++;
  }

  parms->imin=xmin;
  parms->jmin=ymin;
  parms->kmin=zmin;
  parms->imax=xmax;
  parms->jmax=ymax;
  parms->kmax=zmax;

  parms->all_flag=0;
}
#if 1
static int freeTesselationParms(tesselation_parms **parms) {
  int j,k;
  tesselation_parms *tmp=*parms;
  int depth,width,height;
  MRI *mri=tmp->mri;

  depth=mri->depth;
  height=mri->height;
  width=mri->width;

  for (k=0;k<depth;k++) {
    for (j=0;j<height;j++)
      free(tmp->im[k][j]);
    free(tmp->im[k]);
  }
  free(tmp->im);


  free(tmp->imin);
  free(tmp->jmin);
  free(tmp->kmin);
  free(tmp->imax);
  free(tmp->jmax);
  free(tmp->kmax);

  free(*parms);
  return(NO_ERROR) ;
}
#endif
void allocateTesselation(tesselation_parms *parms) {
  int imgsize=parms->imgsize;
  parms->face = (quad_face_type *)lcalloc(MAXFACES,sizeof(quad_face_type));
  parms->maxfaces=MAXFACES;
  parms->face_index_table0 = (int *)lcalloc(6*imgsize,sizeof(int));
  parms->face_index_table1 = (int *)lcalloc(6*imgsize,sizeof(int));

  parms->vertex = 
    (quad_vertex_type *)lcalloc(MAXVERTICES,sizeof(quad_vertex_type));
  parms->maxvertices=MAXVERTICES;
  //  parms->vertex_index_table = (int *)lcalloc(8*imgsize,sizeof(int));

  parms->face_index=0;
  parms->vertex_index=0;

  if ((!parms->face) || (!parms->face_index_table0)
      || (!parms->face_index_table1)
      || (!parms->vertex))
    //     || (!parms->vertex_index_table))
    ErrorExit(ERROR_NO_MEMORY,"MRIStesselate: local tesselation tables");
}

#if 0
static void freeTesselation(tesselation_parms *parms) {
  free(parms->face);
  free(parms->face_index_table0);
  free(parms->face_index_table1);

  free(parms->vertex);
  //  free(parms->vertex_index_table);
}
#endif

MRI* preprocessingStep(tesselation_parms *parms) {
  int i,j,k,width,height,depth,n;
  int xmin,ymin,zmin,xmax,ymax,zmax;
  MRI *mri;

  width=parms->mri->width;
  height=parms->mri->height;
  depth=parms->mri->depth;

  mri=MRIclone(parms->mri,NULL);

  xmin=ymin=zmin=100000;
  xmax=ymax=zmax=0;

  if (!parms->all_flag) {
    for (k=0;k<depth;k++)
      for (j=0;j<height;j++)
        for (i=0;i<width;i++)
          for (n=0;n<parms->number_of_labels;n++)
            if (MRIvox(parms->mri,i,j,k)==parms->label_values[n]) {
              MRIvox(mri,i,j,k)=1;
              if (i<xmin)  xmin=i;
              if (j<ymin)  ymin=j;
              if (k<zmin)  zmin=k;

              if (i>xmax)  xmax=i;
              if (j>ymax)  ymax=j;
              if (k>zmax)  zmax=k;
            }
  } else {
    for (k=0;k<depth;k++)
      for (j=0;j<height;j++)
        for (i=0;i<width;i++)
          if (MRIvox(parms->mri,i,j,k)) {
            MRIvox(mri,i,j,k)=1;
            if (i<xmin)  xmin=i;
            if (j<ymin)  ymin=j;
            if (k<zmin)  zmin=k;

            if (i>xmax)  xmax=i;
            if (j>ymax)  ymax=j;
            if (k>zmax)  zmax=k;
          }
  }

  parms->xmin=xmin-1;
  parms->ymin=ymin-1;
  parms->zmin=zmin-1;

  parms->xmax=xmax+1;
  parms->ymax=ymax+1;
  parms->zmax=zmax+1;
  return mri;

}

void reallocateVertices(tesselation_parms *parms) {
  quad_vertex_type *tmp;
  int k,n;
  int vertex_index=parms->vertex_index;

  parms->maxvertices=(int)(parms->maxvertices*VERTEXINCREASE);
  tmp=(quad_vertex_type*)lcalloc(parms->maxvertices,sizeof(quad_vertex_type));
  if (!tmp)
    ErrorExit(ERROR_NOMEMORY, "%s: max vertices %d exceeded",
              Progname,parms->maxvertices) ;
  for (k=0;k<vertex_index;k++) {
    tmp[k].imnr=parms->vertex[k].imnr;
    tmp[k].i=parms->vertex[k].i;
    tmp[k].j=parms->vertex[k].j;
    tmp[k].num=parms->vertex[k].num;
    for (n=0;n<9;n++)
      tmp[k].f[n]=parms->vertex[k].f[n];
  }
  free(parms->vertex);
  parms->vertex=tmp;
}
void reallocateFaces(tesselation_parms *parms) {
  quad_face_type *tmp;
  int k,n;
  int face_index=parms->face_index;

  parms->maxfaces=(int)(parms->maxfaces*FACEINCREASE);
  tmp=(quad_face_type*)lcalloc(parms->maxfaces,sizeof(quad_face_type));
  if (!tmp)
    ErrorExit(ERROR_NOMEMORY, "%s: max faces %d exceeded",
              Progname,parms->maxfaces) ;
  for (k=0;k<face_index;k++) {
    tmp[k].imnr=parms->face[k].imnr;
    tmp[k].i=parms->face[k].i;
    tmp[k].j=parms->face[k].j;
    tmp[k].f=parms->face[k].f;
    tmp[k].num=parms->face[k].num;
    for (n=0;n<4;n++)
      tmp[k].v[n]=parms->face[k].v[n];
  }

  free(parms->face);
  parms->face=tmp;
}


int saveTesselation2(tesselation_parms *parms) {
  int vno,m,n,fno;
  int nVFMultiplier=1;
  quad_face_type *face2;
  quad_vertex_type *vertex2;
  MRIS* mris;
  FACE *face;
  float x, y, z, xhi, xlo, yhi, ylo, zhi, zlo ;
  float st,ps,xx1,yy0,zz1;
  float j,i,imnr;
  double xw,yw,zw;

  /*necessary for the coord transformation*/
  ps=parms->mri->ps;
  st=parms->mri->thick;

  yy0=parms->mri->ystart;
  xx1=parms->mri->xend;
  zz1=parms->mri->zend;

  mris=MRISoverAlloc(nVFMultiplier*parms->vertex_index,nVFMultiplier*parms->face_index
                     ,parms->vertex_index,parms->face_index);

  MRIScopyVolGeomFromMRI(mris, parms->mri) ;
  fprintf(stderr,"\n(surface with %d faces and %d vertices)...",
          parms->face_index,parms->vertex_index);

  mris->type=MRIS_TRIANGULAR_SURFACE;
  /*first init vertices*/
  for (vno=0;vno<mris->nvertices;vno++) {
    VERTEX_TOPOLOGY* const vertext = &mris->vertices_topology[vno] ;

    vertex2= &parms->vertex[vno];
    i=vertex2->i;
    j=vertex2->j;
    imnr=vertex2->imnr;

    //      x = xx1-(j-0.5)*ps;
    //y = yy0+(imnr-0.5)*st;
    //z = zz1-(i-0.5)*ps;

    MRISsurfaceRASFromVoxelCached(mris, parms->mri, i, j, imnr, &xw, &yw,&zw);
    //    MRIvoxelToWorld(parms->mri,i,j,imnr,&xw,&yw,&zw);
    x=xw;
    y=yw;
    z=zw;

    MRISsetXYZ(mris, vno,
      x,
      y,
      z);
    vertext->num=0;
  }
  /*then init faces*/
  for (m = 0 ; m < parms->face_index ; m ++) {

    face2=&parms->face[m];

    mris->faces[m].v[0] = face2->v[0] ;
    mris->faces[m].v[1] = face2->v[1] ;
    mris->faces[m].v[2] = face2->v[2] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      mris->vertices_topology[mris->faces[m].v[n]].num++;
  }

  /*allocate indices & faces for each vertex*/
  for (vno = 0 ; vno< mris->nvertices ; vno++) {
    // vertex = &mris->vertices[vno] ;
    mris->vertices_topology[vno].f =
      (int *)calloc(mris->vertices_topology[vno].num,sizeof(int));
    if (!mris->vertices_topology[vno].f)
      ErrorExit(ERROR_NOMEMORY,
                "%s: could not allocate %d faces at %dth vertex",
                Progname, mris->vertices_topology[vno].num,vno) ;
    mris->vertices_topology[vno].n =
      (uchar *)calloc(mris->vertices_topology[vno].num,sizeof(uchar));
    if (!mris->vertices_topology[vno].n)
      ErrorExit(ERROR_NOMEMORY,
                "%s: could not allocate %d indices at %dth vertex",
                Progname, mris->vertices_topology[vno].num,vno) ;
    mris->vertices_topology[vno].num = 0 ;
  }

  /*init faces for each vertex*/
  for (fno = 0 ; fno < mris->nfaces ; fno++) {
    face = &mris->faces[fno] ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      mris->vertices_topology[face->v[n]].f[mris->vertices_topology[face->v[n]].num++] = fno;
  }

  mrisCheckVertexFaceTopology(mris);

  /*necessary initialization*/
  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    mris->vertices[vno].curv = 0;
    mris->vertices[vno].origarea = -1;
    mris->vertices[vno].border = 0;

    for (n=0;n<mris->vertices_topology[vno].num;n++) {
      for (m=0;m<VERTICES_PER_FACE;m++) {
        if (mris->faces[mris->vertices_topology[vno].f[n]].v[m] == vno)
          mris->vertices_topology[vno].n[n] = m;
      }
    }
    x = mris->vertices[vno].x;
    y = mris->vertices[vno].y;
    z = mris->vertices[vno].z;
    if (x>xhi) xhi=x;
    if (x<xlo) xlo=x;
    if (y>yhi) yhi=y;
    if (y<ylo) ylo=y;
    if (z>zhi) zhi=z;
    if (z<zlo) zlo=z;
  }

  mris->xlo = xlo ;
  mris->ylo = ylo ;
  mris->zlo = zlo ;
  mris->xhi = xhi ;
  mris->yhi = yhi ;
  mris->zhi = zhi ;
  mris->xctr = (xhi+xlo)/2;
  mris->yctr = (yhi+ylo)/2;
  mris->zctr = (zhi+zlo)/2;

  mrisCompleteTopology(mris);
  MRIScomputeNormals(mris);

  mris->type = MRIS_TRIANGULAR_SURFACE; /*not so sure about that*/
  MRISsetNeighborhoodSizeAndDist(mris, 2) ;
  MRIScomputeSecondFundamentalForm(mris) ;
  MRISuseMeanCurvature(mris) ;
  mris->radius = MRISaverageRadius(mris) ;
  MRIScomputeMetricProperties(mris) ;
  MRISstoreCurrentPositions(mris) ;
  parms->mris_table[parms->ind]=mris;
  //MRISdivideLongEdges(mris, 0.3);

#if 0
  {
    int f;
    float v0[3],v1[3],d1,d2,d3;
    for (f=0;f<mris->nfaces;f++) {
      v0[0] = mris->vertices[mris->faces[f].v[2]].x - 
        mris->vertices[mris->faces[f].v[0]].x;
      v0[1] = mris->vertices[mris->faces[f].v[2]].y - 
        mris->vertices[mris->faces[f].v[0]].y;
      v0[2] = mris->vertices[mris->faces[f].v[2]].z - 
        mris->vertices[mris->faces[f].v[0]].z;
      v1[0] = mris->vertices[mris->faces[f].v[1]].x - 
        mris->vertices[mris->faces[f].v[0]].x;
      v1[1] = mris->vertices[mris->faces[f].v[1]].y - 
        mris->vertices[mris->faces[f].v[0]].y;
      v1[2] = mris->vertices[mris->faces[f].v[1]].z - 
        mris->vertices[mris->faces[f].v[0]].z;
      d1 = -v1[1]*v0[2] + v0[1]*v1[2];
      d2 = v1[0]*v0[2] - v0[0]*v1[2];
      d3 = -v1[0]*v0[1] + v0[0]*v1[1];
      if (sqrt(d1*d1+d2*d2+d3*d3)/2>0.5)
        fprintf(stderr,"\n Face #%d -> %f %d %d %d ",
                f,sqrt(d1*d1+d2*d2+d3*d3)/2,mris->faces[f].v[0],
                mris->faces[f].v[1],mris->faces[f].v[2]);
    }
  }
#endif

  return(NO_ERROR) ;
}

void generateMCtesselation(tesselation_parms * parms) {
  int i,j,k,width,imgsize,*tab1,*tab2,ref,ind,nf,p;
  int xmin,ymin,zmin,xmax,ymax,zmax;
  int vt[12],*vk1,*vk2,*vj1,*vj2,*tmp;
  int f_c[12],vind[12];
  MRI *mri;
  int vertex_index;
  int debug=0;

  fprintf(stderr,"\npreprocessing...");
  mri=preprocessingStep(parms);
  allocateTesselation(parms);
  fprintf(stderr,"done\n");

  width=mri->width;
  imgsize=mri->width*mri->height;


  xmin=parms->xmin;
  ymin=parms->ymin;
  zmin=parms->zmin;
  xmax=parms->xmax;
  ymax=parms->ymax;
  zmax=parms->zmax;

  tab1=(int*)calloc(imgsize,sizeof(int));
  tab2=(int*)calloc(imgsize,sizeof(int));

  //vertices memory
  vk1=(int*)calloc(2*imgsize,sizeof(int));
  vk2=(int*)calloc(2*imgsize,sizeof(int));
  vj1=(int*)calloc(width,sizeof(int));
  vj2=(int*)calloc(width,sizeof(int));

  f_c[0]=0;
  f_c[1]=1;
  f_c[2]=2*width;
  f_c[3]=3;
  f_c[8]=0;
  f_c[9]=1;
  f_c[10]=2*width;
  f_c[11]=3;
  f_c[4]=0;
  f_c[5]=1;
  f_c[6]=0;
  f_c[7]=1;

  fprintf(stderr,"starting generation of surface...");

  for (k=zmin;k<zmax;k++) {
    if (!(k % 10))
      fprintf(stderr,"\n   slice nb %d...",k);
    for (j=ymin;j<ymax;j++) {
      for (i=xmin;i<xmax;i++) {

        ind=i+width*j;
        ref=0;
        if (((k+1)<mri->depth) && MRIvox(mri,i,j,k+1))
          ref+=16;
        if (((k+1)<mri->depth) && ((i+1)<mri->width) && MRIvox(mri,i+1,j,k+1))
          ref+=32;
        if (((k+1)<mri->depth) && ((j+1)<mri->height) && MRIvox(mri,i,j+1,k+1))
          ref+=64;
        if (((k+1)<mri->depth) && 
            ((j+1)<mri->height) && 
            ((i+1)<mri->width) && 
            MRIvox(mri,i+1,j+1,k+1))
          ref+=128;

        tab2[ind]=(ref/16);
        ref+=tab1[ind]; //this is the indice of the cube


        nf=0;
        switch (parms->connectivity) {
        case 1:
          while (MC6p[ref][3*nf]>=0) nf++;
          break;
        case 2:
          while (MC18[ref][3*nf]>=0) nf++;
          break;
        case 3:
          while (MC6[ref][3*nf]>=0) nf++;
          break;
        default:
          while (MC26[ref][3*nf]>=0) nf++;
          break;
        }
        if (nf==0) continue;
        //fprintf(stderr,"\n %d %d. ",j,i);
        //fprintf(stderr," ref=%d, ind=%d,",ref,ind);
        //fprintf(stderr,"%d face(s)\n",nf);

        memset(vt,0,12*sizeof(int));
        memset(vind,0,12*sizeof(int));

        switch (parms->connectivity) {
        case 1:
          for (p=0;p<3*nf;p++) vt[MC6p[ref][p]]++;
          break;
        case 2:
          for (p=0;p<3*nf;p++) vt[MC18[ref][p]]++;
          break;
        case 3:
          for (p=0;p<3*nf;p++) vt[MC6[ref][p]]++;
          break;
        default:
          for (p=0;p<3*nf;p++) vt[MC26[ref][p]]++;
          break;
        }


        //find references of vertices and eventually allocate them!
        for (p=0;p<4;p++) //find the vertex number
          if (vt[p])
            vind[p]=vk1[2*ind+f_c[p]];
        for (p=4;p<6;p++)
          if (vt[p])
            vind[p]=vj1[i+f_c[p]];
        if (vt[6])         //already created
        {
          vind[6]=vj2[i];
        }
        if (vt[7]) //create a new vertex number and save it into v7 and vj2
        {
          vertex_index=parms->vertex_index;
          //fprintf(stderr,"c(7):%d- ",vertex_index);
          if (vertex_index >= parms->maxvertices-1)
            reallocateVertices(parms);
          vind[7]=vertex_index;
          parms->vertex[vertex_index].imnr = k+0.5;
          parms->vertex[vertex_index].i = i+1;
          parms->vertex[vertex_index].j = j+1;
          parms->vertex[vertex_index].num = 0;
          vj2[i+1]=vertex_index;
          parms->vertex_index++;

          //     printf(stderr," !!!!*!!! ");
        }
        if (vt[8]) //already created
        {
          vind[8]=vk2[2*ind+f_c[8]];
        }
        if (vt[9]) //already created
        {
          vind[9]=vk2[2*ind+f_c[9]];
        }
        if (vt[10]) //create a new vertex number and save it into vk2
        {
          vertex_index=parms->vertex_index;
          //fprintf(stderr,"c(10):%d- ",vertex_index);
          if (vertex_index >= parms->maxvertices-1)
            reallocateVertices(parms);
          vind[10]=vertex_index;
          parms->vertex[vertex_index].imnr = k+1;
          parms->vertex[vertex_index].i = i+0.5;
          parms->vertex[vertex_index].j = j+1;
          parms->vertex[vertex_index].num = 0;
          vk2[2*ind+f_c[10]]=vertex_index;
          parms->vertex_index++;
        }
        if (vt[11]) //create a new vertex number and save it into vk2
        {
          vertex_index=parms->vertex_index;
          //fprintf(stderr,"c(11):%d- ",vertex_index);
          if (vertex_index >= parms->maxvertices-1)
            reallocateVertices(parms);
          vind[11]=vertex_index;
          parms->vertex[vertex_index].imnr = k+1;
          parms->vertex[vertex_index].i = i+1;
          parms->vertex[vertex_index].j = j+0.5;
          parms->vertex[vertex_index].num = 0;
          vk2[2*ind+f_c[11]]=vertex_index;
          parms->vertex_index++;
        }
        //now create faces
        for (p=0;p<nf;p++) {
          int face_index=parms->face_index;
          if (face_index >= parms->maxfaces-1)
            reallocateFaces(parms);

          switch (parms->connectivity) {
          case 1:
            parms->face[face_index].v[0] = vind[MC6p[ref][3*p]];
            parms->face[face_index].v[1] = vind[MC6p[ref][3*p+1]];
            parms->face[face_index].v[2] = vind[MC6p[ref][3*p+2]];
            break;
          case 2:
            parms->face[face_index].v[0] = vind[MC18[ref][3*p]];
            parms->face[face_index].v[1] = vind[MC18[ref][3*p+1]];
            parms->face[face_index].v[2] = vind[MC18[ref][3*p+2]];
            break;
          case 3:
            parms->face[face_index].v[0] = vind[MC6[ref][3*p]];
            parms->face[face_index].v[1] = vind[MC6[ref][3*p+1]];
            parms->face[face_index].v[2] = vind[MC6[ref][3*p+2]];
            break;
          default:
            parms->face[face_index].v[0] = vind[MC26[ref][3*p]];
            parms->face[face_index].v[1] = vind[MC26[ref][3*p+1]];
            parms->face[face_index].v[2] = vind[MC26[ref][3*p+2]];
            break;
          }
          parms->face_index++;
        }
        //fprintf(stderr,"%d %d,",parms->face_index,parms->vertex_index);
        //fprintf(stderr,"?");

        //fprintf(stderr,"\n ");
        //for(p=0;p<12;p++) fprintf(stderr,": p=%d vt=%d v=%d ",
        //p,vt[p],vind[p]);

        if (debug)
          pause();

      }
      tmp=vj1;
      vj1=vj2;
      vj2=tmp;
      memset(vj2,-1,width*sizeof(int));
    }
    tmp=tab1;
    tab1=tab2;
    tab2=tmp;
    memset(tab2,-1,imgsize*sizeof(int));

    tmp=vk1;
    vk1=vk2;
    vk2=tmp;
    memset(vk2,-1,2*imgsize*sizeof(int));

  }
  free(tab1);
  free(tab2);
  free(vj1);
  free(vj2);
  free(vk1);
  free(vk2);
  MRIfree(&mri);
  fprintf(stderr,"\nconstructing final surface...");
  saveTesselation2(parms);
  fprintf(stderr,"done\n");
}

int main(int argc, char *argv[]) {
  tesselation_parms *parms;
  MRIS **mris_table, *mris,*mris_corrected;
  MRI *mri, *mri_orig;

  std::string cmdline = getAllInfo(argc, argv, "mri_mc");

  Progname=argv[0];

  if (argc > 1 && (stricmp(argv[1], "-d") == 0)) {
    downsample = atoi(argv[2]) ;
    argc -= 2;
    argv += 2 ;
    printf("downsampling input volume %d times\n", downsample) ;
  }

  if (argc < 4) {
    fprintf(stderr,"\n\nUSAGE: mri_mc input_volume "
            "label_value output_surface [connectivity]");
    fprintf(stderr,
            "\noption connectivity: 1=6+,2=18,3=6,4=26 (default=1)\n\n");
    exit(-1);
  }

  parms=(tesselation_parms*)calloc(1,sizeof(tesselation_parms));
  if (!parms)
    ErrorExit(ERROR_NOMEMORY, "tesselation parms\n") ;
  mri=MRIread(argv[1]);
  if (downsample > 0) {
    MRI *mri_tmp ;
    mri_tmp = MRIdownsample2(mri, NULL) ;
    MRIfree(&mri) ;
    mri = mri_tmp ;
  }
  {
    MRI *mri_tmp ;
    mri_tmp = MRIalloc(mri->width+2, mri->height+2, mri->depth+2, mri->type) ;
    MRIextractInto(mri, mri_tmp, 
                   0, 0, 0, 
                   mri->width, mri->height, mri->depth, 
                   1, 1, 1) ;
    //MRIfree(&mri) ;
    mri_orig = mri;
    mri = mri_tmp ;
  }
  MRIreInitCache(mri);
  if (mri->type != MRI_UCHAR) {
    MRI *mri_tmp ;
    float min_val, max_val ;

    MRIvalRange(mri, &min_val, &max_val) ;
    if (min_val < 0 || max_val > 255)
      ErrorExit
        (ERROR_UNSUPPORTED, 
         "%s: input volume (val range [%2.1f %2.1f]) must be "
         "convertible to UCHAR",
         Progname, min_val, max_val) ;
    printf("changing type of input volume to 8 bits/voxel...\n") ;
    mri_tmp = MRIchangeType(mri, MRI_UCHAR, 0.0, 0.999, TRUE) ;
    MRIfree(&mri) ;
    mri = mri_tmp ;
  }

  parms->mri=mri;

  parms->number_of_labels=1; //only one single label
  parms->label_values=(int*)malloc(sizeof(int));
  parms->label_values[0]=atoi(argv[2]);//label;
  parms->ind=0;
  mris_table=(MRIS**)malloc(sizeof(MRIS*)); //final surface information
  parms->mris_table=mris_table;
  if ((!parms->label_values) || (!mris_table))
    ErrorExit(ERROR_NOMEMORY, "labels/surfaces tables\n") ;

  if (argc==5) parms->connectivity=atoi(argv[4]);//connectivity;
  else parms->connectivity=1;

  initTesselationParms(parms);

  generateMCtesselation(parms);

  free(parms->label_values);
  mris=parms->mris_table[0];
  free(parms->mris_table);
  freeTesselationParms(&parms);

  {
    float dist,max_e=0.0;
    int n,p,vn0,vn2;
    fprintf(stderr,"computing the maximum edge length...");
    for (n = 0 ; n < mris->nvertices ; n++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[n];
      VERTEX          const * const v  = &mris->vertices         [n];
      for (p = 0 ; p < vt->vnum ; p++) {
        VERTEX const * const vp = &mris->vertices[vt->v[p]];
        dist=SQR(vp->x - v->x)+SQR(vp->y - v->y)+SQR(vp->z - v->z);
        if (dist>max_e) max_e=dist;
      }
    }
    fprintf(stderr,"%f mm",sqrt(max_e));
    fprintf(stderr,"\nreversing orientation of faces...");
    for (n = 0 ; n < mris->nfaces ; n++) {
      vn0=mris->faces[n].v[0];
      vn2=mris->faces[n].v[2];
      /* vertex 0 becomes vertex 2 */
      { 
        VERTEX_TOPOLOGY* const v=&mris->vertices_topology[vn0];
        for (p = 0 ; p < v->num ; p++)
          if (v->f[p]==n)
            v->n[p]=2;
        mris->faces[n].v[2]=vn0;
      }
      /* vertex 2 becomes vertex 0 */
      { 
        VERTEX_TOPOLOGY* const v=&mris->vertices_topology[vn2];
        for (p = 0 ; p < v->num ; p++)
          if (v->f[p]==n)
            v->n[p]=0;
        mris->faces[n].v[0]=vn2;
      }
    }
  }

  mrisCheckVertexFaceTopology(mris);

  fprintf(stderr,"\nchecking orientation of surface...");
  MRISmarkOrientationChanges(mris);
  mris_corrected=MRISextractMainComponent(mris,0,1,0);

  MRISfree(&mris);

  fprintf(stderr,"\nwriting out surface...");
  //MRISaddCommandLine(mris_corrected, cmdline);
  strcpy(mris_corrected->fname, argv[1]);
  MRIScopyVolGeomFromMRI(mris_corrected, mri_orig) ;
  //if (mriConformed(mri_orig) == 0) {
  //  printf("input volume is not conformed - using useRealRAS=1\n") ;
  //  mris_corrected->useRealRAS = 1 ;
  //}   // (mr) maybe bad idea to assume this, e.g. in highres stream volume will not be 256 cube
  //  getVolGeom(mri, &mris_corrected->vg);
  MRISwrite(mris_corrected,argv[3]);
  fprintf(stderr,"done\n");

  MRIfree(&mri);
  MRIfree(&mri_orig);
  MRISfree(&mris_corrected);

  return 0;
}


