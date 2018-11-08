#define COMPILING_MRISURF_TOPOLOGY_FRIEND_CHECKED
/**
 * @file  mri_multispectral_segment.c
 * @brief segment tissue classes based on T1/PD volumes
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Florent Segonne
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:23 $
 *    $Revision: 1.7 $
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
#include <ctype.h>

#include "mri.h"
#include "diag.h"
#include "macros.h"
#include "error.h"
#include "cma.h"
#include "histo.h"
#include "mri.h"

#define DEBUG_MODE 1
#define OUTPUT_SURFACES 0

#define FILL_OUTSIDE   1
#define FILL_INSIDE    -1

#define MAXVERTICES 25000
#define MAXFACES 25000

static double field_strength = 3 ;

#define GM_MAX_T1_3T  3000
#define WM_MIN_T1_3T   800

#define GM_MAX_T1_1p5T  1200
#define WM_MIN_T1_1p5T   500

static long PDthresh = 0 ;
static double  GM_max_T1 = GM_MAX_T1_3T ; // used to be 1200
static double  WM_min_T1 = WM_MIN_T1_3T  ;
typedef struct ss_vertex_type_ {
  float x,y,z;
  float nx,ny,nz;
  float xf,yf,zf;
  float nxf,nyf,nzf;
  float ox,oy,oz;
  float mx,my,mz;
  float nc,onc,snc;
  int vnum;
  int v[10];
  int fnum;
  int f[10];
}
ss_vertex_type;

float direction[26][3];

ss_vertex_type vertex[MAXVERTICES];
int face[MAXFACES][3];
int nvertices,nfaces;

double xCOG,yCOG,zCOG,rad_Brain;
double txCOG,tyCOG,tzCOG;

int width,height,depth;

int WM_INTENSITY,WM_VARIANCE,WM_HALF_MAX,WM_HALF_MIN,WM_MAX,WM_MIN;
int CSF_intensity,CSF_HALF_MAX,CSF_MAX,CSF_MIN;
int GM_MIN, GM_INTENSITY,TRANSITION_INTENSITY;

int skull_type=0;
int surf_out=0;

/*new variables*/
int SKULL_PD;
int PEAK1;
int MRI_correction=0;

MRI *mri_T1,*mri_PD,*mri_Err,*mri_dst,*mri_test,*mri_CSF;
MRI *mri_orig = NULL ;
int mriT1=0, mriPD=0, mriErr=0,mriCSF=0,mriOut=0,mriSURF=0;
const char *Progname;
char *T1_fname, *out_fname,*Err_fname,*PD_fname,*tmp_fname,
*CSF_fname,surf_fname[512];

#if 0
static void shrink_Brain(void);
static float rtanh(float x);
#endif
static int Reading(void);

static int get_option(int argc, char *argv[]) ;

static void Error(char *string);

static void lisse(unsigned long *tab,int m);
static void analyse_curve(unsigned long *tab,int,int);
static void PDcurve(MRI* mri_PD);
static void brain_params(void);

static void read_geometry_init(void);
static void write_geometry_init(char *fname) ;
static void read_geometry(char *fname);

#if 0
static void normal_vector(float *nx,float *ny,float *nz,int f);
static void writevertices(char *fname);
static void writefaces(char *fname);
#endif

static void find_normal(float nx,float ny, float nz,float* n1,float *n2);
static void init_direction(void);
static void mean(float tab[4][9],float *moy);
static void normal_face(int f,float *n);
static void compute_normals(void);

static void shrink_Inner_Skull(void);
static void shrink_Outer_Skull(void);
static void shrink_Outer_Skin(void);

static void peel_brain(MRI *mri, float h,int type, short val);

static void write_image(MRI *mri);
static void write_surface(char *fname);
static void GenerateMRI(void);
static void label_voxel(void);
static void intensity_correction(void);


unsigned long Compteur[2000][5000];

/*Error routine*/

static void Error(char *string) {
  fprintf(stderr, "\nError %s\n",string) ;
  exit(1) ;
}


static void brain_params(void) {
  int i,j,k;
  unsigned long n;
  double x,y,z;
  unsigned long masse;
  int threshold=PEAK1*90/100;
  double mean,var,fmasse,fval;
  int a,b,c,d=3,nb=1,nbvox;
  int val;


  txCOG = tyCOG = tzCOG = 0;
  fmasse=0;
  for (k=0;k<depth;k++)
    for (j=0;j<height;j++)
      for (i=0;i<width;i++) {
        txCOG+= i*MRIFvox(mri_dst,i,j,k);
        tyCOG+= j*MRIFvox(mri_dst,i,j,k);
        tzCOG+= k*MRIFvox(mri_dst,i,j,k);

        fmasse+=MRIFvox(mri_dst,i,j,k);
      }

  txCOG/=fmasse;
  tyCOG/=fmasse;
  tzCOG/=fmasse;

  fprintf(stderr,"\n      first estimation of the COG coord: x=%d,y=%d, z=%d",(int)txCOG,(int)tyCOG,(int)tzCOG);


  nbvox=((2*nb+1)*SQR(2*nb+1));

  x=y=z=0;
  n=0;
  fmasse=0;
  for (k=nb*d;k<depth-nb*d;k++)
    for (j=nb*d;j<height-nb*d;j++)
      for (i=nb*d;i<width-nb*d;i++)
        if (MRIFvox(mri_dst,i,j,k)>threshold) {
          mean=0;
          var=0;
          for (c=-nb;c<nb+1;c++)
            for (b=-nb;b<nb+1;b++)
              for (a=-nb;a<nb+1;a++) {
                val=MRIFvox(mri_dst,i+a*d,j+b*d,k+c*d);
                mean+=val;
                var+=SQR(val);
              }
          mean/=nbvox;
          var=var/nbvox-SQR(mean);
          fval=1/(1+SQR(var));
          fmasse+=fval;
          x+=i*fval;
          y+=j*fval;
          z+=k*fval;
          n++;
        }

  xCOG=x/fmasse;
  yCOG=y/fmasse;
  zCOG=z/fmasse;

  masse=0;
  for (k=0;k<depth;k++)
    for (j=0;j<height;j++)
      for (i=0;i<width;i++)
        if (MRIFvox(mri_dst,i,j,k)>threshold) {
          masse+=MRIFvox(mri_dst,i,j,k);
          rad_Brain+=(SQR(i-xCOG)+SQR(j-yCOG)+SQR(k-zCOG))*MRIFvox(mri_dst,i,j,k);
        }

  rad_Brain=sqrt(rad_Brain/masse);

  fprintf(stderr,"\n      second estimation of the COG coord: x=%d,y=%d, z=%d, r=%d",(int)xCOG,(int)yCOG,(int)zCOG,(int)rad_Brain);

  txCOG=xCOG-txCOG;
  tyCOG=yCOG-tyCOG;
  tzCOG=zCOG-tzCOG;

  fprintf(stderr,"\n      direction of the top of the brain: x=%d,y=%d, z=%d",(int)txCOG,(int)tyCOG,(int)tzCOG);

}


static void
normal_face(int f,float *n) {
  float v1[3],v2[3];

  v1[0] = vertex[face[f][0]].x-vertex[face[f][1]].x;
  v1[1] = vertex[face[f][0]].y-vertex[face[f][1]].y;
  v1[2] = vertex[face[f][0]].z-vertex[face[f][1]].z;
  v2[0] = vertex[face[f][2]].x-vertex[face[f][1]].x;
  v2[1] = vertex[face[f][2]].y-vertex[face[f][1]].y;
  v2[2] = vertex[face[f][2]].z-vertex[face[f][1]].z;
  n[0] = v1[1]*v2[2]-v1[2]*v2[1];
  n[1] = v1[2]*v2[0]-v1[0]*v2[2];
  n[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

static void
compute_normals(void) {
  int j,k;
  ss_vertex_type *v;
  float n[3],nt[3],d;

  for (k=0;k<nvertices;k++) {
    v = &vertex[k];
    n[0] = n[1] = n[2] = 0;
    for (j=0;j<v->fnum;j++) {
      normal_face(v->f[j],nt);
      n[0] += nt[0];
      n[1] += nt[1];
      n[2] += nt[2];
    }
    d = sqrt(SQR(n[0])+SQR(n[1])+SQR(n[2]));
    v->nx = n[0]/d;
    v->ny = n[1]/d;
    v->nz = n[2]/d;
  }
}

static void read_geometry_init(void) {
  int i,j,k,n,last,next,skiplast,skipnext;
  char fname[512], *mri_dir;
  FILE *fp;

  mri_dir = getenv("FREESURFER_HOME");
  sprintf(fname,"%s/lib/bem/ic5.tri",mri_dir);

  fp = fopen(fname,"r");
  if (fp==NULL) {
    fprintf(stderr,"\n\nmri_strip_skull: ### cannot open file %s\n",fname);
    return;
  }


  fscanf(fp,"%d",&nvertices);
  for (k=0;k<nvertices;k++) {
    fscanf(fp,"%*d %f %f %f",
           &vertex[k].x,&vertex[k].y,&vertex[k].z);
    vertex[k].mx = vertex[k].my = vertex[k].mz = 0;
    vertex[k].vnum = 0;
    vertex[k].fnum = 0;  /* marty */
    vertex[k].nc = 0;
    vertex[k].snc = 0;
  }
  fscanf(fp,"%d",&nfaces);
  for (k=0;k<nfaces;k++) {
    fscanf(fp,"%*d");
    for (n=0;n<3;n++) {
      fscanf(fp,"%d",&face[k][n]);
      face[k][n]--;
    }
  }

  fclose(fp);
  fprintf(stderr,"mri_strip_skull: triangle file %s read\n",fname);
  fprintf(stderr,"nvertices=%d, nfaces=%d\n",nvertices,nfaces);

  for (k=0;k<nfaces;k++) {
    for (i=0;i<3;i++) {
      vertex[face[k][i]].f[vertex[face[k][i]].fnum++] = k;
      last = (i>0)?i-1:2;
      next = (i<2)?i+1:0;
      skiplast = skipnext = FALSE;
      for (j=0;j<vertex[face[k][i]].vnum;j++) {
        if (vertex[face[k][i]].v[j]==face[k][last])
          skiplast = TRUE;
        if (vertex[face[k][i]].v[j]==face[k][next])
          skipnext = TRUE;
      }
      if (!skiplast)
        vertex[face[k][i]].v[vertex[face[k][i]].vnum++]=face[k][last];
      if (!skipnext)
        vertex[face[k][i]].v[vertex[face[k][i]].vnum++]=face[k][next];
    }
  }
  compute_normals();
}


static void
write_geometry_init(char *fname) {
  FILE *fp;
  int  k,n;

  fp = fopen(fname,"w");
  if (fp==NULL) {
    printf("### File %s not found\n",fname);
    return;
  }
  fprintf(fp,"%5d\n",nvertices);
  for (k=0;k<nvertices;k++) {
    fprintf(fp,"%5d%10.4f%10.4f%10.4f\n",k+1,vertex[k].x,vertex[k].y,vertex[k].z);
  }
  fprintf(fp,"%5d\n",nfaces);
  for (k=0;k<nfaces;k++) {
    fprintf(fp,"%5d",k+1);
    for (n=0;n<3;n++) {
      fprintf(fp,"%5d",face[k][n]+1);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  printf("triangle file %s written\n",fname);
}

static void read_geometry(char *fname) {
  int i,j,k,n,last,next,skiplast,skipnext;
  FILE *fp;



  fp = fopen(fname,"r");
  if (fp==NULL) {
    fprintf(stderr,"\n\nmri_strip_skull: ### cannot open file %s\n",fname);
    return;
  }


  fscanf(fp,"%d",&nvertices);
  for (k=0;k<nvertices;k++) {
    fscanf(fp,"%f %f %f",
           &vertex[k].x,&vertex[k].y,&vertex[k].z);
    vertex[k].mx = vertex[k].my = vertex[k].mz = 0;
    vertex[k].vnum = 0;
    vertex[k].fnum = 0;  /* marty */
    vertex[k].nc = 0;
    vertex[k].snc = 0;
  }
  fscanf(fp,"%d",&nfaces);
  for (k=0;k<nfaces;k++) {
    for (n=0;n<3;n++) {
      fscanf(fp,"%d",&face[k][n]);
      /*face[k][n];*/
    }
  }

  fclose(fp);
  /*    fprintf(stderr,"mri_strip_skull: triangle file %s read\n",fname);
     fprintf(stderr,"nvertices=%d, nfaces=%d\n",nvertices,nfaces);*/

  for (k=0;k<nfaces;k++) {
    for (i=0;i<3;i++) {
      vertex[face[k][i]].f[vertex[face[k][i]].fnum++] = k;
      last = (i>0)?i-1:2;
      next = (i<2)?i+1:0;
      skiplast = skipnext = FALSE;
      for (j=0;j<vertex[face[k][i]].vnum;j++) {
        if (vertex[face[k][i]].v[j]==face[k][last])
          skiplast = TRUE;
        if (vertex[face[k][i]].v[j]==face[k][next])
          skipnext = TRUE;
      }
      if (!skiplast)
        vertex[face[k][i]].v[vertex[face[k][i]].vnum++]=face[k][last];
      if (!skipnext)
        vertex[face[k][i]].v[vertex[face[k][i]].vnum++]=face[k][next];
    }
  }
  compute_normals();
}

static void
init_surf_to_image(float rx, float ry, float rz) {
  int k;
  for (k=0;k<nvertices;k++) {
    vertex[k].x = rx*vertex[k].x+xCOG;
    vertex[k].y = ry*vertex[k].y+yCOG;
    vertex[k].z = rz*vertex[k].z+zCOG;
  }
}

#if 0
static float rtanh(float x) {
  return (x<0.0)?0.0:tanh(x);
}
#endif



/* Find 2 normals to a  vector nx, ny, nz */
static void find_normal(float nx,float ny, float nz,float* n1,float *n2) {
  float ps,ps_buff;
  int p,k;

  k=0;
  ps=10;
  for (p=0;p<26;p++) {
    ps_buff=direction[p][0]*nx+direction[p][1]*ny+direction[p][2]*nz;
    if (ps_buff<0)
      ps_buff=-ps_buff;
    if (ps_buff<ps) {
      ps=ps_buff;
      k=p;
    }
  }
  n1[0]=direction[k][0];
  n1[1]=direction[k][1];
  n1[2]=direction[k][2];

  n2[0]=ny*n1[2]-nz*n1[1];
  n2[1]=nz*n1[0]-nx*n1[2];
  n2[2]=nx*n1[1]-ny*n1[0];

  ps=sqrt(SQR(n2[0])+SQR(n2[1])+SQR(n2[2]));
  n2[0]/=ps;
  n2[1]/=ps;
  n2[2]/=ps;

}

static void init_direction(void) {
  int i,j,k,p=0;
  float norm;

  for (i=-1;i<2;i++)
    for (j=-1;j<2;j++)
      for (k=-1;k<2;k++)
        if (i || j || k) {
          norm=sqrt(SQR(i)+SQR(j)+SQR(k));
          direction[p][0]=i/norm;
          direction[p][1]=j/norm;
          direction[p++][2]=k/norm;
        }
}


/**************main shrinking-expanding routines**************************/
/*************************************************************************/

static void
shrink_Inner_Skull(void) {
  float x,y,z,sx,sy,sz,sd,sxn,syn,szn,sxt,syt,szt,nc;
  float force,force1,force2,force3;
  float d,dx,dy,dz,nx,ny,nz,avgnc;
  ss_vertex_type *v;
  int iter,k,m,n;
  int ninside=30,noutside=10;
  int it,jt,kt,h,niter;
  float r,F,E,rmin=5,rmax=20.;
  float decay=0.8,update=0.9;

  float outsamp[50],insamp[50];
  float fsteepness=0.5,outmax=0,fzero=65*PEAK1/100,valt,force0;

  int MRIflag=1,int_smooth=1;

  double lm,d10m[3],d10,f1m,f2m,dm,dbuff;
  float ***dist;

  float cout,cout_prec,coutbuff,varbuff,mean_sd[10],mean_dist[10];


  dist = (float ***) malloc( nvertices*sizeof(float**) );

  for ( it = 0; it < nvertices; it++ ) {
    dist[it] = (float**) malloc( 4*sizeof(float*) );
    for ( jt = 0; jt < 4; jt++ ) {
      dist[it][jt] = (float*) calloc( 3, sizeof(float));
    }
  }


  E=(1/rmin+1/rmax)/2;
  F=6/(1/rmin-1/rmax);

  for (k=0;k<nvertices;k++)
    for (m=0;m<4;m++)
      for (n=0;n<3;n++)
        dist[k][m][n]=0;

  for (n=0;n<10;n++) {
    mean_sd[n]=0;
    mean_dist[n]=0;
  }

  niter =int_smooth;
  force = 0.0f ;

  cout_prec = 0;

  for (iter=0;niter;iter++) {
    lm = d10 = f1m = f2m = dm = 0;
    for (k=0;k<nvertices;k++) {
      v = &vertex[k];
      v->ox = v->x;
      v->oy = v->y;
      v->oz = v->z;
      v->onc=v->nc;
    }

    for (k=0;k<nvertices;k++) {
      v = &vertex[k];
      x = v->ox;
      y = v->oy;
      z = v->oz;
      nx = v->nx;
      ny = v->ny;
      nz = v->nz;
      sx=sy=sz=sd=0;
      n=0;


      for (m=0;m<v->vnum;m++) {
        sx += dx = vertex[v->v[m]].ox - x;
        sy += dy = vertex[v->v[m]].oy - y;
        sz += dz = vertex[v->v[m]].oz - z;
        sd += sqrt(dx*dx+dy*dy+dz*dz);
        n++;
      }
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
      sd = sd/n;

      lm+=sd;

      nc = sx*nx+sy*ny+sz*nz;

      sxn = nc*nx;
      syn = nc*ny;
      szn = nc*nz;
      sxt=sx-sxn;
      syt=sy-syn;
      szt=sz-szn;

      v->nc=nc;

      force1=0;
      if (nc) {
        r= (nc>0) ? nc : -nc;
        r=SQR(sd)/(2*r);
        force1=(1+tanh(F*(1/r-E)))/2;
      } else
        Error("normal pbm !");

      f1m+=force1;

      if (MRIflag) {


        for (h= -noutside;h<ninside;h++) {
          kt = (int)(z-nz*h+0.5);
          it = (int)(x-nx*h+0.5);
          jt = (int)(y-ny*h+0.5);
          if (kt<0||kt>=depth||it<0||it>=width||jt<0||jt>=height)
            valt = 0;
          else {
            valt = MRIFvox(mri_dst,it,jt,kt);
          }
          if (h<=0)
            outsamp[-h] = valt;
          if (h>=0)
            insamp[h] = valt;

        }
        outmax = 0;
        force2=-1;
        for (h=1;h<noutside;h++) {
          valt = outsamp[h];
          force2*=valt<fzero? 0:1;
        }
        force3 = 1;
        for (h=0;h<ninside;h++) {
          valt = insamp[h];
          force3 *= valt<fzero ? 0:1;
        }

        valt = insamp[0];
        force0 = tanh((valt-fzero)*fsteepness)*force3;

        force = 0.1*force2+1.0*(force3-0.1)+0.5*force0;
      } else
        force=0;


      if (!MRIflag) {
        avgnc = 0;
        for (m=0;m<v->vnum;m++)
          avgnc += vertex[v->v[m]].onc;
        avgnc /= v->vnum;
        force += tanh((nc-avgnc)*0.5);
      } else
        force += tanh(nc*0.1);

      f2m+=force;

      if ((d=sqrt(sxt*sxt+syt*syt+szt*szt))>1.0) {
        sxt /= d;
        syt /= d;
        szt /= d;
      }


      dx = sxt*0.5 + force1*sxn+v->nx*force;
      dy = syt*0.5 + force1*syn+v->ny*force;
      dz = szt*0.5 + force1*szn+v->nz*force;

      if (MRIflag) {
        dx = decay*v->mx+update*dx;
        dy = decay*v->my+update*dy;
        dz = decay*v->mz+update*dz;
      }

      if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0) {
        dx /= d;
        dy /= d;
        dz /= d;
      }

      v->mx = dx;
      v->my = dy;
      v->mz = dz;

      d=sqrt(dx*dx+dy*dy+dz*dz);

      dm+=d;

      dist[k][iter%4][0]=x;
      dist[k][iter%4][1]=y;
      dist[k][iter%4][2]=z;

      d10m[0] = d10m[1] = d10m[2] = 0;

      for (n=0;n<4;n++) {
        d10m[0]+=dist[k][n][0]/4;
        d10m[1]+=dist[k][n][1]/4;
        d10m[2]+=dist[k][n][2]/4;
      }

      dbuff=0;
      for (n=0;n<4;n++)
        dbuff+=SQR(dist[k][n][0]-d10m[0])+SQR(dist[k][n][1]-d10m[1])+SQR(dist[k][n][2]-d10m[2]);

      d10+=dbuff/4;

      v->x += dx;
      v->y += dy;
      v->z += dz;
    }

    lm /=nvertices;
    f1m /=nvertices;
    f2m /=nvertices;
    dm /=nvertices;
    d10 /=nvertices;

    mean_sd[iter%10]=lm;
    mean_dist[iter%10]=d10;

    coutbuff=0;
    for (n=0;n<10;n++)
      coutbuff+=mean_sd[n]/10;

    varbuff=0;
    for (n=0;n<10;n++)
      varbuff+=SQR(mean_sd[n]-coutbuff);


    cout=varbuff;

    coutbuff=0;
    for (n=0;n<10;n++)
      coutbuff+=mean_dist[n]/10;

    varbuff=0;
    for (n=0;n<10;n++)
      varbuff+=SQR(mean_dist[n]-coutbuff);

    cout+=10*varbuff;

    coutbuff=cout;

    cout=(cout_prec+cout)/2;

    cout_prec=coutbuff;

    if (MRIflag)
      compute_normals();

    if ((niter==int_smooth) && !(iter % 5))
      fprintf(stderr,
              "\n%d: lm=%5.3f,f1m=%5.3f,f2m=%5.3f,dm=%5.3f,d10m=%5.3f,c=%5.3f"
              ,iter,lm,f1m,f2m,dm,d10,100*cout);

    if (niter==int_smooth) {
      if (((iter>20)&&(10000*cout<5))||(iter>100)) {
        niter--;
        MRIflag=0;
      };
    } else
      niter--;
  }
  fprintf(stderr,"\n%d iterations",iter);

  compute_normals();
}


static void mean(float tab[4][9],float *moy) {
  int p;
  for (p=0;p<4;p++)
    moy[p]=(2*tab[p][4]+tab[p][1]+tab[p][3]+tab[p][5]+
            tab[p][7])/6;
}


/* to get the Outer Skull: the criterium to find is not so obvious...*/

static void shrink_Outer_Skull(void) {
  float x,y,z,sx,sy,sz,sd,sxn,syn,szn,sxt,syt,szt,nc;
  float force,force1;
  float f0,f1,f2,f3,f4,f5;

  float d,dx,dy,dz,nx,ny,nz;
  ss_vertex_type *v;
  int iter,k,m,n;
  float test_samp[4][9];
  int a,b;
  int it,jt,kt,h,niter;
  float r,F,E,rmin=7,rmax=30.;
  float decay=0.8,update=0.9;
  float fzero,fmax;
  float val;
  int MRIflag=1,int_smooth=1;

  double lm,d10m[3],d10,f1m,f2m,dm,dbuff;
  float ***dist;
  int nb_GM,nb_TR;
  float cout,pcout=0,coutbuff,varbuff,mean_sd[10],mean_dist[10];
  float n1[3],n2[3];

  dist = (float ***) malloc( nvertices*sizeof(float**) );

  for ( it = 0; it < nvertices; it++ ) {
    dist[it] = (float**) malloc( 4*sizeof(float*) );
    for ( jt = 0; jt < 4; jt++ ) {
      dist[it][jt] = (float*) calloc( 3, sizeof(float));
    }
  }

  E=(1/rmin+1/rmax)/2;
  F=6/(1/rmin-1/rmax);

  fzero=PEAK1/10;
  fmax=PEAK1/2;

  for (k=0;k<nvertices;k++)
    for (m=0;m<4;m++)
      for (n=0;n<3;n++)
        dist[k][m][n]=0;

  for (n=0;n<10;n++) {
    mean_sd[n]=0;
    mean_dist[n]=0;
  }

  niter =int_smooth;
  force = 0.0f ;
  pcout=0;

  for (iter=0;niter;iter++) {
    cout = lm = d10 = f1m = f2m = dm = 0;
    for (k=0;k<nvertices;k++) {
      v = &vertex[k];
      v->ox = v->x;
      v->oy = v->y;
      v->oz = v->z;
      v->onc=v->nc;
    }

    for (k=0;k<nvertices;k++) {
      v = &vertex[k];
      x = v->ox;
      y = v->oy;
      z = v->oz;
      nx = v->nx;
      ny = v->ny;
      nz = v->nz;
      sx=sy=sz=sd=0;
      n=0;
      for (m=0;m<v->vnum;m++) {
        sx += dx = vertex[v->v[m]].ox - x;
        sy += dy = vertex[v->v[m]].oy - y;
        sz += dz = vertex[v->v[m]].oz - z;
        sd += sqrt(dx*dx+dy*dy+dz*dz);
        n++;
      }
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
      sd = sd/n;

      lm+=sd;

      nc = sx*nx+sy*ny+sz*nz;

      sxn = nc*nx;
      syn = nc*ny;
      szn = nc*nz;
      sxt=sx-sxn;
      syt=sy-syn;
      szt=sz-szn;

      v->nc=nc;


      force1=0;
      if (nc) {
        r= (nc>0) ? nc : -nc;
        r=SQR(sd)/(2*r);
        force1=(1+tanh(F*(1/r-E)))/2;
      } else
        Error("pbm de normale!");

      f1m+=force1;

      f0 = f1 = f2 = f3 = f4 = f5 =0;


      find_normal(nx,ny,nz,n1,n2);
      for (h=0;h<4;h++)
        for (a=-1;a<2;a++)
          for (b=-1;b<2;b++) {
            kt = (int)(z-nz*h+n1[2]*a+n2[2]*b+0.5);
            it = (int)(x-nx*h+n1[0]*a+n2[0]*b+0.5);
            jt = (int)(y-ny*h+n1[1]*a+n2[1]*b+0.5);
            if ((kt<0||kt>=depth||it<0||it>=width||jt<0||jt>=height))
              val=0;
            else
              val=MRIFvox(mri_dst,it,jt,kt);

            test_samp[h][3*b+a+4] = val;
          }


      val=test_samp[0][4];

      if (!val)
        force=0.25;
      else {
        nb_GM=0;
        nb_TR=0;
        for (h=0;h<4;h++) {
          for (a=0;a<9;a++) {
            if (test_samp[h][a]<=fzero)
              nb_GM++;
            if (test_samp[h][a]>fmax)
              nb_TR++;
          }
        }

        if (nb_TR>=18)
          force=-0.25;
        else if (nb_GM>=15)
          force=0.5;
        else {
          if (nb_GM>9 && nb_TR<9)
            force=0.3;
          else
            force=-0.1;
        }
      }




      force += tanh(nc*0.1);

      f2m+=force;

      if ((d=sqrt(sxt*sxt+syt*syt+szt*szt))>1.0) {
        sxt /= d;
        syt /= d;
        szt /= d;
      }

      force1=force1;
      dx = sxt*0.5 + force1*sxn+v->nx*force;
      dy = syt*0.5 + force1*syn+v->ny*force;
      dz = szt*0.5 + force1*szn+v->nz*force;

      dx = decay*v->mx+update*dx;
      dy = decay*v->my+update*dy;
      dz = decay*v->mz+update*dz;


      if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0) {
        dx /= d;
        dy /= d;
        dz /= d;
      }

      v->mx = dx;
      v->my = dy;
      v->mz = dz;


      d=sqrt(dx*dx+dy*dy+dz*dz);

      dm+=d;

      dist[k][iter%4][0]=x;
      dist[k][iter%4][1]=y;
      dist[k][iter%4][2]=z;

      d10m[0] = d10m[1] = d10m[2] = 0;

      for (n=0;n<4;n++) {
        d10m[0]+=dist[k][n][0]/4;
        d10m[1]+=dist[k][n][1]/4;
        d10m[2]+=dist[k][n][2]/4;
      }

      dbuff=0;
      for (n=0;n<4;n++)
        dbuff+=SQR(dist[k][n][0]-d10m[0])+SQR(dist[k][n][1]-d10m[1])+
               SQR(dist[k][n][2]-d10m[2]);

      d10+=dbuff/4;

      v->x += dx;
      v->y += dy;
      v->z += dz;
    }

    lm /=nvertices;
    f1m /=nvertices;
    f2m /=nvertices;
    dm /=nvertices;
    d10 /=nvertices;

    mean_sd[iter%10]=lm;
    mean_dist[iter%10]=d10;

    coutbuff=0;
    for (n=0;n<10;n++)
      coutbuff+=mean_sd[n]/10;

    varbuff=0;
    for (n=0;n<10;n++)
      varbuff+=SQR(mean_sd[n]-coutbuff);

    cout=varbuff;

    coutbuff=0;
    for (n=0;n<10;n++)
      coutbuff+=mean_dist[n]/10;

    varbuff=0;
    for (n=0;n<10;n++)
      varbuff+=SQR(mean_dist[n]-coutbuff);

    cout+=10*varbuff;

    coutbuff=cout;

    cout=(cout+pcout)/2;

    pcout=coutbuff;

    if (MRIflag)
      compute_normals();

    if ((niter==int_smooth) && !(iter % 5))
      fprintf(stderr,
              "%d: lm=%5.3f,f1m=%5.3f,f2m=%5.3f,dm=%5.3f,d10m=%5.3f,c=%5.3f\n"
              ,iter,lm,f1m,f2m,dm,d10,100*cout);

    if (niter==int_smooth) {
      if (((iter>20)&&(10000*cout<1))||(iter>100))
        niter--;
    } else
      niter--;
  }
  fprintf(stderr,"%d iterations",iter);
  compute_normals();
}


/*to get the Outer Skin*/

static void shrink_Outer_Skin(void) {
  float x,y,z,sx,sy,sz,sd,sxn,syn,szn,sxt,syt,szt,nc;
  float force,force1;
  float f0,f1,f2,f3,f4,f5;

  float d,dx,dy,dz,nx,ny,nz;
  ss_vertex_type *v;
  int iter,k,m,n;
  float samp_mean[4];
  float test_samp[4][9];
  int a,b;
  int it,jt,kt,h,niter;
  float r,F,E,rmin=5,rmax=20.;
  float decay=0.8,update=0.9;
  float fzero;
  float val;
  int MRIflag=1,int_smooth=5;  /*smoothing for 5 iterations*/

  double lm,d10m[3],d10,f1m,f2m,dm,dbuff;
  float ***dist;
  int nb_GM,nb_TR,nb_GTM;
  float cout,pcout=0,coutbuff,varbuff,mean_sd[10],mean_dist[10];
  float n1[3],n2[3];

  dist = (float ***) malloc( nvertices*sizeof(float**) );

  for ( it = 0; it < nvertices; it++ ) {
    dist[it] = (float**) malloc( 4*sizeof(float*) );
    for ( jt = 0; jt < 4; jt++ ) {
      dist[it][jt] = (float*) calloc( 3, sizeof(float));
    }
  }

  E=(1/rmin+1/rmax)/2;
  F=6/(1/rmin-1/rmax);

  fzero=PEAK1/20;
  if (field_strength > 1.5)
    fzero=PEAK1/10;
  if (PDthresh > 0)
    fzero = PDthresh ; // BRF

  for (k=0;k<nvertices;k++)
    for (m=0;m<4;m++)
      for (n=0;n<3;n++)
        dist[k][m][n]=0;

  for (n=0;n<10;n++) {
    mean_sd[n]=0;
    mean_dist[n]=0;
  }

  niter =int_smooth;
  force = 0.0f ;
  pcout=0;

  for (iter=0;niter;iter++) {
    cout = lm = d10 = f1m = f2m = dm = 0;
    for (k=0;k<nvertices;k++) {
      v = &vertex[k];
      v->ox = v->x;
      v->oy = v->y;
      v->oz = v->z;
      v->onc=v->nc;
    }

    for (k=0;k<nvertices;k++) {
      v = &vertex[k];
      x = v->ox;
      y = v->oy;
      z = v->oz;
      nx = v->nx;
      ny = v->ny;
      nz = v->nz;
      sx=sy=sz=sd=0;
      n=0;
      for (m=0;m<v->vnum;m++) {
        sx += dx = vertex[v->v[m]].ox - x;
        sy += dy = vertex[v->v[m]].oy - y;
        sz += dz = vertex[v->v[m]].oz - z;
        sd += sqrt(dx*dx+dy*dy+dz*dz);
        n++;
      }
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
      sd = sd/n;

      lm+=sd;

      nc = sx*nx+sy*ny+sz*nz;

      sxn = nc*nx;
      syn = nc*ny;
      szn = nc*nz;
      sxt=sx-sxn;
      syt=sy-syn;
      szt=sz-szn;

      v->nc=nc;


      force1=0;
      if (nc) {
        r= (nc>0) ? nc : -nc;
        r=SQR(sd)/(2*r);
        force1=(1+tanh(F*(1/r-E)))/2;
      } else
        Error("pbm de normale!");

      f1m+=force1;

      f0 = f1 = f2 = f3 = f4 = f5 =0;



      /******************************/

      find_normal(nx,ny,nz,n1,n2);
      for (h=0;h<4;h++)
        for (a=-1;a<2;a++)
          for (b=-1;b<2;b++) {
            kt = (int)(z-nz*h+n1[2]*a+n2[2]*b+0.5);
            it = (int)(x-nx*h+n1[0]*a+n2[0]*b+0.5);
            jt = (int)(y-ny*h+n1[1]*a+n2[1]*b+0.5);
            if ((kt<0||kt>=depth||it<0||it>=width||jt<0||jt>=height))
              val=0;
            else
              val=MRIFvox(mri_dst,it,jt,kt);

            test_samp[h][3*b+a+4] = val;
          }

#if 0
      val=test_samp[0][4];

      if (!val)                  /*very simple criterium !!!!!!*/
        force=-0.25;
      else {
        nb_GM=0;
        for (h=0;h<3;h++)
          for (a=0;a<9;a++)
            if (test_samp[h][a]>=fzero)
              nb_GM++;
        if (nb_GM>=10)
          force=0.3;
        else
          force=-0.2;
      }

#else                     /*more sophisticated criterium*/
      /*but, the idea is the same*/

      if (!val)
        force=-0.25;
      else if (val<=fzero/2)
        force=-0.1;
      else {
        mean(test_samp,samp_mean);
        nb_GM=0;
        nb_TR=0;
        nb_GTM=0;
        for (h=0;h<3;h++) {
          if (samp_mean[h]>=fzero)
            nb_GM++;
          if (samp_mean[h]<fzero/2)
            nb_TR++;
        }

        if (nb_TR>=2)
          force=-0.2;
        else if (nb_GM==3 )
          force=0.5;
        else if (nb_GM==2 && samp_mean[0]>fzero)
          force=0.5;
        else if (nb_TR==0)
          force=0.3;
        else {
          nb_GM=0;
          nb_TR=0;
          for (h=0;h<3;h++) {
            for (a=0;a<9;a++) {
              if (test_samp[h][a]>=fzero)
                nb_GM++;
              else if (test_samp[h][a]<fzero/2)
                nb_TR++;
              else
                nb_GTM++;
            }
          }

          if (nb_TR>=18)
            force=-0.3;
          else if (nb_GM>=18)
            force=0.5;
          else if (nb_GM>=15)
            force=0.3;
          else {
            if (nb_GM>9 && nb_TR<9)
              force=0.5;
            else if (nb_GTM>30)
              force=0.1;
            else
              force=-0.0;
          }
        }
      }

      force += tanh(nc*0.1);
#endif

      f2m+=force;

      if ((d=sqrt(sxt*sxt+syt*syt+szt*szt))>1.0) {
        sxt /= d;
        syt /= d;
        szt /= d;
      }

      force1=force1;
      dx = sxt*0.5 + force1*sxn+v->nx*force;
      dy = syt*0.5 + force1*syn+v->ny*force;
      dz = szt*0.5 + force1*szn+v->nz*force;

      dx = decay*v->mx+update*dx;
      dy = decay*v->my+update*dy;
      dz = decay*v->mz+update*dz;

      if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0) {
        dx /= d;
        dy /= d;
        dz /= d;
      }

      v->mx = dx;
      v->my = dy;
      v->mz = dz;


      d=sqrt(dx*dx+dy*dy+dz*dz);

      dm+=d;

      dist[k][iter%4][0]=x;
      dist[k][iter%4][1]=y;
      dist[k][iter%4][2]=z;

      d10m[0] = d10m[1] = d10m[2] = 0;

      for (n=0;n<4;n++) {
        d10m[0]+=dist[k][n][0]/4;
        d10m[1]+=dist[k][n][1]/4;
        d10m[2]+=dist[k][n][2]/4;
      }

      dbuff=0;
      for (n=0;n<4;n++)
        dbuff+=SQR(dist[k][n][0]-d10m[0])+SQR(dist[k][n][1]-d10m[1])+
               SQR(dist[k][n][2]-d10m[2]);

      d10+=dbuff/4;

      v->x += dx;
      v->y += dy;
      v->z += dz;
    }

    lm /=nvertices;
    f1m /=nvertices;
    f2m /=nvertices;
    dm /=nvertices;
    d10 /=nvertices;

    mean_sd[iter%10]=lm;
    mean_dist[iter%10]=d10;

    coutbuff=0;
    for (n=0;n<10;n++)
      coutbuff+=mean_sd[n]/10;

    varbuff=0;
    for (n=0;n<10;n++)
      varbuff+=SQR(mean_sd[n]-coutbuff);

    cout=varbuff;

    coutbuff=0;
    for (n=0;n<10;n++)
      coutbuff+=mean_dist[n]/10;

    varbuff=0;
    for (n=0;n<10;n++)
      varbuff+=SQR(mean_dist[n]-coutbuff);

    cout+=10*varbuff;

    coutbuff=cout;

    cout=(cout+pcout)/2;

    pcout=coutbuff;

    if (MRIflag)
      compute_normals();

    if ((niter==int_smooth) && !(iter % 5))
      fprintf(stderr,
              "%d: lm=%5.3f,f1m=%5.3f,f2m=%5.3f,dm=%5.3f,d10m=%5.3f,c=%5.3f\n"
              ,iter,lm,f1m,f2m,dm,d10,100*cout);

    if (niter==int_smooth) {
      if (((iter>20)&&(10000*cout<5))||(iter>230))
        niter--;
    } else {
      niter--;     /*smoothing*/
      rmin=3.33;
      rmax=10;
      E=(1/rmin+1/rmax)/2;
      F=6/(1/rmin-1/rmax);
    }
  }
  fprintf(stderr,"%d iterations",iter);
  compute_normals();
}




/************** routine not used *********************/
/*to be implemented to get the right brain surface (without CSF...)*/

#if 0
static void
shrink_Brain(void) {
  float x,y,z,sx,sy,sz,sd,sxn,syn,szn,sxt,syt,szt,nc;
  float force,force1,force2,force3;
  float d,dx,dy,dz,nx,ny,nz,avgnc;
  ss_vertex_type *v;
  int iter,k,m,n;
  int ninside=30,noutside=10;
  int it,jt,kt,h,niter;
  float r,F,E,rmin=5,rmax=20.;
  float decay=0.8,update=0.9;

  float outsamp[50],insamp[50];
  float fsteepness=0.5,outmax=0,fzero=PEAK1/2,valt,force0;

  int MRIflag=1,int_smooth=1;

  double lm,d10m[3],d10,f1m,f2m,dm,dbuff;
  float ***dist;

  float cout,cout_prec,coutbuff,varbuff,mean_sd[10],mean_dist[10];

  float test_samp[4][9],n1[3],n2[3],samp_mean[4];
  int a,b;
  int nb_BR;


  dist = (float ***) malloc( nvertices*sizeof(float**) );

  for ( it = 0; it < nvertices; it++ ) {
    dist[it] = (float**) malloc( 4*sizeof(float*) );
    for ( jt = 0; jt < 4; jt++ ) {
      dist[it][jt] = (float*) calloc( 3, sizeof(float));
    }
  }

  init_direction();

  E=(1/rmin+1/rmax)/2;
  F=6/(1/rmin-1/rmax);

  for (k=0;k<nvertices;k++)
    for (m=0;m<4;m++)
      for (n=0;n<3;n++)
        dist[k][m][n]=0;

  for (n=0;n<10;n++) {
    mean_sd[n]=0;
    mean_dist[n]=0;
  }

  niter =int_smooth;
  force = 0.0f ;

  cout_prec = 0;

  for (iter=0;niter;iter++) {
    lm = d10 = f1m = f2m = dm = 0;
    for (k=0;k<nvertices;k++) {
      v = &vertex[k];
      v->ox = v->x;
      v->oy = v->y;
      v->oz = v->z;
      v->onc=v->nc;
    }

    for (k=0;k<nvertices;k++) {
      v = &vertex[k];
      x = v->ox;
      y = v->oy;
      z = v->oz;
      nx = v->nx;
      ny = v->ny;
      nz = v->nz;
      sx=sy=sz=sd=0;
      n=0;
      for (m=0;m<v->vnum;m++) {
        sx += dx = vertex[v->v[m]].ox - x;
        sy += dy = vertex[v->v[m]].oy - y;
        sz += dz = vertex[v->v[m]].oz - z;
        sd += sqrt(dx*dx+dy*dy+dz*dz);
        n++;
      }
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
      sd = sd/n;

      lm+=sd;

      nc = sx*nx+sy*ny+sz*nz;

      sxn = nc*nx;
      syn = nc*ny;
      szn = nc*nz;
      sxt=sx-sxn;
      syt=sy-syn;
      szt=sz-szn;

      v->nc=nc;

      force1=0;
      if (nc) {
        r= (nc>0) ? nc : -nc;
        r=SQR(sd)/(2*r);
        force1=(1+tanh(F*(1/r-E)))/2;
      } else
        Error("pbm de normale!");

      f1m+=force1;

      if (MRIflag) {


#if 1
        find_normal(nx,ny,nz,n1,n2);
        for (h=0;h<4;h++)
          for (a=-1;a<2;a++)
            for (b=-1;b<2;b++) {
              kt = (int)(z-nz*h+n1[2]*a+n2[2]*b+0.5);
              it = (int)(x-nx*h+n1[0]*a+n2[0]*b+0.5);
              jt = (int)(y-ny*h+n1[1]*a+n2[1]*b+0.5);
              if ((kt<0||kt>=depth||it<0||it>=width||jt<0||jt>=height))
                valt=0;
              else {
                if (MRIFvox(mri_T1,it,jt,kt)>4000)
                  valt=0;
                else if (MRIFvox(mri_T1,it,jt,kt)>2000)
                  valt=fzero-1;
                else
                  valt = MRIFvox(mri_PD,it,jt,kt);
              }
              test_samp[h][3*b+a+4] = valt;
            }

        valt=test_samp[0][4];

        if (valt==0)
          force=-0.1;
        else if (valt<fzero)
          force=-0.0;
        else
          force=0.2;

        /*        {
              mean(test_samp,samp_mean);
              {
                nb_BR=0;
                for (h=0;h<4;h++)
                {
                  if (samp_mean[h]>=fzero)
                    nb_BR++;
                }

                if (nb_BR<2)
                  force=-0.2;
                else if (nb_BR>=3)
                  force=0.7;
                else
                {
                  nb_BR=0;
                  for (h=0;h<4;h++)
                  {
                    for (a=0;a<9;a++)
                    {
                      if (test_samp[h][a]>=fzero)
                        nb_BR++;
                    }
                  }
                  if (nb_BR<15)
                    force=-0.3;
                  else if (nb_BR>=30)
                    force=0.2;
                  else if (nb_BR>=12)
                    force=0.2;
                }
              }

            }
        */

#endif
      } else
        force=0;


      if (!MRIflag) {
        avgnc = 0;
        for (m=0;m<v->vnum;m++)
          avgnc += vertex[v->v[m]].onc;
        avgnc /= v->vnum;
        force += tanh((nc-avgnc)*0.5);
      } else
        force += tanh(nc*0.1);

      f2m+=force;

      if ((d=sqrt(sxt*sxt+syt*syt+szt*szt))>1.0) {
        sxt /= d;
        syt /= d;
        szt /= d;
      }


      dx = sxt*0.5 + force1*sxn+v->nx*force;
      dy = syt*0.5 + force1*syn+v->ny*force;
      dz = szt*0.5 + force1*szn+v->nz*force;

      if (MRIflag) {
        dx = decay*v->mx+update*dx;
        dy = decay*v->my+update*dy;
        dz = decay*v->mz+update*dz;
      }

      if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0) {
        dx /= d;
        dy /= d;
        dz /= d;
      }

      v->mx = dx;
      v->my = dy;
      v->mz = dz;

      d=sqrt(dx*dx+dy*dy+dz*dz);

      dm+=d;

      dist[k][iter%4][0]=x;
      dist[k][iter%4][1]=y;
      dist[k][iter%4][2]=z;

      d10m[0] = d10m[1] = d10m[2] = 0;

      for (n=0;n<4;n++) {
        d10m[0]+=dist[k][n][0]/4;
        d10m[1]+=dist[k][n][1]/4;
        d10m[2]+=dist[k][n][2]/4;
      }

      dbuff=0;
      for (n=0;n<4;n++)
        dbuff+=SQR(dist[k][n][0]-d10m[0])+SQR(dist[k][n][1]-d10m[1])+SQR(dist[k][n][2]-d10m[2]);

      d10+=dbuff/4;

      v->x += dx;
      v->y += dy;
      v->z += dz;
    }

    lm /=nvertices;
    f1m /=nvertices;
    f2m /=nvertices;
    dm /=nvertices;
    d10 /=nvertices;

    mean_sd[iter%10]=lm;
    mean_dist[iter%10]=d10;

    coutbuff=0;
    for (n=0;n<10;n++)
      coutbuff+=mean_sd[n]/10;

    varbuff=0;
    for (n=0;n<10;n++)
      varbuff+=SQR(mean_sd[n]-coutbuff);


    cout=varbuff;

    coutbuff=0;
    for (n=0;n<10;n++)
      coutbuff+=mean_dist[n]/10;

    varbuff=0;
    for (n=0;n<10;n++)
      varbuff+=SQR(mean_dist[n]-coutbuff);

    cout+=10*varbuff;

    coutbuff=cout;

    cout=(cout_prec+cout)/2;

    cout_prec=coutbuff;

    if (MRIflag)
      compute_normals();

    if ((niter==int_smooth) && !(iter % 5))
      fprintf(stderr,
              "\n%d: lm=%5.3f,f1m=%5.3f,f2m=%5.3f,dm=%5.3f,d10m=%5.3f,c=%5.3f"
              ,iter,lm,f1m,f2m,dm,d10,100*cout);

    if (niter==int_smooth) {
      if (((iter>20)&&(10000*cout<5))||(iter>100)) {
        niter--;
        MRIflag=0;
      };
    } else
      niter--;
  }
  fprintf(stderr,"\n%d iterations",iter);

  compute_normals();
}
#endif


static void write_image(MRI *mri) {
  int i,j,imnr,k,u,v;
  float x0,y0,z0,x1,y1,z1,x2,y2,z2,d0,d1,d2,dmax;
  float px0,py0,pz0,px1,py1,pz1,px,py,pz;
  int numu,numv;


  for (k=0;k<nfaces;k++) {
    x0 = vertex[face[k][0]].x;
    y0 = vertex[face[k][0]].y;
    z0 = vertex[face[k][0]].z;
    x1 = vertex[face[k][1]].x;
    y1 = vertex[face[k][1]].y;
    z1 = vertex[face[k][1]].z;
    x2 = vertex[face[k][2]].x;
    y2 = vertex[face[k][2]].y;
    z2 = vertex[face[k][2]].z;
    d0 = sqrt(SQR(x1-x0)+SQR(y1-y0)+SQR(z1-z0));
    d1 = sqrt(SQR(x2-x1)+SQR(y2-y1)+SQR(z2-z1));
    d2 = sqrt(SQR(x0-x2)+SQR(y0-y2)+SQR(z0-z2));
    dmax = (d0>=d1&&d0>=d2)?d0:(d1>=d0&&d1>=d2)?d1:d2;
    numu = ceil(2*d0);
    numv = ceil(2*dmax);


    for (v=0;v<=numv;v++) {
      px0 = x0 + (x2-x0)*v/numv;
      py0 = y0 + (y2-y0)*v/numv;
      pz0 = z0 + (z2-z0)*v/numv;
      px1 = x1 + (x2-x1)*v/numv;
      py1 = y1 + (y2-y1)*v/numv;
      pz1 = z1 + (z2-z1)*v/numv;
      for (u=0;u<=numu;u++) {
        px = px0 + (px1-px0)*u/numu;
        py = py0 + (py1-py0)*u/numu;
        pz = pz0 + (pz1-pz0)*u/numu;

        i=(int)(px+0.5);
        j = (int)(py+0.5);
        imnr = (int)(pz+0.5);

        if (i>=0 && i<width && j>=0 && j<height && imnr>=0 && imnr<depth)
          MRIFvox(mri,i,j,imnr) = 2500;
      }
    }
  }
  for (k=0;k<nvertices;k++) {
    i=(int)(vertex[k].x+0.5);
    j = (int)(vertex[k].y+0.5);
    imnr = (int)(vertex[k].z+0.5);
    if (i>=0 && i<width && j>=0 && j<height && imnr>=0 && imnr<depth)
      MRIFvox(mri,i,j,imnr) = 1000;
  }
}



static int Reading(void) {
  fprintf(stderr,"\nreading...");

  mri_T1=MRIread(T1_fname);
  if (mri_T1 == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open %s", T1_fname) ;
  mri_PD=MRIread(PD_fname);
  if (mri_PD == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not open %s", PD_fname) ;
  if (mri_T1->type != MRI_FLOAT)
  {
    MRI *mri_tmp = MRIchangeType(mri_T1, MRI_FLOAT, 0, 1, 1) ;
    MRIfree(&mri_T1) ; mri_T1 = mri_tmp ;
  }
  if (mri_PD->type != MRI_FLOAT)
  {
    MRI *mri_tmp = MRIchangeType(mri_PD, MRI_FLOAT, 0, 1, 1) ;
    MRIfree(&mri_PD) ; mri_PD = mri_tmp ;
  }
  {
    MRI *mri_conformed = MRIalloc(256, 256, 256, MRI_FLOAT), *mri_tmp ; 

    mri_tmp = MRIresampleFill(mri_T1, mri_conformed, SAMPLE_TRILINEAR, 0) ;
    MRIfree(&mri_T1) ; mri_T1 = mri_tmp ;
    mri_tmp = MRIresampleFill(mri_PD, mri_conformed, SAMPLE_TRILINEAR, 0) ;
    mri_orig = mri_PD ; mri_PD = mri_tmp ;
    MRIfree(&mri_conformed) ;
  }

  if (mriErr) {
    mri_Err=MRIread(Err_fname);
    if (mri_Err==NULL) {
      fprintf(stderr,"\n      reading of the file %s impossible",Err_fname);
      mriErr=0;
    }
  }

  if (mri_T1==NULL)
    Error("\nT1=NULL");
  if (mri_PD==NULL)
    Error("\nPD=NULL");


  width=mri_T1->width;
  height=mri_T1->height;
  depth=mri_T1->depth;

  fprintf(stderr,"\ndone");
  fprintf(stderr,"\ncloning...");

  mri_dst = MRIclone(mri_T1, NULL) ;
  if (!mri_dst)
    Error("\n mri_dst NUll");

  mri_test=MRIclone(mri_T1, NULL) ;
  if (!mri_test)
    Error("\n mri_test NUll");

  mri_CSF=MRIclone(mri_T1, NULL) ;
  if (!mri_CSF)
    Error("\n mri_CSF NUll");

  fprintf(stderr,"\ndone");

  return 1;
}


/*routine to select an area inside or outside the surface,and
 fill the complement in the MR Image with the short val (usually =0)*/
static void peel_brain(MRI *mri, float height_offset,int type, short val) {
  int i,j,k,imnr;
  float x0,y0,z0,x1,y1,z1,x2,y2,z2,d0,d1,d2,dmax,u,v;
  float px,py,pz,px0,py0,pz0,px1,py1,pz1;
  int numu,numv,totalfilled,newfilled;

  MRI *mri_buff;

  mri_buff= MRIalloc(width, height, depth, MRI_UCHAR) ;

  for (k=0;k<nvertices;k++) {
    vertex[k].xf=vertex[k].x;
    vertex[k].yf=vertex[k].y;
    vertex[k].zf=vertex[k].z;
    vertex[k].nxf=vertex[k].nx;
    vertex[k].nyf=vertex[k].ny;
    vertex[k].nzf=vertex[k].nz;

    vertex[k].x+=height_offset*vertex[k].nx;
    vertex[k].y+=height_offset*vertex[k].ny;
    vertex[k].z+=height_offset*vertex[k].nz;
  }

  for (k=0;k<nfaces;k++) {
    x0 = vertex[face[k][0]].x;
    y0 = vertex[face[k][0]].y;
    z0 = vertex[face[k][0]].z;
    x1 = vertex[face[k][1]].x;
    y1 = vertex[face[k][1]].y;
    z1 = vertex[face[k][1]].z;
    x2 = vertex[face[k][2]].x;
    y2 = vertex[face[k][2]].y;
    z2 = vertex[face[k][2]].z;
    d0 = sqrt(SQR(x1-x0)+SQR(y1-y0)+SQR(z1-z0));
    d1 = sqrt(SQR(x2-x1)+SQR(y2-y1)+SQR(z2-z1));
    d2 = sqrt(SQR(x0-x2)+SQR(y0-y2)+SQR(z0-z2));
    dmax = (d0>=d1&&d0>=d2)?d0:(d1>=d0&&d1>=d2)?d1:d2;
    numu = ceil(2*d0);
    numv = ceil(2*dmax);

    for (v=0;v<=numv;v++) {
      px0 = x0 + (x2-x0)*v/numv;
      py0 = y0 + (y2-y0)*v/numv;
      pz0 = z0 + (z2-z0)*v/numv;
      px1 = x1 + (x2-x1)*v/numv;
      py1 = y1 + (y2-y1)*v/numv;
      pz1 = z1 + (z2-z1)*v/numv;
      for (u=0;u<=numu;u++) {
        px = px0 + (px1-px0)*u/numu;
        py = py0 + (py1-py0)*u/numu;
        pz = pz0 + (pz1-pz0)*u/numu;


        i=(int)(px+0.5);
        j = (int)(py+0.5);
        imnr = (int)(pz+0.5);

        if (i>=0 && i<width && j>=0 && j<height && imnr>=0 && imnr<depth)
          MRIvox(mri_buff,i,j,imnr) = 255;

      }
    }
  }


  MRIvox(mri_buff,1,1,1)= 64;
  totalfilled = newfilled = 1;
  while (newfilled>0) {
    newfilled = 0;
    for (k=1;k<depth-1;k++)
      for (j=1;j<height-1;j++)
        for (i=1;i<width-1;i++)
          if (MRIvox(mri_buff,i,j,k)==0)
            if (MRIvox(mri_buff,i,j,k-1)==64||MRIvox(mri_buff,i,j-1,k)==64||
                MRIvox(mri_buff,i-1,j,k)==64) {
              MRIvox(mri_buff,i,j,k)= 64;
              newfilled++;
            }

    for (k=depth-2;k>=1;k--)
      for (j=height-2;j>=1;j--)
        for (i=width-2;i>=1;i--)
          if (MRIvox(mri_buff,i,j,k)==0)
            if (MRIvox(mri_buff,i,j,k+1)==64||MRIvox(mri_buff,i,j+1,k)==64||
                MRIvox(mri_buff,i+1,j,k)==64) {
              MRIvox(mri_buff,i,j,k) = 64;
              newfilled++;
            }
    totalfilled += newfilled;
    /*    fprintf(stderr,"filled: new=%d, total=%d of %d\n",
          newfilled,totalfilled,(depth-2)*(width-2)*(height-2));*/
  }

  if (type==-1)
    for (k=1;k<depth-1;k++)
      for (j=1;j<height-1;j++)
        for (i=1;i<width-1;i++) {
          if (MRIvox(mri_buff,i,j,k)==64)
            MRIFvox(mri,i,j,k) = val;
        }
  else if (type==1)
    for (k=1;k<depth-1;k++)
      for (j=1;j<height-1;j++)
        for (i=1;i<width-1;i++) {
          if (MRIvox(mri_buff,i,j,k)!=64 && MRIvox(mri_buff,i,j,k)!=255)
            MRIFvox(mri,i,j,k) = val;
        }

  for (k=0;k<nvertices;k++) {
    vertex[k].x=vertex[k].xf;
    vertex[k].y=vertex[k].yf;
    vertex[k].z=vertex[k].zf;
    vertex[k].nx=vertex[k].nxf;
    vertex[k].ny=vertex[k].nyf;
    vertex[k].nz=vertex[k].nzf;
  }

  free(mri_buff);

  fprintf(stderr,"\n      mri_strip_skull: done peeling brain");
}







/* smooth the tab curve*/
static void lisse(unsigned long *tab,int m) {
  int k,n,p;
  unsigned long tmp[3];
  unsigned long buff;
  p=2;
  tmp[0]=0;
  tmp[1]=0;
  for (k=2;k<m-2;k++) {
    buff=0;
    for (n=-2;n<3;n++)
      buff+=tab[k+n]/5;
    tmp[p]=buff;
    p=(p+1)%3;
    tab[k-2]=tmp[p];
  }
  p=(p+1)%3;
  tab[m-4]=tmp[p];
  p=(p+1)%3;
  tab[m-3]=tmp[p];
  tab[m-2]=0;
  tab[m-1]=0;
}


#define HISTO_SIZE 300

static void analyse_curve(unsigned long *tab,int dim,int scale_factor) {
  int n,k,k1,k2;
  unsigned long max;
  unsigned long *tab_buff;
  float a,b,Sxy,Sx,Sy,Sxx;

  tab_buff=(unsigned long *)calloc(dim,sizeof(unsigned long));

  lisse(tab,dim);

#if DEBUG_MODE
  n=0;
  for (k=1;k<(HISTO_SIZE-1);k++)
    if (n<tab[k])
      n=tab[k];

  for (k=1;k<HISTO_SIZE;k++) {
    fprintf(stderr,"\n%3d ",k);
    for (k1=0;k1<(int)((double)tab[k]*65/n);k1++)
      fprintf(stderr,".");
  }
#endif

  max=0;
  for (k=0;k<dim;k++)
    if (max<tab[k])
      max=tab[k];


  k1=0;
  k2=dim;
  for (k=dim-1;k>0;k--) {
    if (!k1) {
      if (tab[k]<max/5)
        k2=k;
      else
        k1=k2;
    } else {
      if (tab[k]<3*max/4)
        k1=k;
      else break;
    }
  }




  n=k2-k1+1;
  Sxy = Sx = Sy = Sxx = 0;
  for (k=k1;k<=k2;k++) {
    Sxy+=(float)k*tab[k];
    Sx+=k;
    Sy+=tab[k];
    Sxx+=k*k;
  }

  if (FZERO((n*Sxx-Sx*Sx)))
    a = -1 ;
  else
    a=(n*Sxy-Sy*Sx)/(n*Sxx-Sx*Sx);
  b=-(a*Sx-Sy)/n;
  SKULL_PD=(int)((-b/a));

  max=0;
  for (k=SKULL_PD/2;k<dim;k++)
    if (tab[k]>max) {
      PEAK1=k;
      max=tab[k];
    }


  SKULL_PD*=scale_factor;
  fprintf(stderr," k1=%d k2=%d SKULL_PD=%d %d", k1, k2,SKULL_PD/scale_factor,PEAK1);

  PEAK1*=scale_factor;


}

static void PDcurve(MRI* mri) {
  int i,j,k,n,val;
  unsigned long countPD[HISTO_SIZE],nb_vox;
  int max,scale;

  for (k=0;k<HISTO_SIZE;k++)
    countPD[k]=0;

  max=0;
  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++) {
        val=MRIFvox(mri,i,j,k);
        if (val>max)
          max=val;
      }

  if (max > 2000)
    max = 2000 ;
#if 0
  {
    long total, count ;
    double max_pct = .99 ;
    HISTOGRAM *h ;

    h = MRIhistogram(mri, 300) ;
    HISTOplot(h, "h.plt") ;
    total = HISTOtotal(h) ;
    for (count = 0L, n = 0 ; n <h->nbins ; n++)
    {
      count += h->counts[n] ;
      if ((double)count / (double)total > max_pct)
        break ;
    }
    max = h->bins[n] ;
  }
#endif

  nb_vox=0;
  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++) {
        val=MRIFvox(mri,i,j,k);
        if (val) {
          n=val*(HISTO_SIZE-1)/max;
          if (n > (HISTO_SIZE-1))
            n = (HISTO_SIZE-1) ;
          countPD[n]++;
          nb_vox++;
        }
      }


  scale=(HISTO_SIZE-1);
  for (k=298;k>=0;k--) {
    countPD[k]+=countPD[k+1];
    if (countPD[k]*1000<nb_vox)
      scale=k;
    else
      break;
  }

  scale=MIN((HISTO_SIZE-1),scale*12/10);

  fprintf(stderr,"max=%d scale=%d  %d ",max,scale*10/15,max*scale/(HISTO_SIZE-1)/(HISTO_SIZE-1));

  for (k=0;k<HISTO_SIZE;k++)
    countPD[k]=0;


  scale=max*scale/(HISTO_SIZE-1)/(HISTO_SIZE-1);

  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++) {
        val=MRIFvox(mri,i,j,k);
        if (val) {
          n=MIN((HISTO_SIZE-1),val/scale);
          countPD[n]++;
        }
      }
  n=0;



  /* for(k=1;k<(HISTO_SIZE-1);k++)
   if (n<countPD[k])
    n=countPD[k];

  for(k=1;k<HISTO_SIZE;k++)
  {
   fprintf(stderr,"\n%3d ",k);
   for(i=0;i<(int)((double)countPD[k]*65/n);i++)
    fprintf(stderr,".");
    }*/


  analyse_curve(countPD,HISTO_SIZE,scale);

}



static void write_surface(char *fname) {
  int  k, n ;
  MRI_SURFACE *mris ;
  double      x, y, z ;
  FILE        *fout ;

  mris = MRISalloc(nvertices, nfaces) ;
  for (k=0;k<nvertices;k++)
  {
    MRISsurfaceRASFromVoxelCached(mris, mri_PD, vertex[k].x, vertex[k].y, vertex[k].z, &x, &y, &z) ;
    MRISsetXYZ(mris,k,x,y,z);
  }
  for (k=0;k<nfaces;k++)
    for (n = 0 ; n < 3 ; n++)
      mris->faces[k].v[n] = face[k][n] ;
  
  mrisCheckVertexFaceTopology(mris);
  
  mris->type = MRIS_TRIANGULAR_SURFACE ;
  MRISwrite(mris, fname) ;
  MRISfree(&mris) ;

  fout=fopen(fname,"w");
  if (fout==NULL) {
    printf(" ### File %s not found\n",fname);
    return;
  }

  fprintf(fout," %d \n",nvertices);
  for (k=0;k<nvertices;k++)
    fprintf(fout, " %f %f %f \n",vertex[k].x,vertex[k].y,vertex[k].z);
  fprintf(fout," %d \n",nfaces);
  for (k=0;k<nfaces;k++) {
    for (n=0;n<3;n++)
      fprintf(fout," %d",face[k][n]);
    fprintf(fout,"\n");
  }

  fclose(fout);
}


void savetab(char *fname) {
  int k,n;
  FILE *fout;

  fout=fopen(fname,"w");

  for (k=0;k<2000;k++) {
    for (n=0;n<5000;n++)
      fprintf(fout, " %ld",Compteur[k][n]);
    fprintf(fout, "\n");
  }

  fclose(fout);
}

static void GenerateMRI(void) {
  int i,j,k;

  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++) 
      {
        if (MRIFvox(mri_T1,i,j,k))
          MRIFvox(mri_dst,i,j,k)=(MRIFvox(mri_PD,i,j,k)*exp(-100/MRIFvox(mri_T1,i,j,k)));
        else
          MRIFvox(mri_dst,i,j,k)=0;

        MRIFvox(mri_test,i,j,k)=MRIFvox(mri_PD,i,j,k);
      }
}

static void intensity_correction(void) {
  int i,j,k,f,nf;

  nf=sqrt(SQR(txCOG)+SQR(tyCOG)+SQR(tzCOG));

  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++) {
        f=((i-xCOG)*txCOG+(j-yCOG)*tyCOG+(k-zCOG)*tzCOG)/nf;
        if (f >= 0) {
          MRIFvox(mri_dst,i,j,k)=((float)MRIFvox(mri_dst,i,j,k)*(1+f/rad_Brain*0.5));
        }
      }
}


static void label_voxel(void) {
  int i,j,k;
  char fname[512];

  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++)
        MRIFvox(mri_CSF,i,j,k)=5;

  strcpy(fname,"surface_Inner_Skull");
  read_geometry(fname);
  peel_brain(mri_CSF,-1,-1,9);

  strcpy(fname,"surface_Outer_Skull");
  read_geometry(fname);
  peel_brain(mri_CSF,1,-1,7);

  strcpy (fname,"surface_Outer_Skin");
  read_geometry(fname);
  peel_brain(mri_CSF,-3,-1,0); /*originally 0*/


  for (i=0;i<width;i++)
    for (j=0;j<height;j++)
      for (k=0;k<depth;k++)
      {
        if (i == Gx && j == Gy && k == Gz)
          DiagBreak() ;
        if (MRIFvox(mri_CSF,i,j,k)==5) {
          if (MRIFvox(mri_T1,i,j,k)>GM_max_T1 && MRIFvox(mri_PD,i,j,k)>80*PEAK1/100)
            MRIFvox(mri_CSF,i,j,k)=CSF; /*CSF*/
          
#if 0   /*could be usefull to determine white matter tissue- gray...
          However, the right threshold have to be found...*/
          
          else if (MRIFvox(mri_T1,i,j,k)>800 && MRIFvox(mri_T1,i,j,k)<1400)
            MRIFvox(mri_CSF,i,j,k)=800;  /*Gray matter*/
          
          else if (MRIFvox(mri_T1,i,j,k)>500 && MRIFvox(mri_T1,i,j,k)<700)
            MRIFvox(mri_T1,i,j,k)=600; /*white matter*/
          
          else
            MRIFvox(mri_T1,i,j,k)=500; /*ambiguous*/
#endif
        } 
        else if (MRIFvox(mri_CSF,i,j,k)==7) 
        {
          if (MRIFvox(mri_T1,i,j,k)<WM_min_T1)
            MRIFvox(mri_CSF,i,j,k)=Bone; /*SKULL*/
        }
      }
  
  strcpy(fname,"surface_Outer_Skin");
  read_geometry(fname);
  peel_brain(mri_CSF,-3,-1,7); /*originally 0*/

  strcpy(fname,"surface_Outer_Skin");
  read_geometry(fname);
  peel_brain(mri_CSF,0,-1,0);
}

static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */

  if (!stricmp(option, "T1")) {
    mriT1=1;
    T1_fname = argv[2];
    nargs = 1 ;
    fprintf(stderr," ") ;
  } else if (!stricmp(option, "PD")) {
    mriPD=1;
    PD_fname = argv[2];
    nargs = 1 ;
    fprintf(stderr," ") ;
  } else if (!stricmp(option, "PDthresh")) {
    PDthresh=atol(argv[2]);
    nargs = 1 ;
    fprintf(stderr,"using intracranial PD thresh %ld\n", PDthresh) ;
  } else if (!stricmp(option, "ERR")) {
    mriErr=1;
    Err_fname = argv[2];
    nargs = 1 ;
    fprintf(stderr," ") ;
  } else if (!stricmp(option, "3T")) {
    GM_max_T1 = GM_MAX_T1_3T ; // used to be 1200
    WM_min_T1 = WM_MIN_T1_3T  ;
    fprintf(stderr,"using 3T T1 thresholds ") ;
  } else if (!stricmp(option, "1.5T")) {
    field_strength = 1.5 ;
    GM_max_T1 = GM_MAX_T1_1p5T ;
    WM_min_T1 = WM_MIN_T1_1p5T  ;
    fprintf(stderr,"using 3T T1 thresholds ") ;
  } else if (!stricmp(option, "max_GM_T1")) {
    GM_max_T1 = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"using max GM T1 threshold = %2.0f\n", GM_max_T1) ;
  } else if (!stricmp(option, "min_WM_T1")) {
    WM_min_T1 = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr,"using min WM T1 threshold = %2.0f\n", WM_min_T1) ;
  } else if (!stricmp(option, "CSF")) {
    mriCSF=1;
    CSF_fname = argv[2];
    nargs = 1 ;
    fprintf(stderr," ") ;
  } else if (!stricmp(option, "OUT")) {
    mriOut=1;
    out_fname = argv[2];
    nargs = 1 ;
    fprintf(stderr," ") ;
  } else if (!stricmp(option, "correction")) {
    MRI_correction=1;
    fprintf(stderr," mode correction of the MR intensity: top of the head...") ;
  } else if (!stricmp(option, "SURF")) {
    mriSURF = 1;
    strcpy(surf_fname,argv[2]);
    nargs=1;
  } else switch (toupper(*option)) {
    default:
      printf("Mode:          unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

int main(int argc, char *argv[]) {
  int nargs;
  char fname[512];

  Progname=argv[0];


  /************* Command line****************/

  fprintf(stderr,"\n");

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }



  if (!(mriT1*mriPD)) {
    fprintf(stderr, "\nUsage: %s -[option] file -[option] file ...\n", Progname);
    printf(" specify -3T to use default 3T thresholds (GM < %d, WM > %d)\n", GM_MAX_T1_3T, WM_MIN_T1_3T) ;
    printf(" specify -1.5T to use default 1.5T thresholds (GM < %d, WM > %d)\n", GM_MAX_T1_1p5T, WM_MIN_T1_1p5T) ;
    printf(" specify -max_GM_T1 <thresh> to set the max GM T1 threshold\n") ;
    printf(" specify -min_WM_T1 <thresh> to set the min WM T1 threshold\n") ;
    exit(1);
  };


  /**************Program*********************/

  Reading(); /*load and clone the MRI structures: T1, PD, Err, etc. */

  init_direction(); /*initialize the variable direction used
                   in the different template processes.
                   not very easy to understand in a first reading*/

  GenerateMRI();   /*generate the mri_dst structure
                   used to find the surfaces:
                   mri_dst=mri_PD*exp(-100/mri_T1).
                   this approach is better since it "removes"
                   some of the errors done during the generation
                   of the T1 and PD maps */


  PDcurve(mri_dst);  /*generate a curve with the mri_dst MR image
                    and find some parameters useful for the next steps:*/
  brain_params();    /*for instance, the COG coordinates,
                    an approximate brain radius, the intensity of the
                    CSF proton density, the main peak (PEAK1) in the
                    previously construted curve (this peak is supposed
                    to correspond to the white-gray matter density,
                    which means about 70% of the CSF PD*/

  if (MRI_correction)
    intensity_correction();


  fprintf(stderr,"\n %d %d \n",SKULL_PD,PEAK1);


  /***********Inner Skull Template Deformation**********/
  /*Try to localize the brain, using the fact that its proton density
   is very high.... threshold about 80% of PEAK1*/
  /*save the surface*/

  /*load the tesselated surface*/
  read_geometry_init();

  /*initialize the surface with the brain parameters*/
  init_surf_to_image(0.7*rad_Brain,0.7*rad_Brain,0.7*rad_Brain);

  /*template deformation process: expanding surface*/
  shrink_Inner_Skull();

  strcpy(fname,"surface_Inner_Skull");
  fprintf(stderr,"\nsaving the file %s ... ",fname);
  write_surface(fname);
  fprintf(stderr,"done");

  peel_brain(mri_dst,2,1,0); /*the skull is assume to have, at least,
                        a 2-mm thickness*/

  /*write surfaces in the MRI sturcture mri_test*/
  write_image(mri_test);


  /*****************Outer Skin Template Deformation***************/

  /*load the tesselated surface*/
  read_geometry_init();

  /*initialize the surface with the brain parameters*/
  init_surf_to_image(2*rad_Brain,2*rad_Brain,2*rad_Brain);

  /*template deformation process: shrinking surface*/
  shrink_Outer_Skin();

  /*save the surface*/
  // strcpy(fname,"./");
  strcpy(fname,"surface_Outer_Skin");
  fprintf(stderr,"\nsaving the file %s ...",fname);
  write_surface(fname);
  fprintf(stderr,"done");

#if OUTPUT_SURFACES
  /* sprintf(surf_fname,"v");
  writevertices(surf_fname);
  sprintf(surf_fname,"f");
  writefaces(surf_fname);*/
#endif

  peel_brain(mri_dst,-1,-1,PEAK1); /*originally -3 mm */

  /*write surfaces in the MRI sturcture mri_test*/
  write_image(mri_test);


  /*******************Outer Skull Template Deformation*************/

  /*load the tesselated surface: load the Inner Skull Surface*/
  strcpy(fname,"surface_Inner_Skull");
  fprintf(stderr,"\nloading the file %s ...",fname);
  read_geometry(fname);
  fprintf(stderr,"done");


  /*template deformation process: expanding surface*/
  shrink_Outer_Skull();

  /*save the surface*/
  strcpy(fname,"surface_Outer_Skull");
  fprintf(stderr,"\nsaving the file %s ...",fname);
  write_surface(fname);
  fprintf(stderr,"done");

  /*write surfaces in the MRI sturcture mri_test*/
  write_image(mri_test);


  /******************Label Process******************************/
  /*simple routine to get the voxelised labels of different tissue*/

  fprintf(stderr,"\nlabelize the volume...");
  label_voxel();
  fprintf(stderr,"\ndone\n");

  if (mri_orig)  // reslice back to original voxel coords
  {
    MRI *mri_tmp ;
    if (mriCSF)
    {
      mri_tmp = MRIresampleFill(mri_CSF, mri_orig, SAMPLE_NEAREST, 0) ;
      MRIfree(&mri_CSF) ; mri_CSF = mri_tmp ;
    }
    mri_tmp = MRIresampleFill(mri_test, mri_orig, SAMPLE_NEAREST, 0) ;
    MRIfree(&mri_test) ; mri_test = mri_tmp ;

    MRIfree(&mri_orig) ;
  }
  if (mriCSF)
    MRIwrite(mri_CSF,CSF_fname);
  if (mriOut)
    MRIwrite(mri_test,out_fname);

  if (mriSURF) {
    strcpy(fname,surf_fname);
    strcat(fname,"surface_Inner_Skull");
    read_geometry(fname);
    write_geometry_init(fname);
    strcpy(fname,surf_fname);
    strcat(fname,"surface_Outer_Skull");
    write_surface(fname);
    write_geometry_init(fname);
    strcpy(fname,surf_fname);
    strcat(fname,"surface_Outer_Skin");
    write_surface(fname);
    write_geometry_init(fname);
  }

  return 0;
}


#if 0
static void
normal_vector(float *nx,float *ny,float *nz,int f) {
  float v1[3],v2[3],d;

  v1[0] = vertex[face[f][0]].x-vertex[face[f][1]].x;
  v1[1] = vertex[face[f][0]].y-vertex[face[f][1]].y;
  v1[2] = vertex[face[f][0]].z-vertex[face[f][1]].z;
  v2[0] = vertex[face[f][2]].x-vertex[face[f][1]].x;
  v2[1] = vertex[face[f][2]].y-vertex[face[f][1]].y;
  v2[2] = vertex[face[f][2]].z-vertex[face[f][1]].z;
  *nx = v1[1]*v2[2]-v1[2]*v2[1];
  *ny = v1[2]*v2[0]-v1[0]*v2[2];
  *nz = v1[0]*v2[1]-v1[1]*v2[0];

  d = sqrt((*nx)*(*nx)+(*ny)*(*ny)+(*nz)*(*nz));
  *nx /= d;
  *ny /= d;
  *nz /= d;
}
static void writevertices(char *fname) {
  int k;
  FILE *fout;
  fout=fopen(fname,"w");
  for (k=0;k<nvertices;k++)
    fprintf(fout, " %f %f %f \n",vertex[k].x,vertex[k].y,vertex[k].z);
  fclose(fout);
}

static void writefaces(char *fname) {
  int k,n;
  FILE *fout;
  fout=fopen(fname,"w");
  for (k=0;k<nfaces;k++) {
    for (n=0;n<3;n++)
      fprintf(fout," %d",face[k][n]);
    fprintf(fout,"\n");
  }
  fclose(fout);
}
#endif

