/*
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
#include <errno.h>

#include "macros.h"
#include "fio.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int spherical_coordinate(double x, double y, double z,double *pphi,
                                double *ptheta) ;
static int parameterization_coordinate(float x,
                                       float y,float z,int *pu,int *pv);
static void AnnotToParameterization(MRI_SURFACE *mris) ;
static void ParameterizationToAnnot(MRI_SURFACE *mris) ;
static void UpdateAnnotHist(void) ;
static void GetAnnotMode(void) ;
static int ReadAnnotFile(MRI_SURFACE *mris, char *fname) ;
static int WriteAnnotFile(MRI_SURFACE *mris, char *fname) ;
static int WriteAnnotFreqFile(MRI_SURFACE *mris,char *fname) ;
static int WriteAnnotHistFile(MRI_SURFACE *mris, char *fname) ;

const char *Progname ;

static int normalize_flag = 0 ;
static int condition_no = 0 ;
static int stat_flag = 0 ;
static char *output_surf_name = NULL ;
static float sigma = 0.0f ;

static int *AnnotLabel, *AnnotCount, **AnnotMapLabel, **AnnotMapCount;
static int ***AnnotHistLabel, ***AnnotHistCount;
static int **AnnotHistLabelVertex, **AnnotHistCountVertex, *AnnotHistNumVertex;
static int udim=256, vdim=128, maxlabels=25, maxvertices=500000, TotalCount;

int
main(int argc, char *argv[]) {
  char         **av, *in_fname, *out_fname, *surf_name, fname[200], *sdir,
  *hemi ;
  int          ac, nargs, i ;
  MRI_SURFACE  *mris ;

  int u,v,index;

  nargs = handleVersionOption(argc, argv, "mris_label_mode");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  sdir = getenv("SUBJECTS_DIR") ;
  if (!sdir)
    ErrorExit(ERROR_BADPARM, "%s: no SUBJECTS_DIR in envoronment.\n",Progname);
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 6)
    usage_exit() ;

  in_fname = argv[1] ;
  hemi = argv[2] ;
  surf_name = argv[3] ;
  out_fname = argv[argc-1] ;

  AnnotLabel = (int *)calloc(maxvertices,sizeof(int));
  AnnotCount = (int *)calloc(maxvertices,sizeof(int));
  AnnotHistNumVertex = (int *)calloc(maxvertices,sizeof(int));
  AnnotHistLabelVertex = (int **)calloc(maxvertices,sizeof(int *));
  AnnotHistCountVertex = (int **)calloc(maxvertices,sizeof(int *));
  AnnotMapLabel = (int **)calloc(udim,sizeof(int *));
  AnnotMapCount = (int **)calloc(udim,sizeof(int *));
  AnnotHistCount = (int ***)calloc(udim,sizeof(int **));
  AnnotHistLabel = (int ***)calloc(udim,sizeof(int **));
  for (u=0;u<udim;u++) {
    AnnotMapLabel[u] = (int *)calloc(vdim,sizeof(int));
    AnnotMapCount[u] = (int *)calloc(vdim,sizeof(int));
    AnnotHistCount[u] = (int **)calloc(udim,sizeof(int *));
    AnnotHistLabel[u] = (int **)calloc(udim,sizeof(int *));
    for (v=0;v<vdim;v++) {
      AnnotMapLabel[u][v] = -1;
      AnnotMapCount[u][v] = 0;
      AnnotHistCount[u][v] = (int *)calloc(maxlabels,sizeof(int));
      AnnotHistLabel[u][v] = (int *)calloc(maxlabels,sizeof(int));
      for (index=0;index<maxlabels;index++) {
        AnnotHistCount[u][v][index] = 0;
        AnnotHistLabel[u][v][index] = -1;
      }
    }
  }

  TotalCount = 0;
  for (i = 4 ; i < argc-1 ; i++) {
    TotalCount++;
    fprintf(stderr, "processing subject %s...\n", argv[i]) ;
    sprintf(fname, "%s/%s/surf/%s.%s", sdir, argv[i], hemi, surf_name) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    sprintf(fname, "%s/%s/label/%s%s.annot", sdir, argv[i], hemi, in_fname) ;
    ReadAnnotFile(mris, fname);
    AnnotToParameterization(mris);
    UpdateAnnotHist();
    if (i < argc-2)
      MRISfree(&mris) ;
  }
  GetAnnotMode() ;
  ParameterizationToAnnot(mris) ;
  sprintf(fname, "%s/%s/label/%s%s.annot", sdir, argv[i-1], hemi, out_fname) ;
  WriteAnnotFile(mris, fname) ;
  sprintf(fname, "%s/%s/label/%s%s.annot.freq", sdir,argv[i-1],hemi,out_fname);
  WriteAnnotFreqFile(mris, fname) ;
  sprintf(fname, "%s/%s/label/%s%s.annot.hist", sdir,argv[i-1],hemi,out_fname);
  WriteAnnotHistFile(mris, fname) ;

  MRISfree(&mris) ;
  exit(0) ;
  return(0) ;  /* for ansi */

}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option)) {
    case 'A':
      sigma = atof(argv[2]) ;
      fprintf(stderr, "blurring thickness measures with sigma=%2.3f\n",sigma);
      nargs = 1 ;
      break ;
    case 'O':
      output_surf_name = argv[2] ;
      nargs = 1 ;
      fprintf(stderr, "painting output onto subject %s.\n", output_surf_name);
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    case 'S':   /* write out stats */
      stat_flag = 1 ;
      condition_no = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "writing out summary statistics as condition %d\n",
              condition_no) ;
      break ;
    case 'N':
      normalize_flag = 1 ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <input curv. file> <hemi> <surface> <subject> ... "
          " <output curv file >\n", Progname) ;
  fprintf(stderr, "the output curvature file will be painted onto the last "
          "subject name specified\non the command line.\n"
          "if the -s flag is specified then the last parameter specifies\n"
          "the directory in which to write the statistical maps.\n"
          "if the -o flag is specified then it overrides the last subject\n"
          "as the output surface.\n") ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will add a template into an average surface.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr, "-s <cond #>     generate summary statistics and write\n"
          "                them into sigavg<cond #>-<hemi>.w and\n"
          "                sigvar<cond #>-<hemi>.w.\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}


static int
parameterization_coordinate(float x, float y, float z, int *pu, int *pv) {
  double phi, theta, uf, vf ;
  int    u, v ;

  spherical_coordinate(x, y, z, &phi, &theta) ;
  uf = udim * phi / PHI_MAX ;
  vf = vdim * theta / THETA_MAX ;
  u = nint(uf) ;
  v = nint(vf) ;
  if (u < 0)  /* enforce spherical topology  */
    u = -u ;
  if (u >= udim)
    u = udim - (u-udim+1) ;
  if (v < 0)  /* enforce spherical topology  */
    v += vdim ;
  if (v >= vdim)
    v -= vdim ;
  *pu = u ;
  *pv = v ;
  return(NO_ERROR) ;
}

static int
spherical_coordinate(double x, double y, double z,double *pphi,double *ptheta) {
  double r, d ;

  r = sqrt(x*x + y*y + z*z) ;
  d = r*r-z*z ;
  if (d < 0.0)
    d = 0.0 ;

  *pphi = atan2(sqrt(d), z) ;
  *ptheta = atan2(y/r, x/r) ;
  if (*ptheta < 0.0f)
    *ptheta = 2 * M_PI + *ptheta ;  /* make it 0 --> 2*PI */
  return(NO_ERROR) ;
}

static void
AnnotToParameterization(MRI_SURFACE *mris) {
  int u,v,du,dv,vertex_index,k,nundef;

  for (u=0;u<udim;u++)
    for (v=0;v<vdim;v++)
      AnnotMapLabel[u][v] = -1;
  vertex_index = mris->nvertices;
  for (k=0;k<vertex_index;k++) {
    parameterization_coordinate(mris->vertices[k].x,mris->vertices[k].y,
                                mris->vertices[k].z,&u,&v);
    AnnotMapLabel[u][v] = AnnotLabel[k];
  }
  nundef = 1;
  while (nundef>0) {
    nundef = 0;
    for (u=0;u<udim;u++)
      for (v=0;v<vdim;v++)
        if (AnnotMapLabel[u][v]==-1) /* undefined AnnotLabel */
        {
          nundef++;
          for (du= -1;du<=1;du++)
            for (dv= -1;dv<=1;dv++)
              if ((u+du>=0)&&(u+du<udim)&&(v+dv>=0)&&(v+dv<vdim)&&
                  (AnnotMapLabel[u+du][v+dv]!=-1))
                AnnotMapLabel[u][v] = AnnotMapLabel[u+du][v+dv];
        }
    printf("AnnotToParameterization: nundef=%d\n",nundef);
  }
}

static void
ParameterizationToAnnot(MRI_SURFACE *mris) {
  int u,v,vertex_index,k,i,num,*list;

  vertex_index = mris->nvertices;
  for (k=0;k<vertex_index;k++) {
    parameterization_coordinate(mris->vertices[k].x,mris->vertices[k].y,
                                mris->vertices[k].z,&u,&v);
    AnnotLabel[k] = AnnotMapLabel[u][v];
    AnnotCount[k] = AnnotMapCount[u][v];
    list = AnnotHistLabel[u][v];
    for (i=0;i<maxlabels&&list[i]!=-1;i++) ;
    num = AnnotHistNumVertex[k] = i;
    AnnotHistLabelVertex[k] = (int *)calloc(num,sizeof(int));
    AnnotHistCountVertex[k] = (int *)calloc(num,sizeof(int));
    for (i=0;i<maxlabels&&list[i]!=-1;i++) {
      AnnotHistLabelVertex[k][i] = AnnotHistLabel[u][v][i];
      AnnotHistCountVertex[k][i] = AnnotHistCount[u][v][i];
    }
  }
}

static void
UpdateAnnotHist(void) {
  int u,v,label,*list,i;

  for (u=0;u<udim;u++)
    for (v=0;v<vdim;v++) {
      label = AnnotMapLabel[u][v];
      list = AnnotHistLabel[u][v];
      for (i=0;i<maxlabels&&list[i]!=-1&&list[i]!=label;i++);
      if (list[i]==label) AnnotHistCount[u][v][i]++;
      else {
        AnnotHistLabel[u][v][i] = label;
        AnnotHistCount[u][v][i] = 1;
      }
    }
}

static void
GetAnnotMode(void) {
  int u,v,*list,i,imaxcnt;

  for (u=0;u<udim;u++)
    for (v=0;v<vdim;v++) {
      list = AnnotHistCount[u][v];
      imaxcnt = 0;
      for (i=0;i<maxlabels&&list[i]>0;i++)
        if (list[i]>list[imaxcnt]) imaxcnt=i;
      AnnotMapLabel[u][v] = AnnotHistLabel[u][v][imaxcnt];
      AnnotMapCount[u][v] = AnnotHistCount[u][v][imaxcnt];
    }
}

static int
ReadAnnotFile(MRI_SURFACE *mris, char *fname) {
  int i,j,k,num,vertex_index;
  FILE *fp;

  vertex_index = mris->nvertices;
  fp = fopen(fname,"r");
  if (fp==NULL)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE, "ReadAnnotFile(%s): file not found",fname));
  for (k=0;k<vertex_index;k++)
    AnnotLabel[k]=0;
  num = freadInt(fp);
  printf("num=%d\n",num);
  for (j=0;j<num;j++) {
    k = freadInt(fp);
    i = freadFloat(fp);
    if (k>=vertex_index||k<0)
      printf("vertex index out of range: %d i=%d\n",k,i);
    else
      AnnotLabel[k] = i;
  }
  fclose(fp);
  return(NO_ERROR) ;
}

static int
WriteAnnotFile(MRI_SURFACE *mris, char *fname) {
  int k,num,i,vertex_index;
  FILE *fp;

  vertex_index = mris->nvertices;
  fp = fopen(fname,"wb");
  if (fp==NULL)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "WriteAnnotFile(%s) - could not open file (%d)",
                 fname, errno)) ;
  for (k=0,num=0;k<vertex_index;k++) if (AnnotLabel[k]!=0) num++;
  printf("num = %d\n",num);
  fwriteFloat(num,fp);
  for (k=0;k<vertex_index;k++) {
    if (AnnotLabel[k]!=0) {
      fwriteInt(k,fp);
      i = AnnotLabel[k];
      fwriteInt(i,fp);
    }
  }
  fclose(fp);
  printf("file %s written\n",fname);
  return(NO_ERROR) ;
}

static int
WriteAnnotFreqFile(MRI_SURFACE *mris,char *fname) {
  int k,num,vertex_index;
  float f;
  FILE *fp;

  vertex_index = mris->nvertices;
  fp = fopen(fname,"wb");
  if (fp==NULL)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "WriteAnnotFreqFile(%s): can't create file (%d)",
                 fname, errno)) ;
  for (k=0,num=0;k<vertex_index;k++) if (AnnotLabel[k]!=0) num++;
  printf("num = %d\n",num);
  fwriteInt(num,fp);
  for (k=0;k<vertex_index;k++) {
    if (AnnotLabel[k]!=0) {
      fwriteInt(k,fp);
      f = ((float)AnnotCount[k])/TotalCount;
      fwriteFloat(f,fp);
    }
  }
  fclose(fp);
  printf("file %s written\n",fname);
  return(NO_ERROR) ;
}

static int
WriteAnnotHistFile(MRI_SURFACE *mris, char *fname) {
  int k,num,vertex_index,i;
  FILE *fp;

  vertex_index = mris->nvertices;
  fp = fopen(fname,"wb");
  if (fp==NULL)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "WriteAnnotHistFile(%s) - could not open file (%d)",
                 fname, errno)) ;
  for (k=0,num=0;k<vertex_index;k++) if (AnnotHistNumVertex[k]>0) num++;
  printf("num = %d\n",num);
  fwriteInt(num,fp);
  for (k=0;k<vertex_index;k++) {
    if (AnnotHistNumVertex[k]>0) {
      fwriteInt(k,fp);
      fwriteInt(AnnotHistNumVertex[k],fp);
      for (i=0;i<AnnotHistNumVertex[k];i++)
        fwriteInt(AnnotHistLabelVertex[k][i],fp);
      for (i=0;i<AnnotHistNumVertex[k];i++)
        fwriteInt(AnnotHistCountVertex[k][i],fp);
    }
  }
  fclose(fp);
  printf("file %s written\n",fname);
  return(NO_ERROR) ;
}
