/**
 * @file  spherical_stats.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.8 $
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
#include <string.h>
#include <math.h>
#include <ctype.h>


#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"


static char vcid[] = "$Id: spherical_stats.c,v 1.8 2011/03/02 00:04:40 nicks Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static float scale = 1 ;
static int all=0;
static int global_stats_per_label=0;
static int plt=1;

static char subjects_dir[STRLEN] ;
static char *hemi;
static char *sphere_fname = "sphere.reg" ;
static char *manual_annotation = "parc.annot";
static char *auto_annotation = "raparc.annot";
static char *man_vs_auto = "man_vs_auto";
static char *distance_error = "dist_error";
static char *mean_fname = "mean_stat";
static char *var_fname ="var_stat";
//static char *prob_fname ="prob_stat";
//static char *orig_name = "smoothwm" ;
static char *alignment_fname = "alignment";
static char *template_fname=NULL;
static char* suffix=NULL;
static char* current_subject=NULL;


static int *collapses[50]; /* not freed !!! */
static int ncollapses[50];

static int labels[100];
static float global_nvertices[100];
static float mean_global_dist[100];
static float mean_global_Efract[100];
static float mean_global_median_25[100];
static float mean_global_median_50[100];
static float mean_global_median_75[100];
static float mean_max[100];
static float mean_area[100];
static float mean_fraction[100];
static FILE* label_file[100];

static int global_count[100];
static int nlabels;

static float Total_Brain_Area;

static int do_combine = 1;

HISTOGRAM *histos[100];
#define MAX_DIST 6.0


#define IMAGES_PER_SURFACE   3   /* mean, variance, and dof */
//#define SURFACES         sizeof(curvature_names) / sizeof(curvature_names[0])/
//#define PARAM_IMAGES         (IMAGES_PER_SURFACE*SURFACES)
#define NIMAGES 3
#define PARAM_IMAGES         (IMAGES_PER_SURFACE*NIMAGES) // 9 -> man_vs_auto / sse / correlation

#define MAX_NUMBER_OF_SUBJECTS 100

static void normalize(HISTO *histo) {
  int n ;
  float total=0.0,scale;
  for (n=0;n<histo->nbins;n++)
    total += histo->counts[n];

  scale = 1/(total*histo->bin_size);
  if (scale)
    for (n=0;n<histo->nbins;n++)
      histo->counts[n] *= scale;
}

static void printGlobalStats(MRIS *mris, char *fname) {
  int m;
  int r,g,b;
  COLOR_TABLE *ct;
  CTE* cte;
  int index;
  FILE *f;
  char histoname[500];

  ct = mris->ct;

  f=fopen(fname,"w+");
  if (!f) {
    fprintf(stderr,"could not open file\n");
    return;
  }

  for ( m = 0 ; m < nlabels ; m++) {
    CTABfindAnnotation(ct,labels[m],&index) ;
    if (index >= 0 || index < ct->nentries) {
      cte = ct->entries[index];
      fprintf(f,"Label %s \n",cte->name);
      sprintf(histoname,"%s/%s.%s.plt", subjects_dir,hemi,cte->name);
    } else {
      MRISAnnotToRGB(labels[m],r,g,b);
      fprintf(f,"Label %d - [ %d , %d , %d ] \n",labels[m],r,g,b);
      sprintf(histoname,"%s/%s.%d_%d_%d.plt", subjects_dir,hemi,r,g,b);
    }
    if (global_count[m]) {
      global_nvertices[m] /= (float)global_count[m];
      mean_global_dist[m] /= (float)global_count[m];
      mean_global_Efract[m] /= (float)global_count[m];
      mean_global_median_25[m] /= (float)global_count[m];
      mean_global_median_50[m] /= (float)global_count[m];
      mean_global_median_75[m] /= (float)global_count[m];
      mean_max[m] /= (float)global_count[m];
      mean_fraction[m] /= (float)global_count[m];
      mean_area[m] /= (float)global_count[m];
    }

    normalize(histos[m]);

    fprintf(f,"Vertices# :  %3.0f - Area : %2.3f%% \n Fraction*sqrt(area) : %2.3f - Fraction : %2.3f - Mean Distance Error : %f \n",global_nvertices[m],100.0*mean_area[m],mean_fraction[m],100.0*mean_global_Efract[m],mean_global_dist[m]);

    fprintf(f,"Median (25,50,75) = ( %f , %f , %f ) - Max : %f\n",mean_global_median_25[m],mean_global_median_50[m],mean_global_median_75[m],mean_max[m]);

    fprintf(f,"\n");
    if (plt)
      HISTOplot(histos[m],histoname);
  }
  fclose(f);

}

static float computeScaleFactor(MRIS *mris,char *fname) {
  int n;
  float area,brain_scale=1.0f;
  MRISsaveVertexPositions(mris,CANONICAL_VERTICES);

  MRISreadVertexPositions(mris,fname);
  MRIScomputeTriangleProperties(mris) ;
  area=0.0f;
  for ( n = 0 ; n < mris->nfaces; n++)
    area += mris->faces[n].area;

  if (!FZERO(area))
    brain_scale = sqrt(area/(4*PI*10000));

  Total_Brain_Area=area;

  fprintf(stderr,"Area = %f (scaling factor = %f )\n",area,brain_scale);

  MRISrestoreVertexPositions(mris,CANONICAL_VERTICES);

  return brain_scale;
}

static float findMedian(HISTO *histo, float med) {
  int n ;
  float total=0.0;
  med *= 0.01;

  for (n=0;n<histo->nbins;n++) {
    total += histo->bin_size*histo->counts[n];
    if (total>med) break;
  }

  return (histo->bin_size*n);
}
static void printStats(MRIS *mris, char *fname) {
  FILE *f;
  int m,n,nvertices;
  int r,g,b,bin;
  float mean_dist,Efract,larea;
  float median_25,median_50,median_75;
  COLOR_TABLE *ct;
  float maxv;
  CTE* cte;
  int index;
  HISTOGRAM *histo;
  VERTEX *v;

  histo=HISTOalloc(100);
  histo->bin_size=MAX_DIST/100.0;

  ct = mris->ct;

  f=fopen(fname,"w+");
  if (!f) {
    fprintf(stderr,"could not open file\n");
    return;
  }

  for ( m = 0 ; m < nlabels ; m++) {

    CTABfindAnnotation(ct,labels[m],&index) ;
    if (index >= 0 || index < ct->nentries) {
      cte = ct->entries[index];
      fprintf(f,"Label %s \n",cte->name);
    } else {
      MRISAnnotToRGB(labels[m],r,g,b);
      fprintf(f,"Label %d - [ %d , %d , %d ] \n",labels[m],r,g,b);
    }

    maxv=0.0f;
    Efract=0.0f;
    mean_dist=0.0f;
    nvertices=0;
    for ( n = 0 ; n < mris->nvertices ; n++) {
      v=&mris->vertices[n];
      if (v->ripflag) continue;

      if (v->val != labels[m]) continue;

      nvertices++;

      if (v->val2==1) { /* mislabeled point */
        mean_dist += v->d;
        if (v->d>maxv)
          maxv=v->d;
        Efract += 1.0f;
        bin=MIN(histo->nbins-1,MAX(0,(int)(v->d/histo->bin_size)));
        histo->counts[bin]++;
      }
    }
    if (nvertices==0)
      fprintf(stderr,"Could not find any common vertices for label %d\n",labels[m]);
    else
      Efract = Efract/(float)nvertices;
    if (Efract) mean_dist = mean_dist/(Efract*nvertices);

    larea=0;
    for (n=0;n<mris->nfaces;n++) {
      if (mris->vertices[mris->faces[n].v[0]].val!=labels[m]) continue;
      if (mris->vertices[mris->faces[n].v[1]].val!=labels[m]) continue;
      if (mris->vertices[mris->faces[n].v[2]].val!=labels[m]) continue;
      larea+=mris->faces[n].area;
    }

    larea /= Total_Brain_Area;

    global_nvertices[m] += (float)nvertices;
    mean_global_dist[m] += mean_dist;
    mean_global_Efract[m] += Efract;
    mean_max[m] += maxv;
    mean_area[m] += larea;
    mean_fraction[m] += 100.0*sqrt(larea)*Efract;
    global_count[m]++;

    normalize(histo);

    for (n=0;n<histo->nbins;n++)
      histos[m]->counts[n]+=histo->counts[n];

    median_25=findMedian(histo,25);
    median_50=findMedian(histo,50);
    median_75=findMedian(histo,75);

    mean_global_median_25[m] += median_25;
    mean_global_median_50[m] += median_50;
    mean_global_median_75[m] += median_75;

    fprintf(f,"Vertices# :  %d - Area : %2.3f%% \nFraction : %2.3f%% - Fraction*sqrt(Area) : %2.3f \n" ,
            nvertices,100.0*larea,100.0*Efract,100.0*Efract*sqrt(larea));
    fprintf(f,"Mean Distance Error : %f - Median (25,50,75) = ( %f , %f , %f ) - Max : %f\n",mean_dist,median_25,median_50,median_75,maxv);
    fprintf(f,"\n");

    if (global_stats_per_label) {
      fprintf(label_file[m],"SUBJECT %s\n",current_subject);
      fprintf(label_file[m],"           FRACTION  :       %2.3f\n",100.0*Efract);
      fprintf(label_file[m],"           MEAN DIST :       %f\n",mean_dist);
      fprintf(label_file[m],"Vertices# :  %d - Area : %2.3f%% \nFraction : %2.3f%% - Fraction*sqrt(Area) : %2.3f \n" ,
              nvertices,100.0*larea,100.0*Efract,100.0*Efract*sqrt(larea));
      fprintf(label_file[m],"Mean Distance Error : %f - Median (25,50,75) = ( %f , %f , %f ) - Max : %f\n",mean_dist,median_25,median_50,median_75,maxv);
      fprintf(label_file[m],"\n");
    }

  }


  HISTOfree(&histo);
  fclose(f);
}



static int checkCollapsing(MRIS *mris) {
  int n,p,q;
  VERTEX *v;

  for (p=0; p < 50 ;p++) {

    if (ncollapses[p] == 0) continue;

    for (n=0; n < mris->nvertices; n++) {
      v=&mris->vertices[n];
      if (v->ripflag) continue;

      for ( q = 0 ; q < ncollapses[p] ; q++)
        if (v->annotation == collapses[p][q])
          v->annotation = collapses[p][0];
    }

  }
  return NO_ERROR;
}

#ifndef TNORM
#define TNORM(a,b,c) ( sqrt( (a)*(a) + (b)*(b) + (c)*(c) ) )
#endif

static float sphericaldistance(float x1,float y1,float z1,float x2,float y2,float z2) {
  float o,coso,dot,n1,n2;
  float dist;

  dot=x1*x2+y1*y2+z1*z2;

  n1=TNORM(x1,y1,z1);
  n2=TNORM(x2,y2,z2);

  coso=dot/(n1*n2);

  o=acos(coso);

  dist=fabs((o*(n1+n2)/2.0f));

  if (isnan(dist)) {
    return TNORM(x1-x2,y1-y2,z1-z2);
  }
  //fprintf(stderr,"\n ISNAN n1=%f n2=%f - coso=%f o=%f ",n1,n2,coso,o);

  // return TNORM(x1-x2,y1-y2,z1-z2); //cartesian distances

  return dist;
}

static void computeDistances(MRIS *mris) {
  int n,m,p,count;
  int annotation;
  float dist,x,y,z;
  VERTEX *v,*vp;
  int nlistvertices,*vertexlist;

  static int first=1; /* first time compute list of labels */

  fprintf(stderr,"compute distances for each label...");

  if (first) {  /* zero labels so far */
    nlabels=0;
    memset(global_nvertices,0,100*sizeof(float));
    memset(mean_global_dist,0,100*sizeof(float));
    memset(mean_global_Efract,0,100*sizeof(float));
    memset(global_count,0,100*sizeof(int));
    memset(mean_global_median_25,0,100*sizeof(float));
    memset(mean_global_median_50,0,100*sizeof(float));
    memset(mean_global_median_75,0,100*sizeof(float));
    memset(mean_max,0,100*sizeof(float));
    memset(mean_area,0,100*sizeof(float));
    memset(mean_fraction,0,100*sizeof(float));
  }

  /* set fieldsign values to -1 */
  for (n = 0 ; n < mris->nvertices ;n++)
    mris->vertices[n].fieldsign=-1;

  for (n=0 ; n < mris->nvertices ; n++) {
    v=&mris->vertices[n];
    if (v->ripflag) continue;
    if (v->fieldsign>=0) continue;
    annotation=v->annotation;

    //fprintf(stderr,"\nANNOTATION %d [%f]",annotation,v->fieldsign);

    if (first) {
      histos[nlabels]=HISTOalloc(100) /* min = 0 - max = 10mm */;
      histos[nlabels]->bin_size=MAX_DIST/100.0; /* 0.1 mm */
      labels[nlabels++]=annotation;
    }

    //  fprintf(stderr,"label = %d [%d]\n",nlabels,annotation);
    /* compute distances for label 'annotation' */

    /* first locate border labels */
    for ( nlistvertices = m = 0 ; m <  mris->nvertices ; m++) {
      v=&mris->vertices[m];
      if (v->ripflag) continue;
      if (v->annotation!=annotation) continue; /* not a label of interest */
      //if(v->fieldsign>0) fprintf(stderr,"#");
      dist=0.0f;
      count=0;
      x=v->x;
      y=v->y;
      z=v->z;
      for ( p=0 ; p < v->vnum; p++) {
        vp=&mris->vertices[v->v[p]];
        if (vp->ripflag) continue;
        if (vp->annotation!=annotation) { /* border vertex */
          dist += sphericaldistance(x,y,z,vp->x,vp->y,vp->z);
          count++;
        }
      }
      if (count>0) {
        v->fieldsign=dist/(float)count;
        nlistvertices++;
      }
    }

    //  fprintf(stderr,"nlist = %d nlistvertices\n",nlistvertices);

    vertexlist=(int*)malloc(nlistvertices*sizeof(int));
    /* init the list */
    for ( nlistvertices = m = 0 ; m <  mris->nvertices ; m++) {
      v=&mris->vertices[m];
      if (v->ripflag) continue;
      if (v->annotation==annotation && v->fieldsign>=0) {
        //fprintf(stderr,"%d-",nlistvertices);
        vertexlist[nlistvertices++]=m;
      }
    }

    //  fprintf(stderr,"/");
    /* compute the inside distance to the border*/
    for ( m = 0 ; m <  mris->nvertices ; m++) {
      v=&mris->vertices[m];
      if (v->ripflag) continue;
      if (v->annotation!=annotation) continue; /* not a label of interest */
      if (v->fieldsign>=0) continue; /*border label */
      /* compute distance to border */
      x=v->x;
      y=v->y;
      z=v->z;
      dist=10000.0f;
      for ( p = 0 ; p < nlistvertices ; p++) {
        vp=&mris->vertices[vertexlist[p]];
        dist=MIN(dist,sphericaldistance(x,y,z,vp->x,vp->y,vp->z));
      }
      v->fieldsign=dist;
    }
    if (vertexlist) free(vertexlist);
  }
  fprintf(stderr,"done\n");

  if (global_stats_per_label && first) {
    COLOR_TABLE *ct;
    CTE* cte;
    int index;
    char fname[500];
    /* allocation */
    ct = mris->ct;
    for ( m = 0 ; m < nlabels ; m++) {
      CTABfindAnnotation(ct,labels[m],&index) ;
      if (index >= 0 || index < ct->nentries) {
        cte = ct->entries[index];
        if (suffix)
          sprintf(fname,"%s/stats/%s.%s.%s.txt", subjects_dir,hemi,suffix,cte->name);
        else
          sprintf(fname,"%s/stats/%s.%s.txt", subjects_dir,hemi,cte->name);

        label_file[m]=fopen(fname,"w+");
      } else {
        if (suffix)
          sprintf(fname,"%s/stats/%s.%s.%d.txt", subjects_dir,hemi,suffix,labels[m]);
        else
          sprintf(fname,"%s/stats/%s.%d.txt", subjects_dir,hemi,labels[m]);
        label_file[m]=fopen(fname,"w+");
      }
      if (label_file[m]==NULL) {
        fprintf(stderr,"could not open file %s\n",fname);
        global_stats_per_label=0;
        break;
      }
    }
  }

  first=0;
}

int main(int argc, char *argv[]) {
  char         **av, *subjects_fname[MAX_NUMBER_OF_SUBJECTS],fname[STRLEN],*cp,*subject_fname;
  int          ac, nargs ;
  float brain_scale;
  MRI_SURFACE  *mris , *mris_atlas;
  MRI_SP *mrisp,*mrisp_template, *mrisp_average=NULL;
  VERTEX *v;
  int n,m,nsubjects;
  // int          msec, minutes, seconds ;
  // struct timeb start;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: spherical_stats.c,v 1.8 2011/03/02 00:04:40 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  memset(ncollapses,0,50*sizeof(int));

  //  TimerStart(&start) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (!strlen(subjects_dir)) /* hasn't been set on command line */
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
                Progname);
    strcpy(subjects_dir, cp) ;
  }

  if (argc < 3)
    usage_exit() ;

  /* hemisphere information */
  hemi = argv[1];

  /* load the list of subjects */
  nsubjects=0;
  for ( n = 2 ; n < argc ; n++) {
    subjects_fname[nsubjects]=argv[n];
    nsubjects++;
  }

  fprintf(stderr,"%d subjects have been selected\n",nsubjects);

  fprintf(stderr, "creating new parameterization...\n") ;
  mrisp_template = MRISPalloc(scale, PARAM_IMAGES);

  if (template_fname) {
    fprintf(stderr,"loading template from file %s\n",template_fname);
    mrisp_average=MRISPread(template_fname);
    if (!mrisp_average) ErrorExit(ERROR_NOFILE, "%s: could not read template file %s",
                                    Progname, template_fname) ;
  }

  for (n=0 ; n < nsubjects ; n++) {
    fprintf(stderr,"\n\nPROCESSING SUBJECT %d out of %d subjects\n",n+1,nsubjects);

    subject_fname=subjects_fname[n];
    current_subject = subject_fname;
    sprintf(fname,"%s/%s/surf/%s.%s", subjects_dir,subject_fname,hemi,sphere_fname);
    fprintf(stderr, "reading surface from %s...\n", fname) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s for %s",
                Progname, fname, subject_fname) ;


    sprintf(fname,"%s/%s/surf/%s.white", subjects_dir,subject_fname,hemi);
    brain_scale=computeScaleFactor(mris,fname);

    sprintf(fname,"%s/%s/label/%s.%s", subjects_dir,subject_fname,hemi,manual_annotation);
    fprintf(stderr, "reading manually labeled file from %s...\n", fname) ;
    MRISreadAnnotation(mris,fname);

    checkCollapsing(mris);
    computeDistances(mris);

    /* save annotation of manual into val and fieldsign into fsmask */
    for (m=0;m<mris->nvertices;m++) {
      v=&mris->vertices[m];
      v->val=v->annotation;
      v->fsmask=v->fieldsign;
    }

    sprintf(fname,"%s/%s/label/%s.%s", subjects_dir,subject_fname,hemi,auto_annotation);
    fprintf(stderr, "reading manually labeled file from %s...\n", fname) ;
    MRISreadAnnotation(mris,fname);
    checkCollapsing(mris);
    computeDistances(mris);
    MRISclearCurvature(mris);

    /* generating maps */
    for (m=0;m<mris->nvertices;m++) {
      v=&mris->vertices[m];
      v->d=0;
      if (v->annotation==(int)v->val)
        v->curv=0;
      else
        v->curv=1;
      v->val2=v->curv; /* temporary saving into v->val2 */
      v->d=0.5*brain_scale*(v->fieldsign+v->fsmask); /* error */
    }

    /* save stats */
    if (suffix)
      sprintf(fname,"%s/%s/surf/%s.%s.stats.txt", subjects_dir,subject_fname,hemi,suffix);
    else
      sprintf(fname,"%s/%s/surf/%s.stats.txt", subjects_dir,subject_fname,hemi);
    fprintf(stderr, "writing out label stats in file %s...\n", fname) ;
    printStats(mris,fname);

    /* save curvature */
    if (suffix)
      sprintf(fname,"%s/%s/surf/%s.%s.%s", subjects_dir,subject_fname,hemi,suffix,man_vs_auto);
    else
      sprintf(fname,"%s/%s/surf/%s.%s", subjects_dir,subject_fname,hemi,man_vs_auto);
    fprintf(stderr, "writing out differences in file %s.%s for subject %s...\n", hemi,man_vs_auto,subject_fname) ;
    MRISwriteCurvature(mris,fname);

    /* combine differences into template */
    fprintf(stderr,"Combining parameterizations 0\n");

    if (do_combine) {
      mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
      MRISPcombine(mrisp, mrisp_template, 0) ;
      MRISPfree(&mrisp) ;
    }

    for (m=0;m<mris->nvertices;m++) {
      v=&mris->vertices[m];
      if (v->val2==1)
        v->curv=v->d;
    }

    /* save distance error */
    if (suffix)
      sprintf(fname,"%s/%s/surf/%s.%s.%s", subjects_dir,subject_fname,hemi,suffix,distance_error);
    else
      sprintf(fname,"%s/%s/surf/%s.%s", subjects_dir,subject_fname,hemi,distance_error);
    fprintf(stderr, "writing out distaces errors in file %s.%s for subject %s...\n", hemi,man_vs_auto,subject_fname) ;
    MRISwriteCurvature(mris,fname);

    if ( n == nsubjects-1 ) { /* last subject */
      if (suffix)
        sprintf(fname,"%s/stats/%s.%s.global_stats_%d.txt", subjects_dir,hemi,suffix,nsubjects);
      else
        sprintf(fname,"%s/stats/%s.global_stats_%d.txt", subjects_dir,hemi,nsubjects);
      fprintf(stderr, "writting global stats in %s...\n", fname) ;
      printGlobalStats(mris,fname);
    }

    /* combine differences into template */
    fprintf(stderr,"Combining parameterizations 3\n");

    if (do_combine) {
      mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
      MRISPcombine(mrisp, mrisp_template, 3) ;
      MRISPfree(&mrisp) ;
    }

    if (mrisp_average) {
      /* loading sulcal information */
      sprintf(fname,"%s/%s/surf/%s.sulc", subjects_dir,subject_fname,hemi);
      fprintf(stderr, "reading sulcal file from %s...\n", fname) ;
      MRISreadCurvature(mris,fname);
      MRISnormalizeCurvature(mris, NORM_MEAN) ;
      for (m=0;m<mris->nvertices;m++) {
        v=&mris->vertices[m];
        v->valbak=v->curv;
      }
      MRISfromParameterization(mrisp_average,mris,6);
      MRISnormalizeCurvature(mris, NORM_MEAN) ;

      for (m=0;m<mris->nvertices;m++) {
        v=&mris->vertices[m];
        v->curv = fabs(v->valbak-v->curv);
      }
      if (suffix)
        sprintf(fname,"%s/%s/surf/%s.%s.%s", subjects_dir,subject_fname,hemi,suffix,alignment_fname);
      else
        sprintf(fname,"%s/%s/surf/%s.%s", subjects_dir,subject_fname,hemi,alignment_fname);
      fprintf(stderr, "writing out sulcal differences in alignment into file %s.%s for subject %s...\n", hemi,alignment_fname,subject_fname) ;
      MRISwriteCurvature(mris,fname);

      fprintf(stderr,"Combining parameterizations 6 \n");
      if (do_combine) {
        mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
        MRISPcombine(mrisp, mrisp_template, 6) ;
        MRISPfree(&mrisp) ;
      }


#if 0
      /* correlation */
      for (m=0;m<mris->nvertices;m++) {
        v=&mris->vertices[m];
        v->curv = v->curv*v->val2;
      }
      mrisp = MRIStoParameterization(mris, NULL, scale, 0) ;
      MRISPcombine(mrisp, mrisp_template, 6) ;
      MRISPfree(&mrisp) ;
#endif
    }

    MRISfree(&mris) ;

  }
  if (global_stats_per_label) { /* closing files */
    for (m=0;m<nlabels;m++)
      fclose(label_file[m]);
  }

  fprintf(stderr,"\nMean and Variance Information...\n");

  if (all && do_combine) {
    for (n=0 ; n < nsubjects ; n++) {
      fprintf(stderr,"\n\nPROCESSING SUBJECT %d out of %d subjects\n",n+1,nsubjects);

      subject_fname=subjects_fname[n];
      current_subject = subject_fname;
      sprintf(fname,"%s/%s/surf/%s.%s", subjects_dir,subject_fname,hemi,sphere_fname);
      fprintf(stderr, "reading surface from %s...\n", fname) ;
      mris_atlas = MRISread(fname) ;

#if 0
      if (n==0) {
        sprintf(fname,"%s/%s/label/%s.%s", subjects_dir,subject_fname,hemi,manual_annotation);
        MRISreadAnnotation(mris_atlas,fname);

        if (suffix)
          sprintf(fname,"%s/stats/%s.%s.global_stats.txt", subjects_dir,hemi,suffix);
        else
          sprintf(fname,"%s/stats/%s.global_stats.txt", subjects_dir,hemi);
        fprintf(stderr, "writting global stats in %s...\n", fname) ;
        printGlobalStats(mris_atlas,fname);
      }
#endif



      fprintf(stderr,"Extracting mean information...\n");
      MRISclearCurvature(mris_atlas);
      MRISfromParameterization(mrisp_template, mris_atlas, 0) ;
      if (suffix)
        sprintf(fname,"%s/%s/surf/%s.%s.%s", subjects_dir,subject_fname,hemi,suffix,mean_fname);
      else
        sprintf(fname,"%s/%s/surf/%s.%s", subjects_dir,subject_fname,hemi,mean_fname);
      fprintf(stderr,"writting mean information into %s of subject %s \n",mean_fname,subject_fname);
      MRISwriteCurvature(mris_atlas,fname);

      fprintf(stderr,"Extracting var information...\n");
      MRISclearCurvature(mris_atlas);
      MRISfromParameterization(mrisp_template, mris_atlas, 1) ;
      if (suffix)
        sprintf(fname,"%s/%s/surf/%s.%s.%s", subjects_dir,subject_fname,hemi,suffix,var_fname);
      else
        sprintf(fname,"%s/%s/surf/%s.%s", subjects_dir,subject_fname,hemi,var_fname);
      fprintf(stderr,"writting mean information into %s of subject %s \n",var_fname,subject_fname);
      MRISwriteCurvature(mris_atlas,fname);

      fprintf(stderr,"Extracting distance error information...\n");
      MRISclearCurvature(mris_atlas);
      MRISfromParameterization(mrisp_template, mris_atlas, 3) ;
      if (suffix)
        sprintf(fname,"%s/%s/surf/%s.%s.mean_distance_error", subjects_dir,subject_fname,hemi,suffix);
      else
        sprintf(fname,"%s/%s/surf/%s.mean_distance_error", subjects_dir,subject_fname,hemi);
      fprintf(stderr,"writting mean distance information into %s of subject %s \n",fname,subject_fname);
      MRISwriteCurvature(mris_atlas,fname);


      if (mrisp_average) {
        MRISclearCurvature(mris_atlas);
        MRISfromParameterization(mrisp_template, mris_atlas, 6) ;
        if (suffix)
          sprintf(fname,"%s/%s/surf/%s.%s.error", subjects_dir,subject_fname,hemi,suffix);
        else
          sprintf(fname,"%s/%s/surf/%s.error", subjects_dir,subject_fname,hemi);
        fprintf(stderr,"writting mean information into %s\n",fname);
        MRISwriteCurvature(mris_atlas,fname);

#if 0
        MRISclearCurvature(mris_atlas);
        MRISfromParameterization(mrisp_template, mris_atlas, 6) ;
        if (suffix)
          sprintf(fname,"%s/%s/surf/%s.%s.corr", subjects_dir,subject_fname,hemi,suffix);
        else
          sprintf(fname,"%s/%s/surf/%s.corr", subjects_dir,subject_fname,hemi);
        fprintf(stderr,"writting mean information into %s \n",fname);
        MRISwriteCurvature(mris_atlas,fname);
#endif
      }
      MRISfree(&mris_atlas);
    }
  } else if (do_combine) {

    subject_fname=subjects_fname[0];
    sprintf(fname,"%s/%s/surf/%s.%s", subjects_dir,subject_fname,hemi,sphere_fname);
    fprintf(stderr, "reading surface from %s...\n", fname) ;
    mris_atlas = MRISread(fname) ;

#if 0
    sprintf(fname,"%s/%s/label/%s.%s", subjects_dir,subject_fname,hemi,manual_annotation);
    fprintf(stderr,"reading annotation file %s\n",fname);
    MRISreadAnnotation(mris_atlas,fname);
    if (suffix)
      sprintf(fname,"%s/stats/%s.%s.global_stats.txt", subjects_dir,hemi,suffix);
    else
      sprintf(fname,"%s/stats/%s.global_stats.txt", subjects_dir,hemi);
    fprintf(stderr, "writting global stats in %s...\n", fname) ;
    printGlobalStats(mris_atlas,fname);
#endif

    fprintf(stderr,"Extracting mean information...\n");
    MRISclearCurvature(mris_atlas);
    MRISfromParameterization(mrisp_template, mris_atlas, 0) ;
    if (suffix)
      sprintf(fname,"%s/%s/surf/%s.%s.%s", subjects_dir,subject_fname,hemi,suffix,mean_fname);
    else
      sprintf(fname,"%s/%s/surf/%s.%s", subjects_dir,subject_fname,hemi,mean_fname);
    fprintf(stderr,"writting mean information into %s of subject %s \n",mean_fname,subject_fname);
    MRISwriteCurvature(mris_atlas,fname);

    fprintf(stderr,"Extracting var information...\n");
    MRISclearCurvature(mris_atlas);
    MRISfromParameterization(mrisp_template, mris_atlas, 1) ;

    if (suffix)
      sprintf(fname,"%s/%s/surf/%s.%s.%s", subjects_dir,subject_fname,hemi,suffix,var_fname);
    else
      sprintf(fname,"%s/%s/surf/%s.%s", subjects_dir,subject_fname,hemi,var_fname);
    fprintf(stderr,"writting mean information into %s of subject %s \n",var_fname,subject_fname);
    MRISwriteCurvature(mris_atlas,fname);

    fprintf(stderr,"Extracting distance error information...\n");
    MRISclearCurvature(mris_atlas);
    MRISfromParameterization(mrisp_template, mris_atlas, 3) ;
    if (suffix)
      sprintf(fname,"%s/%s/surf/%s.%s.mean_distance_error", subjects_dir,subject_fname,hemi,suffix);
    else
      sprintf(fname,"%s/%s/surf/%s.mean_distance_error", subjects_dir,subject_fname,hemi);
    fprintf(stderr,"writting mean distance information into %s of subject %s \n",fname,subject_fname);
    MRISwriteCurvature(mris_atlas,fname);

    if (mrisp_average) {
      MRISclearCurvature(mris_atlas);
      MRISfromParameterization(mrisp_template, mris_atlas, 6) ;
      if (suffix)
        sprintf(fname,"%s/%s/surf/%s.%s.error", subjects_dir,subject_fname,hemi,suffix);
      else
        sprintf(fname,"%s/%s/surf/%s.error", subjects_dir,subject_fname,hemi);
      fprintf(stderr,"writting mean information into %s\n",fname);
      MRISwriteCurvature(mris_atlas,fname);

#if 0
      MRISclearCurvature(mris_atlas);
      MRISfromParameterization(mrisp_template, mris_atlas, 6) ;
      if (suffix)
        sprintf(fname,"%s/%s/surf/%s.%s.corr", subjects_dir,subject_fname,hemi,suffix);
      else
        sprintf(fname,"%s/%s/surf/%s.corr", subjects_dir,subject_fname,hemi);
      fprintf(stderr,"writting mean information into %s \n",fname);
      MRISwriteCurvature(mris_atlas,fname);
#endif
    }
    MRISfree(&mris_atlas);
  }
  if (mrisp_average) MRISPfree(&mrisp_average);
  MRISPfree(&mrisp_template);

  for (n=0;n<50;n++)
    if (ncollapses[n])
      free(collapses[n]);

  for (n=0; n<nlabels;n++)
    HISTOfree(&histos[n]);

  // msec = TimerStop(&start) ;
  //seconds = (int)((float)msec/1000.0f) ;
  //minutes = seconds / 60 ;
  //seconds = seconds % 60 ;
  //printf("spherical statistics took %d minutes and %d seconds.\n", minutes, seconds) ;
  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int    n,m,nargs = 0,r,g,b = 0 ;
  char   *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "SDIR")) {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  } else if (!stricmp(option, "all_stats")) {
    global_stats_per_label=1;
    fprintf(stderr,"generating global stats file per label\n");
  } else if (!stricmp(option, "do_not_combine")) {
    do_combine=0;
    fprintf(stderr,"do not generate maps\n");
  } else if (!stricmp(option, "all")) {
    all=1;
    fprintf(stderr,"projecting mean stats on all brains\n");
  } else if (!stricmp(option, "sphere")) {
    sphere_fname=argv[2];
    fprintf(stderr,"using sphere %s\n",sphere_fname);
    nargs=1;
  } else if (!stricmp(option, "suffix")) {
    suffix=argv[2];
    fprintf(stderr,"using suffix %s in all files\n",suffix);
    nargs=1;
  } else if (!stricmp(option, "auto_annot")) {
    auto_annotation=argv[2];
    fprintf(stderr,"using file %s to read automated annotation\n",auto_annotation);
    nargs=1;
  } else if (!stricmp(option, "plt")) {
    plt=1;
    fprintf(stderr,"generating histograms\n");
  } else if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "template")) {
    template_fname=argv[2];
    fprintf(stderr,"loading brain template from %s\n",template_fname);
    nargs=1;
  } else switch (toupper(*option)) {
    case 'C':
      for ( n = 0 ; n<50;n++)
        if (ncollapses[n]==0)
          break;
      ncollapses[n]=atoi(argv[2]);
      fprintf(stderr,"collapsing %d labels : ",ncollapses[n]);
      collapses[n]=malloc(ncollapses[n]*sizeof(int));
      for (m = 0; m < ncollapses[n] ; m++) {
        r=atoi(argv[3*m+3]);
        g=atoi(argv[3*m+4]);
        b=atoi(argv[3*m+5]);
        MRISRGBToAnnot(r,g,b,collapses[n][m]);
        fprintf(stderr,"%d ",collapses[n][m]);
      }
      fprintf(stderr,"\n");
      nargs = 1+3*ncollapses[n];
      break;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
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

static void print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <hemisphere> <subjects_1 subjects_2 ...>\n",
          Progname) ;
  fprintf(stderr,
          "\nThis program computes statistics from a set of manually and automatically labeld brains.\n"
          "The results are written out as rh.man_vs_auto files for each brain\n"
          "The first brain will be used to write the global statistics\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
}

static void
print_help(void) {
  print_usage() ;

  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

