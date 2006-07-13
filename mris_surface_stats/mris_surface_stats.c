//
// mris_surface_stats.c
// compute the stat-map of scalars defined on the surface 
// original author: Xiao Han
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: fischl $
// Revision Date  : $Date: 2006/07/13 16:29:58 $
// Revision       : $Revision: 1.2 $
//
////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "mri.h"
#include "mri_identify.h"
#include "mrisurf.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "matrix.h"
#include "transform.h"
#include "version.h"
#include "label.h"
#define DEBUG 0

#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
MRI *MyMRISsmoothMRI(MRIS *Surf, MRI *Src, int nSmoothSteps, MRI *Targ);

char *Progname ;

static char *stat_fname = NULL; /* filename for output stat */
static char *surf_fname = NULL; /* filename for underlying surface */

static char *mean_fname = NULL; /* filename for output stat */
static char *absmean_fname = NULL; /* filename for output stat */
static char *absstd_fname = NULL; /* filename for output stat */
static char *mask_fname = NULL; /* filename for subcortical mask */

static char *zscore_fname = NULL; /* filename for output stat */
static int first_group_size = 0;  /* number of subjects for first group */

char *srctypestring = "paint";
int srctype = MRI_VOLUME_TYPE_UNKNOWN;
char *trgtypestring = "paint";
int trgtype = MRI_VOLUME_TYPE_UNKNOWN;

static void usage_exit(int code) ;

int debugflag = 0;
int debugvtx = 0;


static int nSmoothSteps = 0;

#ifdef MAX_SURFACES
#undef MAX_SURFACES
#endif
#define MAX_SURFACES 200

static char vcid[] = "$Id: mris_surface_stats.c,v 1.2 2006/07/13 16:29:58 fischl Exp $";

int
main(int argc, char *argv[])
{
  char   **av, *in_fname;
  int    ac, nargs, i, index, total;
  MRIS    *mri_surf;
  MRI    *surfVal[MAX_SURFACES];
  MRI    *mriMean, *mriStd, *mriAbsMean, *mriAbsStd;
  MRI    *mriMean2, *mriStd2, *mriAbsMean2, *mriAbsStd2;
  MRI    *Zscore;
  double mean, std, scalar, absmean, absstd;
  LABEL *masklabel = NULL;
  VERTEX *v;

  int    msec, minutes, seconds, nsurfaces, nsurfaces_total ;
  struct timeb start ;
  
  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_surface_stats.c,v 1.2 2006/07/13 16:29:58 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;
  
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;
  
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  
  if (argc < 2)
    usage_exit(1) ;

  // printf("command line parsing finished\n");

  if(first_group_size > 0 && (zscore_fname == NULL)){
    printf("Pls use -zscore to specify output Zscore; o.w., nothing will be output. Exit. \n");
    exit(0);
  }
  
  if(argc > 2 && stat_fname == NULL){
    printf("Use -out_name option to specify output stat file name. Exit. \n");
    exit(1);
  }

  if(surf_fname == NULL){
    printf("Use -surf_name option to specify underlying surface file name. Exit. \n");
    exit(1);
  }
    
  mri_surf = MRISread(surf_fname) ;
  if (mri_surf == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s",
	      Progname, surf_fname) ;

  if(mask_fname){
    // printf("read in surface mask ...\n");
    masklabel = LabelRead(NULL, mask_fname);
    if(masklabel == NULL){
      ErrorExit(ERROR_NOFILE, "%s: could not read label file %s",
		Progname, mask_fname) ;
      exit(1);
    }

    for(index = 0; index < mri_surf->nvertices; index++){
      v = &mri_surf->vertices[index] ;
      v->ripflag = 1;
    }

    for(index = 0; index < masklabel->n_points; index++){
      v = &mri_surf->vertices[masklabel->lv[index].vno] ;
      v->ripflag = 0;
    }
    LabelFree(&masklabel) ;
  }
  
  /*** Read in the input data ***/
  nsurfaces = 0 ;
  for (i = 1 ; i < argc; i++){
    in_fname = argv[i] ;
    printf("reading %s...\n", in_fname) ;
    
    if(!strcmp(srctypestring,"curv")){ /* curvature file */
      if(MRISreadCurvatureFile(mri_surf, in_fname) != 0){
	printf("ERROR: reading curvature file\n");
	exit(1);
      }
      surfVal[nsurfaces] = MRIcopyMRIS(NULL, mri_surf, 0, "curv");
    }
    else if(!strcmp(srctypestring,"paint") || !strcmp(srctypestring,"w")){
      MRISreadValues(mri_surf,in_fname);
      surfVal[nsurfaces] = MRIcopyMRIS(NULL, mri_surf, 0, "val");
    }else{
      printf("ERROR: unknown data file format\n");
      exit(1);
    }
    
    if(surfVal[nsurfaces] == NULL){
      fprintf(stderr, "ERROR loading data values from %s\n", in_fname);
    }
    
    if(nSmoothSteps > 0){
      printf("Smooth input data  by %d steps\n", nSmoothSteps);
      MyMRISsmoothMRI(mri_surf, surfVal[nsurfaces], nSmoothSteps, surfVal[nsurfaces]);
      
    }
    nsurfaces++ ;
  }
  
  // printf("All data read in\n");
  
  ///////////////////////////////////////////////////////////////////////////
  nsurfaces_total = nsurfaces ;   /* all surfaces read in */
  // printf("Compute statistics \n");

  if(nsurfaces_total == 1){
    // printf("Only one surface, compute the stats across the whole surface (excluding the masked areas if mask is given) \n");
    mean = 0; absmean = 0; std = 0; total = 0;
    for(index = 0; index < mri_surf->nvertices; index++){
      if(mri_surf->vertices[index].ripflag) continue;
      scalar = MRIgetVoxVal(surfVal[0],index,0,0,0);
      total++;
      mean += scalar; 
      std += scalar*scalar;
      absmean += ((scalar > 0) ? scalar : (-scalar));
    }
    
    mean /= (total + 1e-30);
    absmean /= (total + 1e-30);
    std  /= (total + 1e-30);
    absstd = std - absmean*absmean;
    if(absstd < 0)absstd = 0;
    else
      absstd = sqrt(absstd);
    
    std = std - mean*mean;
    if(std < 0) std = 0;
    else
      std = sqrt(std);
    
    printf("Stats over surface (%d vertices): mean = %g, std = %g, absmean = %g, absstd = %g \n", total, mean, std, absmean, absstd);

    for(i=0; i < nsurfaces_total; i++){
      MRIfree(&surfVal[i]);
    }
    MRISfree(&mri_surf);
    
    exit(0);
    
  }
  
  mriMean = MRIclone(surfVal[0], NULL);
  mriAbsMean = MRIclone(surfVal[0], NULL);
  mriStd = MRIclone(surfVal[0], NULL);
  mriAbsStd = MRIclone(surfVal[0], NULL);

  if(first_group_size <= 0) first_group_size = nsurfaces_total;

  for(index = 0; index < mri_surf->nvertices; index++){
    mean = 0; std = 0; absmean = 0;
    for(i = 0; i < first_group_size; i++){
      scalar = MRIgetVoxVal(surfVal[i],index,0,0,0);
      if(debugflag == 1 && debugvtx == index){
	printf("Data[%d] = %g\n", i+1, scalar);
      }
      mean += scalar; std += scalar*scalar;
      absmean += ((scalar > 0) ? scalar : (-scalar));
    }
    mean /= (first_group_size + 1e-30);
    absmean /= (first_group_size + 1e-30);
    std  /= (first_group_size + 1e-30);
    absstd = std - absmean*absmean;
    if(absstd < 0)absstd = 0;
    else
      absstd = sqrt(absstd);

    std = std - mean*mean;
    if(std < 0) std = 0;
    else
      std = sqrt(std);

    if(debugflag == 1 && debugvtx == index){
      printf("mean of group1 = %g, std of group1 = %g \n", mean, std);
    }
  
    MRIsetVoxVal(mriMean, index, 0, 0, 0, mean);
    MRIsetVoxVal(mriAbsMean, index, 0, 0, 0, absmean);
    MRIsetVoxVal(mriAbsStd, index, 0, 0, 0, absstd);
    MRIsetVoxVal(mriStd, index, 0, 0, 0, std);
  }

  if(first_group_size < nsurfaces_total - 1){
    //compute the stats for the 2nd group; at least two subjects 

    mriMean2 = MRIclone(surfVal[0], NULL);
    mriAbsMean2 = MRIclone(surfVal[0], NULL);
    mriStd2 = MRIclone(surfVal[0], NULL);
    mriAbsStd2 = MRIclone(surfVal[0], NULL);
    Zscore = MRIclone(surfVal[0], NULL);

    for(index = 0; index < mri_surf->nvertices; index++){
      mean = 0; std = 0; absmean = 0;
      for(i = first_group_size; i < nsurfaces_total; i++){
	scalar = MRIgetVoxVal(surfVal[i],index,0,0,0);
	if(debugflag == 1 && debugvtx == index){
	  printf("Data[%d] = %g\n", i+1, scalar);
	}
	mean += scalar; std += scalar*scalar;
	absmean += ((scalar > 0) ? scalar : (-scalar));
      }
      mean /= (nsurfaces_total - first_group_size + 1e-30);
      absmean /= (nsurfaces_total - first_group_size+ 1e-30);
      std  /= (nsurfaces_total -first_group_size + 1e-30);
      absstd = std - absmean*absmean;
      if(absstd < 0)absstd = 0;
      else
	absstd = sqrt(absstd);
      
      std = std - mean*mean;
      if(std < 0) std = 0;
      else
	std = sqrt(std);
      
      MRIsetVoxVal(mriMean2, index, 0, 0, 0, mean);
      MRIsetVoxVal(mriAbsMean2, index, 0, 0, 0, absmean);
      MRIsetVoxVal(mriAbsStd2, index, 0, 0, 0, absstd);
      MRIsetVoxVal(mriStd2, index, 0, 0, 0, std);

      scalar = MRIgetVoxVal(mriMean, index, 0, 0, 0) - mean;
      scalar /= sqrt(0.5*(MRIgetVoxVal(mriStd, index, 0, 0, 0)*MRIgetVoxVal(mriStd, index, 0, 0, 0) + std*std));
  
      if(debugflag == 1 && debugvtx == index){
	printf("mean of group2 = %g, std of group1 = %g \n", mean, std);
	printf("zscore = %g \n", scalar);
      }
      
      if (scalar < 0) scalar = -scalar;
      MRIsetVoxVal(Zscore, index, 0, 0, 0, scalar);
    }
    
    if(zscore_fname != NULL){
      MRIScopyMRI(mri_surf, Zscore, 0, "curv");
      if(!strcmp(trgtypestring,"paint") || !strcmp(trgtypestring,"w")){
	
	/* This function will remove a zero-valued vertices */
	/* Make sense, since default value is considered as zero */
	/* But it will confuse the processing with matlab! */
	/* So I copy the data to the curv field to force every value is 
	 *  written out
	 */
	
	MRISwriteCurvatureToWFile(mri_surf,zscore_fname);
      }
    }else if(!strcmp(trgtypestring,"curv")){
      MRISwriteCurvature(mri_surf,zscore_fname);
    }else{
      fprintf(stderr, "ERROR unknown output file format.\n");      
    }

    // These are computed but never used yet
    MRIfree(&mriAbsMean2);
    MRIfree(&mriMean2);
    MRIfree(&mriStd2);
    MRIfree(&mriAbsStd2);
    
  }
  
  else{
    MRIScopyMRI(mri_surf, mriStd, 0, "curv");
    if(!strcmp(trgtypestring,"paint") || !strcmp(trgtypestring,"w")){
      
      /* This function will remove a zero-valued vertices */
      /* Make sense, since default value is considered as zero */
      /* But it will confuse the processing with matlab! */
      /* So I copy the data to the curv field to force every value is 
       *  written out
       */
      
      MRISwriteCurvatureToWFile(mri_surf,stat_fname);
      if(mean_fname != NULL){
	MRIScopyMRI(mri_surf, mriMean, 0, "curv");
	MRISwriteCurvatureToWFile(mri_surf,mean_fname);
      }
      if(absmean_fname != NULL){
	MRIScopyMRI(mri_surf, mriAbsMean, 0, "curv");
	MRISwriteCurvatureToWFile(mri_surf,absmean_fname);
      }
      if(absstd_fname != NULL){
	MRIScopyMRI(mri_surf, mriAbsStd, 0, "curv");
	MRISwriteCurvatureToWFile(mri_surf,absstd_fname);
      }
    }else if(!strcmp(trgtypestring,"curv")){
      MRISwriteCurvature(mri_surf,stat_fname);
      if(mean_fname != NULL){
	MRIScopyMRI(mri_surf, mriMean, 0, "curv");
	MRISwriteCurvature(mri_surf,mean_fname);
      }
      if(absmean_fname != NULL){
	MRIScopyMRI(mri_surf, mriAbsMean, 0, "curv");
	MRISwriteCurvature(mri_surf,absmean_fname);
      }
      if(absstd_fname != NULL){
      MRIScopyMRI(mri_surf, mriAbsStd, 0, "curv");
      MRISwriteCurvatureToWFile(mri_surf,absstd_fname);
      }
    }else{
      fprintf(stderr, "ERROR unknown output file format.\n");      
    }
  }
  

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("Std computation took %d minutes and %d seconds.\n", minutes, seconds) ;
  
  for(i=0; i < nsurfaces_total; i++){
    MRIfree(&surfVal[i]);
  }
  MRISfree(&mri_surf);
  MRIfree(&mriMean);
  MRIfree(&mriAbsMean);
  MRIfree(&mriStd);
  MRIfree(&mriAbsStd);

  exit(0);
}

/* --------------------------------------------- */
static void print_version(void)
{
  fprintf(stdout, "%s\n", vcid) ;
  exit(1) ;
}


/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;
  
  option = argv[1] + 1 ;            /* past '-' */

  if (!stricmp(option, "out_fname")
      || !stricmp(option, "out_name")
      || !stricmp(option, "out")
      )
    {
      stat_fname = argv[2];
      printf("using %s as output file \n", stat_fname);
      nargs = 1;
    }
  else if (!stricmp(option, "mean_fname")
      || !stricmp(option, "mean_name")
      || !stricmp(option, "mean")
      )
    {
      mean_fname = argv[2];
      printf("using %s as output file for mean \n", mean_fname);
      nargs = 1;
    }
  else if (!stricmp(option, "mask_fname")
      || !stricmp(option, "mask_name")
      || !stricmp(option, "mask")
      )
    {
      mask_fname = argv[2];
      //      printf("using %s as  file for surface mask (vertices within mask be exclued from statistics) \n", mask_fname);
      nargs = 1;
    }
  else if (!stricmp(option, "absmean_fname")
      || !stricmp(option, "absmean_name")
      || !stricmp(option, "absmean")
      )
    {
      absmean_fname = argv[2];
      printf("using %s as output file for abs-mean \n", absmean_fname);
      nargs = 1;
    }
  else if (!stricmp(option, "absstd_fname")
      || !stricmp(option, "absstd_name")
      || !stricmp(option, "absstd")
      )
    {
      absstd_fname = argv[2];
      printf("using %s as output file for std-of-abs-mean \n", absstd_fname);
      nargs = 1;
    }
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "src_type"))
    {
      srctypestring = argv[2];
      srctype = string_to_type(srctypestring);
      nargs = 1 ;
    }
  else if(!stricmp(option, "surf")
	  || !stricmp(option, "surf_name")
	  || !stricmp(option, "surf_file")
	  ){
    surf_fname = argv[2];
    nargs = 1;
    //    printf("underlying surface is %s\n", surf_fname);
  }
  else if (!stricmp(option, "nsmooth"))
  {
    nSmoothSteps = atoi(argv[2]);
    nargs = 1;
    printf("Perform %d steps of smoothing of input data\n", nSmoothSteps);
  }
  else if (!stricmp(option, "first_group_size"))
  {
    first_group_size = atoi(argv[2]);
    nargs = 1;
    printf("The first %d of input data belongs to one group. \n", first_group_size);
  }
  else if (!stricmp(option, "zscore"))
  {
    zscore_fname = argv[2];
    nargs = 1;
    printf("Output zscore to file %s \n", zscore_fname);
  }
  else if (!stricmp(option, "trg_type"))
  {
    trgtypestring = argv[2];
    trgtype = string_to_type(trgtypestring);
    nargs = 1 ;
  }
  else if(!stricmp(option, "debug"))
    {
      debugflag = 1;
      debugvtx = atoi(argv[2]);
      nargs = 1;
    }
  else{
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    usage_exit(0) ;
    exit(1) ;
  }
  
  return(nargs) ;
}

/*----------------------------------------------------------------------
  Parameters:
  
  Description:
  ----------------------------------------------------------------------*/
static void
usage_exit(int code)
{
  printf("usage: %s [options] <data1> ... <dataN>\n", Progname) ;
  printf("This program computes the std of several data maps (paint-type) defined on the same surface. Other stats can be ouput using -absmean or -mean options.\n");
    
  printf("Options includes:\n");
  printf("\t -debug %%d to set the surface vertex to debug \n");
  printf("\t -surf_name %%s to set the surface filename\n");
  printf("\t -mask_name %%s to set the filename for surface mask \n");
  printf("\t -out_name %%s to set the output filename (std of data)\n");
  printf("\t -mean %%s to set the output filename for mean\n");
  printf("\t -absmean %%s to set the output filename for abs-mean\n");
  printf("\t -absstd %%s to set the output filename for std-of-abs-mean \n");
  printf("\t -zscore %%s to set the output filename for zscore (only if first_group_size > 0) \n");
  printf("\t -first_group_size %%d: how many subjects in the beginning belong to first group \n");
  printf("\t -nsmooth %%d number of smoothing steps\n");
  printf("\t -src_type %%s input surface data format (default = paint) \n");
  printf("\t -trg_type  %%s output format (default = paint) \n");

  exit(code) ;
}

/*-------------------------------------------------------------------*/
MRI *MyMRISsmoothMRI(MRIS *Surf, MRI *Src, int nSmoothSteps, MRI *Targ){
  // In this version, vertex flagged won't be included in the smoothing
  int nnbrs, nthstep, frame, vtx, nbrvtx, nthnbr;
  float val;
  MRI *SrcTmp;
  int truenbr;
  
  if(Surf->nvertices != Src->width){
    printf("ERROR: MyMRISsmooth: Surf/Src dimension mismatch\n");
    return(NULL);
  }
  
  if(Targ == NULL){
    Targ = MRIallocSequence(Src->width, Src->height, Src->depth, 
			    MRI_FLOAT, Src->nframes);
    if(Targ==NULL){
      printf("ERROR: MyMRISsmooth: could not alloc\n");
      return(NULL);
    }
  }
  else{
    if(Src->width   != Targ->width  || 
       Src->height  != Targ->height || 
       Src->depth   != Targ->depth  ||
       Src->nframes != Targ->nframes){
      printf("ERROR: MyMRISsmooth: output dimension mismatch\n");
      return(NULL);
    }
    if(Targ->type != MRI_FLOAT){
      printf("ERROR: MyMRISsmooth: structure passed is not MRI_FLOAT\n");
      return(NULL);
    }
  }
  
  SrcTmp = MRIcopy(Src,NULL);
  for(nthstep = 0; nthstep < nSmoothSteps; nthstep ++){
    //printf("Step = %d\n",nthstep); fflush(stdout);
    
    for(vtx = 0; vtx < Surf->nvertices; vtx++){
      nnbrs = Surf->vertices[vtx].vnum;
      
      for(frame = 0; frame < Targ->nframes; frame ++){
	val = MRIFseq_vox(SrcTmp,vtx,0,0,frame);
	truenbr = 0;
	for(nthnbr = 0; nthnbr < nnbrs; nthnbr++){
	  nbrvtx = Surf->vertices[vtx].v[nthnbr];
	  if(Surf->vertices[nbrvtx].ripflag) continue;
	  truenbr++;
	  val += MRIFseq_vox(SrcTmp,nbrvtx,0,0,frame) ;
	}/* end loop over neighbor */
	//average of itself and all its neighbors
	MRIFseq_vox(Targ,vtx,0,0,frame) = (val/(truenbr+1));
      }/* end loop over frame */
      
    } /* end loop over vertex */
    
    MRIcopy(Targ,SrcTmp);
  }/* end loop over smooth step */
  
  MRIfree(&SrcTmp);
  
  return(Targ);
}
