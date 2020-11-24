/**
 * @brief Maps aparc labels to aseg
 *
 * Maps the cortical labels from the automatic cortical parcellation (aparc)
 * to the automatic segmentation volume (aseg). The result can be used as
 * the aseg would. The algorithm is to find each aseg voxel labeled as
 * cortex (3 and 42) and assign it the label of the closest cortical vertex.
 * If the voxel is not in the ribbon (as defined by mri/lh.ribbon and
 * rh.ribbon), then the voxel is marked as unknown (0). This can be turned
 * off with --noribbon. The cortical parcellation is obtained from
 * subject/label/hemi.aparc.annot which should be based on the
 * curvature.buckner40.filled.desikan_killiany.gcs atlas
 * The aseg is obtained from subject/mri/aseg.mgz and should be based on
 * the RB40_talairach_2005-07-20.gca atlas. If these atlases are used, then
 * the segmentations can be viewed with tkmedit and the FreeSurferColorLUT.txt
 * color table found in $FREESURFER_HOME. These are the default atlases
 * used by recon-all.
 */
/*
 * Original Author: Doug Greve
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
#include "macros.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "annotation.h"
#include "version.h"
#include "mrisegment.h"
#include "cma.h"
#include "gca.h"
#include "cmdargs.h"
#ifdef _OPENMP
#include "romp_support.h"
#endif

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
int FindClosestLRWPVertexNo(int c, int r, int s,
                            int *lhwvtx, int *lhpvtx,
                            int *rhwvtx, int *rhpvtx,
                            MATRIX *Vox2RAS,
                            MRIS *lhwite,  MRIS *lhpial,
                            MRIS *rhwhite, MRIS *rhpial,
                            MHT *lhwhite_hash, MHT *lhpial_hash,
                            MHT *rhwhite_hash, MHT *rhpial_hash);
int CCSegment(MRI *seg, int segid, int segidunknown);

int main(int argc, char *argv[]) ;

const char *Progname = NULL;
static char *relabel_gca_name = NULL ;
static char *relabel_norm_name = NULL ;
static char *relabel_xform_name = NULL ;
static char *relabel_label_intensities_name = NULL ;

static char *SUBJECTS_DIR = NULL;
static char *subject = NULL;
static char *OutASegFile = NULL;
static char *OutAParcFile = NULL;
static char *OutDistFile = NULL;
static int debug = 0;
static int UseRibbon = 0;
static int UseNewRibbon = 1;
static MRI *ASeg, *filled, *mritmp;
static MRI *AParc;
static MRI *Dist;
static MRI *lhRibbon=NULL,*rhRibbon=NULL,*RibbonSeg;
static MRIS *lhwhite, *rhwhite;
static MRIS *lhpial, *rhpial;
static MHT *lhwhite_hash, *rhwhite_hash;
static MHT *lhpial_hash, *rhpial_hash;
static int  lhwvtx, lhpvtx, rhwvtx, rhpvtx;
static MATRIX *Vox2RAS;
static float dmaxctx = 5.0;
static int LabelWM=0;
static int LabelHypoAsWM=0;
static int RipUnknown = 0;

static char tmpstr[2000];
static char annotfile[1000];
static const char *annotname = "aparc";
static const char *asegname = "aseg";
static int baseoffset = 0;
static float hashres = 16;

static int normal_smoothing_iterations = 10 ;
int crsTest = 0, ctest=0, rtest=0, stest=0;
int UseHash = 1;
int DoLH=1, DoRH=1, LHOnly=0, RHOnly=0;

char *CtxSegFile = NULL;
MRI *CtxSeg = NULL;

int FixParaHipWM = 1;
double BRFdotCheck(MRIS *surf, int vtxno, int c, int r, int s, MRI *AParc);
int nthreads=1;

/*--------------------------------------------------*/
int main(int argc, char **argv)
{
  int nargs, err, c, nctx, annot,vtxno,nripped;
  int annotid;
  int nbrute=0;
  MRI    *mri_fixed = NULL, *mri_lh_dist, *mri_rh_dist, *mri_dist=NULL;
  TRANSFORM *xform  = NULL;
  GCA *gca = NULL ;

  nargs = handleVersionOption(argc, argv, "mri_aparc2aseg");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0)
  {
    usage_exit();
  }

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR==NULL)
  {
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

#ifdef _OPENMP
  printf("%d avail.processors, using %d\n",omp_get_num_procs(),omp_get_max_threads());
#endif

  if(DoLH){
    /* ------ Load subject's lh white surface ------ */
    sprintf(tmpstr,"%s/%s/surf/lh.white",SUBJECTS_DIR,subject);
    printf("\nReading lh white surface \n %s\n",tmpstr);
    lhwhite = MRISread(tmpstr);
    if (lhwhite == NULL)  {
      fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
      exit(1);
    }
    /* ------ Load subject's lh pial surface ------ */
    sprintf(tmpstr,"%s/%s/surf/lh.pial",SUBJECTS_DIR,subject);
    printf("\nReading lh pial surface \n %s\n",tmpstr);
    lhpial = MRISread(tmpstr);
    if (lhpial == NULL)  {
      fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
      exit(1);
    }
    if (lhwhite->nvertices != lhpial->nvertices)  {
      printf("ERROR: lh white and pial have a different number of "
	     "vertices (%d,%d)\n",
	     lhwhite->nvertices,lhpial->nvertices);
      exit(1);
    }

    /* ------ Load lh annotation ------ */
    sprintf(annotfile,"%s/%s/label/lh.%s.annot",SUBJECTS_DIR,subject,annotname);
    printf("\nLoading lh annotations from %s\n",annotfile);
    err = MRISreadAnnotation(lhwhite, annotfile);
    if (err)  {
      printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile);
      exit(1);
    }
    if(lhwhite->ct) printf("Have color table for lh white annotation\n");
    if (UseRibbon)  {
      sprintf(tmpstr,"%s/%s/mri/lh.ribbon.mgz",SUBJECTS_DIR,subject);
      printf("Loading lh ribbon mask from %s\n",tmpstr);
      lhRibbon = MRIread(tmpstr);
      if (lhRibbon == NULL)    {
	printf("ERROR: loading %s\n",tmpstr);
	exit(1);
      }
    }
    if(RipUnknown)  {
      printf("Ripping vertices labeled as unkown\n");
      nripped = 0;
      for (vtxno = 0; vtxno < lhwhite->nvertices; vtxno++) {
	annot = lhwhite->vertices[vtxno].annotation;
	CTABfindAnnotation(lhwhite->ct, annot, &annotid);
	// Sometimes the annotation will be "none" indicated by
	// annotid = -1. We interpret this as "unknown".
	if (annotid == 0 || annotid == -1) {
	  lhwhite->vertices[vtxno].ripflag = 1;
	  lhpial->vertices[vtxno].ripflag = 1;
	  nripped++;
	}
      }
      printf("Ripped %d vertices from left hemi\n",nripped);
    }
    printf("\n");
    printf("Building hash of lh white\n");
    lhwhite_hash = MHTcreateVertexTable_Resolution(lhwhite, CURRENT_VERTICES,hashres);
    printf("\n");
    printf("Building hash of lh pial\n");
    lhpial_hash = MHTcreateVertexTable_Resolution(lhpial, CURRENT_VERTICES,hashres);
  }

  if(DoRH){
    /* ------ Load subject's rh white surface ------ */
    sprintf(tmpstr,"%s/%s/surf/rh.white",SUBJECTS_DIR,subject);
    printf("\nReading rh white surface \n %s\n",tmpstr);
    rhwhite = MRISread(tmpstr);
    if (rhwhite == NULL)      {
      fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
      exit(1);
    }
    /* ------ Load subject's rh pial surface ------ */
    sprintf(tmpstr,"%s/%s/surf/rh.pial",SUBJECTS_DIR,subject);
    printf("\nReading rh pial surface \n %s\n",tmpstr);
    rhpial = MRISread(tmpstr);
    if (rhpial == NULL)  {
      fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
      exit(1);
    }
    if (rhwhite->nvertices != rhpial->nvertices)  {
      printf("ERROR: rh white and pial have a different "
	     "number of vertices (%d,%d)\n",
	     rhwhite->nvertices,rhpial->nvertices);
      exit(1);
    }

    /* ------ Load rh annotation ------ */
    sprintf(annotfile,"%s/%s/label/rh.%s.annot",SUBJECTS_DIR,subject,annotname);
    printf("\nLoading rh annotations from %s\n",annotfile);
    err = MRISreadAnnotation(rhwhite, annotfile);
    if (err)  {
      printf("ERROR: MRISreadAnnotation() failed %s\n",annotfile);
      exit(1);
    }
    if(rhwhite->ct)    printf("Have color table for rh white annotation\n");
    if (UseRibbon)  {
      sprintf(tmpstr,"%s/%s/mri/rh.ribbon.mgz",SUBJECTS_DIR,subject);
      printf("Loading rh ribbon mask from %s\n",tmpstr);
      rhRibbon = MRIread(tmpstr);
      if (rhRibbon == NULL)	{
	printf("ERROR: loading  %s\n",tmpstr);
	exit(1);
      }
    }
    if (RipUnknown)  {
      printf("Ripping vertices labeled as unkown\n");
      nripped = 0;
      for (vtxno = 0; vtxno < rhwhite->nvertices; vtxno++)	{
	annot = rhwhite->vertices[vtxno].annotation;
	CTABfindAnnotation(rhwhite->ct, annot, &annotid);
	if (annotid == 0 || annotid == -1)      {
	  rhwhite->vertices[vtxno].ripflag = 1;
	  rhpial->vertices[vtxno].ripflag = 1;
	  nripped++;
	}
      }
      printf("Ripped %d vertices from right hemi\n",nripped);
    }
    printf("\n");
    printf("Building hash of rh white\n");
    rhwhite_hash = MHTcreateVertexTable_Resolution(rhwhite, CURRENT_VERTICES,hashres);
    printf("\n");
    printf("Building hash of rh pial\n");
    rhpial_hash = MHTcreateVertexTable_Resolution(rhpial, CURRENT_VERTICES,hashres);
  }

  if(UseNewRibbon){
    sprintf(tmpstr,"%s/%s/mri/ribbon.mgz",SUBJECTS_DIR,subject);
    printf("Loading ribbon segmentation from %s\n",tmpstr);
    RibbonSeg = MRIread(tmpstr);
    if (RibbonSeg == NULL){
      printf("ERROR: loading %s\n",tmpstr);
      exit(1);
    }
  }

  if(LabelHypoAsWM){
    /* This section used to use filled.mgz instead of ribbon.mgz. This is used
     to determine whether a hypointensity is WM. filled.mgz is not guaranteed
     to have a non-zero value in a hypointensity.  Not sure why I used
     filled in the first place, esp since ribbon is used above. */
    sprintf(tmpstr,"%s/%s/mri/ribbon.mgz",SUBJECTS_DIR,subject);
    printf("Loading filled from %s\n",tmpstr);
    filled = MRIread(tmpstr);
    if (filled == NULL){
      printf("ERROR: loading filled %s\n",tmpstr);
      exit(1);
    }
    int MatchList[2];
    MatchList[0] = 2;
    MatchList[1] = 41;
    mritmp = MRIbinarizeMatch(filled, MatchList, 2, 0, NULL);
    MRIfree(&filled);
    filled = mritmp;
  }

  /* ------ Load ASeg ------ */
  if(!fio_FileExistsReadable(asegname)){
    sprintf(tmpstr,"%s/%s/mri/%s.mgz",SUBJECTS_DIR,subject,asegname);
    if (!fio_FileExistsReadable(tmpstr))      {
      sprintf(tmpstr,"%s/%s/mri/%s.mgh",SUBJECTS_DIR,subject,asegname);
      if (!fio_FileExistsReadable(tmpstr))	{
	sprintf(tmpstr,"%s/%s/mri/aseg/COR-.info",SUBJECTS_DIR,subject);
	if (!fio_FileExistsReadable(tmpstr))	      {
	  printf("ERROR: cannot find aseg %s\n",asegname);
	  exit(1);
	}
	else {
	  sprintf(tmpstr,"%s/%s/mri/aseg/",SUBJECTS_DIR,subject);
	}
      }
    }
  }
  else {
    strcpy(tmpstr,asegname);
  }
  
  printf("\nLoading aseg from %s\n",tmpstr);
  ASeg = MRIread(tmpstr);
  if (ASeg == NULL)  {
    printf("ERROR: loading aseg %s\n",tmpstr);
    exit(1);
  }
  mritmp = MRIchangeType(ASeg,MRI_INT,0,0,1);
  MRIfree(&ASeg);
  ASeg = mritmp;

  if(DoLH){
    mri_lh_dist = MRIcloneDifferentType(ASeg, MRI_FLOAT) ;
    MRIScomputeDistanceToSurface(lhwhite, mri_lh_dist, mri_lh_dist->xsize) ;
    if(LHOnly) mri_dist = mri_lh_dist;
  }
  if(DoRH){
    mri_rh_dist = MRIcloneDifferentType(ASeg, MRI_FLOAT) ;
    MRIScomputeDistanceToSurface(rhwhite, mri_rh_dist, mri_rh_dist->xsize) ;
    if(RHOnly) mri_dist = mri_rh_dist;
  }
  if(DoLH && DoRH){
    mri_dist = MRImin(mri_lh_dist, mri_rh_dist, NULL) ;
    MRIfree(&mri_lh_dist) ; 
    MRIfree(&mri_rh_dist) ;
  }

  if (relabel_norm_name)  {
    mri_fixed = MRIcloneDifferentType(ASeg, MRI_UCHAR) ;
    MRIsetValues(mri_fixed, 1) ;
  }
  if (CtxSegFile)  {
    printf("Loading Ctx Seg File %s\n",CtxSegFile);
    CtxSeg = MRIread(CtxSegFile);
    if (CtxSeg == NULL)
    {
      exit(1);
    }
  }

  AParc = MRIclone(ASeg,NULL);
  if (OutDistFile != NULL)  {
    Dist = MRIclone(ASeg,NULL);
    mritmp = MRIchangeType(Dist,MRI_FLOAT,0,0,0);
    if (mritmp == NULL)
    {
      printf("ERROR: could change type\n");
      exit(1);
    }
    MRIfree(&Dist);
    Dist = mritmp;
  }

  Vox2RAS = MRIxfmCRS2XYZtkreg(ASeg);
  printf("ASeg Vox2RAS: -----------\n");
  MatrixPrint(stdout,Vox2RAS);
  printf("-------------------------\n");

  if (crsTest)  {
    printf("Testing point %d %d %d\n",ctest,rtest,stest);
    err = FindClosestLRWPVertexNo(ctest,rtest,stest,
                                  &lhwvtx, &lhpvtx,
                                  &rhwvtx, &rhpvtx, Vox2RAS,
                                  lhwhite,  lhpial,
                                  rhwhite, rhpial,
                                  lhwhite_hash, lhpial_hash,
                                  rhwhite_hash, rhpial_hash);

    printf("Result: err = %d\n",err);
    exit(err);
  }

  nctx = 0;
  annot = 0;
  annotid = 0;
  nbrute = 0;

  if(DoLH){
    MRISsmoothSurfaceNormals(lhwhite, normal_smoothing_iterations) ;
    MRISsmoothSurfaceNormals(lhpial, normal_smoothing_iterations) ;
  }
  if(DoRH){
    MRISsmoothSurfaceNormals(rhpial, normal_smoothing_iterations) ;
    MRISsmoothSurfaceNormals(rhwhite, normal_smoothing_iterations) ;
  }

  if (relabel_gca_name != NULL)    // reclassify voxels interior to white that are likely to be something else
  {
    MRI    *mri_norm;
    FILE   *fp ;
    int    *labels, nlines, i, mean, label;
    float  *intensities ;
    char   *cp, line[STRLEN], label_name[STRLEN] ;

    printf("relabeling unlikely voxels in interior of white matter\n") ;

    mri_norm = MRIread(relabel_norm_name) ;
    if (mri_norm == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load norm volume from %s\n", relabel_norm_name) ;

#if 0
    // BF included this in his relabel code, but not clear what it does
    // variables don't seem to be used and shadow variables in the main code
    MRI    *mri_rh_dist, *mri_lh_dist, *mri_dist;
    if(DoLH){
      mri_lh_dist = MRIcloneDifferentType(mri_norm, MRI_FLOAT) ;
      MRIScomputeDistanceToSurface(lhwhite, mri_lh_dist, mri_lh_dist->xsize) ; 
      if(LHOnly) mri_dist = mri_lh_dist;
    }
    if(DoRH){
      mri_rh_dist = MRIcloneDifferentType(mri_norm, MRI_FLOAT) ;
      MRIScomputeDistanceToSurface(rhwhite, mri_rh_dist, mri_rh_dist->xsize) ; 
      if(RHOnly) mri_dist = mri_rh_dist;
    }
    if(DoLH && DoRH){
      mri_dist = MRIcombineDistanceTransforms(mri_lh_dist, mri_rh_dist, NULL) ;
      MRIfree(&mri_lh_dist) ; MRIfree(&mri_rh_dist) ;
    }
#endif
    
    xform = TransformRead(relabel_xform_name) ;
    if (xform == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load transform from %s\n", relabel_xform_name) ;

    gca = GCAread(relabel_gca_name) ;
    if (gca == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load gca from %s\n", relabel_gca_name) ;

    fp = fopen(relabel_label_intensities_name, "r") ;
    if (!fp)
      ErrorExit(ERROR_NOFILE, "%s: could not read %s",
                Progname, relabel_label_intensities_name) ;
    cp = fgetl(line, 199, fp) ;
    nlines = 0 ;
    while (cp)
    {
      nlines++ ;
      cp = fgetl(line, 199, fp) ;
    }
    rewind(fp) ;
    printf("reading %d labels from %s\n", nlines,relabel_label_intensities_name) ;
    labels = (int *)calloc(nlines, sizeof(int)) ;
    intensities = (float *)calloc(nlines, sizeof(float)) ;
    cp = fgetl(line, 199, fp) ;
    for (i = 0 ; i < nlines ; i++)
    {
      sscanf(cp, "%d %s %*f %*f %d", &label, label_name, &mean) ;
      labels[i] = label ;
      intensities[i] = mean ;
      if (labels[i] == Left_Cerebral_White_Matter)
      {
        DiagBreak() ;
      }
      cp = fgetl(line, 199, fp) ;
    }
    GCArenormalizeIntensities(gca, labels, intensities, nlines) ;
    free(labels) ;
    free(intensities) ;

    TransformInvert(xform, mri_norm) ;
    // edit GCA to disallow cortical labels at points interior to white that we are relabeling
    {
      int        x, y, z, n ;
      GCA_PRIOR *gcap ;

      for (x = 0 ; x < mri_norm->width ; x++)
	for (y = 0 ; y < mri_norm->height ; y++)
	  for (z = 0 ; z < mri_norm->depth ; z++)
	  {
	    if (x  == Gx && y == Gy && z == Gz)
	      DiagBreak() ;
	    if (MRIgetVoxVal(mri_fixed, x, y, z, 0) > 0)
	      continue ;
	    gcap = getGCAP(gca, mri_norm, xform, x, y, z) ;
	    if (gcap == NULL)
	      continue ;
	    for (n = 0 ; n < gcap->nlabels ; n++)
	      if (IS_CORTEX(gcap->labels[n]))
	      { 
		int n2 ;
		double total_p ;
		gcap->priors[n] = 1e-5 ;
		for (total_p = 0.0, n2 = 0 ; n2 < gcap->nlabels ; n2++)
		  total_p += gcap->priors[n2] ;
		for (n2 = 0 ; n2 < gcap->nlabels ; n2++)
		  gcap->priors[n2] /= total_p ;
	      }
	  }
    }
    Ggca_x = Gx ; Ggca_y = Gy ; Ggca_z = Gz ; // diagnostics
  }

  // Go through each voxel in the aseg
  printf("\nLabeling Slice (%d)\n",ASeg->width);
  
  MHT_maybeParallel_begin();
  #ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+ : nbrute, nctx)
  #endif
  for (c=0; c < ASeg->width; c++){
    int r,s,asegid,IsWM,IsCblumCtx,IsCortex,IsHypo, RibbonVal,lhRibbonVal,rhRibbonVal,lhwvtx,rhwvtx,lhpvtx,rhpvtx;
    int annot=0,annotid,hemi=0,segval=0;
    double dthresh,dist,dot;
    float dlhw,drhw,dlhp,drhp,dmin=1e7;
    struct { float x,y,z; } vtx;
    MATRIX *CRS, *RAS;

    printf("%3d ",c);
    if (c%20 ==19) printf("\n");
    fflush(stdout);
    for (r=0; r < ASeg->height; r++)    {
      for (s=0; s < ASeg->depth; s++)      {

	if (c == Gx && r == Gy && s == Gz)
	  DiagBreak() ;

        asegid = MRIgetVoxVal(ASeg,c,r,s,0);
	if(LHOnly && (asegid == Right_Cerebral_Cortex || asegid == Right_Cerebral_White_Matter)) continue;
	if(RHOnly && (asegid ==  Left_Cerebral_Cortex || asegid ==  Left_Cerebral_White_Matter))continue;
	IsCortex = IS_CORTEX(asegid) ;
	IsHypo = IS_HYPO(asegid) ;
        if(asegid == Left_Cerebral_White_Matter || asegid == Right_Cerebral_White_Matter)  IsWM = 1;
        else                                                                               IsWM = 0;
	//  what is 172??? (BRF)
        if(asegid == Left_Cerebellum_Cortex || asegid == Right_Cerebellum_Cortex || asegid == 172)  IsCblumCtx = 1;
        else                                                                                        IsCblumCtx = 0;
        if(IsHypo && LabelHypoAsWM && MRIgetVoxVal(filled,c,r,s,0)) IsWM = 1;

        // integrate surface information
        //
        // Only Do This for GM,WM or Unknown labels in the ASEG !!!
        // priority is given to the ribbon computed from the surface namely
	//  aseg=SubCortGM => keep GM (unless possibly CblumCtx)
        //  ribbon=GM => GM
        //  aseg=GM AND ribbon=WM => WM
        //  ribbon=UNKNOWN => UNKNOWN
	// Note: have to be careful when setting values that are "Unknown" in the ribbon
	// because the aseg.presurf is often not very accurate for cortex, eg, sometimes an entire
	// gyrus (including sulcal CSF) can be labeled as WM, obviously, you don't want to
	// just transfer the label. 
        if(UseNewRibbon){
	  if(IsCortex || IsWM || (asegid==Unknown || asegid == CSF) || IsCblumCtx) {
	    RibbonVal = MRIgetVoxVal(RibbonSeg,c,r,s,0);
	    if(!IsCblumCtx && asegid != CSF) MRIsetVoxVal(ASeg,c,r,s,0, RibbonVal);
	    if(RibbonVal==Left_Cerebral_White_Matter || RibbonVal==Right_Cerebral_White_Matter) {
	      // Ribbon says it is WM
	      IsWM = 1;
	      IsCortex = 0;
	    }
	    else if(RibbonVal==Left_Cerebral_Cortex || RibbonVal==Right_Cerebral_Cortex) {
	      // Ribbon says it is Ctx
	      IsWM = 0;
	      IsCortex = 1;
	      if(IsCblumCtx) MRIsetVoxVal(ASeg,c,r,s,0, RibbonVal);
	    }
	    if(RibbonVal==Unknown)  {
	      // Ribbon says it is unknown
	      IsWM = 0;
	      IsCortex = 0;
	    }
	  }
        }

	if (c == Gx && r == Gy && s == Gz)
	  DiagBreak() ;
	if (gca &&
	    (GCAisPossible(gca, ASeg, Left_Hippocampus, xform, c, r, s, 0)  ||
	     GCAisPossible(gca, ASeg, Right_Hippocampus, xform, c, r, s, 0)  ||
	     GCAisPossible(gca, ASeg, Left_Amygdala, xform, c, r, s, 0)  ||
	     GCAisPossible(gca, ASeg, Right_Amygdala, xform, c, r, s, 0)))
	  dthresh = -1.5 ;  // don't trust surfaces much in MTL
	else
	  dthresh = 0.5 ;
	dist = MRIgetVoxVal(mri_dist, c, r, s, 0) ;
	if (IsWM && (asegid == 0 || asegid == CSF) && mri_fixed != NULL && dist < dthresh)  // interior to white matter but labeled unknown
	  MRIsetVoxVal(mri_fixed, c, r, s, 0, 0) ;     // allow it to be relabeled below

        // If it's not labeled as cortex or wm in the aseg, skip
        if(!IsCortex && !IsWM)continue;

        // If it's wm but not labeling wm, skip
        if(IsWM && !LabelWM) continue;

        // Check whether this point is in the ribbon
        if(UseRibbon) {
	  lhRibbonVal = 0;
	  rhRibbonVal = 0;
	  if(DoLH) lhRibbonVal = MRIgetVoxVal(lhRibbon,c,r,s,0);
          if(DoRH) rhRibbonVal = MRIgetVoxVal(rhRibbon,c,r,s,0);
          if(IsCortex) {
            // ASeg says it's in cortex, or other logic says so
            if (lhRibbonVal < 0.5 && rhRibbonVal < 0.5) {
              // but it is not part of the ribbon,
              // so set it to unknown (0) and go to the next voxel.
              MRIsetVoxVal(ASeg,c,r,s,0,0);
              continue;
            }
          }
        }

        // Convert the CRS to RAS
	CRS = MatrixAlloc(4,1,MATRIX_REAL);
	CRS->rptr[4][1] = 1;
	RAS = MatrixAlloc(4,1,MATRIX_REAL);
	RAS->rptr[4][1] = 1;
        CRS->rptr[1][1] = c;
        CRS->rptr[2][1] = r;
        CRS->rptr[3][1] = s;
        RAS = MatrixMultiply(Vox2RAS,CRS,RAS);
        vtx.x = RAS->rptr[1][1];
        vtx.y = RAS->rptr[2][1];
        vtx.z = RAS->rptr[3][1];

        // Get the index of the closest vertex in the
        // lh.white, lh.pial, rh.white, rh.pial
        if(UseHash) {
	  if(DoLH){
	    lhwvtx = MHTfindClosestVertexNoXYZ(lhwhite_hash,lhwhite,vtx.x,vtx.y,vtx.z,&dlhw);
	    lhpvtx = MHTfindClosestVertexNoXYZ(lhpial_hash, lhpial, vtx.x,vtx.y,vtx.z,&dlhp);
	  } else {
	    lhwvtx = -1;
	    lhpvtx = -1;
	  }
	  if(DoRH){
	    rhwvtx = MHTfindClosestVertexNoXYZ(rhwhite_hash,rhwhite,vtx.x,vtx.y,vtx.z,&drhw);
	    rhpvtx = MHTfindClosestVertexNoXYZ(rhpial_hash, rhpial, vtx.x,vtx.y,vtx.z,&drhp);
	  } else {
	    rhwvtx = -1;
	    rhpvtx = -1;
	  }
          if (lhwvtx < 0 && lhpvtx < 0 && rhwvtx < 0 && rhpvtx < 0) {
            /*
            printf("  Could not map to any surface with hash table:\n");
            printf("  crs = %d %d %d, ras = %6.4f %6.4f %6.4f \n",
            c,r,s,vtx.x,vtx.y,vtx.z);
            printf("  Using brute force search %d ... \n",nbrute);
            fflush(stdout);
            */
	    if(DoLH){
	      lhwvtx = MRISfindClosestVertex(lhwhite,vtx.x,vtx.y,vtx.z,&dlhw, CURRENT_VERTICES);
	      lhpvtx = MRISfindClosestVertex(lhpial,vtx.x,vtx.y,vtx.z,&dlhp, CURRENT_VERTICES);
	    }
	    if(DoRH){
	      rhwvtx = MRISfindClosestVertex(rhwhite,vtx.x,vtx.y,vtx.z,&drhw, CURRENT_VERTICES);
	      rhpvtx = MRISfindClosestVertex(rhpial,vtx.x,vtx.y,vtx.z,&drhp, CURRENT_VERTICES);
	    }
            nbrute ++;
            //exit(1);
          }
        }
        else
        {
	  if(DoLH){
	    lhwvtx = MRISfindClosestVertex(lhwhite,vtx.x,vtx.y,vtx.z,&dlhw, CURRENT_VERTICES);
	    lhpvtx = MRISfindClosestVertex(lhpial,vtx.x,vtx.y,vtx.z,&dlhp, CURRENT_VERTICES);
	  } else {
	    lhwvtx = -1;
	    lhpvtx = -1;
	  }
	  if(DoRH){
	    rhwvtx = MRISfindClosestVertex(rhwhite,vtx.x,vtx.y,vtx.z,&drhw, CURRENT_VERTICES);
	    rhpvtx = MRISfindClosestVertex(rhpial,vtx.x,vtx.y,vtx.z,&drhp, CURRENT_VERTICES);
	  } else {
	    rhwvtx = -1;
	    rhpvtx = -1;
	  }
        }

	/* added some checks here to make sure closest vertex (usually pial but can be white) isn't on
	   the other bank of a sulcus or through a thin white matter strand. This removes inaccurate voxels
	   that used to speckle the aparc+aseg
	*/
        if (lhwvtx < 0)       dlhw = 1000000000000000.0;
	else if(!LabelWM){
	  dot = BRFdotCheck(lhwhite,lhwvtx,c,r,s,AParc);
	  if (dot < 0) 
	  {
	    if (MRIneighbors(ASeg, c, r, s, Left_Cerebral_Cortex) > 0) // only do expensive check if it is possible
	      dlhw = MRISfindMinDistanceVertexWithDotCheck(lhwhite, c, r, s, AParc, 1, &lhwvtx) ;
	    else
	      dlhw = 1000000000000000.0;
	  }
	}

        if (lhpvtx < 0) dlhp = 1000000000000000.0;
	else if(!LabelWM){
	  dot = BRFdotCheck(lhpial,lhpvtx,c,r,s,AParc);
	  if (dot > 0)   // pial surface normal should point in same direction as vector from voxel to vertex
	  {
	    if (MRIneighbors(ASeg, c, r, s, Left_Cerebral_Cortex) > 0) // only do expensive check if it is possible
	      dlhp = MRISfindMinDistanceVertexWithDotCheck(lhpial, c, r, s, AParc, -1, &lhpvtx) ;
	    else
	      dlhp = 1000000000000000.0;
	  }
	}

        if (rhwvtx < 0) drhw = 1000000000000000.0;
	else if(!LabelWM){
	  dot = BRFdotCheck(rhwhite,rhwvtx,c,r,s,AParc);
	  if (dot < 0)
	  {
	    if (MRIneighbors(ASeg, c, r, s, Right_Cerebral_Cortex) > 0) // only do expensive check if it is possible
	      drhw = MRISfindMinDistanceVertexWithDotCheck(rhwhite, c, r, s, AParc, 1, &rhwvtx) ;
	    else
	      drhw = 1000000000000000.0;
	  }
	}
        if (rhpvtx < 0) drhp = 1000000000000000.0;
	else if(!LabelWM){
	  dot = BRFdotCheck(rhpial,rhpvtx,c,r,s,AParc);
	  if (dot > 0) 
	  {
	    if (MRIneighbors(ASeg, c, r, s, Right_Cerebral_Cortex) > 0) // only do expensive check if it is possible
	      drhp = MRISfindMinDistanceVertexWithDotCheck(rhpial, c, r, s, AParc, -1, &rhpvtx) ;
	    else
	      drhp = 1000000000000000.0;
	  }
	}

        if (dlhw <= dlhp && dlhw < drhw && dlhw < drhp && lhwvtx >= 0) {
          annot = lhwhite->vertices[lhwvtx].annotation;
          hemi = 1;
          if (lhwhite->ct) CTABfindAnnotation(lhwhite->ct, annot, &annotid);
          else annotid = annotation_to_index(annot);
          dmin = dlhw;
        }
        if (dlhp < dlhw && dlhp < drhw && dlhp < drhp && lhpvtx >= 0) {
          annot = lhwhite->vertices[lhpvtx].annotation;
          hemi = 1;
          if (lhwhite->ct) CTABfindAnnotation(lhwhite->ct, annot, &annotid);
          else annotid = annotation_to_index(annot);
          dmin = dlhp;
        }

        if (drhw < dlhp && drhw < dlhw && drhw <= drhp && rhwvtx >= 0)
        {
          annot = rhwhite->vertices[rhwvtx].annotation;
          hemi = 2;
          if (rhwhite->ct) CTABfindAnnotation(rhwhite->ct, annot, &annotid);
          else annotid = annotation_to_index(annot);
          dmin = drhw;
        }
        if (drhp < dlhp && drhp < drhw && drhp < dlhw && rhpvtx >= 0)
        {
          annot = rhwhite->vertices[rhpvtx].annotation;
          hemi = 2;
          if (rhwhite->ct) CTABfindAnnotation(rhwhite->ct, annot, &annotid);
          else annotid = annotation_to_index(annot);
          dmin = drhp;
        }

        // Sometimes the annotation will be "none" indicated by
        // annotid = -1. We interpret this as "unknown".
        if (annotid == -1) annotid = 0;

	/* If the cortical label is "unkown", it is difficult to
	   determine what to put here. If the aseg says it is WM, then
	   that is kept. If the aseg says it is GM, then it is given
	   "ctx-?h-unknown". These voxels can show up in funny places
	   (eg, between hippo and amyg), so this is really just a
	   hack. The real fix should be the surface creation or the
	   aseg. */
	if(annotid == 0 && !LabelWM){
	  if(asegid == Left_Cerebral_Cortex)  MRIsetVoxVal(ASeg,c,r,s,0,1000);
	  else if(asegid == Right_Cerebral_Cortex) MRIsetVoxVal(ASeg,c,r,s,0,2000);
	  else MRIsetVoxVal(ASeg,c,r,s,0,asegid);
	  MatrixFree(&CRS);
	  MatrixFree(&RAS);
	  continue;
	  //{if(hemi == 1) MRIsetVoxVal(ASeg,c,r,s,0,Left_Cerebral_White_Matter);
	  //if(hemi == 2) MRIsetVoxVal(ASeg,c,r,s,0,Right_Cerebral_White_Matter);
	  //continue;
	}

        // why was this here in the first place?
        /*
               if (annotid == 0 &&
                   lhwvtx >= 0 &&
                   lhpvtx >= 0 &&
                   rhwvtx >= 0 &&
                   rhpvtx >= 0) {
                 printf("%d %d %d %d\n",
                        lhwhite->vertices[lhwvtx].ripflag,
                        lhpial->vertices[lhpvtx].ripflag,
                        rhwhite->vertices[rhwvtx].ripflag,
                        rhpial->vertices[rhpvtx].ripflag);
          } */

        if ( IsCortex && hemi == 1) segval = annotid+1000 + baseoffset;  //ctx-lh
        if ( IsCortex && hemi == 2) segval = annotid+2000 + baseoffset;  //ctx-rh
        if (!IsCortex && hemi == 1) segval = annotid+3000 + baseoffset;  // wm-lh
        if (!IsCortex && hemi == 2) segval = annotid+4000 + baseoffset;  // wm-rh
        if (!IsCortex && dmin > dmaxctx && hemi == 1) 
	{
	  if (dmin > 2*fabs(dist))  // in medial wall (dist is to a ripped vertex so not reflected in dmin)
	    segval = asegid ;
	  else
	    segval = Left_Unsegmented_WM ;
	}
        if (!IsCortex && dmin > dmaxctx && hemi == 2) 
	{
	  if (dmin > 2*fabs(dist))  // in medial wall (dist is to a ripped vertex so not reflected in dmin)
	    segval = asegid ;
	  else
	    segval = Right_Unsegmented_WM;
	}

	if (LabelWM) // unsegmented wm shouldn't be labeled as such
	{
	  switch (segval)
	  {
	  case Left_Cerebral_White_Matter:
	  case wm_lh_unknown:
	    segval = Left_Unsegmented_WM ;
	    break ;
	  case Right_Cerebral_White_Matter:
	  case wm_rh_unknown:
	    segval = Right_Unsegmented_WM ;
	    break ;
	  default:
	    break ;
	  }
	}

        // This is a hack for getting the right cortical seg with --rip-unknown
        // The aparc+aseg should be passed as CtxSeg. Used with WMParc
        if (IsCortex && CtxSeg) segval = MRIgetVoxVal(CtxSeg,c,r,s,0);

        MRIsetVoxVal(ASeg,c,r,s,0,segval);
        MRIsetVoxVal(AParc,c,r,s,0,annot);
        if (OutDistFile != NULL) MRIsetVoxVal(Dist,c,r,s,0,dmin);

        if (debug || annotid == -1) {
          // Gets here when there is no label at the found vertex.
          // This is different than having a vertex labeled as "unknown"
          if (!debug)
          {
	    MatrixFree(&CRS);
	    MatrixFree(&RAS);
            continue;
          }
          printf("\n");
          printf("Found closest vertex, but it has no label.\n");
          printf("aseg id = %d\n",asegid);
          printf("crs = %d %d %d, ras = %6.4f %6.4f %6.4f \n",
                 c,r,s,vtx.x,vtx.y,vtx.z);
          if (lhwvtx > 0) printf("lhw  %d  %7.5f     %6.4f  %6.4f  %6.4f\n",
                                   lhwvtx, dlhw,
                                   lhwhite->vertices[lhwvtx].x,
                                   lhwhite->vertices[lhwvtx].y,
                                   lhwhite->vertices[lhwvtx].z);
          if (lhpvtx > 0) printf("lhp  %d  %7.5f     %6.4f  %6.4f  %6.4f\n",
                                   lhpvtx, dlhp,
                                   lhpial->vertices[lhpvtx].x,
                                   lhpial->vertices[lhpvtx].y,
                                   lhpial->vertices[lhpvtx].z);
          if (rhwvtx > 0) printf("rhw  %d  %7.5f     %6.4f  %6.4f  %6.4f\n",
                                   rhwvtx, drhw,
                                   rhwhite->vertices[rhwvtx].x,
                                   rhwhite->vertices[rhwvtx].y,
                                   rhwhite->vertices[rhwvtx].z);
          if (rhpvtx > 0) printf("rhp  %d  %7.5f     %6.4f  %6.4f  %6.4f\n",
                                   rhpvtx, drhp,
                                   rhpial->vertices[rhpvtx].x,
                                   rhpial->vertices[rhpvtx].y,
                                   rhpial->vertices[rhpvtx].z);
          printf("annot = %d, annotid = %d\n",annot,annotid);
          CTABprintASCII(lhwhite->ct,stdout);
	  MatrixFree(&CRS);
	  MatrixFree(&RAS);
          continue;
        }

	MatrixFree(&CRS);
	MatrixFree(&RAS);
        nctx++;
      } // slice
    } // row
  } // col
  MHT_maybeParallel_end();
  printf("nctx = %d\n",nctx);
  printf("Used brute-force search on %d voxels\n",nbrute);

  if (relabel_gca_name != NULL)    // reclassify voxels interior to white that are likely to be something else
  {
    MRI *mri_norm ;
    int nchanged = 0 ;

    printf("relabeling unlikely voxels in interior of white matter\n") ;

    mri_norm = MRIread(relabel_norm_name) ;
    if (mri_norm == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not load norm volume from %s\n", relabel_norm_name) ;

    {
      int        x, y, z, i ;
      MRI        *mri_tmp  = NULL, *mri_aseg_orig = NULL ;
      

      mri_aseg_orig = MRIcopy(ASeg, NULL) ;
      Ggca_x = Gx ; Ggca_y = Gy ; Ggca_z = Gz ; // diagnostics
      GCAregularizeCovariance(gca,1.0);   // don't use covariances for this classification
      GCAreclassifyUsingGibbsPriors(mri_norm, gca, ASeg, xform, 10, mri_fixed, GCA_RESTART_ONCE, NULL, 0.5, 0.5);
      for (i = 0 ; i < 2 ; i++)
      {
	mri_tmp = MRIcopy(ASeg, mri_tmp);
	  for (x = 0 ; x < mri_norm->width ; x++)
	    for (y = 0 ; y < mri_norm->height ; y++)
	      for (z = 0 ; z < mri_norm->depth ; z++)
	      {
		if (x  == Gx && y == Gy && z == Gz)
		  DiagBreak() ;
		if (MRIgetVoxVal(mri_fixed, x, y, z, 0) > 0)
		  continue ;
		if ((int)MRIgetVoxVal(ASeg, x, y, z, 0) > 0)  // only process voxels that are interior to the ribbon and unknown - shouldn't be
		  continue ;
		if ((MRIlabelsInNbhd(mri_tmp,  x,  y,  z, 1, Left_Lateral_Ventricle)  > 0 ) ||
		    (MRIlabelsInNbhd(mri_tmp,  x,  y,  z, 1, Right_Lateral_Ventricle)  > 0))   // neighbors a ventricular label
		{
		  GCA_NODE *gcan ;
		  GCA_PRIOR *gcap ;
		  int      xn, yn, zn, n, max_label = -1, label, xp, yp, zp ;
		  double   mah_dist, min_mah_dist, prior ;
		  float    vals[MAX_GCA_INPUTS] ;

#define REGION_WSIZE 3
		  if (x  == Gx && y == Gy && z == Gz)
		    DiagBreak() ;
		  GCAsourceVoxelToNode( gca, mri_norm, xform, x,  y, z, &xn, &yn, &zn) ;
		  gcan = GCAbuildRegionalGCAN(gca, xn, yn, zn, REGION_WSIZE);
		  GCAsourceVoxelToPrior( gca, mri_norm, xform, x,  y, z, &xp, &yp, &zp) ;
		  gcap = GCAbuildRegionalGCAP(gca, xp, yp, zp, REGION_WSIZE*gca->node_spacing/gca->prior_spacing-1);
		  load_vals(mri_norm, x, y, z, vals, mri_norm->nframes) ;
		  min_mah_dist = 1e10 ;
		  for (n = 0 ; n < gcan->nlabels ; n++)
		  {
		    label = gcan->labels[n] ;
		    if (IS_CORTEX(label) || label == 0 || IS_CEREBELLAR_GM(label) || IS_CEREBELLAR_WM(label))
		      continue ;   // prohibited interior to white
		    if ((MRIlabelsInNbhd6(mri_tmp,  x,  y,  z,  label)  == 0) ||
			(MRIlabelsInNbhd(mri_aseg_orig,  x,  y,  z, 2,  label) == 0))
		      continue ; // only if another voxel with this label exists nearby
		    mah_dist = GCAmahDist( &gcan->gcs[n], vals, mri_norm->nframes);
		    prior = getPrior(gcap, label) ;
		    if (mah_dist+log(prior) < min_mah_dist)
		    {
		      min_mah_dist = mah_dist+log(prior) ;
		      max_label = gcan->labels[n] ;
		    }
		  }
		  if (max_label < 0)
		  {
		    max_label = 0 ;
		    for (n = 0 ; n < gcan->nlabels ; n++)
		    {
		      label = gcan->labels[n] ;
		      if (IS_CORTEX(label) || label == 0 || IS_CEREBELLAR_GM(label) || IS_CEREBELLAR_WM(label))
			continue ;   // prohibited interior to white
		      mah_dist = GCAmahDist( &gcan->gcs[n], vals, mri_norm->nframes);
		      prior = getPrior(gcap, label) ;
		      if (mah_dist+log(prior) < min_mah_dist)
		      {
			min_mah_dist = mah_dist+log(prior) ;
			max_label = gcan->labels[n] ;
		      }
		    }
		  }

		  if (x  == Gx && y == Gy && z == Gz)
		    printf("reclassifying unknown voxel (%d, %d, %d) that neighbors ventricle %s --> %s (%d)\n", 
			   x, y, z, cma_label_to_name(MRIgetVoxVal(mri_tmp, x, y, z, 0)), cma_label_to_name(max_label), max_label) ;
		  MRIsetVoxVal(ASeg, x, y, z, 0, max_label) ;
		  GCAfreeRegionalGCAN(&gcan) ;
		}
	      }
      
	  if (Gdiag & DIAG_WRITE)
	  {
	    char fname[STRLEN], fonly[STRLEN] ;
	    FileNameRemoveExtension(OutASegFile, fonly) ;
	    int req = snprintf(fname, STRLEN, "%s.%3.3d.mgz", fonly, i) ;  
	    if( req >= STRLEN ) {
	      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	    }
	    printf("writing iter %d to %s\n", i, fname) ;
	    MRIwrite(ASeg, fname) ;
	  }
      }

      // expand into voxels that are adjacent to lots of voxels that are ventricle
      for (i = 2 ; i < 10 ; i++)
      {
	int xi, yi, zi, xk, yk, zk ;
	nchanged = 0 ;
	mri_tmp = MRIcopy(ASeg, mri_tmp);
	  for (x = 0 ; x < mri_norm->width ; x++)
	    for (y = 0 ; y < mri_norm->height ; y++)
	      for (z = 0 ; z < mri_norm->depth ; z++)
	      {
		if (x  == Gx && y == Gy && z == Gz)
		  DiagBreak() ;
		if (MRIgetVoxVal(mri_fixed, x, y, z, 0) > 0)
		  continue ;
		if ((int)MRIgetVoxVal(ASeg, x, y, z, 0) > 0)  // only process voxels that are interior to the ribbon and unknown - shouldn't be
		  continue ;
		if (MRIgetVoxVal(mri_dist, x, y, z, 0) > -2)
		  continue ;  // only if it is pretty far interior

		if ((MRIlabelsInNbhd(mri_tmp,  x,  y,  z, 1, Left_Lateral_Ventricle)  > 4 ) ||
		    (MRIlabelsInNbhd(mri_tmp,  x,  y,  z, 1, Right_Lateral_Ventricle)  > 4))   // neighbors a bunch of ventricular labels
		{
		  GCA_NODE *gcan ;
		  GCA_PRIOR *gcap ;
		  int      xn, yn, zn, n, max_label = -1, label, xp, yp, zp ;
		  double   mah_dist, min_mah_dist, prior ;
		  float    vals[MAX_GCA_INPUTS] ;

#define REGION_WSIZE 3
		  if (x  == Gx && y == Gy && z == Gz)
		    DiagBreak() ;
		  GCAsourceVoxelToNode( gca, mri_norm, xform, x,  y, z, &xn, &yn, &zn) ;
		  gcan = GCAbuildRegionalGCAN(gca, xn, yn, zn, REGION_WSIZE);
		  GCAsourceVoxelToPrior( gca, mri_norm, xform, x,  y, z, &xp, &yp, &zp) ;
		  gcap = GCAbuildRegionalGCAP(gca, xp, yp, zp, REGION_WSIZE*gca->node_spacing/gca->prior_spacing-1);
		  load_vals(mri_norm, x, y, z, vals, mri_norm->nframes) ;
		  min_mah_dist = 1e10 ;
		  for (n = 0 ; n < gcan->nlabels ; n++)
		  {
		    label = gcan->labels[n] ;
		    if (IS_CORTEX(label) || label == 0 || IS_CEREBELLAR_GM(label) || IS_CEREBELLAR_WM(label))
		      continue ;   // prohibited interior to white
		    mah_dist = GCAmahDist( &gcan->gcs[n], vals, mri_norm->nframes);
		    prior = getPrior(gcap, label) ;
		    if (mah_dist+log(prior) < min_mah_dist)
		    {
		      min_mah_dist = mah_dist+log(prior) ;
		      max_label = gcan->labels[n] ;
		    }
		  }
		  // if the max label is vent AND every neighboring vent is itself neighbored by lots of vent, relabel it
		  if (IS_VENTRICLE(max_label))  
		  {
		    int is_vent = 1, vlabels, olabel ;

		    for (xk = -1 ; xk <= 1 ; xk++)
		      for (yk = -1 ; yk <= 1 ; yk++)
			for (zk = -1 ; zk <= 1 ; zk++)
			{
			  xi = mri_tmp->xi[x+xk] ; yi = mri_tmp->yi[y+yk] ;  zi = mri_tmp->zi[z+zk] ;
			  olabel = MRIgetVoxVal(mri_tmp, xi, yi, zi, 0) ;
			  if (IS_VENTRICLE(olabel))
			  {
			    if (olabel == Right_Lateral_Ventricle)
			      vlabels = MRIlabelsInNbhd(mri_tmp,  xi,  yi,  zi, 1, Right_Lateral_Ventricle) ;
			    else
			      vlabels = MRIlabelsInNbhd(mri_tmp,  xi,  yi,  zi, 1, Left_Lateral_Ventricle) ;
			    if (vlabels < 9)
			      is_vent = 0 ;
			  }
			}
		    if (is_vent)
		    {
		      if (x  == Gx && y == Gy && z == Gz)
			printf("reclassifying unknown voxel (%d, %d, %d) that neighbors ventricle %s --> %s (%d)\n", 
			       x, y, z, cma_label_to_name(MRIgetVoxVal(mri_tmp, x, y, z, 0)), cma_label_to_name(max_label), max_label) ;
		      MRIsetVoxVal(ASeg, x, y, z, 0, max_label) ;
		      nchanged++ ;
		    }
		  }
		  GCAfreeRegionalGCAN(&gcan) ;
		}
	      }
	  
	  if (Gdiag & DIAG_WRITE)
	  {
	    char fname[STRLEN], fonly[STRLEN] ;
	    FileNameRemoveExtension(OutASegFile, fonly) ;
	    int req = snprintf(fname, STRLEN, "%s.%3.3d.mgz", fonly, i) ; 
	    if( req >= STRLEN ) {
	      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	    }
	    printf("writing iter %d to %s\n", i, fname) ;
	    MRIwrite(ASeg, fname) ;
	  }
	  printf("nchanged = %d\n", nchanged) ;
	  if (nchanged == 0)
	    break ;
      }

      // remove singleton voxels
      for (x = 0 ; x < mri_norm->width ; x++)
	for (y = 0 ; y < mri_norm->height ; y++)
	  for (z = 0 ; z < mri_norm->depth ; z++)
	  {
	    int label_orig, label_new ;
	    
	    if (x  == Gx && y == Gy && z == Gz)
	      DiagBreak() ;
	    if (MRIgetVoxVal(mri_fixed, x, y, z,0 ) > 0)
	      continue ;
	    
	    label_orig = MRIgetVoxVal(mri_aseg_orig, x, y, z, 0) ;
	    label_new = MRIgetVoxVal(ASeg, x, y, z, 0) ;
	    if (label_orig == label_new)
	      continue ;
	    if  (MRIlabelsInNbhd(ASeg,  x,  y,  z,  1, label_new) <= 1)
	    {
	      if (x  == Gx && y == Gy && z == Gz)
		printf("voxel(%d, %d, %d): reverting label back to %s (%d), was %s (%d)\n",
		       x, y, z, cma_label_to_name(label_orig), label_orig, cma_label_to_name(label_new), label_new) ;
	      MRIsetVoxVal(ASeg, x, y, z, 0, label_orig) ;
	    }
	  }


      MRIfree(&mri_tmp) ; MRIfree(&mri_aseg_orig);
    }

    MRIfree(&mri_norm) ; MRIfree(&mri_fixed) ; GCAfree(&gca) ; TransformFree(&xform) ;
  }
  if (mri_dist)
    MRIfree(&mri_dist) ;
  if (FixParaHipWM)
  {
    /* This is a bit of a hack. There are some vertices that have been
       ripped because they are "unkown". When the above alorithm finds
       these, it searches for the closest known vertex. If this is
       less than dmax away, then the wm voxel gets labeled
       accordingly.  However, there are often some voxels near
       ventralDC that are just close enough in 3d space to parahip to
       get labeled even though they are very far away along the
       surface. These voxels end up forming an island. CCSegment()
       will eliminate any islands. Unforunately, CCSegment() uses
       6-neighbor (face) definition of connectedness, so some voxels
       may be eliminated.
     */
    printf("Fixing Parahip LH WM\n");
    CCSegment(ASeg, 3016, Left_Unsegmented_WM); //3016 = lhphwm, 5001 = unsegmented WM left
    printf("Fixing Parahip RH WM\n");
    CCSegment(ASeg, 4016, Right_Unsegmented_WM); //4016 = rhphwm, 5002 = unsegmented WM right
  }

  // embed color lookup table
  if (!ASeg->ct) ASeg->ct = CTABreadDefault();
  if (!AParc->ct) AParc->ct = CTABreadDefault();

  printf("Writing output aseg to %s\n",OutASegFile);
  MRIwrite(ASeg,OutASegFile);

  if (OutAParcFile != NULL)
  {
    printf("Writing output aparc to %s\n",OutAParcFile);
    MRIwrite(AParc,OutAParcFile);
  }
  if (OutDistFile != NULL)
  {
    printf("Writing output dist file to %s\n",OutDistFile);
    MRIwrite(Dist,OutDistFile);
  }

  printf("#VMPC# mri_aparc2aseg VmPeak  %d\n",GetVmPeak());
  printf("mri_aparc2aseg done\n");

  return(0);
}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1)
  {
    usage_exit();
  }

  nargc   = argc;
  pargv = argv;
  while (nargc > 0)
  {
    option = pargv[0];
    if (debug)
    {
      printf("%d %s\n",nargc,option);
    }
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help")||
        !strcasecmp(option, "-h")||
        !strcasecmp(option, "--usage")||
        !strcasecmp(option, "-u")) print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug")) debug = 1;
    else if (!strcasecmp(option, "--lh")) {
      LHOnly = 1;
      DoLH = 1;
      RHOnly = 0;
      DoRH = 0;
    }
    else if (!strcasecmp(option, "--rh")) {
      LHOnly = 0;
      DoLH = 0;
      RHOnly = 1;
      DoRH = 1;
    }
    else if (!strcasecmp(option, "--relabel")){
      relabel_norm_name = pargv[0] ;
      relabel_xform_name = pargv[1] ;
      relabel_gca_name = pargv[2] ;
      relabel_label_intensities_name = pargv[3] ;
      printf("relabeling unlikely voxels interior to white matter surface:\n\tnorm: %s\n\t XFORM: %s\n\tGCA: %s\n\tlabel intensities: %s\n",
	     relabel_norm_name, relabel_xform_name, relabel_gca_name, relabel_label_intensities_name) ;
      nargsused = 4;
    }
    else if (!strcasecmp(option, "--no-relabel")){
      relabel_norm_name = NULL;
      relabel_xform_name =  NULL;
      relabel_gca_name =  NULL;
      relabel_label_intensities_name =  NULL;
    }
    else if (!strcasecmp(option, "--debug_voxel"))
    {
      if (nargc < 3)
      {
        argnerr(option,3);
      }
      Gx = atoi(pargv[0]) ;
      Gy = atoi(pargv[1]) ;
      Gz = atoi(pargv[2]) ;
      nargsused = 3;
    }
    else if (!strcasecmp(option, "--smooth_normals"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      normal_smoothing_iterations = atoi(pargv[0]) ;
      nargsused = 1 ;
    }
    // This was --ribbon, but changed to --old-ribbon 4/17/08 DNG
    else if (!strcasecmp(option, "--old-ribbon"))
    {
      UseRibbon = 1;
      UseNewRibbon = 0;
    }
    else if (!strcasecmp(option, "--volmask") ||
             !strcasecmp(option, "--new-ribbon"))
    {
      UseNewRibbon = 1;
    }
    else if (!strcasecmp(option, "--noribbon"))
    {
      UseRibbon = 0;
      UseNewRibbon = 0;
    }
    else if (!strcasecmp(option, "--labelwm"))
    {
      LabelWM = 1;
    }
    else if (!strcasecmp(option, "--fix-parahipwm"))
    {
      FixParaHipWM = 1;
    }
    else if (!strcasecmp(option, "--no-fix-parahipwm"))
    {
      FixParaHipWM = 0;
    }
    else if (!strcasecmp(option, "--hypo-as-wm"))
    {
      LabelHypoAsWM = 1;
    }
    else if (!strcasecmp(option, "--rip-unknown"))
    {
      RipUnknown = 1;
    }
    else if (!strcasecmp(option, "--no-hash"))
    {
      UseHash = 0;
    }
    else if (!strcmp(option, "--sd"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      SUBJECTS_DIR = pargv[0];
      setenv("SUBJECTS_DIR",SUBJECTS_DIR,1);
      nargsused = 1;
    }
    else if (!strcmp(option, "--s"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      subject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--oaseg") || !strcmp(option, "--o"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      OutASegFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--a2005s"))
    {
      annotname = "aparc.a2005s";
      baseoffset = 100;
    }
    else if (!strcmp(option, "--a2009s"))
    {
      annotname = "aparc.a2009s";
      baseoffset = 10100;
    }
    else if (!strcmp(option, "--base-offset")){
      if (nargc < 1)argnerr(option,1);
      sscanf(pargv[0],"%d",&baseoffset);
      nargsused = 1;
    }
   else if (!strcmp(option, "--aseg"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      asegname = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--annot"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      annotname = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--annot-table"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      // annotation_table_file is declared in annotation.h
      // default is $FREESURFER_HOME/Simple_surface_labels2009.txt
      annotation_table_file = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--oaparc"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      OutAParcFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--ctxseg"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      CtxSegFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--dist"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      OutDistFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--hashres"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%f",&hashres);
      nargsused = 1;
    }
    else if (!strcmp(option, "--wmparc-dmax"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      sscanf(pargv[0],"%f",&dmaxctx);
      nargsused = 1;
    }
    else if (!strcmp(option, "--crs-test"))
    {
      if (nargc < 3)
      {
        argnerr(option,3);
      }
      sscanf(pargv[0],"%d",&ctest);
      sscanf(pargv[1],"%d",&rtest);
      sscanf(pargv[2],"%d",&stest);
      crsTest = 1;
      nargsused = 3;
    }
    else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nthreads);
      #ifdef _OPENMP
      omp_set_num_threads(nthreads);
      #endif
      nargsused = 1;
    } 
    else
    {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option))
      {
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      }
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
#include "mri_aparc2aseg.help.xml.h"
static void print_usage(void)
{
  outputHelpXml(mri_aparc2aseg_help_xml,mri_aparc2aseg_help_xml_len);
}
/* --------------------------------------------- */
static void print_help(void)
{
  outputHelpXml(mri_aparc2aseg_help_xml,mri_aparc2aseg_help_xml_len);
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if (n==1)
  {
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  }
  else
  {
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  }
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void)
{
  if (subject == NULL)
  {
    printf("ERROR: must specify a subject\n");
    exit(1);
  }
  if (OutASegFile == NULL)
  {
    sprintf(tmpstr,"%s/%s/mri/%s+aseg.mgz",SUBJECTS_DIR,subject,annotname);
    OutASegFile = strcpyalloc(tmpstr);
  }
  if (UseRibbon && UseNewRibbon)
  {
    printf("ERROR: cannot --old-ribbon and --new-ribbon\n");
    exit(1);
  }
  if (CtxSegFile && ! RipUnknown)
  {
    printf("ERROR: can only use --ctxseg with --rip-unknown\n");
    exit(1);
  }
  if (CtxSegFile)
  {
    if (!fio_FileExistsReadable(CtxSegFile))
    {
      sprintf(tmpstr,"%s/%s/mri/%s",SUBJECTS_DIR,subject,CtxSegFile);
      if (! fio_FileExistsReadable(tmpstr))
      {
        printf("ERROR: cannot find %s or %s\n",CtxSegFile,tmpstr);
        exit(1);
      }
      CtxSegFile = strcpyalloc(tmpstr);
    }
  }
  if (FixParaHipWM && ! LabelWM)
  {
    FixParaHipWM  = 0;
  }

  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"SUBJECTS_DIR %s\n",SUBJECTS_DIR);
  fprintf(fp,"subject %s\n",subject);
  fprintf(fp,"outvol %s\n",OutASegFile);
  fprintf(fp,"useribbon %d\n",UseRibbon);
  fprintf(fp,"baseoffset %d\n",baseoffset);
  if (LabelWM)
  {
    printf("labeling wm\n");
    if (LabelHypoAsWM)
    {
      printf("labeling hypo-intensities as wm\n");
    }
    printf("dmaxctx %f\n",dmaxctx);
  }
  fprintf(fp,"RipUnknown %d\n",RipUnknown);
  if (CtxSegFile)
  {
    fprintf(fp,"CtxSeg %s\n",CtxSegFile);
  }
  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if (len < 2)
  {
    return(0);
  }

  if (flag[0] == '-' && flag[1] != '-')
  {
    return(1);
  }
  return(0);
}

/*---------------------------------------------------------------*/
int FindClosestLRWPVertexNo(int c, int r, int s,
                            int *lhwvtx, int *lhpvtx,
                            int *rhwvtx, int *rhpvtx,
                            MATRIX *Vox2RAS,
                            MRIS *lhwite,  MRIS *lhpial,
                            MRIS *rhwhite, MRIS *rhpial,
                            MHT *lhwhite_hash, MHT *lhpial_hash,
                            MHT *rhwhite_hash, MHT *rhpial_hash)
{
  static MATRIX *CRS = NULL;
  static MATRIX *RAS = NULL;
  static struct { float x,y,z; } vtx;
  static float dlhw, dlhp, drhw, drhp,dmin;
  int annot, hemi, annotid;

  if (CRS == NULL)
  {
    CRS = MatrixAlloc(4,1,MATRIX_REAL);
    CRS->rptr[4][1] = 1;
    RAS = MatrixAlloc(4,1,MATRIX_REAL);
    RAS->rptr[4][1] = 1;
  }

  CRS->rptr[1][1] = c;
  CRS->rptr[2][1] = r;
  CRS->rptr[3][1] = s;
  RAS = MatrixMultiply(Vox2RAS,CRS,RAS);
  vtx.x = RAS->rptr[1][1];
  vtx.y = RAS->rptr[2][1];
  vtx.z = RAS->rptr[3][1];

  *lhwvtx = MHTfindClosestVertexNoXYZ(lhwhite_hash,lhwhite,vtx.x,vtx.y,vtx.z,&dlhw);
  *lhpvtx = MHTfindClosestVertexNoXYZ(lhpial_hash, lhpial, vtx.x,vtx.y,vtx.z,&dlhp);
  *rhwvtx = MHTfindClosestVertexNoXYZ(rhwhite_hash,rhwhite,vtx.x,vtx.y,vtx.z,&drhw);
  *rhpvtx = MHTfindClosestVertexNoXYZ(rhpial_hash, rhpial, vtx.x,vtx.y,vtx.z,&drhp);

  printf("lh white: %d %g\n",*lhwvtx,dlhw);
  printf("lh pial:  %d %g\n",*lhpvtx,dlhp);
  printf("rh white: %d %g\n",*rhwvtx,drhw);
  printf("rh pial:  %d %g\n",*rhpvtx,drhp);

  hemi = 0;
  dmin = -1;
  if (*lhwvtx < 0 && *lhpvtx < 0 && *rhwvtx < 0 && *rhpvtx < 0)
  {
    printf("ERROR2: could not map to any surface.\n");
    printf("crs = %d %d %d, ras = %6.4f %6.4f %6.4f \n",
           c,r,s,vtx.x,vtx.y,vtx.z);
    printf("Using Bruce Force\n");
    *lhwvtx = MRISfindClosestVertex(lhwhite,vtx.x,vtx.y,vtx.z,&dlhw, CURRENT_VERTICES);
    *lhpvtx = MRISfindClosestVertex(lhpial, vtx.x,vtx.y,vtx.z,&dlhp, CURRENT_VERTICES);
    *rhwvtx = MRISfindClosestVertex(rhwhite,vtx.x,vtx.y,vtx.z,&drhw, CURRENT_VERTICES);
    *rhpvtx = MRISfindClosestVertex(rhpial, vtx.x,vtx.y,vtx.z,&drhp, CURRENT_VERTICES);
    printf("lh white: %d %g\n",*lhwvtx,dlhw);
    printf("lh pial:  %d %g\n",*lhpvtx,dlhp);
    printf("rh white: %d %g\n",*rhwvtx,drhw);
    printf("rh pial:  %d %g\n",*rhpvtx,drhp);
    return(1);
  }
  if (dlhw <= dlhp && dlhw < drhw && dlhw < drhp && (*lhwvtx >= 0))
  {
    annot = lhwhite->vertices[*lhwvtx].annotation;
    hemi = 1;
    if (lhwhite->ct)
    {
      CTABfindAnnotation(lhwhite->ct, annot, &annotid);
    }
    else
    {
      annotid = annotation_to_index(annot);
    }
    dmin = dlhw;
  }
  if (dlhp < dlhw && dlhp < drhw && dlhp < drhp && (*lhpvtx >= 0))
  {
    annot = lhwhite->vertices[*lhpvtx].annotation;
    hemi = 1;
    if (lhwhite->ct)
    {
      CTABfindAnnotation(lhwhite->ct, annot, &annotid);
    }
    else
    {
      annotid = annotation_to_index(annot);
    }
    dmin = dlhp;
  }

  if (drhw < dlhp && drhw < dlhw && drhw <= drhp && (*rhwvtx >= 0))
  {
    annot = rhwhite->vertices[*rhwvtx].annotation;
    hemi = 2;
    if (rhwhite->ct)
    {
      CTABfindAnnotation(lhwhite->ct, annot, &annotid);
    }
    else
    {
      annotid = annotation_to_index(annot);
    }
    dmin = drhw;
  }
  if (drhp < dlhp && drhp < drhw && drhp < dlhw && (*rhpvtx >= 0))
  {
    annot = rhwhite->vertices[*rhpvtx].annotation;
    hemi = 2;
    if (rhwhite->ct)
    {
      CTABfindAnnotation(lhwhite->ct, annot, &annotid);
    }
    else
    {
      annotid = annotation_to_index(annot);
    }
    dmin = drhp;
  }

  printf("hemi = %d, annotid = %d, dist = %g\n",hemi,annotid,dmin);

  return(0);
}

/*!
  \fn int CCSegment(MRI *seg, int segid, int segidunknown)
  Constraints a sementation ID to consist of voxels that
  are spatially contiguous (6 face neighbors, not edge
  or corner). The voxels in the largest cluster are not
  changed. The voxels in the other clusters are set to
  segidunknown.
*/
int CCSegment(MRI *seg, int segid, int segidunknown)
{
  MRI_SEGMENTATION *sgmnt;
  int k,kmax,index,c,r,s;

  sgmnt = MRIsegment(seg,segid-.5,segid+.5);
  printf("  Found %d clusters\n",sgmnt->nsegments);

  kmax = 0;
  for (k=0; k < sgmnt->nsegments; k++)
    if (sgmnt->segments[k].nvoxels > sgmnt->segments[kmax].nvoxels)
    {
      kmax = k;
    }

  for (k=0; k < sgmnt->nsegments; k++)
  {
    printf("     %d k %f\n",k,sgmnt->segments[k].area);
    if (k==kmax)
    {
      continue;
    }
    for (index = 0; index < sgmnt->segments[k].nvoxels; index++)
    {
      c = sgmnt->segments[k].voxels[index].x;
      r = sgmnt->segments[k].voxels[index].y;
      s = sgmnt->segments[k].voxels[index].z;
      MRIsetVoxVal(seg,c,r,s,0,segidunknown);
    }
  }
  MRIsegmentFree(&sgmnt);
  return(0);
}

/* BRF added some checks here to make sure closest vertex (usually
   pial but can be white) isn't on the other bank of a sulcus or
   through a thin white matter strand. This removes inaccurate voxels
   that used to speckle the aparc+aseg */
double BRFdotCheck(MRIS *surf, int vtxno, int c, int r, int s, MRI *AParc)
{
  double dx, dy, dz, nx, ny, nz, xv, yv, zv, x1, y1, z1, dot ;
  VERTEX *v ;
  v = &surf->vertices[vtxno] ;
  MRISvertexToVoxel(surf, v, AParc, &xv, &yv, &zv) ;
  x1 = v->x + v->nx ;  y1 = v->y + v->ny ; z1 = v->z + v->nz ; 
  MRISsurfaceRASToVoxel(surf, AParc, x1, y1, z1, &nx, &ny, &nz) ;
  nx -= xv ; ny -= yv ; nz -= zv ;  // normal in voxel coords
  dx = c-xv ; dy = r-yv ; dz = s-zv ;
  dot = dx*nx + dy*ny + dz*nz ;
  return(dot);
}


    
  
