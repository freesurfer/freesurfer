/**
 * @brief Computes the desired target location of a surface, mainly for exploring target placement strategies
 *
 */
/*
 * Original Author: Douglas N Greve
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
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <float.h>
#include <errno.h>

#include "macros.h"
#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "cmdargs.h"
#include "dmatrix.h"
#include "matrix.h"
#include "cma.h"
#include "romp_support.h"
#include "version.h"

#undef private
int main(int argc, char *argv[]) ;
const char *Progname = "mris_target_pos";
static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int debug = 0, checkoptsonly = 0;

// A Simple class to sort floats and keep an index
class FloatInt {
public:
  float f;
  int i;
  static int compare(const void *v1, const void *v2)
  {
    FloatInt i1, i2;
    i1 = *((FloatInt *)v1);
    i2 = *((FloatInt *)v2);
    if (fabs(i1.f) < fabs(i2.f)) return (-1);
    if (fabs(i1.f) > fabs(i2.f)) return (+1);
    return(0);  // equal
  }
  int sort(float *flist, int *ilist, int nlist)
  {
    FloatInt *filist = (FloatInt*)calloc(sizeof(FloatInt),nlist);
    int n;
    for(n=0; n < nlist; n++) {
      filist[n].f = flist[n];
      filist[n].i = ilist[n];
    }
    qsort(filist, nlist, sizeof(FloatInt), FloatInt::compare);
    for(n=0; n < nlist; n++) {
      flist[n] = filist[n].f;
      ilist[n]=filist[n].i;
    }
    free(filist);
    filist=NULL;
    return(0);
  }
};

class MRIS_TARGET_POS
{
public:
  MRI *vol=NULL; // intensity volume to place surfaces on
  MRI *seg=NULL; // segmentation to stop from going across hemi
  MRIS *surf=NULL; // surface
  MRI  *proj=NULL; // surf intensity projection along the normal
  MRI  *projsm=NULL; // surf intensity projection along the normal, smoothed
  MRI  *projsmgrad=NULL; // gradient of smoothed surf intensity projection along the normal
  MRI *mask=NULL;
  double MinSampleDist,MaxSampleDist,DeltaSampleDist,*SampleDist=NULL;
  int NSampleDist=0,nSampleDist0;
  double MinSearchDist,MaxSearchDist;
  int interptype=SAMPLE_TRILINEAR;
  double Contrast = -1; // -1 for T1 +1 for T2
  double sigma; // guass std to smoothing along projetion, in mm
  double *sigmalist;
  int nsigma;
  MATRIX **Glist;
  double InwardThresh=0, OutwardThresh=0;
  MRIS_SurfRAS2VoxelMap* sras2v_map;
  AutoDetGWStats adgws;
  int InitSampleDist(void);
  int alloc(void);
  int GetVertexTarget(const int vtxno);
};

int MRIS_TARGET_POS::InitSampleDist(void)
{
  NSampleDist = round((MaxSampleDist-MinSampleDist)/DeltaSampleDist) + 1;
  double MaxSampleDist2 = MinSampleDist + (NSampleDist-1)*DeltaSampleDist;
  if(MaxSampleDist2 != MaxSampleDist){
    printf("INFO: resetting MaxSampleDist from %g to %g\n",MaxSampleDist,MaxSampleDist2);
    MaxSampleDist = MaxSampleDist2;
  }
  SampleDist = (double *) calloc(sizeof(double),NSampleDist);
  int n;
  nSampleDist0 = 0;
  double minSD = 10e10;
  for(n=0; n < NSampleDist;n++){
    SampleDist[n] = MinSampleDist + n*DeltaSampleDist;
    //printf("%3d %6.4f\n",n,SampleDist[n]);
    if(fabs(SampleDist[n])<minSD){
      minSD = fabs(SampleDist[n]);
      nSampleDist0 = n;
    }
  }
  printf("%d %g %g %g %g %d\n",NSampleDist,MinSampleDist,MaxSampleDist,DeltaSampleDist,minSD,nSampleDist0);
  return(0);
}

int MRIS_TARGET_POS::alloc(void)
{
  proj = MRIallocSequence(surf->nvertices,1,1,MRI_FLOAT,NSampleDist);
  MRIcopyHeader(vol,proj);
  MRIcopyPulseParameters(vol,proj);

  projsm = MRIallocSequence(surf->nvertices,1,1,MRI_FLOAT,NSampleDist);
  MRIcopyHeader(vol,projsm);
  MRIcopyPulseParameters(vol,projsm);

  projsmgrad = MRIallocSequence(surf->nvertices,1,1,MRI_FLOAT,NSampleDist);
  MRIcopyHeader(vol,projsmgrad);
  MRIcopyPulseParameters(vol,projsmgrad);

  sras2v_map = MRIS_makeRAS2VoxelMap(vol, surf);

  MRISsaveVertexPositions(surf, ORIGINAL_VERTICES) ; 

  nsigma = 4;
  sigmalist = (double*)calloc(sizeof(double),nsigma);
  Glist = (MATRIX**) calloc(sizeof(MATRIX*),nsigma);
  int n;
  sigmalist[0] = sigma;
  for(n=1; n < nsigma; n++)  sigmalist[n] = 2*sigmalist[n-1];
  for(n=0; n < nsigma; n++) {
    MATRIX *G = GaussianMatrix(NSampleDist, sigmalist[n]/DeltaSampleDist, 0, NULL);
    // Force the sum of each row to be 1
    int i;
    for(i=0; i < NSampleDist; i++){
      int j;
      double gsum=0;
      for(j=0; j < NSampleDist; j++)gsum = gsum + G->rptr[i+1][j+1];
      for(j=0; j < NSampleDist; j++) G->rptr[i+1][j+1] /= gsum;
    }
    Glist[n] = G;
  }
  return(0);
}

int MRIS_TARGET_POS::GetVertexTarget(const int vtxno)
{
  VERTEX *v = &(surf->vertices[vtxno]);

  // No matter what, set the target to itself. Set value as well?
  v->targx = v->x;
  v->targy = v->y;
  v->targz = v->z;
  v->d = 0;
  v->val = 0; // should prob be set
  v->marked = 0; // assume the worst

  if(v->ripflag) return(1);

  // Get intensity projection
  MATRIX *ip = MatrixAlloc(NSampleDist,1,MATRIX_REAL);
  MATRIX *ipvalid = MatrixAlloc(NSampleDist,1,MATRIX_REAL);
  int n;
  for(n=0; n < NSampleDist; n++){
    double x,y,z, c,r,s, val;
    // project xyz to this distance
    x = v->x + v->nx * SampleDist[n];
    y = v->y + v->ny * SampleDist[n];
    z = v->z + v->nz * SampleDist[n];
    // Compute the col, row, slice at this point
    MRIS_useRAS2VoxelMap(sras2v_map, vol, x, y, z, &c, &r, &s);
    // Skip if out of the volume
    if(MRIindexNotInVolume(vol, c,r,s)) continue;
    // Skip if not in the mask
    if(mask){
       double const mval = MRIgetVoxVal(mask, nint(c), nint(r), nint(s), 0);
       if(mval < 0.5) continue;
    }
    // Skip if in the other hemi
    if(seg){
       double const label = MRIgetVoxVal(seg, nint(c), nint(r), nint(s), 0);
       if((surf->hemisphere == LEFT_HEMISPHERE  && IS_RH_CLASS(label)) ||
	  (surf->hemisphere == RIGHT_HEMISPHERE && IS_LH_CLASS(label))) {
	 continue;
       }
    }
    // OK, now sample the value at this distance
    MRIsampleVolumeType(vol, c, r, s, &val, interptype);
    ip->rptr[n+1][1] = val;
    MRIsetVoxVal(proj,vtxno,0,0,n,val);
    ipvalid->rptr[n+1][1] = 1; // set as valid
  }

  /* Below does not match CBV. It should go through each sigma looking for
     a valid bracket. If none is found, then it should increase sigma.
     Once it has a valid bracket, then it looks for the gradmax.
   */
  /* Go through each smoothing level starting at low smoothing. If a
     valid maxgrad is found, then stop*/
  MATRIX *G=NULL, *ipsm=NULL, *ipsmgrad=NULL;
  double tsigma=-1; 
  int nOuterLimit = -1;
  int nInnerLimit = -1;
  int ValidBracket = 0;
  int nthsigma;
  for (nthsigma = 0; nthsigma < nsigma; nthsigma++){
    
    // Smooth the projection (what about invalid points?)
    tsigma = sigmalist[nthsigma];
    G = Glist[nthsigma];
    ipsm = MatrixMultiplyD(G,ip,ipsm);

    // Compute the gradient at each distance
    if(ipsmgrad == NULL) ipsmgrad = MatrixAlloc(NSampleDist,1,MATRIX_REAL);
    for(n=0; n < NSampleDist; n++){
      int k = n + 1;
      int kp1 = k + 1;
      int km1 = k - 1;
      double scale = (2*DeltaSampleDist);
      if(n==0)	           {km1 = k; scale = DeltaSampleDist;}
      if(n==NSampleDist-1) {kp1 = k; scale = DeltaSampleDist;}
      ipsmgrad->rptr[k][1] = Contrast*(ipsm->rptr[kp1][1] - ipsm->rptr[km1][1])/scale;
    }

    if(Gdiag_no == vtxno && nthsigma == 0){
      MatrixWriteTxt("ipsmgrad.mtx",ipsmgrad);
    }

    // Look outward
    int breakcode = 0;
    for(n=nSampleDist0; n < NSampleDist-1; n++){
      int k = n + 1;
      double const val = ip->rptr[k][1];
      double const grad = ipsmgrad->rptr[k][1];
      if(SampleDist[n] > MaxSearchDist) breakcode = 1;
      else if(Contrast < 0 && val < OutwardThresh) breakcode = 2;
      else if(Contrast > 0 && val > OutwardThresh) breakcode = 3;
      else if(grad < 0) breakcode = 4;
      else if(!ipvalid->rptr[k][1]) breakcode = 5;
      if(Gdiag_no == vtxno){
	printf("%g %g %g\n",SampleDist[n],val,grad);
	fflush(stdout);
      }
      if(breakcode){
	v->marked3 = breakcode;
	if(Gdiag_no == vtxno){
	  printf("vno=%d s=%4.1f breaking outward loop code=%d, d=%g n=%d val=%g thresh=%g grad=%g\n",
		 vtxno,tsigma,breakcode,SampleDist[n],n,val,OutwardThresh,grad);
	  fflush(stdout);
	}
	break;
      }
    }
    nOuterLimit = n-1;

    // Now look inward
    breakcode = 0;
    for(n=nSampleDist0; n > 0; n--){
      int k = n + 1;
      double const val = ip->rptr[k][1];
      double const grad = ipsmgrad->rptr[k][1];
      if(SampleDist[n] < MinSearchDist) breakcode = 1;
      else if(Contrast < 0 && val > InwardThresh) breakcode = 2;
      else if(Contrast > 0 && val < InwardThresh) breakcode = 3;
      else if(grad < 0) breakcode = 4;
      else if(!ipvalid->rptr[k][1]) breakcode = 4;
      if(breakcode){
	v->undefval = breakcode;
	if(Gdiag_no == vtxno){
	  printf("vno=%d s=%4.1f breaking inward loop code=%d, d=%g n=%d val=%g thresh=%g grad=%g\n",
		 vtxno,tsigma,breakcode,SampleDist[n],n,val,InwardThresh,grad);
	  fflush(stdout);
	}
	break;
      }
    }
    nInnerLimit = n+1;

    if(nOuterLimit > nInnerLimit){
      // Found a valid bracket
      ValidBracket = 1;
      break;
    }
    if(Gdiag_no == vtxno){
      printf("vno=%d s=%4.1f could not find a valid bracket for this sigma, increasing\n",vtxno,tsigma);
      fflush(stdout);
    }
  } // sigma loop

  // Load smoothed projections and gradient into an MRI struct
  for(n=0; n < NSampleDist; n++) {
    MRIsetVoxVal(projsm,vtxno,0,0,n,ipsm->rptr[n+1][1]);
    MRIsetVoxVal(projsmgrad,vtxno,0,0,n,ipsmgrad->rptr[n+1][1]);
  }

  if(Gdiag_no == vtxno){
    printf("vno=%d ValidBracket=%d   %d %d  ",vtxno,ValidBracket,nInnerLimit,nOuterLimit);
    if(ValidBracket) printf("%g %g",SampleDist[nInnerLimit],SampleDist[nOuterLimit]);
    printf("\n");
    fflush(stdout);
  }

  // Now, at this sigma, fine the peaks
  double maxgrad=-1;
  int nmaxgrad=-1;
  int nlocalpeaks = 0;
  if(ValidBracket){
    // Find best maximum in the gradient within the bracket
    for(n=nInnerLimit; n <= nOuterLimit; n++){
      int k = n + 1;
      double const grad = ipsmgrad->rptr[k][1];
      if(grad > ipsmgrad->rptr[k-1][1] && grad > ipsmgrad->rptr[k+1][1]){
	// This is a local peak
	// Need to account for plateaus
	nlocalpeaks++;
	if(Gdiag_no == vtxno){
	  printf("vno=%d Found local peak at %d %g grad=%g\n",vtxno,n,SampleDist[n],grad);
	  fflush(stdout);
	}
	if(maxgrad < grad){
	  // And this is the best peak so far
	  maxgrad = grad;
	  nmaxgrad = n;
	}
      }
    }
  }

  v->val2 = tsigma; // smoothing level along gradient used to find the target
  v->mean = maxgrad;     // derivative at target intensity
  v->marked2 = nlocalpeaks; 
  if(nmaxgrad >= 0) {
    v->marked = 1;         // vertex has good data
    v->val  = ip->rptr[nmaxgrad+1][1]; // target intensity
    v->d    = SampleDist[nmaxgrad];   // dist to target intensity along normal
    v->targx = v->x + v->nx * SampleDist[nmaxgrad];
    v->targy = v->y + v->ny * SampleDist[nmaxgrad];
    v->targz = v->z + v->nz * SampleDist[nmaxgrad];
  }
  else {
    if(Gdiag_no == vtxno){
      printf("vno=%d could not find a valid peak\n",vtxno);
      fflush(stdout);
    }
  }

  if(Gdiag_no == vtxno){
    printf("vno=%d marked %d  d=%g sigma=%g val=%g grad=%g nlocal=%d nI=%d nO=%d\n",
	   vtxno,v->marked,v->d,tsigma,v->val,v->mean,nlocalpeaks,nInnerLimit,nOuterLimit);
    MatrixWriteTxt("ip.dat", ip);
    fflush(stdout);
  }

  MatrixFree(&ip);
  MatrixFree(&ipsm);
  MatrixFree(&ipsmgrad);
  MatrixFree(&ipvalid);
  return(0);
}

const char *insurfname=NULL, *outsurfname=NULL, *involname=NULL;
const char *adgwsinfile=NULL;
const char *format = "mgh";
const char *interpmethod = "trilinear";
char tmpstr[2000];
const char *dumpdir = NULL;
struct utsname uts;
char *cmdline, cwd[2000];
int   sinchw;
MRIS_TARGET_POS mtp;
int DoCBV = 0;
const char *riplabelfile=NULL;
int nthreads = 1;
int SavePointSet(const char *fname, MRIS *surf, int npoints);
int npointset=50;

/*--------------------------------------------------*/
int main(int argc, char **argv) {
  FILE *fp;
  int nargs;

  mtp.MinSampleDist = -4;
  mtp.MaxSampleDist = +6;
  mtp.DeltaSampleDist = 0.1;
  mtp.MinSearchDist = -3.0;
  mtp.MaxSearchDist = +6.0;
  //mtp.interptype=SAMPLE_TRILINEAR;
  mtp.sigma=0.2; 

  nargs = handleVersionOption(argc, argv, "mris_target_pos");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);
#ifdef _OPENMP
  printf("%d avail.processors, using %d\n",omp_get_num_procs(),omp_get_max_threads());
#endif

  mtp.vol = MRIread(involname);
  mtp.surf = MRISread(insurfname);
  mtp.InitSampleDist();
  mtp.alloc();
  printf("interp %d\n",mtp.interptype);

  if(riplabelfile){
    LABEL *riplabel = LabelRead("",riplabelfile);
    if(riplabel == NULL) exit(1);
    printf("Ripping based on label %s\n",riplabelfile);
    MRISripNotLabel(mtp.surf,riplabel);
  }

  const char *field;

  MRISclearMarks(mtp.surf); /* for soap bubble smoothing later */

  MRIS* surftarget = MRISclone(mtp.surf);
  double border_hi = mtp.adgws.white_border_hi;
  double outside_low = mtp.adgws.white_outside_low;

  MRISclearMarks(mtp.surf); /* for soap bubble smoothing later */

  if(DoCBV){
    double inside_hi = mtp.adgws.white_inside_hi;
    double outside_hi = mtp.adgws.white_outside_hi;
    double border_low = 0;
    double sigma = 1.0;
    
    MRIScomputeBorderValues(mtp.surf, mtp.vol, NULL, inside_hi,border_hi,border_low,outside_low,outside_hi,
			    sigma, 10, NULL, GRAY_WHITE, NULL, 0, 0, NULL,-1,-1) ;
    
  }
  else {
    mtp.InwardThresh = border_hi;
    mtp.OutwardThresh = outside_low;
    printf("InwardThresh  = %g\n",mtp.InwardThresh);
    printf("OutwardThresh = %g\n",mtp.OutwardThresh);
    
    Timer mytimer ;
    int msec;
    mytimer.reset();
    int vtxno;
#ifdef HAVE_OPENMP
#pragma omp parallel for 
#endif
    for(vtxno=0; vtxno < mtp.surf->nvertices; vtxno++){
      mtp.GetVertexTarget(vtxno);
    }
    msec = mytimer.milliseconds() ;
    printf("TargetPos finished in %6.4f min\n",(float)msec/(60*1000.0f)); fflush(stdout);
  }

  double drms = 0, dmean=0, dl1=0, dmax=0;
  int nnotrip = 0, ndmax=0;
  for(int v=0; v < mtp.surf->nvertices;v++){
    surftarget->vertices[v].x = mtp.surf->vertices[v].targx;	
    surftarget->vertices[v].y = mtp.surf->vertices[v].targy;	
    surftarget->vertices[v].z = mtp.surf->vertices[v].targz;	
    if(!mtp.surf->vertices[v].ripflag) {
      drms += (mtp.surf->vertices[v].d*mtp.surf->vertices[v].d);
      dmean += mtp.surf->vertices[v].d;
      dl1 += fabs(mtp.surf->vertices[v].d);
      if(fabs(dmax) < fabs(mtp.surf->vertices[v].d)) {
        dmax = mtp.surf->vertices[v].d;
        ndmax = v;
      }
      nnotrip ++;
    }
  }
  MRISwrite(surftarget,outsurfname);
  drms = sqrt(drms)/nnotrip;
  dmean = dmean/nnotrip;
  dl1 = dl1/nnotrip;
  printf("RMS %g  L1 %g Mean %g Max %g %d n=%d\n",drms,dl1,dmean,dmax,ndmax,nnotrip);

  if(dumpdir){
    int err = mkdir(dumpdir,0777);
    if(err != 0 && errno != EEXIST) {
      printf("ERROR: creating directory %s\n",dumpdir);
      perror(NULL);
      return(1);
    }
    sprintf(tmpstr,"%s/marked.%s",dumpdir,format);
    field = "marked"; MRISwriteField(mtp.surf, &field, 1, tmpstr);
    sprintf(tmpstr,"%s/sigma.%s",dumpdir,format);
    field = "val2"; MRISwriteField(mtp.surf, &field, 1, tmpstr);
    sprintf(tmpstr,"%s/val.%s",dumpdir,format);
    field = "val"; MRISwriteField(mtp.surf, &field, 1, tmpstr);
    sprintf(tmpstr,"%s/opt.dist.%s",dumpdir,format);
    field = "d"; MRISwriteField(mtp.surf, &field, 1, tmpstr);
    sprintf(tmpstr,"%s/maxgrad.%s",dumpdir,format);
    field = "mean"; MRISwriteField(mtp.surf, &field, 1, tmpstr);
    sprintf(tmpstr,"%s/rms.dat",dumpdir);
    fp = fopen(tmpstr,"w");
    fprintf(fp,"%g %g %g %g\n",drms,dmean,dl1,dmax);
    fclose(fp);
    sprintf(tmpstr,"%s/dist.ps",dumpdir);
    SavePointSet(tmpstr, mtp.surf, npointset);
    if(DoCBV == 0){
      sprintf(tmpstr,"%s/thresh.dat",dumpdir);
      fp = fopen(tmpstr,"w");
      fprintf(fp,"%g %g %d %g %g\n",mtp.InwardThresh,mtp.OutwardThresh,
	      mtp.NSampleDist,mtp.MinSampleDist,mtp.DeltaSampleDist);
      fclose(fp);
      sprintf(tmpstr,"%s/nlocalpeaks.%s",dumpdir,format);
      field = "marked2"; MRISwriteField(mtp.surf, &field, 1, tmpstr);
      sprintf(tmpstr,"%s/outer-break.%s",dumpdir,format);
      field = "marked3"; MRISwriteField(mtp.surf, &field, 1, tmpstr);
      sprintf(tmpstr,"%s/inner-break.%s",dumpdir,format);
      field = "undefval"; MRISwriteField(mtp.surf, &field, 1, tmpstr);
      sprintf(tmpstr,"%s/proj.%s",dumpdir,format);
      MRIwrite(mtp.proj,tmpstr);
      sprintf(tmpstr,"%s/projsm.%s",dumpdir,format);
      MRIwrite(mtp.projsm,tmpstr);
      sprintf(tmpstr,"%s/projsmgrad.%s",dumpdir,format);
      MRIwrite(mtp.projsmgrad,tmpstr);
      for(int n=0; n < mtp.nsigma; n++){
        sprintf(tmpstr,"%s/G%03d.mtx",dumpdir,n);
        MatrixWriteTxt(tmpstr, mtp.Glist[0]);
      }
    }
  }

  printf("mris_target_pos done\n");
  return(0);
  exit(0);

}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv;
  char *option ;
  //FILE *fp;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (CMDstringMatch(option, "--help"))  print_help() ;
    else if (CMDstringMatch(option, "--version")) print_version() ;
    else if (CMDstringMatch(option, "--debug"))   debug = 1;
    else if (CMDstringMatch(option, "--checkopts"))   checkoptsonly = 1;
    else if (CMDstringMatch(option, "--cbv"))    DoCBV = 1;
    else if (CMDstringMatch(option, "--no-cbv")) DoCBV = 0;

    else if (CMDstringMatch(option, "--dump")) {
      if (nargc < 1) CMDargNErr(option,1);
      dumpdir = pargv[0];
      nargsused = 1;
    }
    else if (CMDstringMatch(option, "--v")) {
      if (nargc < 1) CMDargNErr(option,1);
      involname = pargv[0];
      nargsused = 1;
    }
    else if (CMDstringMatch(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);
      insurfname = pargv[0];
      nargsused = 1;
    }
    else if (CMDstringMatch(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      outsurfname = pargv[0];
      nargsused = 1;
    }
    else if (CMDstringMatch(option, "--l")) {
      if (nargc < 1) CMDargNErr(option,1);
      riplabelfile = pargv[0];
      nargsused = 1;
    }
    else if (CMDstringMatch(option, "--npointset")) {
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&npointset);
      nargsused = 1;
    }
    else if(CMDstringMatch(option, "--adgws")){
      if(nargc < 1) CMDargNErr(option,1);
      adgwsinfile = pargv[0];
      int err = mtp.adgws.Read(pargv[0]);
      if(err){
	printf("ERROR: reading %s\n",pargv[0]);
	exit(1);
      }
      nargsused = 1;
    } 
    else if (CMDstringMatch(option, "--sigma")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf( pargv[0],"%lf",&mtp.sigma);
      nargsused = 1;
    }
    else if (CMDstringMatch(option, "--interp")) {
      if (nargc < 1) CMDargNErr(option,1);
      interpmethod = pargv[0];
      nargsused = 1;
      if (!strcmp(interpmethod,"sinc") && CMDnthIsArg(nargc, pargv, 1)) {
        sscanf(pargv[1],"%d",&sinchw);
        nargsused ++;
      }
    } 
    else if (CMDstringMatch(option, "--trilin")) {
      interpmethod = "trilinear";
    } 
    else if (CMDstringMatch(option, "--nearest")) {
      interpmethod = "nearest";
    } 
    else if (CMDstringMatch(option, "--cubic")) {
      interpmethod = "cubic";
    } 
    else if (CMDstringMatch(option, "--debug-vertex")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--threads") || !strcasecmp(option, "--nthreads") ){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nthreads);
      #ifdef _OPENMP
      omp_set_num_threads(nthreads);
      #endif
      nargsused = 1;
    } 
    else if (CMDstringMatch(option, "--nii")) format = "nii";
    else if (CMDstringMatch(option, "--nii.gz")) format = "nii.gz";
    else if (CMDstringMatch(option, "--mgh")) format = "mgh";
    else if (CMDstringMatch(option, "--mgz")) format = "mgz";
    else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (CMDsingleDash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}


/* --------------------------------------------- */
static void print_usage(void) 
{
  printf("\n");
  printf("USAGE: mris_target_pos\n");
  printf("\n");
  printf("   --v inputvolume\n");
  printf("   --i inputsurf\n");
  printf("   --o outputsurf\n");
  printf("   --adgws adgwsfile\n");
  printf("   --l label\n");
  printf("   --interp interpname\n");
  printf("   --debug-vertex vtxno\n");
  printf("   --cbv\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  return;
}

/* --------------------------------------------- */
static void check_options(void) {
  if(insurfname == NULL){
    printf("ERROR: must specify an input surface\n");
    exit(1);
  }
  if(outsurfname == NULL){
    printf("ERROR: must specify an output surface\n");
    exit(1);
  }
  if(involname == NULL){
    printf("ERROR: must specify an input volume\n");
    exit(1);
  }
  if(adgwsinfile == NULL){
    printf("ERROR: must specify an input adgws with --adgws\n");
    exit(1);
  }
  mtp.interptype = MRIinterpCode(interpmethod);
  if(mtp.interptype < 0) {
    printf("ERROR: interpolation method %s unrecognized\n",interpmethod);
    printf("       legal values are nearest, trilin, and sinc\n");
    exit(1);
  }
  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fprintf(fp,"inputsurf %s\n",insurfname);
  fprintf(fp,"outputsurf %s\n",outsurfname);
  fprintf(fp,"interp %s\n",interpmethod);

  return;
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n");
  exit(1) ;
}

/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}

int SavePointSet(const char *fname, MRIS *surf, int npoints)
{
  int n,vtxno;
  VERTEX *v;
  FILE *fp;

  FloatInt *flist = (FloatInt*)calloc(sizeof(FloatInt),surf->nvertices);
  int nlist = 0;
  for (n=0; n < surf->nvertices; n++) {
    vtxno = n;
    v = &(surf->vertices[vtxno]);
    if(v->ripflag) continue;
    flist[nlist].f = v->d;
    flist[nlist].i = n;
    nlist++;
  }
  qsort(flist, nlist, sizeof(FloatInt), FloatInt::compare);

  int nsave = MIN(npoints,nlist);

  fp = fopen(fname,"w");
  for (n=0; n < nsave; n++) {
    vtxno = flist[nlist-n-1].i;
    v = &(surf->vertices[vtxno]);
    fprintf(fp,"%g %g %g\n",v->x,v->y,v->z);
    //printf("%d %g %g %g  %g\n",n,v->x,v->y,v->z,v->d);
  }
  fprintf(fp,"info\n");
  fprintf(fp,"numpoints %d\n",npoints);
  fprintf(fp,"useRealRAS 0\n");
  fclose(fp);

  free(flist);
  flist = NULL;

  return(0);
}
