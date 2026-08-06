/**
 * @brief Places surface based on an intensity input image. This is meant to provide
 * replacement functionality for mris_make_surfaces in a form that is easier to 
 * maintain.
 */
/*
 * Original Author: Douglas N Greve (but basically a rewrite of mris_make_surfaces by BF)
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

#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "surfgrad.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "surfcluster.h"
#include "fio.h"
#include "version.h"
//#include "label.h"
//#include "annotation.h"
#include "cmdargs.h"
#include "romp_support.h"
#include <string>
#include <iostream>
#include <fstream>

class MRISsphere2mri {
public:
  MRIS *sphere=NULL; // just sphere, eg lh.sphere
  MRIS *surfreg=NULL; // sphere.reg of this subject
  LTA *lta = NULL;
  //LABEL *label=NULL;
  double lhfsasph1xyz[3] = {11.40, 15.42, -98.14}; //vno 136722 
  double lhfsasph2xyz[3] = {-4.52, 30.08, -95.26}; //vno 52903
  double lhfsasph3xyz[3] = {2.92, 36.70, -92.98}; //vno 136674
  double rhfsasph1xyz[3] = {-32.1000, 23.9500, -91.6300}; //90206
  double rhfsasph2xyz[3] = {-10.3700, 40.0200, -91.0500}; //89549
  double rhfsasph3xyz[3] = {-17.4600, 42.9400, -88.6100}; //13969
  //90206 89549 13969
  int fsavno1=-1, fsavno2=-1;
  int nvox[2] = {-1,-1};
  double P0[2] = {-180,-90}, P1[2] = {180,90};
  double dthresh = 2;
  MRI *sphere2mri(MRI *ov, int direction=1);
  void xyz2rpt(double x, double y, double z, double *r, double *phi, double *theta){
    *r = sqrt(x*x+y*y+z*z);
    *theta = 90 - acos(z/(*r))*180/M_PI; // -90 to +90
    *phi = atan2(y,x)*180/M_PI;     // -180 to +180
  }
  void rpt2xyz(double r, double phi, double theta, double *x, double *y, double *z){
    double th = (90-theta)*M_PI/180;
    *x = r*sin(th)*cos(phi*M_PI/180);
    *y = r*sin(th)*sin(phi*M_PI/180);
    *z = r*cos(th);
  }
  void rpt2xyz(double r, double phi, double theta, float *x, float *y, float *z){
    double xx,yy,zz;
    this->rpt2xyz(r,phi,theta,&xx,&yy,&zz);
    *x = xx; *y = yy; *z = zz;
  }
};
MRI *MRISsphere2mri::sphere2mri(MRI *ov, int direction)
{
  // direction=1=sphere2mri
  // direction=2=mri2sphere
  if(direction != 1 && direction != 2){
    printf("ERROR: diretion=%d must be 1 or 2\n",direction);
    return(NULL);
  }

  // Three points to define the common coordinate system. Map them to the current surface via sphere.reg=surfreg
  float dmin;
  int vno1, vno2, vno3;
  if(surfreg->hemisphere == 0){
    vno1 = MRISfindClosestVertex(surfreg, lhfsasph1xyz[0], lhfsasph1xyz[1], lhfsasph1xyz[2], &dmin, CURRENT_VERTICES);
    vno2 = MRISfindClosestVertex(surfreg, lhfsasph2xyz[0], lhfsasph2xyz[1], lhfsasph2xyz[2], &dmin, CURRENT_VERTICES);
    vno3 = MRISfindClosestVertex(surfreg, lhfsasph3xyz[0], lhfsasph3xyz[1], lhfsasph3xyz[2], &dmin, CURRENT_VERTICES);
  }
  else{
    vno1 = MRISfindClosestVertex(surfreg, rhfsasph1xyz[0], rhfsasph1xyz[1], rhfsasph1xyz[2], &dmin, CURRENT_VERTICES);
    vno2 = MRISfindClosestVertex(surfreg, rhfsasph2xyz[0], rhfsasph2xyz[1], rhfsasph2xyz[2], &dmin, CURRENT_VERTICES);
    vno3 = MRISfindClosestVertex(surfreg, rhfsasph3xyz[0], rhfsasph3xyz[1], rhfsasph3xyz[2], &dmin, CURRENT_VERTICES);
  }

  // From here on out we use the unregistered ?h.sphere (less distortion)

  // vno1 and vno2 are two vertices that dictate the direction of the columns/width
  VERTEX *v1 = &(sphere->vertices[vno1]);
  VERTEX *v2 = &(sphere->vertices[vno2]);
  VERTEX *v3 = &(sphere->vertices[vno3]);
  //double d12 = sqrt( pow(v1->x-v2->x,2.0) + pow(v1->y-v2->y,2.0) + pow(v1->z-v2->z,2.0));
  //double radius = sqrt( pow(v1->x,2.0) + pow(v1->y,2.0) + pow(v1->z,2.0));
  double radius, phi1, theta1, phi2, theta2, phi3, theta3;
  this->xyz2rpt(v1->x,v1->y,v1->z,&radius,&phi1,&theta1);
  this->xyz2rpt(v2->x,v2->y,v2->z,&radius,&phi2,&theta2);
  this->xyz2rpt(v3->x,v3->y,v3->z,&radius,&phi3,&theta3);
  printf("Prior to any +=360 phi %6.4lf %6.4lf %6.4lf\n",phi1,phi2,phi3);
  // This can happen on rh
  if(phi1 < -90) phi1 += 360;
  if(phi2 < -90) phi2 += 360;
  if(phi3 < -90) phi3 += 360;
  double cphi = (phi1+phi2)/2;
  double ctheta = (theta1+theta2)/2;
  double calpha = atan2((theta2-theta1),(phi2-phi1))*180/M_PI;
  if(surfreg->hemisphere != 0) calpha = 180-calpha;
  double cx,cy,cz;
  this->rpt2xyz(radius,cphi,ctheta,&cx,&cy,&cz);
  int cvno = MRISfindClosestVertex(sphere, cx, cy, cz, &dmin, CURRENT_VERTICES);
  VERTEX *cv = &(sphere->vertices[cvno]);
  printf("vno1=%d (%g,%g,%g) vno2=%d (%g,%g,%g) vno3=%d (%g,%g,%g) \n",vno1,
	 v1->x,v1->y,v1->z, vno2,v2->x,v2->y,v2->z, vno3,v3->x,v3->y,v3->z);
  printf("phi = (%g,%g,%g) theta=(%g,%g,%g)\n",phi1,phi2,phi3,theta1,theta2,theta3);
  printf("center phi=%g theta=%g alpha=%g %4.1f %4.1f %4.1f  cvno=%d  %4.1f %4.1f %4.1f  d=%g\n",
	 cphi,ctheta,calpha,cx,cy,cz,cvno,cv->x,cv->y,cv->z,dmin);
  printf("radius %g\n",radius);
  printf("v1 %d %g %g %g\n",vno1,v1->x,v1->y,v1->z);
  printf("v2 %d %g %g %g\n",vno2,v2->x,v2->y,v2->z);

  double par[12];
  for(int k=0; k<12; k++) par[k]=0;
  par[3] = calpha; // rotation about x
  par[4] = ctheta; // rotation about y
  par[5] = cphi; // rotation about z (applied first)
  par[6] = 1; par[7] = 1; par[8] = 1;
  MATRIX *T = TranformAffineParams2Matrix(par, NULL);

  //VERTEX *v = &(sphere->vertices[cvno]);
  //printf("before %g %g %g\n",v->x,v->y,v->z);
  MRISmatrixMultiply(sphere, T);
  //MRISwrite(sphere,"lh.sphere.T");
  //printf("after  %g %g %g\n",v->x,v->y,v->z);
  MatrixFree(&T);
  MHT *Hash = MHTcreateVertexTable_Resolution(sphere, CURRENT_VERTICES, 16);

  MRI *ovsurf=NULL, *ovslice=NULL;
  if(direction == 1){
    ovsurf = ov;
    ovslice = MRIallocSequence(nvox[0],nvox[1],1,ovsurf->type,ovsurf->nframes);
    MRIcopyPulseParameters(ovsurf,ovslice);
    if(ovsurf->ct) ovslice->ct = CTABdeepCopy(ovsurf->ct);  
    ovslice->valid = 1;
    
    // The coordinates of the MRI structure will be in phi/theta degrees
    ovslice->xsize = fabs((P1[0]-P0[0])/nvox[0]);
    ovslice->ysize = fabs((P1[1]-P0[1])/nvox[1]);
    ovslice->zsize = 1; // flat in this dim

    // The center
    ovslice->c_r  = 0.5;
    ovslice->c_a  = (P1[0]+P0[0])/2;
    ovslice->c_s  = (P1[1]+P0[1])/2;
    // Set col axis to A->P 
    ovslice->x_r = 0;
    ovslice->x_a = 1;
    ovslice->x_s = 0;
    // Set row axis to be orthog to col axis
    ovslice->y_r = 0;
    ovslice->y_a = 0;
    ovslice->y_s = 1;
    // Set slice axis to point in x
    ovslice->z_r = 1;
    ovslice->z_a = 0;
    ovslice->z_s = 0;
  } else {
    ovslice = ov;
    ovsurf = MRIallocSequence(sphere->nvertices,1,1,ovslice->type,ovslice->nframes);
    MRIcopyPulseParameters(ovslice,ovsurf);
    if(ovslice->ct) ovsurf->ct = CTABdeepCopy(ovslice->ct);  
  }
  MATRIX *vox2ras = ovslice->get_Vox2RAS(0);
  printf("x_r %g x_a %g\n",ovslice->x_r,ovslice->x_a);
  printf("center %g %g %g\n",ovslice->c_r,ovslice->c_a,ovslice->c_s);
  //printf("vox2ras = [\n");
  //MatrixPrint(stdout,vox2ras);
  //printf("];\n");

  // Go through each point in the mri and find closest vertex in patch
  // Probably a more efficient way to do this with hash, but the data
  // are pretty small. Could probably improve the interp
  MATRIX *vox = MatrixAlloc(4,1,MATRIX_REAL);
  vox->rptr[4][1] = 1;
  MATRIX *phitheta=NULL;
  for(int c=0; c < ovslice->width; c++){
    for(int r=0; r < ovslice->height; r++){
      if(surfreg->hemisphere == 0) vox->rptr[1][1] = c; //lh 
      if(surfreg->hemisphere == 1) vox->rptr[1][1] = ovslice->width - c - 1; //rh
      vox->rptr[2][1] = r;
      vox->rptr[3][1] = 0;
      phitheta = MatrixMultiplyD(vox2ras,vox,phitheta);
      double p = phitheta->rptr[2][1];
      double t = phitheta->rptr[3][1];
      double x,y,z;
      this->rpt2xyz(radius,p,t,&x,&y,&z);
      int vno = Hash->findClosestVertexNoXYZ(x,y,z, &dmin);
      for(int f=0; f < ovslice->nframes; f++){
	double val;
	if(direction == 1){
	  val = MRIgetVoxVal(ovsurf,vno,0,0,f);
	  MRIsetVoxVal(ovslice,c,r,0,f,val);
	} else {
	  val = MRIgetVoxVal(ovslice,c,r,0,f);
	  // If the slice is 0 here and if the surf already has a
	  // non-zero value here, then don't overwrite it. This is to
	  // try to prevent loss of non-zero voxels at the edge. These
	  // may be lost in the forward (sphere2mri) loop.
	  if(val == 0 && MRIgetVoxVal(ovsurf,vno,0,0,f)!=0) continue;
	  MRIsetVoxVal(ovsurf,vno,0,0,f,val);
	}
      }
    }
  }

  if(direction == 2){ // should only be done with seg
    // Fill holes: Bin. Invert. Cluster. Loop through voxels in clusters 2-N
    MATRIX *ras2vox = MatrixInverse(vox2ras,NULL);
    SURFCLUSTERSUM *scs;
    int NClusters=0;
    MRIScopyMRI(sphere, ovsurf, 0, "val");
    // Set marked by binarizing map as it is
    for(int vno=0; vno < sphere->nvertices; vno++) {
      VERTEX *v = &(sphere->vertices[vno]);
      if(v->val > .001) v->marked = 1;
      else              v->marked = 0;
    }
    // Close using marked by dilating by 2 and then eroding by 1. The
    // goal here is only to make sure that any voids are actual holes
    // and not inlets.
    MRISdilateMarked(sphere,2);
    MRISerodeMarked(sphere,1);
    // Reset val to be the inverse of the closed marked
    for(int vno=0; vno < sphere->nvertices; vno++) {
      VERTEX *v = &(sphere->vertices[vno]);
      v->val = !v->marked;
    }
    // Cluster to find any holes
    scs = sclustMapSurfClusters(sphere,.0001,-1,0,0,&NClusters,NULL,NULL);
    printf("nc = %d\n",NClusters);
    for(int vno=0; vno < sphere->nvertices; vno++){
      VERTEX *v = &(sphere->vertices[vno]);
      int cno = v->undefval;
      // If it is not marked and not in a hole then skip vertex
      if(! v->marked && cno < 2) continue;
      double phi, theta;
      phitheta->rptr[4][1] = 1; // just to make sure
      this->xyz2rpt(v->x,v->y,v->z,&radius,&phi,&theta);
      phitheta->rptr[1][1] = 0;
      phitheta->rptr[2][1] = phi; 
      phitheta->rptr[3][1] = theta;
      vox = MatrixMultiplyD(ras2vox,phitheta,vox);
      int c = round(vox->rptr[1][1]);
      int r = round(vox->rptr[2][1]);
      if(c < 0 || c >= ovslice->width) continue;
      if(r < 0 || r >= ovslice->height) continue;
      if(surfreg->hemisphere == 1) c = ovslice->width - c  - 1; //rh
      double val = MRIgetVoxVal(ovslice,c,r,0,0);
      //printf("cno=%d vno = %d  pt=(%g,%g) cr = (%d %d)  %g\n",cno,vno,phi,theta,c,r,val);
      MRIsetVoxVal(ovsurf,vno,0,0,0,val);
    }
    MatrixFree(&ras2vox);
  }

  MatrixFree(&vox2ras);
  MatrixFree(&vox);
  MatrixFree(&phitheta);
  if(direction == 1) return(ovslice);
  else return(ovsurf);
}

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);

struct utsname uts;
char *cmdline, cwd[2000];
int debug = 0, checkoptsonly = 0;

int main(int argc, char *argv[]) ;

const char *Progname = "mris_sphere2mris";
char *spherefile=NULL, *sphereregfile=NULL, *ltafile=NULL, *ctabfile=NULL;
char *ovsurffile=NULL, *ovslicefile=NULL;
int direction = 1;

/*--------------------------------------------------*/
int main(int argc, char **argv) 
{
  int nargs, i, msec;
  double        spring_scale = 1;
  MRI *invol, *seg=NULL, *wm, *involCBV, *involPS;
  Timer timer ;
  char *cmdline2, cwd[2000];
  //char *field=NULL;

  nargs = handleVersionOption(argc, argv, "mris_sphere2mri");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);
  cmdline2 = argv2cmdline(argc,argv);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  Gdiag |= DIAG_SHOW ;

  if(argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if(checkoptsonly) return(0);

  MRI *mri=NULL;
  MRISsphere2mri s2m;
  s2m.sphere = MRISread(spherefile);
  s2m.surfreg = MRISread(sphereregfile);
  s2m.lta = LTAread(ltafile);
  if(s2m.lta->type != REGISTER_DAT){
    printf("Changing LTA to TKR\n");
    LTAchangeType(s2m.lta, REGISTER_DAT);
  }
  LTAmat2RotMat(s2m.lta);
  MRISltaMultiply(s2m.sphere,s2m.lta);
  if(direction == 1){
    s2m.nvox[0] = 400;
    s2m.nvox[1] = 200;
    s2m.P0[0] = -90; s2m.P0[1] = -45;
    s2m.P1[0] = +90; s2m.P1[1] = +45;
    mri = MRIread(ovsurffile); //ovsurf is input, ovslice output
  } else {
    mri = MRIread(ovslicefile); //ovslice is input, ovsurf output
  }
  if(ctabfile) mri->ct = CTABreadASCII(ctabfile);
  MRI *ss = s2m.sphere2mri(mri,direction);
  if(direction == 1)  MRIwrite(ss,ovslicefile);
  else                MRIwrite(ss,ovsurffile);

  printf("#VMPC# mris_sphere2mris VmPeak  %d\n",GetVmPeak());
  printf("mris_sphere2mri done\n");

  return(0);

}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if(!strcasecmp(option, "--version")) print_version() ;
    else if(!strcasecmp(option, "--debug"))   debug = 1;
    else if(!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if(!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if(!strcmp(option, "--surf-to-slice")) direction = 1;
    else if(!strcmp(option, "--slice-to-surf")) direction = 2;
    else if(!strcmp(option, "--sphere")){
      if(nargc < 1) CMDargNErr(option,1);
      spherefile = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--sphere.reg")){
      if(nargc < 1) CMDargNErr(option,1);
      sphereregfile = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--lta")){
      if(nargc < 1) CMDargNErr(option,1);
      ltafile = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--ov-surf")){
      if(nargc < 1) CMDargNErr(option,1);
      ovsurffile = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--ov-slice")){
      if(nargc < 1) CMDargNErr(option,1);
      ovslicefile = pargv[0];
      nargsused = 1;
    }
    else if(!strcmp(option, "--ctab")){
      if(nargc < 1) CMDargNErr(option,1);
      ctabfile = pargv[0];
      nargsused = 1;
    }

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
static void check_options(void) {
  if(spherefile == NULL){
    printf("ERROR: must specify --sphere\n");
    exit(1);
  }
  if(sphereregfile == NULL){
    printf("ERROR: must specify --sphere.reg\n");
    exit(1);
  }
  if(ovsurffile == NULL){
    printf("ERROR: must specify --ov-surf\n");
    exit(1);
  }
  if(ovslicefile == NULL){
    printf("ERROR: must specify --ov-slice\n");
    exit(1);
  }
  if(ltafile == NULL){
    printf("ERROR: must specify --lta\n");
    exit(1);
  }
  return;
}



#include "mris_sphere2mri.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mris_sphere2mri_help_xml,mris_sphere2mri_help_xml_len);
}

static void
print_help(void)
{
  print_usage() ;
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
  return;
}

