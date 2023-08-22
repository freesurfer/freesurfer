/**
 * @name mri_si_prep
 * @brief prepares a volume for smart interpol
 *
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
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <float.h>
#include <vector>
#include <array>
#include <sys/utsname.h>
#include <sys/stat.h>
#include "mrisurf.h"
#include "mrisurf_topology.h"
#include "mrisutils.h"
#include "timer.h"
#include "utils.h"
#include "fio.h"
#include "mri2.h"
#include "region.h"
#include "error.h"
#include "dmatrix.h"
#include "matrix.h"
#include "diag.h"
#include "colortab.h"
#include "transform.h"
#include "cmdargs.h"
#include "version.h"
#ifdef _OPENMP
#include "romp_support.h"
#endif

static int  parse_commandline(int argc, char **argv);
static void check_options();
static void print_usage();
static void usage_exit();
static void print_help();
static void print_version();
static void dump_options(FILE *fp);
int SIprep(MRI *mri, MRI *seg, std::vector<int> segnos, int nskip, int npad, int interpdim, MRI **pmriout, MRI **psegout);

/*---------------------------------------------------------------------*/
struct utsname uts;
char *cmdline, cwd[2000];
int checkoptsonly = 0;
const char *Progname = NULL;
char *involpath=NULL,*insegpath=NULL, *outvolpath=NULL, *outsegpath=NULL;
std::vector<int> segnos;
int nskip=0, npad = -1, interpdim = -1, nthreads = 1;

int main(int argc, char** argv)
{
  int nargs,err;

  nargs = handleVersionOption(argc, argv, "mri_si_prep");
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
  if(checkoptsonly) return(0);

  printf("Reading %s\n",involpath);
  MRI *invol = MRIread(involpath);
  if(invol == NULL) exit(1);
  printf("Reading %s\n",insegpath);
  MRI *inseg = MRIread(insegpath);
  if(inseg == NULL) exit(1);

  printf("Running prep\n"); fflush(stdout);
  MRI *outvol=NULL, *outseg=NULL;
  err =  SIprep(invol, inseg, segnos, nskip, npad, interpdim, &outvol, &outseg);
  if(err) exit(1);

  printf("Writing %s\n",outvolpath);  
  err = MRIwrite(outvol,outvolpath);
  if(err) exit(1);

  printf("Writing %s\n",outsegpath);  
  err = MRIwrite(outseg,outsegpath);
  if(err) exit(1);

  printf("#VMPC# mri_si_prep VmPeak  %d\n",GetVmPeak());
  printf("mri_si_prep done\n");

  return 0;
}

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))             print_help() ;
    else if (!strcasecmp(option, "--version"))     print_version() ;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if(!strcasecmp(option, "--i")){
      if(nargc < 2) CMDargNErr(option,2);
      involpath = pargv[0];
      insegpath = pargv[1];
      nargsused = 2;
    }
    else if(!strcasecmp(option, "--o")){
      if(nargc < 2) CMDargNErr(option,2);
      outvolpath = pargv[0];
      outsegpath = pargv[1];
      nargsused = 2;
    }
    else if (!strcasecmp(option, "--segno")){
      if(nargc < 1) CMDargNErr(option,1);
      int segno;
      sscanf(pargv[0],"%d",&segno);
      segnos.push_back(segno);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--nskip")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nskip);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--npad")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&npad);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--dim")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&interpdim);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--ax"))  interpdim=4;
    else if (!strcasecmp(option, "--cor")) interpdim=5;
    else if (!strcasecmp(option, "--sag")) interpdim=6;
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
      printf("ERROR: Option %s unknown\n", option);
      if (CMDsingleDash(option))
        printf("       Did you really mean -%s ?\n\n", option);
      print_help();
      exit(1);
    }

    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* --------------------------------------------- */
static void check_options()
{
  dump_options(stdout);

  if(involpath == NULL){
    printf("ERROR: must spec --i\n");
    exit(1);
  }
  if(outvolpath == NULL){
    printf("ERROR: must spec --o\n");
    exit(1);
  }
  if(segnos.size() == 0){
    printf("ERROR: must spec --segno\n");
    exit(1);
  }
  if(interpdim < 0){
    printf("ERROR: must specify a dimension with --dim or --ax, --cor, --sag\n");
    exit(1);
  }
  if(nskip < 1){
    printf("ERROR: must specify --nskip\n");
    exit(1);
  }

  return;
}
/* ------------------------------------------------------ */
#include "mri_si_prep.help.xml.h"
static void print_usage()
{
  outputHelpXml(mri_si_prep_help_xml, mri_si_prep_help_xml_len);
}
/* ------------------------------------------------------ */
static void print_help()
{
  print_usage();
  exit(1) ;
}
/* ------------------------------------------------------ */
static void usage_exit() {
  print_usage();
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version() {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n", cwd);
  fprintf(fp,"cmdline %s\n", cmdline);
  fprintf(fp,"sysname  %s\n", uts.sysname);
  fprintf(fp,"hostname %s\n", uts.nodename);
  fprintf(fp,"machine  %s\n", uts.machine);
  fprintf(fp,"\n");
}

int SIprep(MRI *mri, MRI *seg, std::vector<int> segnos, int nskip, int npad, int interpdim, MRI **pmriout, MRI **psegout)
{
  int c,cmin=10000, cmax=-1, rmin=10000, rmax=-1, smin=10000, smax=-1;

  printf("SIprep(): nskip=%d npad=%d interpdim=%d segnos=",nskip,npad,interpdim);
  std::vector<int>::iterator itr = segnos.begin();
  while(itr != segnos.end()) {printf("%d ",*itr); itr++;}
  printf("\n");
  printf("mridim %d %d %d\n",mri->width,mri->height,mri->depth);
  printf("segdim %d %d %d\n",seg->width,seg->height,seg->depth);

  // check for dim mismatches
  int err = MRIdimMismatch(mri, seg, 0);
  if(err){
    printf("ERROR: SIprep(): dimension mismatch\n");
    return(1);
  }

  // Determin which dim (1,2,3) to skip
  if(interpdim < 1 || interpdim > 6){
    printf("ERROR: SIprep(): interpdim=%d out of range",interpdim);    
    return(1);
  }
  char ostr[5];
  ostr[4] = '\0';
  MRIdircosToOrientationString(mri,ostr);
  printf("ostring %s\n",ostr);
  int skipdim = -1;
  if(interpdim <= 3) skipdim = interpdim;
  else {
    int n;
    switch(interpdim){
    case 4: // axial
      for(n=0; n < 3; n++)	if(ostr[n] == 'S' || ostr[n] == 'I') break;
      break;
    case 5: // cor
      for(n=0; n < 3; n++)	if(ostr[n] == 'A' || ostr[n] == 'P') break;
      break;
    case 6: // sagittal
      for(n=0; n < 3; n++)	if(ostr[n] == 'L' || ostr[n] == 'R') break;
      break;
    }	  
    skipdim = n+1;
  }
  printf("skipdim = %d\n",skipdim);

  // Get the bounding box. Could have used REGIONgetBoundingBoxM() but not parallel
  int nhits = 0;
#ifdef HAVE_OPENMP
#pragma omp parallel for reduction(max: cmax,rmax,smax) reduction(min: cmin,rmin,smin) reduction(+: nhits) 
#endif
  for(c=0; c < mri->width; c++){
    for(int r=0; r < mri->height; r++){
      for(int s=0; s < mri->depth; s++){
        int segno = MRIgetVoxVal(seg, c, r, s, 0);
	int hit=0;
	std::vector<int>::iterator itr = segnos.begin();
	while(itr != segnos.end()){
	  if(*itr == segno) {
	    hit = 1;
	    break;
	  }
	  itr++;
	}
	if(!hit) continue;
	nhits ++;
        if(cmin > c) cmin = c;
        if(rmin > r) rmin = r;
        if(smin > s) smin = s;
        if(cmax < c) cmax = c;
        if(rmax < r) rmax = r;
        if(smax < s) smax = s;
      }
    }
  }
  printf("Init BB  %3d %3d %3d  %3d %3d %3d\n",cmin,rmin,smin,cmax,rmax,smax);
  printf("nhits=%d\n",nhits);
  if(nhits == 0){
    printf("ERROR: cannot find any matching voxels\n");
    return(1);
  }

  // Apply the padding to the non-skip dims. If npad < 0, then use
  // the full FoV in the non-skip dims.
  if(skipdim != 1){
    if(npad >= 0){
      cmin = MAX(cmin-npad,0);
      cmax = MIN(cmax+npad,mri->width-1);
    } else {
      cmin = 0;
      cmax = mri->width-1;
    }
  }
  if(skipdim != 2){
    if(npad >= 0){
      rmin = MAX(rmin-npad,0);
      rmax = MIN(rmax+npad,mri->height-1);
    } else {
      rmin = 0;
      rmax = mri->height-1;
    }
  }
  if(skipdim != 3){
    if(npad >= 0){
      smin = MAX(smin-npad,0);
      smax = MIN(smax+npad,mri->depth-1);
    } else {
      smin = 0;
      smax = mri->depth-1;
    }
  }
  printf("Final BB %3d %3d %3d  %3d %3d %3d\n",cmin,rmin,smin,cmax,rmax,smax);

  // Crop the MRI and the seg
  MRI_REGION region;
  region.x = cmin;
  region.y = rmin;
  region.z = smin;
  region.dx = cmax-cmin+1;
  region.dy = rmax-rmin+1;
  region.dz = smax-smin+1;
  printf("region: "); REGIONprint(stdout, &region);
  MRI *mriout = MRIextractRegion(mri, NULL, &region);
  if(mriout == NULL) return(1);
  MRIcopyPulseParameters(mri, mriout); // don't copy header
  MRI *segout = MRIextractRegion(seg, NULL, &region);
  if(segout == NULL) return(1);
  MRIcopyPulseParameters(seg, segout); // don't copy header
  if(seg->ct) segout->ct = CTABdeepCopy(seg->ct);

  // Go through the outseg and zero skipped slices and non-segno
  // voxels Should check whether each nskip slice has some label in
  // it. Hard to do in paralle. Does the last slice have to have
  // label? The first slice is guaranteed to have label. 
  nhits = 0;
#ifdef HAVE_OPENMP
#pragma omp parallel for reduction(+: nhits) 
#endif
  for(c=0; c < mriout->width; c++){
    for(int r=0; r < mriout->height; r++){
      for(int s=0; s < mriout->depth; s++){
	int zslice=0;
	if(skipdim == 1) zslice = (c%nskip);
	if(skipdim == 2) zslice = (r%nskip);
	if(skipdim == 3) zslice = (s%nskip);
        int segno = MRIgetVoxVal(segout, c, r, s, 0);
	int hit=0;
	std::vector<int>::iterator itr = segnos.begin();
	while(itr != segnos.end()){
	  if(*itr == segno) {
	    hit = 1;
	    break;
	  }
	  itr++;
	}
	if(!hit || zslice){
	  MRIsetVoxVal(segout,c,r,s,0, 0);
	  continue;
	}
	MRIsetVoxVal(segout,c,r,s,0, segnos[0]); // have to choose one segno
	nhits++;
      }
    }
  }
  printf("nhits %d\n",nhits);

  if(*pmriout) MRIfree(pmriout);
  if(*psegout) MRIfree(psegout);
  *pmriout = mriout;
  *psegout = segout;

  return(0);
}
