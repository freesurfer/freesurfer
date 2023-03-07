/**
 * @name mris_skeletonize
 * @brief skeletonizes the curvature to create labels of the gyral
	crowns or the sulcal fundi.
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
#include "annotation.h"
#include "error.h"
#include "dmatrix.h"
#include "matrix.h"
#include "diag.h"
#include "pointset.h"
#include "dtk.fs.h"
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
static void dump_options();

class SurfSkeleton 
{
public:
  MRIS *m_surf=NULL,*m_sphere=NULL;
  MRI *m_mask=NULL;
  std::vector<std::vector<int>> vtxnbrlist; // nbrs in spatially contiguous order
  int NeighborhoodSize = 2;
  double m_threshold = 0.3;
  double m_fwhm = 0;
  std::vector<std::pair<double,int>> edgelist; // val,vno
  std::vector<double> surfvals;
  int build_vtxnbrlist(void){
    // For each vertex, get the list of neighboring vertices in the 
    // proper spatial order
    vtxnbrlist.clear();
    for(int vno=0; vno < m_surf->nvertices; vno++){
      std::vector<int> nbrlist = MRISgetSpatOrderedNeighbors(m_surf,vno);
      vtxnbrlist.push_back(nbrlist);
    }
    return(0);
  };
  int isedge(int vno){
    // Test whether a given vertex is on the edge of the mask
    if(!MRIgetVoxVal(m_mask,vno,0,0,0)) return(0); // do this here?
    VERTEX_TOPOLOGY const * const vt = &m_surf->vertices_topology[vno];
    // count the number of neighbors that are 0 and the number that are 1
    // if both are greater than 0, then it has to be an edge
    int n0=0, n1=0;
    for(int n = 0; n < vt->vnum; n++) {
      int vnvno = vt->v[n];
      int m = MRIgetVoxVal(m_mask,vnvno,0,0,0);
      if(m) n1++;
      else  n0++;
      if(n0>0 && n1>0) return(1);
    }
    return(0);
  };
  int erodeable(int vno) {
    // Test whether a given vertex can be eroded without causing
    // the skeleton to become disconnected
    if(!MRIgetVoxVal(m_mask,vno,0,0,0)) return(0);
    // count the number of 0 islands (clusters)
    std::vector<int> nbrlist = vtxnbrlist[vno];
    int nzc=0,mprev=0,nnz=0,first0=0,last0=0;
    for(int n = 0; n < nbrlist.size(); n++) {
      int vnvno = nbrlist[n];
      int m = MRIgetVoxVal(m_mask,vnvno,0,0,0);
      nnz += m; // keep track of the number of non-zero nbrs
      if(n==0 && m==0) first0=1;
      if(n==nbrlist.size()-1 && m==0) last0=1;
      if(m!=0) {
	mprev = m;
	continue;
      }
      if(n==0 || mprev) nzc++; // found a potentially new cluster (m=0 here)
      mprev = m;
    }
    if(first0 && last0 && nzc > 1) nzc--; // cluster was not new, wrap-around
    if(nzc > 1) return(0); // more than 1 cluster/island, so can't erode
    if(nnz == 1) return(0); // this is a terminal vertex, can't erode
    return(1);
  };
  int set_edgelist(void){
    // Get a list of the current edge vertices
    edgelist.clear();
    for(int vno=0; vno < m_surf->nvertices; vno++){
      if(!isedge(vno)) continue;
      edgelist.push_back(std::make_pair(surfvals[vno],vno));
    }
    return((int)edgelist.size());
  }
  void load_k1(double scale){
    // Compute k1 and load it into surfvals; apply the given scale (-1 or gyrus, +1 for sulcus)
    printf("k1 scale %g\n",scale);
    surfvals.clear();
    MRISresetNeighborhoodSize(m_surf, NeighborhoodSize) ; 
    MRIScomputeSecondFundamentalForm(m_surf) ;
    for(int vno=0; vno < m_surf->nvertices; vno++) surfvals.push_back(scale*m_surf->vertices[vno].k1);
  }
  void load_k1_gyrus(void) {load_k1(-1);}
  void load_k1_sulcus(void){load_k1(+1);}
  void load_curv_nonmaxsup(double scale){
    // Compute non-max suppressed curv and load it into surfvals; apply the given scale (-1 or gyrus, +1 for sulcus)
    printf("nonmaxsup scale = %g\n",scale);
    surfvals.clear();
    MRISresetNeighborhoodSize(m_surf, NeighborhoodSize) ; 
    MRIScomputeSecondFundamentalForm(m_surf) ; 
    //MRISresetNeighborhoodSize(m_sphere, NeighborhoodSize) ; 
    for(int vno=0; vno < m_surf->nvertices; vno++) m_sphere->vertices[vno].curv = m_surf->vertices[vno].H;
    //for(int vno=0; vno < 100; vno++) printf("%d %g %g\n",vno,m_surf->vertices[vno].H,m_sphere->vertices[vno].curv);
    MRISnonmaxSuppress(m_sphere) ;
    for(int vno=0; vno < m_surf->nvertices; vno++) surfvals.push_back(scale*m_sphere->vertices[vno].curv);
  }
  void load_curv_nonmaxsup_gyrus(void)  {load_curv_nonmaxsup(-1);}
  void load_curv_nonmaxsup_sulcus(void) {load_curv_nonmaxsup(+1);}
  MRI *surfvals2mri(void){
    // Load surfvals into an MRI structure
    MRI *mri = MRIalloc(m_surf->nvertices,1,1,MRI_FLOAT);
    for(int vno=0; vno < m_surf->nvertices; vno++) MRIsetVoxVal(mri,vno,0,0,0,surfvals[vno]);
    return(mri);
  }
  void mri2surfvals(MRI *mri){
    // Copy MRI structure into surfvals 
    surfvals.clear();
    for(int vno=0; vno < m_surf->nvertices; vno++) surfvals.push_back(MRIgetVoxVal(mri,vno,0,0,0));
  }
  int smooth_surfvals(double fwhm){
    // Smooth surfvals by the given FWHM
    int niters = MRISfwhm2niters(fwhm,m_surf);
    printf("Smoothing input by fwhm=%lf, niters=%d \n",fwhm,niters);
    MRI *mri = surfvals2mri();
    MRISsmoothMRI(m_surf, mri, niters, NULL, mri);
    mri2surfvals(mri);
    MRIfree(&mri);
    return(0);
  }
  int write_surfvals(char *fname){
    // Write surfvals as an MRI file
    MRI *mri = surfvals2mri();
    int err = MRIwrite(mri,fname);
    MRIfree(&mri);
    return(err);
  }
  int write_pointset(char *fname){
    // Write current edge list vertices as a point set. When the
    // procedure is done, this will be the skeleton.
    fsPointSet ps;
    MATRIX *tkras2scannerras = m_surf->vg.get_TkregRAS2RAS();
    MATRIX *tkxyz = MatrixAlloc(4,1,MATRIX_REAL);
    tkxyz->rptr[4][1] = 1;
    MATRIX *scannerxyz = NULL;
    for(int n=0; n < edgelist.size(); n++){
      int vno = edgelist[n].second;
      fsPointSet::Point p;
      tkxyz->rptr[1][1] = m_surf->vertices[vno].x;
      tkxyz->rptr[2][1] = m_surf->vertices[vno].y;
      tkxyz->rptr[3][1] = m_surf->vertices[vno].z;
      scannerxyz = MatrixMultiplyD(tkras2scannerras,tkxyz,scannerxyz);
      p.x = scannerxyz->rptr[1][1];
      p.y = scannerxyz->rptr[2][1];
      p.z = scannerxyz->rptr[3][1];
      ps.add(p);
    }
    ps.vox2ras = "scanner_ras";
    int err = ps.save(fname);
    MatrixFree(&tkras2scannerras);
    MatrixFree(&tkxyz);
    MatrixFree(&scannerxyz);
    return(err);
  }

  int skeletonize(void){
    // When this is done, the m_mask will be a binary mask of the skeleton and edgelist
    // will be a list of the vertices in the skeleton.
    if(surfvals.size()==0){
      printf("ERROR: skeletonize(): no input\n");
      return(1);
    }
    if(m_fwhm > 0) smooth_surfvals(m_fwhm);
    if(m_mask) MRIfree(&m_mask);
    m_mask = MRIalloc(m_surf->nvertices,1,1,MRI_INT);
    for(int vno=0; vno < m_surf->nvertices; vno++) if(surfvals[vno] > m_threshold) MRIsetVoxVal(m_mask,vno,0,0,0, 1);
    build_vtxnbrlist();

    int iter=0, nerodedtot=0;
    while(1){
      int neroded=0;
      iter++;
      set_edgelist();
      sort(edgelist.begin(),edgelist.end());// sorts lowest to highest
      for(int eno=0; eno < edgelist.size(); eno++){
	int vno = edgelist[eno].second;
	if(!erodeable(vno)) continue;
	MRIsetVoxVal(m_mask,vno,0,0,0, 0);
	neroded++;
      }
      printf("i=%d nedge=%d neroded=%d\n",iter,(int)edgelist.size(),neroded);
      nerodedtot += neroded;
      if(neroded == 0) break;
    }
    printf("nedges = %d, nerodedtot = %d\n",(int)edgelist.size(),nerodedtot);
    return(0);
  }
};


/*---------------------------------------------------------------------*/
struct utsname uts;
char *cmdline, cwd[2000];
int checkoptsonly = 0;
const char *Progname = NULL;
char *surfpath=NULL, *spherepath=NULL, *outmaskpath=NULL,*outsurfvalspath=NULL,*pointsetpath=NULL,*labelpath=NULL,*tractpath=NULL;
SurfSkeleton sk;
double scale = -1; // gyrus, crowns
char *subject=NULL, *hemi=NULL, *surfname=NULL;

int main(int argc, char** argv)
{
  int nargs,err;

  nargs = handleVersionOption(argc, argv, "mris_skeletonize");
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

  sk.m_surf = MRISread(surfpath);
  if(!sk.m_surf) exit(1);
  sk.m_sphere = MRISread(spherepath);
  if(!sk.m_sphere) exit(1);

  sk.load_curv_nonmaxsup(scale);

  sk.skeletonize();
  if(outmaskpath) {
    err = MRIwrite(sk.m_mask,outmaskpath);
    if(err) exit(err);
  }
  if(outsurfvalspath){
    err = sk.write_surfvals(outsurfvalspath);
    if(err) exit(err);
  }
  if(pointsetpath){
    err = sk.write_pointset(pointsetpath);
    if(err) exit(err);
  }
  printf("mris_skeletonize done\n");

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
    else if (!strcasecmp(option, "--gyrus"))  scale = -1; 
    else if (!strcasecmp(option, "--sulcus")) scale = +1; 

    else if (!strcasecmp(option, "--s")){
      if(nargc < 3) CMDargNErr(option,3);
      subject = pargv[0];
      hemi = pargv[1];
      surfname = pargv[2];
      nargsused = 3;
    }
    else if (!strcasecmp(option, "--surf")){
      if(nargc < 1) CMDargNErr(option,1);
      surfpath = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--sphere")){
      if(nargc < 1) CMDargNErr(option,1);
      spherepath = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--mask")){
      if(nargc < 1) CMDargNErr(option,1);
      outmaskpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--nbrsize")){
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&sk.NeighborhoodSize);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--threshold") || !strcasecmp(option, "--thresh")){
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&sk.m_threshold);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--fwhm")){
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&sk.m_fwhm);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--surfvals")){
      if (nargc < 1) CMDargNErr(option,1);
      outsurfvalspath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--ps")){
      if (nargc < 1) CMDargNErr(option,1);
      pointsetpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--l")){
      if (nargc < 1) CMDargNErr(option,1);
      labelpath = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--tract")){
      if (nargc < 1) CMDargNErr(option,1);
      tractpath = pargv[0];
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
  dump_options();

  // atlasMesh is required
  if(surfpath == NULL){
    printf("ERROR: must spec --surf\n");
    exit(1);
  }
  if(spherepath == NULL){
    printf("ERROR: must spec --sphere\n");
    exit(1);
  }
  if(outmaskpath==NULL && outsurfvalspath==NULL && pointsetpath==NULL && labelpath==NULL && tractpath==NULL){
    printf("ERROR: not outputf specified\n");
    exit(1);
  }

  return;
}
/* ------------------------------------------------------ */
#include "mris_skeletonize.help.xml.h"
static void print_usage()
{
  outputHelpXml(mris_skeletonize_help_xml, mris_skeletonize_help_xml_len);
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
static void dump_options() {
  printf("\n");
  printf("%s\n", getVersion().c_str());
  printf("cwd %s\n", cwd);
  printf("cmdline %s\n", cmdline);
  printf("sysname  %s\n", uts.sysname);
  printf("hostname %s\n", uts.nodename);
  printf("machine  %s\n", uts.machine);
  printf("\n");
}
