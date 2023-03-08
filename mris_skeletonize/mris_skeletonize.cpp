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
#include "surfcluster.h"
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
  std::vector<int> terminals, branchpoints;
  class branchstruct {public: std::vector<int> vtxlist; double length=0;};
  std::vector<branchstruct> branches;  
  int NeighborhoodSize = 2;
  double m_threshold = 0.3;
  int nkeep=0; // number of clusters to keep, 0 means don't do clustering
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
  int nmasknbrs(int vno){
    // Count the number of neighbors in the mask
    std::vector<int> nbrlist = vtxnbrlist[vno];
    int nnbrs = 0;
    for(int n=0; n < nbrlist.size(); n++) if(MRIgetVoxVal(m_mask,nbrlist[n],0,0,0)>0) nnbrs++;
    return(nnbrs);
  }
  void get_termbranch(void){
    // Get lists of terminals and branchpoints
    terminals.clear();
    branchpoints.clear();
    for(int n=0; n < edgelist.size(); n++){
      int vno = edgelist[n].second;
      int nnbrs = nmasknbrs(vno);
      if(nnbrs == 1) terminals.push_back(vno);
      if(nnbrs  > 2) branchpoints.push_back(vno);
      // nnbrs==0 is possible, maybe?
    }
    printf("Found %d terminals and %d branchpoints\n",(int)terminals.size(),(int)branchpoints.size());
  };
  int get_branches(void){
    // This was meant as a segmentation of the skeleton into individual segments ("branches"). This 
    // code works but only on spurs. More needs to be done, but had to work on other things.
    branches.clear();
    std::vector<int> hitmap(m_surf->nvertices,0);// only makes a difference for terminals
    for(int n=0; n < terminals.size(); n++){
      int vno = terminals[n];
      if(hitmap[vno]) continue;
      branchstruct branch;
      branch.vtxlist.push_back(vno);
      hitmap[vno] = 1;
      printf("terminal %d %d %d ",n,(int)branches.size()+1,vno);
      fflush(stdout);
      // Move out from this vertex, finding the next active but unhit neighbor 
      while(1){
	// If this vertex is a terminal or branchpoint, so stop
	// searching.  Note that the same terminals/branchpoints will
	// be represented across multiple branches
	int nmnbrs = nmasknbrs(vno);
	if(branch.vtxlist.size()>1 && nmnbrs != 2) break;
	// Otherwise, look through the neighbors for the next one (has to be there)
	std::vector<int> nbrlist = vtxnbrlist[vno];
	int found = 0;
	for(int k=0; k < nbrlist.size(); k++){
	  int nbrvno = nbrlist[k];
	  if(!MRIgetVoxVal(m_mask,nbrvno,0,0,0)) continue;
	  if(hitmap[nbrvno]) continue;
	  // If it gets here, then this vetex is the next step (there
	  // can be one and only one).  Add it to the list and update
	  // the hit map.
	  branch.vtxlist.push_back(nbrvno);
	  hitmap[nbrvno] = 1;
	  //int len = (int)branch.vtxlist.size();
	  //printf("  Terminal %d %6d nm=%d nnbrs=%d k=%d nbrvno=%6d m=%d len=%d\n",
	  //	 n,vno,nmnbrs,(int)nbrlist.size(),k,nbrvno,(int)MRIgetVoxVal(m_mask,nbrvno,0,0,0),len);
	  //fflush(stdout);
	  found = 1;
	  vno = nbrvno;
	  break;
	} // loop over neighbors
	if(!found) break;
      }// while(!done) getting members of the branch
      printf(" len = %5d\n",(int)branch.vtxlist.size()); fflush(stdout);
      branches.push_back(branch);
    }
    printf("Found %d branches\n",(int)branches.size());
    return((int)branches.size());
  }
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
  int write_label(char *fname, char *subject){
    // Write current edge list vertices as a surface label. When the
    // procedure is done, this will be the skeleton.
    LABEL *label = LabelAlloc(edgelist.size()+1, subject, "skeleton");
    label->n_points = edgelist.size();
    for(int n=0; n < edgelist.size(); n++){
      int vno = edgelist[n].second;
      label->lv[n].vno = vno;
      label->lv[n].x = m_surf->vertices[vno].x;
      label->lv[n].y = m_surf->vertices[vno].y;
      label->lv[n].z = m_surf->vertices[vno].z;
      label->lv[n].stat = surfvals[vno];
    }
    int err = LabelWrite(label,fname);
    return(err);
  }
  int cluster(void){
    printf("Clusterizing %d\n",nkeep);
    int NClusters=0;
    MRIScopyMRI(m_surf,m_mask,0,"val"); // overwrites the val field
    SURFCLUSTERSUM *scs = sclustMapSurfClusters(m_surf,0.5,-1,+1,0,&NClusters,NULL,NULL);
    printf("Found %d clusters\n",NClusters);
    double cmaxsize = sclustMaxClusterArea(scs, NClusters);
    printf("Max cluster size %lf\n",cmaxsize);
    SURFCLUSTERSUM *scs2 = SortSurfClusterSum(scs, NClusters);
    free(scs);
    scs = scs2;
    sclustAnnot(m_surf, NClusters); // create an annotation for convenience
    for(int vno=0; vno < m_surf->nvertices; vno++){
      int cno = m_surf->vertices[vno].undefval;
      for(int k=0; k < nkeep; k++){
	double v=0;
	if(scs[k].clusterno == cno) v=1;
	MRIsetVoxVal(m_mask,vno,0,0,0,v);
      }
    }
    return(0);
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
    if(nkeep>0) cluster();

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
char *surfvalspath=NULL;
char *surfpath=NULL, *spherepath=NULL, *outmaskpath=NULL,*outsurfvalspath=NULL,*pointsetpath=NULL,*labelpath=NULL,*tractpath=NULL;
char *outdir=NULL;
SurfSkeleton sk;
double scale = -1; // gyrus, crowns
char *subject=NULL, *hemi=NULL, *surfname=NULL;
int inputtype=0;

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

  if(outdir){
    err = mkdir(outdir,0777);
    if (err != 0 && errno != EEXIST) {
      printf("ERROR: creating directory %s\n",outdir);
      perror(NULL);
      return(1);
    }
  }

  sk.m_surf = MRISread(surfpath);
  if(!sk.m_surf) exit(1);

  if(inputtype == 1){
    MRI *mri = MRIread(surfvalspath);
    if(mri==NULL) exit(1);
    sk.mri2surfvals(mri);
    MRIfree(&mri);
  }
  if(inputtype == 2){  
    sk.m_sphere = MRISread(spherepath);
    sk.load_curv_nonmaxsup(scale);
    if(!sk.m_sphere) exit(1);
  }
  if(inputtype == 3){  
    sk.load_k1(scale);
  }

  sk.skeletonize();

  if(outmaskpath) {
    err = MRIwrite(sk.m_mask,outmaskpath);
    if(err) exit(err);
  }
  if(outsurfvalspath){
    err = sk.write_surfvals(outsurfvalspath);
    if(err) exit(err);
  }
  if(labelpath){
    err = sk.write_label(labelpath,subject);
    if(err) exit(err);
  }
  if(pointsetpath){
    // returns 1 always
    sk.write_pointset(pointsetpath);
  }
  if(outdir){
    char tmpstr[2000];
    sprintf(tmpstr,"%s/cluster.mgz",outdir);
    printf("Writing annotation %s\n",tmpstr);
    MRISwriteAnnotation(sk.m_surf, tmpstr);
  }

  //sk.get_termbranch();
  //sk.get_branches();

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
    else if (!strcasecmp(option, "--gyrus")  || !strcasecmp(option, "--crown"))  scale = -1; 
    else if (!strcasecmp(option, "--sulcus") || !strcasecmp(option, "--fundus")) scale = +1; 

    else if(!strcasecmp(option, "--surfvals")){
      if (nargc < 1) CMDargNErr(option,1);
      surfvalspath = pargv[0];
      inputtype = 1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--curv-nonmaxsup")) inputtype = 2;
    else if (!strcasecmp(option, "--k1"))             inputtype = 3;

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
    else if (!strcasecmp(option, "--outdir")){
      if(nargc < 1) CMDargNErr(option,1);
      outdir = pargv[0];
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
    else if(!strcasecmp(option, "--cluster")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&sk.nkeep);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--fwhm")){
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&sk.m_fwhm);
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--out-surfvals")){
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

  if(inputtype == 0){
    printf("ERROR: input type not specified, use --surfvals, --curv-nonmaxsup, or --k1\n");
    exit(1);
  }

  if(outdir){
    char tmpstr[2000];
    if(!outmaskpath){
      sprintf(tmpstr,"%s/skeleton.mgz",outdir);
      outmaskpath = strcpyalloc(tmpstr);
    }
    if(!pointsetpath){
      sprintf(tmpstr,"%s/skeleton.json",outdir);
      pointsetpath = strcpyalloc(tmpstr);
    }
    if(!outsurfvalspath){
      sprintf(tmpstr,"%s/surfvals.mgz",outdir);
      outsurfvalspath = strcpyalloc(tmpstr);
    }
    if(!labelpath){
      sprintf(tmpstr,"%s/skeleton.label",outdir);
      labelpath = strcpyalloc(tmpstr);
    }
    if(!tractpath){
      sprintf(tmpstr,"%s/skeleton.tract",outdir);
      tractpath = strcpyalloc(tmpstr);
    }
  }

  // atlasMesh is required
  if(surfpath == NULL){
    printf("ERROR: must spec --surf\n");
    exit(1);
  }
  if(inputtype == 2 && spherepath == NULL){
    printf("ERROR: must spec --sphere with --curv-nonmaxsup\n");
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
