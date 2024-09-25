#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include "mri.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "timer.h"
#include "utils.h"
#include "error.h"
#include "diag.h"
//#include "dmatrix.h"
//#include "surfgrad.h"
//#include "fsglm.h"
//#include "mrisurf_metricProperties.h"

#define export
#define ENS_PRINT_INFO
#define ENS_PRINT_WARN
#include <armadillo>
//#include <ensmallen.hpp>
#include "romp_support.h"

//using namespace ens;
//using namespace arma;

int GetP(MRIS *surf, int vno)
{
  VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vno]);
  VERTEX *vtx = &(surf->vertices[vno]);

  FILE *fp = fopen("P.txt","w");
  fprintf(fp,"%20.18f %20.18f %20.18f\n ",vtx->x,vtx->y,vtx->z);
  int vprev=vno,vnbr=0;
  vnbr = vtop->v[0];
  std::vector<int> vlist;
  vlist.push_back(vno);
  for(int n=0; n < vtop->vnum; n++){
    vlist.push_back(vnbr);
    VERTEX *nbrvtx = &surf->vertices[vnbr];
    fprintf(fp,"%20.18f %20.18f %20.18f\n",nbrvtx->x,nbrvtx->y,nbrvtx->z);
    if(n == vtop->vnum-1) break;
    for(int k=0; k < vtop->vnum; k++){    
      int fno = vtop->f[k];
      FACE *f = &surf->faces[fno];
      int nhits = 0, cnothit = -1;
      for(int c=0; c < 3; c++){
	if(f->v[c] == vno || f->v[c] == vnbr) nhits++;
	else cnothit = c;
      }
      if(nhits == 2){
	int vnohit = f->v[cnothit];
	int mhits = 0;
	for(int m=0; m < vlist.size(); m++){
	  if(vlist[m] == vnohit) mhits++;
	}
	if(mhits == 0){
	  vprev=vnbr;
	  vnbr = vnohit;
	  break;
	}
      }
    }
  }
  for(int n=vtop->num; n < vtop->vtotal; n++){
    vnbr = vtop->v[n];
    VERTEX *nbrvtx = &surf->vertices[vnbr];
    fprintf(fp,"%20.18f %20.18f %20.18f\n",nbrvtx->x,nbrvtx->y,nbrvtx->z);
  }
  fclose(fp);
  return(0);
}

//======================================================================
arma::colvec FSdx = {1,0,0};
arma::colvec FSdy = {0,1,0};
arma::colvec FSdz = {0,0,1};
arma::mat FS_J_v1_nA = {{0,-1,0},{1,0,0},{0,0,0}};
arma::mat FS_J_v1_nB = {{0,0,-1},{0,0,0},{1,0,0}};

class FSface{ //----------------
public:
  MRIS *surf;
  int fno;
  FACE *f;
  arma::mat P = arma::mat(3,3);
  arma::colvec norm, v0, v1, m;
  double L;
  int RecomputeJ=1;
  std::array<arma::mat,3> J_n_p;
  int DoDelta=0, nDelta=0, cDelta=0;
  double Delta = 0.0001;
  arma::colvec computeNorm(arma::mat &vxyz);
  int computeNormJ(void);
};

class FSvertexcell { 
public:
  MRIS *surf;
  std::vector<FSface> faces;
  int vno;
  VERTEX_TOPOLOGY *vtop=NULL;
  int npoints=0, nnnbrs=0;
  arma::colvec norm, t, v1, v2, e1, e2;
  double L, L1, L2;
  int cuse=0;
  std::vector<arma::mat> J_n_p, J_e1_p, J_e2_p;
  int computeNorm(std::vector<FSface> const &faces);
  int computeNormJ(std::vector<FSface> const &faces);
};

class FSmeancurv  { //----------------
public:
  double H;
  arma::mat J_H_p;
  FSvertexcell vc;
  arma::mat P,M,Q,Qt,D,F,beta;
  arma::colvec  u, v, w;
  int compute(arma::mat &vxyz);
  int computeJ(void);
}; 
int FSmeancurv::computeJ(void)
{
  //arma::sp_mat dMdx,dMdy,dMdz;
 arma::mat dMdx,dMdy,dMdz;
  J_H_p.zeros(vc.npoints,3);
  for(int n=0; n < vc.npoints; n++){
    dMdx.zeros(vc.npoints,3);
    dMdy.zeros(vc.npoints,3);
    dMdz.zeros(vc.npoints,3);
    if(n==0){
      for(int k=1; k < vc.npoints; k++){
	dMdx(k,0) = -1;
	dMdy(k,1) = -1;
	dMdz(k,2) = -1;
      }
    }
    else{
      dMdx(n,0) = 1;
      dMdy(n,1) = 1;
      dMdz(n,2) = 1;
    }

    arma::colvec dudx, dvdx, dwdx;
    dudx = dMdx*vc.e1;
    dvdx = dMdx*vc.e2;
    dwdx = dMdx*vc.norm;
    if(n < vc.nnnbrs+1){
      dudx += M*vc.J_e1_p[n].col(0);
      dvdx += M*vc.J_e2_p[n].col(0);
      dwdx += M*vc.J_n_p[n].col(0);
    }
    arma::mat dQdx = arma::join_horiz(arma::join_horiz(2*u%dudx,(2*dudx%v+2*dvdx%u)),2*v%dvdx);
    arma::mat dQdxt = arma::trans(dQdx);
    arma::mat dDdx = dQdxt*Q + Qt*dQdx;
    arma::mat dFdx = -F*dDdx*F;
    arma::mat dbetadx = dFdx*(Qt*w) + F*(dQdxt*w) + F*(Qt*dwdx);
    double J_H_x = dbetadx(0)+dbetadx(2);

    arma::colvec dudy, dvdy, dwdy;
    dudy = dMdy*vc.e1;
    dvdy = dMdy*vc.e2;
    dwdy = dMdy*vc.norm;
    if(n < vc.nnnbrs+1){
      dudy += M*vc.J_e1_p[n].col(1);
      dvdy += M*vc.J_e2_p[n].col(1);
      dwdy += M*vc.J_n_p[n].col(1);
    }
    arma::mat dQdy = arma::join_horiz(arma::join_horiz(2*u%dudy,(2*dudy%v+2*dvdy%u)),2*v%dvdy);
    arma::mat dQdyt = arma::trans(dQdy);
    arma::mat dDdy = dQdyt*Q + Qt*dQdy;
    arma::mat dFdy = -F*dDdy*F;
    arma::mat dbetady = dFdy*(Qt*w) + F*(dQdyt*w) + F*(Qt*dwdy);
    double J_H_y = dbetady(0)+dbetady(2);

    arma::colvec dudz, dvdz, dwdz;
    dudz = dMdz*vc.e1;
    dvdz = dMdz*vc.e2;
    dwdz = dMdz*vc.norm;
    if(n < vc.nnnbrs+1){
      dudz += M*vc.J_e1_p[n].col(2);
      dvdz += M*vc.J_e2_p[n].col(2);
      dwdz += M*vc.J_n_p[n].col(2);
    }
    arma::mat dQdz = arma::join_horiz(arma::join_horiz(2*u%dudz,(2*dudz%v+2*dvdz%u)),2*v%dvdz);
    arma::mat dQdzt = arma::trans(dQdz);
    arma::mat dDdz = dQdzt*Q + Qt*dQdz;
    arma::mat dFdz = -F*dDdz*F;
    arma::mat dbetadz = dFdz*(Qt*w) + F*(dQdzt*w) + F*(Qt*dwdz);
    double J_H_z = dbetadz(0)+dbetadz(2);

    J_H_p(n,0) = J_H_x;
    J_H_p(n,1) = J_H_y;
    J_H_p(n,2) = J_H_z;

  } // end loop over points
  return(0);
}
int FSmeancurv::compute(arma::mat &vxyz)
{
  int vno = vc.vno;
  VERTEX_TOPOLOGY *vtop = vc.vtop;
  VERTEX *vtx = &(vc.surf->vertices[vno]);
  int npoints = vtop->vtotal+1;

  P.zeros(npoints,3);
  //P(0,0) = vtx->x;
  //P(0,1) = vtx->y;
  //P(0,2) = vtx->z;
  P(0,0) = vxyz(vno,0);
  P(0,1) = vxyz(vno,1);
  P(0,2) = vxyz(vno,2);
  for(int n=1; n < npoints; n++){
    int nbrvno = vtop->v[n-1];
    P(n,0) = vxyz(nbrvno,0);
    P(n,1) = vxyz(nbrvno,1);
    P(n,2) = vxyz(nbrvno,2);
    //VERTEX *nbrvtx = &(vc.surf->vertices[nbrvno]);
    //P(n,0) = nbrvtx->x;
    //P(n,1) = nbrvtx->y;
    //P(n,2) = nbrvtx->z;
  }
  M = P - arma::repmat(P.row(0),npoints,1);

  u = M*vc.e1;
  v = M*vc.e2;
  w = M*vc.norm;

  Q = arma::join_horiz(arma::join_horiz(u%u,2*u%v),v%v);
  Qt = arma::trans(Q);
  D = Qt*Q;
  F = arma::inv(D);
  beta = F*Qt*w;
  H = (beta(0)+beta(2));
  vtx->curv = H;
  //printf("meancurv %d %g\n",vno,km);

  return(0);
}


class FSsurf //####################################################
{
public:
  MRIS *surf;
  std::vector<FSface> faces;
  std::vector<FSvertexcell> vertices;
  std::vector<FSmeancurv> mcurvs;
  arma::mat vxyz;
  int CopyVXYZ(void){
    vxyz.clear();
    vxyz.set_size(surf->nvertices,3);
    for(int vno=0; vno < surf->nvertices; vno++){
      VERTEX *vtx = &surf->vertices[vno];
      vxyz(vno,0) = vtx->x;
      vxyz(vno,1) = vtx->y;
      vxyz(vno,2) = vtx->z;
    }
    return(0);
  }
  int ComputeFaces(int DoJ, int vno0=-1){
    if(faces.size() != surf->nfaces){
      faces.clear();
      faces.resize(surf->nfaces);
    }
    if(vno0 > -1){
      VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vno0]);
      for(int n=0; n < vtop->vnum; n++){
	int fno = vtop->f[n];
	faces[fno].surf = surf;
	faces[fno].fno = fno;
	faces[fno].computeNorm(vxyz);
	if(DoJ) faces[fno].computeNormJ();
      }
      return(0);
    }
    #ifdef HAVE_OPENMP
     #pragma omp parallel for 
    #endif
    for(int fno=0; fno < surf->nfaces; fno++){
      faces[fno].surf = surf;
      faces[fno].fno = fno;
      faces[fno].computeNorm(vxyz);
      if(DoJ) faces[fno].computeNormJ();
    }
    return(0);
  }
  int ComputeVertices(int DoJ, int vno0=-1){
    if(vertices.size() != surf->nvertices){
      vertices.clear(); vertices.resize(surf->nvertices);
    }
    int vnostart=0;
    int vnostop=surf->nvertices;
    if(vno0 > -1){vnostart=vno0; vnostop=vno0+1;}
    #ifdef HAVE_OPENMP
     #pragma omp parallel for 
    #endif
    for(int vno=vnostart; vno < vnostop; vno++){
      vertices[vno].surf = surf;
      vertices[vno].vno = vno;
      vertices[vno].computeNorm(faces);
      if(DoJ) vertices[vno].computeNormJ(faces);
    }
    return(0);
  }
  int ComputeMeanCurvs(int DoJ, int vno0=-1){
    if(mcurvs.size() != surf->nvertices){
      mcurvs.clear();
      mcurvs.resize(surf->nvertices);
    }
    int vnostart=0;
    int vnostop=surf->nvertices;
    if(vno0 > -1){vnostart=vno0; vnostop=vno0+1;}
    #ifdef HAVE_OPENMP
     #pragma omp parallel for 
    #endif
    for(int vno=vnostart; vno < vnostop; vno++){
      mcurvs[vno].vc = vertices[vno]; // be good to have a pointer
      mcurvs[vno].compute(vxyz);
      if(DoJ) mcurvs[vno].computeJ();
    }
    return(0);
  }
};
arma::colvec FSface::computeNorm(arma::mat &vxyz)
{
  f = &surf->faces[fno]; 
  for(int r=0; r<3; r++){
    int vno = f->v[r];
    for(int k=0; k<3; k++){
      P(k,r) = vxyz(vno,k);
    }
    //VERTEX *vtx = &surf->vertices[vno];
    //P(0,r)= vtx->x;
    //P(1,r)= vtx->y;
    //P(2,r)= vtx->z;
  }
  double psave=0;
  if(DoDelta){
    psave = P(cDelta,nDelta);
    P(cDelta,nDelta) += Delta;
  }
  v0 = P.col(2)-P.col(0);
  v1 = P.col(1)-P.col(2);
  m = arma::cross(v1,v0);
  L = arma::norm(m);
  norm = m/L;
  //float snorm[3];
  //mrisNormalFace(surf, fno, 1, snorm);
  RecomputeJ=1;
  if(DoDelta) P(cDelta,nDelta)=psave; 
  return(norm);
}

int FSface::computeNormJ(void)
{
  if(!RecomputeJ) return(0);

  arma::mat J_m_v0x = arma::cross(v1,FSdx);
  arma::mat J_m_v0y = arma::cross(v1,FSdy);
  arma::mat J_m_v0z = arma::cross(v1,FSdz);
  arma::mat J_m_v0 = arma::join_horiz(arma::join_horiz(J_m_v0x,J_m_v0y),J_m_v0z);

  arma::mat J_m_v1x = arma::cross(FSdx,v0);
  arma::mat J_m_v1y = arma::cross(FSdy,v0);
  arma::mat J_m_v1z = arma::cross(FSdz,v0);
  arma::mat J_m_v1 = arma::join_horiz(arma::join_horiz(J_m_v1x,J_m_v1y),J_m_v1z);

  arma::mat J_m_p0 = -J_m_v0;
  arma::mat J_m_p1 =  J_m_v1;
  arma::mat J_m_p2 =  J_m_v0 - J_m_v1;

  // J_v0_p0=-I; J_v1_p0=0; J_v0_p1=0; J_v1_p1=I; J_v0_p2=I; J_v1_p2=-I;
  // J_L_pi = (m'*J_m_pi)/(L*L)
  // J_n_pi = J_m_pi/L - (m*J_L_pi)/L

  // Change in face normal wrt a change in corner position
  arma::mat J_L_p0 = (arma::trans(m)*J_m_p0)/L;
  J_n_p[0] = J_m_p0/L - (m*J_L_p0)/(L*L);
  arma::mat J_L_p1 = (arma::trans(m)*J_m_p1)/L;
  J_n_p[1] = J_m_p1/L - (m*J_L_p1)/(L*L);
  arma::mat J_L_p2 = (arma::trans(m)*J_m_p2)/L;
  J_n_p[2] = J_m_p2/L - (m*J_L_p2)/(L*L);

  return(0);
}

int FSvertexcell::computeNorm(std::vector<FSface> const &faces)
{
  vtop = &(surf->vertices_topology[vno]);
  npoints = vtop->vtotal+1;
  nnnbrs = vtop->vnum;
  t.zeros(3);
  for(int k=0; k<vtop->vnum; k++){
    int fno = vtop->f[k];
    t += faces[fno].norm;
  }
  L = arma::norm(t);
  norm = t/L;

  v1.zeros(3);
  if(fabs(norm(1))>fabs(norm(2))){
    // v1 = [-y x 0]
    v1(0) = -norm(1);
    v1(1) =  norm(0);
    v1(2) =  0;
    cuse = 1;
  }
  else{
    // v1 = [-z 0 x]
    v1(0) = -norm(2);
    v1(1) =  0;      
    v1(2) =  norm(0);
    cuse = 2;
  }
  L1 = arma::norm(v1);
  e1 = v1/L1;

  v2 = arma::cross(norm,e1);
  L2 = arma::norm(v2);
  e2 = v2/L2;

  //printf("e1 %g %g %g  e2 %g%g %g\n",e1(0),e1(1),e1(2),e1(0),e2(1),e2(2));

  //VERTEX *vtx = &(surf->vertices[vno]);
  //if(fabs(vtx->nx-norm(0))>.0001 || fabs(vtx->ny-norm(1))>.0001 || fabs(vtx->nz-norm(2))>.0001)
  //printf("%d  %g %g %g   %g %g %g\n",vno,vtx->nx,vtx->ny,vtx->nz,norm(0),norm(1),norm(2));
  return(0);
}

int FSvertexcell::computeNormJ(std::vector<FSface> const &faces)
{
  if(J_n_p.size() != nnnbrs+1){
    J_n_p.clear();
    J_n_p.resize(nnnbrs+1);
  }
  if(J_e1_p.size() != nnnbrs+1){
    J_e1_p.clear();
    J_e1_p.resize(nnnbrs+1);
  }
  if(J_e2_p.size() != nnnbrs+1){
    J_e2_p.clear();
    J_e2_p.resize(nnnbrs+1);
  }

  // Compute J_e?_n - grad  of e? wrt the normal as this 
  // does not depend on the point
  arma::mat J_e1_n;
  arma::rowvec J_L1_n;
  if(cuse == 1){
    J_L1_n = (arma::trans(v1)*FS_J_v1_nA)/L1;
    J_e1_n = FS_J_v1_nA/L1 - (v1*J_L1_n)/(L1*L1);
  } else {
    J_L1_n = (arma::trans(v1)*FS_J_v1_nB)/L1;
    J_e1_n = FS_J_v1_nB/L1 - (v1*J_L1_n)/(L1*L1);
  }

  arma::mat J_e2_n;
  arma::mat J_v2_nx = arma::cross(FSdx,e1) + arma::cross(norm,J_e1_n.col(0));
  arma::mat J_v2_ny = arma::cross(FSdy,e1) + arma::cross(norm,J_e1_n.col(1));
  arma::mat J_v2_nz = arma::cross(FSdz,e1) + arma::cross(norm,J_e1_n.col(2));
  arma::mat J_v2_n = arma::join_horiz(arma::join_horiz(J_v2_nx,J_v2_ny),J_v2_nz);
  arma::rowvec J_L2_n = (arma::trans(v2)*J_v2_n)/L2;
  J_e2_n = J_v2_n/L2 - (v2*J_L2_n)/(L2*L2);


  for(int n=0; n < nnnbrs+1; n++){
    int nvno;
    if(n == 0) nvno=vno;
    else       nvno = vtop->v[n-1];
    arma::mat J_t_p;
    J_t_p.zeros(3,3);
    int nhits = 0;
    for(int k=0; k<vtop->vnum; k++){
      int fno = vtop->f[k];
      FACE *vface = &surf->faces[fno];
      int m;
      for(m=0; m<3; m++) if(vface->v[m] == nvno) break;
      if(m==3) continue; // this vertex not in this face
      nhits++;
      J_t_p += faces[fno].J_n_p[m];
    }
    // nhits should be 2 or vnum
    arma::mat J_L_p = (arma::trans(t)*J_t_p)/L;
    J_n_p[n]  = J_t_p/L - (t*J_L_p)/(L*L);
    J_e1_p[n] = J_e1_n*J_n_p[n];
    J_e2_p[n] = J_e2_n*J_n_p[n];
  }

  if(0){
  printf("cuse = %d  L1 = %g\n",cuse,L1);
  FS_J_v1_nA.print("J_v1_nA = [");
  FS_J_v1_nB.print("J_v1_nB = [");
  J_L1_n.print("J_L1_n = [");
  J_e1_n.print("J_e1_n = [");
  J_e2_n.print("J_e2_n = [");
  J_n_p[1].print("J_n_p2=[");
  J_e1_p[1].print("J_e1_p2=[");
  J_e2_p[1].print("J_e2_p2=[");
  norm.print("norm");
  e1.print("e1");
  e2.print("e2");
  }
  //printf("%d %d %d %d\n",n,nvno,vtop->vnum,nhits);
  //printf("J_t_p\n");
  //J_t_p.print();
  //printf("J_n_p\n");
  //J_n_p.print();

  return(0);
}
#if 0

double FSsurf::test_face_NormJ(double thresh)
{
  // Note that this needs to be done at double precision. Applying
  // the delta to v->{xyz} does not work beause v->{xyz} is single precision
  FSsurf::face myface;
  FSsurf::face myface2;
  myface.surf = surf;
  myface2.surf = surf;
  myface2.DoDelta = 1;
  // Tricky to set the delta
  myface2.Delta = 1e6*arma::datum::eps; // 10e-8 might be ok
  double maxd = 0;
  int fnomax=0;
  for(int fno=0; fno < surf->nfaces; fno++){
    myface.fno = fno;
    myface.computeNorm();
    myface.computeNormJ();
    myface2.fno = fno;
    for(int n=0; n<3; n++){
      //int vno = myface.f->v[n];
      myface2.nDelta = n;
      for(int c=0; c<3; c++){
	myface2.cDelta = c;
	//double saveval;
	//if(c==0) {saveval = vtx->x; vtx->x = (vtx->x + Delta);}
	//if(c==1) {saveval = vtx->y; vtx->y = (vtx->y + Delta);}
	//if(c==2) {saveval = vtx->z; vtx->z = (vtx->z + Delta);}
	myface2.computeNorm();
	//if(c==0) vtx->x = saveval;
	//if(c==1) vtx->y = saveval;
	//if(c==2) vtx->z = saveval;
	arma::mat d = (myface2.norm-myface.norm)/myface2.Delta;
	arma::mat e = arma::abs(d-myface.J_n_p[n].col(c));
	double maxevnc = e.max();
	if(maxd < maxevnc) {maxd = maxevnc;fnomax=fno;}
	if(maxevnc > thresh){ 
	  printf("%d %d %d maxvnc %g %g %g\n",fno,n,c,maxevnc,myface.L,myface2.Delta);
	  myface.P.print("P---------------");
	  myface.J_n_p[n].col(c).print("J_n_p---------------");
	  d.print("numJ-------------");
	  e.print("e-------------");
	  myface.norm.print("norm-----------");
	  myface.P.print("P-------------");
	  myface.v0.print("v0-------------");
	  myface.v1.print("v1-------------");
	}
      }
    }
  }
  printf("maxd %g fnomax %d\n",maxd,fnomax);
  return(maxd);
}

#endif

//MAIN ----------------------------------------------------
int main(int argc, char **argv) 
{
#ifdef HAVE_OPENMP
  omp_set_num_threads(20);
#endif
  //Gdiag_no = 72533;
  MRIS *surf;
  //MRI *mri=NULL; // *mri2;
  int vno;

  Timer mytimer;

  sscanf(argv[3],"%d\n",&vno);

  surf = MRISread(argv[1]);
  MRISsetNeighborhoodSizeAndDist(surf,2);

  //for(int n=0; n < 5; n++){
  //  double t0 = mytimer.seconds();
  //  mrisComputeTangentPlanes(surf);
  //  MRIScomputeSecondFundamentalFormThresholded(surf,-1);
  //  printf("dt  %6.4f\n",mytimer.seconds()-t0);
  //}

  VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vno]);
  VERTEX *vtx = &(surf->vertices[vno]);
  printf("vno %d vnum %d vtotal %d\n",vno,vtop->vnum,vtop->vtotal);

  GetP(surf,vno);// saves P.txt
  FILE *fp = fopen("nnnbrs.dat","w");
  fprintf(fp,"%d\n",vtop->vnum);
  fclose(fp);
  //printf("FS: k1 %g k2 %g H %g\n",vtx->k1,vtx->k2,vtx->H);

  FSsurf fs;
  fs.surf = surf;
  fs.CopyVXYZ();
  fs.ComputeFaces(1,vno);
  fs.ComputeVertices(1,vno);
  fs.ComputeMeanCurvs(1,vno);
  //fs.mcurvs[vno].J_H_p.print("J_H_p = [");
  //printf("];\n");
  //fs.mcurvs[vno].J_H_p.save("J_H_p.txt",arma::raw_ascii);

  FSsurf fs2;
  fs2.surf = surf;
  fs2.CopyVXYZ();
  fs2.ComputeFaces(0,vno);
  fs2.ComputeVertices(0,vno);
  fs2.ComputeMeanCurvs(0,vno);

  double t = mytimer.seconds();
  double Delta = 1e-10;
  for(vno=0; vno < 100; vno++){
    fs.ComputeFaces(1,vno);
    fs.ComputeVertices(1,vno);
    fs.ComputeMeanCurvs(1,vno);
    VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vno]);
    arma::mat J(vtop->vtotal+1,3);
    J.zeros();
    for(int n=0; n < vtop->vtotal+1; n++){
      int vnbrno=0;
      if(n==0) vnbrno = vno;
      else     vnbrno = vtop->v[n-1];
      for(int c=0; c < 3; c++){
	double psave = fs2.vxyz(vnbrno,c);
	fs2.vxyz(vnbrno,c) += Delta;
	fs2.ComputeFaces(0,vno);
	fs2.ComputeVertices(0,vno);
	fs2.ComputeMeanCurvs(0,vno);
	fs2.vxyz(vnbrno,c) = psave;
	J(n,c) = (fs2.mcurvs[vno].H-fs.mcurvs[vno].H)/Delta;
      }
    }
    //J.print("Jnum = [");
    //printf("];\n");
    //J.save("Jnum.txt",arma::raw_ascii);
    arma::mat Jdiff = J-fs.mcurvs[vno].J_H_p;
    arma::mat Jdiffabs = arma::abs(Jdiff);
    //Jdiff.print("Jdiff = [");
    printf("%d max Jdiff %g\n",vno,Jdiffabs.max());
  }
  printf("dt = %6.4f\n",mytimer.seconds()-t);

  //mri = MRIcopyMRIS(NULL, surf, 0, "curv"); // start at z to autoalloc
  //MRIwrite(mri,argv[2]);
  //exit(0);

  for(int n=0; n < 5; n++){
    double t = mytimer.seconds();
    fs.CopyVXYZ();
    fs.ComputeFaces(0);
    fs.ComputeVertices(0);
    fs.ComputeMeanCurvs(0);
    printf("H %d %6.4f\n",n,mytimer.seconds()-t);
  }
  for(int n=0; n < 5; n++){
    double t = mytimer.seconds();
    fs.CopyVXYZ();
    fs.ComputeFaces(1);
    fs.ComputeVertices(1);
    fs.ComputeMeanCurvs(1);
    printf("dH %d %6.4f\n",n,mytimer.seconds()-t);
  }

  printf("%g %g %g\n",vtx->nx,vtx->ny,vtx->nz);
  printf("%g %g %g\n",fs.vertices[vno].norm(0),fs.vertices[vno].norm(1),fs.vertices[vno].norm(2));
  printf("%g %g %g\n",vtx->nx-fs.vertices[vno].norm(0),vtx->ny-fs.vertices[vno].norm(1),vtx->nz-fs.vertices[vno].norm(2));


  //Curv2(surf, vno);

  exit(0);

  printf("vno = %d  %7.4f %7.4f %7.4f\n",vno,vtx->x,vtx->y,vtx->z);

  fp = fopen("faces.dat","w");
  for(int n=0; n<vtop->vnum; n++){
    //printf("%2d %6d %6d  ",n,vtop->v[n],vtop->f[n]);
    FACE *f = &surf->faces[vtop->f[n]];
    //for(int k=0; k < 3; k++)  printf("%6d ",f->v[k]);
    for(int k=0; k < 3; k++){
      VERTEX *kvtx = &surf->vertices[f->v[k]];
      fprintf(fp," %7.4f %7.4f %7.4f ",kvtx->x,kvtx->y,kvtx->z);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);


  printf("Nbr %6.4f\n",mytimer.seconds());

  printf("Tan %6.4f sec\n",mytimer.seconds());

  //Curv2(surf, vno);
  //Curv3(surf, vno);


  exit(0);
}

