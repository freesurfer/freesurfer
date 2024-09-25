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
//#define ENS_PRINT_INFO
//#define ENS_PRINT_WARN
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
  arma::mat P = arma::mat(3,3); //col vects
  arma::colvec norm, v0, v1, m;
  double L;
  int RecomputeJ=1;
  std::array<arma::mat,3> J_n_p;
  arma::colvec computeNorm(arma::mat &vxyz);
  int computeNormJ(void);
  int computeNormJslow(void);
};
arma::mat computeJ_v_cross_xyz(arma::colvec const &v);

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
  int computeMaybeFaster(arma::mat &vxyz);
  int computeJ(void);
  int computeJslow(void);
}; 

int FSmeancurv::computeJ(void)
{
  J_H_p.zeros(vc.npoints,3);
  for(int n=0; n < vc.npoints; n++){
    arma::colvec dudx, dudy, dudz,  dvdx, dvdy, dvdz,  dwdx, dwdy, dwdz;
    if(n==0){
      dudx.zeros(vc.npoints); dudx.fill(-vc.e1(0)); dudx(0) = 0;
      dudy.zeros(vc.npoints); dudy.fill(-vc.e1(1)); dudy(0) = 0;
      dudz.zeros(vc.npoints); dudz.fill(-vc.e1(2)); dudz(0) = 0;
      // Note: either e1(1) or e1(2) will be 0 so dont' have to fill
      dvdx.zeros(vc.npoints); dvdx.fill(-vc.e2(0)); dvdx(0) = 0;
      dvdy.zeros(vc.npoints); dvdy.fill(-vc.e2(1)); dvdy(0) = 0;
      dvdz.zeros(vc.npoints); dvdz.fill(-vc.e2(2)); dvdz(0) = 0;
      dwdx.zeros(vc.npoints); dwdx.fill(-vc.norm(0)); dwdx(0) = 0;
      dwdy.zeros(vc.npoints); dwdy.fill(-vc.norm(1)); dwdy(0) = 0;
      dwdz.zeros(vc.npoints); dwdz.fill(-vc.norm(2)); dwdz(0) = 0;
    }
    else{
      dudx.zeros(vc.npoints); dudx(n) = vc.e1(0);
      dudy.zeros(vc.npoints); dudy(n) = vc.e1(1);
      dudz.zeros(vc.npoints); dudz(n) = vc.e1(2);
      dvdx.zeros(vc.npoints); dvdx(n) = vc.e2(0);
      dvdy.zeros(vc.npoints); dvdy(n) = vc.e2(1);
      dvdz.zeros(vc.npoints); dvdz(n) = vc.e2(2);
      dwdx.zeros(vc.npoints); dwdx(n) = vc.norm(0);
      dwdy.zeros(vc.npoints); dwdy(n) = vc.norm(1);
      dwdz.zeros(vc.npoints); dwdz(n) = vc.norm(2);
    }
    if(n < vc.nnnbrs+1){
      // one of the rows of J_e1_p will be 0
      dudx += M*vc.J_e1_p[n].col(0);
      dvdx += M*vc.J_e2_p[n].col(0);
      dwdx += M*vc.J_n_p[n].col(0);
      dudy += M*vc.J_e1_p[n].col(1);
      dvdy += M*vc.J_e2_p[n].col(1);
      dwdy += M*vc.J_n_p[n].col(1);
      dudz += M*vc.J_e1_p[n].col(2);
      dvdz += M*vc.J_e2_p[n].col(2);
      dwdz += M*vc.J_n_p[n].col(2);
    }
    arma::mat dQdx = arma::join_horiz(arma::join_horiz(2*u%dudx,(2*dudx%v+2*dvdx%u)),2*v%dvdx);
    arma::mat dQdxt = arma::trans(dQdx);
    arma::mat dDdx = dQdxt*Q + Qt*dQdx;
    arma::mat dFdx = -F*dDdx*F;
    arma::mat dbetadx = dFdx*(Qt*w) + F*(dQdxt*w) + F*(Qt*dwdx);
    double J_H_x = dbetadx(0)+dbetadx(2);

    arma::mat dQdy = arma::join_horiz(arma::join_horiz(2*u%dudy,(2*dudy%v+2*dvdy%u)),2*v%dvdy);
    arma::mat dQdyt = arma::trans(dQdy);
    arma::mat dDdy = dQdyt*Q + Qt*dQdy;
    arma::mat dFdy = -F*dDdy*F;
    arma::mat dbetady = dFdy*(Qt*w) + F*(dQdyt*w) + F*(Qt*dwdy);
    double J_H_y = dbetady(0)+dbetady(2);

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
int FSmeancurv::computeJslow(void)
{
  arma::mat dMdx,dMdy,dMdz;// sparse matrix slows down even more
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

  // xyz are row vectors in P
  // First is always self/center vertex
  // Next nnnbrs are the nearest neighbors
  P.zeros(npoints,3);
  P(0,0) = vxyz(vno,0);
  P(0,1) = vxyz(vno,1);
  P(0,2) = vxyz(vno,2);
  for(int n=1; n < npoints; n++){
    int nbrvno = vtop->v[n-1];
    P(n,0) = vxyz(nbrvno,0);
    P(n,1) = vxyz(nbrvno,1);
    P(n,2) = vxyz(nbrvno,2);
  }
  // Subtract the center vertex from all rows
  M = P - arma::repmat(P.row(0),npoints,1);

  // Rotate into the norm/tang space
  u = M*vc.e1;
  v = M*vc.e2;
  w = M*vc.norm;

  // Construct design matrix and solve. Keeping the intermediate steps
  // is helpful for computeing J
  Q = arma::join_horiz(arma::join_horiz(u%u,2*u%v),v%v);
  Qt = arma::trans(Q);
  D = Qt*Q;
  F = arma::inv(D);
  beta = F*Qt*w;
  H = (beta(0)+beta(2));
  vtx->curv = H;

  return(0);
}

int FSmeancurv::computeMaybeFaster(arma::mat &vxyz)
{
  // This does not appear to run faster (maybe slower)
  int vno = vc.vno;
  VERTEX_TOPOLOGY *vtop = vc.vtop;
  VERTEX *vtx = &(vc.surf->vertices[vno]);
  int npoints = vtop->vtotal+1;

  // xyz are row vectors in P
  P.zeros(npoints,3);
  P(0,0) = vxyz(vno,0);
  P(0,1) = vxyz(vno,1);
  P(0,2) = vxyz(vno,2);
  for(int n=1; n < npoints; n++){
    int nbrvno = vtop->v[n-1];
    P(n,0) = vxyz(nbrvno,0);
    P(n,1) = vxyz(nbrvno,1);
    P(n,2) = vxyz(nbrvno,2);
  }
  M = P - arma::repmat(P.row(0),npoints,1);

  // u  = M*e1, but one element of e1 is 0
  u = M.col(0)*vc.e1(0)+M.col(vc.cuse)*vc.e1(vc.cuse);
  v = M*vc.e2;
  w = M*vc.norm;

  arma::colvec u2 = u%u;
  arma::colvec uv = 2*u%v;
  arma::colvec v2 = v%v;
  Q = arma::join_horiz(arma::join_horiz(u2,uv),v2);
  Qt = arma::trans(Q);

  double d11 = arma::accu(u2%u2);
  double d12 = arma::accu(u2%uv);//=d21
  double d13 = arma::accu(u2%v2);//=d31
  double d22 = arma::accu(uv%uv);
  double d23 = arma::accu(uv%v2);//=d32
  double d33 = arma::accu(v2%v2);
  arma::mat DD = {{d11,d12,d13},{d12,d22,d23},{d13,d23,d33}};
  D = DD;
  //D = Qt*Q;
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
    // If vno0 != -1, then only apply to faces around vno0, otherwise
    // all faces
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
    // If vno0 != -1, then only apply to vno0, otherwise all
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
    // If vno0 != -1, then only apply to vno0, otherwise all
    if(mcurvs.size() != surf->nvertices){
      mcurvs.clear();
      mcurvs.resize(surf->nvertices);
    }
    if(vno0 > -1){
      mcurvs[vno0].vc = vertices[vno0]; // be good to have a pointer
      mcurvs[vno0].compute(vxyz);
      if(DoJ) mcurvs[vno0].computeJ();
    }
    else {
      #ifdef HAVE_OPENMP
       #pragma omp parallel for 
      #endif
      for(int vno=0; vno < surf->nvertices; vno++){
        mcurvs[vno].vc = vertices[vno]; // be good to have a pointer
	mcurvs[vno].compute(vxyz);
	if(DoJ) mcurvs[vno].computeJ();
      }
    }
    return(0);
  }
};
arma::colvec FSface::computeNorm(arma::mat &vxyz)
{
  // Load P. Note xyz is a column in P, not a row as it is in vxyz
  f = &surf->faces[fno]; 
  for(int r=0; r<3; r++){
    int vno = f->v[r];
    for(int k=0; k<3; k++){
      P(k,r) = vxyz(vno,k);
    }
  }
  v0 = P.col(2)-P.col(0);
  v1 = P.col(1)-P.col(2);
  m = arma::cross(v1,v0);
  L = arma::norm(m);
  norm = m/L;
  RecomputeJ=1;
  return(norm);
}
arma::mat computeJ_v_cross_xyz(arma::colvec const &v)
{
  // This is the Jacobian of the cross of v with dx,dy,dz
  // To get the cross of dx,dy,dz with v, just the negative
  arma::mat J_m_v = {{0,-v(2),v(1)},{v(2),0,-v(0)},{-v(1),v(0),0}};
  return(J_m_v);
}
int FSface::computeNormJslow(void)
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

int FSface::computeNormJ(void)
{
  // fast version
  if(!RecomputeJ) return(0);

  arma::mat J_m_v0 =  computeJ_v_cross_xyz(v1);
  arma::mat J_m_v1 = -computeJ_v_cross_xyz(v0);

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
  // Loop over neighboring triangles of the center vertex
  // and sum up, then average the normals to get t
  t.zeros(3);
  for(int k=0; k<vtop->vnum; k++){
    int fno = vtop->f[k];
    t += faces[fno].norm;
  }
  // Now normalize because the sum of normals is not unit 1
  L = arma::norm(t);
  norm = t/L;

  // To get e1, need a vector that is normal to the vertex
  // normal. This can be constructed directly from the normal.  Take
  // the largest of y and z to reduce the likelihood of numerical
  // problems if both (x,y) or (x,z) are both small. Note that
  // rotations of e1/e2 about the normal do not matter.
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
  // Now normalize
  L1 = arma::norm(v1);
  e1 = v1/L1;

  // e2 is the normal to both norm and e1
  v2 = arma::cross(norm,e1);
  L2 = arma::norm(v2);
  e2 = v2/L2;
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

  // Compute J_e?_n - grad of e? wrt the normal. Note: these do not
  // depend on the point, just the normal. 

  // Compute J_e1_n. This could be sped up a little because nA and nB
  // are sparse, but probably won't make much of a diff because
  // outside the point loop.
  arma::mat J_e1_n;
  arma::rowvec J_L1_n;
  if(cuse == 1){
    J_L1_n = (arma::trans(v1)*FS_J_v1_nA)/L1;
    J_e1_n = FS_J_v1_nA/L1 - (v1*J_L1_n)/(L1*L1);
  } else {
    J_L1_n = (arma::trans(v1)*FS_J_v1_nB)/L1;
    J_e1_n = FS_J_v1_nB/L1 - (v1*J_L1_n)/(L1*L1);
  }

  // Compute J_e2_n.  This could be sped up a little using
  // computeJ_v_cross_xyz, but probably does not seem to make much of
  // a diff.  Keeping the original slice it is simpler/clearer.
  arma::mat J_v2_nx = arma::cross(FSdx,e1) + arma::cross(norm,J_e1_n.col(0));
  arma::mat J_v2_ny = arma::cross(FSdy,e1) + arma::cross(norm,J_e1_n.col(1));
  arma::mat J_v2_nz = arma::cross(FSdz,e1) + arma::cross(norm,J_e1_n.col(2));
  arma::mat J_v2_n = arma::join_horiz(arma::join_horiz(J_v2_nx,J_v2_ny),J_v2_nz);
  //arma::colvec ax = arma::cross(norm,J_e1_n.col(0));
  //arma::colvec ay = arma::cross(norm,J_e1_n.col(1));
  //arma::colvec az = arma::cross(norm,J_e1_n.col(2));
  //arma::mat J_v2_n = -computeJ_v_cross_xyz(e1) + arma::join_horiz(arma::join_horiz(ax,ay),az);
  arma::rowvec J_L2_n = (arma::trans(v2)*J_v2_n)/L2;
  arma::mat J_e2_n;
  J_e2_n = J_v2_n/L2 - (v2*J_L2_n)/(L2*L2);

  // Now compute the J of norm, e1, and e2 wrt a given point.

  // Go through the first nnnbrs+1 points (and not all npoints)
  // because only changes at the first nnnbrs+1 points will affect
  // norm, e1, e2 at the center
  for(int n=0; n < nnnbrs+1; n++){
    int nvno;
    if(n == 0) nvno=vno; // first is always center vertex
    else       nvno = vtop->v[n-1]; // rest are nearest neighbors
    arma::mat J_t_p; 
    J_t_p.zeros(3,3);
    // Go through the all the neighboring faces to see which ones
    // contain both the center vertex and the current neighbor vertex
    int nhits = 0;
    for(int k=0; k<vtop->vnum; k++){
      int fno = vtop->f[k];
      FACE *vface = &surf->faces[fno];
      int m;
      for(m=0; m<3; m++) if(vface->v[m] == nvno) break;
      if(m==3) continue; // this nbr vertex not in this face
      nhits++;
      // This face is shared by both center and neighbor so add its J to the total
      J_t_p += faces[fno].J_n_p[m];
    }
    // Note: nhits should be 2 or vnum
    arma::mat J_L_p = (arma::trans(t)*J_t_p)/L;
    J_n_p[n]  = J_t_p/L - (t*J_L_p)/(L*L);
    // Apply the chain rule to get J_e?_p
    J_e1_p[n] = J_e1_n*J_n_p[n];
    J_e2_p[n] = J_e2_n*J_n_p[n];
  }

  return(0);
}

double FStestMeanCurvGrad(MRIS *surf, double Delta, int vno0=-1);
double FStestMeanCurvGrad(MRIS *surf, double Delta, int vno0)
{
  int vnostart=0;
  int vnostop=surf->nvertices;
  if(vno0 > -1){vnostart=vno0; vnostop=vno0+1;}
  int nverts = vnostop-vnostart;

  int threads = 1;
  #ifdef HAVE_OPENMP
  threads = omp_get_max_threads();
  #endif

  Timer mytimer;
  double t0,dt; 

  printf("FStestMeanCurvGrad() %d %d threads=%d Delta=%g\n",vnostart,vnostop,threads,Delta); 
  fflush(stdout);

  // The profiling might not be accurate. It seems like the 2nd time
  // it gets run, it runs much faster. Kind of weird.
  printf("Computing H and gradient ..."); fflush(stdout);
  FSsurf fs;
  t0 = mytimer.seconds();
  fs.surf = surf;
  fs.CopyVXYZ();
  fs.ComputeFaces(1,vno0);
  fs.ComputeVertices(1,vno0);
  fs.ComputeMeanCurvs(1,vno0);
  dt  = mytimer.seconds()-t0;
  printf("  dt = %g sec\n",dt);fflush(stdout);

  printf("Computing H  ..."); fflush(stdout);
  t0 = mytimer.seconds();
  FSsurf fs2;
  fs2.surf = surf;
  fs2.CopyVXYZ();
  // Might not need to do these, but it is fast
  fs2.ComputeFaces(0,vno0);
  fs2.ComputeVertices(0,vno0);
  fs2.ComputeMeanCurvs(0,vno0);
  dt  = mytimer.seconds()-t0;
  printf("  dt = %g sec\n",dt);fflush(stdout);

  printf("Computing numerical gradient\n");fflush(stdout);
  t0 = mytimer.seconds();
  double maxmax = 0;
  double JdiffRabsMax=0;
  int vnomax=0;
  for(int vno=vnostart; vno < vnostop; vno++){
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
    } // Loop over point set
    arma::mat Jdiff = J-fs.mcurvs[vno].J_H_p;
    arma::mat Jdiffabs = arma::abs(Jdiff);
    double Jdiffabsmax = Jdiffabs.max();
    if(maxmax < Jdiffabsmax) {
      maxmax = Jdiffabsmax;
      vnomax = vno;
      arma::mat JdiffR = Jdiff/((J+fs.mcurvs[vno].J_H_p)/2);
      arma::mat JdiffRabs = arma::abs(JdiffR);
      JdiffRabsMax = 100*JdiffRabs.max(); // make a percent
    }
    if(vno0 > -1){
      printf("vno = %d H=%g\n",vno,fs.mcurvs[vno].H);
      fs.mcurvs[vno].J_H_p.print("J_H_p = [");printf("];\n");
      J.print("Jnum = ["); printf("];\n");
      fs.mcurvs[vno].J_H_p.save("J_H_p.txt",arma::raw_ascii);
      J.save("Jnum.txt",arma::raw_ascii);
      Jdiff.print("Jdiff = ["); printf("];\n");
      printf("%d max Jdiff %g\n",vno,Jdiffabsmax);
      fs.mcurvs[vno].vc.norm.print("norm");
      fs.mcurvs[vno].vc.e1.print("e1");
      fs.mcurvs[vno].vc.e2.print("e2");
    }
  }
  dt  = mytimer.seconds()-t0;
  printf("maxmax= %6.4lf, rmax %g, vnomax = %d, dt=%g sec, nv=%d, dtp = %g sec\n",
	 maxmax,JdiffRabsMax,vnomax,dt,nverts,dt/nverts);
  fflush(stdout);
  return(maxmax);
}

//MAIN ----------------------------------------------------
int main(int argc, char **argv) 
{
  //Gdiag_no = 72533;
  MRIS *surf;
  //MRI *mri=NULL; // *mri2;
  int vno, threads=1;
  Timer mytimer;

  sscanf(argv[1],"%d\n",&threads);
  sscanf(argv[2],"%d\n",&vno);
  surf = MRISread(argv[3]);
  MRISsetNeighborhoodSizeAndDist(surf,2);

#ifdef HAVE_OPENMP
  omp_set_num_threads(threads);
#endif


  // Profiling. Runs faster the 2nd time around. Weird
  FSsurf fs;
  fs.surf = surf;
  for(int n=0; n < 2; n++){
    double t = mytimer.seconds();
    fs.CopyVXYZ();
    fs.ComputeFaces(0);
    fs.ComputeVertices(0);
    fs.ComputeMeanCurvs(0);
    printf("H %d %6.4f\n",n,mytimer.seconds()-t);
  }
  for(int n=0; n < 2; n++){
    double t = mytimer.seconds();
    fs.CopyVXYZ();
    fs.ComputeFaces(1);
    fs.ComputeVertices(1);
    fs.ComputeMeanCurvs(1);
    printf("dH %d %6.4f\n",n,mytimer.seconds()-t);
  }

  double Delta = 1e-10;
  FStestMeanCurvGrad(surf, Delta, vno);

  FILE *fp=NULL;
  VERTEX_TOPOLOGY *vtop;
  VERTEX *vtx;
  if(vno > -1){
    vtop = &(surf->vertices_topology[vno]);
    vtx = &(surf->vertices[vno]);
    printf("vno %d vnum %d vtotal %d\n",vno,vtop->vnum,vtop->vtotal);
    GetP(surf,vno);// saves P.txt
    fp = fopen("nnnbrs.dat","w");
    fprintf(fp,"%d\n",vtop->vnum);
    fclose(fp);
  }

  exit(0);

  //for(int n=0; n < 5; n++){ need #include to do this now
  //  double t0 = mytimer.seconds();
  //  mrisComputeTangentPlanes(surf);
  //  MRIScomputeSecondFundamentalFormThresholded(surf,-1);
  //  printf("dt  %6.4f\n",mytimer.seconds()-t0);
  //}


  //printf("FS: k1 %g k2 %g H %g\n",vtx->k1,vtx->k2,vtx->H);


  //mri = MRIcopyMRIS(NULL, surf, 0, "curv"); // start at z to autoalloc
  //MRIwrite(mri,argv[2]);
  //exit(0);


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

