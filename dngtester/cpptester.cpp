#include <sys/utsname.h>
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
#include "surfgrad.h"
#include "mrishash.h"
#include "mrisurf_metricProperties.h"
#include "mrishash_internals.h"
#include "cmdargs.h"
#include "version.h"
#include "label.h"

#define export
#define ENS_PRINT_INFO
#define ENS_PRINT_WARN
//#include <armadillo>
#include <ensmallen.hpp>
#include "romp_support.h"

//using namespace ens;
//using namespace arma;

int GetP(MRIS *surf, int vno)
{
  // This gets the xyz of vertices in the extended neighborhood (P)
  // around the the center vertex where the nearest neighbor vertices
  // are in a contiguous order. This allows for the contruction of the
  // triangles as the first points plus any two other contiguous
  // points. Eg, if there are 6 nearest neighbors then there are 6
  // nearest triangles. The xyz of the 1st triangle is P0+P1+P2. The
  // next one is P0+P2+P3, etc (P will be more than 6 because it
  // contains all the vertices from the extended neighborhood). Doing
  // it in this way allows for determining the xyz of the triangles
  // without having to have a separate structure from the xyz in P.
  // One of my matlab programs assumes this format.
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
#if 1
//======================================================================
int CBV(char *adgwsfile, MRI *mri_brain, MRIS *surf, int surftype, float max_cbv_dist=5.0);

arma::colvec FSdx = {1,0,0};
arma::colvec FSdy = {0,1,0};
arma::colvec FSdz = {0,0,1};
arma::mat FS_J_v1_nA = {{0,-1,0},{1,0,0},{0,0,0}};
arma::mat FS_J_v1_nB = {{0,0,-1},{0,0,0},{1,0,0}};
arma::mat computeJ_v_cross_xyz(arma::colvec const &v);
arma::mat computeJ_v_cross_xyz(arma::colvec const &v)
{
  // This is the Jacobian of the cross product of v with some other
  // vector (ie, cross(v,v2)) wrt to the xyz of the other vector.  To
  // get the Jacobian of cross of v2 with v (cross(v2,v)), just
  // multiply by -1. This is the same as cross(v,d{xyz})
  arma::mat J_m_v = {{0,-v(2),v(1)},{v(2),0,-v(0)},{-v(1),v(0),0}};
  return(J_m_v);
}

// New face structure to accomodate Armadillo and the computatio of
// the Jacobian of the normal wrt each point
class FSface{
public:
  MRIS *surf;
  int fno;
  FACE *f;
  arma::mat P = arma::mat(3,3); //col vects
  arma::colvec norm, v0, v1, m;
  double L;
  int RecomputeJ=1;
  // J_n_p[3] is the jacobian, each J_n_p[i] is a 3x3 matrix for
  // corner i. A row of J is the gradient wrt p.
  std::array<arma::mat,3> J_n_p;
  arma::colvec computeNorm(arma::mat &vxyz);
  int computeNormJ(void);
  int computeNormJslow(void);
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

  // J_n_p[i] is a 3x3 jacobian matrix. 

  return(0);
}

// This is a new vertex class to accomodate Armodillo and to assist in
// the computation of the Jacobian of the normal and tangent vectors
// wrt each point in the cell.
class FSvertexcell { 
public:
  MRIS *surf;
  int vno; // center vertex of the cell
  std::vector<FSface> faces; // faces in the cell
  VERTEX_TOPOLOGY *vtop=NULL;//topology of the center
  int nnnbrs=0; // number of nearest neighbors = vtop->vnum
  int npoints=0; // number of extended neighbors, including self = vtop->vtotal+1
  arma::colvec norm, e1, e2; // normal, 1st tangent vector, 2nd tangent vector
  arma::colvec t, v1, v2;
  double L, L1, L2;
  int cuse=0;
  // Jacobians. Each is a 3x3 matrix. Each row is the gradient of n or e1/2 wrt pxyz
  std::vector<arma::mat> J_n_p, J_e1_p, J_e2_p;
  int computeNorm(std::vector<FSface> const &faces);
  int computeNormJ(std::vector<FSface> const &faces);
};
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
  // This function computes the Jacobian of the vertex normal/e1/e2
  // vectors with respect to changes in the xyz of the center point or
  // nearest neighbors of the center point (points in the extended
  // neighborhood out side of the nearest do not affect the normal
  // calcuation). This is used in the computation of the mean curv
  // Jacobian to adjust/correct for changes in the norm/e1/e2 in
  // response to changes in these vertices above and beyond how such
  // changes would affect mean curv directly.
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

  // First, compute J_e?_n - grad of e? wrt the normal. Note: these do not
  // depend on the point directly, just the normal (as indicated by
  // the _n).

  // Compute J_e1_n. This could be sped up a little because nA and nB
  // are sparse, but probably won't make much of a diff because
  // outside the point loop.
  arma::mat J_e1_n;
  arma::rowvec J_L1_n;
  // cuse was set in FSvertexcell::computeNorm() and indicates whether e1
  // was computed from xy or xz
  if(cuse == 1){ // xy
    J_L1_n = (arma::trans(v1)*FS_J_v1_nA)/L1;
    J_e1_n = FS_J_v1_nA/L1 - (v1*J_L1_n)/(L1*L1);
  } else { // xz
    J_L1_n = (arma::trans(v1)*FS_J_v1_nB)/L1;
    J_e1_n = FS_J_v1_nB/L1 - (v1*J_L1_n)/(L1*L1);
  }

  // Compute J_e2_n.  This could be sped up a little using
  // computeJ_v_cross_xyz, but does not seem to make much of
  // a diff.  Keeping the original since it is simpler/clearer.
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


// Mean curvature structure for a vertex
class FSmeancurv{ 
public:
  FSvertexcell vc; // vertex info for this vertex
  double H; // mean curvature for this vertex
  // Jacobian wrt each point in the extended neighborhood
  arma::mat J_H_p; // npointsx3
  arma::mat P,M,Q,Qt,D,F,beta;
  arma::colvec  u, v, w;
  // For when optimizing: HO=target H, esp=H-H0, cost=eps^2
  double H0=0,eps=0,cost=0; 
  int compute(arma::mat &vxyz);
  int computeMaybeFaster(arma::mat &vxyz);
  int computeJ(void);
  int computeJslow(void);
}; 

int FSmeancurv::computeJ(void)
{
  // Compute Jacobian of H wrt each point in the extended neighborhood.
  // See the "slow" version for more matrix detail
  J_H_p.zeros(vc.npoints,3);
  for(int n=0; n < vc.npoints; n++){
    arma::colvec dudx, dudy, dudz,  dvdx, dvdy, dvdz,  dwdx, dwdy, dwdz;
    if(n==0){
      dudx.zeros(vc.npoints); dudx.fill(-vc.e1(0)); dudx(0) = 0;
      dudy.zeros(vc.npoints); dudy.fill(-vc.e1(1)); dudy(0) = 0;
      dudz.zeros(vc.npoints); dudz.fill(-vc.e1(2)); dudz(0) = 0;
      // Note: either e1(1) or e1(2) will be 0 so don't need to fill
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
      // This is the adjustment for changes in e1/e2/norm when the
      // center or nearest neighbor changes position.
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

  return(0); //  fast
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
      // This is the adjustment/correction for changes in e1/e2/norm
      // when the center or nearest neighbor changes position
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

  return(0); // slow
}
int FSmeancurv::compute(arma::mat &vxyz)
{
  int vno = vc.vno;
  VERTEX_TOPOLOGY *vtop = vc.vtop;
  VERTEX *vtx = &(vc.surf->vertices[vno]);
  int npoints = vtop->vtotal+1;

  // xyz are row vectors in P
  // First (P(0,:)) is always self/center vertex
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

  eps = H-H0; // can be done outside of class
  cost = eps*eps; // can be done outside of class

  return(0);
}

int FSmeancurv::computeMaybeFaster(arma::mat &vxyz)
{
  // This does not appear to run faster (maybe even slower)
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

  return(0);
}

// --------------------------------------------------------------
// Surface class to assist in the computing of mean curvature and
// its gradient.
class FSsurf
{
public:
  MRIS *surf=NULL; // Classic FS surface
  MRI *mask=NULL; // not used yet
  std::vector<FSface> faces;
  std::vector<FSvertexcell> vertices;
  std::vector<FSmeancurv> mcurvs;
  arma::mat vxyz; // xyz of vertices copied from surf
  int surftype = GRAY_WHITE; //GRAY_CSF
  int Label2Mask(LABEL *label);

  double long HCost(int vno0=-1);
  int HCostJ(int vno0=-1);
  arma::mat J_cH_p; // cost of H wrt p = nvertices x 3

  double long IntensityCostAndGrad(int DoJ, int vno0=-1);
  char *adgwsfile=NULL;
  MRI *mri_brain=NULL;
  double sigma_global=2.0;
  float max_cbv_dist = 5.0 ; // same as max_thickness in MMS
  arma::mat J_cI_p; // cost of intensity wrt p = nvertices x 3

  double long TangentialSpringCG(void);
  arma::mat J_cTS_p; // cost of intensity wrt p = nvertices x 3

  long double NTSpringCG(int SpringType);
  arma::mat J_cNS_p; // cost of intensity wrt p = nvertices x 3

  double SurfRepulsionCG(void);
  arma::mat J_cRep_p; 
  MHT *mht=NULL;

  double TargetSurfCG(void);
  arma::mat J_cTarg_p; 

  // Each element of VnoInNbrhdOf corresponds to a vertex. Each
  // element has list of center vertices in which the element vertex
  // is in the extended neighborhood of (includes self). The pair is
  // the vertex number of the center vertex and the nth point within
  // the extended neighborhood that the element vertex corresponds to.
  // Note that VnoInNbrhdOf does NOT indicate the extended neighbors
  // of the the element vertex; that info is already stored in vtop->v[]
  std::vector<std::vector<std::pair<int,int>>> VnoInNbrhdOf;
  int SetVnoInNbrhdOf(void){
    VnoInNbrhdOf.clear();
    VnoInNbrhdOf.resize(surf->nvertices);
    // First make vno in the nbr of itself
    for(int i=0; i < surf->nvertices; i++) {
      std::pair<int,int> p(i,0);
      VnoInNbrhdOf[i].push_back(p);
    }
    for(int i=0; i < surf->nvertices; i++){
      VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[i]);
      for(int n=0; n < vtop->vtotal; n++){
	int vno = vtop->v[n];
	// Indicate that vno is in the nbrhd of i at the (n+1)th position
	std::pair<int,int> p(i,n+1);
	VnoInNbrhdOf[vno].push_back(p);
      }
    }
    return(0);
  }
  // Copy the XYZ from the surface v->{x,y,z} into vxyz
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
  // Copy XYZ from vxyz to the surface v->{x,y,z} 
  int CopyVXYZtoSurf(const arma::mat &vxyz){
    for(int vno=0; vno < surf->nvertices; vno++){
      VERTEX *vtx = &surf->vertices[vno];
      vtx->x=vxyz(vno,0);
      vtx->y=vxyz(vno,1);
      vtx->z=vxyz(vno,2);
    }
    return(0);
  }
  int ComputeFaces(int DoJ, int vno0=-1){
    // DoJ=1 will cause the Jacobian to be computed
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
    // DoJ=1 will cause the Jacobian to be computed
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
    // DoJ=1 will cause the Jacobian to be computed
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
  void SetH0(double H0){// Set H0 (target H value) to a constant
    for(int vno=0; vno < surf->nvertices; vno++) mcurvs[vno].H0 = H0;
  }
  void SetH0ToCurv(void){// Set H0 (target H value) to v->curv
    for(int vno=0; vno < surf->nvertices; vno++) mcurvs[vno].H0 = surf->vertices[vno].curv;
  }
  void SetH0ToMRI(MRI *mri){// Set H0 (target H value) to values in MRI struct
    for(int vno=0; vno < surf->nvertices; vno++) mcurvs[vno].H0 = MRIgetVoxVal(mri,vno,0,0,0);
  }
};

int FSsurf::Label2Mask(LABEL *label)
{
  if(mask) MRIfree(&mask);
  mask = MRIalloc(surf->nvertices,1,1,MRI_UCHAR);
  // Create a mask=1 of vertices in the label.
  for(int n = 0; n < label->n_points; n++){
    int vno = label->lv[n].vno;
    MRIsetVoxVal(mask,vno,0,0,0,1);
  }
  // Extend the mask outside of the label to those vertices that are
  // in the extended nbhd of each vertex. Set mask=2 in these
  // vertices.
  for(int vno=0; vno < surf->nvertices; vno++){
    int m = MRIgetVoxVal(mask,vno,0,0,0);
    if(m != 1) continue;
    VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vno]);
    for(int nbr=0; nbr < vtop->vtotal; nbr++){
      int nbrvno = vtop->v[nbr];
      m = MRIgetVoxVal(mask,nbrvno,0,0,0);
      if(m) continue;
      MRIsetVoxVal(mask,nbrvno,0,0,0,2);
    }
  }
  // Extend the mask again to those vertices that are in the extended
  // nbhd of each marked vertex. Set mask=3 in these vertices.
  // This will assure that all the vertices/faces that could 
  // contribute to the cost. 
  for(int vno=0; vno < surf->nvertices; vno++){
    int m = MRIgetVoxVal(mask,vno,0,0,0);
    if(m != 1 && m != 2) continue;
    VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vno]);
    for(int nbr=0; nbr < vtop->vtotal; nbr++){
      int nbrvno = vtop->v[nbr];
      m = MRIgetVoxVal(mask,nbrvno,0,0,0);
      if(m) continue;
      MRIsetVoxVal(mask,nbrvno,0,0,0,3);
    }
  }

  return(0);
}

double FSsurf::TargetSurfCG(void)
{

  // Have to run this to get position
  CBV(adgwsfile, mri_brain, surf, surftype, max_cbv_dist);

  J_cTarg_p.zeros(surf->nvertices,3);
  double cost = 0;
  for(int vno = 0; vno < surf->nvertices; vno++) {
    VERTEX *v = &surf->vertices[vno];
    double x,y,z;
    x = vxyz(vno,0);
    y = vxyz(vno,1);
    z = vxyz(vno,2);
    double dx = x - v->targx;
    double dy = y - v->targy;
    double dz = z - v->targz;
    double d2 = dx*dx + dy*dy + dz*dz;
    double d = sqrt(d2);
    cost += d2;
    J_cTarg_p(vno,0) = dx/(d+FLT_EPSILON);
    J_cTarg_p(vno,1) = dy/(d+FLT_EPSILON);
    J_cTarg_p(vno,2) = dz/(d+FLT_EPSILON);
  }
  return(cost);
}

double FSsurf::SurfRepulsionCG(void)
{
  int vno, max_vno, i;
  float dx, dy, dz, x, y, z, sx, sy, sz, norm[3], dot=0;
  float max_scale, max_dot;
  double scale;

  J_cRep_p.zeros(surf->nvertices,3);
  double cost = 0;
  for (vno = 0; vno < surf->nvertices; vno++) {
    VERTEX *v = &surf->vertices[vno];

    x = vxyz(vno,0);
    y = vxyz(vno,1);
    z = vxyz(vno,2);

    MHBT *bucket = MHTacqBucket(mht, x, y, z);
    if (!bucket) continue;
    sx = sy = sz = 0.0;
    max_dot = max_scale = 0.0;
    max_vno = 0;

    MHB *bin = bucket->bins;
    for (i = 0; i < bucket->nused; i++, bin++) {
      VERTEX *vn = &(surf->vertices[bin->fno]);
      dx = x - vn->origx;
      dy = y - vn->origy;
      dz = z - vn->origz;
      mrisComputeOrigNormal(surf, bin->fno, norm);
      dot = dx * norm[0] + dy * norm[1] + dz * norm[2];
      if(dot > 1) continue;
      //scale = l_repulse * pow(1.0 - (double)dot, 4.0);
      scale = pow(1.0 - (double)dot, 4.0);
      sx += (scale * v->nx);
      sy += (scale * v->ny);
      sz += (scale * v->nz);
    }

    J_cRep_p(vno,0) = sx;
    J_cRep_p(vno,1) = sy;
    J_cRep_p(vno,2) = sz;

    MHTrelBucket(&bucket);
  }
  return(cost);
}


long double FSsurf::NTSpringCG(int SpringType)
{
  if(SpringType == 1) J_cNS_p.zeros(surf->nvertices,3);
  if(SpringType == 2) J_cTS_p.zeros(surf->nvertices,3);

  long double cost=0;
  for(int vno = 0; vno < surf->nvertices; vno++){
    VERTEX * const vertex = &surf->vertices[vno];
    VERTEX_TOPOLOGY const * const vertext = &surf->vertices_topology[vno];
    float nx = vertex->nx;
    float ny = vertex->ny;
    float nz = vertex->nz;

    float x = vxyz(vno,0);
    float y = vxyz(vno,1);
    float z = vxyz(vno,2);

    // compute the average distance between center and nbrs
    float sx = 0, sy = 0, sz = 0;
    for(int m = 0; m < vertext->vnum; m++) {
      int vnonbr = vertext->v[m];
      sx += vxyz(vnonbr,0)-x;
      sy += vxyz(vnonbr,1)-y;
      sz += vxyz(vnonbr,2)-z;
    }
    sx /= vertext->vnum;
    sy /= vertext->vnum;
    sz /= vertext->vnum;

    // project onto normal
    sx *= nx;
    sy *= ny;
    sz *= nz;

    if(SpringType == 1){
      J_cNS_p(vno,0) = sx;
      J_cNS_p(vno,1) = sy;
      J_cNS_p(vno,2) = sz;
    }
    if(SpringType == 2){
      // remove normal
      double nc = sx + sy + sz;
      sx -= (nc*nx);
      sy -= (nc*ny);
      sz -= (nc*nz);
      J_cTS_p(vno,0) = sx;
      J_cTS_p(vno,1) = sy;
      J_cTS_p(vno,2) = sz;
    }
    cost += (sx*sx + sy*sy + sz*sz);

  }
  return(cost);
}

double long  FSsurf::TangentialSpringCG(void)
{// replaced by NTSpringCG()
  J_cTS_p.zeros(surf->nvertices,3);

  long double cost=0;
  for(int vno = 0; vno < surf->nvertices; vno++){
    VERTEX * const vertex = &surf->vertices[vno];
    VERTEX_TOPOLOGY const * const vertext = &surf->vertices_topology[vno];
    
    float x = vxyz(vno,0);
    float y = vxyz(vno,1);
    float z = vxyz(vno,2);
    
    float sx = 0, sy = 0, sz = 0;
    for(int m = 0; m < vertext->vnum; m++) {
      int vnonbr = vertext->v[m];
      sx += vxyz(vnonbr,0)-x;
      sy += vxyz(vnonbr,1)-y;
      sz += vxyz(vnonbr,2)-z;
    }
    sx /= vertext->vnum;
    sy /= vertext->vnum;
    sz /= vertext->vnum;
    cost += (sx*sx + sy*sy + sz*sz);
    
    // Would normals need to be updated?
    float nx = vertex->nx;
    float ny = vertex->ny;
    float nz = vertex->nz;
    
    // project onto normal
    float nc = sx * nx + sy * ny + sz * nz;
    
    // remove normal component and scale
    J_cTS_p(vno,0) = (sx - nc * nx);
    J_cTS_p(vno,0) = (sy - nc * ny);
    J_cTS_p(vno,0) = (sz - nc * nz);
  }
  return(cost);
}


int CBV(char *adgwsfile, MRI *mri_brain, MRIS *surf, int surftype, float max_cbv_dist)
{
  double sigma_global=2.0;
  //float max_cbv_dist = 5.0 ; // same as max_thickness in MMS
  AutoDetGWStats adgws;
  int err = adgws.Read(adgwsfile);
  if(err){
    printf("ERROR: reading %s\n",adgwsfile);
    exit(1);
  }
  double inside_hi=0, border_hi=0, border_low=0, outside_low=0, outside_hi=0;
  if(surftype == GRAY_WHITE){
    inside_hi = adgws.white_inside_hi;
    border_hi = adgws.white_border_hi;
    if(adgws.white_border_low_factor > -9){
      double f = adgws.white_border_low_factor;
      border_low = f*adgws.gray_mean + (1-f)*adgws.white_mean;
    }
    outside_low = adgws.white_outside_low;
    outside_hi = adgws.white_outside_hi;
  }
  if(surftype == GRAY_CSF){
    inside_hi = adgws.pial_inside_hi;
    border_hi = adgws.pial_border_hi;
    border_low = adgws.pial_border_low;
    outside_low = adgws.pial_outside_low;
    outside_hi = adgws.pial_outside_hi;
  }
  printf("CBV ih=%g bh=%g bl=%g ol=%g oh=%g\n",inside_hi, border_hi, border_low, outside_low, outside_hi);
  printf("    sigma=%g  maxdist=%g  2maxdist=%g surftype=%d\n",sigma_global, max_cbv_dist, 2*max_cbv_dist, surftype);
  MRIScomputeBorderValues(surf, mri_brain, NULL, inside_hi, border_hi, border_low, outside_low, outside_hi,
			  sigma_global, 2*max_cbv_dist, NULL, surftype, NULL, 0.5, 0, NULL,-1,-1,0) ;
  return(0);
}

double long FSsurf::IntensityCostAndGrad(int DoJ, int vno0)
{
  int vno;
  VERTEX *v;
  float x, y, z, nx, ny, nz, dx, dy, dz;
  double xw, yw, zw, delI, delV;
  AutoDetGWStats adgws;
  int err = adgws.Read(adgwsfile);
  if(err) exit(1);
  double inside_hi=0, border_hi=0, border_low=0, outside_low=0, outside_hi=0;
  if(surftype == GRAY_WHITE){
    inside_hi = adgws.white_inside_hi;
    border_hi = adgws.white_border_hi;
    if(adgws.white_border_low_factor > -9){
      double f = adgws.white_border_low_factor;
      border_low = f*adgws.gray_mean + (1-f)*adgws.white_mean;
    }
    outside_low = adgws.white_outside_low;
    outside_hi = adgws.white_outside_hi;
  }
  if(surftype == GRAY_CSF){
    inside_hi = adgws.pial_inside_hi;
    border_hi = adgws.pial_border_hi;
    border_low = adgws.pial_border_low;
    outside_low = adgws.pial_outside_low;
    outside_hi = adgws.pial_outside_hi;
  }
  MRI *seg=NULL;

  if(DoJ){
    // This sets target intensity value v->val
    MRIScomputeBorderValues(surf, mri_brain, NULL, inside_hi,border_hi,border_low,outside_low,outside_hi,
			    sigma_global, 2*max_cbv_dist, NULL, surftype, NULL, 0.5, 0, seg,-1,-1,0) ;
    int vavgs=2; // spatial smoothing of the target value
    if(vavgs > 0) MRISaverageMarkedVals(surf, vavgs) ;
    J_cI_p.zeros(surf->nvertices,3);
  }

  long double cost = 0;
  for (vno = 0; vno < surf->nvertices; vno++) {
    v = &surf->vertices[vno];
    if (v->ripflag || v->val < 0) continue;
    x = vxyz(vno,0);
    y = vxyz(vno,1);
    z = vxyz(vno,2);

    // Sample the actual intensity at this vertex
    double valActual=0;
    MRISsurfaceRASToVoxelCached(surf, mri_brain, x, y, z, &xw, &yw, &zw);
    MRIsampleVolume(mri_brain, xw, yw, zw, &valActual);
    //if(vno == Gdiag_no) printf("vno =%d xyz = %g %g %g   %g %g %g  valActual=%g valTarg=%g\n",vno,x,y,z,xw,yw,zw,valActual,v->val);

    // epsilon = Difference between target intensity and actual intensity
    delV = valActual - v->val; // was v->val - valActual; 
    // Dont allow the difference to be greater than 5 or less than -5
    // Hidden parameter 5
    //if (delV > 5)        delV = 5;
    //else if (delV < -5)  delV = -5;
    cost += (delV*delV);

    if(!DoJ) continue;

    FSvertexcell vc = vertices[vno];
    nx = vc.norm(0);
    ny = vc.norm(1);
    nz = vc.norm(2);

    /* Numerically compute intensity gradient along the normal. Only used to get the right sign */
#if 0
    double val_outside, val_inside , k, ktotal_outside, ktotal_inside;
    double sigma = v->val2; // smoothing level for this vertex 
    if (FZERO(sigma)) sigma = sigma_global;
    if (FZERO(sigma)) sigma = 0.25;
    // Hidden parameter used to compute the step size
    double step_size;
    step_size = MIN(sigma / 2, MIN(mri_brain->xsize, MIN(mri_brain->ysize, mri_brain->zsize)) * 0.5);
    double dist, val;
    int n;
    ktotal_inside = ktotal_outside = 0.0;
    n = 0, val_outside = val_inside = 0.0;
    for(dist = step_size; dist <= 2 * sigma; dist += step_size, n++) {
      k = exp(-dist * dist / (2 * sigma * sigma));
      xw = x + dist * nx;
      yw = y + dist * ny;
      zw = z + dist * nz;
      ktotal_outside += k;
      MRISsurfaceRASToVoxelCached(surf, mri_brain, xw, yw, zw, &xw, &yw, &zw);
      MRIsampleVolume(mri_brain, xw, yw, zw, &val);
      val_outside += k * val;

      xw = x - dist * nx;
      yw = y - dist * ny;
      zw = z - dist * nz;
      MRISsurfaceRASToVoxelCached(surf, mri_brain, xw, yw, zw, &xw, &yw, &zw);
      MRIsampleVolume(mri_brain, xw, yw, zw, &val);
      val_inside += k * val;
      ktotal_inside += k;
    }
    if(ktotal_inside  > 0) val_inside  /= (double)ktotal_inside;
    if(ktotal_outside > 0) val_outside /= (double)ktotal_outside;

    // Gradient of the intensity at this location wrt a change along the normal
    delI = (val_outside - val_inside) / 2.0;
    // Change delI into +1 or -1
    //if (!FZERO(delI))  delI /= fabs(delI);
    //else               delI = -1; /* intensities tend to increase inwards */
#endif

    double Delta = 1.0, valin, valout, c, r, s;
    xw = x - nx * Delta/2;
    yw = y - ny * Delta/2;
    zw = z - nz * Delta/2;
    MRISsurfaceRASToVoxelCached(surf, mri_brain, xw, yw, zw, &c,&r,&s);
    MRIsampleVolume(mri_brain, c, r, s, &valin);
    xw = x + nx * Delta/2;
    yw = y + ny * Delta/2;
    zw = z + nz * Delta/2;
    MRISsurfaceRASToVoxelCached(surf, mri_brain, xw, yw, zw, &c,&r,&s);
    MRIsampleVolume(mri_brain, c, r, s, &valout);
    delI = (valout-valin)/Delta;

    // Weight intensity error by cost weighting
    double del = delV * delI;

    // Set to push vertex in the normal direction by this amount
    dx = nx * del;
    dy = ny * del;
    dz = nz * del;

    if(vno == Gdiag_no){
      printf("vno=%d  valActual=%g val=%g delV=%g delI=%g valin=%g valout=%g (%g,%g,%g) n=(%g,%g,%g)\n",
        vno,valActual,v->val,delV,delI,valin,valout,dx,dy,dz,nx,ny,nz);
      printf("   xyz = %g %g %g   crs = %g %g %g  valActual=%g valTarg=%g delV=%g\n",xw,yw,zw,c,r,s,valActual,v->val,delV);
      fflush(stdout);
    }

    J_cI_p(vno,0) = dx;
    J_cI_p(vno,1) = dy;
    J_cI_p(vno,2) = dz;
  } // vertices
  return(cost);
}

double long FSsurf::HCost(int vno0)
{
  // could go over xnbhrd of vno0, but this should be pretty fast
  double long cost = 0;
  #ifdef HAVE_OPENMP
  #pragma omp parallel for reduction(+:cost) 
  #endif
  for(int vno=0; vno < surf->nvertices; vno++) {
    mcurvs[vno].eps = mcurvs[vno].H-mcurvs[vno].H0; // epsilon=H-H0
    mcurvs[vno].cost = pow(mcurvs[vno].eps,2.0); // cost = epsilon^2
    cost += mcurvs[vno].cost;
  }
  return(cost);
}

int FSsurf::HCostJ(int vno0)
{
  // Compute J_cH_p, the Jacobian of the cost of H wrt the xyz of the vertices
  // J_cH_p is nvertices x 3
  int vnostart=0;
  int vnostop=surf->nvertices;
  if(vno0 > -1){vnostart=vno0; vnostop=vno0+1;}

  J_cH_p.zeros(surf->nvertices,3);
  for(int vno=vnostart; vno < vnostop; vno++) {
    for(int k=0; k < VnoInNbrhdOf[vno].size(); k++){
      std::pair<int,int> p = VnoInNbrhdOf[vno][k];
      int j = p.first;
      int n = p.second;
      J_cH_p.row(vno) += (mcurvs[j].eps*mcurvs[j].J_H_p.row(n));
    }
  }
  J_cH_p *= 2;
  return(0);
}
#endif

// Class to interface with ensmallen to optimize a cost function
// Before running the class, make sure to set
//  myopt.fs.surf = surf;
//  myopt.fs.SetVnoInNbrhdOf();
//  myopt.fs.CopyVXYZ();
// Do not pass fs.vxyz as the argument to the Optimizer, eg,
//    double cost = optimizer.Optimize(moyopt, vxyz);
//  NOT -->  double cost = optimizer.Optimize(moyopt, myopt.fs.vxyz); NOT!!
class MyOpt
{
public:
  int ncalls=0;
  int nevcalls=0;
  int GradOnly=0;
  FSsurf fs;
  double wI = 1;
  double wS = 3;
  double wRep = 0;
  double wH = 10000;
  double wTarg = 0;
  // It appears that it will use EvalWithGrad() only if it is there,
  // which is strange because I would have through that the line
  // search would have used Eval. If EvalWithGrad() is removed
  // leaving only Eval and Grad, then it will call them separately
  // and a differnt number of times, but it ends up being slower. 
  double Evaluate(const arma::mat& vxyz){
    Timer mytimer; double t0 = mytimer.seconds();
    double c = 0, cH=0, cI=0, cNS=0, cTS=0, cRep=0,cTarg=0;
    fs.vxyz = vxyz;
    fs.CopyVXYZtoSurf(vxyz);
    fs.ComputeFaces(0);
    fs.ComputeVertices(0);
    if(wH>0){
      fs.ComputeMeanCurvs(0);
      cH = fs.HCost();
      c += wH*cH;
    }
    if(wI > 0){
      cI = fs.IntensityCostAndGrad(0);
      c += wI*cI;
    }
    if(wRep > 0){
      cRep = fs.SurfRepulsionCG();
      c += wRep*cRep;
    }
    if(wS > 0){
      cNS = fs.NTSpringCG(1);
      cTS = fs.NTSpringCG(2);
      c += (wS*(cNS+cTS));
    }
    if(wTarg > 0){
      cTarg = wTarg*fs.TargetSurfCG();
      c += cTarg;
    }
    nevcalls ++;
    printf("  Eval %4d tot=%7.5lf  I=%7.5lf NS=%7.5lf TS=%7.5lf H=%7.5lf Targ=%7.5lf  dt=%4.3lf\n",
	   nevcalls,c,cI,cNS,cTS,cH,cTarg,mytimer.seconds()-t0);
    fflush(stdout);
    return(c);
  }
  double EvaluateWithGradient(const arma::mat& vxyz, arma::mat& g){
    Timer mytimer; double t0 = mytimer.seconds();
    double c = 0, cH=0, cI=0, cNS=0, cTS=0, cRep=0, cTarg=0;
    fs.vxyz = vxyz;
    fs.CopyVXYZtoSurf(vxyz);
    fs.ComputeFaces(1);
    fs.ComputeVertices(1);
    if(wH>0){
      fs.ComputeMeanCurvs(1);
      cH = fs.HCost();
      fs.HCostJ();
      c += wH*cH;
      g += wH*fs.J_cH_p;
    }
    if(wI > 0){
      cI = fs.IntensityCostAndGrad(1);
      c += wI*cI;
      g += wI*fs.J_cI_p;
    }
    if(wRep > 0){
      cRep = fs.SurfRepulsionCG();
      c += wRep*cRep;
      g += wRep*fs.J_cRep_p;
    }
    if(wS > 0){
      cNS = fs.NTSpringCG(1);
      cTS = fs.NTSpringCG(2);
      c += (wS*(cNS+cTS));
      g += (wS*(fs.J_cNS_p+fs.J_cTS_p));
    }
    if(wTarg > 0){
      cTarg = wTarg*fs.TargetSurfCG();
      c += cTarg;
      g += (wTarg*fs.J_cTarg_p);
    }
    ncalls ++;
    printf("  EWG%d %4d tot=%7.5lf  I=%7.5lf NS=%7.5lf TS=%7.5lf H=%7.5lf cTarg=%7.5lf  dt=%4.3lf\n",
	   GradOnly,ncalls,c,cI,cNS,cTS,cH,cTarg,mytimer.seconds()-t0);
    if(Gdiag_no > 0){
      VERTEX *vtx = &(fs.surf->vertices[Gdiag_no]);
      printf("  EWGvno %d  (%g,%g,%g) (%g,%g,%g)  (%g,%g,%g)\n",Gdiag_no,vtx->x,vtx->y,vtx->z,
	     vxyz(Gdiag_no,0),vxyz(Gdiag_no,1),vxyz(Gdiag_no,2),
	     g(Gdiag_no,0),g(Gdiag_no,1),g(Gdiag_no,2));
    }
    fflush(stdout);
    return(c);
  }
  void Gradient(const arma::mat& vxyz, arma::mat& g){
    GradOnly=1;
    EvaluateWithGradient(vxyz, g);
    GradOnly=0;
  }
  double LinMin(arma::mat& vxyz, double stepsize, int nsteps, int niters){
    arma::mat vxyz0 = vxyz, vxyzmin;
    arma::mat g;

    int vno = 175;
    double cmin = 10e10;
    for(int k=0; k < 1000; k++){
      g.zeros(vxyz.n_rows,vxyz.n_cols);
      double c0 = EvaluateWithGradient(vxyz, g);
      if(k==0) cmin = c0;
      vxyzmin = vxyz;
      double cprev = c0;
      printf("#$$ k = %d =============================================\n",k);
      printf("LM 00 %6.4lf   %6.4lf %6.4lf %6.4lf \n",cmin,vxyz(vno,0),vxyz(vno,1),vxyz(vno,2));
      printf(" Jvno %6.4lf %6.4lf %6.4lf eps=%g\n",g(vno,0),g(vno,1),g(vno,2),fs.mcurvs[vno].eps);
      printf(" vno H0 %6.4lf H %6.4lf \n",fs.mcurvs[vno].H0,fs.mcurvs[vno].H);
      
      int minhit = 0;
      for(int n = 0; n < nsteps; n++){
	double d = n; //-nsteps/2.0 + n;
	arma::mat vxyzstep = vxyz - g*stepsize*d; // should check neg as well
	double cstep = Evaluate(vxyzstep);
	printf("LM %d %12.9lf   %6.4lf %6.4lf %6.4lf ",n,cstep,vxyz(vno,0),vxyz(vno,1),vxyz(vno,2));
	if(cmin > cstep) {
	  cmin = cstep;
	  vxyzmin = vxyzstep;
	  minhit = 1;
	  printf("min\n");
	}
	else if((cstep > cprev)){ // && minhit
	  printf("increase\n");
	  //printf(" Cost increaed, breaking\n");
	  break;
	}
	else printf("\n");
	fflush(stdout);
	cprev = cstep;
      } // step
      vxyz = vxyzmin;
    }
    return(cmin);
  }

};

// Numeric test of Jacobian of mean curv wrt P
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
  fs.SetVnoInNbrhdOf();
  fs.CopyVXYZ();
  fs.ComputeFaces(1,vno0);
  fs.ComputeVertices(1,vno0);
  fs.ComputeMeanCurvs(1,vno0);
  dt  = mytimer.seconds()-t0;
  printf("  dt = %g sec\n",dt);fflush(stdout);
  t0 = mytimer.seconds();
  printf("H and J  dt = %g sec\n",dt);fflush(stdout);

  printf("Computing H  ..."); fflush(stdout);
  t0 = mytimer.seconds();
  FSsurf fs2;
  fs2.surf = surf;
  fs2.SetVnoInNbrhdOf();
  fs2.CopyVXYZ();
  // Might not need to do these, but it is fast
  fs2.ComputeFaces(0,vno0);
  fs2.ComputeVertices(0,vno0);
  fs2.ComputeMeanCurvs(0,vno0);
  dt  = mytimer.seconds()-t0;
  printf("  dt = %g sec\n",dt);fflush(stdout);

  printf("Computing numerical gradient\n");fflush(stdout);
  double maxmax = 0;
  double JdiffRabsMax=0;
  int vnomax=0;
  t0 = mytimer.seconds();
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
	J(n,c) = (fs2.mcurvs[vno].H-fs.mcurvs[vno].H)/Delta;
	fs2.vxyz(vnbrno,c) = psave;
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
  printf("maxmax= %6.4lf, rmax %g, vnomax = %d, dt=%g sec, nv=%d, dtp = %g sec  %d\n",
	 maxmax,JdiffRabsMax,vnomax,dt,nverts,dt/nverts,GetVmPeak());
  fflush(stdout);

  return(maxmax);
}

// Numeric test of Jacobian of mean curv cost wrt vertex xyz
double FStestMeanCurvCostGrad(MRIS *surf, double Delta, int vno0=-1);
double FStestMeanCurvCostGrad(MRIS *surf, double Delta, int vno0)
{
  int vnostart=0;
  int vnostop=surf->nvertices;
  if(vno0 > -1){vnostart=vno0; vnostop=vno0+1;}

  int threads = 1;
  #ifdef HAVE_OPENMP
  threads = omp_get_max_threads();
  #endif

  Timer mytimer;
  double t0,dt; 

  printf("FStestMeanCurvCostGrad() %d %d threads=%d Delta=%g\n",vnostart,vnostop,threads,Delta); 
  fflush(stdout);

  // The profiling might not be accurate. It seems like the 2nd time
  // it gets run, it runs much faster. Kind of weird.
  printf("Computing H and gradient ..."); fflush(stdout);
  FSsurf fs;
  fs.surf = surf;
  fs.SetVnoInNbrhdOf();
  fs.CopyVXYZ();
  fs.ComputeFaces(1);
  fs.ComputeVertices(1);
  fs.ComputeMeanCurvs(1);
  double long hc = fs.HCost();
  fs.HCostJ(vno0);

  printf("H cost = %g\n",(double)hc);
  fs.J_cH_p.save("J_cH_p.dat",arma::raw_ascii);

  printf("Computing H  ..."); fflush(stdout);
  FSsurf fs2;
  fs2.surf = surf;
  fs2.SetVnoInNbrhdOf();
  fs2.CopyVXYZ();
  // Might not need to do these, but it is fast
  fs2.ComputeFaces(0);
  fs2.ComputeVertices(0);
  fs2.ComputeMeanCurvs(0);
  fs2.HCost();

  printf("Computing numerical gradient\n");fflush(stdout);
  arma::mat Jc(surf->nvertices,3);
  Jc.zeros();
  t0 = mytimer.seconds();
  for(int vno=vnostart; vno < vnostop; vno++){
    //VERTEX_TOPOLOGY *vtop = &(surf->vertices_topology[vno]);
    for(int c=0; c < 3; c++){
      double psave = fs2.vxyz(vno,c);
      fs2.vxyz(vno,c) += Delta;
      // to use vno in compute, must be able to go over all nbrhds
      fs2.ComputeFaces(0);
      fs2.ComputeVertices(0);
      fs2.ComputeMeanCurvs(0);
      double long hc2 = fs2.HCost();
      double long dd = (hc2-hc)/Delta;
      Jc(vno,c) = dd;
      if(0){
      for(int k=0; k < fs.VnoInNbrhdOf[vno].size(); k++){
	std::pair<int,int> p = fs.VnoInNbrhdOf[vno][k];
	int j = p.first;
	//int n = p.second;
	double dd = (fs2.mcurvs[j].cost-fs.mcurvs[j].cost)/Delta;
	// Use vno because this is the change in the total cost caused by
	// vno regardless of what vertex j it shows up in
	Jc(vno,c) += dd;
      }
      }
      fs2.vxyz(vno,c) = psave;
    } // Loop over coords
  } // Loop over vertices

  dt  = mytimer.seconds()-t0;
  printf("dt = %g\n",dt);
  fflush(stdout);
  Jc.save("J_cH_p.num.txt",arma::raw_ascii);
  arma::mat Jdiff = Jc-fs.J_cH_p;
  arma::mat Jdiffabs = arma::abs(Jdiff);
  double Jdiffabsmax = Jdiffabs.max();
  printf("HCostJ diff %g\n",Jdiffabsmax);

  return(0);
}

int SurfDiffStats(MRIS *surf1, MRIS *surf2)
{
  double rmsmax=0, rmssum=0, rmssum2=0;
  int vnomax = -1;
  for(int vno=0; vno < surf1->nvertices; vno++){
    VERTEX *vtx1 = &surf1->vertices[vno];
    VERTEX *vtx2 = &surf2->vertices[vno];
    double dx = vtx2->x - vtx1->x;
    double dy = vtx2->y - vtx1->y;
    double dz = vtx2->z - vtx1->z;
    double rms = sqrt(dx*dx + dy*dy + dz*dz);
    if(rmsmax < rms) {
      rmsmax = rms;
      vnomax = vno;
    }
    rmssum += rms;
    rmssum2 += (rms*rms);
  }
  printf("RMS: max=%g imax=%d mean=%g  %g\n",rmsmax,vnomax,rmssum/surf1->nvertices, rmssum2);
  return(0);
}

MRI *AramaMat2MRI(arma::mat mat)
{
  int nframes = mat.n_cols;
  int nvox = mat.n_rows;
  MRI *mri = MRIallocSequence(nvox,1,1,MRI_FLOAT,nframes);

  //printf("nvox = %d nframes = %d\n",nvox,nframes);
  for(int r=0; r<nvox; r++){
    for(int c=0; c<nframes; c++){
      //printf("r=%d c=%d\n",r,c); fflush(stdout);
      MRIsetVoxVal(mri,r,0,0,c,mat(r,c));
      fflush(stderr);
    }
  }
  return(mri);
}



int  parse_commandline(int argc, char **argv);
void check_options();
void print_usage();
void usage_exit();
void print_help();
void dump_options(FILE *fp);

char *surffile = NULL;
char *adgwsfile = NULL;
char *imrifile = NULL;
char *outsurffile = NULL;
double wI = 1;
double wS = 30;
double wRep = 0;
double wH = 0;
double wTarg = 0;
int threads=1;
int surftype=1;
int iters = 1;
int saveiters = 0;
int debug = 0, checkoptsonly = 0;
char *cmdline2, cwd[2000];
char *outcurvfile=NULL;
char *curvfile=NULL;
char *initcurvfile=NULL;
LABEL *label=NULL;
char *targetsurffile=NULL;
double cbvdist=5.0;
double NoiseLevel=0;
int DoBFGS = 1;

//MAIN ----------------------------------------------------
int main(int argc, char **argv) 
{
  int nargs = handleVersionOption(argc, argv, "cpptester");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  char *cmdline = argv2cmdline(argc,argv);
  //uname(&uts);
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

  printf("\ncd %s\n",cwd);
  printf("%s\n\n",cmdline);
  dump_options(stdout);

#ifdef HAVE_OPENMP
  omp_set_num_threads(threads);
#endif

  check_options();
  //if(checkoptsonly) return(0);

  Timer mytimer;
  //double tstart = mytimer.seconds();

  MyOpt mo;

  // Read in the surface and set up topo
  mo.fs.surf = MRISread(surffile); // surface to optimize
  if(!mo.fs.surf) exit(1);
  MRIS *surf0 = MRISread(surffile);
  MRISsaveVertexPositions(mo.fs.surf, ORIGINAL_VERTICES) ; // This is used for repulsion and CBV
  MRISsetNeighborhoodSizeAndDist(mo.fs.surf,2);

  if(label) mo.fs.Label2Mask(label);

  mo.fs.adgwsfile = adgwsfile;
  if(imrifile)  mo.fs.mri_brain = MRIread(imrifile);
  mo.fs.surftype = surftype;
  mo.wI = wI;
  mo.wS = wS;
  mo.wH = wH;
  mo.wTarg = wTarg;
  mo.fs.max_cbv_dist = cbvdist;

  // Set up the optimization structure
  mo.fs.SetVnoInNbrhdOf();
  mo.fs.CopyVXYZ();
  mo.fs.ComputeFaces(0);
  mo.fs.ComputeVertices(0);
  mo.fs.ComputeMeanCurvs(0);

  if(curvfile){
    printf("Setting target curv from %s\n",curvfile);
    MRI *curv0 = MRIread(curvfile);
    if(curv0==NULL) exit(1);
    mo.fs.SetH0ToMRI(curv0);
  } 
  else  mo.fs.SetH0ToCurv();

  if(initcurvfile){
    // Write out the initial H
    MRI *initcurv = MRIcopyMRIS(NULL, mo.fs.surf, 0, "curv");
    int err = MRIwrite(initcurv,initcurvfile);
    if(err) exit(1);
    MRIfree(&initcurv);
  }

  if(NoiseLevel > 0){
    printf("Adding noise to the surface %g\n",NoiseLevel);
    for(int vno=0; vno < mo.fs.surf->nvertices; vno++){
      VERTEX *vtx = &(mo.fs.surf->vertices[vno]);
      vtx->x += (2*NoiseLevel*vtx->nx)*(drand48()-0.5);
      vtx->y += (2*NoiseLevel*vtx->ny)*(drand48()-0.5);
      vtx->z += (2*NoiseLevel*vtx->nz)*(drand48()-0.5);
    }
    MRISwrite( mo.fs.surf,"noisy.srf");
    // Recompute
    mo.fs.CopyVXYZ();
    mo.fs.ComputeFaces(0);
    mo.fs.ComputeVertices(0);
    mo.fs.ComputeMeanCurvs(0);
    printf("Difference between true input and noisy surfs\n");
    MRISdiffSimple(mo.fs.surf, surf0, 0, .00000001, 0);
  }

  // Get the initial cost
  double c0 = mo.Evaluate(mo.fs.vxyz);
  printf("Init cost = %g\n",c0);fflush(stdout);

  printf("Starting optimization\n"); fflush(stdout);
  double t0 = mytimer.seconds();


  if(DoBFGS){
    double armijoConst = 1e-4; // default is 1e-4
    double wolfe = .9; // default is 0.9
    ens::L_BFGS optimizer(10, 100, armijoConst, wolfe, 1e-6, 1e-15, 50, 1e-20);
    arma::mat vxyz = mo.fs.vxyz;
    for(int k=0; k < iters; k++){
      double cost = optimizer.Optimize(mo, vxyz); //ens::PrintLoss(stdout)
      printf("k= %2d cost = %g   %g  %d %d\n",k,cost,mytimer.seconds()-t0,mo.ncalls,mo.nevcalls);fflush(stdout);
      //printf("#VMPC#k%d  VmPeak  %d\n",k,GetVmPeak());
      mo.fs.CopyVXYZtoSurf(vxyz);
      if(saveiters){
	char tmpstr[1000];
	sprintf(tmpstr,"%s.k%02d",outsurffile,k);
	MRISwrite(mo.fs.surf,tmpstr);
      }
    }
  }
  else  mo.LinMin(mo.fs.vxyz,.1,20,iters);

  if(label){
    for(int vno=0; vno < mo.fs.surf->nvertices; vno++){
      int m = MRIgetVoxVal(mo.fs.mask,vno,0,0,0);
      if(m) continue;
      VERTEX *vtx = &mo.fs.surf->vertices[vno];
      VERTEX *vtx0 = &surf0->vertices[vno];
      vtx->x = vtx0->x;
      vtx->y = vtx0->y;
      vtx->z = vtx0->z;
    }
  }
  if(!saveiters) MRISwrite(mo.fs.surf,outsurffile);

  mo.fs.CopyVXYZ();// from surface
  mo.fs.ComputeFaces(0);
  mo.fs.ComputeVertices(0);
  mo.fs.ComputeMeanCurvs(0);
  if(outcurvfile){
    // Write out the final H
    MRI *finalcurv = MRIcopyMRIS(NULL, mo.fs.surf, 0, "curv");
    int err = MRIwrite(finalcurv,outcurvfile);
    if(err) exit(1);
    MRIfree(&finalcurv);
  }
  if(targetsurffile){
    MRISsaveVertexPositions(mo.fs.surf, TMP2_VERTICES) ; 
    for(int v=0; v < mo.fs.surf->nvertices;v++){
      mo.fs.surf->vertices[v].x = mo.fs.surf->vertices[v].targx;	
      mo.fs.surf->vertices[v].y = mo.fs.surf->vertices[v].targy;	
      mo.fs.surf->vertices[v].z = mo.fs.surf->vertices[v].targz;	
    }
    int err = MRISwrite(mo.fs.surf,targetsurffile);
    if(err) exit(1);
    MRISrestoreVertexPositions(mo.fs.surf, TMP2_VERTICES);
  }

  if(0 && wH>0){
    printf("Saving J\n");
    MRI *J = AramaMat2MRI(mo.fs.J_cH_p);
    MRIwrite(J,"j.mgz");
    MRIfree(&J);
  }

  printf("Difference between true input and output surfs\n");
  MRISdiffSimple(mo.fs.surf, surf0, 0, .00000001, 0);

  printf("Init cost = %g\n",c0);fflush(stdout);
  printf("cpptester-runtime   %g %d %d\n",mytimer.seconds()-t0,mo.ncalls,mo.nevcalls);fflush(stdout);
  printf("#VMPC# cpptester VmPeak  %d\n",GetVmPeak());

  exit(0);

  #if 0
  //double Delta = 1e-10;
  //FStestMeanCurvCostGrad(surf, Delta, vno);
  //printf("==============================================\n");
  //FStestMeanCurvGrad(surf, Delta, vno);
  //printf("==============================================\n");
  //SurfDiffStats(surf,surf0);
  #endif


}

int parse_commandline(int argc, char **argv) {
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
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--white")) surftype = GRAY_WHITE; 
    else if (!strcasecmp(option, "--pial"))  surftype = GRAY_CSF; 
    else if (!strcasecmp(option, "--save-iters"))  saveiters = 1;

    else if(!strcasecmp(option, "--surf")){
      if (nargc < 1) CMDargNErr(option,1);
      surffile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--curv")){
      if (nargc < 1) CMDargNErr(option,1);
      curvfile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--out-curv")){
      if (nargc < 1) CMDargNErr(option,1);
      outcurvfile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--init-curv")){
      // for writing out the starting curv file
      if (nargc < 1) CMDargNErr(option,1);
      initcurvfile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--o")){
      if (nargc < 1) CMDargNErr(option,1);
      outsurffile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--target")){
      if (nargc < 1) CMDargNErr(option,1);
      targetsurffile = pargv[0];
      nargsused = 1;
    }
    else if(!strcasecmp(option, "--vol")){
      if(nargc < 2) CMDargNErr(option,1);
      imrifile = pargv[0];
      adgwsfile = pargv[1];
      nargsused = 2;
    }
    else if (!strcasecmp(option, "--wI")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&wI);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--wS")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&wS);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--wH")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&wH);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--wT")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&wTarg);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--cbv-dist")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&cbvdist);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--noise")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&NoiseLevel);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--threads")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&threads);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--iters")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&iters);
      nargsused = 1;
    }
    else if(!strcmp(option, "--debug-vertex")){
      if(nargc < 1) CMDargNErr(option,1);
      Gdiag_no = atoi(pargv[0]) ;
      printf("Gdiag_no set to %d\n",Gdiag_no);
      nargsused = 1;
    }
    else if(!strcmp(option, "--label")){
      if(nargc < 1) CMDargNErr(option,1);
      label = LabelRead(NULL,pargv[0]) ;
      if(!label) exit(1);
      nargsused = 1;
    }
    else if(!strcmp(option, "--cbv-test")){
      // 0=insurf 1=mri 2=agws 3=maxdist 4=outsurf 5=outval
      // Uses surftype from above
      if(nargc < 5) CMDargNErr(option,3);
      MRIS *surf = MRISread(pargv[0]);
      MRISsaveVertexPositions(surf, ORIGINAL_VERTICES);
      MRI *mri = MRIread(pargv[1]);
      float maxdist;
      sscanf(pargv[3],"%f",&maxdist);
      CBV(pargv[2],mri,surf,surftype,maxdist);
      for(int v=0; v<surf->nvertices;v++){
	surf->vertices[v].x = surf->vertices[v].targx;	
	surf->vertices[v].y = surf->vertices[v].targy;	
	surf->vertices[v].z = surf->vertices[v].targz;	
      }
      int err = MRISwrite(surf,pargv[4]);
      if(err) exit(1);
      MRI *val = MRIcopyMRIS(NULL, surf, 0, "val");
      err = MRIwrite(val,pargv[5]);
      if(err) exit(1);
      MRIfree(&val);
      printf("\nmyfv %s -f %s:overlay=%s -f %s:overlay=%s:edgecolor=green\n\n",
	     pargv[1],pargv[0],pargv[5],pargv[4],pargv[5]);
      exit(0);
      nargsused = 3;
    }
    else  {
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

void print_usage(void)
{
  printf("cpptester \n");
  printf("  --surf inputsurf \n");
  printf("  --curv curvfile  : target curv (will get from surf if not spec)\n");
  printf("  --o outputsurf \n");
  printf("  --out-curv outcurvfile   : write out the curvature of the final surface\n");
  printf("  --init-curv initcurvfile : write out the curvature of the input surface\n");
  printf("  --vol intensityvol adgws \n");
  printf("  --wI wI : intensity weight (default is %lf)\n",wI);
  printf("  --wS wS : spring weight (default is %lf)\n",wS);
  printf("  --wH wH : meancurv weight (default is %lf)\n",wH);
  printf("  --wT wTarg : weight for target surf (default is %lf)\n",wTarg);
  printf("  --cbv-dist dist : (default is %lf)\n",cbvdist);
  printf("  --white or --pial \n");
  printf("  --target targetsurf : CBV target surface\n");
  printf("  --iters iters (defualt %d) \n",iters);
  printf("  --save-iters \n");
  printf("  --threads threads \n");
  printf("  --label labelfile\n");
  printf("  --cbv-test 0=insurf 1=mri 2=agws 3=maxdist 4=targetsurf 5=targval\n");
  printf("  --noise NoiseLevel : add uniform noise +/-Level in normal\n");
  printf("  --debug-vertex vno \n");
}

void print_help(void)
{
  print_usage() ;
  exit(1) ;
}
/* ------------------------------------------------------ */
void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  //fprintf(fp,"cwd %s\n",cwd);
  //fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"surffile  %s\n",surffile);
  if(curvfile) fprintf(fp,"curvfile  %s\n",curvfile);
  fprintf(fp,"adwgs %s\n",adgwsfile);
  fprintf(fp,"mri  %s\n",imrifile);
  fprintf(fp,"outsurf     %s\n",outsurffile);
  fprintf(fp,"surftype  %d (w=%d, p=%d)\n",surftype,GRAY_WHITE,GRAY_CSF);
  fprintf(fp,"wI  %g\n",wI);
  fprintf(fp,"wS  %g\n",wS);
  fprintf(fp,"wH  %g\n",wH);
  fprintf(fp,"wT  %g\n",wTarg);
  fprintf(fp,"cbvdist  %g\n",cbvdist);
  fprintf(fp,"NoiseLevel  %g\n",NoiseLevel);
  fprintf(fp,"threads  %d\n",threads);
  return;
}
void check_options(void){
  if(imrifile == NULL && (wI > 0 || wTarg > 0)){
    printf("ERROR: must specify an input mri\n");
    exit(1);
  }
  if(surffile == NULL){
    printf("ERROR: must specify an input surface\n");
    exit(1);
  }
  if(outsurffile == NULL){
    printf("ERROR: must specify an output surface\n");
    exit(1);
  }
}
