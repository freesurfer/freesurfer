#ifndef __kvlTetrahedronAspectRatio_h
#define __kvlTetrahedronAspectRatio_h

#include "kvlAtlasMesh.h"


namespace kvl
{


//
// From VTK's Graphics/vtkMeshQuality.cxx
// and stuff in ./Common/vtkMath.h
//

static double Dot(const double x[3], const double y[3])
{
  return (x[0]*y[0] + x[1]*y[1] + x[2]*y[2]);
};


static void Cross(const double x[3], const double y[3], double z[3])
{
  double Zx = x[1]*y[2] - x[2]*y[1];
  double Zy = x[2]*y[0] - x[0]*y[2];
  double Zz = x[0]*y[1] - x[1]*y[0];
  z[0] = Zx; z[1] = Zy; z[2] = Zz;
};


static double Norm(const double x[3])
{
  return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
};

// 
static double Determinant3x3(const double c1[3],
                             const double c2[3],
                             const double c3[3])
{
  return c1[0]*c2[1]*c3[2] + c2[0]*c3[1]*c1[2] + c3[0]*c1[1]*c2[2] -
         c1[0]*c3[1]*c2[2] - c2[0]*c1[1]*c3[2] - c3[0]*c2[1]*c1[2];
};


#if 0
static double TetrahedronAspectRatio( const AtlasMesh::PointType&  p0,
                                      const AtlasMesh::PointType&  p1,
                                      const AtlasMesh::PointType&  p2,
                                      const AtlasMesh::PointType&  p3 )
{

  double ab[3],bc[3],ac[3],ad[3],bd[3],cd[3];
  double t0,t1,t2,t3,t4,t5;
  double ma,mb,mc,hm;
  const double normal_coeff = sqrt(6.) / 12.;

  ab[0] = p1[0]-p0[0];
  ab[1] = p1[1]-p0[1];
  ab[2] = p1[2]-p0[2];

  bc[0] = p2[0]-p1[0];
  bc[1] = p2[1]-p1[1];
  bc[2] = p2[2]-p1[2];

  ac[0] = p2[0]-p0[0];
  ac[1] = p2[1]-p0[1];
  ac[2] = p2[2]-p0[2];

  ad[0] = p3[0]-p0[0];
  ad[1] = p3[1]-p0[1];
  ad[2] = p3[2]-p0[2];

  bd[0] = p3[0]-p1[0];
  bd[1] = p3[1]-p1[1];
  bd[2] = p3[2]-p1[2];

  cd[0] = p3[0]-p2[0];
  cd[1] = p3[1]-p2[1];
  cd[2] = p3[2]-p2[2];

  t0 = Dot(ab,ab);
  t1 = Dot(bc,bc);
  t2 = Dot(ac,ac);
  t3 = Dot(ad,ad);
  t4 = Dot(bd,bd);
  t5 = Dot(cd,cd);

  ma = t0 > t1 ? t0 : t1;
  mb = t2 > t3 ? t2 : t3;
  mc = t4 > t5 ? t4 : t5;
  hm = ma > mb ? ma : mb;
  hm = hm > mc ? sqrt(hm) : sqrt(mc);

  Cross(ab,bc,bd);
  Norm(bd);
  Cross(ab,ad,bd);
  Norm(bd);
  Cross(ac,ad,bd);
  Norm(bd);
  Cross(bc,cd,bd);
  Norm(bd);

  t4 = fabs(Determinant3x3(ab,ac,ad));

  return normal_coeff * hm * (t0 + t1 + t2 + t3) / t4;
};
#endif

static double TetrahedronRadiusRatio( const AtlasMesh::PointType&  p0,
                                      const AtlasMesh::PointType&  p1,
                                      const AtlasMesh::PointType&  p2,
                                      const AtlasMesh::PointType&  p3 )
{
  double ab[3],bc[3],ac[3],ad[3],bd[3],cd[3],u[3];
  double abc,abd,acd,bcd,a,b,c,det;
  const double normal_coeff = 1. / 12.;

  ab[0] = p1[0]-p0[0];
  ab[1] = p1[1]-p0[1];
  ab[2] = p1[2]-p0[2];

  bc[0] = p2[0]-p1[0];
  bc[1] = p2[1]-p1[1];
  bc[2] = p2[2]-p1[2];

  ac[0] = p2[0]-p0[0];
  ac[1] = p2[1]-p0[1];
  ac[2] = p2[2]-p0[2];

  ad[0] = p3[0]-p0[0];
  ad[1] = p3[1]-p0[1];
  ad[2] = p3[2]-p0[2];

  bd[0] = p3[0]-p1[0];
  bd[1] = p3[1]-p1[1];
  bd[2] = p3[2]-p1[2];

  cd[0] = p3[0]-p2[0];
  cd[1] = p3[1]-p2[1];
  cd[2] = p3[2]-p2[2];

  a = sqrt(Dot(ab,ab) * Dot(cd,cd));
  b = sqrt(Dot(ac,ac) * Dot(bd,bd));
  c = sqrt(Dot(ad,ad) * Dot(bc,bc));

  Cross(ab,bc,u);
  abc = Norm(u);
  Cross(ab,ad,u);
  abd = Norm(u);
  Cross(ac,ad,u);
  acd = Norm(u);
  Cross(bc,cd,u);
  bcd = Norm(u);

  det = Determinant3x3(ab,ac,ad);

  return normal_coeff * sqrt((a+b+c) * (a+b-c) * (a+c-b) * (b+c-a)) * (abc + abd + acd + bcd) \
    / (det * det);
};






} // end namespace kvl


#endif
