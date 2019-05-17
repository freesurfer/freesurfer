
#include <unistd.h>

#include <math.h>
#include <map>
#include <string>
#include <sstream>

#include "misc.h"
#include "timer.h"

#include "fem_3d.h"

using namespace std;

#define CCS(str) (const_cast<char*>(str.c_str()))

//-------------------------------
//
// Misc
//
//-------------------------------

double sign(double x)
{
  if ( x>1e-13)
    return 1.0;
  else if ( x<-1e-13)
    return -1.0;
  return .0;
}

tIntCoords icoords(int x, int y, int z)
{
  tIntCoords retval;
  retval(0) =x;
  retval(1) =y;
  retval(2) =z;
  return retval;
}

//---------------------------------------
//
// class Material3d
//
//--------------------------------------

//----------
// Definitions
//------------

Material3d::Material3d()
    : VMaterial()
{}

Material3d::Material3d(double E,
                       double nu)
    : VMaterial(E,nu)
{}

SmallMatrix
Material3d::get_matrix() const
{
  SmallMatrix m(6,6);
  m.set(0);

  double dfactor = (get_E()/(1+get_nu()))/(1-2.0*get_nu());
  double alpha = (1.0-get_nu())*dfactor;
  double beta  = (1-2*get_nu())/2 * dfactor;
  double dnu   = dfactor * get_nu();

  m(0,0) = alpha;
  m(0,1) = dnu;
  m(0,2) = dnu;
  m(1,0) = dnu;
  m(1,1) = alpha;
  m(1,2) = dnu;
  m(2,0) = dnu;
  m(2,2) = alpha;

  m(3,3) = beta;
  m(4,4) = beta;
  m(5,5) = beta;

  return m;
}

//-----------------------------------------------------
//
// class Element3d
//
// implements virtual members of the TElement class
// - get_matrix
// - dst_contains
// - src_contains
// - dir_img
// - src_vol
// - orientation_pb
//
//------------------------------------------------------

//----------
// Declarations
//-------------


//------
// Definitions
//------------

Element3d::Element3d()
    : TElement<3>(),
    m_isShapeUpdated(false),
    m_isInterpolUpdated(false)
{}

SmallMatrix
Element3d::get_matrix() const
{
  assert( m_vpNodes.size() == size_t(4) );

  if ( !m_isShapeUpdated )
    update_coefs();

  double dbuf;
  double dfactor = m_cpMaterial->get_E() /
                   (1+m_cpMaterial->get_nu()) /
                   (1 - 2* m_cpMaterial->get_nu()); // used in front of the elasticity matrix
  double dalpha =
    1- m_cpMaterial->get_nu(); // first 3 diagonal elements
  // of the elasticity matrix
  double dbeta =
    .5 * ( 1 - 2* m_cpMaterial->get_nu() ); // next 3 diagonal elements
  // of the elasticity matrix
  double dnu = m_cpMaterial->get_nu();

//  double dvol = abs( vol( m_vpNodes[0]->coords(),
  //   m_vpNodes[1]->coords(),
  //   m_vpNodes[2]->coords(),
  //   m_vpNodes[3]->coords() ) );
  //dfactor *= dvol;

  SmallMatrix mret(12,12);
  for (int i=0; i<4; i++)
    for (int j=0; j<4; j++)
    {
      // compute local matrix

      SmallMatrix m(3,3);

      dbuf = dfactor * (dalpha * m_cb(i)*m_cb(j) + 
                        dbeta * m_cc(i)*m_cc(j) + dbeta * m_cd(i)*m_cd(j) );
      m(0,0) = dbuf;

      dbuf = dfactor * ( dnu * m_cb(i)*m_cc(j) + dbeta * m_cc(i)*m_cb(j) );
      m(0,1) = dbuf;

      dbuf = dfactor * ( dnu * m_cb(i)*m_cd(j) + dbeta *m_cd(i)*m_cb(j) );
      m(0,2) = dbuf;

      dbuf = dfactor * ( dnu* m_cb(j)*m_cc(i) + dbeta * m_cb(i)*m_cc(j) );
      m(1,0) = dbuf;

      dbuf = dfactor * ( dalpha * m_cc(i)*m_cc(j) + 
                         dbeta * m_cb(i)*m_cb(j) + dbeta * m_cd(i)*m_cd(j) );
      m(1,1) = dbuf;

      dbuf = dfactor * ( dnu *m_cc(i)*m_cd(j) + dbeta * m_cd(i)*m_cc(j) );
      m(1,2) = dbuf;

      dbuf = dfactor * ( dnu *m_cd(i)*m_cb(j) + dbeta * m_cb(i)*m_cd(j) );
      m(2,0) = dbuf;

      dbuf = dfactor * ( dnu * m_cd(i)*m_cc(j) + dbeta * m_cc(i)*m_cd(j) );
      m(2,1) = dbuf;

      dbuf = dfactor * ( dalpha * m_cd(i)*m_cd(j) + dbeta * m_cc(i)*m_cc(j) + 
                         dbeta * m_cb(i)*m_cb(j) );
      m(2,2) = dbuf;

      mret.set_block( m,
                      i*3,
                      j*3);
    }

  return mret;

}

double
Element3d::det( const tDblCoords& c0,
                const tDblCoords& c1,
                const tDblCoords& c2) const
{
  return c0(0) * ( c1(1)*c2(2) - c1(2)*c2(1) )
         - c1(0) * ( c0(1)*c2(2) - c2(1)*c0(2) )
         + c2(0) * ( c0(1)*c1(2) - c1(1)*c0(2) );
}

double
Element3d::det( double m00, double m01, double m02,
                double m10, double m11, double m12,
                double m20, double m21, double m22) const
{
  return m00 * ( m11*m22 - m12*m21 )
         - m10 * ( m01*m22 - m21*m02 )
         + m20 * ( m01*m12 - m11*m02 );
}
#if 0
void
Element3d::update_coefs() const
{
  assert(m_vpNodes.size()==size_t(4));

  tDblCoords c[4];
  for (int i=0; i<4; ++i) c[i] = m_vpNodes[i]->coords();

  double dvol = det(c[1],c[2],c[3])
                - det(c[0],c[2],c[3])
                + det(c[0],c[1],c[3])
                - det(c[0],c[1],c[2]);

  int alpha = -1;
  int arg_1, arg_2, arg_3;
  tDblCoords c_1(1.0), c_2(1.0), c_3(1.0);
  // linear algebra - solve the system for the shape fct coefs
  for ( int i=0; i<4; ++i)
  {
    arg_1 = i>0?0:1;
    arg_2 = i>1?1:2;
    arg_3 = i>2?2:3;

    c_1.set(1.0);
    c_2.set(1.0);
    c_3.set(1.0);
    c_1(1) = c[arg_1](1);
    c_1(2) = c[arg_1](2);
    c_2(1) = c[arg_2](1);
    c_1(2) = c[arg_2](2);
    c_3(1) = c[arg_3](1);
    c_3(2) = c[arg_3](2);
    m_cb(i) = alpha * det(c_1,c_2,c_3) / dvol;
    alpha *= -1;

    c_1.set(1.0);
    c_2.set(1.0);
    c_3.set(1.0);
    c_1(1) = c[arg_1](0);
    c_1(2) = c[arg_1](2);
    c_2(1) = c[arg_2](0);
    c_2(2) = c[arg_2](2);
    c_3(1) = c[arg_3](0);
    c_3(2) = c[arg_3](2);
    m_cc(i) = alpha * det(c_1,c_2,c_3) / dvol;
    alpha *= -1;

    c_1.set(1.0);
    c_2.set(1.0);
    c_3.set(1.0);
    c_1(1) = c[arg_1](0);
    c_1(2) = c[arg_1](1);
    c_2(1) = c[arg_2](0);
    c_2(2) = c[arg_2](1);
    c_3(1) = c[arg_3](0);
    c_3(2) = c[arg_3](1);
    m_cd(i) = alpha * det(c_1,c_2,c_3) / dvol;

    m_ca(i) = alpha * det( c[arg_1], c[arg_2], c[arg_3]) / dvol;
  }
  m_updated = true;
}
#else
void
Element3d::update_coefs() const
{
  assert(m_vpNodes.size()==size_t(4));

  double dvol = det( m_vpNodes[1]->coords(),
                     m_vpNodes[2]->coords(),
                     m_vpNodes[3]->coords() )
                - det( m_vpNodes[0]->coords(),
                       m_vpNodes[2]->coords(),
                       m_vpNodes[3]->coords() )
                + det( m_vpNodes[0]->coords(),
                       m_vpNodes[1]->coords(),
                       m_vpNodes[3]->coords() )
                - det( m_vpNodes[0]->coords(),
                       m_vpNodes[1]->coords(),
                       m_vpNodes[2]->coords() );

  int alpha = -1;
  int arg_1, arg_2, arg_3;
  for ( int i=0; i<4; ++i)
  {
    arg_1 = i>0?0:1;
    arg_2 = i>1?1:2;
    arg_3 = i>2?2:3;

    m_cb(i) = alpha * det( 1,
                           m_vpNodes[arg_1]->coords()(1),
                           m_vpNodes[arg_1]->coords()(2),
                           1,
                           m_vpNodes[arg_2]->coords()(1),
                           m_vpNodes[arg_2]->coords()(2),
                           1,
                           m_vpNodes[arg_3]->coords()(1),
                           m_vpNodes[arg_3]->coords()(2) )
      / dvol;
    alpha *= -1;

    m_cc(i) = alpha * det( 1,
                           m_vpNodes[arg_1]->coords()(0),
                           m_vpNodes[arg_1]->coords()(2),
                           1,
                           m_vpNodes[arg_2]->coords()(0),
                           m_vpNodes[arg_2]->coords()(2),
                           1,
                           m_vpNodes[arg_3]->coords()(0),
                           m_vpNodes[arg_3]->coords()(2) )
      / dvol;
    alpha *= -1;

    m_cd(i) = alpha * det( 1,
                           m_vpNodes[arg_1]->coords()(0),
                           m_vpNodes[arg_1]->coords()(1),
                           1,
                           m_vpNodes[arg_2]->coords()(0),
                           m_vpNodes[arg_2]->coords()(1),
                           1,
                           m_vpNodes[arg_3]->coords()(0),
                           m_vpNodes[arg_3]->coords()(1) )
      / dvol;
    alpha *= -1;

    m_ca(i) = alpha * det( m_vpNodes[arg_1]->coords(),
                           m_vpNodes[arg_2]->coords(),
                           m_vpNodes[arg_3]->coords() ) / dvol;

  }
  m_isShapeUpdated = true;

}
#endif

void
Element3d::update_interpol_coefs() const
{
  assert(m_vpNodes.size()==size_t(4));

  if ( !m_isShapeUpdated )
    update_coefs();

  m_int_a.set(0.0);
  m_int_b.set(0.0);
  m_int_c.set(0.0);
  m_int_d.set(0.0);

  for (int i=0; i<4; ++i)
  {
    const Coords3d& delta = m_vpNodes[i]->delta();
    m_int_a += m_ca(i) * delta;
    m_int_b += m_cb(i) * delta;
    m_int_c += m_cc(i) * delta;
    m_int_d += m_cd(i) * delta;
  }
  m_isInterpolUpdated = true;
}

void
Element3d::print(std::ostream& os) const
{
  TElement<3>::print(os);
  std::cout << " a = " << m_ca << std::endl
  << " b = " << m_cb << std::endl
  << " c = " << m_cc << std::endl
  << " d = " << m_cd << std::endl;
}

bool
Element3d::dst_contains(const tDblCoords& c) const
{
  assert(m_vpNodes.size()==size_t(4));

  return contains( m_vpNodes[0]->dst_coords(),
                   m_vpNodes[1]->dst_coords(),
                   m_vpNodes[2]->dst_coords(),
                   m_vpNodes[3]->dst_coords(),
                   c);
}

tDblCoords
Element3d::inv_img(const tDblCoords& c) const
{
  if ( !m_isInterpolUpdated )
    update_interpol_coefs();

  m_int_b(0) += 1.0;
  m_int_c(1) += 1.0;
  m_int_d(2) += 1.0;

  double inv_sys_det = 1.0 / det( m_int_b, m_int_c, m_int_d);

  if ( abs(inv_sys_det) > 1.0e+4 )
  {
    cerr << " Element3d::inv_img -> tiny det\n";
    exit(1);
  }
  tDblCoords csrc;
  tDblCoords cbuf = c - m_int_a;

  csrc(0) = inv_sys_det * det(cbuf, m_int_c, m_int_d);
  csrc(1) = inv_sys_det * det(m_int_b, cbuf, m_int_d);
  csrc(2) = inv_sys_det * det(m_int_b, m_int_c, cbuf);

  // reset vars to init vals at the end
  m_int_b(0) -= 1.0;
  m_int_c(1) -= 1.0;
  m_int_d(2) -= 1.0;

  return csrc;
}

double
Element3d::dst_volume() const
{
  assert( m_vpNodes.size() == size_t(4) );

  return abs( vol( m_vpNodes[0]->dst_coords(),
                   m_vpNodes[1]->dst_coords(),
                   m_vpNodes[2]->dst_coords(),
                   m_vpNodes[3]->dst_coords() ) );
}

bool
Element3d::src_contains(const tDblCoords& c) const
{
  assert( m_vpNodes.size() == size_t(4) );

  return contains( m_vpNodes[0]->coords(),
                   m_vpNodes[1]->coords(),
                   m_vpNodes[2]->coords(),
                   m_vpNodes[3]->coords(),
                   c);
}

bool
Element3d::contains( const tDblCoords& c1,
                     const tDblCoords& c2,
                     const tDblCoords& c3,
                     const tDblCoords& c4,
                     const tDblCoords& c) const
{
  tDblCoords cmin, cmax;
  cmin = c1;
  cmax = c1;

  cmin = min( cmin, c2 );
  cmin = min( cmin, c3 );
  cmin = min( cmin, c4 );

  cmax = max( cmax, c2 );
  cmax = max( cmax, c3 );
  cmax = max( cmax, c4 );

  for (unsigned int ui=0; ui<3; ++ui)
    if ( c(ui) < cmin(ui) || c(ui) > cmax(ui) )
      return false;

  double dvol_total = abs( vol( c1,c2,c3,c4));
  double dvol_1 = abs( vol( c,c1,c2,c3));
  double dvol_2 = abs( vol( c,c1,c2,c4));
  double dvol_3 = abs( vol( c,c1,c3,c4));
  double dvol_4 = abs( vol( c,c2,c3,c4));

  return ( abs( dvol_total - dvol_1 - dvol_2 - dvol_3 - dvol_4 ) < 1.0e-5 );
}

tDblCoords
Element3d::dir_img(const tDblCoords& src_coords) const
{
  if ( !m_isInterpolUpdated )
    update_interpol_coefs();

  return src_coords + m_int_a
         + src_coords(0) * m_int_b
         + src_coords(1) * m_int_c
         + src_coords(2) * m_int_d;

}

double
Element3d::shape_fct(int node_id,
                     const tCoords& pt) const
{
  if ( !m_isShapeUpdated )
    update_coefs();

  return m_ca(node_id) +
         pt(0) * m_cb(node_id) +
         pt(1) * m_cc(node_id) +
         pt(2) * m_cd(node_id);
}

double
Element3d::src_volume() const
{
  assert( m_vpNodes.size() == size_t(4) );

  return abs( vol( m_vpNodes[0]->coords(),
                   m_vpNodes[1]->coords(),
                   m_vpNodes[2]->coords(),
                   m_vpNodes[3]->coords() ) );
}

bool
Element3d::orientation_pb(Frame f) const
{
  assert ( m_vpNodes.size() == size_t(4) );

  switch (f)
  {
  case both:
  {
    double dvol_ref = vol( m_vpNodes[0]->coords(),
                           m_vpNodes[1]->coords(),
                           m_vpNodes[2]->coords(),
                           m_vpNodes[3]->coords() );
    double dvol_def = vol( m_vpNodes[0]->dst_coords(),
                           m_vpNodes[1]->dst_coords(),
                           m_vpNodes[2]->dst_coords(),
                           m_vpNodes[3]->dst_coords() );

    if ( sign(dvol_def) != sign(dvol_ref) )
      return true;
  }
  break;
  case src:
    if ( sign( vol(m_vpNodes[0]->coords(),
                   m_vpNodes[1]->coords(),
                   m_vpNodes[2]->coords(),
                   m_vpNodes[3]->coords() ) ) <= 0)
      return true;
    break;
  case dst:
    if ( sign( vol(m_vpNodes[0]->dst_coords(),
                   m_vpNodes[1]->dst_coords(),
                   m_vpNodes[2]->dst_coords(),
                   m_vpNodes[3]->dst_coords() ) ) <= 0)
      return true;
    break;
  default:
    ;
  }

  return false;
}

bool
Element3d::orientation_test(double dalpha) const
{
  double dvol_ref = vol( m_vpNodes[0]->coords(),
                         m_vpNodes[1]->coords(),
                         m_vpNodes[2]->coords(),
                         m_vpNodes[3]->coords()
                       );
  double dvol_def =
    vol( m_vpNodes[0]->coords() + m_vpNodes[0]->delta() * dalpha,
         m_vpNodes[1]->coords() + m_vpNodes[1]->delta() * dalpha,
         m_vpNodes[2]->coords() + m_vpNodes[2]->delta() * dalpha,
         m_vpNodes[3]->coords() + m_vpNodes[3]->delta() * dalpha
      );
  if ( sign(dvol_def) != sign(dvol_ref) )
    return true;

  return false;
}

double
Element3d::vol( const tDblCoords& c1, const tDblCoords& c2,
                const tDblCoords& c3, const tDblCoords& c4) const
{
  return det( c2-c1,
              c3-c1,
              c4-c1);
}

//------------------------------------------
//
// CMesh3d
//
//------------------------------------------



class ElementProxy : public toct::TCubicRegion<3>
{
public:
  typedef TElement<3> Element;
  typedef TCoords<double,3> Coords3d;

  ElementProxy(CMesh3d* pmesh, unsigned int elt_idx)
  {
    m_pelt = pmesh->fetch_elt(elt_idx);
    m_pelt->src_box( this->m_cmin, this->m_cmax );
  }

  bool contains(const Coords3d& pt) const
  {
    return m_pelt->src_contains(pt);
  }
  Element* elt() const
  {
    return m_pelt;
  }
private:
  Element* m_pelt;
};


//---------------------------------------
//
// CMesh3d class implementation
//
//---------------------

CMesh3d::CMesh3d()
    : TMesh3d(),
    m_maxNodes(20),
    m_poctree(NULL)
{}

CMesh3d::CMesh3d(const CMesh3d& cmesh)
    : TMesh3d(cmesh), m_maxNodes(cmesh.m_maxNodes)
{
  if ( cmesh.m_poctree )
  {
    m_poctree = new OctreeType( *cmesh.m_poctree );
  }
}

CMesh3d::~CMesh3d()
{
  if ( m_poctree ) delete m_poctree;
}

int
CMesh3d::build_index_src()
{

  std::cout << " building index src\n";
  m_vpEltBlock.clear();

  if (m_poctree) delete m_poctree;
  Coords3d cbuf = m_cmax - m_cmin;
  //  double ddelta = std::max( cbuf(0), std::max( cbuf(1), cbuf(2)));
  cbuf = m_cmin + cbuf;
  m_poctree = new OctreeType( m_cmin, cbuf,6, m_maxNodes );
  m_vpEltBlock.reserve( this->get_no_elts() );
  for (unsigned int ui=0, noItems = this->get_no_elts(); ui < noItems; ++ui)
    m_vpEltBlock.push_back( ElementProxy(this, ui) );

  std::cout << " done building the list\n";
  unsigned int count = 0;//, oldPercentage = 0, percentage;
  Timer timer;
  for (std::vector<ElementProxy>::const_iterator cit = m_vpEltBlock.begin();
       cit != m_vpEltBlock.end(); ++cit, ++count )
  {
    m_poctree->insertItem( &*cit );
    if ( !(count % 100000) )
    {
      std::cout << "\t count inserted = " << count
      << " elapsed = " << timer.seconds() << " seconds "
      << " element count = " << m_poctree->getElementCount()
      << std::endl;
      timer.reset();
    }
  }

  std::cout << " done building octree - total elements = "
  << m_poctree->getElementCount() << std::endl;

  return 0;
}

const CMesh3d::tElement*
CMesh3d::element_at_point(const tCoords& c) const
{
  if ( const ElementProxy* cep = m_poctree->element_at_point(c) )
    return cep->elt();
  return NULL;
}

CMesh3d::tElement*
CMesh3d::element_at_point(const tCoords& c)
{
  if (const ElementProxy* cep = m_poctree->element_at_point(c) )
    return cep->elt();
  return NULL;
}



const CMesh3d::tNode*
CMesh3d::closest_node(const tCoords& c) const
{
  // get element at point
  const tElement* cpelt = this->element_at_point(c);

  if ( !cpelt ) return NULL;

  double dbuf, dMin = ( m_cmax - m_cmin ).norm();
  tNode* pnode = NULL;
  tNode* pargmin = NULL;
  for ( unsigned int ui=0, nnodes = cpelt->no_nodes();
        ui < nnodes; ++ui )
  {
    cpelt->get_node(ui, &pnode);
    dbuf = (c - pnode->coords()).norm();
    if ( dbuf < dMin )
    {
      pargmin = pnode;
      dMin = dbuf;
    }
  } // next ui

  return pargmin;
}


CMesh3d::tNode*
CMesh3d::closest_node(const tCoords& c)
{
  // get element at point
  tElement* pelt = this->element_at_point(c);

  if ( !pelt ) return NULL;

  double dbuf, dMin = ( m_cmax - m_cmin ).norm();
  tNode* pnode = NULL;
  tNode* pargmin = NULL;
  for ( unsigned int ui=0, nnodes = pelt->no_nodes();
        ui < nnodes; ++ui )
  {
    pelt->get_node(ui, &pnode);
    dbuf = (c - pnode->coords()).norm();
    if ( dbuf < dMin )
    {
      pargmin = pnode;
      dMin = dbuf;
    }
  } // next ui

  return pargmin;
}


//-------------------------------------------------------
//
// Function that creates a mesh
//
//-------------------------------------------------------


DelaunayMesh::DelaunayMesh(PointsListType& sp,
                           tDblCoords cmin,
                           tDblCoords cmax,
                           double vol,
                           double de, double dnu)
    : surfPoints(sp),dEltVol(vol), m_cmin(cmin), m_cmax(cmax),
    m_de(de), m_dnu(dnu)
{}

tetgenio* DelaunayMesh::createDelaunay()
{
  tetgenio in;
  tetgenio* out = new tetgenio;

  // All indices start from 1
  in.firstnumber = 1;

  in.numberofpoints = 8 + surfPoints.size();
  in.pointlist = new REAL[in.numberofpoints * 3];

  // fill in the points
  for (unsigned int pt=0; pt<4; ++pt)
    for (unsigned int ui=0; ui<3; ++ui)
      in.pointlist[ui+3*pt] = m_cmin(ui);

  in.pointlist[3] = m_cmax(0);

  in.pointlist[6] = m_cmax(0);
  in.pointlist[7] = m_cmax(1);

  in.pointlist[10] = m_cmax(1);

  for (unsigned int ui=4; ui<8; ++ui)
  {
    in.pointlist[ ui*3 ]    = in.pointlist[ (ui-4) * 3];
    in.pointlist[ ui*3 + 1] = in.pointlist[ (ui-4) * 3 + 1 ];
    in.pointlist[ ui*3 + 2] = m_cmax(2);
  }

  // add the points from the container
  unsigned int ui = 8;
  for ( PointsListType::const_iterator cit = surfPoints.begin();
        cit != surfPoints.end(); ++cit, ++ui)
  {
    for (unsigned int a = 0; a<3; ++a)
      in.pointlist[ui*3 + a] = (*cit)(a);
  }


  // setup facets
  in.numberoffacets = 6;
  in.facetlist = new tetgenio::facet[ in.numberoffacets ];
  in.facetmarkerlist = new int[in.numberoffacets];

  // Facet 1. The leftmost facet.
  setupFacet(&in.facetlist[0],
             1,2,3,4);

  setupFacet(&in.facetlist[1],
             5,6,7,8);

  setupFacet(&in.facetlist[2],
             1,5,6,2);

  setupFacet(&in.facetlist[3],
             2,6,7,3);

  setupFacet(&in.facetlist[4],
             3,7,8,4);

  setupFacet(&in.facetlist[5],
             4,8,5,1);

  for (int i=0; i<in.numberoffacets; ++i)
    in.facetmarkerlist[i] = 0;

  // tetrahedralize
  double volConstraint = 3.0f;
  std::ostringstream os;
  os << "a" << volConstraint;
  os.flush();

  char* pchBuf = new char[100];

  sprintf(pchBuf, "pq1.414a%f",
          dEltVol);

  tetrahedralize( pchBuf,
                  &in, out);

  delete[] pchBuf;

  return out;

}

/*

the input of this function is a Delaunay tetrahedral mesh

- recover nodes
- recover tetrahedra

*/
void
DelaunayMesh::convertFormat(tetgenio* p)
{
  // if destination mesh doesn't exist, exit error
  if ( !m_pmesh )
    throw " DelaunayMesh - NULL mesh";

  // setup nodes
  int firstnumber = p->firstnumber;
  int noPoints = p->numberofpoints;

  tDblCoords tc, dcZero(.0);
  REAL* pdbl = & p->pointlist[0];
  bool is_active[3];
  std::fill_n(is_active,3, false);

  for (int index=0; index<noPoints; ++index, pdbl+=3)
  {
    tc(0) = *pdbl;
    tc(1) = *(pdbl+1);
    tc(2) = *(pdbl+2);

    TNode<3>* pnode = Constructor::node(index,
                                        tc,
                                        dcZero,
                                        is_active);
    m_pmesh->add_node(pnode);
  } // next index

  // setup elements
  int numberoftetrahedra = p->numberoftetrahedra;
  int* pint = &p->tetrahedronlist[0];
  for (int index=0; index< numberoftetrahedra; ++index, pint+=4)
  {
    Element3d* pelt = Constructor::elt( index,
                                        m_pmesh->node(*pint - firstnumber),
                                        m_pmesh->node(*(pint+1) - firstnumber),
                                        m_pmesh->node(*(pint+2) - firstnumber),
                                        m_pmesh->node(*(pint+3) - firstnumber) );
    m_pmesh->add_elt( pelt );
  } // next index, pint
}

CMesh3d* DelaunayMesh::get()
{
  m_pmesh = new CMesh3d; // leak?

  tetgenio* out = this->createDelaunay();
  this->convertFormat( out );

  // dbg - save mesh
  // out->save_nodes((char*)"iteration");
  // out->save_elements((char*)"iteration");
  // out->save_faces((char*)"iteration");

  delete out;

  m_pmesh->set_constants(m_de, m_dnu);
  m_pmesh->m_cmin = m_cmin;
  m_pmesh->m_cmax = m_cmax;

  return m_pmesh;
}

void
DelaunayMesh::setupFacet(tetgenio::facet* f,
                         int a, int b, int c, int d)
{
  tetgenio::polygon* p;

  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;

  p = &f->polygonlist[0];
  p->numberofvertices = 4;

  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = a;
  p->vertexlist[1] = b;
  p->vertexlist[2] = c;
  p->vertexlist[3] = d;
}

//-------------------------------------
//
// Class Delaunay2 - uses the mtr file to have the mesh be refined
// near the surface
//
//-------------------------------------

#if 0
AdaptiveDelaunay::AdaptiveDelaunay(PointsListType& sp,
                                   tDblCoords cmin,
                                   tDblCoords cmax,
                                   double dEltVol,
                                   double de, double dnu,
                                   const char* fname)
    : DelaunayMesh(sp,cmin,cmax,dEltVol,de,dnu), m_baseName(fname),
    m_volBoundary(10.0)
{}

// this function uses files on the hard-disk
// that is why a baseName is needed
//
//
tetgenio*
AdaptiveDelaunay::createDelaunay()
{
  std::cout << "AdaptiveDelaunay - createDelaunay\n";
  tetgenio* box= this->create_bkg_mesh();

  tetgenio* mesh = new tetgenio;
  {
    char* pchBuf = new char[100];
    sprintf(pchBuf, "pVYC");
    tetrahedralize(pchBuf, box,  mesh);
    delete[] pchBuf;
  }

  std::string strMesh = m_baseName;
  std::string strBkgMesh = strMesh + ".b";

  char* meshName = const_cast<char*>(m_baseName);
  std::cout << " saving box mesh - " << meshName << std::endl;
  box->save_nodes( meshName );
  box->save_elements( meshName );
  box->save_faces( meshName );
  box->save_poly( meshName );

  meshName = const_cast<char*>( strBkgMesh.c_str() );
  std::cout << " saving bkg mesh - " << meshName << std::endl;
  mesh->save_nodes( meshName );
  mesh->save_elements( meshName );
  mesh->save_faces( meshName );
  mesh->save_poly( meshName );

  // create the MTR file
  //
  // as an internal convention, the first 8 vertices
  // are the edges of the cube
  std::ofstream ofs( (strBkgMesh+".mtr").c_str() );
  ofs << surfPoints.size() + 8 << " 1\n";
  for (unsigned int ui(0), npts(8); ui<npts; ++ui)
    ofs << m_volBoundary << std::endl;
  for (unsigned int ui(0), npts(surfPoints.size());
       ui < npts; ++ui)
    ofs << dEltVol << std::endl;
  ofs.close();

#if 0
  char* sysCall[4];
  sysCall[0] = "tetgen"; // tetgen needs to be in the path
  char* pchBuf = new char[100];
  sprintf(pchBuf, "tetgen -pqb %s", m_baseName);
  sysCall[1] = pchBuf;
  sysCall[2] = const_cast<char*>(m_baseName);
  sysCall[3] = (char*)0;
#endif

  char* pchBuf = new char[100];
  sprintf(pchBuf, "tetgen -pqb %s", m_baseName);

  int retVal = system( pchBuf );
  std::cout << " ret val from system = " << retVal << std::endl;
  if ( retVal )
    throw std::string(" The call to tetgen failed");
  delete[] pchBuf;

  //------------------------------------

  // 1. cleanup previous meshes
  delete box;
  delete mesh;

  // 2. load created mesh
  pchBuf = new char[100];
  sprintf(pchBuf, "%s.1", m_baseName);
  std::cout << "loading mesh 1304 " << pchBuf << std::endl;
  mesh = new tetgenio;
  mesh->load_tetmesh( pchBuf );
  delete[] pchBuf;

  return mesh;
}

tetgenio*
AdaptiveDelaunay::create_bkg_mesh()
{
  tetgenio *box = new tetgenio;

  box->firstnumber = 1;
  box->numberofpoints = 8 + surfPoints.size();
  box->pointlist = new REAL[ box->numberofpoints * 3 ];

  // fill the first 8 points
  for (unsigned int pt=0; pt<4; ++pt)
    for (unsigned int ui=0; ui<3; ++ui)
      box->pointlist[ui+3*pt] = m_cmin(ui);

  box->pointlist[3] = m_cmax(0);
  box->pointlist[6] = m_cmax(0);
  box->pointlist[7] = m_cmax(1);
  box->pointlist[10]= m_cmax(1);

  for (unsigned int ui=4; ui<8; ++ui)
  {
    box->pointlist[ ui*3 ] = box->pointlist[ (ui-4) * 3] ;
    box->pointlist[ ui*3 + 1] = box->pointlist[ (ui-4) * 3 + 1];
    box->pointlist[ ui*3 + 2] = m_cmax(2);
  }

  // add the points from the surface
  unsigned int ui=8;
  unsigned int counter = 0;
  for ( PointsListType::const_iterator cit = surfPoints.begin();
        cit != surfPoints.end();
        ++cit , ++ui)
  {
    for (unsigned int a=0; a<3; ++a)
      box->pointlist[ui*3+a] = (*cit)(a);
  } // next cit

  // setup facets
  box->numberoffacets = 6;
  box->facetlist = new tetgenio::facet[ box->numberoffacets ];
  box->facetmarkerlist = new int[box->numberoffacets];

  int indices[] = { 1, 2, 3, 4,
                    5, 6, 7, 8,
                    1, 5, 6, 2,
                    2, 6, 7, 3,
                    3, 7, 8, 4,
                    4, 8, 5, 1 };
  for (unsigned int a(0), ai(0); a<6; ++a, ai+=4)
    this->setupFacet(&box->facetlist[a],
                     indices[ai], indices[ai+1],
                     indices[ai+2], indices[ai+3]
                    );

  return box;
}
#endif

MeshCreator::MeshCreator(tDblCoords cmin,
                         tDblCoords cmax,
                         tIntCoords ticks,
                         double de, double dnu)
    : m_cmin(cmin), m_cmax(cmax),
    m_mnodes(),
    m_ticks(ticks),
    m_step(-1.0),
    m_d3_zero(0.0),
    m_de(de), m_dnu(dnu)
{
  std::cout << " Mesh creator ticks = " << ticks << std::endl;

  std::fill_n(m_is_active,3,false);

  for (int i=0; i<3; ++i)
    m_step(i) = (m_cmax(i)-m_cmin(i)) / (m_ticks(i)-1.0) *.5;

  m_pmesh = new CMesh3d;
  m_pmesh->m_cmin = m_cmin;
  m_pmesh->m_cmax = m_cmax;
}

CMesh3d*
MeshCreator::get()
{
  setup_nodes();
  setup_elts();

  return m_pmesh;
}

void
MeshCreator::add_node(int id,
                      int x, int y, int z)
{
  tIntCoords index;
  index(0) = x;
  index(1) = y;
  index(2) = z;
  tDblCoords tc;
  tc(0) = double(x) * m_step(0);
  tc(1) = double(y) * m_step(1);
  tc(2) = double(z) * m_step(2);
  tc += m_cmin;

  TNode<3>* pnode = Constructor::node(id,
                                      tc,
                                      m_d3_zero,
                                      m_is_active);
  m_pmesh->add_node(pnode);
  m_mnodes[index] = pnode;
}

void
MeshCreator::setup_nodes()
{
  int id = 0;

  for ( int i=0; i<m_ticks(0); ++i)
    for ( int j=0; j<m_ticks(1); ++j)
      for ( int k=0; k<m_ticks(2); ++k)
      {
        add_node(id++, 2*i, 2*j, 2*k);

        if ( i<m_ticks(0)-1 && j<m_ticks(1)-1 )
          add_node(id++, 2*i+1, 2*j+1, 2*k);
        if ( i<m_ticks(0)-1 && k<m_ticks(2)-1 )
          add_node(id++, 2*i+1, 2*j, 2*k+1);
        if ( j<m_ticks(1)-1 && k<m_ticks(2)-1 )
          add_node(id++, 2*i, 2*j+1, 2*k+1);
        if ( i<m_ticks(0)-1 && j<m_ticks(1)-1 && k<m_ticks(2)-1 )
          add_node(id++, 2*i+1, 2*j+1, 2*k+1);
      }
}

void
MeshCreator::setup_elt_group(int& id,
                             TNode<3>* pcubeCenter,
                             TNode<3>* pfaceCenter,
                             TNode<3>* p0,
                             TNode<3>* p1,
                             TNode<3>* p2,
                             TNode<3>* p3)
{
  Element3d* pelt = NULL;

  pelt = Constructor::elt( id++,
                           pcubeCenter, pfaceCenter,
                           p0, p1 );
  m_pmesh->add_elt(pelt);

  pelt = Constructor::elt( id++,
                           pcubeCenter, pfaceCenter,
                           p1, p2 );
  m_pmesh->add_elt(pelt);

  pelt = Constructor::elt( id++,
                           pcubeCenter, pfaceCenter,
                           p2, p3 );
  m_pmesh->add_elt(pelt);

  pelt = Constructor::elt( id++,
                           pcubeCenter, pfaceCenter,
                           p3, p0);
  m_pmesh->add_elt(pelt);

}

void
MeshCreator::setup_elts()
{
  int id=0;
  NodeMapType::iterator it;
  TNode<3> *pn0, *pn1, *pn2, *pn3,
  *pn4, *pn5, *pn6, *pn7, *pnij0,
  *pnik0, *pnjk0, *pnij1, *pnik1, *pnjk1, *pncenter;

  vector<TNode<3>*> vpnodes;

  for ( int k=0; k<m_ticks(2)-1; ++k)
    for ( int j=0; j<m_ticks(1)-1; ++j)
      for ( int i=0; i<m_ticks(0)-1; ++i)
      {
        // 15 nodes

        pn0 = get_node( 2*i,2*j,2*k);
        pn1 = get_node( 2*(i+1), 2*j, 2*k);
        pn2 = get_node( 2*(i+1), 2*(j+1), 2*k);
        pn3 = get_node( 2*i, 2*(j+1), 2*k);

        pn4 = get_node( 2*i, 2*j, 2*(k+1));
        pn5 = get_node( 2*(i+1), 2*j, 2*(k+1) );
        pn6 = get_node( 2*(i+1), 2*(j+1), 2*(k+1) );
        pn7 = get_node( 2*i, 2*(j+1), 2*(k+1) );

        pnij0 = get_node(2*i+1, 2*j+1, 2*k);
        pnik0 = get_node(2*i+1, 2*j, 2*k+1);
        pnjk0 = get_node(2*i, 2*j+1, 2*k+1);

        pnij1 = get_node( 2*i+1, 2*j+1, 2*(k+1) );
        pnik1 = get_node( 2*i+1, 2*(j+1), 2*k+1 );
        pnjk1 = get_node( 2*(i+1), 2*j+1, 2*k+1 );

        pncenter = get_node( 2*i+1, 2*j+1, 2*k+1);

        // 24 elements
        setup_elt_group(id, pncenter, pnij0,
                        pn3, pn2, pn1, pn0);

        setup_elt_group(id, pncenter, pnij1,
                        pn4, pn5, pn6, pn7);

        setup_elt_group(id, pncenter, pnjk1,
                        pn1, pn2, pn6, pn5 );

        setup_elt_group(id, pncenter, pnik0,
                        pn0, pn1, pn5, pn4 );

        setup_elt_group(id, pncenter, pnjk0,
                        pn3, pn0, pn4, pn7);

        setup_elt_group(id, pncenter, pnik1,
                        pn2, pn3, pn7, pn6);

      }

  m_pmesh->set_constants(m_de, m_dnu);
}

