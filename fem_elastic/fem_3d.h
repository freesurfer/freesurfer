
#ifndef H_FEM_3D_H
#define H_FEM_3D_H

#include <list>

#include "mesh.h"

// Delaunay dependence
#include "tetgen.h"

#include "toctree.hpp"

typedef TCoords<int,3> tIntCoords;
typedef TCoords<double,3> tDblCoords;

tIntCoords icoords(int x, int y, int z);

class Material3d : public VMaterial
{
public:
  Material3d();
  Material3d(double E, double nu);

  virtual SmallMatrix get_matrix() const;
};

class Element3d : public TElement<3>
{
public:

  Element3d();

  virtual SmallMatrix get_matrix() const;
  virtual bool dst_contains(const tCoords& c) const;
  virtual tCoords inv_img(const tDblCoords& dst_coords) const;
  virtual double dst_volume() const;

  virtual bool src_contains(const tDblCoords& c) const;
  virtual tCoords dir_img(const tDblCoords& src_coords) const;
  virtual double src_volume() const;

  virtual bool orientation_pb(Frame f=both) const;
  virtual bool orientation_test(double dalpha) const;

  // this next function will update coeffs for the shape function
  void  update_coefs() const;

  virtual void print(std::ostream& os) const;

  virtual double shape_fct(int node_id, const tCoords& pt) const;
private:
  bool contains(const tDblCoords& c1, const tDblCoords& c2,
                const tDblCoords& c3, const tDblCoords& c4,
                const tCoords& ci) const;
  double vol(const tDblCoords& c1, const tDblCoords& c2,
             const tDblCoords& c3, const tDblCoords& c4) const;

  //
  // the shape function will be
  // a + bx + cy + dz
  mutable TCoords<double,4> m_ca, m_cb, m_cc, m_cd;
  mutable bool m_isShapeUpdated;

  //
  // interpol coefs
  // based on the previous, this assembles the shape functions of the elements
  mutable tCoords m_int_a, m_int_b, m_int_c, m_int_d;
  void update_interpol_coefs() const;
  mutable bool m_isInterpolUpdated;

  double det( const tCoords& c1,
              const tCoords& c2,
              const tCoords& c3) const;
  double det(double, double, double,
             double, double, double,
             double, double, double) const;
};


//------------------------------------------------
//
// class that returns pointers to implemented objects
//
//------------------------------------------------

class Constructor
{
public:
  Constructor()
  {}

  static Material3d* material(std::string strLabel,
                              double d_e,
                              double d_nu)
  {
    Material3d* pmat = new Material3d(d_e, d_nu);
    pmat->set_label(strLabel);
    return pmat;
  }
  static TNode<3>* node( int id,
                         tDblCoords coords,
                         tDblCoords delta,
                         const bool *is_active)
  {
    TNode<3>* pnode = new TNode<3>;
    pnode->set_id(id);
    pnode->set_coords(coords);
    for (int i=0; i<3; ++i)
    {
      pnode->set_dof_val( size_t(i), delta(i));
      pnode->set_active( size_t(i), is_active[i]);
    }
    return pnode;
  }
  static Element3d* elt( int id,
                         const VMaterial* cpmat,
                         std::vector<TNode<3>*>& vpnodes)
  {
    Element3d* pelt = new Element3d;
    pelt->set_id(id);
    pelt->set_material(cpmat);
    for ( std::vector<TNode<3>*>::iterator it = vpnodes.begin();
          it != vpnodes.end();
          ++it )
      pelt->add_node( *it );

    return pelt;
  }
  static Element3d* elt( int id,
                         TNode<3>* pn1,
                         TNode<3>* pn2,
                         TNode<3>* pn3,
                         TNode<3>* pn4)
  {
    Element3d* pelt = new Element3d;
    pelt->set_id(id);
    pelt->add_node(pn1);
    pelt->add_node(pn2);
    pelt->add_node(pn3);
    pelt->add_node(pn4);
    return pelt;
  }
};

typedef TMesh<Constructor, 3> TMesh3d;

// forward declaration
// definition in fem_3d.cpp
class ElementProxy;

class CMesh3d : public TMesh3d
{
public:
  typedef TMesh3d Superclass;
  typedef Superclass::tNode tNode;
  typedef Superclass::tElement tElement;
  typedef Superclass::tCoords tCoords;

  CMesh3d();
  CMesh3d(const CMesh3d&);
  CMesh3d& operator=(const CMesh3d&);

  ~CMesh3d();

  const tNode* closest_node(const tCoords& c) const;
  const tElement* element_at_point(const tCoords& c) const;

  tNode* closest_node(const tCoords& c);
  tElement* element_at_point(const tCoords& c);

  // this function's implementation actually uses an octree
  //
  //      since this is a virtual function, the argument will
  //      be set via a member variable
  unsigned int m_maxNodes;
  int build_index_src();

protected:

private:
  typedef toct::Octree<ElementProxy,3> OctreeType;
  OctreeType* m_poctree;
  std::vector<ElementProxy> m_vpEltBlock;
};


class DelaunayMesh
{
public:
  typedef std::list<tDblCoords> PointsListType;

  DelaunayMesh(PointsListType& sp,
               tDblCoords cmin,
               tDblCoords cmax,
               double dEltVol,
               double de, double dnu);
  virtual ~DelaunayMesh()
  {};
  CMesh3d* get();

protected:
  virtual tetgenio* createDelaunay();

  const PointsListType& surfPoints;

  double dEltVol;
  tDblCoords m_cmin, m_cmax;
  double m_de, m_dnu;

  void setupFacet(tetgenio::facet*,
                  int, int, int, int);

private:
  void convertFormat(tetgenio*);

  CMesh3d* m_pmesh;
};

#if 0
class AdaptiveDelaunay : public DelaunayMesh
{
public:
  typedef DelaunayMesh Superclass;
  typedef Superclass::PointsListType PointsListType;

  AdaptiveDelaunay(PointsListType& sp,
                   tDblCoords cmin,
                   tDblCoords cmax,
                   double dEltVol,
                   double de, double dnu,
                   const char* fname);
protected:
  virtual tetgenio* createDelaunay();
private:
  const char* m_baseName;
  double m_volBoundary;

  tetgenio* create_bkg_mesh();
};

#endif

class MeshCreator
{
public:
  MeshCreator(tDblCoords cmin,
              tDblCoords cmax,
              tIntCoords ticks,
              double de, double dnu);

  CMesh3d* get(); // main method

private:
  tDblCoords m_cmin, m_cmax;

  typedef std::map<tIntCoords, TNode<3>*, cless<int,3> > NodeMapType;
  NodeMapType m_mnodes;
  tIntCoords m_ticks;
  tDblCoords m_step;
  CMesh3d* m_pmesh;
  bool m_is_active[3];
  tDblCoords m_d3_zero;

  double m_de, m_dnu; // material constants -> to be set after
  //elements are created

  void add_node(int id, int x, int y, int z); // %2==1 => centroid

  void setup_nodes();
  void setup_elts();

  void setup_elt_group(int& id,
                       TNode<3>* pcubeCenter,
                       TNode<3>* pfaceCenter,
                       TNode<3>* p0,
                       TNode<3>* p1,
                       TNode<3>* p2,
                       TNode<3>* p3);

  TNode<3>* get_node(int x,int y, int z)
  {
    NodeMapType::iterator it
    = m_mnodes.find( icoords(x,y,z) );
    if ( it==m_mnodes.end() )
    {
      std::cerr << " MeshCreator::setup_elts -> could not locate index "
      << icoords(x,y,z) << std::endl;
      exit(1);
    }
    return it->second;
  }

};

#endif // H_FEM_3D_H
