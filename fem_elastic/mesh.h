
#ifndef H_MESH_H
#define H_MESH_H

#include <cmath>
#include <fstream>
#include <sstream>
#include <list>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>
#include <stack>
#include <stdexcept>

#include "coords.h"
#include "material.h"
#include "element.h"

#include "gmpError.h"

#include "ZlibStringCompressor.h"
#include "streamIoMisc.h"

#include "stl_utils.hpp"

#if 0


#include "fio.h"
;
#endif

typedef std::stack<int> tStack;


//----------------------------------------------------------------------
//
//   Class definitions
//
//----------------------


//--------------------------------------------------------------
//
// will use the 3 creators to create pointers to the specific
// classes used in this case. Since the interface are already
// implemented, each of these classes is required to
// provide a static operator() that will return a pointer to objects.
template<class Cstr,
int n>
class TMesh
{
public:
  typedef TCoords<int,n> tIntCoords;
  typedef TNode<n> tNode;
  typedef TElement<n> tElement;
  typedef TCoords<double, n> tCoords;

  typedef std::vector<tNode*> NodeContainerType;
  typedef std::vector<tElement*> ElementContainerType;
  typedef std::map<std::string, VMaterial*> MaterialContainerType;


  TMesh();
  TMesh(const char*);
  TMesh(const TMesh&);
  TMesh& operator=(const TMesh&);

  virtual ~TMesh();
  void load(const char*);
  void load(std::istream& is);
  void save(const char*) const;
  void save(std::ostream& os) const;


  std::pair<TCoords<double,n>,
  TCoords<double,n> > invert(); // in-place inversion

  // material management

  // label used by default
  void set_constants(double _e,
                     double _nu); // will set constants for all the elements

  void clear_materials(); // will clear the map and reinit it
  void add_material(std::string strLabel, VMaterial* vm);
  void assign_material(unsigned int i, std::string strLabel);

  void set_default_material(VMaterial* pmaterial);

  typedef typename MaterialContainerType::iterator
  MaterialIteratorType;
  typedef typename MaterialContainerType::const_iterator
  MaterialConstIteratorType;
  void get_material_iterators(MaterialIteratorType& begin,
                              MaterialIteratorType& end)
  {
    begin = m_mpMaterial.begin();
    end   = m_mpMaterial.end();
  }

  void get_material_iterators(MaterialConstIteratorType& begin,
                              MaterialConstIteratorType& end) const
  {
    begin = m_mpMaterial.begin();
    end   = m_mpMaterial.end();
  }
  unsigned int get_no_materials() const
  {
    return m_mpMaterial.size();
  }

  //--------------------------

  // nodes management
  typedef typename NodeContainerType::iterator NodeIterator;
  typedef typename NodeContainerType::const_iterator NodeConstIterator;
  void get_node_citer(NodeConstIterator& begin, NodeConstIterator& end)
  {
    begin = m_vpNodes.begin();
    end = m_vpNodes.end();
  }
  void get_node_iter(NodeIterator& begin, NodeIterator& end)
  {
    begin = m_vpNodes.begin();
    end = m_vpNodes.end();
  }

  unsigned int get_no_nodes() const
  {
    return m_vpNodes.size();
  }
  tNode* node(unsigned int i) const
  {
    return m_vpNodes[i];
  }
  bool get_node(unsigned int i, tNode** ppNode) const;
  virtual const tNode* closest_node(const tCoords& c) const = 0;
  virtual const tElement* element_at_point(const tCoords& c) const = 0;

  virtual tNode* closest_node(const tCoords& c) = 0;
  virtual tElement* element_at_point(const tCoords& c) = 0;

  unsigned int add_node(tNode* pnode);
  //---------------------------

  typedef typename ElementContainerType::iterator ElementIterator;
  typedef typename ElementContainerType::const_iterator
  ElementConstIterator;


  // elements management
  int add_elt(tElement* elt);
  unsigned int get_no_elts() const
  {
    return m_vpElements.size();
  }
  const tElement* get_elt(unsigned int i) const
  {
    assert(i<m_vpElements.size());
    return m_vpElements[i];
  }
  tElement* fetch_elt(unsigned int i)
  {
    assert(i<m_vpElements.size());
    return m_vpElements[i];
  }

  // recover all the elements that have a common NODE
  // in common with current one
  typedef std::list<unsigned int> ElementIndexContainer;
  int get_neighboring_elts(unsigned int i, int radius,
                           ElementIndexContainer& eltIndex) const;

  template<class In> void remove_elts(In begin, In end) /*throw(gmpErr)*/;
  bool check_elt_id() const;
  //---------------------------

  // Mapping functionality
  // if signalTopology is true and the element has a topological defect,
  //    an invalid point will be returned
  tCoords        dir_img(const tCoords&,
                         bool signalTopology = false) const /*throw(gmpErr) */;
  void           get_dst_box(tCoords& cmin,
                             tCoords& cmax) const;
  //----------------------------

  void           get_src_box(tCoords& cmin,
                             tCoords& cmax) const
  {
    cmin = m_cmin;
    cmax = m_cmax;
  }

  //----------------------------
  // Some summary orientation checking

  // fast orientation problem checking - stops at the first true
  bool orientation_pb(Frame f=both) const;

  // for this version, the template argument should be an
  // output_iterator to the container at hand
  template <class InsertIterator> void orientation_pb(Frame f,
      InsertIterator ii) const;

  //----------------------------


  virtual int  build_index_src() = 0; // negative return => error

  double src_volume() const;
  double dst_volume() const;

  tCoords m_cmin, m_cmax;

protected:
  void init();
  void clone(const TMesh& mesh); // copy constructor and = operator
  void free();

  NodeContainerType m_vpNodes;
  ElementContainerType m_vpElements;

  MaterialContainerType m_mpMaterial;

};

//---------------------------------------------------------------
//
// class implementations

//--------------------------------------------
//
//   TMesh implementation
//
//--------------------------------------------

template<class Cstr,
int n>
TMesh<Cstr,n>::TMesh()
    :     m_cmin(), m_cmax(),
    m_vpNodes(),
    m_vpElements(),
    m_mpMaterial()
{
  m_cmin.invalidate();
  m_cmax.invalidate();
}

template<class Cstr,
int n>
TMesh<Cstr,n>::TMesh(const char* fname)
    : m_cmin(), m_cmax(),
    m_vpNodes(),
    m_vpElements(),
    m_mpMaterial()
{
  m_cmin.invalidate();
  m_cmax.invalidate();
  load(fname);
}

template<class Cstr,
int n>
TMesh<Cstr,n>
::TMesh(const TMesh<Cstr,n>& mesh)
    :    m_cmin(), m_cmax(),
    m_vpNodes(),
    m_vpElements(),
    m_mpMaterial()
{
  m_cmin.invalidate();
  m_cmax.invalidate();
  clone(mesh);
}

template<class Cstr,
int n>
TMesh<Cstr,n>&
TMesh<Cstr,n>::operator=(const TMesh<Cstr,n>& mesh)
{
  free();
  clone(mesh);
  return *this;
}

template<class Cstr,
int n>
TMesh<Cstr,n>::~TMesh()
{
  free();
}

template<class Cstr,
int n>
void
TMesh<Cstr,n>::free()
{
  // nodes
  std::transform( m_vpNodes.begin(), m_vpNodes.end(),
                  m_vpNodes.begin(),
                  FunctorDeletePointer()
                );
#if 0
  for ( typename NodeContainerType::const_iterator cit = m_vpNodes.begin();
        cit != m_vpNodes.end();
        ++cit )
    delete (*cit);
  m_vpNodes.clear();
#endif

  // elements
  std::transform( m_vpElements.begin(), m_vpElements.end(),
                  m_vpElements.begin(),
                  FunctorDeletePointer()
                );
#if 0
  for ( typename ElementContainerType::const_iterator cit = m_vpElements.begin();
        cit != m_vpElements.end();
        ++cit )
    delete (*cit);
  m_vpElements.clear();
#endif

  // materials
  for ( MaterialConstIteratorType cit = m_mpMaterial.begin();
        cit != m_mpMaterial.end();
        ++cit )
    delete cit->second;
  m_mpMaterial.clear();
}

template<class Cstr,
int n>
void
TMesh<Cstr,n>::load(const char* fname)
{
  // open stream
  std::ifstream ifs(fname, std::ios::binary);
  if ( !ifs )
    throw "TMesh::load - failed opening input stream";

  ifs.seekg(0, std::ios::end);
  unsigned int size = ifs.tellg();
  char* dataBuffer = new char[size];
  ifs.seekg(0, std::ios::beg);
  ifs.read(dataBuffer, size);

  const std::string strCompressed(dataBuffer, size);

  ZlibStringCompressor compressor;

  const std::string strInflated = compressor.inflate( strCompressed );

  // feed the inflated stream to the load function
  std::istringstream is(strInflated, std::ios::binary);

  this->load(is);
  delete[] dataBuffer;
}

template<class Cstr,
int n>
void
TMesh<Cstr,n>::load(std::istream& is)
{
  unsigned int noItems, uilen;

  // assume stream is open in binary mode
  if ( !is )
  {
    throw "TMesh::load - stream isn't open";
  }

  // no materials
  noItems = TRead<unsigned int>(is);

  // read actual materials
  for (unsigned int ui=0; ui<noItems; ++ui)
  {
    double de, dnu;
    uilen = TRead<unsigned int>(is);
    char* buffer = new char[uilen];

    is.read(buffer, uilen);
    de = TRead<double>(is);
    dnu= TRead<double>(is);

    std::string label(buffer, uilen);
    VMaterial *pmaterial = Cstr::material(label, de, dnu);

    m_mpMaterial[ label ] = pmaterial;
    delete[] buffer;
  } // next ui

  // bounding box
  bool flag = TRead<bool>(is);

  if ( flag )
  {
    m_cmin.validate();
    m_cmax.validate();

    for (unsigned int ui=0; ui<n; ++ui)
      m_cmin(ui) = TRead<double>(is);
    for (unsigned int ui=0; ui<n; ++ui)
      m_cmax(ui) = TRead<double>(is);
  }
  else
  {
    m_cmin.invalidate();
    m_cmax.invalidate();
  }

  // nodes
  bool is_active[n];

  typedef typename std::map<int,tNode*> NodeMapType;
//  typedef typename NodeMapType::iterator NodeMapIterator;
  typedef typename NodeMapType::const_iterator NodeMapConstIterator;
  NodeMapType node_map;

  noItems = TRead<unsigned int>(is);
  tCoords coords, delta;
  for (unsigned int ui=0;
       ui < noItems; ++ui )
  {
    int id = TRead<int>(is);

    for (unsigned int uj=0; uj<n; ++uj)
      coords(uj) = TRead<double>(is);

    for (unsigned int uj=0; uj<n; ++uj)
      delta(uj) = TRead<double>(is);

    for (unsigned int uj=0; uj<n; ++uj)
      is_active[uj] = TRead<bool>(is);

    tNode* pnode = Cstr::node(id, coords,
                              delta, is_active);
    node_map[id] = pnode;
    m_vpNodes.push_back(pnode);
  } // next ui

  // elements
  noItems = TRead<unsigned int>(is);
  int node[n+1];
  MaterialConstIteratorType materialIter;


  for (unsigned int ui=0;
       ui < noItems;
       ++ui )
  {
    int id = TRead<int>(is);
    // material label
    uilen = TRead<unsigned int>(is);
    char* buffer = new char[uilen];

    is.read(buffer, uilen);
    std::string materialLabel(buffer, uilen);
    delete[] buffer;

    materialIter = m_mpMaterial.find(materialLabel);
    if ( materialIter == m_mpMaterial.end() )
    {
      std::ostringstream excOs;
      excOs << "TMesh load - invalid material label "
      << materialLabel;
      throw excOs.str().c_str();
    }
    for (unsigned int uj=0; uj<=n; ++uj)
      node[uj] = TRead<int>(is);

    std::vector<tNode*> vpNodes;
    for (unsigned int uj=0; uj<=n; ++uj)
    {
      NodeMapConstIterator  cit
      = node_map.find( node[uj] );
      if ( cit == node_map.end() )
        throw "TMesh::load - invalid element node";

      vpNodes.push_back( cit->second );
    } // next j

    tElement* pelt = Cstr::elt( id, materialIter->second, vpNodes );
    m_vpElements.push_back( pelt );

  } // next ui


}

template<class Cstr,
int n>
void
TMesh<Cstr,n>::clone(const TMesh& mesh)
{
  // this method will assume it receives the object in a virgin state

  int counter = 0;
  // 1. start by adding materials
  MaterialConstIteratorType cbegin, cend;
  mesh.get_material_iterators(cbegin, cend);
  for ( MaterialConstIteratorType cit = cbegin; cit!=cend; ++cit)
  {
    ++counter;
    VMaterial *material = Cstr::material(cit->first, cit->second->get_E(),
                                         cit->second->get_nu() );
    this->add_material( cit->first, material );
  }

  // 3. process nodes
  unsigned int no_of_nodes = mesh.get_no_nodes();
  tNode* pnode;
  tNode* cpnode;
  bool isActive[n];
  for ( unsigned int i=0; i<no_of_nodes; ++i)
  {
    if ( !mesh.get_node(i, &cpnode) )
      std::cerr << "TMesh::clone - didn't find node\n";
    for (int u=0; u<n; ++u) isActive[u] = cpnode->is_dof_active(u);
    pnode = Cstr::node( cpnode->get_id(),
                        cpnode->coords(),
                        cpnode->delta(),
                        isActive );
    m_vpNodes.push_back(pnode);
  }

  // 4. process elements
  unsigned int no_of_elts = mesh.get_no_elts();
  std::vector<tNode*> vpNodes;
  const tElement* cpelt;
  tElement* pelt;
  for (unsigned int i=0; i<no_of_elts; ++i)
  {
    cpelt = mesh.get_elt(i);
    vpNodes.clear();
    for ( int u=0; u<cpelt->no_nodes(); ++u)
    {
      cpelt->get_node(u, &pnode);
      vpNodes.push_back( m_vpNodes[pnode->get_id()] );
    }

    VMaterial* pmat = m_mpMaterial.find( cpelt->material()->label() )->second;
    if ( !pmat )
    {
      std::cerr << " TMesh::clone -> didn't find material\n";
    }
    pelt = Cstr::elt( cpelt->get_id(),
                      pmat,
                      vpNodes );
    m_vpElements.push_back( pelt );
  }

  m_cmin = mesh.m_cmin;
  m_cmax = mesh.m_cmax;
}

template<class Cstr,
int n>
void
TMesh<Cstr,n>::save(const char* fname) const
{
  // open stream
  std::ostringstream os( std::ios::binary);
  if ( !os )
    throw "TMesh::save - failed opening output stream";

  this->save(os);

  ZlibStringCompressor compressor;
  std::string compressed = compressor.compress( os.str(),
                           Z_BEST_COMPRESSION);

  std::ofstream ofs(fname, std::ios::binary);

  ofs.write( compressed.c_str(),
             sizeof(char) * compressed.size() );

  ofs.close();

}

template<class Cstr,
int n>
void
TMesh<Cstr,n>::save(std::ostream& os) const
{

  // materials
  // no items
  TWrite(os, (unsigned int)m_mpMaterial.size() );
  // items
  for ( MaterialConstIteratorType cit = m_mpMaterial.begin();
        cit != m_mpMaterial.end(); ++cit )
  {
    TWrite(os, (unsigned int)cit->first.size());
    os.write( cit->first.c_str(), sizeof(char)* cit->first.size() );
    TWrite(os, cit->second->get_E());
    TWrite(os, cit->second->get_nu() );
  } // next cit

  // bbox
  bool flag = m_cmin.isValid() && m_cmax.isValid();
  TWrite<bool>(os, flag );

  if ( flag )
  {
    for (unsigned int ui=0; ui<n; ++ui)
      TWrite( os, m_cmin(ui) );
    for (unsigned int ui=0; ui<n; ++ui)
      TWrite( os, m_cmax(ui) );
  }

  // nodes
  TWrite(os, (unsigned int)m_vpNodes.size() );
  for ( typename NodeContainerType::const_iterator
        cit = m_vpNodes.begin();
        cit != m_vpNodes.end(); ++cit )
  {
    TWrite( os,  (int)(*cit)->get_id() );

    for (unsigned int ui=0; ui<n; ++ui)
      TWrite(os, (*cit)->coords()(ui));

    for (unsigned int ui=0; ui<n; ++ui)
      TWrite(os, (*cit)->delta()(ui) );

    for (unsigned int ui=0; ui<n; ++ui)
      TWrite(os, (bool)(*cit)->is_dof_active(ui) );
  } // next cit


  // elements
  // no items
  TWrite(os, (unsigned int)m_vpElements.size() );
  tNode* pnode = NULL;
  for ( ElementConstIterator cit = m_vpElements.begin();
        cit != m_vpElements.end(); ++cit )
  {
    TWrite( os, (*cit)->get_id() );

    TWrite( os, (unsigned int)(*cit)->material()->label().size() );
    os.write( (*cit)->material()->label().c_str(),
              sizeof(char) * (*cit)->material()->label().size() );

    for (unsigned int ui=0; ui<=n; ++ui)
    {
      if ( !(*cit)->get_node(ui, &pnode) )
        throw "TMesh::save error getting element node";
      TWrite(os, pnode->get_id() );
    } // next ui
  } // next cit
}

template<class Cstr, int n>
void
TMesh<Cstr,n>::set_constants(double _e,
                             double _nu)
{
  clear_materials();

  VMaterial* pmat = Cstr::material("default",_e,_nu);

  add_material("default",pmat);

  for ( typename std::vector<tElement*>::iterator it = m_vpElements.begin();
        it != m_vpElements.end();
        ++it )
    (*it)->set_material(pmat);
}

struct FunctorAddMaterial
{
  const VMaterial* cp_material;
  FunctorAddMaterial(const VMaterial* input)
      : cp_material(input)
  {}
  template<class T> T operator()(T v) const
  {
    v->set_material(cp_material);
  }
};

template<class Cstr, int n>
void
TMesh<Cstr,n>::set_default_material(VMaterial* pmaterial)
{
  clear_materials();

  if ( pmaterial->label().empty() ) pmaterial->set_label("default");

  this->add_material(pmaterial->label(), pmaterial);
  std::transform( m_vpElements.begin(), m_vpElements.end(),
                  m_vpElements.begin(),
                  FunctorAddMaterial(pmaterial)
                );
#if 0
  for ( typename std::vector<tElement*>::iterator it = m_vpElements.begin();
        it != m_vpElements.end();
        ++it )
    (*it)->set_material(pmaterial);
#endif
}

template<class Cstr,
int n>
void
TMesh<Cstr,n>::clear_materials()
{
  if ( m_mpMaterial.size() )
  {
    for ( MaterialIteratorType it = m_mpMaterial.begin();
          it != m_mpMaterial.end();
          ++it )
      delete it->second;
    m_mpMaterial.clear();
  }
}

template<class Cstr,
int n>
void
TMesh<Cstr,n>::add_material(std::string strLabel,
                            VMaterial* mat)
{
  if ( !strLabel.size() )
  {
    std::cerr << "TMesh::add_material -> void label for constants "
    << mat->get_E() << " " << mat->get_nu() << std::endl;
    return;
  }

  // check if the label exists already
  MaterialConstIteratorType mit = m_mpMaterial.find(strLabel);
  if ( mit != m_mpMaterial.end() ) // it's there
  {
    *(mit->second) = *mat; // by copy, the elements' address is not invalidated
  }
  else
    m_mpMaterial[strLabel] = mat;

  // make sure the label is also set in the material
  if ( mat->label() != strLabel )
    mat->set_label(strLabel);
}

template<class Cstr, int n>
void
TMesh<Cstr,n>::assign_material(unsigned i,
                               std::string strLabel)
{
  if ( i > m_vpElements.size()-1 )
  {
    std::cerr << "TMesh::assign_material -> index "
    << i << " beyond max val \n";
    return;
  }
  MaterialConstIteratorType  mit = m_mpMaterial.find(strLabel);
  if ( mit == m_mpMaterial.end() )
  {
    std::cerr << "TMesh::assign_material -> cannot find label "
    << strLabel << " in list of materials\n";
    return;
  }
  m_vpElements[i]->set_material(mit->second);
}

template<class Cstr, int n>
bool
TMesh<Cstr,n>::get_node(unsigned int i,
                        tNode** ppnode) const
{
  if ( i>=m_vpNodes.size() )
  {
    std::cerr << "TMesh::get_node -> incorrect index "
    << i << " >= " << m_vpNodes.size() << std::endl;
    *ppnode = NULL;
    return false;
  }
  *ppnode = m_vpNodes[i];
  return true;
}

#if 0
#undef __FUNCT__
#define __FUNCT__ "TMesh::closest_node"
template<class Cstr, int n>
TNode<n>*
TMesh<Cstr,n>::closest_node(tCoords c) const
{
  bool bdbg = false;
  bdbg = (c(0)==51)&&(c(1)==96)&&(c(2)==203);

  tIntCoords ic = m_hash.query_index(c);
  if (bdbg) std::cout << " cnode c = " << c << " index = " << ic << std::endl;

  // proceed in a fast-marching-like manner
  std::set<tIntCoords, cless<int,n> > set_visited;
  std::list< tIntCoords > lst_buffer;

  set_visited.insert(ic);
  lst_buffer.push_back(ic);

  double dist = 100.0;
  tNode* pargmin = NULL;
  tNode* pnode = NULL;
  do
  {
    for ( typename std::list<tIntCoords >::iterator it = lst_buffer.begin();
          it != lst_buffer.end();
          ++it )
    {
      double buf_dist = local_closest_node( m_hash.query_pos(*it), c, &pnode);
      if (pnode)
      {
        if (buf_dist<dist)
        {
          dist = buf_dist;
          pargmin = pnode;
        }
      }
    }
    if (!pargmin)
    {
      std::list<tIntCoords > lst_swap;
      for ( typename std::list<tIntCoords >::iterator it = lst_buffer.begin();
            it != lst_buffer.end();
            ++it )
      {
        std::stack<tIntCoords >* pstack = neighbors(*it);
        while ( !pstack->empty() )
        {
          bool bContinue = false;
          for (int alpha=0; alpha<n; ++alpha)
          {
            if ( pstack->top()(alpha) < 0 ||
                 pstack->top()(alpha) >= m_hash.ticks()(alpha) )
              bContinue = true;
          }
          if ( !bContinue )
          {
            typename std::set<tIntCoords >::const_iterator
              csit = set_visited.find( pstack->top() );
            if ( csit == set_visited.end() )
            {
              set_visited.insert( pstack->top() );
              lst_swap.push_back( pstack->top() );
            }
          }
          pstack->pop();
        }
        delete pstack;
      }
      lst_buffer = lst_swap;
    }
  }
  while ( !pargmin && !lst_buffer.empty() );

  return pargmin;
}

#undef __FUNCT__
#define __FUNCT__ "TMesh::element_at_point"
template<class Cstr, int n>
TElement<n>*
TMesh<Cstr,n>::element_at_point(const tCoords& c) const
{
  std::ostringstream os;

  bool bdbg=true;
  tNode* pnode = this->closest_node(c);
  if (bdbg)
  {
    if (pnode) os << "\t valid pnode\n";
    else os << "\t invalid pnode\n";
  }

  if ( pnode )
  {
    if (bdbg) os << " node id = " << pnode->get_id() << std::endl;
    // go through all the elements which contain the node
    // and check if any belongs
    typedef typename tNode::tElt_citer NodeEltsConstIterator;
    NodeEltsConstIterator cit, end;
    pnode->get_elt_citer(cit, end);

    for (; cit != end; ++cit )
    {
      if (bdbg) os << "\t elt = " << *cit << std::endl;
      if ( m_vpElements[ *cit ]->src_contains(c) )
        return m_vpElements[ *cit ];
    } // next cit
  }
  if (bdbg) std::cout << os.str() << std::endl;

  return NULL;
}


template<class Cstr,int n>
double
TMesh<Cstr,n>::local_closest_node(tStack nodes,
                                  tCoords c,
                                  tNode** ppnode) const throw(gmpErr)
{
  bool bdbg = false;
  bdbg = (c(0)==51)&&(c(1)==96)&&(c(2)==203);

  if (bdbg) std::cout << "\t\t stack size = " << nodes.size() << std::endl;
  tNode* pnode = NULL;
  tNode* pargmin = NULL;
  double dist = 10000.0;

  while ( !nodes.empty() )
  {
    pnode = m_vpNodes[ nodes.top() ];
    nodes.pop();
    if ( !pnode )
      throw gmpErr("TMesh::local_closest_node -> error from point ");
    double buf_dist = ( c - pnode->coords() ).norm();
    if (bdbg) std::cout << " node id = " << pnode->get_id()
      << " dist = " << buf_dist << std::endl;
    if ( buf_dist<dist)
    {
      dist = buf_dist;
      pargmin = pnode;
    }
  }
  *ppnode = pargmin;
  return dist;
}

#endif

template<class Cstr,int n>
void
TMesh<Cstr,n>::get_dst_box(tCoords& cmin,
                           tCoords& cmax) const
{
  typename std::vector<tElement*>::const_iterator cit
  = m_vpElements.begin();

  (*cit)->dst_box(cmin, cmax);
  tCoords elt_min, elt_max;
  for ( ++cit; cit != m_vpElements.end(); ++cit )
  {
    (*cit)->dst_box(elt_min, elt_max);
    cmin = min(elt_min, cmin);
    cmax = max(elt_max, cmax);
  }
}

template<class Cstr,int n>
unsigned int
TMesh<Cstr,n>::add_node(tNode* pnode)
{
  m_vpNodes.push_back(pnode);
  return m_vpNodes.size();
}

template<class Cstr, int n>
int
TMesh<Cstr,n>::add_elt(tElement* pelt)
{
  m_vpElements.push_back(pelt);
  return (int)m_vpElements.size() - 1;
}

#if 0
template<class Cstr, int n>
int
TMesh<Cstr,n>::build_index_src()
{
  m_hash.reset();
  std::cout << " build index src \n\t cmin = " << m_cmin
  << " cmax = " << m_cmax << std::endl;
  m_hash.init(m_cmin, m_cmax);

  for ( unsigned int i= unsigned int(0); i<m_vpNodes.size();
        ++i)
    m_hash.insert_pt( m_vpNodes[i]->coords(), i);
  return 1;
}
#endif

struct FunctorSourceVolumeAccumulator
{
  template<class In, class T> T operator()(T  init,
      In data)
  {
    return init+data->src_volume();
  }
};

struct FunctorDestinationVolumeAccumulator
{
  template<class In, class T> T operator()(T  init,
      In data)
  {
    return init+data->dst_volume();
  }
};

template<class Cstr,int n>
double
TMesh<Cstr,n>::src_volume() const
{
#if 0
  double dvol = 0.0;
  for ( typename std::vector<tElement*>::const_iterator cit
        = m_vpElements.begin();
        cit != m_vpElements.end();
        ++cit )
    dvol += (*cit)->src_volume();

  return dvol;
#endif
  return std::accumulate( m_vpElements.begin(),
                          m_vpElements.end(),
                          0.0,
                          FunctorSourceVolumeAccumulator()
                        );
}

template<class Cstr,int n>
double
TMesh<Cstr,n>::dst_volume() const
{
#if 0
  double dvol = 0.0;
  for (typename std::vector<tElement*>::const_iterator cit
       = m_vpElements.begin();
       cit != m_vpElements.end();
       ++cit )
    dvol += (*cit)->dst_volume();

  return dvol;
#endif
  return std::accumulate( m_vpElements.begin(),
                          m_vpElements.end(),
                          0.0,
                          FunctorDestinationVolumeAccumulator()
                        );
}

struct FunctorIndexGrabber
{
  template<class In, class Out>  Out operator()(In data) const
  {
    return data->get_id();
  }
};

struct FunctorClearElementTies
{
  template<class T>  T operator()(T v) const
  {
    v->clear_elts();
    return v;
  }
};

// use a mutable data member to avoid the overhead of
// constantly re-declaring this
template<class Cstr, int n>
struct FunctorUpdateNodesFromElements
{
  typedef TNode<n> tNode;
  mutable tNode* pnode;
  template<class T> T operator()(T v) const
  {
    for (unsigned int ui=0, nnodes(v->no_nodes());
         ui < nnodes; ++ui)
    {
      if ( !v->get_node(ui, &pnode) )
        throw std::logic_error("TMesh::remove_elts -> invalid node call");
      pnode->add_elt( v->get_id() );
    }
  }
};

template<class Cstr, int n>
template<class In>
void
TMesh<Cstr,n>::remove_elts(In begin, In end) /*throw(gmpErr)*/
{
  bool bdbg = false;
  if (bdbg) std::cout << " TMesh::remove_elts\n";

  // do a set operation
  //
  // copy the structure in a vector
  // sort it
  // create a container of the id's
  // sort it
  // perform a set difference
  // replace vector container with remaining ids
  //
  // for simplicity, I am assuming element ids
  //     also serve as random iterators
  //
  // after the set difference, update the ids of the elements...
  //     how to do this better?
  //
  // finally, update the nodes

  typedef std::vector<int> IntVector;
  IntVector vecIdx;
  std::copy( begin, end, std::back_inserter(vecIdx) );

  std::sort( vecIdx.begin(), vecIdx.end() );
  std::unique(vecIdx.begin(), vecIdx.end() );

  IntVector vecInitialMesh;
  vecInitialMesh.reserve( m_vpElements.size() );
  std::transform( m_vpElements.begin(),
                  m_vpElements.end(),
                  std::back_inserter(vecInitialMesh),
                  FunctorIndexGrabber()
                );

  IntVector vecDifference;
  std::set_difference( vecInitialMesh.begin(), vecInitialMesh.end(),
                       vecIdx.begin(), vecIdx.end(),
                       std::back_inserter( vecDifference )
                     );

  // finally re-create the container
  ElementContainerType newElementContainer;
  newElementContainer.reserve( vecDifference.size() );

  for ( ElementConstIterator cit = m_vpElements.begin();
        cit != m_vpElements.end();
        ++cit )
  {
    if ( std::binary_search( vecDifference.begin(),
                             vecDifference.end(),
                             (*cit)->get_id() ) )
      newElementContainer.push_back( *cit );
  } // next cit

  m_vpElements = newElementContainer;

#if 0
  // go through elements and replace element ids
  unsigned int count = 0;
  for ( ElementIterator it = m_vpElements.begin();
        it != m_vpElements.end(); ++it, ++count )
  {
    (*it)->set_id( count );
  } // next it

  std::transform( m_vpNodes.begin(), m_vpNodes.end(),
                  m_vpNodes.begin(),
                  FunctorClearElementTies() );
  std::transform( m_vpElements.begin(), m_vpElements.end(),
                  m_vpElements.begin(),
                  FunctorUpdateNodesFromElements() );

  lst_elts.sort();
  lst_elts.unique();

  In cit = begin;
  std::vector<tElement*> vbuf_elts;

  // assumption = the element id's are increasing
  // in the vector structure
  for ( typename std::vector<tElement*>::const_iterator vit =
          m_vpElements.begin();
        vit != m_vpElements.end();
        ++vit)
  {
    if ( (*vit)->get_id() != *cit )
      vbuf_elts.push_back( *vit );
    else
      ++cit;
  }
  m_vpElements = vbuf_elts;
  vbuf_elts.clear();

  // re-parse el <-> node relations
  for ( typename std::vector<tNode*>::iterator it = m_vpNodes.begin();
        it != m_vpNodes.end();
        ++it )
    (*it)->clear_elts();

  tNode* pnode = NULL;
  for ( unsigned int i = 0;
        i < m_vpElements.size();
        ++i )
  {
    tElement* pelt = m_vpElements[i];
    pelt->set_id(i);
    for ( unsigned int eid = 0; eid < (unsigned int)(pelt->no_nodes()); ++eid )
    {
      if ( !pelt->get_node(eid, &pnode ) )
        throw gmpErr("TMesh::remove_elts -> invalid node call ");
      pnode->add_elt(i);
    }
  }

  // go through the nodes and remove those which belong to no elts
  std::vector<tNode*> vbuf_pnodes;
  for ( unsigned int i = 0;
        i < m_vpNodes.size();
        ++i )
  {
    pnode = m_vpNodes[i];
    if ( pnode->get_nelts() )
    {
      vbuf_pnodes.push_back( pnode );
      pnode->set_id( vbuf_pnodes.size() -1 );
    }
    else delete pnode;
  }
  m_vpNodes = vbuf_pnodes;
#endif

}

template<int n>
bool
positive( const TCoords<int,n>& pt)
{
  for (int i=0; i<n; ++i)
    if ( pt(i) <= 0 ) return false;

  return true;

}

// proposed strategy for the direct image
//
// use the closest_node method -> it should do the trick most of the times
template<class Cstr, int n>
TCoords<double,n>
TMesh<Cstr,n>::dir_img(const tCoords& c_src,
                       bool signalTopology) const /*throw(gmpErr) */
{
  tCoords img;
  const tElement* pelt = this->element_at_point( c_src );

  if ( !pelt )
    img.status() = cOutOfBounds;
  else if ( pelt->orientation_pb() )
    img.invalidate();
  else
    img = pelt->dir_img( c_src );

  return img;
}

template<class Cstr, int n>
bool
TMesh<Cstr,n>::check_elt_id() const
{
  int index = 0;
  for (typename std::vector<tElement*>::const_iterator cit = m_vpElements.begin();
       cit != m_vpElements.end(); ++cit, ++index )
  {
    if ( index != (*cit)->get_id())
    {
      std::cerr << " id position mismatch -> " << index << " <> " <<
      (*cit)->get_id() << std::endl;
    }
  } // next cit

  // got through all the elts and check the back-cor
  typedef typename tNode::tElt_citer EltConstIterator;
  EltConstIterator cit_begin, cit_end;
  for (typename std::vector<tNode*>::const_iterator cit = m_vpNodes.begin();
       cit != m_vpNodes.end(); ++cit)
  {
    (*cit)->get_elt_citer( cit_begin, cit_end );
    while ( cit_begin != cit_end )
    {
      if ( m_vpElements[*cit_begin] == 0 )
      {
        std::cerr << " Invalid elt for node " << *(*cit) << std::endl;
        std::cerr << " No of elts = " << (*cit)->get_nelts() << std::endl;
        std::cerr << " elts = \n";
        while ( cit_begin != cit_end )
        {
          std::cerr << *cit_begin << "|";
          ++cit_begin;
        }
        std::cout << std::endl;
        exit(1);
      }
      ++cit_begin;
    }
  }
  return true;
}

// faster than the full test - stops at the first error found
template<class Cstr, int n>
bool
TMesh<Cstr,n>::orientation_pb(Frame f) const
{
  typename std::vector<tElement*>::const_iterator cit;
  for ( cit = m_vpElements.begin();
        cit != m_vpElements.end();
        ++cit )
    if ( (*cit)->orientation_pb(f) ) return true;

  return false;
}

template<class Cstr, int n>
template<class InsertIterator>
void
TMesh<Cstr,n>::orientation_pb(Frame f,
                              InsertIterator ii) const
{
  typename std::vector<tElement*>::const_iterator cit;
  for ( cit = m_vpElements.begin();
        cit != m_vpElements.end();
        ++cit )
    if ( (*cit)->orientation_pb(f) ) ii = (*cit)->get_id();

}

template<class Cstr, int n>
std::pair<TCoords<double,n>, TCoords<double,n> >
TMesh<Cstr,n>::invert()
{
  tCoords cmin, cmax;
  cmin.set(10.0);
  cmax.set(-10.0);

  for ( typename NodeContainerType::iterator it = m_vpNodes.begin();
        it != m_vpNodes.end(); ++it )
  {
    (*it)->invert();

    cmin = min( cmin, (*it)->coords() );
    cmax = max( cmax, (*it)->coords() );
  } // next it

  return std::make_pair(cmin, cmax);
}

template<class Cstr, int n>
int
TMesh<Cstr,n>::get_neighboring_elts(unsigned int i,
                                    int radius,
                                    ElementIndexContainer& eltIndex) const
{
  typedef std::set<unsigned int> IndexSetType;
  typedef std::stack<unsigned int> IndexStackType;
  IndexStackType nodeStack; // used to find new elements
  IndexStackType eltStack; // used to find new nodes
  IndexSetType nodeSet;
  IndexSetType eltSet;

  eltStack.push(i);
  eltSet.insert(i);

  IndexSetType::iterator it;

  tNode* pnode;
  for (int step = 0;
       step<radius && !eltStack.empty();
       ++step)
  {
    // populate the node stack and set
    while (!eltStack.empty())
    {
      unsigned int crt = eltStack.top();
      eltStack.pop();
      const tElement* cpelt = this->get_elt(crt);

      // insert the nodes in the set if they are new
      for (unsigned int ui(0), nnodes(cpelt->no_nodes());
           ui < nnodes; ++ui)
      {
        this->get_node(ui, &pnode);
        // see if this is a new node
        it = nodeSet.find( pnode->id() );
        if ( it == nodeSet.end() )
        {
          nodeStack.push( pnode->id() );
          nodeSet.insert( pnode->id() );
        }
      }
    }

    // use the node stack to add new elements
    while (!nodeStack.empty())
    {
      unsigned int crt = nodeStack.top();
      nodeStack.pop();
      this->get_node(crt, &pnode);

      typename tNode::tElt_citer begin, end;
      pnode->get_elt_citer(begin, end);
      for (; begin!=end; ++begin)
      {
        it = eltSet.find(*begin);
        if ( it == eltSet.end() )
        {
          eltStack.push(*begin);
          eltSet.insert(*begin);
        }
      }
    }

  } // next step

  eltIndex.clear();
  for ( IndexSetType::const_iterator cit = eltSet.begin();
        cit != eltSet.end(); ++cit)
    eltIndex.push_back( *cit );

  return eltIndex.size();
}

#endif

