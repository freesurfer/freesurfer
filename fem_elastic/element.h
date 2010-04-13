
#ifndef H_ELEMENT_H
#define H_ELEMENT_H

#include <iostream>
#include <vector>

#include "material.h"
#include "node.h"

//-----------
//
// Class definition
//

typedef enum
{
  src,
  dst,
  both
} Frame;


template<int n>
class TElement
{
public:
  typedef TCoords<double,n> tCoords;
  TElement();
  virtual ~TElement()
  {}

  void set_material(const VMaterial* cpMaterial)
  {
    m_cpMaterial=cpMaterial;
  }
  const VMaterial* material() const
  {
    return m_cpMaterial;
  }

  int   get_id() const
  {
    return m_id;
  }
  int   id() const
  {
    return m_id;
  }
  void  set_id(int id);

  bool get_node(int index, TNode<n>** ppNode) const;
  int  no_nodes() const
  {
    return (int)m_vpNodes.size();
  }

  void add_node(TNode<n>* pNode);
  virtual SmallMatrix get_matrix() const=0; // compute the stiffness matrix

  int no_dofs() const
  {
    return n*m_vpNodes.size();
  }

  //--
  virtual bool dst_contains(const tCoords& c) const =0;
  void dst_box(tCoords& c_min, tCoords& c_max) const;
  virtual tCoords inv_img(const tCoords& dst_coords) const =0;
  virtual double dst_volume() const=0; // returns the volume
  // element (=area in 2D)
  //--

  //---
  virtual bool src_contains(const tCoords& c) const =0;
  void src_box(tCoords& c_min, tCoords& c_max) const;
  virtual tCoords dir_img(const tCoords& src_coords) const=0;
  virtual double src_volume() const=0; //  volume element (=area in 2D)
  //---

  virtual bool orientation_pb(Frame f=both) const=0; // returns 1 if
  // the determinant changes signs
  virtual bool orientation_test(double dalpha) const=0;

  // a valid return value will be comprised between 0 and 1
  virtual double shape_fct(int node_id, const tCoords& pt) const = 0;

  virtual void print(std::ostream& os) const;

protected:
  int      m_id; // -1 = not set

  typedef typename std::vector<TNode<n>*> NodeContainerType;
  NodeContainerType m_vpNodes;  // no ownership over nodes ->
  //shared between elements
  const VMaterial*       m_cpMaterial; // material is shared
};

template<int n>
std::ostream& operator<<(std::ostream& os, const TElement<n>& elt);

//-------------------------------------------
//
// class implementation
//

template<int n>
TElement<n>::TElement()
    : m_id(-1)
{
  m_cpMaterial = NULL;
}

template<int n>
void
TElement<n>::set_id(int id)
{
  m_id = id;
}

template<int n>
void
TElement<n>::add_node(TNode<n>* pNode)
{
  m_vpNodes.push_back(pNode);
  pNode->add_elt(m_id);
}

template<int n>
bool
TElement<n>::get_node( int index,
                       TNode<n>** ppNode) const
{
  if ( index>=(int)m_vpNodes.size() )
  {
    std::cerr << " TElement::get_node -> improper index value = "
    << index << std::endl;
    *ppNode = NULL;
    return false;
  }
  *ppNode = m_vpNodes[index];
  return true;
}

template<int n>
void
TElement<n>::dst_box(tCoords& c_min,
                     tCoords& c_max) const
{
  if ( m_vpNodes.empty() )
  {
    c_min = tCoords(-1.0);
    c_max = tCoords(-1.0);
    return;
  }

  typename std::vector<TNode<n>*>::const_iterator cit = m_vpNodes.begin();
  c_min = (*cit)->dst_coords();
  c_max = c_min;

  tCoords tc_buf;
  for ( ++cit; cit != m_vpNodes.end(); ++cit )
  {
    tc_buf = (*cit)->dst_coords();
    c_min = min( c_min, tc_buf );
    c_max = max( c_max, tc_buf );
  }

}

template<int n>
void
TElement<n>::src_box(tCoords& c_min,
                     tCoords& c_max) const
{
  if ( m_vpNodes.empty() )
  {
    c_min = tCoords(-1.0);
    c_max = tCoords(-1.0);
    return;
  }
  typename std::vector<TNode<n>*>::const_iterator cit = m_vpNodes.begin();
  c_min = (*cit)->coords();
  c_max = c_min;

  tCoords tc_buf;
  for ( ++cit;
        cit != m_vpNodes.end();
        ++cit )
  {
    tc_buf = (*cit)->coords();
    c_min = min( c_min, tc_buf);
    c_max = max( c_max, tc_buf);
  }
}

template<int n>
void
TElement<n>::print(std::ostream& os) const
{
  os << " elt id = " << m_id << std::endl;
  os << " nodes \n";
  for ( typename NodeContainerType::const_iterator cit = m_vpNodes.begin();
        cit != m_vpNodes.end(); ++cit )
    os << **cit << std::endl;
  // todo add material
}

template<int n>
std::ostream& operator<<(std::ostream& os,
                         const TElement<n>& elt)
{
  elt.print(os);
  return os;
}


#endif // H_ELEMENT_H
