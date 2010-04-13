
#ifndef H_NODE_H
#define H_NODE_H

#include <set>

#include "coords.h"

template<int n>
class TNode
{
public:
  TNode();
  TNode(const TNode& tn);
  typedef TCoords<double, n> tCoords;

  inline void set_coords(tCoords c)
  {
    m_coords = c;
  }
  inline void set_id(int id)
  {
    m_id = id;
  }
  void set_dof_val(size_t i, double val);
  inline void set_delta(const tCoords& d)
  {
    m_delta = d;
  }
  void set_active(size_t i, bool flag);

  void invert();

  //---
  // interface to access elements which contain this node
  void add_elt(int elt_id);
  void clear_elts();
  void remove_elt(int elt_id);
  // get const iterators to the indices of the elements that contain this node
  typedef typename std::set<unsigned int>::const_iterator tElt_citer;
  void get_elt_citer( tElt_citer& cit_begin, tElt_citer& cit_end) const;
  int  get_nelts() const;

  int     get_id() const
  {
    return m_id;
  }
  int     id() const
  {
    return m_id;
  }
  int     get_no_dofs() const;
  bool    is_dof_active(int no) const
  {
    return m_bdof_active[no];
  }
  const bool* is_active() const
  {
    return &m_bdof_active[0];
  }
  double  get_dof(int no) const
  {
    return m_delta(no);
  }
  double& dof(int no);

  const tCoords& coords() const
  {
    return m_coords;
  }
  const tCoords& delta() const
  {
    return m_delta;
  }
  void setDst(const tCoords& dst)
  {
    m_delta = dst - m_coords;
  }
  tCoords dst_coords() const
  {
    return m_coords+m_delta;
  }

  void set_bc(const tCoords& c); // activates constraints on all directions
  void set_bc(int no, double val); // activates a constraint in one direction -> DISPLACEMENT

  void print(std::ostream& os) const;

protected:

private:
  int m_id;
  tCoords m_coords;
  tCoords m_delta;
  bool    m_bdof_active[n]; // TRUE if value prescribed

  std::set<unsigned int> m_elts; // set of elements which contain the node

  void clone(const TNode& tn);
};

template<int n>
std::ostream& operator<<(std::ostream& os, const TNode<n>& node);

//------------------------------------------------------
//
// Definitions
//

template<int n>
TNode<n>::TNode()
    : m_coords(),
    m_delta(),
    m_elts()
{
  std::fill_n(m_bdof_active, n, false);
  m_id = -1; // not set
}

template<int n>
TNode<n>::TNode(const TNode<n>& tn)
    : m_coords(),
    m_delta()
{
  clone();
}

template<int n>
void
TNode<n>::clone(const TNode<n>& tn)
{
  m_coords = tn.m_coords;
  m_delta  = tn.m_delta;
  m_id     = tn.m_id;
  std::copy( (tn.m_bdof_active), (tn.m_bdof_active)+n,
             m_bdof_active );
  m_elts   = tn.m_elts;
}

template<int n>
void
TNode<n>::set_dof_val(size_t i,
                      double val)
{
  m_delta(i) = val;
}

template<int n>
void
TNode<n>::set_active(size_t i,
                     bool flag)
{
  m_bdof_active[i] = flag;
}

template<int n>
void
TNode<n>::add_elt(int elt_id)
{
  m_elts.insert(elt_id);
}

template<int n>
void
TNode<n>::remove_elt(int elt_id)
{
  m_elts.erase(elt_id);
}

template<int n>
void
TNode<n>::clear_elts()
{
  m_elts.clear();
}

template<int n>
void
TNode<n>::get_elt_citer( tElt_citer& cit_begin,
                         tElt_citer& cit_end) const
{
  cit_begin = m_elts.begin();
  cit_end   = m_elts.end();
}

template<int n>
int
TNode<n>::get_nelts() const
{
  return (int)m_elts.size();
}

template<int n>
int
TNode<n>::get_no_dofs() const
{
  return n;
}

template<int n>
void
TNode<n>::set_bc(const tCoords& c)
{
  m_delta = c;
  std::fill_n(m_bdof_active, n, true);
}

template<int n>
void
TNode<n>::set_bc(int no,
                 double val)
{
  if ( no<0 || no>n )
  {
    std::cerr << " TNode::set_bc -> trying to access invalid entry "
    << no << std::endl;
    return;
  }
#if 0
  if ( m_bdof_active[no] )
  {
    std::cerr << " TNode::set_bc WARNING!!! -> BC already set for "
    << no << "\n\t imposing value " << val
    << " over " << m_delta[no] << std::endl;
  }
#endif
  m_bdof_active[no] = true;
  m_delta[no] = val;
}

template<int n>
void
TNode<n>::print(std::ostream& os) const
{
  os << " node [[ id=" << m_id << " coords=" << m_coords << " delta=" << m_delta << " ]]";
}

template<int n>
std::ostream& operator<<(std::ostream& os,
                         const TNode<n>& node)
{
  node.print(os);
  return os;
}

template<int n>
void
TNode<n>::invert()
{
  m_coords += m_delta;
  m_delta  = m_delta * -1.0;
}

#endif
