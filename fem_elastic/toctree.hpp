
#ifndef _h_own_octree_hpp
#define _h_own_octree_hpp

#include <vector>

#include "stl_utils.hpp"
#include "coords.h"

#define POW2(N) 1<<N

namespace toct
{

template<unsigned int N>
class TCubicRegion
{
public:
  typedef TCoords<double,N> CoordsType;
  TCubicRegion() : m_cmin(), m_cmax()
  {}
  TCubicRegion(const CoordsType& cmin, const CoordsType& cmax) : m_cmin(cmin), m_cmax(cmax)
  {}

  bool isInside(const CoordsType& pt) const
  {
    return ( (m_cmin<=pt) && (pt<=m_cmax) );
  }

  /*

  decides if 2 cubic regions are overlapping or not
  no informative priors are used, so there are a lot of comparisons necessary for this

  */
  bool isOverlapping(const TCubicRegion& other) const
  {
    // one of the points of the "other" needs to be inside "this" one
    CoordsType pt;
    static const unsigned int npoints =  1<<N;

    for (unsigned int ptIdx=0; ptIdx<npoints; ++ptIdx)
    {
      for (unsigned int ui=0; ui<N; ++ui)
      {
        pt(ui) = (ptIdx & (1<<ui))?other.m_cmin(ui):other.m_cmax(ui);
      } // next ui
      if ( this->isInside(pt) ) return true;
    } // next ptIdx

    return false;
  }

  const CoordsType& min() const
  {
    return m_cmin;
  }
  const CoordsType& max() const
  {
    return m_cmax;
  }

  CoordsType& min()
  {
    return m_cmin;
  }
  CoordsType& max()
  {
    return m_cmax;
  }

protected:
  CoordsType m_cmin, m_cmax;
};

struct OctreeData
{
  unsigned int maxLevels;
  unsigned int maxElementsInLeaf;
  unsigned int currentLevel;
  OctreeData() : maxLevels(), maxElementsInLeaf(), currentLevel()
  {}
  // provide constructor specifically for the most current operation
  OctreeData( const OctreeData& other, bool incrementLevel =false )
      : maxLevels(other.maxLevels), maxElementsInLeaf(other.maxElementsInLeaf),
      currentLevel(other.currentLevel)
  {
    if (incrementLevel) ++this->currentLevel;
  }
};

/*---------------------------------------------------------------------------

Generic Node class - this is inspired by a Composite pattern (see GoF)

However, in this case it seemed better to separate the leaves and the intermediate nodes

This class is supposed to support both dimensions 2 and 3 by default

*/
template<class ElementProxy, unsigned int N>
class TNode : public TCubicRegion<N>
{
public:
  typedef typename std::array<unsigned char, POW2(N)> DeterminantType;
protected:
  /*
    since the underlying type is unsigned char,
    this will only work in dimension 2, 3 or 4

    !!! being static, for this to link, it needs to be initialized
  */
  static const DeterminantType& GetDeterminant()
  {
    static DeterminantType s_octantDeterminant
    = InitializeOctantDeterminant();

    return s_octantDeterminant;
  }
public:
  typedef TCubicRegion<N> Superclass;
  typedef typename Superclass::CoordsType CoordsType;

  TNode(const CoordsType& cmin=CoordsType(),
        const CoordsType& cmax=CoordsType() )
      : TCubicRegion<N>(cmin,cmax)
  {}
  virtual ~TNode()
  {}

  virtual bool insertItem(const ElementProxy* ep,
                          TNode*& pCaller,
                          const OctreeData& data) = 0;
  virtual const ElementProxy* element_at_point
  (const CoordsType& pt) const = 0;
  virtual unsigned int getElementCount() const = 0;


private:
  // rule => abc --> 1<<(4+a) + 1<<(2+b) + 1<<c
  //
  // USE 2 bits per dimension
  static DeterminantType InitializeOctantDeterminant()
  {
    DeterminantType det;
    std::cout << " InitializeOctantDeterminant\n";
    unsigned int p2, res;
    for (unsigned int idx=0; idx<POW2(N); ++idx)
    {
      p2 = 1;
      res = 0;
      for (unsigned int ui=0; ui<N; ++ui, p2<<=1)
      {
        res += 1<<( 2*ui + ( (p2&idx) ? 1: 0 ) );
      }
      det[idx] = res;
    }
    return det;
  }

};

template<class ElementProxy, unsigned int N> class TerminalNode;

/*-----------------------------------------------------------------------------

Specific implementation of intermediate nodes. These will not contain elements themselves,
only route requests toward specific branches.

Octree implementation => each node will hold exactly 8 branches.
The separation happens according to the orthogonal directions - the most likely point to do the separation is the center.

The most time-critical operation occurs when some element data is inserted. Hence, it makes sense to optimize for this.

In this case, a natural idea is to minimize the number of comparisons necessary for the insertion. Namely, IF
at the moment of the insertion it has already been checked that the element overlaps with this branch, then all that is
actually needed is to decide which of the octants are overlapping with the element....

*/

// rule => abc --> 1<<(4+a) + 1<<(2+b) + 1<<c

template<class ElementProxy, unsigned int N>
class TIntermediateNode : public TNode<ElementProxy,N>
{
  typedef TNode<ElementProxy,N> Superclass;

  // I know that each intermediate node will have exactly 8 children
  //       so enforce it
  typedef std::array<Superclass*,POW2(N)> NodeContainer;

  // the next enum is used for bitwise flag operations on each enumeration when deciding
  //     which octant is overlapping with the current region.
  enum OctantOverlapPerDimension
  {
    Lower = 1,
    Upper,    // will automatically be 2
    Both      // will automatically be 3
  };

public:
  typedef typename Superclass::CoordsType CoordsType;

  // define the "kids" here
  TIntermediateNode(const CoordsType& cmin=CoordsType(),
                    const CoordsType& cmax=CoordsType() ) : Superclass(cmin,cmax),
      m_cmidPoint(),
      m_items()
  {
    typedef TerminalNode<ElementProxy,N> TerminalNode;

    m_cmidPoint = .5 * ( this->m_cmin + this->m_cmax );

    // create the new terminal nodes
    // --> implicitly, decide on the ordering - essential in the following to optimize for speed

    CoordsType bufMin, bufMax;
    unsigned int p2;
    for (unsigned int ptIndex =0; ptIndex < POW2(N); ++ptIndex)
    {
      bufMin = this->m_cmin;
      bufMax = this->m_cmidPoint;
      p2 = 1;
      // this update rule is tightly connected with the numbering of the octants
      // if n = abc (base 2), then a~x, b~y, c~z
      for (int ui=N-1; ui>=0; --ui, p2<<=1)
      {
        if ( ptIndex & p2 )
        {
          bufMin(ui) = this->m_cmidPoint(ui);
          bufMax(ui) = this->m_cmax(ui);
        }
      } // next ui
      m_items[ptIndex] = new TerminalNode( bufMin, bufMax );
    } // next ptIndex

  }
  ~TIntermediateNode()
  {
    std::transform( m_items.begin(), m_items.end(),
                    m_items.begin(),
                    FunctorDeletePointer()
                  );
  }

  bool insertItem(const ElementProxy* ep,
                  Superclass*& pCaller,
                  const OctreeData& data)
  {
    unsigned char overlappingOctants = this->FindOverlappingOctants(*ep);
    OctreeData newData(data,true);
    bool inserted = false;
    for (unsigned int ui=0; ui< POW2(N); ++ui)
      if ( (Superclass::GetDeterminant()[ui] & overlappingOctants)
           == Superclass::GetDeterminant()[ui] )
      {
        inserted = true;
        m_items[ui]->insertItem(ep, m_items[ui], newData);
      }
    if ( !inserted )
      std::cerr << " element NOT inserted - overlappingOctants = " << (int)overlappingOctants << std::endl;


    return inserted;
  }

  /*
    apply the same strategy as for the insertion, namely minimize the number of comparisons
    needed to decide which is the appropriate bucket

    hence assume a test has already determined the point fits in this branch
    the only thing remaining then is to decide which branch to continue on

    use the fact that branches are encoded in the array as abc (binary decomposition),
    where a->x , b->y and c->z
  */
  const ElementProxy* element_at_point(const CoordsType& pt) const
  {

    unsigned int branch = 0, exponent;
    for (unsigned int ui=0; ui<N; ++ui)
    {
      exponent = N-1-ui;
      if ( pt(ui) > m_cmidPoint(ui) )
        branch += (1<<exponent);
    }

    return m_items[branch]->element_at_point(pt);

  }

  unsigned int getElementCount() const
  {
    unsigned int count = 0;
    for ( typename NodeContainer::const_iterator cit = m_items.begin();
          cit != m_items.end(); ++cit )
      count += (*cit)->getElementCount();
    return count;
  }

private:
  // the separation into octants is based on the comparison with the mid point.
  CoordsType  m_cmidPoint;
  NodeContainer m_items;

  /*
    when this function is called, IT IS KNOWN there is an overlap with THIS region
    the purpose is to decide which of the sub-regions overlaps with the argument
  */
  unsigned char FindOverlappingOctants(const TCubicRegion<N>& region)
  {
    const double deps = 1.0e-5;
    unsigned char ret = 0;
    int exponent;
    for (unsigned int ui=0; ui<N; ++ui)
    {
      exponent = 2*(N-1-ui);
      if ( region.min()(ui) < this->m_cmidPoint(ui) + deps )
      {
        if ( region.max()(ui) > this->m_cmidPoint(ui) - deps )
          ret += ( Both << exponent );
        else
          ret += ( Lower << exponent );
      }
      else
        ret += ( Upper << exponent );
    } // next ui
    return ret;
  }
};

/*--------------------------------------------------

*/


template<class ElementProxy, unsigned int N>
class TerminalNode : public TNode<ElementProxy,N>
{
  typedef TNode<ElementProxy,N> Superclass;
  typedef std::vector<const ElementProxy*> ProxyContainer;
public:
  typedef typename Superclass::CoordsType CoordsType;
  TerminalNode(const CoordsType& cmin,
               const CoordsType& cmax) : Superclass(cmin, cmax),
      m_elements()
  {}

  bool insertItem(const ElementProxy* ep,
                  Superclass*& pCaller,
                  const OctreeData& data)
  {
    // at this point, already know it is overlapping

    this->m_elements.push_back(ep);

    // decide if element goes in the current pool,
    //     or a split occurs
    if ( data.currentLevel < data.maxLevels &&
         data.maxElementsInLeaf < m_elements.size() )
    {
      // do a split

      // create an intermediate node with the same geometry
      //
      // when the new node is created, it readily creates terminal nodes for coming elements
      typedef TIntermediateNode<ElementProxy,N> IntermediateNodeType;
      IntermediateNodeType* pnewNode =
        new IntermediateNodeType(this->m_cmin,
                                 this->m_cmax);

      OctreeData newData(data, true);

      Superclass* pnodeReplacement = static_cast<Superclass*>(pnewNode);
      // after creating the new nodes, re-insert all the current element proxies in the new structure
      for ( typename ProxyContainer::iterator cit = this->m_elements.begin();
            cit != this->m_elements.end(); ++cit )
      {
        pnodeReplacement->insertItem( &(*(*cit)), pnodeReplacement, newData );
      } // next cit

      // finally, replace the current node in the parent
      delete pCaller;
      pCaller = pnodeReplacement;
    }

    return true;
  }

  const ElementProxy* element_at_point(const CoordsType& pt) const
  {
    for ( typename ProxyContainer::const_iterator cit = m_elements.begin();
          cit != m_elements.end(); ++cit )
      if ( (*cit)->contains(pt) ) return *cit; //else ++g_missCount;

    return NULL;
  }

  unsigned int getElementCount() const
  {
    return static_cast<unsigned int>( m_elements.size() );
  }

private:
  ProxyContainer m_elements;
};

template<class ElementProxy, unsigned int N>
class Octree
{
  typedef TNode<ElementProxy,N> NodeType;
public:
  typedef TCoords<double,N> CoordsType;
  Octree( const CoordsType& cmin,
          const CoordsType& cmax,
          unsigned int maxLevels,
          unsigned int maxElementsInLeaf )
      : m_octreeData()
  {
    typedef TerminalNode<ElementProxy,N>
    TerminalNodeType;
    m_pnode = new TerminalNodeType(cmin,
                                   cmax);
    m_octreeData.maxLevels = maxLevels;
    m_octreeData.maxElementsInLeaf = maxElementsInLeaf;
    m_octreeData.currentLevel = 0;
  }

  ~Octree()
  {
    delete m_pnode;
  }

  bool insertItem(const ElementProxy* ep)
  {
    // check if overlapping
    if ( !this->m_pnode->isOverlapping( *ep ) ) return false;

    return m_pnode->insertItem(ep, m_pnode, m_octreeData);
  }

  const ElementProxy* element_at_point(const CoordsType& pt) const
  {
    // perform initial test on the whole grid bounding box
    if ( !m_pnode->isInside(pt) ) return NULL;

    return m_pnode->element_at_point(pt);
  }

  unsigned int getElementCount() const
  {
    return m_pnode->getElementCount();
  }

private:
  NodeType* m_pnode;
  OctreeData m_octreeData;
};

};

#endif
