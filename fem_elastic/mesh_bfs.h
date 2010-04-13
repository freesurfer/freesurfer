/*

Gheorghe Postelnicu, 2007

Mesh BFS - will do a breadth-first-search around a given seed

*/

#ifndef _h_mesh_bfs_h
#define _h_mesh_bfs_h

#include <functional>
#include <queue>
#include <set>
#include <vector>

#include "mesh.h"

template<class TMesh>
class TMeshBfs
{
public:
  typedef TMesh MeshType;

  typedef typename MeshType::tNode      NodeType;
  typedef typename MeshType::tElement   ElementType;
  typedef typename MeshType::tCoords    CoordsType;

  typedef std::vector<unsigned int>     IndexContainerType;

  TMeshBfs(const MeshType* pmesh) : m_pmesh(pmesh)
  {}

  // will go through the element nodes and stop at a given radius
  //
  // why constrain the output data type
  // just templatize this over an std::output_iterator
  //
  // what is the best way of constraining it?
  // 1. set S (template argument) = std::output_iterator<Container,e.g. list>
  // 2. set S (template argument) = Container and
  //force std::output_iterator<Container>
  //
  // will choose 1. for now
  template<class OutputIterator>
  void Visit(unsigned int seed,
             unsigned int radius,
             OutputIterator oit);

protected:
  const MeshType* m_pmesh;
};

template<class T, unsigned int n>
class CircularArray
{
public:
  CircularArray() : m_counter(0)
  {}

  inline T& current()
  {
    return m_data[ m_counter % n ];
  }
  inline const T& current() const
  {
    return m_data[ m_counter % n ];
  }

  inline T& getIncrement(int i)
  {
    return m_data[ (m_counter+i) % n ];
  }
  inline const T& getIncrement(int i) const
  {
    return m_data[ (m_counter+i) % n ];
  }
  void shift(int value=0)
  {
    m_counter += value;
  }
protected:
  int m_counter;
  T           m_data[n];

};

template<class T>
class BiCircularArray : public CircularArray<T,2>
{
public:
  BiCircularArray() : CircularArray<T,2>()
  {}
  void next()
  {
    if (this->m_counter) this->shift(-1);
    else this->shift(1);
  }
  inline T& other()
  {
    return this->getIncrement(1);
  }
  inline const T& other() const
  {
    return this->getIncrement(1);
  }
};

/*
  the seed and the indices refer to ELEMENTS.

  The return will be channeled through the output iterator
*/
template<class TMesh>
template<class OutputIterator>
void
TMeshBfs<TMesh>::Visit(unsigned int seed,
                       unsigned int radius,
                       OutputIterator oit)
{
  typedef std::set<unsigned int>   IndexSetType;
  typedef std::vector<unsigned int> IndexVectorType;
  typedef std::queue<unsigned int> IndexQueueType;

  typename NodeType::tElt_citer eltBegin, eltEnd;
  IndexSetType nodeSet; // contains the nodes which have already been visited
  IndexSetType elementSet;

  BiCircularArray<IndexQueueType> queueSwitch;

  const ElementType* cpElt;
  NodeType*    pnode;

  if ( !m_pmesh ) throw std::string (" NULL mesh in TMeshBfs");

  queueSwitch.current().push(seed);
  int eltIndex;

  // the radius refers to the ELEMENT graph
  //
  while ( radius-- && !queueSwitch.current().empty() )
  {
    IndexSetType eltWave;
    std::insert_iterator<IndexSetType> iiSet(eltWave, eltWave.begin());

    while (!queueSwitch.current().empty() )
    {
      // if element has already been visited, continue
      eltIndex = queueSwitch.current().front();
      queueSwitch.current().pop();

      if ( elementSet.find( eltIndex ) != elementSet.end() ) continue;
      elementSet.insert( eltIndex );
      cpElt = m_pmesh->get_elt( eltIndex );
      if ( !cpElt ) throw std::string( " NULL element in TMeshBfs::Visit");

      // get nodes belonging to the element
      for (unsigned int ui(0), nnodes(cpElt->no_nodes());
           ui < nnodes; ++ui)
      {
        if ( !cpElt->get_node(ui, &pnode) )
          throw std::string( " Failed to get node in TMeshBfs::Visit");

        // continue if node was already visited
        // this speeds up by preventing backward propagation
        if ( nodeSet.find( pnode->get_id() ) != nodeSet.end() ) continue;

        nodeSet.insert( pnode->get_id() );

        // get all elements that contain node
        pnode->get_elt_citer(eltBegin, eltEnd);
        std::copy( eltBegin, eltEnd, iiSet );
      } // next ui
    } // end while queueSwitch

    // the next queue will be made up from the set difference
    // between the nodes already visited and the current candidates
    IndexVectorType vec;
    std::insert_iterator<IndexVectorType> iiVec(vec, vec.begin());
    std::set_difference( eltWave.begin(), eltWave.end(),
                         elementSet.begin(), elementSet.end(),
                         iiVec );
    // the other queue is always empty
    for ( typename IndexVectorType::const_iterator cit = vec.begin();
          cit != vec.end(); ++cit )
    {
      queueSwitch.other().push(*cit);
    } // next cit
    queueSwitch.next();

  } // end while radius

  std::copy( elementSet.begin(), elementSet.end(), oit );
}

#endif
