/*

Gheorghe Postelnicu, 2006

Problems Mesh Crop Filter - crops a certain radius around problem elements

DFS with given radius - the radius is reinitialized at each problem encountered

*/

#ifndef _h_pbmesh_crop_h
#define _h_pbmesh_crop_h

#include <functional>
#include <map>
#include <set>
#include <stack>

#include "mesh.h"

/*****************************************

        Class declaration

 boundaryNodes - indices of the boundary nodes in the NEW mesh
 nodeMap - map NEW -> OLD
 elementMap - map NEW -> OLD

******************************************/

template<class Cstr, int n>
class TPbMeshCrop
{
public:
  typedef TMesh<Cstr,n> MeshType;

  typedef typename MeshType::tNode NodeType;
  typedef typename MeshType::tElement ElementType;
  typedef typename MeshType::tCoords CoordsType;
  typedef std::set<unsigned int> IndexSetType;

  typedef std::map<unsigned int, unsigned int> MapType;

  TPbMeshCrop(const MeshType* pmesh);

  void Crop(MeshType* outMesh,
            MapType&  nodeMap,
            MapType&  elementMap,
            IndexSetType& boundaryNodes,
            unsigned int seed,
            unsigned int radius);

private:
  TPbMeshCrop(const TPbMeshCrop&); // purposely not implemented
  void operator=(const TPbMeshCrop&); // purposely not implemented

  IndexSetType m_leafElements;

  const MeshType* m_pmesh;

  void bfs(unsigned int seed,
           unsigned int radius,
           unsigned int maxRadius,
           IndexSetType& eltVisited,
           IndexSetType& nodeVisited);

  std::pair<CoordsType,CoordsType> ComputeBoundingBox(MeshType* pmesh);
};

/**********************************************

         Class implementation

***********************************************/

template<class Cstr, int n>
TPbMeshCrop<Cstr,n>::TPbMeshCrop(const MeshType* pmesh)
{
  m_pmesh = pmesh;
}

template<class Cstr, int n>
void
TPbMeshCrop<Cstr,n>::Crop(MeshType* outMesh,
                          MapType& nodeMap,
                          MapType& elementMap,
                          IndexSetType& boundaryNodes,
                          unsigned int seed,
                          unsigned int radius)
{
  m_leafElements.clear();

  // get elements
  IndexSetType eltVisited; // will contain all the elements
  // used to populate the mesh
  IndexSetType nodeVisited; // will contain all
  // INTERNAL nodes -> used to determine boundary nodes
  this->bfs( seed,
             radius,
             radius,
             eltVisited,
             nodeVisited );

  // get all the nodes for the elements in the set
  IndexSetType allNodes;
  for ( typename IndexSetType::const_iterator cit = eltVisited.begin();
        cit != eltVisited.end();
        ++cit )
  {
    NodeType* pnode = NULL;
    const ElementType* cpelt = m_pmesh->get_elt(*cit);
    for (unsigned int ui(0), nnodes(cpelt->no_nodes());
         ui < nnodes; ++ui)
    {
      cpelt->get_node(ui, &pnode);
      allNodes.insert( pnode->get_id() );
    } // next ui
  } // next cit

  // populate nodes
  NodeType* pnodeIn;
  NodeType* pnodeOut;
  unsigned int index = 0;
  MapType internalMap;
  for ( typename IndexSetType::const_iterator cit = allNodes.begin();
        cit != allNodes.end();
        ++cit, ++index )
  {
    m_pmesh->get_node(*cit, &pnodeIn);

    pnodeOut = Cstr::node(index,
                          pnodeIn->coords(),
                          pnodeIn->delta(),
                          pnodeIn->is_active() );
    outMesh->add_node(pnodeOut);
    internalMap[ *cit ] = index;
    nodeMap[ index ] = *cit;
  } // next cit, index

  // setup elements
  index = 0;
  typename MapType::const_iterator mapCiter;
  const ElementType* cpelt;
  std::vector<NodeType*> vpnodes;
  for ( typename IndexSetType::const_iterator cit = eltVisited.begin();
        cit != eltVisited.end(); ++cit, ++index )
  {
    vpnodes.clear();
    cpelt = m_pmesh->get_elt(*cit);
    for (unsigned int ui(0), nnodes(cpelt->no_nodes());
         ui < nnodes; ++ui)
    {
      cpelt->get_node(ui, &pnodeIn);
      outMesh->get_node( internalMap[pnodeIn->id()], &pnodeOut );
      vpnodes.push_back( pnodeOut );
    }
    ElementType* pelt = Cstr::elt(index,
                                  NULL,
                                  vpnodes );
    outMesh->add_elt( pelt );
    elementMap[ index ] = cpelt->id() ;
  } // next cit, index

  // compute the bounding box
  std::pair<CoordsType, CoordsType> bbox;
  bbox = this->ComputeBoundingBox( outMesh );

  outMesh->m_cmin = bbox.first;
  outMesh->m_cmax = bbox.second;

  // fill the boundary nodes structure
  typedef std::vector<unsigned int> VectorType;
  VectorType boundaryInitial;
  std::insert_iterator<VectorType> ii(boundaryInitial, boundaryInitial.begin());
  std::set_difference( allNodes.begin(), allNodes.end(),
                       nodeVisited.begin(), nodeVisited.end(),
                       ii );

  // convert the node indices to the new mesh
  for ( typename VectorType::const_iterator cit = boundaryInitial.begin();
        cit != boundaryInitial.end(); ++cit)
    boundaryNodes.insert( internalMap[*cit] );
}

/*

2 tetrahedra are considered neighbors if they share a node

*/
template<class Cstr, int n>
void
TPbMeshCrop<Cstr,n>::bfs(unsigned int seed,
                         unsigned int radius,
                         unsigned int maxRadius,
                         IndexSetType& eltVisited,
                         IndexSetType& nodeVisited)
{
  IndexSetType tetSet; // used to populate the neighbor list
  NodeType* pnode = NULL;
  const ElementType* cpseed = NULL;
  //unused: const ElementType* cpelt = NULL;

  typename IndexSetType::const_iterator cit;
  typename NodeType::tElt_citer eltBegin, eltEnd;

  // if element was already visited, return
  if ( eltVisited.find(seed) != eltVisited.end() )
  {
    // if not a leaf, truly visited
    if ( m_leafElements.find(seed) == m_leafElements.end() )
      return;
  }

  cpseed = m_pmesh->get_elt(seed);
  if ( !cpseed  )
  {
    std::cerr << " Failed to extract element with seed " << seed << std::endl;
    exit(1);
  }

  // criterium for seed-reinit
  if ( cpseed->orientation_pb( dst ) )
    radius = maxRadius;

  eltVisited.insert( seed );

  if ( !radius )
  {
    m_leafElements.insert( seed );
    return;
  }
  else
  {
    // if previously a leaf, erase it
    cit = m_leafElements.find( seed );
    if ( cit != m_leafElements.end() )
      m_leafElements.erase( cit );
  }
  radius -= 1;

  // populate the set of neighbors
  for (unsigned int ui(0), nnodes(cpseed->no_nodes());
       ui < nnodes; ++ui)
  {
    if ( !cpseed->get_node(ui, &pnode))
    {
      std::cerr << " Failed to get node for element " << seed << std::endl;
      exit(1);
    }
    // continue if node was already visited
    // no backward propagation
    if ( nodeVisited.find( pnode->get_id() ) != nodeVisited.end() ) continue;

    nodeVisited.insert( pnode->get_id() );

    for ( pnode->get_elt_citer(eltBegin, eltEnd);
          eltBegin != eltEnd;
          ++eltBegin )
      tetSet.insert( *eltBegin );
  }

  for ( typename IndexSetType::const_iterator tetCiter = tetSet.begin();
        tetCiter != tetSet.end(); ++tetCiter)
    this->bfs( *tetCiter,
               radius,
               maxRadius,
               eltVisited,
               nodeVisited );

}


template<class Cstr, int n>
std::pair<typename TPbMeshCrop<Cstr,n>::CoordsType,
typename TPbMeshCrop<Cstr,n>::CoordsType>
TPbMeshCrop<Cstr,n>::ComputeBoundingBox(MeshType* pmesh)
{
  typename MeshType::NodeConstIterator nodeCiter, nodeEnd;

  pmesh->get_node_citer(nodeCiter, nodeEnd);
  CoordsType cmin, cmax;
  cmin = (*nodeCiter)->coords();
  cmax = (*nodeCiter)->coords();

  for ( ;
        nodeCiter != nodeEnd; ++nodeCiter )
  {
    cmin = min(cmin, (*nodeCiter)->coords());
    cmax = max(cmax, (*nodeCiter)->coords());
  } // next nodeCiter

  return std::make_pair(cmin, cmax);
}

#endif
