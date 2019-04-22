/*

Gheorghe Postelnicu, 2007

*/

#include <algorithm>
#include <iterator>
#include <numeric>

#include "mesh_bfs.h"

#include "pbCluster_mesh_crop.h"

#include "fem_3d.h"
#include "solver.h"

//--------------------------------------------------
//
// MeshCrop class
//
// takes a set of element indices as input
// (as well as a mesh)
//
//
//  and produces a cropped mesh
//
//  as well as a correspondence table for the elements
//  and nodes which are common
//----------------

template<class Cstr, int n>
class TMeshCrop
{
public:
  typedef TMesh<Cstr,n> MeshType;

  typedef typename MeshType::tNode NodeType;
  typedef typename MeshType::tElement ElementType;
  typedef typename MeshType::tCoords CoordsType;

  typedef std::vector<unsigned int> IndexContainerType;
  typedef std::set<unsigned int> IndexSetType;
  typedef std::map<unsigned int, unsigned int> MapType;

  TMeshCrop(const MeshType* pmesh) : m_pmesh(pmesh)
{}

  void Crop(MeshType* outMesh,
            const IndexContainerType& i_vidxElements,
            MapType& o_nodeMap,
            MapType& o_eltMap,
            IndexContainerType& o_boundaryNodes
           );

private:

  const MeshType* m_pmesh;

  std::pair<CoordsType,CoordsType>
  ComputeBoundingBox(MeshType* pmesh);

  // assembles all the node indices (in the master array)
  // that will be part of the cropped mesh
  typedef std::back_insert_iterator<IndexContainerType> IndexContainerII;

  void PopulateNodes( const IndexContainerType& velts,
                      IndexContainerII ii_nodes
                    );
  void PopulateBoundaryNodes( const IndexContainerType& velts,
                              const IndexContainerType& vnodes,
                              IndexContainerII  ii_boundary);
};



//--------------------------------------------------



TopologySolver::TopologySolver(TMesh3d& mesh,
                               unsigned int radius,
                               double de,
                               double dnu) :
    m_mesh(mesh), m_radius(radius), m_dE(de), m_dnu(dnu)
{
  this->BuildClusters();
  this->MergeIntersectingClusters();
  this->DoLocalSmoothing();
}

/*

performs a radius BFS around each element with topology problems.
The indices of the resulting elements are then stored in a list as
smart pointers.


!!!! Each vector of element indices is sorted - this is important
for the set operations to follow

*/
void
TopologySolver::BuildClusters()
{
  TMeshBfs<TMesh3d> bfs(&m_mesh);
  IndexVectorType vecOrigProblems;

  try
  {
    m_mesh.orientation_pb(dst, std::back_inserter(vecOrigProblems) );
  }
  catch (const std::exception& e)
  {
    std::cout
    << " Exception caught while gathering info on topology problems\n"
    << e.what() << std::endl;
  }

  // build container of topology problem indices
  for ( IndexVectorType::const_iterator cit = vecOrigProblems.begin();
        cit != vecOrigProblems.end(); ++cit)
  {
    m_clusterContainer.push_back(
      std::shared_ptr<IndexVectorType>(new IndexVectorType ) );
    bfs.Visit( *cit,
               m_radius,
               std::back_inserter<IndexVectorType>(
                 *(m_clusterContainer.back()) )
             );
    std::sort( m_clusterContainer.back()->begin(),
               m_clusterContainer.back()->end() );
  } // next cit
}

typedef std::vector<unsigned int> VectorType;
class IntersectsWith : public std::unary_function<VectorType, bool>
{
  TopologySolver::IndexVectorPointer& m_container;
public:
  IntersectsWith(TopologySolver::IndexVectorPointer& c) : m_container(c)
  {}

  bool operator()(const TopologySolver::IndexVectorPointer other) const
  {
    TopologySolver::IndexVectorType vecIntersection;
    std::set_intersection( other->begin(), other->end(),
                           m_container->begin(), m_container->end(),
                           std::back_inserter(vecIntersection)
                         );
    return (!vecIntersection.empty());

  }
};

template<class Container>
struct RemoveIfAndCopy
{
  RemoveIfAndCopy(Container& c) : container(c)
  {}
  Container& container;
  template<typename _ForwardIterator, typename _OutputIterator,
  typename _Predicate>
  void Execute( _ForwardIterator __first, _ForwardIterator __last,
                _OutputIterator __out, _Predicate __pred)
  {
    while ( __first != __last )
    {
      if ( __pred(*__first) )
      {
        *__out = *__first;
        __first = container.erase(__first);
      }
      else ++__first;
    } // next it
  }
};

TopologySolver::IndexVectorPointer
do_union(TopologySolver::IndexVectorPointer& v1,
         TopologySolver::IndexVectorPointer& v2)
{
  TopologySolver::IndexVectorPointer v3(new TopologySolver::IndexVectorType);

  std::set_union( v1->begin(), v1->end(),
                  v2->begin(), v2->end(),
                  std::back_inserter(*v3)
                );
  return v3;
}

void
TopologySolver::MergeIntersectingClusters()
{
  RemoveIfAndCopy<ClusterContainerType> listFilter(m_clusterContainer);

  unsigned int counter = 0;
  for (ClusterContainerType::iterator it = m_clusterContainer.begin();
       it != m_clusterContainer.end();
       ++it, ++counter )
  {
    bool bAgain = true;
    while ( bAgain )
    {
      ClusterContainerType::iterator it2 = it;
      ClusterContainerType tmpContainer;

      listFilter.Execute
      ( ++it2, m_clusterContainer.end(),
        std::back_inserter<ClusterContainerType>(tmpContainer),
        IntersectsWith(*it)
      );

      if ( !tmpContainer.empty() )
      {
        IndexVectorPointer newCluster
        = std::accumulate( tmpContainer.begin(), tmpContainer.end(),
                           IndexVectorPointer(new IndexVectorType),
                           do_union );
        *it = newCluster;

      }
      else bAgain = false;
    } // end while bAgain
  } // next it

  std::cout << " TopologyClusters = "
  << m_clusterContainer.size() << std::endl;
}

void
TopologySolver::DoLocalSmoothing()
{
  typedef TMeshCrop<Constructor,3> MeshCropFilter;
  typedef MeshCropFilter::MapType MapType;
  typedef MeshCropFilter::IndexContainerType  Vector;

  typedef TDirectSolver<Constructor,3> SolverType;

  MeshCropFilter mesh_crop(&m_mesh);


  for ( ClusterContainerType::const_iterator cit = m_clusterContainer.begin();
        cit != m_clusterContainer.end();
        ++cit )
  {
    // use the cluster elements to create a cropped mesh
    CMesh3d cropped;
    MapType nodeMap, elementMap;
    Vector boundary;

    mesh_crop.Crop( &cropped,
                    **cit,
                    nodeMap,
                    elementMap,
                    boundary );

    cropped.build_index_src();
    cropped.set_constants( m_dE, m_dnu );

    SolverType solver;
    solver.set_mesh( &cropped );

    {
      CMesh3d::tNode* pnode = NULL;
      for ( Vector::const_iterator  cit = boundary.begin();
            cit != boundary.end(); ++cit )
      {
        cropped.get_node( *cit, &pnode );
        if ( !pnode )
          throw " TopologySolver - NULL node ";
        solver.add_bc_natural( pnode,
                               pnode->delta() );
      }
    }

    solver.solve();

    // pass on results to the original mesh
    {
      CMesh3d::tNode* pnode = NULL, *pnode_orig = NULL;
      for (unsigned int ui(0), nnodes(cropped.get_no_nodes());
           ui < nnodes; ++ui )
      {
        cropped.get_node(ui, &pnode);
        m_mesh.get_node( nodeMap[ui], &pnode_orig );
        pnode_orig->set_delta( pnode->delta() );
      } // next ui
    }

  } // next cit
}

//--------------------------------------------------------------
//
//--------------------------------------------------------------


template<class Cstr, int n>
void
TMeshCrop<Cstr,n>::Crop(MeshType* outMesh,
                        const IndexContainerType& i_vidxElements,
                        MapType& o_nodeMap,
                        MapType& o_eltMap,
                        IndexContainerType& o_boundaryNodes
                       )
{
  NodeType* pnodeIn;
  NodeType* pnodeOut;

  IndexContainerType vnodes; // indices of the concerned nodes in
  // the master mesh

  this->PopulateNodes( i_vidxElements,
                       std::back_inserter(vnodes)
                     );

  IndexContainerType origBoundaryNodes;
  this->PopulateBoundaryNodes( i_vidxElements,
                               vnodes,
                               std::back_inserter(origBoundaryNodes)
                             );
  origBoundaryNodes.erase
  (
    std::unique( origBoundaryNodes.begin(), origBoundaryNodes.end() ),
    origBoundaryNodes.end()
  );

  // populate the new mesh

  // keep track of the correspondence
  MapType internalMap; //  [ orig Mesh index ] = new Mesh index

  // populate nodes
  {
    unsigned int index = 0;
    for (typename IndexContainerType::const_iterator cit
         = vnodes.begin();
         cit != vnodes.end();
         ++cit, ++index )
    {
      m_pmesh->get_node( *cit, &pnodeIn );
      if ( !pnodeIn )
        throw " TMeshCrop - NULL node while cropping ";
      pnodeOut = Cstr::node(index,
                            pnodeIn->coords(),
                            pnodeIn->delta(),
                            pnodeIn->is_active() );
      outMesh->add_node(pnodeOut);

      internalMap[ *cit ] = index;
      o_nodeMap[ index ] = *cit;
    } // next cit, index
  }

  // populate elements
  {
    unsigned int index = 0;
    const ElementType* cpelt;
    std::vector<NodeType*> vpnodes;
    typename MapType::const_iterator map_cit;

    for ( typename IndexContainerType::const_iterator cit
          = i_vidxElements.begin();
          cit != i_vidxElements.end();
          ++cit, ++index )
    {
      vpnodes.clear();
      cpelt = m_pmesh->get_elt(*cit);
      if ( !cpelt )
        throw " TMeshCrop - NULL cpelt while cropping ";
      for (unsigned int ui(0), nnodes(cpelt->no_nodes());
           ui < nnodes; ++ui)
      {
        cpelt->get_node(ui, &pnodeIn);
        if ( !pnodeIn ) throw " TMeshCrop - NULL nodeIn ";
        map_cit = internalMap.find( pnodeIn->id() );
        if ( map_cit == internalMap.end() )
          throw " TMeshCrop - NULL extraction from map ";
        outMesh->get_node( map_cit->second, &pnodeOut );
        vpnodes.push_back( pnodeOut );
      }

      ElementType* pelt = Cstr::elt(index,
                                    NULL,
                                    vpnodes);
      outMesh->add_elt( pelt );
      o_eltMap[index] = cpelt->id();
    } // next cit, index
  }

  // compute the bounding box
  {
    std::pair<CoordsType,CoordsType> bbox;
    bbox = this->ComputeBoundingBox( outMesh );

    outMesh->m_cmin = bbox.first;
    outMesh->m_cmax = bbox.second;
  }

  for ( typename IndexContainerType::const_iterator cit =
          origBoundaryNodes.begin();
        cit != origBoundaryNodes.end(); ++cit )
  {
    typename MapType::const_iterator map_cit = internalMap.find( *cit );
    if ( map_cit == internalMap.end() )
      throw " Crop - exception during query";
    o_boundaryNodes.push_back( map_cit->second );
  }

}

template<class Cstr, int n>
void
TMeshCrop<Cstr,n>::PopulateNodes( const IndexContainerType& velts,
                                  IndexContainerII ii_nodes
                                )
{
  for ( typename IndexContainerType::const_iterator cit
        = velts.begin();
        cit != velts.end();
        ++cit )
  {
    NodeType* pnode = NULL;
    const ElementType* cpelt
    = m_pmesh->get_elt(*cit);

    for (unsigned int ui(0), nnodes(cpelt->no_nodes());
         ui < nnodes; ++ui)
    {
      cpelt->get_node(ui, &pnode);
      if ( !pnode )
        throw " TMeshCrop::PopulateNodes - NULL node ";
      ii_nodes = pnode->get_id();
    }
  }

}

/*

the boundary nodes are those which also belong
to elements other than the ones present in the container

*/
template<class Cstr, int n>
void
TMeshCrop<Cstr,n>::PopulateBoundaryNodes( const IndexContainerType& velts,
    const IndexContainerType& vnodes,
    IndexContainerII  ii_boundary)
{
  typename NodeType::tElt_citer elt_cit, elt_end;

  for ( typename IndexContainerType::const_iterator node_cit = vnodes.begin();
        node_cit != vnodes.end(); ++node_cit )
  {
    NodeType* pnode = NULL;
    m_pmesh->get_node( *node_cit, &pnode );
    if ( !pnode )
      throw " TMeshCrop::PopulateBoundaryNodes - NULL node ";
    for ( pnode->get_elt_citer(elt_cit, elt_end);
          elt_cit != elt_end;
          ++elt_cit )
    {
      unsigned int eltCit = (unsigned int)*elt_cit;
      if ( !std::binary_search( velts.begin(), velts.end(), eltCit ) )
        ii_boundary = *node_cit;
    }
  } // next cit
}

template<class Cstr, int n>
std::pair<typename TMeshCrop<Cstr,n>::CoordsType,
typename TMeshCrop<Cstr,n>::CoordsType>
TMeshCrop<Cstr,n>::ComputeBoundingBox(MeshType* pmesh)
{
  typename MeshType::NodeConstIterator node_cit, node_end;
  pmesh->get_node_citer( node_cit, node_end );
  CoordsType cmin, cmax;

  cmin = (*node_cit)->coords();
  cmax = (*node_cit)->coords();

  for ( ; node_cit != node_end; ++node_cit )
  {
    cmin = min( cmin, (*node_cit)->coords() );
    cmax = max( cmax, (*node_cit)->coords() );
  }

  return std::make_pair(cmin, cmax);
}
