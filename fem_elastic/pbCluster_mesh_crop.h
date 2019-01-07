/*

Gheorghe Postelnicu, 2007

*/
#ifndef _h_pbCluster_mesh_crop_h_
#define _h_pbCluster_mesh_crop_h_

#include <list>
#include <vector>
#include <memory>

#include "fem_3d.h"

class TopologySolver
{
public:
  typedef std::vector<unsigned int> IndexVectorType;
  typedef std::shared_ptr<IndexVectorType> IndexVectorPointer;
  typedef std::list<IndexVectorPointer> ClusterContainerType;

  TopologySolver(TMesh3d& mesh, unsigned int radius=2,
                 double de=1.0, double dnu=0.3);

private:
  void BuildClusters();
  void MergeIntersectingClusters();

  void DoLocalSmoothing();

  TMesh3d& m_mesh;
  unsigned int m_radius;
  double m_dE, m_dnu;

  ClusterContainerType m_clusterContainer;
};

#endif
