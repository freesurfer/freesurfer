#ifndef VENTFIX_H
#define VENTFIX_H

#include "mri.h"
#include "volcluster.h"

#include "pointset.h"

class VentFix
{
public:
  struct RELABELCRITERIA {
    int segA;            // clusters that we are interested
    int *adjSegs;        // clusters that are neighboring segA with topo constraint
    int nAdjSegs;        // number of elements in adjSegs and newsegids
    int topo;            // topology constraint, 
                         // 1=face neighbors only, 2=face neighbors + edge neighbors only, 3=face neighbors + edge neighbors + corner neighbors
    int *newsegids;      // optional, ndwsegids needs to have the same size as inclSegs, new segids that will be used to relabel clusters segA if conditions are met; 
    int *exclSegs;       // clusters that can not neighbor segA with topo constraint
    int nExclSegs;       // number of elements in exclSegs
    int nIters;          // number of iterations to look for cluster members
  };

  static MRI *fixasegps(MRI *asegps, MRI *brainmask, char *segids, float threshmin, int niters, int nmax, int topo);

  static MRI *relabelSegANeighboringSegB(MRI *asegps, RELABELCRITERIA *criteria, fsPointSet *centroid);

private:
  static MRI *ExpandSegIndices(MRI *seg, int segid, MRI *mask, int niters, int nmax, int topo, int *nexpansions, fsPointSet &centroid, MRI* newseg);

  static int hasneighbor(MRI* asegps, VOLCLUSTER *vc, RELABELCRITERIA*criteria, int *col, int *row, int *sli);
  static void relabelCluster(MRI *asegps, VOLCLUSTER *vc, int newsegid, MATRIX *vox2ras, fsPointSet *centroid);
  static int relabelClusterNiters(int c0, int r0, int s0, MRI *asegVol, int clustid, int allowdiag, int newsegid, int nIters, MATRIX *vox2ras, fsPointSet *centroid);
};

#endif
