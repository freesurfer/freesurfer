#ifndef VENTFIX_H
#define VENTFIX_H

#include "mri.h"
#include "volcluster.h"

#include "pointset.h"

class VentFix
{
public:
  static MRI *fixasegps(MRI *asegps, MRI *brainmask, char *segids, float threshmin, int niters, int nmax, int topo);
  static MRI *ExpandSegIndices(MRI *seg, int segid, MRI *mask, int niters, int nmax, int topo, int *nexpansions, fsPointSet &centroid, MRI* newseg);
};

#endif
