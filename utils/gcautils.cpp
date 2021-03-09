/**
 * @brief C++ GCA utilities
 *
 */
/*
 * Original Authors: Richard Edgar
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include <iostream>

#include "minmax.hpp"

#include "gcautils.hpp"

// ============================================================

namespace Freesurfer
{
void GetGCAstats(const GCA *const src)
{
  std::cout << "GCA Vital Statistics" << std::endl;

  std::cout << "ninputs    : " << src->ninputs << std::endl;

  GetGCAnodeStats(src);
  std::cout << std::endl;
  GetGCApriorStats(src);

  std::cout << "Statistics complete =======" << std::endl;
}

// --------------------------------------

void GetGCAnodeStats(const GCA *const src)
{
  std::cout << "Stats from nodes:" << std::endl;

  unsigned int nx, ny, nz;

  nx = src->node_width;
  ny = src->node_height;
  nz = src->node_depth;

  std::cout << "Dimensions : "
            << "( " << nx << ", " << ny << ", " << nz << " )" << std::endl;
  // ---

  MinMax< int > nLabelsNode, nMaxLabelsNode;
  MinMax< short > nLabelsGC1D[GIBBS_NEIGHBORHOOD];

  for (unsigned int ix = 0; ix < nx; ix++) {
    for (unsigned int iy = 0; iy < ny; iy++) {
      for (unsigned int iz = 0; iz < nz; iz++) {
        const GCA_NODE *const gcan = &(src->nodes[ix][iy][iz]);

        nLabelsNode.Accumulate(gcan->nlabels);
        nLabelsNode.Accumulate(gcan->max_labels);

        for (int iGC1D = 0; iGC1D < gcan->nlabels; iGC1D++) {
          const GC1D *const gc1d = &(gcan->gcs[iGC1D]);

          for (unsigned int iNeighbour = 0; iNeighbour < GIBBS_NEIGHBORHOOD; iNeighbour++) {
            nLabelsGC1D[iNeighbour].Accumulate(gc1d->nlabels[iNeighbour]);
          }
        }
      }
    }
  }

  std::cout << "GCA_NODE :" << std::endl;
  std::cout << "  nLabels   : " << nLabelsNode << std::endl;
  std::cout << "  maxLabels : " << nMaxLabelsNode << std::endl;

  std::cout << "  GC1D : " << std::endl;
  for (unsigned int i = 0; i < GIBBS_NEIGHBORHOOD; i++) {
    std::cout << "    nlabels[" << i << "] : " << nLabelsGC1D[i] << std::endl;
  }
}

// --------------------------------------

void GetGCApriorStats(const GCA *const src)
{
  std::cout << "Stats from priors:" << std::endl;

  unsigned int nx, ny, nz;

  nx = src->prior_width;
  ny = src->prior_height;
  nz = src->prior_depth;

  std::cout << "Dimensions : "
            << "( " << nx << ", " << ny << ", " << nz << " )" << std::endl;

  // ---

  MinMax< short > nLabelsPrior;
  MinMax< short > maxLabelsPrior;

  for (unsigned int ix = 0; ix < nx; ix++) {
    for (unsigned int iy = 0; iy < ny; iy++) {
      for (unsigned int iz = 0; iz < nz; iz++) {
        const GCA_PRIOR *const gcap = &(src->priors[ix][iy][iz]);

        nLabelsPrior.Accumulate(gcap->nlabels);
        maxLabelsPrior.Accumulate(gcap->max_labels);
      }
    }
  }

  std::cout << "GCA_PRIOR :" << std::endl;
  std::cout << "  nLabels   : " << nLabelsPrior << std::endl;
  std::cout << "  maxLabels : " << maxLabelsPrior << std::endl;
}
}  // namespace Freesurfer
