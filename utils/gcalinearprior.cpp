/**
 * @brief Class to hold a volume of GCA priors in linear memory
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

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
using namespace std;

#include "minmax.hpp"

#include "gcalinearprior.hpp"

namespace Freesurfer
{
// ==========================================
void GCAlinearPrior::PrintStats(ostream &os) const
{
  os << "Stats for GCAlinearPrior" << endl;
  os << "  Exhumation time = " << exhumeTime << " ms" << endl;
  os << "  Inhumation time = " << inhumeTime << " ms" << endl;
}

// ==========================================
void GCAlinearPrior::ExtractDims(const GCA *const src)
{
  /*!
    Fills in the dimensions required from the given GCA.
    Does this by looping over all voxels, and finding the
    maximum number of elements allocated for each.
  */

  this->xDim = src->prior_width;
  this->yDim = src->prior_height;
  this->zDim = src->prior_depth;

  this->n4D = 0;

  for (int ix = 0; ix < this->xDim; ix++) {
    for (int iy = 0; iy < this->yDim; iy++) {
      for (int iz = 0; iz < this->zDim; iz++) {
        this->n4D += src->priors[ix][iy][iz].nlabels;
      }
    }
  }
}

// ==========================================
void GCAlinearPrior::Allocate(void)
{
  /*!
    Allocates the main arrays according to the 'dimension'
    members.
    All vectors are cleared first, and filled with defaults
  */

  this->bytes = 0;

  const size_t nVoxels = this->xDim * this->yDim * this->zDim;

  //! Space for the offsets4D array
  this->offsets4D.clear();
  this->offsets4D.resize(nVoxels + 1, numeric_limits< size_t >::max());
  this->bytes += this->offsets4D.size() * sizeof(unsigned int);

  //! Space for maxLabels
  this->maxLabels.clear();
  this->maxLabels.resize(nVoxels, 0);
  this->bytes += this->maxLabels.size() * sizeof(short);

  //! Space for the labels
  this->labels.clear();
  this->labels.resize(this->n4D, numeric_limits< unsigned short >::max());
  this->bytes += this->labels.size() * sizeof(unsigned short);

  //! Space for the priors
  this->priors.clear();
  this->priors.resize(this->n4D, numeric_limits< float >::quiet_NaN());
  this->bytes += this->priors.size() * sizeof(float);

  //! Space for totTraining
  this->totTraining.clear();
  this->totTraining.resize(nVoxels, numeric_limits< int >::max());
  this->bytes += this->totTraining.size() * sizeof(int);
}

// ==========================================

void GCAlinearPrior::Exhume(const GCA *const src)
{
  /*!
    This method is responsible for extracting GCA_PRIOR data
    from the source GCA, and packing into the linear arrays
  */

  Timer tExhume;

  this->ExtractDims(src);
  this->Allocate();

  std::vector< unsigned int > labelCounts(this->totTraining.size(), numeric_limits< unsigned int >::max());

  // Handle the 3D data
  for (int ix = 0; ix < this->xDim; ix++) {
    for (int iy = 0; iy < this->yDim; iy++) {
      for (int iz = 0; iz < this->zDim; iz++) {
        const GCA_PRIOR *const gcap = &(src->priors[ix][iy][iz]);

        this->maxVoxelLabel(ix, iy, iz) = gcap->max_labels;
        this->totalTraining(ix, iy, iz) = gcap->total_training;

        labelCounts.at(this->index3D(ix, iy, iz)) = gcap->nlabels;
      }
    }
  }

  // Compute 4D offsets
  this->offsets4D.at(0) = 0;
  std::partial_sum(labelCounts.begin(), labelCounts.end(), ++(this->offsets4D.begin()));

  // Handle the 4D data
  for (int ix = 0; ix < this->xDim; ix++) {
    for (int iy = 0; iy < this->yDim; iy++) {
      for (int iz = 0; iz < this->zDim; iz++) {
        const GCA_PRIOR *const gcap = &(src->priors[ix][iy][iz]);

        for (int iLabel = 0; iLabel < this->voxelLabelCount(ix, iy, iz); iLabel++) {
          this->voxelLabel(ix, iy, iz, iLabel) = gcap->labels[iLabel];
          this->voxelPrior(ix, iy, iz, iLabel) = gcap->priors[iLabel];
        }
      }
    }
  }

  this->exhumeTime = tExhume.milliseconds();
}

// ==========================================
const_GCAprior GCAlinearPrior::GetConstPrior(const int ix, const int iy, const int iz) const
{
  return (const_GCAprior(ix, iy, iz, *this));
}

// ==========================================
void GCAlinearPrior::Inhume(GCA *dst) const
{
  /*!
    Stores data about the priors back into the target
    GCA.
    Starts by destroying all of the GCA_PRIOR data
    which is there.
    One might prefer the subsequent reallocation
    to use the routines in gca.c, but that would
    probably involve writing them first.
  */

  Timer tInhume;

  // Dump the old data
  this->ScorchPriors(dst);

  // Set dimensions
  dst->prior_width = this->xDim;
  dst->prior_height = this->yDim;
  dst->prior_depth = this->zDim;

  // Start allocating
  dst->priors = (GCA_PRIOR ***)calloc(this->xDim, sizeof(GCA_PRIOR **));
  if (!(dst->priors)) {
    cerr << __FUNCTION__ << ": dst->priors allocation failed" << endl;
    exit(EXIT_FAILURE);
  }

  // Loop
  for (int ix = 0; ix < this->xDim; ix++) {
    // Allocate pointer block
    dst->priors[ix] = (GCA_PRIOR **)calloc(this->yDim, sizeof(GCA_PRIOR *));
    if (!(dst->priors[ix])) {
      cerr << __FUNCTION__ << ": dst->priors[ix] allocation failed" << endl;
      exit(EXIT_FAILURE);
    }

    for (int iy = 0; iy < this->yDim; iy++) {
      // Allocate pointer block
      dst->priors[ix][iy] = (GCA_PRIOR *)calloc(this->zDim, sizeof(GCA_PRIOR));
      if (!(dst->priors[ix][iy])) {
        cerr << __FUNCTION__ << ": dst->priors[ix][iy] allocation failed" << endl;
        exit(EXIT_FAILURE);
      }

      for (int iz = 0; iz < this->zDim; iz++) {
        GCA_PRIOR *gcap = &(dst->priors[ix][iy][iz]);

        const_GCAprior cGCAp = this->GetConstPrior(ix, iy, iz);

        gcap->nlabels = cGCAp.labelCount();
        gcap->max_labels = cGCAp.maxLabel();

        gcap->total_training = cGCAp.totalTraining();

        gcap->labels = (unsigned short *)calloc(gcap->nlabels, sizeof(unsigned short));
        gcap->priors = (float *)calloc(gcap->nlabels, sizeof(float));
        for (int iLabel = 0; iLabel < gcap->nlabels; iLabel++) {
          gcap->labels[iLabel] = cGCAp.labels(iLabel);
          gcap->priors[iLabel] = cGCAp.priors(iLabel);
        }
      }
    }
  }

  this->inhumeTime = tInhume.milliseconds();
}

// ==========================================

void GCAlinearPrior::ScorchPriors(GCA *targ) const
{
  /*!
    This method destroys the priors structure of a GCA,
    prior to inhumation of new data
  */
  for (int ix = 0; ix < targ->prior_width; ix++) {
    for (int iy = 0; iy < targ->prior_height; iy++) {
      for (int iz = 0; iz < targ->prior_depth; iz++) {
        free(targ->priors[ix][iy][iz].labels);
        free(targ->priors[ix][iy][iz].priors);
      }
      free(targ->priors[ix][iy]);
    }
    free(targ->priors[ix]);
  }
  free(targ->priors);
}
}  // namespace Freesurfer
