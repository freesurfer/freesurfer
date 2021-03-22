/**
 * @brief Class to hold a volume of GCA nodes in linear memory
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

#include "gcalinearnode.hpp"

namespace Freesurfer
{
// ====================================================
void GCAlinearNode::PrintStats(ostream &os) const
{
  os << "Stats for GCAlinearNode" << endl;
  os << "  Exhumation time = " << exhumeTime << " ms" << endl;
  os << "  Inhumation time = " << inhumeTime << " ms" << endl;
}

// ==========================================
void GCAlinearNode::Exhume(const GCA *const src)
{
  /*!
    This method is responsible for extracting data from
    the source GCA, and packing into the linear arrays
  */

  Timer tExhume;

  this->ExtractDims(src);
  this->Allocate();

  // Array to hold label count at each node
  vector< unsigned int > nLabelsNode(this->nodeTotalTraining.size(), numeric_limits< unsigned int >::max());

  // Handle the 3D data
  for (int ix = 0; ix < this->xDim; ix++) {
    for (int iy = 0; iy < this->yDim; iy++) {
      for (int iz = 0; iz < this->zDim; iz++) {
        const GCA_NODE *const gcan = &(src->nodes[ix][iy][iz]);

        nLabelsNode.at(this->index3D(ix, iy, iz)) = gcan->nlabels;

        this->totalTraining(ix, iy, iz) = gcan->total_training;
        this->maxLabels(ix, iy, iz) = gcan->max_labels;
      }
    }
  }

  // Compute 4D offsets
  this->offsets4D.at(0) = 0;
  std::partial_sum(nLabelsNode.begin(), nLabelsNode.end(), ++(this->offsets4D.begin()));

  // Array to hold directions for each GC1D
  vector< unsigned int > nDirecAtGC1D(this->n4D, numeric_limits< unsigned int >::max());

  // Handle the 4D data
  for (int ix = 0; ix < this->xDim; ix++) {
    for (int iy = 0; iy < this->yDim; iy++) {
      for (int iz = 0; iz < this->zDim; iz++) {
        const GCA_NODE *const gcan = &(src->nodes[ix][iy][iz]);

        for (int iGC1D = 0; iGC1D < this->gc1dCount(ix, iy, iz); iGC1D++) {
          nDirecAtGC1D.at(this->index4D(ix, iy, iz, iGC1D)) = this->gc1dNeighbourDim;

          // Deal with the label held inside each node
          this->labelsAtNode(ix, iy, iz, iGC1D) = gcan->labels[iGC1D];

          // Grab hold of the current GC1D
          const GC1D *const gc1d = &(gcan->gcs[iGC1D]);

          // Deal with the means and variances
          // Recall that ninputs==1
          this->meansAtNodeGC1D(ix, iy, iz, iGC1D) = gc1d->means[0];
          this->variancesAtNodeGC1D(ix, iy, iz, iGC1D) = gc1d->covars[0];

          this->nJustPriorsAtNodeGC1D(ix, iy, iz, iGC1D) = gc1d->n_just_priors;
          this->nTrainingAtNodeGC1D(ix, iy, iz, iGC1D) = gc1d->ntraining;
          this->regularisedAtNodeGC1D(ix, iy, iz, iGC1D) = gc1d->regularized;
        }
      }
    }
  }

  if (this->hasGibbsNeighbourhood) {
    // Compute 5D offsets
    this->offsets5D.at(0) = 0;
    std::partial_sum(nDirecAtGC1D.begin(), nDirecAtGC1D.end(), ++(this->offsets5D.begin()));

    // Array to hold labels for each direction at each GC1D
    vector< unsigned int > nLabelsDirecAtGC1D(this->n5D, numeric_limits< unsigned int >::max());

    // Handle 5D data
    for (int ix = 0; ix < this->xDim; ix++) {
      for (int iy = 0; iy < this->yDim; iy++) {
        for (int iz = 0; iz < this->zDim; iz++) {
          const GCA_NODE *const gcan = &(src->nodes[ix][iy][iz]);

          for (int iGC1D = 0; iGC1D < this->gc1dCount(ix, iy, iz); iGC1D++) {
            const GC1D *const gc1d = &(gcan->gcs[iGC1D]);

            for (int iDirec = 0; iDirec < this->gc1dNeighbourDim; iDirec++) {
              nLabelsDirecAtGC1D.at(this->index5D(ix, iy, iz, iGC1D, iDirec)) = gc1d->nlabels[iDirec];
            }
          }
        }
      }
    }

    // Compute 6D offsets
    this->offsets6D.at(0) = 0;
    std::partial_sum(nLabelsDirecAtGC1D.begin(), nLabelsDirecAtGC1D.end(), ++(this->offsets6D.begin()));

    // Handle 6D data
    for (int ix = 0; ix < this->xDim; ix++) {
      for (int iy = 0; iy < this->yDim; iy++) {
        for (int iz = 0; iz < this->zDim; iz++) {
          const GCA_NODE *const gcan = &(src->nodes[ix][iy][iz]);

          for (int iGC1D = 0; iGC1D < this->gc1dCount(ix, iy, iz); iGC1D++) {
            const GC1D *const gc1d = &(gcan->gcs[iGC1D]);
            for (int iDirec = 0; iDirec < this->gc1dNeighbourDim; iDirec++) {
              for (int iLabel = 0; iLabel < this->nLabelsAtNodeGC1Ddirection(ix, iy, iz, iGC1D, iDirec); iLabel++) {
                this->labelsAtNodeGC1Ddirection(ix, iy, iz, iGC1D, iDirec, iLabel) = gc1d->labels[iDirec][iLabel];
                this->labelPriorsAtNodeGC1Ddirection(ix, iy, iz, iGC1D, iDirec, iLabel) =
                    gc1d->label_priors[iDirec][iLabel];
              }
            }
          }
        }
        // End of per-voxel loop
      }
    }
  }

  this->exhumeTime = tExhume.milliseconds();
}

// ====================================================

const_GCAnode GCAlinearNode::GetConstNode(const int ix, const int iy, const int iz) const
{
  return (const_GCAnode(ix, iy, iz, *this));
}

// ====================================================

void GCAlinearNode::Inhume(GCA *dst) const
{
  /*!
    Stores data back from the linear arrays into a GCA.
    This removes the existing node structure, and reallocates
    it iself.
    One might prefer this not to happen, and to use the
    routines from gca.c.
    However, that would involve writing some of them first.
  */

  Timer tInhume;

  // Dispose of the old node data
  this->ScorchNodes(dst);

  // Set dimensions
  dst->ninputs = 1;
  dst->node_width = this->xDim;
  dst->node_height = this->yDim;
  dst->node_depth = this->zDim;

  // Start allocating
  dst->nodes = (GCA_NODE ***)calloc(this->xDim, sizeof(GCA_NODE **));
  if (!(dst->nodes)) {
    cerr << __FUNCTION__ << ": dst->nodes allocation failed" << endl;
    exit(EXIT_FAILURE);
  }

  for (int ix = 0; ix < this->xDim; ix++) {
    // Allocate pointer block
    dst->nodes[ix] = (GCA_NODE **)calloc(this->yDim, sizeof(GCA_NODE *));
    if (!(dst->nodes[ix])) {
      cerr << __FUNCTION__ << ": dst->nodes[" << ix << "] allocation failed" << endl;
      exit(EXIT_FAILURE);
    }

    for (int iy = 0; iy < this->yDim; iy++) {
      // Allocate pointer block
      dst->nodes[ix][iy] = (GCA_NODE *)calloc(this->zDim, sizeof(GCA_NODE));
      if (!(dst->nodes[ix][iy])) {
        cerr << __FUNCTION__ << ": dst->nodes"
             << "[" << ix << "]"
             << "[" << iy << "]"
             << " allocation failed" << endl;
        exit(EXIT_FAILURE);
      }

      for (int iz = 0; iz < this->zDim; iz++) {
        // Allocate pointer block
        GCA_NODE *const gcan = &(dst->nodes[ix][iy][iz]);
        const const_GCAnode gln = this->GetConstNode(ix, iy, iz);

        gcan->nlabels = gln.gc1dCount();
        gcan->max_labels = gln.maxLabels();
        gcan->total_training = gln.totalTraining();

        // Allocate labels array
        gcan->labels = (unsigned short *)calloc(gln.gc1dCount(), sizeof(unsigned short));
        if (!(gcan->labels)) {
          cerr << __FUNCTION__ << ": dst->nodes"
               << "[" << ix << "]"
               << "[" << iy << "]"
               << "[" << iz << "].labels"
               << " allocation failed" << endl;
          exit(EXIT_FAILURE);
        }

        // Allocate GC1D array
        gcan->gcs = (GC1D *)calloc(this->gc1dCount(ix, iy, iz), sizeof(GC1D));
        if (!(gcan->gcs)) {
          cerr << __FUNCTION__ << ": dst->nodes"
               << "[" << ix << "]"
               << "[" << iy << "]"
               << "[" << iz << "].gcs"
               << " allocation failed" << endl;
          exit(EXIT_FAILURE);
        }

        // Loop over the GC1Ds
        for (int iGC1D = 0; iGC1D < gln.gc1dCount(); iGC1D++) {
          // Do the labels on the side
          gcan->labels[iGC1D] = gln.labels(iGC1D);

          GC1D *const gc1d = &(gcan->gcs[iGC1D]);
          const const_GCAnode_GC1D g1d = gln.GetConstGC1D(iGC1D);

          gc1d->means = (float *)calloc(dst->ninputs,  // Always 1
                                        sizeof(float));
          if (!(gc1d->means)) {
            cerr << __FUNCTION__ << ": Allocation failure of means" << endl;
            exit(EXIT_FAILURE);
          }

          gc1d->covars = (float *)calloc(dst->ninputs,  // Always 1
                                         sizeof(float));
          if (!(gc1d->covars)) {
            cerr << __FUNCTION__ << ": Allocation failure of covars" << endl;
            exit(EXIT_FAILURE);
          }

          // Do the mean and variance (recall ninputs==1)
          gc1d->means[0] = g1d.mean();
          gc1d->covars[0] = g1d.variance();

          gc1d->n_just_priors = g1d.nJustPriors();
          gc1d->ntraining = g1d.nTraining();
          gc1d->regularized = g1d.regularised();

          if (this->hasGibbsNeighbourhood) {
            // Allocate the nlabels array
            gc1d->nlabels = (short *)calloc(this->gc1dNeighbourDim,  // Always 6/GIBBS_NEIGHBORHOOD
                                            sizeof(short));
            if (!(gc1d->nlabels)) {
              cerr << __FUNCTION__ << ": Allocation failure of nlabels" << endl;
              exit(EXIT_FAILURE);
            }

            // Allocate pointers for label_priors
            gc1d->label_priors = (float **)calloc(this->gc1dNeighbourDim, sizeof(float *));
            if (!(gc1d->label_priors)) {
              cerr << __FUNCTION__ << ": Allocation failure of label_priors" << endl;
              exit(EXIT_FAILURE);
            }

            // Allocate pointers for labels
            gc1d->labels = (unsigned short **)calloc(this->gc1dNeighbourDim, sizeof(unsigned short *));
            if (!(gc1d->labels)) {
              cerr << __FUNCTION__ << ": Allocation failure of labels" << endl;
              exit(EXIT_FAILURE);
            }

            for (int iDirec = 0; iDirec < this->gc1dNeighbourDim; iDirec++) {
              // Set the number
              gc1d->nlabels[iDirec] = g1d.nLabels(iDirec);

              // Allocate the memory
              gc1d->label_priors[iDirec] = (float *)calloc(gc1d->nlabels[iDirec], sizeof(float));
              if (!(gc1d->label_priors[iDirec])) {
                cerr << __FUNCTION__ << ": Allocation failure of label_priors" << endl;
                exit(EXIT_FAILURE);
              }

              gc1d->labels[iDirec] = (unsigned short *)calloc(gc1d->nlabels[iDirec], sizeof(unsigned short));
              if (!(gc1d->labels[iDirec])) {
                cerr << __FUNCTION__ << ": Allocation failure of labels" << endl;
                exit(EXIT_FAILURE);
              }

              for (int iLabel = 0; iLabel < g1d.nLabels(iDirec); iLabel++) {
                gc1d->labels[iDirec][iLabel] = g1d.labels(iDirec, iLabel);
                gc1d->label_priors[iDirec][iLabel] = g1d.labelPriors(iDirec, iLabel);
              }
            }
          }
        }
      }
    }
  }

  this->inhumeTime = tInhume.milliseconds();
}

// ====================================================

void GCAlinearNode::Allocate(void)
{
  /*!
    Allocates the main arrays according to
    the values stored in the dimension member
    variables.
    All vectors are cleared first, in order to
    fill with defaults
  */

  const size_t nVoxels = this->xDim * this->yDim * this->zDim;

  // Allocate the offset arrays
  this->offsets4D.resize(nVoxels + 1, numeric_limits< size_t >::max());
  this->offsets5D.resize(this->n4D + 1, numeric_limits< size_t >::max());
  this->offsets6D.resize(this->n5D + 1, numeric_limits< size_t >::max());

  // Allocate the 3D arrays
  this->nodeMaxLabels.resize(nVoxels, numeric_limits< int >::max());
  this->nodeTotalTraining.resize(nVoxels, numeric_limits< int >::max());

  // Allocate the 4D arrays
  this->nodeLabels.resize(this->n4D, numeric_limits< unsigned short >::max());
  this->means.resize(this->n4D, numeric_limits< float >::quiet_NaN());
  this->variances.resize(this->n4D, numeric_limits< float >::quiet_NaN());
  this->nJustPriors.resize(this->n4D, numeric_limits< short >::max());
  this->nTraining.resize(this->n4D, numeric_limits< int >::max());
  this->regularised.resize(this->n4D, numeric_limits< char >::max());

  // Allocate the 6D arrays
  this->gc1dDirecLabelPriors.resize(this->n6D, numeric_limits< float >::quiet_NaN());
  this->gc1dDirecLabels.resize(this->n6D, numeric_limits< unsigned short >::max());
}

// ====================================================

void GCAlinearNode::ExtractDims(const GCA *const src)
{
  /*!
    Fills in the dimensions required from the given GCA.
    Does this by looping over all voxels, and finding the
    maximum number of elements allocated for each.
  */
  if (src->ninputs != 1) {
    cerr << __FUNCTION__ << ": Must have ninputs==1!" << endl;
    exit(EXIT_FAILURE);
  }

  this->xDim = src->node_width;
  this->yDim = src->node_height;
  this->zDim = src->node_depth;

  this->n4D = 0;
  this->n5D = 0;
  this->n6D = 0;

  for (int iz = 0; iz < this->zDim; iz++) {
    for (int iy = 0; iy < this->yDim; iy++) {
      for (int ix = 0; ix < this->xDim; ix++) {
        const GCA_NODE *const gcan = &(src->nodes[ix][iy][iz]);
        this->n4D += gcan->nlabels;

        for (int iGC1D = 0; iGC1D < gcan->nlabels; iGC1D++) {
          const GC1D *const gc1d = &(gcan->gcs[iGC1D]);
          this->n5D += GIBBS_NEIGHBORHOOD;

          if (!(src->flags & GCA_NO_MRF)) {
            for (unsigned int iNeighbour = 0; iNeighbour < GIBBS_NEIGHBORHOOD; iNeighbour++) {
              this->n6D += gc1d->nlabels[iNeighbour];
            }
          }
          else {
            this->hasGibbsNeighbourhood = false;
          }
        }
      }
    }
  }
}

// ====================================================

void GCAlinearNode::ScorchNodes(GCA *targ) const
{
  /*!
    Deletes all of the node related things from
    a GCA, prior to inhumation of new data
  */

  for (int ix = 0; ix < targ->node_width; ix++) {
    for (int iy = 0; iy < targ->node_height; iy++) {
      for (int iz = 0; iz < targ->node_depth; iz++) {
        GCANfree(&(targ->nodes[ix][iy][iz]), targ->ninputs);
      }
      free(targ->nodes[ix][iy]);
    }
    free(targ->nodes[ix]);
  }
  free(targ->nodes);
}

// ###############################################################

const_GCAnode_GC1D const_GCAnode::GetConstGC1D(const int iGC1D) const
{
  if ((iGC1D < 0) || (iGC1D >= this->myGC1Dcount)) {
    std::cerr << __FUNCTION__ << ": Out of range " << iGC1D << std::endl;
    abort();
  }

  const_GCAnode_GC1D gc1d(this->ix, this->iy, this->iz, iGC1D, this->gcaln);

  return (gc1d);
}
}  // namespace Freesurfer
