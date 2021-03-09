/**
 * @brief functions for finding the best rotation for aligning a MRIS with a target
 *
 */
/*
 * Original Author: Bevin Brett
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
#include "MRISrigidBodyAlignGlobal.h"
#include "romp_support.h"
#include "vertexRotator.h"

static float* getFloats(size_t capacity) {
  void* ptr = NULL;
  int status = posix_memalign(&ptr, 64, capacity*sizeof(float));
  if (status) {
    fprintf(stderr, "%s:%d could not posix_memalign %ld aligned floats, status:%d\n", __FILE__, __LINE__, capacity, status);
    ptr = (float*)malloc(capacity*sizeof(float));
    if (!ptr) exit(1);
    fprintf(stderr, "%s:%d but malloc could\n", __FILE__, __LINE__);
  }
  return (float*)ptr;
}

void MRISrigidBodyAlignGlobal_findMinSSE(
  double* new_mina, double* new_minb, double* new_ming, double* new_sse,  // outputs
  MRI_SURFACE*       mris,
  INTEGRATION_PARMS* parms,
  float              min_radians,
  float              max_radians,
  double             ext_sse,
  int                nangles) {

  bool const tracing     = false;
  bool const spreadsheet = false;
  
  // Get all the non-ripped vertices
  // This yields cache-aligned dense coordinates to work with.
  //
  size_t verticesCapacity = mris->nvertices;

  float* const curv                = getFloats(verticesCapacity);

  float* const xv                  = getFloats(verticesCapacity);
  float* const yv                  = getFloats(verticesCapacity);
  float* const zv                  = getFloats(verticesCapacity);

  float* const gammaRotated_xv     = getFloats(verticesCapacity);
  float* const gammaRotated_yv     = getFloats(verticesCapacity);
  float* const gammaRotated_zv     = getFloats(verticesCapacity);

  float* const betaGammaRotated_xv = getFloats(verticesCapacity);
  float* const betaGammaRotated_yv = getFloats(verticesCapacity);
  float* const betaGammaRotated_zv = getFloats(verticesCapacity);

  size_t verticesSize = 0;
  int vno;
  for (vno = 0; vno < mris->nvertices; vno++) {
    VERTEX const * v = &mris->vertices[vno];
    if (v->ripflag) continue;
    xv[verticesSize] = v->x;
    yv[verticesSize] = v->y;
    zv[verticesSize] = v->z;
    curv[verticesSize] = v->curv;
    verticesSize++;
  }

  // Project these coordinates to the surface of a sphere of radius mris->radius
  // This assumes that the vertices are already well away from (0,0,0)
  //
  {
    for (unsigned int i = 0; i < verticesSize; i++) {
      float invLength = mris->radius/sqrtf(squaref(xv[i])+squaref(yv[i])+squaref(zv[i]));
      xv[i] *= invLength;
      yv[i] *= invLength;
      zv[i] *= invLength;
    }
  }

  // TODO sort on z
  // to maximize cache hit rate in the MRISPfunctionVal_radiusR function
  
  // Decide on the number of points along each edge of the grid
  // and what that is in radians
  //
  int gridSize_init = 1;
  while (max_radians / gridSize_init > min_radians) gridSize_init *= 2;
  int const gridSize = gridSize_init;
  int const gridCenterI = gridSize/2;
  float const radiansPerGridCell = max_radians / gridSize;
  
  #define iToRadians(I) (((I) - gridCenterI)*radiansPerGridCell)
  
  // Allocate enough done flags for a cube of this many points
  // Zero'ed since none visited yet
  //
  size_t const doneFlagsSize_inBits = gridSize*gridSize*gridSize;
  size_t const doneFlagsPerElt      = sizeof(long)*8;
  size_t const doneFlagsMask        = doneFlagsPerElt-1;
  size_t const doneFlagsSize        = (doneFlagsSize_inBits + doneFlagsPerElt - 1) / doneFlagsPerElt;
  long*  const doneFlags = (long*)calloc(doneFlagsSize, sizeof(long));
  
  // Support inverting the alpha and the vertices loops
  // The vertices are broken into partitions
  // so the inner loops can be threaded
  //
  int const numberOfVerticesPartitions     = 16;
  int const forAlphasCapacity              = nangles + 1;

  int*      ajsForAlphas                   = (int*)
                                             malloc(forAlphasCapacity*sizeof(int));
  float*    alphasForAlphas                = (float*)
                                             malloc(forAlphasCapacity*sizeof(float));

  int const forAlphasForPartitionsCapacity = forAlphasCapacity * numberOfVerticesPartitions;
  
  double*   ssesForAlphasForPartitions     = (double*)
                                             malloc(forAlphasForPartitionsCapacity*sizeof(double));
  MRISPfunctionValResultForAlpha* 
            fvsForAlphasForPartitions      = (MRISPfunctionValResultForAlpha*)
                                             malloc(forAlphasForPartitionsCapacity*sizeof(MRISPfunctionValResultForAlpha));

  // The old code rotates this sphere is rotated around the x- y- and z- axis and compares this resulting rotated vertices with 
  // the corresponding points in an image, first mapping the points into the image.  The mapping is very expensive.
  // .
  // Rather than do that, this code rotates the vertices around the x and y axes, and then rotates the image around the z axes
  // which avoids a lot of the rotation work
  //
  // The old code converged on a single min, which meant it was very fragile.
  // This code converges on several mins, provided they are far enough apart.
  //
  struct Center { 
    double center_sse;
    int    center_ai, center_bi, center_gi;
    bool   center_sse_known;
  };
  
  // Initialize the centers
  //
  #define centersCapacity 1
  typedef struct Center Centers[centersCapacity];

  Centers outCenters;
  outCenters[0].center_ai = gridSize/2;
  outCenters[0].center_bi = gridSize/2;
  outCenters[0].center_gi = gridSize/2;
  outCenters[0].center_sse = -1.0;  // not known
  outCenters[0].center_sse_known = false;
  int outCentersSize = 1;

  // Walk near the centers, stepping gridSize rather than 1, and going 
  //        center_i + for (j=0 ; j < nangles+1 ; j++) gridStride*(j - nangles/2)       // nangles+1 == forAlphaCapacity because of this
  // 
  int gridStride = (gridSize + nangles - 1)/nangles - 1;
  if (spreadsheet) {
    // format suitable for spreadsheet
    fprintf(stdout, "gridStride %d\n", gridStride); 
  }

  bool changed = true;
  for (;;) {

    // Change to a finer stride when nothing changes
    //  
    if (!changed) {
      if (gridStride == 1) break;
      gridStride /= 2;
      if (spreadsheet) {
        // format suitable for spreadsheet
        fprintf(stdout, "gridStride %d\n", gridStride); 
      }
    }
    changed = false;
    
    // Use the outputs from the previous iteration as the best known centers so far
    //
    Centers inpCenters;
    memcpy(inpCenters, outCenters, outCentersSize*sizeof(struct Center));
    int const inpCentersSize = outCentersSize;
    
    // Search near each to improve them
    //
    int ici;
    for (ici = 0; ici < inpCentersSize; ici++) {
    
      // search around each of the current centers
      //
      struct Center* icp = &inpCenters[ici];
      
      int    const center_ai        = icp->center_ai; 
      int    const center_bi        = icp->center_bi; 
      int    const center_gi        = icp->center_gi;
      double const center_sse       = icp->center_sse;
      bool   const center_sse_known = icp->center_sse_known;

      bool trace = tracing;

      if (tracing) {
        fprintf(stdout, 
          "%s:%d scanning %2.2f degree nbhd of centers[%d] (%2.2f, %2.2f, %2.2f), min sse = %2.2f\n", __FILE__, __LINE__,
          (float)DEGREES(nangles   * radiansPerGridCell), 
          ici,
          (float)DEGREES(iToRadians(center_ai)),  
          (float)DEGREES(iToRadians(center_bi)), 
          (float)DEGREES(iToRadians(center_gi)),
          center_sse_known ? (float)(center_sse + ext_sse) : -666.0);
      }

      int gj;
      for (gj=0; gj < nangles + 1 ; gj++) {
        int   const gi    = center_gi + gridStride*(gj - nangles/2);
        if (gi < 0 || gridSize <= gi) continue;
        float const gamma = iToRadians(gi);
        
        rotateVertices(gammaRotated_xv, gammaRotated_yv, gammaRotated_zv, xv, yv, zv, verticesSize,
          0.0,      // rotate around z axis - last rotation
          0.0,      // rotate around y axis - middle rotation
          gamma);   // rotate around x axis - first rotation

        int bj;
        for (bj=0; bj < nangles + 1 ; bj++) {
          int   const bi   = center_bi + gridStride*(bj - nangles/2);
          if (bi < 0 || gridSize <= bi) continue;
          float const beta = iToRadians(bi);
         
          rotateVertices(betaGammaRotated_xv, betaGammaRotated_yv, betaGammaRotated_zv, gammaRotated_xv, gammaRotated_yv, gammaRotated_zv, verticesSize,
            0.0,    // rotate around z axis - last rotation
            beta,   // rotate around y axis - middle rotation
            0.0);   // rotate around x axis - first rotation

          // select those not already done
          //
          size_t ajsSize = 0;

          int aj;
          for (aj = 0; aj < nangles + 1 ; aj++) {
            int   const ai    = center_ai + gridStride*(aj - nangles/2);
            if (ai < 0 || gridSize <= ai) continue;
            float const alpha = iToRadians(ai);

            int    const doneIndex = gi*gridSize*gridSize + bi*gridSize + ai;
            int    const doneElt   = doneIndex / doneFlagsPerElt;
            size_t const doneFlag  = 1L << (doneIndex&doneFlagsMask);

            bool done = doneFlags[doneElt] & doneFlag;
            if (!done) {
              doneFlags[doneElt] |= doneFlag;
              ajsForAlphas   [ajsSize  ] = aj;
              alphasForAlphas[ajsSize++] = alpha;
            }
          }
              
          // Partition the vertices so can spread across any threads in a deterministic manner
          //
          int const verticesPerPartition = (verticesSize + numberOfVerticesPartitions - 1)/numberOfVerticesPartitions;

          ROMP_PF_begin        
          int partition;
  #ifdef HAVE_OPENMP
          #pragma omp parallel for if_ROMP(assume_reproducible)
  #endif
          for (partition = 0; partition < numberOfVerticesPartitions; partition++) {
            ROMP_PFLB_begin

            int const viLo = partition*verticesPerPartition;
            int const viHi = MIN(verticesSize, (unsigned)viLo + verticesPerPartition);

            double*                         const ssesForAlphas   = &  ssesForAlphasForPartitions[partition*forAlphasCapacity];
            MRISPfunctionValResultForAlpha* const fvsForAlphas    = &   fvsForAlphasForPartitions[partition*forAlphasCapacity];

            { int aj;
              for (aj = 0; aj < forAlphasCapacity; aj++) {
                ssesForAlphas[aj] = 0.0;
              }
            }

            int vi;
            for (vi = viLo; vi < viHi; vi++) {

              // alpha rotates around the z axis
              //      which keeps the z coordinate constant
              //      which can be used to speed up the MRISPfunctionVal_radiusR calls for each vertex
              // hence the movement of this loop inside the vertices loop
              //
              bool const vertexTrace = trace && (vi == 0);

              // process those selected
              //
              MRISPfunctionVal_radiusR(
                  parms->mrisp_template, 
                  fvsForAlphas,             // output values
                  mris->radius, betaGammaRotated_xv[vi], betaGammaRotated_yv[vi], betaGammaRotated_zv[vi], 
                  parms->frame_no, true,
                  alphasForAlphas, ajsSize, // input requests
                  vertexTrace);

              for (unsigned int ajsI = 0; ajsI < ajsSize ; ajsI++) {
                int    const aj     = ajsForAlphas[ajsI];
                double const target = fvsForAlphas[ajsI].curr;
                double const std    = fvsForAlphas[ajsI].next;
                
                int    const ai     = center_ai + gridStride*(aj - nangles/2);
                float  const alpha  = iToRadians(ai);

                if (0) {

                  // This is temporary, until it is built into MRISPfunctionVal_radiusR
                  //
                  float x = 0, y = 0, z = 0;

                  rotateVertices(&x, &y, &z, &betaGammaRotated_xv[vi], &betaGammaRotated_yv[vi], &betaGammaRotated_zv[vi], 1,
                    alpha,    // rotate around z axis - last rotation
                    0.0,      // rotate around y axis - middle rotation
                    0.0);     // rotate around x axis - first rotation

                  bool const vertexAlphaTrace = trace && (vi == 0) && (aj == 0);
                  if (vertexAlphaTrace) {
                    fprintf(stdout, "%s:%d rotated (%g,%g,%g) by (a:%g, b:%g, g:%g) to (%g,%g,%g)\n", __FILE__, __LINE__, 
                      xv[vi],yv[vi],zv[vi],
                      alpha, beta, gamma,
                      x,y,z); 
                  }

                  MRISPfunctionValResultForAlpha targetAndStd;
                  MRISPfunctionVal_radiusR(
                      parms->mrisp_template, &targetAndStd,
                      mris->radius, x, y, z, 
                      parms->frame_no, true,
                      &alpha, 1,
                      vertexAlphaTrace);

                  if (vertexTrace && (aj < 2 || forAlphasCapacity-2 <= aj)) {
                      fprintf(stdout, "%s:%d predicted target:%g target:%g   std:%g v %g\n", __FILE__, __LINE__, 
                        target, targetAndStd.curr,
                        std,    targetAndStd.next);
                  }
                }

                double sqrt_std = sqrt(std);
                if (FZERO(sqrt_std)) {
                  #define DEFAULT_STD 4.0f
                  sqrt_std = DEFAULT_STD /*FSMALL*/;
                }

                float const src = curv[vi];

                double const delta = (src - target) / sqrt_std;
                if (parms->geometry_error) {
                  parms->geometry_error[vno] = (delta * delta);
                }
                if (parms->abs_norm) {
                  ssesForAlphas[aj] += fabs(delta);
                } else {
                  ssesForAlphas[aj] += delta * delta;
                }

              }     // alpha
            }       // vertices

            ROMP_PFLB_end
          }         // partitions
          ROMP_PF_end

          // Combine all the partitions into the 0'th
          //
          {
            int partition;
            for (partition = 1; partition < numberOfVerticesPartitions; partition++) {

              double* const ssesForAlphas = &ssesForAlphasForPartitions[partition*forAlphasCapacity];

              int aj;
              for (aj = 0; aj < forAlphasCapacity; aj++) {
                  ssesForAlphasForPartitions[aj] += ssesForAlphas[aj];
              }
            }
          }

          // Add to the output centers
          //
          for (unsigned int ajsI = 0; ajsI < ajsSize ; ajsI++) {
            int const aj = ajsForAlphas[ajsI];
            int const ai = center_ai + gridStride*(aj - nangles/2);

            double const sse = ssesForAlphasForPartitions[aj];

            if (spreadsheet) {
              int    const ai    = center_ai + gridStride*(aj - nangles/2);
              float  const alpha = iToRadians(ai);
              // format suitable for spreadsheet
              fprintf(stdout, "abgi, %d,%d,%d,  abg, %g,%g,%g, sse,%g\n", ai,bi,gi, alpha,beta,gamma, sse); 
            }

            if (trace) fprintf(stdout, "%s:%d sse:%g after vno:%d\n", __FILE__, __LINE__, sse, vno); 
            trace = false;

            // Consider adding this center to the outCenters, either as a new or as a replacement for an old
            //
            float const radius = (gridStride*nangles) / 3.0f;
                // Anything within this grid index is consider to be 'nearby' and hence replaces the center
            
            int oci;
            for (oci = 0; oci < outCentersSize; oci++) {
              struct Center* oc = &outCenters[oci];
              bool nearBy = 
                  std::abs(oc->center_ai - ai) <= radius
               && std::abs(oc->center_bi - bi) <= radius
               && std::abs(oc->center_gi - gi) <= radius;
              if (nearBy) {
                if (oc->center_sse_known && oc->center_sse < sse) {
                  oci = -1;                                 // Don't replace any
                }
                break;                                      // Was nearBy one,so don't look further
              }
            }
            if (oci == centersCapacity) {                   // A distant local minimum, and centers is full
              oci = 0;                                      // Select the largest existing minimum to be replaced
              int i;
              for (i = 1; i < outCentersSize; i++) {
                if (outCenters[oci].center_sse < outCenters[i].center_sse) oci = i;
              }
              if (outCenters[oci].center_sse <= sse) oci = -1;  // The largest existing min is better than this
            }
            
            // Keep this one
            //
            if (0 <= oci) {
              changed = true;
              struct Center* oc = &outCenters[oci];
              oc->center_ai = ai;
              oc->center_bi = bi;
              oc->center_gi = gi;
              oc->center_sse = sse;
              oc->center_sse_known = true;
              if (oci == outCentersSize) outCentersSize++;
            }
            
          }     // alpha
        }       // beta
      }         // gamma
    }           // centers
    
    if (false || tracing) {
      int oci;
      for (oci = 0; oci < outCentersSize; oci++) {
        struct Center* oc = &outCenters[oci];
        fprintf(stdout, 
          "%s:%d best fit %d at (%2.2f, %2.2f, %2.2f), min sse = %2.2f\n", __FILE__, __LINE__,
          oci, 
          (float)DEGREES(iToRadians(oc->center_ai)),  
          (float)DEGREES(iToRadians(oc->center_bi)),
          (float)DEGREES(iToRadians(oc->center_gi)),
          oc->center_sse_known ? (float)(oc->center_sse + ext_sse) : -666.0);
      }
    }
  }             // gridStride


  // Free the temps
  //
  free(fvsForAlphasForPartitions);  free(ssesForAlphasForPartitions);
  free(alphasForAlphas);            free(ajsForAlphas); 
  free(doneFlags);
  free(betaGammaRotated_zv);        free(betaGammaRotated_yv);          free(betaGammaRotated_xv);
  free(gammaRotated_zv);            free(gammaRotated_yv);              free(gammaRotated_xv);
  free(zv);                         free(yv);                           free(xv);
  free(curv);

  // Return the min center
  //
  { int oci = 0;
    int i;
    for (i = 1; i < outCentersSize; i++) {
      if (outCenters[oci].center_sse > outCenters[i].center_sse) oci = i;
    }
    struct Center* oc = &outCenters[oci];
    *new_mina = iToRadians(oc->center_ai);
    *new_minb = iToRadians(oc->center_bi);
    *new_ming = iToRadians(oc->center_gi);
    *new_sse  = oc->center_sse + ext_sse;
  }
  
#undef iToRadians
}


