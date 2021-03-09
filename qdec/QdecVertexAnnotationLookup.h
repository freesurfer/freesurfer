/**
 * @brief Main QDEC logic
 *
 * A tiny abstract interface class that the view will use to get a
 * annotation string for a vertex number. This is implemented by the
 * window.
 */
/*
 * Original Author: Kevin Teich
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

#ifndef QdecVertexAnnotationLookup_h
#define QdecVertexAnnotationLookup_h

class QdecVertexAnnotationLookup {

 public:

  virtual ~QdecVertexAnnotationLookup(){};

  // Get the annotation string for a given vertex number.
  virtual const char* GetAnnotationForVertex ( int inVertex ) = 0;
};

#endif
