/**
 * @file  QdecVertexAnnotationLookup.h
 * @brief Main QDEC logic
 *
 * A tiny abstract interface class that the view will use to get a
 * annotation string for a vertex number. This is implemented by the
 * window.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/12/13 23:58:14 $
 *    $Revision: 1.1.2.1 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA).
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
