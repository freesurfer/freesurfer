/**
 * @file  ScubaCollectionProperties.h
 * @brief The common properties available to layers
 *
 * An abstract interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/04/09 22:38:09 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2007,
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

#ifndef ScubaCollectionProperties_h
#define ScubaCollectionProperties_h

class ScubaCollectionProperties {

 public:

  // Desription:
  // Return a human-readable label for this collection to be used in
  // the GUI.
  virtual const char* GetLabel () const = 0;
  
  // Desription:
  // Return the opacity for this layer from 0-1. Changes will be
  // broadcast with the "OpacityChanged" message.
  virtual float GetOpacity () const = 0;

 protected:

  ScubaCollectionProperties () {};
  virtual ~ScubaCollectionProperties () {};
};

#endif
