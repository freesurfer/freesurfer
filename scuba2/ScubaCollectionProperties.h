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
 *    $Date: 2011/03/02 00:04:39 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
