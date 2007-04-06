/**
 * @file  ScubaCollectionProperties.h
 * @brief The common view properties
 *
 * An abstract interface implemented by a view. Layers will get
 * a pointer to an object of this type so they can get access to
 * their view settings.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:04 $
 *    $Revision: 1.1 $
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

#ifndef ScubaViewProperties_h
#define ScubaViewProperties_h

class ScubaViewProperties {

 public:

  // Desription:
  // Return the fast mode flag for this layer. If true, the layer
  // should draw as fast as possible and imprecisely if necessary. If
  // false, draw normally. Changes will be broadcast with the
  // "FastModeChanged" message.
  virtual bool GetFastMode () const = 0;
  
  // Desription:
  // Get the value of the Z coordinate (through plane) in 2D
  // mode. Changes will be broadcast with the "Layer2DInfoChanged"
  // message.
  virtual float Get2DRASZ () const = 0;

  // Desription:
  // Get the in plane in 2D mode. Changes will be broadcast with the
  // "Layer2DInfoChanged" message.
  virtual int Get2DInPlane () const = 0;

  // Desription: 
  // Get the value of the coordinates for the planes to display in 3D
  // mode. Changes will be broadcast with the "Layer3DInfoChanged"
  // message. This is only used by 3D layers who display slices of
  // images.
  virtual float Get3DRASX () const = 0;
  virtual float Get3DRASY () const = 0;
  virtual float Get3DRASZ () const = 0;
};

#endif
