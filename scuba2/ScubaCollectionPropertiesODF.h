/**
 * @file  ScubaCollectionPropertiesODF.h
 * @brief The common properties available to path layers
 *
 * An abstract interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 * 
 * This is only for accessing data.
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author$
 *    $Date$
 *    $Revision$
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

#ifndef ScubaCollectionPropertiesODF_h
#define ScubaCollectionPropertiesODF_h

class vtkFSVolumeSource;

class ScubaCollectionPropertiesODF {

 public:

  // Description:
  // Get a pointer to the source volume.
  virtual vtkFSVolumeSource* GetODFVolumeSource () const = 0;
  
 protected:

  ScubaCollectionPropertiesODF () {};
  virtual ~ScubaCollectionPropertiesODF () {};
};

#endif
