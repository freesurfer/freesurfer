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
