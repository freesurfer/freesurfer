/**
 * @file  vtkKWScubaLayerCollectionODF.h
 * @brief Implementation for ODF viewers.
 *
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

#ifndef vtkKWScubaLayerCollectionODF_h
#define vtkKWScubaLayerCollectionODF_h

#include "vtkKWScubaLayerCollection.h"
#include "ScubaCollectionPropertiesODF.h"


class vtkKWScubaLayerCollectionODF : public vtkKWScubaLayerCollection
                                     //BTX
                                     , public ScubaCollectionPropertiesODF
                                     //ETX
{

 public:

  static vtkKWScubaLayerCollectionODF* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayerCollectionODF, 
      vtkKWScubaLayerCollection );

  // Description:
  // Populate a UI page with common controls for this layer type.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();

  // Description:
  // Set the file name for the ODF volume. We'll use this name to
  // infer the names of the other volumes.
  void SetODFVolumeFileName ( const char* ifnVolume );

  // Settings that are shared by multiple layer types.
  vtkFSVolumeSource* GetODFVolumeSource() const;

 protected:

  vtkKWScubaLayerCollectionODF ();
  ~vtkKWScubaLayerCollectionODF ();

  //BTX
  virtual vtkKWScubaLayer*
    MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode );
  //ETX

  void LoadVolumesFromFileName ();

 private:

  //BTX

  // Common VTK pipeline elements ----------------------------------------
  vtkFSVolumeSource* mODFVolumeSource;
  // ---------------------------------------------------------------------

  std::string mfnODFVolume;
  //ETX

};

#endif
