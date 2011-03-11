/**
 * @file  LayerDTI.h
 * @brief Layer class for DTI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/11 23:27:38 $
 *    $Revision: 1.1 $
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

#ifndef LayerDTI_h
#define LayerDTI_h

#include "LayerMRI.h"
#include "vtkSmartPointer.h"
#include <string>

class wxWindow;
class wxCommandEvent;
class LayerPropertiesDTI;

class LayerDTI : public LayerMRI
{
public:
  LayerDTI( LayerMRI* ref );
  virtual ~LayerDTI();

  bool LoadDTIFromFile( wxWindow* wnd, wxCommandEvent& event );

  void SetVectorFileName( const char* filename )
  {
    m_sVectorFileName = filename;
  }

  const char* GetVectorFileName()
  {
    return m_sVectorFileName.c_str();
  }

  inline LayerPropertiesDTI* GetProperties()
  {
    return (LayerPropertiesDTI*)mProperties;
  }

  bool GetVectorValue( double* pos_in, double* v_out );

protected:
  bool DoRotate( std::vector<RotationElement>& rotations, wxWindow* wnd, wxCommandEvent& event );
  void DoRestore();
  void UpdateColorMap();
  void InitializeDTIColorMap( wxWindow* wnd, wxCommandEvent& event );
  
  virtual void UpdateVectorActor( int nPlane );

  FSVolume*  m_vectorSource;
  std::string  m_sVectorFileName;
};

#endif


