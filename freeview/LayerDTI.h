/**
 * @file  LayerDTI.h
 * @brief Layer class for DTI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2010/07/13 20:43:41 $
 *    $Revision: 1.14 $
 *
 * Copyright (C) 2008-2009,
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


