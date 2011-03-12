/**
 * @file  LayerDTI.h
 * @brief Layer class for DTI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: krish $
 *    $Date: 2011/03/12 00:28:49 $
 *    $Revision: 1.16 $
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

class LayerPropertyDTI;

class LayerDTI : public LayerMRI
{
  Q_OBJECT
public:
  LayerDTI( LayerMRI* ref, QObject* parent = NULL );
  virtual ~LayerDTI();

  bool LoadDTIFromFile();

  void SetVectorFileName( const QString& filename )
  {
    m_sVectorFileName = filename;
  }

  QString GetVectorFileName()
  {
    return m_sVectorFileName;
  }

  inline LayerPropertyDTI* GetProperty()
  {
    return (LayerPropertyDTI*)mProperty;
  }

  bool GetVectorValue( double* pos_in, double* v_out );

protected:
  bool DoRotate( std::vector<RotationElement>& rotations );
  void DoRestore();
  void UpdateColorMap();
  void InitializeDTIColorMap();
  
  virtual void UpdateVectorActor( int nPlane );

  FSVolume*  m_vectorSource;
  QString  m_sVectorFileName;
};

#endif


