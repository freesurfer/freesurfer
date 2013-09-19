/**
 * @file  LayerDTI.h
 * @brief Layer class for DTI volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2013/09/19 19:00:50 $
 *    $Revision: 1.20 $
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

  void SetEigenvalueFileName( const QString& filename )
  {
    m_sEigenvalueFileName = filename;
  }

  QString GetEigenvalueFileName()
  {
    return m_sEigenvalueFileName;
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

  FSVolume*  m_eigenvalueSource;   // eigen values
  QString m_sEigenvalueFileName;
};

#endif


