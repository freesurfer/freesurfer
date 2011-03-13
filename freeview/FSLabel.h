/**
 * @file  FSLabel.h
 * @brief Base label class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/13 23:04:17 $
 *    $Revision: 1.9 $
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

#ifndef FSLabel_h
#define FSLabel_h

#include <QObject>
#include "vtkImageData.h"
#include "vtkMatrix4x4.h"

extern "C"
{
#include "label.h"
}

class FSVolume;
class wxWindow;
class wxCommandEvent;

class FSLabel : public QObject
{
  Q_OBJECT
public:
  FSLabel( QObject* parent );
  virtual ~FSLabel();

  bool LabelRead( const QString& filename );
  bool LabelWrite( const QString& filename );

  void UpdateLabelFromImage( vtkImageData* rasImage_in, FSVolume* ref_vol );
  void UpdateRASImage( vtkImageData* rasImage_out, FSVolume* ref_vol );

protected:
  LABEL*   m_label;
  bool   m_bTkReg;
};

#endif


