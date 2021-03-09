/**
 * @brief Base label class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <QList>



#include "label.h"
#undef uchar  // conflicts with Qt
#include "mrisutils.h"


class FSVolume;
class FSSurface;

class FSLabel : public QObject
{
  Q_OBJECT
public:
  FSLabel( QObject* parent, FSVolume* mri_template );
  virtual ~FSLabel();

  bool LabelRead( const QString& filename );
  bool LabelWrite( const QString& filename );

  void Initialize(FSVolume* ref_vol, FSSurface* surf, int coords);

  void UpdateLabelFromImage( vtkImageData* rasImage_in, FSVolume* ref_vol );
  void UpdateRASImage( vtkImageData* rasImage_out, FSVolume* ref_vol, double threshold = -1e10 );
  void FillUnassignedVertices(FSSurface* surf, FSVolume* mri_template, int coords);
  void EditVoxel(int nx, int ny, int nz, int coords, bool bAdd, int* vertices = NULL, int* pnum = NULL);

  bool GetCentroidRASPosition(double* pos, FSVolume* ref_vol);

  void GetStatsRange(double* range);

  bool UpdateStatsRange(double val);

  LABEL* GetRawLabel()
  {
    return m_label;
  }

  bool HasUndo();
  bool HasRedo();

  void Undo();
  void Redo();
  void SaveForUndo();
  void Clear();

protected:
  LABEL*   m_label;
  QList<LABEL*> m_undoBuffer;
  QList<LABEL*> m_redoBuffer;
  double   m_dStatsRange[2];
  LABEL2SURF* m_l2s;
  FSVolume*   m_mri_template;
};

#endif


