/**
 * @brief The common properties available to surface label
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
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

#ifndef SurfaceLabel_h
#define SurfaceLabel_h

#include <QObject>
#include <vtkSmartPointer.h>
#include <QList>



#include "label.h"


class LayerSurface;
class vtkRGBAColorTransferFunction;
class LayerMRI;

class SurfaceLabel  : public QObject
{
  Q_OBJECT
public:
  SurfaceLabel ( LayerSurface* surf, bool bInitializeLabel = false );
  ~SurfaceLabel ();

  enum ColorCode { SolidColor = 0, Heatscale };
//  void SetSurface( LayerSurface* surf );

  QString GetName();

  void SetName( const QString& name );

  bool LoadLabel( const QString& filename );

  void SetColor( double r, double g, double b );

  double* GetColor()
  {
    return m_rgbColor;
  }

  void MapLabel( unsigned char* colordata, int nVertexCount );

  bool GetShowOutline()
  {
    return m_bShowOutline;
  }

  void SetShowOutline(bool bOutline);

  bool IsVisible()
  {
    return m_bVisible;
  }

  double GetThreshold()
  {
    return m_dThreshold;
  }

  double GetHeatscaleMin()
  {
    return m_dHeatscaleMin;
  }

  double GetHeatscaleMax()
  {
    return m_dHeatscaleMax;
  }

  double GetOpacity()
  {
    return m_dOpacity;
  }

  int GetColorCode()
  {
    return m_nColorCode;
  }

  bool GetCentroid(double* x, double* y, double* z, int* nvo);

  LABEL* GetLabelData()
  {
    return m_label;
  }

  void Resample(LayerMRI* mri);
  void Dilate(int nTimes = 1);
  void Erode(int nTimes = 1);
  void Open(int nTimes = 1);
  void Close(int nTimes = 1);

  QString GetFileName()
  {
    return m_strFilename;
  }

  bool HasVertex(int nvo);

  void EditVertices(const QVector<int>& verts, bool bAdd = true);

  bool SaveToFile(const QString& filename = "");

  bool HasUndo()
  {
    return !m_undoBuffer.isEmpty();
  }

  bool HasRedo()
  {
    return !m_redoBuffer.isEmpty();
  }

Q_SIGNALS:
  void SurfaceLabelChanged();
  void SurfaceLabelVisibilityChanged();

public slots:
  void SetVisible(bool flag);
  void SetThreshold(double th);
  void SetColorCode(int nCode);
  void SetHeatscaleMin(double dval);
  void SetHeatscaleMax(double dval);
  void Undo();
  void Redo();
  void SaveForUndo();
  void SetOpacity(double dval);
  void MaskOverlay();

private:
  void UpdateOutline();
  void UpdateLut();

  LABEL*        m_label;
  QString       m_strName;
  LayerSurface* m_surface;
  double        m_rgbColor[3];
  bool          m_bTkReg;
  bool          m_bShowOutline;
  bool          m_bVisible;
  int*          m_nOutlineIndices;
  double        m_dThreshold;
  int           m_nColorCode;
  double        m_dHeatscaleMin;
  double        m_dHeatscaleMax;
  double        m_dOpacity;
  QString       m_strFilename;
  bool          m_bModified;
  QList<LABEL*> m_undoBuffer;
  QList<LABEL*> m_redoBuffer;

  vtkSmartPointer<vtkRGBAColorTransferFunction> m_lut;
};

#endif
