/**
 * @brief Data handler for surface annotation
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

#ifndef SurfaceAnnotation_h
#define SurfaceAnnotation_h

#include <QObject>
#include <vtkSmartPointer.h>
#include <QMap>
#include <QColor>
#include <QVector>

#include "colortab.h"

#ifndef UNASSIGNED_ANNOT_BASE
#define UNASSIGNED_ANNOT_BASE 1000000
#endif

class vtkLookupTable;
class vtkRGBAColorTransferFunction;
class LayerSurface;
class vtkActor;
class vtkPolyData;

struct NewAnnotationLabel
{
  QString name;
  QColor color;
  int id;
};

struct AnnotUndoRedoBufferItem
{
  AnnotUndoRedoBufferItem()
  {
    m_data = NULL;
    m_ctab = NULL;
  }

  void Free()
  {
    if (m_data)
      delete[] m_data;
    m_data = NULL;
    CTABfree(&m_ctab);
  }

  int*          m_data;
  QList<int>    m_listVisibleLabels;
  QMap<int, NewAnnotationLabel> m_mapNewLabels;
  COLOR_TABLE*  m_ctab;
};

class SurfaceAnnotation  : public QObject
{
  Q_OBJECT
public:
  SurfaceAnnotation ( LayerSurface* surf );
  ~SurfaceAnnotation ();

  void SetSurface( LayerSurface* surf );

  QString GetName();

  void SetName( const QString& name );

  bool LoadAnnotation( const QString& fn);

  bool LoadFromSegmentation( const QString& fn);

  bool LoadColorTable( const QString& fn );

  bool InitializeNewAnnotation(const QString& ctab_fn);

  int* GetIndices()
  {
    return (m_bShowOutline ? m_nOutlineIndices : m_nIndices);
  }

  int GetIndexSize()
  {
    return m_nIndexSize;
  }

  COLOR_TABLE* GetColorTable()
  {
    return m_lut;
  }

  int GetIndexAtVertex( int nVertex );

  void GetAnnotationPoint( int nIndex, double* pt_out );

  QString GetAnnotationNameAtIndex( int nIndex );

  QString GetAnnotationNameAtVertex( int nVertex );

  void GetAnnotationColorAtIndex( int nIndex, int* rgb );

  bool GetShowOutline()
  {
    return m_bShowOutline;
  }

  void SetShowOutline(bool bOutline);

  void MapAnnotationColor( unsigned char* colordata );

  QString GetFilename()
  {
    return m_strFilename;
  }

  void SetFilename(const QString& fn)
  {
    m_strFilename = fn;
  }

  QList<int> GetVisibleLabels()
  {
    return m_listVisibleLabels;
  }

  void SetSelectLabel(int nVal, bool bSelected);
  void SetSelectAllLabels();
  void SetUnselectAllLabels();

  void SetHighlightedLabel(int n);

  int GetHighlightedLabel()
  {
    return m_nHighlightedLabel;
  }

  int* GetAnnotationData()
  {
    return m_data;
  }

  QList<int> GetExistingAnnotations()
  {
    return m_listAnnotations;
  }

  void EditLabel(const QVector<int>& verts, int fill_index, const QVariantMap& options);

  void ReassignNewLabel(int old_id, int new_id, const QString& name = "", const QColor& color = QColor());

  QMap<int, NewAnnotationLabel> GetNewLabels()
  {
    return m_mapNewLabels;
  }

  void UpdateLabelInfo(int i, const QString& name, const QColor& color = QColor());

  bool HasColor(const QColor& color)
  {
    return m_listColors.contains(color);
  }

  bool HasUnassignedLabels()
  {
    return !m_mapNewLabels.isEmpty();
  }

  bool SaveToFile(const QString& fn);

  bool HasUndo();
  bool HasRedo();

  void DeleteLabel(int nIndex);

  void CleanUpColorTable();

signals:
  void Modified();

public slots:
  void SetModified()
  {
    emit Modified();
  }

  void Undo();
  void Redo();
  void SaveForUndo();

protected:
  void Reset();
  void UpdateData();

private:
  QColor GenerateNewColor();
  void UpdateColorList();
  int ColorToAnnotation(const QColor& c);
  AnnotUndoRedoBufferItem SaveCurrentUndoRedoBuffer();
  void RestoreFromUndoRedoBuffer(const AnnotUndoRedoBufferItem& item);
  void UpdateColorTable(int nIndex, const QString& name, const QColor& color);

  int*          m_nIndices;
  int*          m_nOutlineIndices;
  int           m_nIndexSize;
  int*          m_nCenterVertices;  // center vertex of each annotation
  int*          m_data;
  QList<int>    m_listAnnotations;

  QString       m_strName;
  COLOR_TABLE*  m_lut;
  LayerSurface* m_surface;
  bool          m_bShowOutline;
  double        m_dOpacity;
  QString       m_strFilename;
  QList<int>    m_listVisibleLabels;
  int           m_nHighlightedLabel;
  QMap<int, NewAnnotationLabel> m_mapNewLabels;
  QVector<QColor> m_listColors;

  QVector<AnnotUndoRedoBufferItem>  m_bufferUndo;
  QVector<AnnotUndoRedoBufferItem>  m_bufferRedo;
};

#endif
