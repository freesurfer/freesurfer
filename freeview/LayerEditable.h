/**
 * @brief Base Layer class for editable volume.
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

#ifndef LayerEditable_h
#define LayerEditable_h

#include "Layer.h"
#include <QString>

class LayerEditable : public Layer
{
  Q_OBJECT
public:
  LayerEditable( QObject* parent = NULL );
  virtual ~LayerEditable();

  virtual bool HasUndo()
  {
    return false;
  }
  virtual bool HasRedo()
  {
    return false;
  }

  virtual void Undo()
  {}
  virtual void Redo()
  {}

  void ResetModified()
  {
    m_bModified = false;
  }

  bool IsModified()
  {
    return m_bModified;
  }

  void SetEditable( bool bEditable = true )
  {
    m_bEditable = bEditable;
  }

  bool IsEditable()
  {
    return m_bEditable;
  }

  QString GetRegFileName()
  {
    return m_sRegFilename;
  }

  void SetRegFileName( const QString& fn )
  {
    m_sRegFilename = fn;
  }


Q_SIGNALS:
  void Modified();

public slots:
  virtual void SetModified();

protected:

  int   m_nMaxUndoSteps;
  bool  m_bModified;

  bool   m_bEditable;

  QString m_sRegFilename;
};

#endif


