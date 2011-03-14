/**
 * @file  LayerEditable.h
 * @brief Base Layer class for editable volume.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:58 $
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

  virtual void SetModified();
  
Q_SIGNALS:
  void Modified();

protected:

  int   m_nMaxUndoSteps;
  bool  m_bModified;

  bool   m_bEditable;

  QString m_sRegFilename;
};

#endif


