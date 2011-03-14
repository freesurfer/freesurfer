/**
 * @file  LayerProperty.h
 * @brief The common properties
 *
 * An interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2011/03/14 21:20:58 $
 *    $Revision: 1.3 $
 *
 * Copyright (C) 2007-2009,
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

#ifndef LayerProperty_h
#define LayerProperty_h

#include <QObject>

class LayerProperty : public QObject
{
    Q_OBJECT
public:
  LayerProperty ( QObject* parent = 0 );
  ~LayerProperty ();

  bool GetShowInfo()
  {
      return m_bShowInfo;
  }

public slots:
  void SetShowInfo( bool bShow );

Q_SIGNALS:
  void ShowInfoChanged( bool bShow );
  void PropertyChanged();
  void DisplayModeChanged();

private:
  bool m_bShowInfo;
};

#endif
