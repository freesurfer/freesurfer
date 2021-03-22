/**
 * @brief The common properties
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
  void EmitChangeSignal()
  {
    emit PropertyChanged();
  }

Q_SIGNALS:
  void ShowInfoChanged( bool bShow );
  void PropertyChanged();
  void DisplayModeChanged();

private:
  bool m_bShowInfo;
};

#endif
