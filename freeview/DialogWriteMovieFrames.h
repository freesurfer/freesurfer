/**
 * @file  DialogWriteMovieFrames.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/13 23:04:17 $
 *    $Revision: 1.6 $
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
 */

#ifndef DIALOGWRITEMOVIEFRAMES_H
#define DIALOGWRITEMOVIEFRAMES_H

#include <QDialog>

namespace Ui
{
class DialogWriteMovieFrames;
}

class DialogWriteMovieFrames : public QDialog
{
  Q_OBJECT

public:
  explicit DialogWriteMovieFrames(QWidget *parent = 0);
  ~DialogWriteMovieFrames();

private:
  Ui::DialogWriteMovieFrames *ui;
};

#endif // DIALOGWRITEMOVIEFRAMES_H
