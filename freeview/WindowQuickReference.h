/**
 * @file  WindowQuickReference.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/13 23:04:18 $
 *    $Revision: 1.8 $
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

#ifndef WINDOWQUICKREFERENCE_H
#define WINDOWQUICKREFERENCE_H

#include <QWidget>

namespace Ui
{
class WindowQuickReference;
}

class WindowQuickReference : public QWidget
{
  Q_OBJECT

public:
  explicit WindowQuickReference(QWidget *parent = 0);
  ~WindowQuickReference();

private:
  Ui::WindowQuickReference *ui;
};

#endif // WINDOWQUICKREFERENCE_H
