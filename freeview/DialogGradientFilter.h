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
 */
#ifndef DIALOGGRADIENTFILTER_H
#define DIALOGGRADIENTFILTER_H

#include <QDialog>

namespace Ui
{
class DialogGradientFilter;
}

class DialogGradientFilter : public QDialog
{
  Q_OBJECT

public:
  explicit DialogGradientFilter(QWidget *parent = 0);
  ~DialogGradientFilter();

  void SetSmoothing(bool smooth);
  bool GetSmoothing();

  void SetSD(double val);
  double GetSD();

private:
  Ui::DialogGradientFilter *ui;
};

#endif // DIALOGGRADIENTFILTER_H
