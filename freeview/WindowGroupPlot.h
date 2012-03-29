/**
 * @file  WindowGroupPlot.h
 * @brief Tool window to plot group data
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/03/29 20:35:50 $
 *    $Revision: 1.1 $
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

#ifndef WINDOWGROUPPLOT_H
#define WINDOWGROUPPLOT_H

#include <QWidget>

namespace Ui {
    class WindowGroupPlot;
}

class FSGroupDescriptor;

class WindowGroupPlot : public QWidget
{
    Q_OBJECT

public:
    explicit WindowGroupPlot(QWidget *parent = 0);
    ~WindowGroupPlot();

    void SetFsgdData(FSGroupDescriptor* fsgd);

private:
    Ui::WindowGroupPlot *ui;
};

#endif // WINDOWGROUPPLOT_H
