/**
 * @file  WindowTimeCourse.h
 * @brief Tool window to display time course data
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2013/05/03 17:52:39 $
 *    $Revision: 1.2.2.7 $
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

#ifndef WINDOWTIMECOURSE_H
#define WINDOWTIMECOURSE_H

#include <QWidget>

namespace Ui {
    class WindowTimeCourse;
}

class WindowTimeCourse : public QWidget
{
    Q_OBJECT

public:
    explicit WindowTimeCourse(QWidget *parent = 0);
    ~WindowTimeCourse();

public slots:
    void UpdateData();
    void OnFrameChanged(int n);
    void SetCurrentFrame(int n);

signals:
    void FrameChanged(int frame);

private:
    Ui::WindowTimeCourse *ui;
};

#endif // WINDOWTIMECOURSE_H
