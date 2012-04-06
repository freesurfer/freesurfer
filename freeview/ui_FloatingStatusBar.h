/********************************************************************************
** Form generated from reading UI file 'FloatingStatusBar.ui'
**
** Created: Fri Apr 6 15:12:22 2012
**      by: Qt User Interface Compiler version 4.7.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QFrame>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QProgressBar>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_FloatingStatusBar
{
public:
    QHBoxLayout *horizontalLayout;
    QFrame *frame;
    QHBoxLayout *horizontalLayout_2;
    QProgressBar *progressBar;

    void setupUi(QWidget *FloatingStatusBar)
    {
        if (FloatingStatusBar->objectName().isEmpty())
            FloatingStatusBar->setObjectName(QString::fromUtf8("FloatingStatusBar"));
        FloatingStatusBar->resize(387, 25);
        FloatingStatusBar->setMinimumSize(QSize(200, 0));
        FloatingStatusBar->setFocusPolicy(Qt::ClickFocus);
        horizontalLayout = new QHBoxLayout(FloatingStatusBar);
        horizontalLayout->setSpacing(0);
        horizontalLayout->setContentsMargins(0, 0, 0, 0);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setSizeConstraint(QLayout::SetFixedSize);
        frame = new QFrame(FloatingStatusBar);
        frame->setObjectName(QString::fromUtf8("frame"));
        frame->setFrameShape(QFrame::NoFrame);
        frame->setFrameShadow(QFrame::Raised);
        horizontalLayout_2 = new QHBoxLayout(frame);
        horizontalLayout_2->setSpacing(5);
        horizontalLayout_2->setContentsMargins(5, 5, 5, 5);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalLayout_2->setSizeConstraint(QLayout::SetFixedSize);
        progressBar = new QProgressBar(frame);
        progressBar->setObjectName(QString::fromUtf8("progressBar"));
        progressBar->setMinimumSize(QSize(170, 0));
        progressBar->setMaximumSize(QSize(16777215, 15));
        progressBar->setFocusPolicy(Qt::StrongFocus);
        progressBar->setValue(24);

        horizontalLayout_2->addWidget(progressBar);


        horizontalLayout->addWidget(frame);


        retranslateUi(FloatingStatusBar);

        QMetaObject::connectSlotsByName(FloatingStatusBar);
    } // setupUi

    void retranslateUi(QWidget *FloatingStatusBar)
    {
        FloatingStatusBar->setWindowTitle(QString());
    } // retranslateUi

};

namespace Ui {
    class FloatingStatusBar: public Ui_FloatingStatusBar {};
} // namespace Ui

QT_END_NAMESPACE

