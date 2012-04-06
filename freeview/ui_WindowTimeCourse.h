/********************************************************************************
** Form generated from reading UI file 'WindowTimeCourse.ui'
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
#include <QtGui/QCheckBox>
#include <QtGui/QFrame>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "WidgetTimeCoursePlot.h"

QT_BEGIN_NAMESPACE

class Ui_WindowTimeCourse
{
public:
    QVBoxLayout *verticalLayout_2;
    QFrame *frame;
    QVBoxLayout *verticalLayout;
    WidgetTimeCoursePlot *widgetPlot;
    QHBoxLayout *horizontalLayout;
    QSpacerItem *horizontalSpacer;
    QCheckBox *checkBoxAutoScale;

    void setupUi(QWidget *WindowTimeCourse)
    {
        if (WindowTimeCourse->objectName().isEmpty())
            WindowTimeCourse->setObjectName(QString::fromUtf8("WindowTimeCourse"));
        WindowTimeCourse->resize(400, 300);
        verticalLayout_2 = new QVBoxLayout(WindowTimeCourse);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        frame = new QFrame(WindowTimeCourse);
        frame->setObjectName(QString::fromUtf8("frame"));
        frame->setFrameShape(QFrame::Panel);
        frame->setFrameShadow(QFrame::Sunken);
        verticalLayout = new QVBoxLayout(frame);
        verticalLayout->setContentsMargins(5, 5, 5, 5);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        widgetPlot = new WidgetTimeCoursePlot(frame);
        widgetPlot->setObjectName(QString::fromUtf8("widgetPlot"));
        widgetPlot->setMinimumSize(QSize(350, 250));
        widgetPlot->setFocusPolicy(Qt::StrongFocus);

        verticalLayout->addWidget(widgetPlot);


        verticalLayout_2->addWidget(frame);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalSpacer = new QSpacerItem(5, 5, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);

        checkBoxAutoScale = new QCheckBox(WindowTimeCourse);
        checkBoxAutoScale->setObjectName(QString::fromUtf8("checkBoxAutoScale"));
        checkBoxAutoScale->setChecked(true);

        horizontalLayout->addWidget(checkBoxAutoScale);


        verticalLayout_2->addLayout(horizontalLayout);


        retranslateUi(WindowTimeCourse);
        QObject::connect(checkBoxAutoScale, SIGNAL(toggled(bool)), widgetPlot, SLOT(SetAutoScale(bool)));

        QMetaObject::connectSlotsByName(WindowTimeCourse);
    } // setupUi

    void retranslateUi(QWidget *WindowTimeCourse)
    {
        WindowTimeCourse->setWindowTitle(QApplication::translate("WindowTimeCourse", "Form", 0, QApplication::UnicodeUTF8));
        checkBoxAutoScale->setText(QApplication::translate("WindowTimeCourse", "Auto Scale", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class WindowTimeCourse: public Ui_WindowTimeCourse {};
} // namespace Ui

QT_END_NAMESPACE

