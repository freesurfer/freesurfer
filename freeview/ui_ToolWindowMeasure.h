/********************************************************************************
** Form generated from reading UI file 'ToolWindowMeasure.ui'
**
** Created: Fri Apr 6 15:12:21 2012
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
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QTextBrowser>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "qtcolorpicker.h"

QT_BEGIN_NAMESPACE

class Ui_ToolWindowMeasure
{
public:
    QAction *actionLine;
    QAction *actionPolyLine;
    QAction *actionSpline;
    QAction *actionRectangle;
    QAction *actionContour;
    QAction *actionLabel;
    QVBoxLayout *verticalLayout_main;
    QToolBar *toolBar;
    QFrame *lineSeparator;
    QHBoxLayout *horizontalLayout;
    QLabel *labelId;
    QSpinBox *spinBoxId;
    QLabel *labelGroup;
    QSpinBox *spinBoxGroup;
    QtColorPicker *colorPickerGroup;
    QTextBrowser *textBrowserInfo;
    QHBoxLayout *horizontalLayout_3;
    QSpacerItem *horizontalSpacer_3;
    QPushButton *pushButtonLoad;
    QPushButton *pushButtonSave;
    QPushButton *pushButtonSaveAll;
    QSpacerItem *horizontalSpacer_4;
    QHBoxLayout *horizontalLayout_2;
    QSpacerItem *horizontalSpacer_2;
    QPushButton *pushButtonUpdate;
    QPushButton *pushButtonCopy;
    QPushButton *pushButtonExport;
    QSpacerItem *horizontalSpacer;

    void setupUi(QWidget *ToolWindowMeasure)
    {
        if (ToolWindowMeasure->objectName().isEmpty())
            ToolWindowMeasure->setObjectName(QString::fromUtf8("ToolWindowMeasure"));
        ToolWindowMeasure->resize(285, 400);
        actionLine = new QAction(ToolWindowMeasure);
        actionLine->setObjectName(QString::fromUtf8("actionLine"));
        actionLine->setCheckable(true);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/measure_line.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLine->setIcon(icon);
        actionPolyLine = new QAction(ToolWindowMeasure);
        actionPolyLine->setObjectName(QString::fromUtf8("actionPolyLine"));
        actionPolyLine->setCheckable(true);
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/resource/icons/measure_polyline.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionPolyLine->setIcon(icon1);
        actionSpline = new QAction(ToolWindowMeasure);
        actionSpline->setObjectName(QString::fromUtf8("actionSpline"));
        actionSpline->setCheckable(true);
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/resource/icons/measure_spline.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSpline->setIcon(icon2);
        actionRectangle = new QAction(ToolWindowMeasure);
        actionRectangle->setObjectName(QString::fromUtf8("actionRectangle"));
        actionRectangle->setCheckable(true);
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/resource/icons/measure_rectangle.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionRectangle->setIcon(icon3);
        actionContour = new QAction(ToolWindowMeasure);
        actionContour->setObjectName(QString::fromUtf8("actionContour"));
        actionContour->setCheckable(true);
        QIcon icon4;
        icon4.addFile(QString::fromUtf8(":/resource/icons/measure_surface_region.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionContour->setIcon(icon4);
        actionLabel = new QAction(ToolWindowMeasure);
        actionLabel->setObjectName(QString::fromUtf8("actionLabel"));
        actionLabel->setCheckable(true);
        QIcon icon5;
        icon5.addFile(QString::fromUtf8(":/resource/icons/measure_label.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLabel->setIcon(icon5);
        verticalLayout_main = new QVBoxLayout(ToolWindowMeasure);
        verticalLayout_main->setSpacing(0);
        verticalLayout_main->setContentsMargins(0, 0, 0, 0);
        verticalLayout_main->setObjectName(QString::fromUtf8("verticalLayout_main"));
        verticalLayout_main->setSizeConstraint(QLayout::SetMinimumSize);
        toolBar = new QToolBar(ToolWindowMeasure);
        toolBar->setObjectName(QString::fromUtf8("toolBar"));
        toolBar->setIconSize(QSize(24, 24));

        verticalLayout_main->addWidget(toolBar);

        lineSeparator = new QFrame(ToolWindowMeasure);
        lineSeparator->setObjectName(QString::fromUtf8("lineSeparator"));
        lineSeparator->setFrameShape(QFrame::HLine);
        lineSeparator->setFrameShadow(QFrame::Sunken);

        verticalLayout_main->addWidget(lineSeparator);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalLayout->setContentsMargins(5, 0, 5, 0);
        labelId = new QLabel(ToolWindowMeasure);
        labelId->setObjectName(QString::fromUtf8("labelId"));

        horizontalLayout->addWidget(labelId);

        spinBoxId = new QSpinBox(ToolWindowMeasure);
        spinBoxId->setObjectName(QString::fromUtf8("spinBoxId"));

        horizontalLayout->addWidget(spinBoxId);

        labelGroup = new QLabel(ToolWindowMeasure);
        labelGroup->setObjectName(QString::fromUtf8("labelGroup"));

        horizontalLayout->addWidget(labelGroup);

        spinBoxGroup = new QSpinBox(ToolWindowMeasure);
        spinBoxGroup->setObjectName(QString::fromUtf8("spinBoxGroup"));

        horizontalLayout->addWidget(spinBoxGroup);

        colorPickerGroup = new QtColorPicker(ToolWindowMeasure);
        colorPickerGroup->setObjectName(QString::fromUtf8("colorPickerGroup"));
        colorPickerGroup->setMinimumSize(QSize(40, 0));

        horizontalLayout->addWidget(colorPickerGroup);

        horizontalLayout->setStretch(4, 1);

        verticalLayout_main->addLayout(horizontalLayout);

        textBrowserInfo = new QTextBrowser(ToolWindowMeasure);
        textBrowserInfo->setObjectName(QString::fromUtf8("textBrowserInfo"));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(1);
        sizePolicy.setHeightForWidth(textBrowserInfo->sizePolicy().hasHeightForWidth());
        textBrowserInfo->setSizePolicy(sizePolicy);
        textBrowserInfo->setMinimumSize(QSize(0, 270));
        textBrowserInfo->setTabStopWidth(50);

        verticalLayout_main->addWidget(textBrowserInfo);

        horizontalLayout_3 = new QHBoxLayout();
#ifndef Q_OS_MAC
        horizontalLayout_3->setSpacing(6);
#endif
        horizontalLayout_3->setContentsMargins(0, 0, 0, 0);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        horizontalSpacer_3 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_3);

        pushButtonLoad = new QPushButton(ToolWindowMeasure);
        pushButtonLoad->setObjectName(QString::fromUtf8("pushButtonLoad"));
        pushButtonLoad->setIconSize(QSize(16, 16));

        horizontalLayout_3->addWidget(pushButtonLoad);

        pushButtonSave = new QPushButton(ToolWindowMeasure);
        pushButtonSave->setObjectName(QString::fromUtf8("pushButtonSave"));

        horizontalLayout_3->addWidget(pushButtonSave);

        pushButtonSaveAll = new QPushButton(ToolWindowMeasure);
        pushButtonSaveAll->setObjectName(QString::fromUtf8("pushButtonSaveAll"));

        horizontalLayout_3->addWidget(pushButtonSaveAll);

        horizontalSpacer_4 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_4);


        verticalLayout_main->addLayout(horizontalLayout_3);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalSpacer_2 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_2);

        pushButtonUpdate = new QPushButton(ToolWindowMeasure);
        pushButtonUpdate->setObjectName(QString::fromUtf8("pushButtonUpdate"));

        horizontalLayout_2->addWidget(pushButtonUpdate);

        pushButtonCopy = new QPushButton(ToolWindowMeasure);
        pushButtonCopy->setObjectName(QString::fromUtf8("pushButtonCopy"));
        pushButtonCopy->setMinimumSize(QSize(0, 0));

        horizontalLayout_2->addWidget(pushButtonCopy);

        pushButtonExport = new QPushButton(ToolWindowMeasure);
        pushButtonExport->setObjectName(QString::fromUtf8("pushButtonExport"));

        horizontalLayout_2->addWidget(pushButtonExport);

        horizontalSpacer = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer);


        verticalLayout_main->addLayout(horizontalLayout_2);

        verticalLayout_main->setStretch(3, 1);

        toolBar->addAction(actionLine);
        toolBar->addAction(actionPolyLine);
        toolBar->addAction(actionSpline);
        toolBar->addAction(actionRectangle);
        toolBar->addAction(actionLabel);
        toolBar->addAction(actionContour);

        retranslateUi(ToolWindowMeasure);
        QObject::connect(spinBoxId, SIGNAL(valueChanged(int)), ToolWindowMeasure, SLOT(OnSpinBoxId(int)));
        QObject::connect(spinBoxGroup, SIGNAL(valueChanged(int)), ToolWindowMeasure, SLOT(OnSpinBoxGroup(int)));
        QObject::connect(pushButtonLoad, SIGNAL(clicked()), ToolWindowMeasure, SLOT(OnLoad()));
        QObject::connect(pushButtonSave, SIGNAL(clicked()), ToolWindowMeasure, SLOT(OnSave()));
        QObject::connect(pushButtonSaveAll, SIGNAL(clicked()), ToolWindowMeasure, SLOT(OnSaveAll()));
        QObject::connect(pushButtonUpdate, SIGNAL(clicked()), ToolWindowMeasure, SLOT(OnUpdate()));
        QObject::connect(pushButtonCopy, SIGNAL(clicked()), ToolWindowMeasure, SLOT(OnCopy()));
        QObject::connect(pushButtonExport, SIGNAL(clicked()), ToolWindowMeasure, SLOT(OnExport()));
        QObject::connect(colorPickerGroup, SIGNAL(colorChanged(QColor)), ToolWindowMeasure, SLOT(OnColorGroup(QColor)));

        QMetaObject::connectSlotsByName(ToolWindowMeasure);
    } // setupUi

    void retranslateUi(QWidget *ToolWindowMeasure)
    {
        ToolWindowMeasure->setWindowTitle(QApplication::translate("ToolWindowMeasure", "Measure", 0, QApplication::UnicodeUTF8));
        actionLine->setText(QApplication::translate("ToolWindowMeasure", "Line", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionLine->setToolTip(QApplication::translate("ToolWindowMeasure", "Line", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionPolyLine->setText(QApplication::translate("ToolWindowMeasure", "PolyLine", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionPolyLine->setToolTip(QApplication::translate("ToolWindowMeasure", "Polyline", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionSpline->setText(QApplication::translate("ToolWindowMeasure", "Spline", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionSpline->setToolTip(QApplication::translate("ToolWindowMeasure", "Spline", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionRectangle->setText(QApplication::translate("ToolWindowMeasure", "Rectangle", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionRectangle->setToolTip(QApplication::translate("ToolWindowMeasure", "Rectangle", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionContour->setText(QApplication::translate("ToolWindowMeasure", "Contour", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionContour->setToolTip(QApplication::translate("ToolWindowMeasure", "IsoSurface Region", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionLabel->setText(QApplication::translate("ToolWindowMeasure", "Label", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionLabel->setToolTip(QApplication::translate("ToolWindowMeasure", "Label", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        labelId->setText(QApplication::translate("ToolWindowMeasure", "Id", 0, QApplication::UnicodeUTF8));
        labelGroup->setText(QApplication::translate("ToolWindowMeasure", "Group", 0, QApplication::UnicodeUTF8));
        pushButtonLoad->setText(QApplication::translate("ToolWindowMeasure", "Load", 0, QApplication::UnicodeUTF8));
        pushButtonSave->setText(QApplication::translate("ToolWindowMeasure", "Save", 0, QApplication::UnicodeUTF8));
        pushButtonSaveAll->setText(QApplication::translate("ToolWindowMeasure", "Save All", 0, QApplication::UnicodeUTF8));
        pushButtonUpdate->setText(QApplication::translate("ToolWindowMeasure", "Update", 0, QApplication::UnicodeUTF8));
        pushButtonCopy->setText(QApplication::translate("ToolWindowMeasure", "Copy", 0, QApplication::UnicodeUTF8));
        pushButtonExport->setText(QApplication::translate("ToolWindowMeasure", "Export", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class ToolWindowMeasure: public Ui_ToolWindowMeasure {};
} // namespace Ui

QT_END_NAMESPACE

