/********************************************************************************
** Form generated from reading UI file 'ToolWindowEdit.ui'
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
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "qtcolorpicker.h"

QT_BEGIN_NAMESPACE

class Ui_ToolWindowEdit
{
public:
    QAction *actionFreeHand;
    QAction *actionPolyLine;
    QAction *actionLiveWire;
    QAction *actionFill;
    QAction *actionContour;
    QAction *actionColorPicker;
    QAction *actionReplaceLabel;
    QVBoxLayout *verticalLayoutMain;
    QToolBar *toolBar;
    QFrame *line;
    QVBoxLayout *verticalLayout;
    QGridLayout *gridLayout;
    QLabel *labelBrushSize;
    QSpinBox *spinBoxBrushSize;
    QSpacerItem *horizontalSpacer;
    QLabel *labelReference;
    QComboBox *comboBoxReference;
    QLabel *labelTolerance;
    QLabel *labelContourValue;
    QLineEdit *lineEditContourValue;
    QSpinBox *spinBoxTolerance;
    QSpacerItem *horizontalSpacer_2;
    QSpacerItem *horizontalSpacer_3;
    QCheckBox *checkBoxConstrain;
    QCheckBox *checkBoxDrawRange;
    QHBoxLayout *horizontalLayout;
    QLabel *labelDrawRangeLow;
    QLineEdit *lineEditDrawRangeLow;
    QSpacerItem *horizontalSpacer_8;
    QLabel *labelDrawRangeHigh;
    QLineEdit *lineEditDrawRangeHigh;
    QSpacerItem *horizontalSpacer_5;
    QCheckBox *checkBoxExcludeRange;
    QHBoxLayout *horizontalLayout_2;
    QLabel *labelExcludeRangeLow;
    QLineEdit *lineEditExcludeRangeLow;
    QSpacerItem *horizontalSpacer_9;
    QLabel *labelExcludeRangeHigh;
    QLineEdit *lineEditExcludeRangeHigh;
    QSpacerItem *horizontalSpacer_6;
    QCheckBox *checkBoxSmooth;
    QHBoxLayout *horizontalLayout_3;
    QLabel *labelSD;
    QLineEdit *lineEditSD;
    QSpacerItem *horizontalSpacer_7;
    QHBoxLayout *horizontalLayout_4;
    QLabel *labelContourColor;
    QtColorPicker *colorPickerContour;
    QSpacerItem *horizontalSpacer_4;
    QCheckBox *checkBoxFill3D;
    QSpacerItem *verticalSpacer_2;
    QLabel *labelTips;
    QSpacerItem *verticalSpacer;

    void setupUi(QWidget *ToolWindowEdit)
    {
        if (ToolWindowEdit->objectName().isEmpty())
            ToolWindowEdit->setObjectName(QString::fromUtf8("ToolWindowEdit"));
        ToolWindowEdit->resize(326, 629);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(ToolWindowEdit->sizePolicy().hasHeightForWidth());
        ToolWindowEdit->setSizePolicy(sizePolicy);
        actionFreeHand = new QAction(ToolWindowEdit);
        actionFreeHand->setObjectName(QString::fromUtf8("actionFreeHand"));
        actionFreeHand->setCheckable(true);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/voxel_draw_freehand.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionFreeHand->setIcon(icon);
        actionPolyLine = new QAction(ToolWindowEdit);
        actionPolyLine->setObjectName(QString::fromUtf8("actionPolyLine"));
        actionPolyLine->setCheckable(true);
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/resource/icons/voxel_draw_polyline.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionPolyLine->setIcon(icon1);
        actionLiveWire = new QAction(ToolWindowEdit);
        actionLiveWire->setObjectName(QString::fromUtf8("actionLiveWire"));
        actionLiveWire->setCheckable(true);
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/resource/icons/voxel_draw_livewire.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLiveWire->setIcon(icon2);
        actionFill = new QAction(ToolWindowEdit);
        actionFill->setObjectName(QString::fromUtf8("actionFill"));
        actionFill->setCheckable(true);
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/resource/icons/voxel_draw_fill.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionFill->setIcon(icon3);
        actionContour = new QAction(ToolWindowEdit);
        actionContour->setObjectName(QString::fromUtf8("actionContour"));
        actionContour->setCheckable(true);
        QIcon icon4;
        icon4.addFile(QString::fromUtf8(":/resource/icons/voxel_draw_contour.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionContour->setIcon(icon4);
        actionColorPicker = new QAction(ToolWindowEdit);
        actionColorPicker->setObjectName(QString::fromUtf8("actionColorPicker"));
        actionColorPicker->setCheckable(true);
        QIcon icon5;
        icon5.addFile(QString::fromUtf8(":/resource/icons/voxel_draw_colorpicker.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionColorPicker->setIcon(icon5);
        actionReplaceLabel = new QAction(ToolWindowEdit);
        actionReplaceLabel->setObjectName(QString::fromUtf8("actionReplaceLabel"));
        QIcon icon6;
        icon6.addFile(QString::fromUtf8(":/resource/icons/voxel_draw_replace_color.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionReplaceLabel->setIcon(icon6);
        verticalLayoutMain = new QVBoxLayout(ToolWindowEdit);
        verticalLayoutMain->setSpacing(0);
        verticalLayoutMain->setContentsMargins(0, 0, 0, 0);
        verticalLayoutMain->setObjectName(QString::fromUtf8("verticalLayoutMain"));
        verticalLayoutMain->setSizeConstraint(QLayout::SetFixedSize);
        toolBar = new QToolBar(ToolWindowEdit);
        toolBar->setObjectName(QString::fromUtf8("toolBar"));
        toolBar->setIconSize(QSize(24, 24));

        verticalLayoutMain->addWidget(toolBar);

        line = new QFrame(ToolWindowEdit);
        line->setObjectName(QString::fromUtf8("line"));
        line->setFrameShape(QFrame::HLine);
        line->setFrameShadow(QFrame::Sunken);

        verticalLayoutMain->addWidget(line);

        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(7);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        verticalLayout->setContentsMargins(9, 5, 9, 5);
        gridLayout = new QGridLayout();
        gridLayout->setSpacing(6);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        labelBrushSize = new QLabel(ToolWindowEdit);
        labelBrushSize->setObjectName(QString::fromUtf8("labelBrushSize"));

        gridLayout->addWidget(labelBrushSize, 0, 0, 1, 1);

        spinBoxBrushSize = new QSpinBox(ToolWindowEdit);
        spinBoxBrushSize->setObjectName(QString::fromUtf8("spinBoxBrushSize"));
        spinBoxBrushSize->setMinimum(1);
        spinBoxBrushSize->setMaximum(512);

        gridLayout->addWidget(spinBoxBrushSize, 0, 1, 1, 1);

        horizontalSpacer = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer, 0, 2, 1, 3);

        labelReference = new QLabel(ToolWindowEdit);
        labelReference->setObjectName(QString::fromUtf8("labelReference"));

        gridLayout->addWidget(labelReference, 1, 0, 1, 1);

        comboBoxReference = new QComboBox(ToolWindowEdit);
        comboBoxReference->setObjectName(QString::fromUtf8("comboBoxReference"));
        QSizePolicy sizePolicy1(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(1);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(comboBoxReference->sizePolicy().hasHeightForWidth());
        comboBoxReference->setSizePolicy(sizePolicy1);
        comboBoxReference->setMaximumSize(QSize(200, 16777215));

        gridLayout->addWidget(comboBoxReference, 1, 1, 1, 3);

        labelTolerance = new QLabel(ToolWindowEdit);
        labelTolerance->setObjectName(QString::fromUtf8("labelTolerance"));

        gridLayout->addWidget(labelTolerance, 2, 0, 1, 1);

        labelContourValue = new QLabel(ToolWindowEdit);
        labelContourValue->setObjectName(QString::fromUtf8("labelContourValue"));

        gridLayout->addWidget(labelContourValue, 3, 0, 1, 1);

        lineEditContourValue = new QLineEdit(ToolWindowEdit);
        lineEditContourValue->setObjectName(QString::fromUtf8("lineEditContourValue"));
        lineEditContourValue->setMinimumSize(QSize(0, 0));
        lineEditContourValue->setMaximumSize(QSize(80, 16777215));
        lineEditContourValue->setBaseSize(QSize(0, 0));

        gridLayout->addWidget(lineEditContourValue, 3, 1, 1, 1);

        spinBoxTolerance = new QSpinBox(ToolWindowEdit);
        spinBoxTolerance->setObjectName(QString::fromUtf8("spinBoxTolerance"));

        gridLayout->addWidget(spinBoxTolerance, 2, 1, 1, 1);

        horizontalSpacer_2 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_2, 2, 3, 1, 1);

        horizontalSpacer_3 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_3, 3, 3, 1, 1);


        verticalLayout->addLayout(gridLayout);

        checkBoxConstrain = new QCheckBox(ToolWindowEdit);
        checkBoxConstrain->setObjectName(QString::fromUtf8("checkBoxConstrain"));

        verticalLayout->addWidget(checkBoxConstrain);

        checkBoxDrawRange = new QCheckBox(ToolWindowEdit);
        checkBoxDrawRange->setObjectName(QString::fromUtf8("checkBoxDrawRange"));

        verticalLayout->addWidget(checkBoxDrawRange);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        labelDrawRangeLow = new QLabel(ToolWindowEdit);
        labelDrawRangeLow->setObjectName(QString::fromUtf8("labelDrawRangeLow"));
        labelDrawRangeLow->setEnabled(false);

        horizontalLayout->addWidget(labelDrawRangeLow);

        lineEditDrawRangeLow = new QLineEdit(ToolWindowEdit);
        lineEditDrawRangeLow->setObjectName(QString::fromUtf8("lineEditDrawRangeLow"));
        lineEditDrawRangeLow->setEnabled(false);
        lineEditDrawRangeLow->setMaximumSize(QSize(60, 16777215));

        horizontalLayout->addWidget(lineEditDrawRangeLow);

        horizontalSpacer_8 = new QSpacerItem(10, 0, QSizePolicy::Fixed, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_8);

        labelDrawRangeHigh = new QLabel(ToolWindowEdit);
        labelDrawRangeHigh->setObjectName(QString::fromUtf8("labelDrawRangeHigh"));
        labelDrawRangeHigh->setEnabled(false);

        horizontalLayout->addWidget(labelDrawRangeHigh);

        lineEditDrawRangeHigh = new QLineEdit(ToolWindowEdit);
        lineEditDrawRangeHigh->setObjectName(QString::fromUtf8("lineEditDrawRangeHigh"));
        lineEditDrawRangeHigh->setEnabled(false);
        lineEditDrawRangeHigh->setMaximumSize(QSize(60, 16777215));

        horizontalLayout->addWidget(lineEditDrawRangeHigh);

        horizontalSpacer_5 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_5);


        verticalLayout->addLayout(horizontalLayout);

        checkBoxExcludeRange = new QCheckBox(ToolWindowEdit);
        checkBoxExcludeRange->setObjectName(QString::fromUtf8("checkBoxExcludeRange"));

        verticalLayout->addWidget(checkBoxExcludeRange);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        labelExcludeRangeLow = new QLabel(ToolWindowEdit);
        labelExcludeRangeLow->setObjectName(QString::fromUtf8("labelExcludeRangeLow"));
        labelExcludeRangeLow->setEnabled(false);

        horizontalLayout_2->addWidget(labelExcludeRangeLow);

        lineEditExcludeRangeLow = new QLineEdit(ToolWindowEdit);
        lineEditExcludeRangeLow->setObjectName(QString::fromUtf8("lineEditExcludeRangeLow"));
        lineEditExcludeRangeLow->setEnabled(false);
        lineEditExcludeRangeLow->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_2->addWidget(lineEditExcludeRangeLow);

        horizontalSpacer_9 = new QSpacerItem(10, 0, QSizePolicy::Fixed, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_9);

        labelExcludeRangeHigh = new QLabel(ToolWindowEdit);
        labelExcludeRangeHigh->setObjectName(QString::fromUtf8("labelExcludeRangeHigh"));
        labelExcludeRangeHigh->setEnabled(false);

        horizontalLayout_2->addWidget(labelExcludeRangeHigh);

        lineEditExcludeRangeHigh = new QLineEdit(ToolWindowEdit);
        lineEditExcludeRangeHigh->setObjectName(QString::fromUtf8("lineEditExcludeRangeHigh"));
        lineEditExcludeRangeHigh->setEnabled(false);
        lineEditExcludeRangeHigh->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_2->addWidget(lineEditExcludeRangeHigh);

        horizontalSpacer_6 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_6);


        verticalLayout->addLayout(horizontalLayout_2);

        checkBoxSmooth = new QCheckBox(ToolWindowEdit);
        checkBoxSmooth->setObjectName(QString::fromUtf8("checkBoxSmooth"));

        verticalLayout->addWidget(checkBoxSmooth);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        labelSD = new QLabel(ToolWindowEdit);
        labelSD->setObjectName(QString::fromUtf8("labelSD"));
        labelSD->setEnabled(false);

        horizontalLayout_3->addWidget(labelSD);

        lineEditSD = new QLineEdit(ToolWindowEdit);
        lineEditSD->setObjectName(QString::fromUtf8("lineEditSD"));
        lineEditSD->setEnabled(false);
        lineEditSD->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_3->addWidget(lineEditSD);

        horizontalSpacer_7 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_7);


        verticalLayout->addLayout(horizontalLayout_3);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        labelContourColor = new QLabel(ToolWindowEdit);
        labelContourColor->setObjectName(QString::fromUtf8("labelContourColor"));

        horizontalLayout_4->addWidget(labelContourColor);

        colorPickerContour = new QtColorPicker(ToolWindowEdit);
        colorPickerContour->setObjectName(QString::fromUtf8("colorPickerContour"));
        colorPickerContour->setMinimumSize(QSize(60, 0));

        horizontalLayout_4->addWidget(colorPickerContour);

        horizontalSpacer_4 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_4);


        verticalLayout->addLayout(horizontalLayout_4);

        checkBoxFill3D = new QCheckBox(ToolWindowEdit);
        checkBoxFill3D->setObjectName(QString::fromUtf8("checkBoxFill3D"));

        verticalLayout->addWidget(checkBoxFill3D);

        verticalSpacer_2 = new QSpacerItem(5, 3, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer_2);

        labelTips = new QLabel(ToolWindowEdit);
        labelTips->setObjectName(QString::fromUtf8("labelTips"));

        verticalLayout->addWidget(labelTips);

        verticalSpacer = new QSpacerItem(0, 5, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);


        verticalLayoutMain->addLayout(verticalLayout);


        toolBar->addAction(actionFreeHand);
        toolBar->addAction(actionPolyLine);
        toolBar->addAction(actionLiveWire);
        toolBar->addAction(actionFill);
        toolBar->addAction(actionContour);
        toolBar->addAction(actionColorPicker);
        toolBar->addSeparator();
        toolBar->addAction(actionReplaceLabel);

        retranslateUi(ToolWindowEdit);
        QObject::connect(checkBoxDrawRange, SIGNAL(toggled(bool)), labelDrawRangeLow, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxDrawRange, SIGNAL(toggled(bool)), lineEditDrawRangeHigh, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxDrawRange, SIGNAL(toggled(bool)), lineEditDrawRangeLow, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxDrawRange, SIGNAL(toggled(bool)), labelDrawRangeHigh, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxExcludeRange, SIGNAL(toggled(bool)), lineEditExcludeRangeLow, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxExcludeRange, SIGNAL(toggled(bool)), labelExcludeRangeLow, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxExcludeRange, SIGNAL(toggled(bool)), labelExcludeRangeHigh, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxSmooth, SIGNAL(toggled(bool)), labelSD, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxExcludeRange, SIGNAL(toggled(bool)), lineEditExcludeRangeHigh, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxSmooth, SIGNAL(toggled(bool)), lineEditSD, SLOT(setEnabled(bool)));
        QObject::connect(comboBoxReference, SIGNAL(currentIndexChanged(int)), ToolWindowEdit, SLOT(OnComboReference(int)));
        QObject::connect(lineEditContourValue, SIGNAL(textEdited(QString)), ToolWindowEdit, SLOT(OnLineEditContourValue(QString)));
        QObject::connect(lineEditDrawRangeLow, SIGNAL(textEdited(QString)), ToolWindowEdit, SLOT(OnDrawRangeChanged(QString)));
        QObject::connect(lineEditDrawRangeHigh, SIGNAL(textEdited(QString)), ToolWindowEdit, SLOT(OnDrawRangeChanged(QString)));
        QObject::connect(lineEditExcludeRangeLow, SIGNAL(textEdited(QString)), ToolWindowEdit, SLOT(OnExcludeRangeChanged(QString)));
        QObject::connect(lineEditExcludeRangeHigh, SIGNAL(textEdited(QString)), ToolWindowEdit, SLOT(OnExcludeRangeChanged(QString)));
        QObject::connect(lineEditSD, SIGNAL(textEdited(QString)), ToolWindowEdit, SLOT(OnLineEditSmoothSD(QString)));
        QObject::connect(actionReplaceLabel, SIGNAL(triggered()), ToolWindowEdit, SLOT(OnReplaceLabel()));

        QMetaObject::connectSlotsByName(ToolWindowEdit);
    } // setupUi

    void retranslateUi(QWidget *ToolWindowEdit)
    {
        ToolWindowEdit->setWindowTitle(QApplication::translate("ToolWindowEdit", "Voxel Edit", 0, QApplication::UnicodeUTF8));
        actionFreeHand->setText(QApplication::translate("ToolWindowEdit", "FreeHand", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionFreeHand->setToolTip(QApplication::translate("ToolWindowEdit", "Freehand", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionPolyLine->setText(QApplication::translate("ToolWindowEdit", "PolyLine", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionPolyLine->setToolTip(QApplication::translate("ToolWindowEdit", "Polyline", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionLiveWire->setText(QApplication::translate("ToolWindowEdit", "LiveWire", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionLiveWire->setToolTip(QApplication::translate("ToolWindowEdit", "Live wire", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionFill->setText(QApplication::translate("ToolWindowEdit", "Fill", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionFill->setToolTip(QApplication::translate("ToolWindowEdit", "Flood Fill", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionContour->setText(QApplication::translate("ToolWindowEdit", "Contour", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionContour->setToolTip(QApplication::translate("ToolWindowEdit", "Contour", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionColorPicker->setText(QApplication::translate("ToolWindowEdit", "ColorPicker", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionColorPicker->setToolTip(QApplication::translate("ToolWindowEdit", "Color picker", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionReplaceLabel->setText(QApplication::translate("ToolWindowEdit", "Replace Label", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionReplaceLabel->setToolTip(QApplication::translate("ToolWindowEdit", "Replace Label", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        labelBrushSize->setText(QApplication::translate("ToolWindowEdit", "Brush size", 0, QApplication::UnicodeUTF8));
        labelReference->setText(QApplication::translate("ToolWindowEdit", "Reference", 0, QApplication::UnicodeUTF8));
        labelTolerance->setText(QApplication::translate("ToolWindowEdit", "Tolerance", 0, QApplication::UnicodeUTF8));
        labelContourValue->setText(QApplication::translate("ToolWindowEdit", "Contour value", 0, QApplication::UnicodeUTF8));
        spinBoxTolerance->setSuffix(QApplication::translate("ToolWindowEdit", "%", 0, QApplication::UnicodeUTF8));
        checkBoxConstrain->setText(QApplication::translate("ToolWindowEdit", "Constrain drawings to only voxels that are \n"
"connected to old labels", 0, QApplication::UnicodeUTF8));
        checkBoxDrawRange->setText(QApplication::translate("ToolWindowEdit", "Only draw on pixels in the range of", 0, QApplication::UnicodeUTF8));
        labelDrawRangeLow->setText(QApplication::translate("ToolWindowEdit", "Low", 0, QApplication::UnicodeUTF8));
        labelDrawRangeHigh->setText(QApplication::translate("ToolWindowEdit", "high", 0, QApplication::UnicodeUTF8));
        checkBoxExcludeRange->setText(QApplication::translate("ToolWindowEdit", "Do not draw on pixels in the range of", 0, QApplication::UnicodeUTF8));
        labelExcludeRangeLow->setText(QApplication::translate("ToolWindowEdit", "Low", 0, QApplication::UnicodeUTF8));
        labelExcludeRangeHigh->setText(QApplication::translate("ToolWindowEdit", "high", 0, QApplication::UnicodeUTF8));
        checkBoxSmooth->setText(QApplication::translate("ToolWindowEdit", "Apply Gaussian smoothing", 0, QApplication::UnicodeUTF8));
        labelSD->setText(QApplication::translate("ToolWindowEdit", "Standard deviation", 0, QApplication::UnicodeUTF8));
        labelContourColor->setText(QApplication::translate("ToolWindowEdit", "Contour color", 0, QApplication::UnicodeUTF8));
        checkBoxFill3D->setText(QApplication::translate("ToolWindowEdit", "Flood fill multiple slices at a time", 0, QApplication::UnicodeUTF8));
        labelTips->setText(QApplication::translate("ToolWindowEdit", "Tips:\n"
"Ctrl + Alt + Left button to adjust contour.\n"
"Left button to add pixels to contour.\n"
"Shift + Left button to break contour.\n"
"Ctrl + Left button to fill isolated contour.\n"
"Shift + Ctrl + Left button to move cursor.\n"
"Alt + H to show/hide contour.\n"
"Alt + L to turn on/off label outline display mode.", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class ToolWindowEdit: public Ui_ToolWindowEdit {};
} // namespace Ui

QT_END_NAMESPACE

