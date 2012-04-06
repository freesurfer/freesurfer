/********************************************************************************
** Form generated from reading UI file 'PanelPointSet.ui'
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
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QScrollArea>
#include <QtGui/QSlider>
#include <QtGui/QSpacerItem>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "LayerTreeWidget.h"
#include "qtcolorpicker.h"

QT_BEGIN_NAMESPACE

class Ui_PanelPointSet
{
public:
    QWidget *scrollAreaWidgetContents;
    QVBoxLayout *verticalLayout;
    QVBoxLayout *verticalLayoutToolbar;
    QToolBar *toolbar;
    LayerTreeWidget *treeWidgetLayers;
    QGridLayout *gridLayout;
    QLabel *labelFileName;
    QLineEdit *lineEditFileName;
    QLabel *labelOpacity;
    QHBoxLayout *horizontalLayout_2;
    QSlider *sliderOpacity;
    QDoubleSpinBox *doubleSpinBoxOpacity;
    QLabel *labelSplineColor;
    QComboBox *comboBoxSplineColor;
    QLabel *labelScalarMap;
    QComboBox *comboBoxScalarMap;
    QLabel *labelMin;
    QHBoxLayout *horizontalLayout_9;
    QSlider *sliderMin;
    QLineEdit *lineEditMin;
    QLabel *labelMid;
    QHBoxLayout *horizontalLayout_7;
    QSlider *sliderMid;
    QLineEdit *lineEditMid;
    QLabel *labelMax;
    QHBoxLayout *horizontalLayout_11;
    QSlider *sliderMax;
    QLineEdit *lineEditMax;
    QLabel *labelOffset;
    QHBoxLayout *horizontalLayout_14;
    QSlider *sliderOffset;
    QLineEdit *lineEditOffset;
    QLabel *labelSplineRadius;
    QHBoxLayout *horizontalLayout_17;
    QLineEdit *lineEditSplineRadius;
    QSpacerItem *horizontalSpacer;
    QLabel *labelRadius;
    QHBoxLayout *horizontalLayout;
    QLineEdit *lineEditRadius;
    QSpacerItem *horizontalSpacer_2;
    QCheckBox *checkBoxSnapToCenter;
    QCheckBox *checkBoxShowSpline;
    QHBoxLayout *horizontalLayout_3;
    QtColorPicker *colorpickerPointColor;
    QSpacerItem *horizontalSpacer_3;
    QLabel *labelColor_2;
    QHBoxLayout *horizontalLayout_4;
    QtColorPicker *colorpickerSplineColor;
    QSpacerItem *horizontalSpacer_4;
    QSpacerItem *verticalSpacer;

    void setupUi(QScrollArea *PanelPointSet)
    {
        if (PanelPointSet->objectName().isEmpty())
            PanelPointSet->setObjectName(QString::fromUtf8("PanelPointSet"));
        PanelPointSet->resize(314, 816);
        PanelPointSet->setFrameShape(QFrame::NoFrame);
        PanelPointSet->setWidgetResizable(true);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QString::fromUtf8("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 314, 816));
        verticalLayout = new QVBoxLayout(scrollAreaWidgetContents);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(5, 0, 5, 5);
        verticalLayoutToolbar = new QVBoxLayout();
        verticalLayoutToolbar->setSpacing(0);
        verticalLayoutToolbar->setObjectName(QString::fromUtf8("verticalLayoutToolbar"));
        verticalLayoutToolbar->setContentsMargins(-1, 0, 0, -1);
        toolbar = new QToolBar(scrollAreaWidgetContents);
        toolbar->setObjectName(QString::fromUtf8("toolbar"));
        toolbar->setIconSize(QSize(24, 24));
        toolbar->setFloatable(false);

        verticalLayoutToolbar->addWidget(toolbar);

        treeWidgetLayers = new LayerTreeWidget(scrollAreaWidgetContents);
        QTreeWidgetItem *__qtreewidgetitem = new QTreeWidgetItem();
        __qtreewidgetitem->setText(0, QString::fromUtf8("1"));
        treeWidgetLayers->setHeaderItem(__qtreewidgetitem);
        treeWidgetLayers->setObjectName(QString::fromUtf8("treeWidgetLayers"));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(treeWidgetLayers->sizePolicy().hasHeightForWidth());
        treeWidgetLayers->setSizePolicy(sizePolicy);
        treeWidgetLayers->setMaximumSize(QSize(16777215, 145));
        treeWidgetLayers->setRootIsDecorated(false);
        treeWidgetLayers->setUniformRowHeights(true);
        treeWidgetLayers->header()->setVisible(false);

        verticalLayoutToolbar->addWidget(treeWidgetLayers);


        verticalLayout->addLayout(verticalLayoutToolbar);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        labelFileName = new QLabel(scrollAreaWidgetContents);
        labelFileName->setObjectName(QString::fromUtf8("labelFileName"));

        gridLayout->addWidget(labelFileName, 0, 0, 1, 1);

        lineEditFileName = new QLineEdit(scrollAreaWidgetContents);
        lineEditFileName->setObjectName(QString::fromUtf8("lineEditFileName"));

        gridLayout->addWidget(lineEditFileName, 0, 1, 1, 1);

        labelOpacity = new QLabel(scrollAreaWidgetContents);
        labelOpacity->setObjectName(QString::fromUtf8("labelOpacity"));

        gridLayout->addWidget(labelOpacity, 1, 0, 1, 1);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        sliderOpacity = new QSlider(scrollAreaWidgetContents);
        sliderOpacity->setObjectName(QString::fromUtf8("sliderOpacity"));
        sliderOpacity->setMaximum(100);
        sliderOpacity->setValue(100);
        sliderOpacity->setOrientation(Qt::Horizontal);

        horizontalLayout_2->addWidget(sliderOpacity);

        doubleSpinBoxOpacity = new QDoubleSpinBox(scrollAreaWidgetContents);
        doubleSpinBoxOpacity->setObjectName(QString::fromUtf8("doubleSpinBoxOpacity"));
        doubleSpinBoxOpacity->setMaximum(1);
        doubleSpinBoxOpacity->setSingleStep(0.1);
        doubleSpinBoxOpacity->setValue(1);

        horizontalLayout_2->addWidget(doubleSpinBoxOpacity);


        gridLayout->addLayout(horizontalLayout_2, 1, 1, 1, 1);

        labelSplineColor = new QLabel(scrollAreaWidgetContents);
        labelSplineColor->setObjectName(QString::fromUtf8("labelSplineColor"));

        gridLayout->addWidget(labelSplineColor, 6, 0, 1, 1);

        comboBoxSplineColor = new QComboBox(scrollAreaWidgetContents);
        comboBoxSplineColor->setObjectName(QString::fromUtf8("comboBoxSplineColor"));

        gridLayout->addWidget(comboBoxSplineColor, 6, 1, 1, 1);

        labelScalarMap = new QLabel(scrollAreaWidgetContents);
        labelScalarMap->setObjectName(QString::fromUtf8("labelScalarMap"));

        gridLayout->addWidget(labelScalarMap, 8, 0, 1, 1);

        comboBoxScalarMap = new QComboBox(scrollAreaWidgetContents);
        comboBoxScalarMap->setObjectName(QString::fromUtf8("comboBoxScalarMap"));

        gridLayout->addWidget(comboBoxScalarMap, 8, 1, 1, 1);

        labelMin = new QLabel(scrollAreaWidgetContents);
        labelMin->setObjectName(QString::fromUtf8("labelMin"));

        gridLayout->addWidget(labelMin, 9, 0, 1, 1);

        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setObjectName(QString::fromUtf8("horizontalLayout_9"));
        sliderMin = new QSlider(scrollAreaWidgetContents);
        sliderMin->setObjectName(QString::fromUtf8("sliderMin"));
        sliderMin->setMaximum(100);
        sliderMin->setOrientation(Qt::Horizontal);

        horizontalLayout_9->addWidget(sliderMin);

        lineEditMin = new QLineEdit(scrollAreaWidgetContents);
        lineEditMin->setObjectName(QString::fromUtf8("lineEditMin"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(lineEditMin->sizePolicy().hasHeightForWidth());
        lineEditMin->setSizePolicy(sizePolicy1);
        lineEditMin->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_9->addWidget(lineEditMin);


        gridLayout->addLayout(horizontalLayout_9, 9, 1, 1, 1);

        labelMid = new QLabel(scrollAreaWidgetContents);
        labelMid->setObjectName(QString::fromUtf8("labelMid"));

        gridLayout->addWidget(labelMid, 10, 0, 1, 1);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setObjectName(QString::fromUtf8("horizontalLayout_7"));
        sliderMid = new QSlider(scrollAreaWidgetContents);
        sliderMid->setObjectName(QString::fromUtf8("sliderMid"));
        sliderMid->setMaximum(100);
        sliderMid->setOrientation(Qt::Horizontal);

        horizontalLayout_7->addWidget(sliderMid);

        lineEditMid = new QLineEdit(scrollAreaWidgetContents);
        lineEditMid->setObjectName(QString::fromUtf8("lineEditMid"));
        sizePolicy1.setHeightForWidth(lineEditMid->sizePolicy().hasHeightForWidth());
        lineEditMid->setSizePolicy(sizePolicy1);
        lineEditMid->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_7->addWidget(lineEditMid);


        gridLayout->addLayout(horizontalLayout_7, 10, 1, 1, 1);

        labelMax = new QLabel(scrollAreaWidgetContents);
        labelMax->setObjectName(QString::fromUtf8("labelMax"));

        gridLayout->addWidget(labelMax, 11, 0, 1, 1);

        horizontalLayout_11 = new QHBoxLayout();
        horizontalLayout_11->setObjectName(QString::fromUtf8("horizontalLayout_11"));
        sliderMax = new QSlider(scrollAreaWidgetContents);
        sliderMax->setObjectName(QString::fromUtf8("sliderMax"));
        sliderMax->setMaximum(100);
        sliderMax->setOrientation(Qt::Horizontal);

        horizontalLayout_11->addWidget(sliderMax);

        lineEditMax = new QLineEdit(scrollAreaWidgetContents);
        lineEditMax->setObjectName(QString::fromUtf8("lineEditMax"));
        sizePolicy1.setHeightForWidth(lineEditMax->sizePolicy().hasHeightForWidth());
        lineEditMax->setSizePolicy(sizePolicy1);
        lineEditMax->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_11->addWidget(lineEditMax);


        gridLayout->addLayout(horizontalLayout_11, 11, 1, 1, 1);

        labelOffset = new QLabel(scrollAreaWidgetContents);
        labelOffset->setObjectName(QString::fromUtf8("labelOffset"));

        gridLayout->addWidget(labelOffset, 12, 0, 1, 1);

        horizontalLayout_14 = new QHBoxLayout();
        horizontalLayout_14->setObjectName(QString::fromUtf8("horizontalLayout_14"));
        sliderOffset = new QSlider(scrollAreaWidgetContents);
        sliderOffset->setObjectName(QString::fromUtf8("sliderOffset"));
        sliderOffset->setMaximum(100);
        sliderOffset->setOrientation(Qt::Horizontal);

        horizontalLayout_14->addWidget(sliderOffset);

        lineEditOffset = new QLineEdit(scrollAreaWidgetContents);
        lineEditOffset->setObjectName(QString::fromUtf8("lineEditOffset"));
        sizePolicy1.setHeightForWidth(lineEditOffset->sizePolicy().hasHeightForWidth());
        lineEditOffset->setSizePolicy(sizePolicy1);
        lineEditOffset->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_14->addWidget(lineEditOffset);


        gridLayout->addLayout(horizontalLayout_14, 12, 1, 1, 1);

        labelSplineRadius = new QLabel(scrollAreaWidgetContents);
        labelSplineRadius->setObjectName(QString::fromUtf8("labelSplineRadius"));

        gridLayout->addWidget(labelSplineRadius, 13, 0, 1, 1);

        horizontalLayout_17 = new QHBoxLayout();
        horizontalLayout_17->setObjectName(QString::fromUtf8("horizontalLayout_17"));
        lineEditSplineRadius = new QLineEdit(scrollAreaWidgetContents);
        lineEditSplineRadius->setObjectName(QString::fromUtf8("lineEditSplineRadius"));
        lineEditSplineRadius->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_17->addWidget(lineEditSplineRadius);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_17->addItem(horizontalSpacer);


        gridLayout->addLayout(horizontalLayout_17, 13, 1, 1, 1);

        labelRadius = new QLabel(scrollAreaWidgetContents);
        labelRadius->setObjectName(QString::fromUtf8("labelRadius"));

        gridLayout->addWidget(labelRadius, 3, 0, 1, 1);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        lineEditRadius = new QLineEdit(scrollAreaWidgetContents);
        lineEditRadius->setObjectName(QString::fromUtf8("lineEditRadius"));
        lineEditRadius->setMaximumSize(QSize(60, 16777215));

        horizontalLayout->addWidget(lineEditRadius);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_2);


        gridLayout->addLayout(horizontalLayout, 3, 1, 1, 1);

        checkBoxSnapToCenter = new QCheckBox(scrollAreaWidgetContents);
        checkBoxSnapToCenter->setObjectName(QString::fromUtf8("checkBoxSnapToCenter"));

        gridLayout->addWidget(checkBoxSnapToCenter, 4, 1, 1, 1);

        checkBoxShowSpline = new QCheckBox(scrollAreaWidgetContents);
        checkBoxShowSpline->setObjectName(QString::fromUtf8("checkBoxShowSpline"));
        checkBoxShowSpline->setChecked(true);

        gridLayout->addWidget(checkBoxShowSpline, 5, 1, 1, 1);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        colorpickerPointColor = new QtColorPicker(scrollAreaWidgetContents);
        colorpickerPointColor->setObjectName(QString::fromUtf8("colorpickerPointColor"));

        horizontalLayout_3->addWidget(colorpickerPointColor);

        horizontalSpacer_3 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_3);


        gridLayout->addLayout(horizontalLayout_3, 2, 1, 1, 1);

        labelColor_2 = new QLabel(scrollAreaWidgetContents);
        labelColor_2->setObjectName(QString::fromUtf8("labelColor_2"));

        gridLayout->addWidget(labelColor_2, 2, 0, 1, 1);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        colorpickerSplineColor = new QtColorPicker(scrollAreaWidgetContents);
        colorpickerSplineColor->setObjectName(QString::fromUtf8("colorpickerSplineColor"));

        horizontalLayout_4->addWidget(colorpickerSplineColor);

        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_4);


        gridLayout->addLayout(horizontalLayout_4, 7, 1, 1, 1);


        verticalLayout->addLayout(gridLayout);

        verticalSpacer = new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        PanelPointSet->setWidget(scrollAreaWidgetContents);

        retranslateUi(PanelPointSet);
        QObject::connect(sliderOpacity, SIGNAL(valueChanged(int)), PanelPointSet, SLOT(OnSliderOpacity(int)));
        QObject::connect(lineEditRadius, SIGNAL(textEdited(QString)), PanelPointSet, SLOT(OnLineEditRadius(QString)));
        QObject::connect(sliderMin, SIGNAL(valueChanged(int)), PanelPointSet, SLOT(OnSliderMin(int)));
        QObject::connect(sliderMid, SIGNAL(valueChanged(int)), PanelPointSet, SLOT(OnSliderMid(int)));
        QObject::connect(sliderMax, SIGNAL(valueChanged(int)), PanelPointSet, SLOT(OnSliderMax(int)));
        QObject::connect(sliderOffset, SIGNAL(valueChanged(int)), PanelPointSet, SLOT(OnSliderOffset(int)));
        QObject::connect(lineEditSplineRadius, SIGNAL(textEdited(QString)), PanelPointSet, SLOT(OnLineEditSplineRadius(QString)));
        QObject::connect(lineEditMin, SIGNAL(textEdited(QString)), PanelPointSet, SLOT(OnLineEditMin(QString)));
        QObject::connect(lineEditMid, SIGNAL(textEdited(QString)), PanelPointSet, SLOT(OnLineEditMid(QString)));
        QObject::connect(lineEditMax, SIGNAL(textEdited(QString)), PanelPointSet, SLOT(OnLineEditMax(QString)));
        QObject::connect(lineEditOffset, SIGNAL(textEdited(QString)), PanelPointSet, SLOT(OnLineEditOffset(QString)));
        QObject::connect(comboBoxScalarMap, SIGNAL(currentIndexChanged(int)), PanelPointSet, SLOT(OnComboScalarMap(int)));

        QMetaObject::connectSlotsByName(PanelPointSet);
    } // setupUi

    void retranslateUi(QScrollArea *PanelPointSet)
    {
        PanelPointSet->setWindowTitle(QApplication::translate("PanelPointSet", "ScrollArea", 0, QApplication::UnicodeUTF8));
        labelFileName->setText(QApplication::translate("PanelPointSet", "File name", 0, QApplication::UnicodeUTF8));
        labelOpacity->setText(QApplication::translate("PanelPointSet", "Opacity", 0, QApplication::UnicodeUTF8));
        labelSplineColor->setText(QApplication::translate("PanelPointSet", "Spline color", 0, QApplication::UnicodeUTF8));
        comboBoxSplineColor->clear();
        comboBoxSplineColor->insertItems(0, QStringList()
         << QApplication::translate("PanelPointSet", "Solid Color", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelPointSet", "Heatscale", 0, QApplication::UnicodeUTF8)
        );
        labelScalarMap->setText(QApplication::translate("PanelPointSet", "Scalar map", 0, QApplication::UnicodeUTF8));
        labelMin->setText(QApplication::translate("PanelPointSet", "Min", 0, QApplication::UnicodeUTF8));
        labelMid->setText(QApplication::translate("PanelPointSet", "Mid", 0, QApplication::UnicodeUTF8));
        labelMax->setText(QApplication::translate("PanelPointSet", "Max", 0, QApplication::UnicodeUTF8));
        labelOffset->setText(QApplication::translate("PanelPointSet", "Offset", 0, QApplication::UnicodeUTF8));
        labelSplineRadius->setText(QApplication::translate("PanelPointSet", "Spline radius", 0, QApplication::UnicodeUTF8));
        labelRadius->setText(QApplication::translate("PanelPointSet", "Radius", 0, QApplication::UnicodeUTF8));
        checkBoxSnapToCenter->setText(QApplication::translate("PanelPointSet", "Snap to voxel center", 0, QApplication::UnicodeUTF8));
        checkBoxShowSpline->setText(QApplication::translate("PanelPointSet", "Show spline", 0, QApplication::UnicodeUTF8));
        labelColor_2->setText(QApplication::translate("PanelPointSet", "Color", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class PanelPointSet: public Ui_PanelPointSet {};
} // namespace Ui

QT_END_NAMESPACE

