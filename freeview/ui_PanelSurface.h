/********************************************************************************
** Form generated from reading UI file 'PanelSurface.ui'
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
#include <QtGui/QPushButton>
#include <QtGui/QScrollArea>
#include <QtGui/QSlider>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "LayerTreeWidget.h"
#include "qtcolorpicker.h"

QT_BEGIN_NAMESPACE

class Ui_PanelSurface
{
public:
    QAction *actionSurfaceMain;
    QAction *actionSurfaceInflated;
    QAction *actionSurfaceWhite;
    QAction *actionSurfacePial;
    QAction *actionSurfaceOriginal;
    QAction *actionLockLayer;
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
    QLabel *labelSurfaceColor;
    QHBoxLayout *horizontalLayout21;
    QtColorPicker *colorpickerSurfaceColor;
    QSpacerItem *horizontalSpacer_21;
    QLabel *labelRender;
    QComboBox *comboBoxRender;
    QLabel *labelMeshColor;
    QComboBox *comboBoxMeshColor;
    QCheckBox *checkBoxShowVertices;
    QLabel *labelVertexColor;
    QHBoxLayout *horizontalLayout20;
    QtColorPicker *colorpickerVertexColor;
    QSpacerItem *horizontalSpacer_20;
    QLabel *labelVertexPointSize;
    QHBoxLayout *horizontalLayout;
    QSpinBox *spinBoxVertexPointSize;
    QSpacerItem *horizontalSpacer_2;
    QLabel *labelCurvature;
    QComboBox *comboBoxCurvature;
    QLabel *labelMidPoint;
    QHBoxLayout *horizontalLayout_6;
    QSlider *sliderMidPoint;
    QLineEdit *lineEditMidPoint;
    QLabel *labelSlope;
    QHBoxLayout *horizontalLayout_10;
    QSlider *sliderSlope;
    QLineEdit *lineEditSlope;
    QLabel *labelOverlay;
    QComboBox *comboBoxOverlay;
    QPushButton *pushButtonConfigureOverlay;
    QLabel *labelAnnotation;
    QComboBox *comboBoxAnnotation;
    QLabel *labelLabel;
    QComboBox *comboBoxLabel;
    QLabel *labelLabelColor;
    QHBoxLayout *horizontalLayout22;
    QtColorPicker *colorpickerLabelColor;
    QCheckBox *checkBoxLabelOutline;
    QSpacerItem *horizontalSpacer_22;
    QLabel *labelEdgeColor;
    QHBoxLayout *horizontalLayout23;
    QtColorPicker *colorpickerEdgeColor;
    QSpacerItem *horizontalSpacer_23;
    QLabel *labelEdgeThickness;
    QHBoxLayout *horizontalLayout_17;
    QSpinBox *spinBoxEdgeThickness;
    QSpacerItem *horizontalSpacer;
    QLabel *labelVectorDisplay;
    QComboBox *comboBoxVectorDisplay;
    QLabel *labelVectorColor;
    QHBoxLayout *horizontalLayout24;
    QtColorPicker *colorpickerVectorColor;
    QSpacerItem *horizontalSpacer_24;
    QLabel *labelVectorPointSize;
    QHBoxLayout *horizontalLayout_3;
    QSpinBox *spinBoxVectorPointSize;
    QSpacerItem *horizontalSpacer_3;
    QLabel *labelPositionOffset;
    QLineEdit *lineEditPositionOffset;
    QHBoxLayout *horizontalLayout22_3;
    QCheckBox *checkBoxAnnotationOutline;
    QSpacerItem *horizontalSpacer_30;
    QCheckBox *checkBoxShowInfo;
    QSpacerItem *verticalSpacer;

    void setupUi(QScrollArea *PanelSurface)
    {
        if (PanelSurface->objectName().isEmpty())
            PanelSurface->setObjectName(QString::fromUtf8("PanelSurface"));
        PanelSurface->resize(332, 1030);
        PanelSurface->setFrameShape(QFrame::NoFrame);
        PanelSurface->setWidgetResizable(true);
        actionSurfaceMain = new QAction(PanelSurface);
        actionSurfaceMain->setObjectName(QString::fromUtf8("actionSurfaceMain"));
        actionSurfaceMain->setCheckable(true);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/surface_main.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSurfaceMain->setIcon(icon);
        actionSurfaceInflated = new QAction(PanelSurface);
        actionSurfaceInflated->setObjectName(QString::fromUtf8("actionSurfaceInflated"));
        actionSurfaceInflated->setCheckable(true);
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/resource/icons/surface_inflated.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSurfaceInflated->setIcon(icon1);
        actionSurfaceWhite = new QAction(PanelSurface);
        actionSurfaceWhite->setObjectName(QString::fromUtf8("actionSurfaceWhite"));
        actionSurfaceWhite->setCheckable(true);
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/resource/icons/surface_white.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSurfaceWhite->setIcon(icon2);
        actionSurfacePial = new QAction(PanelSurface);
        actionSurfacePial->setObjectName(QString::fromUtf8("actionSurfacePial"));
        actionSurfacePial->setCheckable(true);
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/resource/icons/surface_pial.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSurfacePial->setIcon(icon3);
        actionSurfaceOriginal = new QAction(PanelSurface);
        actionSurfaceOriginal->setObjectName(QString::fromUtf8("actionSurfaceOriginal"));
        actionSurfaceOriginal->setCheckable(true);
        QIcon icon4;
        icon4.addFile(QString::fromUtf8(":/resource/icons/surface_original.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionSurfaceOriginal->setIcon(icon4);
        actionLockLayer = new QAction(PanelSurface);
        actionLockLayer->setObjectName(QString::fromUtf8("actionLockLayer"));
        actionLockLayer->setCheckable(true);
        QIcon icon5;
        icon5.addFile(QString::fromUtf8(":/resource/icons/volume_lock.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLockLayer->setIcon(icon5);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QString::fromUtf8("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 332, 1030));
        verticalLayout = new QVBoxLayout(scrollAreaWidgetContents);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(5, 0, 5, 5);
        verticalLayoutToolbar = new QVBoxLayout();
        verticalLayoutToolbar->setSpacing(0);
        verticalLayoutToolbar->setObjectName(QString::fromUtf8("verticalLayoutToolbar"));
        verticalLayoutToolbar->setContentsMargins(-1, -1, -1, 0);
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

        labelSurfaceColor = new QLabel(scrollAreaWidgetContents);
        labelSurfaceColor->setObjectName(QString::fromUtf8("labelSurfaceColor"));

        gridLayout->addWidget(labelSurfaceColor, 2, 0, 1, 1);

        horizontalLayout21 = new QHBoxLayout();
        horizontalLayout21->setObjectName(QString::fromUtf8("horizontalLayout21"));
        colorpickerSurfaceColor = new QtColorPicker(scrollAreaWidgetContents);
        colorpickerSurfaceColor->setObjectName(QString::fromUtf8("colorpickerSurfaceColor"));

        horizontalLayout21->addWidget(colorpickerSurfaceColor);

        horizontalSpacer_21 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout21->addItem(horizontalSpacer_21);


        gridLayout->addLayout(horizontalLayout21, 2, 1, 1, 1);

        labelRender = new QLabel(scrollAreaWidgetContents);
        labelRender->setObjectName(QString::fromUtf8("labelRender"));

        gridLayout->addWidget(labelRender, 3, 0, 1, 1);

        comboBoxRender = new QComboBox(scrollAreaWidgetContents);
        comboBoxRender->setObjectName(QString::fromUtf8("comboBoxRender"));

        gridLayout->addWidget(comboBoxRender, 3, 1, 1, 1);

        labelMeshColor = new QLabel(scrollAreaWidgetContents);
        labelMeshColor->setObjectName(QString::fromUtf8("labelMeshColor"));
        labelMeshColor->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelMeshColor, 4, 0, 1, 1);

        comboBoxMeshColor = new QComboBox(scrollAreaWidgetContents);
        comboBoxMeshColor->setObjectName(QString::fromUtf8("comboBoxMeshColor"));

        gridLayout->addWidget(comboBoxMeshColor, 4, 1, 1, 1);

        checkBoxShowVertices = new QCheckBox(scrollAreaWidgetContents);
        checkBoxShowVertices->setObjectName(QString::fromUtf8("checkBoxShowVertices"));

        gridLayout->addWidget(checkBoxShowVertices, 5, 1, 1, 1);

        labelVertexColor = new QLabel(scrollAreaWidgetContents);
        labelVertexColor->setObjectName(QString::fromUtf8("labelVertexColor"));
        labelVertexColor->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelVertexColor, 6, 0, 1, 1);

        horizontalLayout20 = new QHBoxLayout();
        horizontalLayout20->setObjectName(QString::fromUtf8("horizontalLayout20"));
        colorpickerVertexColor = new QtColorPicker(scrollAreaWidgetContents);
        colorpickerVertexColor->setObjectName(QString::fromUtf8("colorpickerVertexColor"));

        horizontalLayout20->addWidget(colorpickerVertexColor);

        horizontalSpacer_20 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout20->addItem(horizontalSpacer_20);


        gridLayout->addLayout(horizontalLayout20, 6, 1, 1, 1);

        labelVertexPointSize = new QLabel(scrollAreaWidgetContents);
        labelVertexPointSize->setObjectName(QString::fromUtf8("labelVertexPointSize"));
        labelVertexPointSize->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelVertexPointSize, 7, 0, 1, 1);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        spinBoxVertexPointSize = new QSpinBox(scrollAreaWidgetContents);
        spinBoxVertexPointSize->setObjectName(QString::fromUtf8("spinBoxVertexPointSize"));
        spinBoxVertexPointSize->setMinimum(1);
        spinBoxVertexPointSize->setMaximum(100);

        horizontalLayout->addWidget(spinBoxVertexPointSize);

        horizontalSpacer_2 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_2);


        gridLayout->addLayout(horizontalLayout, 7, 1, 1, 1);

        labelCurvature = new QLabel(scrollAreaWidgetContents);
        labelCurvature->setObjectName(QString::fromUtf8("labelCurvature"));

        gridLayout->addWidget(labelCurvature, 8, 0, 1, 1);

        comboBoxCurvature = new QComboBox(scrollAreaWidgetContents);
        comboBoxCurvature->setObjectName(QString::fromUtf8("comboBoxCurvature"));

        gridLayout->addWidget(comboBoxCurvature, 8, 1, 1, 1);

        labelMidPoint = new QLabel(scrollAreaWidgetContents);
        labelMidPoint->setObjectName(QString::fromUtf8("labelMidPoint"));
        labelMidPoint->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelMidPoint, 9, 0, 1, 1);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        sliderMidPoint = new QSlider(scrollAreaWidgetContents);
        sliderMidPoint->setObjectName(QString::fromUtf8("sliderMidPoint"));
        sliderMidPoint->setMaximum(100);
        sliderMidPoint->setOrientation(Qt::Horizontal);

        horizontalLayout_6->addWidget(sliderMidPoint);

        lineEditMidPoint = new QLineEdit(scrollAreaWidgetContents);
        lineEditMidPoint->setObjectName(QString::fromUtf8("lineEditMidPoint"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(lineEditMidPoint->sizePolicy().hasHeightForWidth());
        lineEditMidPoint->setSizePolicy(sizePolicy1);
        lineEditMidPoint->setMaximumSize(QSize(75, 16777215));

        horizontalLayout_6->addWidget(lineEditMidPoint);


        gridLayout->addLayout(horizontalLayout_6, 9, 1, 1, 1);

        labelSlope = new QLabel(scrollAreaWidgetContents);
        labelSlope->setObjectName(QString::fromUtf8("labelSlope"));
        labelSlope->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelSlope, 10, 0, 1, 1);

        horizontalLayout_10 = new QHBoxLayout();
        horizontalLayout_10->setObjectName(QString::fromUtf8("horizontalLayout_10"));
        sliderSlope = new QSlider(scrollAreaWidgetContents);
        sliderSlope->setObjectName(QString::fromUtf8("sliderSlope"));
        sliderSlope->setMaximum(100);
        sliderSlope->setOrientation(Qt::Horizontal);

        horizontalLayout_10->addWidget(sliderSlope);

        lineEditSlope = new QLineEdit(scrollAreaWidgetContents);
        lineEditSlope->setObjectName(QString::fromUtf8("lineEditSlope"));
        sizePolicy1.setHeightForWidth(lineEditSlope->sizePolicy().hasHeightForWidth());
        lineEditSlope->setSizePolicy(sizePolicy1);
        lineEditSlope->setMaximumSize(QSize(75, 16777215));

        horizontalLayout_10->addWidget(lineEditSlope);


        gridLayout->addLayout(horizontalLayout_10, 10, 1, 1, 1);

        labelOverlay = new QLabel(scrollAreaWidgetContents);
        labelOverlay->setObjectName(QString::fromUtf8("labelOverlay"));

        gridLayout->addWidget(labelOverlay, 11, 0, 1, 1);

        comboBoxOverlay = new QComboBox(scrollAreaWidgetContents);
        comboBoxOverlay->setObjectName(QString::fromUtf8("comboBoxOverlay"));

        gridLayout->addWidget(comboBoxOverlay, 11, 1, 1, 1);

        pushButtonConfigureOverlay = new QPushButton(scrollAreaWidgetContents);
        pushButtonConfigureOverlay->setObjectName(QString::fromUtf8("pushButtonConfigureOverlay"));

        gridLayout->addWidget(pushButtonConfigureOverlay, 12, 1, 1, 1);

        labelAnnotation = new QLabel(scrollAreaWidgetContents);
        labelAnnotation->setObjectName(QString::fromUtf8("labelAnnotation"));

        gridLayout->addWidget(labelAnnotation, 13, 0, 1, 1);

        comboBoxAnnotation = new QComboBox(scrollAreaWidgetContents);
        comboBoxAnnotation->setObjectName(QString::fromUtf8("comboBoxAnnotation"));

        gridLayout->addWidget(comboBoxAnnotation, 13, 1, 1, 1);

        labelLabel = new QLabel(scrollAreaWidgetContents);
        labelLabel->setObjectName(QString::fromUtf8("labelLabel"));

        gridLayout->addWidget(labelLabel, 15, 0, 1, 1);

        comboBoxLabel = new QComboBox(scrollAreaWidgetContents);
        comboBoxLabel->setObjectName(QString::fromUtf8("comboBoxLabel"));

        gridLayout->addWidget(comboBoxLabel, 15, 1, 1, 1);

        labelLabelColor = new QLabel(scrollAreaWidgetContents);
        labelLabelColor->setObjectName(QString::fromUtf8("labelLabelColor"));
        labelLabelColor->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelLabelColor, 16, 0, 1, 1);

        horizontalLayout22 = new QHBoxLayout();
        horizontalLayout22->setObjectName(QString::fromUtf8("horizontalLayout22"));
        colorpickerLabelColor = new QtColorPicker(scrollAreaWidgetContents);
        colorpickerLabelColor->setObjectName(QString::fromUtf8("colorpickerLabelColor"));

        horizontalLayout22->addWidget(colorpickerLabelColor);

        checkBoxLabelOutline = new QCheckBox(scrollAreaWidgetContents);
        checkBoxLabelOutline->setObjectName(QString::fromUtf8("checkBoxLabelOutline"));

        horizontalLayout22->addWidget(checkBoxLabelOutline);

        horizontalSpacer_22 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout22->addItem(horizontalSpacer_22);


        gridLayout->addLayout(horizontalLayout22, 16, 1, 1, 1);

        labelEdgeColor = new QLabel(scrollAreaWidgetContents);
        labelEdgeColor->setObjectName(QString::fromUtf8("labelEdgeColor"));

        gridLayout->addWidget(labelEdgeColor, 17, 0, 1, 1);

        horizontalLayout23 = new QHBoxLayout();
        horizontalLayout23->setObjectName(QString::fromUtf8("horizontalLayout23"));
        colorpickerEdgeColor = new QtColorPicker(scrollAreaWidgetContents);
        colorpickerEdgeColor->setObjectName(QString::fromUtf8("colorpickerEdgeColor"));

        horizontalLayout23->addWidget(colorpickerEdgeColor);

        horizontalSpacer_23 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout23->addItem(horizontalSpacer_23);


        gridLayout->addLayout(horizontalLayout23, 17, 1, 1, 1);

        labelEdgeThickness = new QLabel(scrollAreaWidgetContents);
        labelEdgeThickness->setObjectName(QString::fromUtf8("labelEdgeThickness"));

        gridLayout->addWidget(labelEdgeThickness, 18, 0, 1, 1);

        horizontalLayout_17 = new QHBoxLayout();
        horizontalLayout_17->setObjectName(QString::fromUtf8("horizontalLayout_17"));
        spinBoxEdgeThickness = new QSpinBox(scrollAreaWidgetContents);
        spinBoxEdgeThickness->setObjectName(QString::fromUtf8("spinBoxEdgeThickness"));

        horizontalLayout_17->addWidget(spinBoxEdgeThickness);

        horizontalSpacer = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_17->addItem(horizontalSpacer);


        gridLayout->addLayout(horizontalLayout_17, 18, 1, 1, 1);

        labelVectorDisplay = new QLabel(scrollAreaWidgetContents);
        labelVectorDisplay->setObjectName(QString::fromUtf8("labelVectorDisplay"));

        gridLayout->addWidget(labelVectorDisplay, 19, 0, 1, 1);

        comboBoxVectorDisplay = new QComboBox(scrollAreaWidgetContents);
        comboBoxVectorDisplay->setObjectName(QString::fromUtf8("comboBoxVectorDisplay"));

        gridLayout->addWidget(comboBoxVectorDisplay, 19, 1, 1, 1);

        labelVectorColor = new QLabel(scrollAreaWidgetContents);
        labelVectorColor->setObjectName(QString::fromUtf8("labelVectorColor"));
        labelVectorColor->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelVectorColor, 20, 0, 1, 1);

        horizontalLayout24 = new QHBoxLayout();
        horizontalLayout24->setObjectName(QString::fromUtf8("horizontalLayout24"));
        colorpickerVectorColor = new QtColorPicker(scrollAreaWidgetContents);
        colorpickerVectorColor->setObjectName(QString::fromUtf8("colorpickerVectorColor"));

        horizontalLayout24->addWidget(colorpickerVectorColor);

        horizontalSpacer_24 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout24->addItem(horizontalSpacer_24);


        gridLayout->addLayout(horizontalLayout24, 20, 1, 1, 1);

        labelVectorPointSize = new QLabel(scrollAreaWidgetContents);
        labelVectorPointSize->setObjectName(QString::fromUtf8("labelVectorPointSize"));
        labelVectorPointSize->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelVectorPointSize, 21, 0, 1, 1);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        spinBoxVectorPointSize = new QSpinBox(scrollAreaWidgetContents);
        spinBoxVectorPointSize->setObjectName(QString::fromUtf8("spinBoxVectorPointSize"));
        spinBoxVectorPointSize->setMinimum(1);
        spinBoxVectorPointSize->setMaximum(100);

        horizontalLayout_3->addWidget(spinBoxVectorPointSize);

        horizontalSpacer_3 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_3);


        gridLayout->addLayout(horizontalLayout_3, 21, 1, 1, 1);

        labelPositionOffset = new QLabel(scrollAreaWidgetContents);
        labelPositionOffset->setObjectName(QString::fromUtf8("labelPositionOffset"));

        gridLayout->addWidget(labelPositionOffset, 22, 0, 1, 1);

        lineEditPositionOffset = new QLineEdit(scrollAreaWidgetContents);
        lineEditPositionOffset->setObjectName(QString::fromUtf8("lineEditPositionOffset"));

        gridLayout->addWidget(lineEditPositionOffset, 22, 1, 1, 1);

        horizontalLayout22_3 = new QHBoxLayout();
        horizontalLayout22_3->setObjectName(QString::fromUtf8("horizontalLayout22_3"));
        checkBoxAnnotationOutline = new QCheckBox(scrollAreaWidgetContents);
        checkBoxAnnotationOutline->setObjectName(QString::fromUtf8("checkBoxAnnotationOutline"));

        horizontalLayout22_3->addWidget(checkBoxAnnotationOutline);

        horizontalSpacer_30 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout22_3->addItem(horizontalSpacer_30);


        gridLayout->addLayout(horizontalLayout22_3, 14, 1, 1, 1);


        verticalLayout->addLayout(gridLayout);

        checkBoxShowInfo = new QCheckBox(scrollAreaWidgetContents);
        checkBoxShowInfo->setObjectName(QString::fromUtf8("checkBoxShowInfo"));

        verticalLayout->addWidget(checkBoxShowInfo);

        verticalSpacer = new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        PanelSurface->setWidget(scrollAreaWidgetContents);

        toolbar->addAction(actionSurfaceMain);
        toolbar->addAction(actionSurfaceInflated);
        toolbar->addAction(actionSurfaceWhite);
        toolbar->addAction(actionSurfacePial);
        toolbar->addAction(actionSurfaceOriginal);
        toolbar->addSeparator();
        toolbar->addAction(actionLockLayer);

        retranslateUi(PanelSurface);
        QObject::connect(sliderOpacity, SIGNAL(valueChanged(int)), PanelSurface, SLOT(OnSliderOpacity(int)));
        QObject::connect(comboBoxCurvature, SIGNAL(currentIndexChanged(int)), PanelSurface, SLOT(OnComboCurvature(int)));
        QObject::connect(sliderMidPoint, SIGNAL(valueChanged(int)), PanelSurface, SLOT(OnSliderMidPoint(int)));
        QObject::connect(sliderSlope, SIGNAL(valueChanged(int)), PanelSurface, SLOT(OnSliderSlope(int)));
        QObject::connect(lineEditMidPoint, SIGNAL(textEdited(QString)), PanelSurface, SLOT(OnLineEditMidPoint(QString)));
        QObject::connect(lineEditSlope, SIGNAL(textEdited(QString)), PanelSurface, SLOT(OnLineEditSlope(QString)));
        QObject::connect(comboBoxOverlay, SIGNAL(currentIndexChanged(int)), PanelSurface, SLOT(OnComboOverlay(int)));
        QObject::connect(pushButtonConfigureOverlay, SIGNAL(clicked()), PanelSurface, SLOT(OnButtonConfigureOverlay()));
        QObject::connect(comboBoxAnnotation, SIGNAL(currentIndexChanged(int)), PanelSurface, SLOT(OnComboAnnotation(int)));
        QObject::connect(comboBoxLabel, SIGNAL(currentIndexChanged(int)), PanelSurface, SLOT(OnComboLabel(int)));
        QObject::connect(comboBoxVectorDisplay, SIGNAL(currentIndexChanged(int)), PanelSurface, SLOT(OnComboVector(int)));
        QObject::connect(lineEditPositionOffset, SIGNAL(returnPressed()), PanelSurface, SLOT(OnEditPositionOffset()));

        QMetaObject::connectSlotsByName(PanelSurface);
    } // setupUi

    void retranslateUi(QScrollArea *PanelSurface)
    {
        PanelSurface->setWindowTitle(QApplication::translate("PanelSurface", "ScrollArea", 0, QApplication::UnicodeUTF8));
        actionSurfaceMain->setText(QApplication::translate("PanelSurface", "Main", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionSurfaceMain->setToolTip(QApplication::translate("PanelSurface", "Main", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionSurfaceInflated->setText(QApplication::translate("PanelSurface", "Inflated", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionSurfaceInflated->setToolTip(QApplication::translate("PanelSurface", "Inflated", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionSurfaceWhite->setText(QApplication::translate("PanelSurface", "White", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionSurfaceWhite->setToolTip(QApplication::translate("PanelSurface", "White", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionSurfacePial->setText(QApplication::translate("PanelSurface", "Pial", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionSurfacePial->setToolTip(QApplication::translate("PanelSurface", "Pial", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionSurfaceOriginal->setText(QApplication::translate("PanelSurface", "Original", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionSurfaceOriginal->setToolTip(QApplication::translate("PanelSurface", "Original", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionLockLayer->setText(QApplication::translate("PanelSurface", "Lock Layer", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionLockLayer->setToolTip(QApplication::translate("PanelSurface", "Lock", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        labelFileName->setText(QApplication::translate("PanelSurface", "File name", 0, QApplication::UnicodeUTF8));
        labelOpacity->setText(QApplication::translate("PanelSurface", "Opacity", 0, QApplication::UnicodeUTF8));
        labelSurfaceColor->setText(QApplication::translate("PanelSurface", "Color", 0, QApplication::UnicodeUTF8));
        labelRender->setText(QApplication::translate("PanelSurface", "Render", 0, QApplication::UnicodeUTF8));
        comboBoxRender->clear();
        comboBoxRender->insertItems(0, QStringList()
         << QApplication::translate("PanelSurface", "Surface", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelSurface", "Mesh", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelSurface", "Surface & mesh", 0, QApplication::UnicodeUTF8)
        );
        labelMeshColor->setText(QApplication::translate("PanelSurface", "Mesh color", 0, QApplication::UnicodeUTF8));
        comboBoxMeshColor->clear();
        comboBoxMeshColor->insertItems(0, QStringList()
         << QApplication::translate("PanelSurface", "Same as surface", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelSurface", "Curvature", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelSurface", "Overlay", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelSurface", "Solid color", 0, QApplication::UnicodeUTF8)
        );
        checkBoxShowVertices->setText(QApplication::translate("PanelSurface", "Show vertices", 0, QApplication::UnicodeUTF8));
        labelVertexColor->setText(QApplication::translate("PanelSurface", "Color", 0, QApplication::UnicodeUTF8));
        labelVertexPointSize->setText(QApplication::translate("PanelSurface", "Point size", 0, QApplication::UnicodeUTF8));
        labelCurvature->setText(QApplication::translate("PanelSurface", "Curvature", 0, QApplication::UnicodeUTF8));
        comboBoxCurvature->clear();
        comboBoxCurvature->insertItems(0, QStringList()
         << QApplication::translate("PanelSurface", "Off", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelSurface", "Threshold", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelSurface", "Binary", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelSurface", "Reload...", 0, QApplication::UnicodeUTF8)
        );
        labelMidPoint->setText(QApplication::translate("PanelSurface", "Mid point", 0, QApplication::UnicodeUTF8));
        labelSlope->setText(QApplication::translate("PanelSurface", "Slope", 0, QApplication::UnicodeUTF8));
        labelOverlay->setText(QApplication::translate("PanelSurface", "Overlay", 0, QApplication::UnicodeUTF8));
        pushButtonConfigureOverlay->setText(QApplication::translate("PanelSurface", "Configure Overlay", 0, QApplication::UnicodeUTF8));
        labelAnnotation->setText(QApplication::translate("PanelSurface", "Annotation", 0, QApplication::UnicodeUTF8));
        comboBoxAnnotation->clear();
        comboBoxAnnotation->insertItems(0, QStringList()
         << QApplication::translate("PanelSurface", "Off", 0, QApplication::UnicodeUTF8)
        );
        labelLabel->setText(QApplication::translate("PanelSurface", "Label", 0, QApplication::UnicodeUTF8));
        labelLabelColor->setText(QApplication::translate("PanelSurface", "Color", 0, QApplication::UnicodeUTF8));
        checkBoxLabelOutline->setText(QApplication::translate("PanelSurface", "Show outline only", 0, QApplication::UnicodeUTF8));
        labelEdgeColor->setText(QApplication::translate("PanelSurface", "Edge color", 0, QApplication::UnicodeUTF8));
        labelEdgeThickness->setText(QApplication::translate("PanelSurface", "Edge thickness", 0, QApplication::UnicodeUTF8));
        labelVectorDisplay->setText(QApplication::translate("PanelSurface", "Vector display", 0, QApplication::UnicodeUTF8));
        labelVectorColor->setText(QApplication::translate("PanelSurface", "Color", 0, QApplication::UnicodeUTF8));
        labelVectorPointSize->setText(QApplication::translate("PanelSurface", "Point size", 0, QApplication::UnicodeUTF8));
        labelPositionOffset->setText(QApplication::translate("PanelSurface", "Position offset", 0, QApplication::UnicodeUTF8));
        checkBoxAnnotationOutline->setText(QApplication::translate("PanelSurface", "Show outline only", 0, QApplication::UnicodeUTF8));
        checkBoxShowInfo->setText(QApplication::translate("PanelSurface", "Show in Info Panel", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class PanelSurface: public Ui_PanelSurface {};
} // namespace Ui

QT_END_NAMESPACE

