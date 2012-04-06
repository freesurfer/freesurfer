/********************************************************************************
** Form generated from reading UI file 'PanelVolume.ui'
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
#include <QtGui/QTreeWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "LayerTreeWidget.h"
#include "qtcolorpicker.h"

QT_BEGIN_NAMESPACE

class Ui_PanelVolume
{
public:
    QAction *actionMoveLayerUp;
    QAction *actionMoveLayerDown;
    QAction *actionLockLayer;
    QAction *actionCopySetting;
    QAction *actionPasteSetting;
    QAction *actionPasteSettingToAll;
    QWidget *scrollAreaWidgetContents;
    QVBoxLayout *verticalLayout;
    QVBoxLayout *verticalLayoutToolbar;
    QToolBar *toolbar;
    LayerTreeWidget *treeWidgetLayers;
    QToolBar *toolbar2;
    QGridLayout *gridLayout_2;
    QLabel *labelFileName;
    QLineEdit *lineEditFileName;
    QLabel *labelFrame;
    QHBoxLayout *horizontalLayout;
    QSlider *sliderFrame;
    QSpinBox *spinBoxFrame;
    QCheckBox *checkBoxDisplayVector;
    QCheckBox *checkBoxDisplayTensor;
    QLabel *labelRenderObject;
    QComboBox *comboBoxRenderObject;
    QLabel *labelInversion;
    QComboBox *comboBoxInversion;
    QLabel *labelMask;
    QComboBox *comboBoxMask;
    QLabel *labelOpacity;
    QHBoxLayout *horizontalLayout_2;
    QSlider *sliderOpacity;
    QDoubleSpinBox *doubleSpinBoxOpacity;
    QCheckBox *checkBoxSmooth;
    QCheckBox *checkBoxUpsample;
    QLabel *labelColorMap;
    QComboBox *comboBoxColorMap;
    QLabel *labelLookUpTable;
    QComboBox *comboBoxLookUpTable;
    QLabel *labelDirectionCode;
    QComboBox *comboBoxDirectionCode;
    QLabel *labelWindow;
    QHBoxLayout *horizontalLayout_6;
    QSlider *sliderWindow;
    QLineEdit *lineEditWindow;
    QLabel *labelLevel;
    QHBoxLayout *horizontalLayout_10;
    QSlider *sliderLevel;
    QLineEdit *lineEditLevel;
    QLabel *labelMin;
    QHBoxLayout *horizontalLayout_9;
    QSlider *sliderMin;
    QLineEdit *lineEditMin;
    QLabel *labelMid;
    QHBoxLayout *horizontalLayout_7;
    QSlider *sliderMid;
    QLineEdit *lineEditMid;
    QHBoxLayout *horizontalLayout_11;
    QSlider *sliderMax;
    QLineEdit *lineEditMax;
    QLabel *labelOffset;
    QHBoxLayout *horizontalLayout_14;
    QSlider *sliderOffset;
    QLineEdit *lineEditOffset;
    QGridLayout *gridLayout;
    QCheckBox *checkBoxTruncate;
    QSpacerItem *horizontalSpacer_2;
    QCheckBox *checkBoxClearBackground;
    QCheckBox *checkBoxInvert;
    QCheckBox *checkBoxClearHigher;
    QLabel *labelBrushValue;
    QHBoxLayout *horizontalLayout_17;
    QLineEdit *lineEditBrushValue;
    QLabel *colorLabelBrushValue;
    QSpacerItem *horizontalSpacer;
    QLabel *labelMax;
    QVBoxLayout *verticalLayout_2;
    QCheckBox *checkBoxRememberFrame;
    QLabel *labelRememberFrame;
    QCheckBox *checkBoxShowExistingLabels;
    QTreeWidget *treeWidgetColorTable;
    QGridLayout *gridLayout_3;
    QLabel *labelTrackVolumeThreshold;
    QHBoxLayout *horizontalLayout_19;
    QSlider *sliderTrackVolumeThresholdLow;
    QLineEdit *lineEditTrackVolumeThresholdLow;
    QCheckBox *checkBoxProjectionMap;
    QCheckBox *checkBoxShowOutline;
    QCheckBox *checkBoxShowContour;
    QGridLayout *gridLayoutContour;
    QLabel *labelContourThresholdLow;
    QCheckBox *checkBoxContourExtractAll;
    QLabel *labelContourColor;
    QtColorPicker *colorPickerContour;
    QCheckBox *checkBoxUseColorMap;
    QLabel *labelContourThresholdHigh;
    QHBoxLayout *horizontalLayout_13;
    QSlider *sliderContourThresholdHigh;
    QLineEdit *lineEditContourThresholdHigh;
    QLabel *labelSmoothIteration;
    QHBoxLayout *horizontalLayout_15;
    QSlider *sliderContourSmoothIteration;
    QLineEdit *lineEditContourSmoothIteration;
    QHBoxLayout *horizontalLayout_16;
    QPushButton *pushButtonContourSave;
    QSpacerItem *horizontalSpacer_3;
    QHBoxLayout *horizontalLayout_12;
    QSlider *sliderContourThresholdLow;
    QLineEdit *lineEditContourThresholdLow;
    QCheckBox *checkBoxShowInfo;
    QSpacerItem *verticalSpacer;

    void setupUi(QScrollArea *PanelVolume)
    {
        if (PanelVolume->objectName().isEmpty())
            PanelVolume->setObjectName(QString::fromUtf8("PanelVolume"));
        PanelVolume->resize(326, 1398);
        PanelVolume->setFrameShape(QFrame::NoFrame);
        PanelVolume->setFrameShadow(QFrame::Plain);
        PanelVolume->setWidgetResizable(true);
        actionMoveLayerUp = new QAction(PanelVolume);
        actionMoveLayerUp->setObjectName(QString::fromUtf8("actionMoveLayerUp"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/move_up.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionMoveLayerUp->setIcon(icon);
        actionMoveLayerDown = new QAction(PanelVolume);
        actionMoveLayerDown->setObjectName(QString::fromUtf8("actionMoveLayerDown"));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/resource/icons/move_down.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionMoveLayerDown->setIcon(icon1);
        actionLockLayer = new QAction(PanelVolume);
        actionLockLayer->setObjectName(QString::fromUtf8("actionLockLayer"));
        actionLockLayer->setCheckable(true);
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/resource/icons/volume_lock.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLockLayer->setIcon(icon2);
        actionCopySetting = new QAction(PanelVolume);
        actionCopySetting->setObjectName(QString::fromUtf8("actionCopySetting"));
        actionCopySetting->setEnabled(false);
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/resource/icons/volume_copy_setting.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionCopySetting->setIcon(icon3);
        actionPasteSetting = new QAction(PanelVolume);
        actionPasteSetting->setObjectName(QString::fromUtf8("actionPasteSetting"));
        actionPasteSetting->setEnabled(false);
        QIcon icon4;
        icon4.addFile(QString::fromUtf8(":/resource/icons/volume_paste_setting.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionPasteSetting->setIcon(icon4);
        actionPasteSettingToAll = new QAction(PanelVolume);
        actionPasteSettingToAll->setObjectName(QString::fromUtf8("actionPasteSettingToAll"));
        actionPasteSettingToAll->setEnabled(false);
        QIcon icon5;
        icon5.addFile(QString::fromUtf8(":/resource/icons/volume_paste_setting_all.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionPasteSettingToAll->setIcon(icon5);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QString::fromUtf8("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 308, 1632));
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
        treeWidgetLayers->setEditTriggers(QAbstractItemView::EditKeyPressed|QAbstractItemView::SelectedClicked);
        treeWidgetLayers->setSelectionMode(QAbstractItemView::SingleSelection);
        treeWidgetLayers->setRootIsDecorated(false);
        treeWidgetLayers->setUniformRowHeights(true);
        treeWidgetLayers->header()->setVisible(false);

        verticalLayoutToolbar->addWidget(treeWidgetLayers);

        toolbar2 = new QToolBar(scrollAreaWidgetContents);
        toolbar2->setObjectName(QString::fromUtf8("toolbar2"));
        toolbar2->setIconSize(QSize(24, 24));
        toolbar2->setFloatable(false);

        verticalLayoutToolbar->addWidget(toolbar2);


        verticalLayout->addLayout(verticalLayoutToolbar);

        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        labelFileName = new QLabel(scrollAreaWidgetContents);
        labelFileName->setObjectName(QString::fromUtf8("labelFileName"));

        gridLayout_2->addWidget(labelFileName, 0, 0, 1, 1);

        lineEditFileName = new QLineEdit(scrollAreaWidgetContents);
        lineEditFileName->setObjectName(QString::fromUtf8("lineEditFileName"));
        lineEditFileName->setReadOnly(true);

        gridLayout_2->addWidget(lineEditFileName, 0, 1, 1, 1);

        labelFrame = new QLabel(scrollAreaWidgetContents);
        labelFrame->setObjectName(QString::fromUtf8("labelFrame"));

        gridLayout_2->addWidget(labelFrame, 1, 0, 1, 1);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        sliderFrame = new QSlider(scrollAreaWidgetContents);
        sliderFrame->setObjectName(QString::fromUtf8("sliderFrame"));
        sliderFrame->setMaximum(100);
        sliderFrame->setOrientation(Qt::Horizontal);

        horizontalLayout->addWidget(sliderFrame);

        spinBoxFrame = new QSpinBox(scrollAreaWidgetContents);
        spinBoxFrame->setObjectName(QString::fromUtf8("spinBoxFrame"));
        spinBoxFrame->setMinimum(1);

        horizontalLayout->addWidget(spinBoxFrame);


        gridLayout_2->addLayout(horizontalLayout, 1, 1, 1, 1);

        checkBoxDisplayVector = new QCheckBox(scrollAreaWidgetContents);
        checkBoxDisplayVector->setObjectName(QString::fromUtf8("checkBoxDisplayVector"));

        gridLayout_2->addWidget(checkBoxDisplayVector, 3, 1, 1, 1);

        checkBoxDisplayTensor = new QCheckBox(scrollAreaWidgetContents);
        checkBoxDisplayTensor->setObjectName(QString::fromUtf8("checkBoxDisplayTensor"));

        gridLayout_2->addWidget(checkBoxDisplayTensor, 4, 1, 1, 1);

        labelRenderObject = new QLabel(scrollAreaWidgetContents);
        labelRenderObject->setObjectName(QString::fromUtf8("labelRenderObject"));

        gridLayout_2->addWidget(labelRenderObject, 5, 0, 1, 1);

        comboBoxRenderObject = new QComboBox(scrollAreaWidgetContents);
        comboBoxRenderObject->setObjectName(QString::fromUtf8("comboBoxRenderObject"));

        gridLayout_2->addWidget(comboBoxRenderObject, 5, 1, 1, 1);

        labelInversion = new QLabel(scrollAreaWidgetContents);
        labelInversion->setObjectName(QString::fromUtf8("labelInversion"));

        gridLayout_2->addWidget(labelInversion, 6, 0, 1, 1);

        comboBoxInversion = new QComboBox(scrollAreaWidgetContents);
        comboBoxInversion->setObjectName(QString::fromUtf8("comboBoxInversion"));

        gridLayout_2->addWidget(comboBoxInversion, 6, 1, 1, 1);

        labelMask = new QLabel(scrollAreaWidgetContents);
        labelMask->setObjectName(QString::fromUtf8("labelMask"));

        gridLayout_2->addWidget(labelMask, 7, 0, 1, 1);

        comboBoxMask = new QComboBox(scrollAreaWidgetContents);
        comboBoxMask->setObjectName(QString::fromUtf8("comboBoxMask"));

        gridLayout_2->addWidget(comboBoxMask, 7, 1, 1, 1);

        labelOpacity = new QLabel(scrollAreaWidgetContents);
        labelOpacity->setObjectName(QString::fromUtf8("labelOpacity"));

        gridLayout_2->addWidget(labelOpacity, 8, 0, 1, 1);

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


        gridLayout_2->addLayout(horizontalLayout_2, 8, 1, 1, 1);

        checkBoxSmooth = new QCheckBox(scrollAreaWidgetContents);
        checkBoxSmooth->setObjectName(QString::fromUtf8("checkBoxSmooth"));

        gridLayout_2->addWidget(checkBoxSmooth, 9, 1, 1, 1);

        checkBoxUpsample = new QCheckBox(scrollAreaWidgetContents);
        checkBoxUpsample->setObjectName(QString::fromUtf8("checkBoxUpsample"));

        gridLayout_2->addWidget(checkBoxUpsample, 10, 1, 1, 1);

        labelColorMap = new QLabel(scrollAreaWidgetContents);
        labelColorMap->setObjectName(QString::fromUtf8("labelColorMap"));

        gridLayout_2->addWidget(labelColorMap, 11, 0, 1, 1);

        comboBoxColorMap = new QComboBox(scrollAreaWidgetContents);
        comboBoxColorMap->setObjectName(QString::fromUtf8("comboBoxColorMap"));

        gridLayout_2->addWidget(comboBoxColorMap, 11, 1, 1, 1);

        labelLookUpTable = new QLabel(scrollAreaWidgetContents);
        labelLookUpTable->setObjectName(QString::fromUtf8("labelLookUpTable"));

        gridLayout_2->addWidget(labelLookUpTable, 12, 0, 1, 1);

        comboBoxLookUpTable = new QComboBox(scrollAreaWidgetContents);
        comboBoxLookUpTable->setObjectName(QString::fromUtf8("comboBoxLookUpTable"));

        gridLayout_2->addWidget(comboBoxLookUpTable, 12, 1, 1, 1);

        labelDirectionCode = new QLabel(scrollAreaWidgetContents);
        labelDirectionCode->setObjectName(QString::fromUtf8("labelDirectionCode"));

        gridLayout_2->addWidget(labelDirectionCode, 13, 0, 1, 1);

        comboBoxDirectionCode = new QComboBox(scrollAreaWidgetContents);
        comboBoxDirectionCode->setObjectName(QString::fromUtf8("comboBoxDirectionCode"));

        gridLayout_2->addWidget(comboBoxDirectionCode, 13, 1, 1, 1);

        labelWindow = new QLabel(scrollAreaWidgetContents);
        labelWindow->setObjectName(QString::fromUtf8("labelWindow"));

        gridLayout_2->addWidget(labelWindow, 14, 0, 1, 1);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        sliderWindow = new QSlider(scrollAreaWidgetContents);
        sliderWindow->setObjectName(QString::fromUtf8("sliderWindow"));
        sliderWindow->setMaximum(100);
        sliderWindow->setOrientation(Qt::Horizontal);

        horizontalLayout_6->addWidget(sliderWindow);

        lineEditWindow = new QLineEdit(scrollAreaWidgetContents);
        lineEditWindow->setObjectName(QString::fromUtf8("lineEditWindow"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(lineEditWindow->sizePolicy().hasHeightForWidth());
        lineEditWindow->setSizePolicy(sizePolicy1);
        lineEditWindow->setMaximumSize(QSize(75, 16777215));

        horizontalLayout_6->addWidget(lineEditWindow);


        gridLayout_2->addLayout(horizontalLayout_6, 14, 1, 1, 1);

        labelLevel = new QLabel(scrollAreaWidgetContents);
        labelLevel->setObjectName(QString::fromUtf8("labelLevel"));

        gridLayout_2->addWidget(labelLevel, 15, 0, 1, 1);

        horizontalLayout_10 = new QHBoxLayout();
        horizontalLayout_10->setObjectName(QString::fromUtf8("horizontalLayout_10"));
        sliderLevel = new QSlider(scrollAreaWidgetContents);
        sliderLevel->setObjectName(QString::fromUtf8("sliderLevel"));
        sliderLevel->setMaximum(100);
        sliderLevel->setOrientation(Qt::Horizontal);

        horizontalLayout_10->addWidget(sliderLevel);

        lineEditLevel = new QLineEdit(scrollAreaWidgetContents);
        lineEditLevel->setObjectName(QString::fromUtf8("lineEditLevel"));
        sizePolicy1.setHeightForWidth(lineEditLevel->sizePolicy().hasHeightForWidth());
        lineEditLevel->setSizePolicy(sizePolicy1);
        lineEditLevel->setMaximumSize(QSize(75, 16777215));

        horizontalLayout_10->addWidget(lineEditLevel);


        gridLayout_2->addLayout(horizontalLayout_10, 15, 1, 1, 1);

        labelMin = new QLabel(scrollAreaWidgetContents);
        labelMin->setObjectName(QString::fromUtf8("labelMin"));

        gridLayout_2->addWidget(labelMin, 16, 0, 1, 1);

        horizontalLayout_9 = new QHBoxLayout();
        horizontalLayout_9->setObjectName(QString::fromUtf8("horizontalLayout_9"));
        sliderMin = new QSlider(scrollAreaWidgetContents);
        sliderMin->setObjectName(QString::fromUtf8("sliderMin"));
        sliderMin->setMaximum(100);
        sliderMin->setOrientation(Qt::Horizontal);

        horizontalLayout_9->addWidget(sliderMin);

        lineEditMin = new QLineEdit(scrollAreaWidgetContents);
        lineEditMin->setObjectName(QString::fromUtf8("lineEditMin"));
        sizePolicy1.setHeightForWidth(lineEditMin->sizePolicy().hasHeightForWidth());
        lineEditMin->setSizePolicy(sizePolicy1);
        lineEditMin->setMaximumSize(QSize(75, 16777215));

        horizontalLayout_9->addWidget(lineEditMin);


        gridLayout_2->addLayout(horizontalLayout_9, 16, 1, 1, 1);

        labelMid = new QLabel(scrollAreaWidgetContents);
        labelMid->setObjectName(QString::fromUtf8("labelMid"));

        gridLayout_2->addWidget(labelMid, 17, 0, 1, 1);

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
        lineEditMid->setMaximumSize(QSize(75, 16777215));

        horizontalLayout_7->addWidget(lineEditMid);


        gridLayout_2->addLayout(horizontalLayout_7, 17, 1, 1, 1);

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
        lineEditMax->setMaximumSize(QSize(75, 16777215));

        horizontalLayout_11->addWidget(lineEditMax);


        gridLayout_2->addLayout(horizontalLayout_11, 18, 1, 1, 1);

        labelOffset = new QLabel(scrollAreaWidgetContents);
        labelOffset->setObjectName(QString::fromUtf8("labelOffset"));

        gridLayout_2->addWidget(labelOffset, 19, 0, 1, 1);

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
        lineEditOffset->setMaximumSize(QSize(75, 16777215));

        horizontalLayout_14->addWidget(lineEditOffset);


        gridLayout_2->addLayout(horizontalLayout_14, 19, 1, 1, 1);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        checkBoxTruncate = new QCheckBox(scrollAreaWidgetContents);
        checkBoxTruncate->setObjectName(QString::fromUtf8("checkBoxTruncate"));

        gridLayout->addWidget(checkBoxTruncate, 0, 0, 1, 1);

        horizontalSpacer_2 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_2, 0, 2, 1, 1);

        checkBoxClearBackground = new QCheckBox(scrollAreaWidgetContents);
        checkBoxClearBackground->setObjectName(QString::fromUtf8("checkBoxClearBackground"));

        gridLayout->addWidget(checkBoxClearBackground, 2, 0, 1, 2);

        checkBoxInvert = new QCheckBox(scrollAreaWidgetContents);
        checkBoxInvert->setObjectName(QString::fromUtf8("checkBoxInvert"));

        gridLayout->addWidget(checkBoxInvert, 0, 1, 1, 1);

        checkBoxClearHigher = new QCheckBox(scrollAreaWidgetContents);
        checkBoxClearHigher->setObjectName(QString::fromUtf8("checkBoxClearHigher"));

        gridLayout->addWidget(checkBoxClearHigher, 1, 0, 1, 2);


        gridLayout_2->addLayout(gridLayout, 20, 1, 1, 1);

        labelBrushValue = new QLabel(scrollAreaWidgetContents);
        labelBrushValue->setObjectName(QString::fromUtf8("labelBrushValue"));

        gridLayout_2->addWidget(labelBrushValue, 21, 0, 1, 1);

        horizontalLayout_17 = new QHBoxLayout();
        horizontalLayout_17->setObjectName(QString::fromUtf8("horizontalLayout_17"));
        lineEditBrushValue = new QLineEdit(scrollAreaWidgetContents);
        lineEditBrushValue->setObjectName(QString::fromUtf8("lineEditBrushValue"));
        lineEditBrushValue->setMaximumSize(QSize(60, 16777215));
        lineEditBrushValue->setFocusPolicy(Qt::WheelFocus);

        horizontalLayout_17->addWidget(lineEditBrushValue);

        colorLabelBrushValue = new QLabel(scrollAreaWidgetContents);
        colorLabelBrushValue->setObjectName(QString::fromUtf8("colorLabelBrushValue"));
        colorLabelBrushValue->setMinimumSize(QSize(40, 0));
        colorLabelBrushValue->setMargin(3);

        horizontalLayout_17->addWidget(colorLabelBrushValue);

        horizontalSpacer = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_17->addItem(horizontalSpacer);


        gridLayout_2->addLayout(horizontalLayout_17, 21, 1, 1, 1);

        labelMax = new QLabel(scrollAreaWidgetContents);
        labelMax->setObjectName(QString::fromUtf8("labelMax"));

        gridLayout_2->addWidget(labelMax, 18, 0, 1, 1);

        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setSpacing(0);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        checkBoxRememberFrame = new QCheckBox(scrollAreaWidgetContents);
        checkBoxRememberFrame->setObjectName(QString::fromUtf8("checkBoxRememberFrame"));

        verticalLayout_2->addWidget(checkBoxRememberFrame);

        labelRememberFrame = new QLabel(scrollAreaWidgetContents);
        labelRememberFrame->setObjectName(QString::fromUtf8("labelRememberFrame"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy2.setHorizontalStretch(1);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(labelRememberFrame->sizePolicy().hasHeightForWidth());
        labelRememberFrame->setSizePolicy(sizePolicy2);
        labelRememberFrame->setWordWrap(true);
        labelRememberFrame->setIndent(21);

        verticalLayout_2->addWidget(labelRememberFrame);


        gridLayout_2->addLayout(verticalLayout_2, 2, 1, 1, 1);


        verticalLayout->addLayout(gridLayout_2);

        checkBoxShowExistingLabels = new QCheckBox(scrollAreaWidgetContents);
        checkBoxShowExistingLabels->setObjectName(QString::fromUtf8("checkBoxShowExistingLabels"));

        verticalLayout->addWidget(checkBoxShowExistingLabels);

        treeWidgetColorTable = new QTreeWidget(scrollAreaWidgetContents);
        QTreeWidgetItem *__qtreewidgetitem1 = new QTreeWidgetItem();
        __qtreewidgetitem1->setText(0, QString::fromUtf8("1"));
        treeWidgetColorTable->setHeaderItem(__qtreewidgetitem1);
        treeWidgetColorTable->setObjectName(QString::fromUtf8("treeWidgetColorTable"));
        QSizePolicy sizePolicy3(QSizePolicy::Expanding, QSizePolicy::Expanding);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(1);
        sizePolicy3.setHeightForWidth(treeWidgetColorTable->sizePolicy().hasHeightForWidth());
        treeWidgetColorTable->setSizePolicy(sizePolicy3);
        treeWidgetColorTable->setMinimumSize(QSize(0, 150));
        treeWidgetColorTable->setAlternatingRowColors(false);
        treeWidgetColorTable->setRootIsDecorated(false);
        treeWidgetColorTable->setUniformRowHeights(true);
        treeWidgetColorTable->header()->setVisible(false);

        verticalLayout->addWidget(treeWidgetColorTable);

        gridLayout_3 = new QGridLayout();
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        gridLayout_3->setContentsMargins(-1, 0, -1, -1);
        labelTrackVolumeThreshold = new QLabel(scrollAreaWidgetContents);
        labelTrackVolumeThreshold->setObjectName(QString::fromUtf8("labelTrackVolumeThreshold"));

        gridLayout_3->addWidget(labelTrackVolumeThreshold, 0, 0, 1, 1);

        horizontalLayout_19 = new QHBoxLayout();
        horizontalLayout_19->setObjectName(QString::fromUtf8("horizontalLayout_19"));
        sliderTrackVolumeThresholdLow = new QSlider(scrollAreaWidgetContents);
        sliderTrackVolumeThresholdLow->setObjectName(QString::fromUtf8("sliderTrackVolumeThresholdLow"));
        sliderTrackVolumeThresholdLow->setMaximum(100);
        sliderTrackVolumeThresholdLow->setOrientation(Qt::Horizontal);

        horizontalLayout_19->addWidget(sliderTrackVolumeThresholdLow);

        lineEditTrackVolumeThresholdLow = new QLineEdit(scrollAreaWidgetContents);
        lineEditTrackVolumeThresholdLow->setObjectName(QString::fromUtf8("lineEditTrackVolumeThresholdLow"));
        sizePolicy1.setHeightForWidth(lineEditTrackVolumeThresholdLow->sizePolicy().hasHeightForWidth());
        lineEditTrackVolumeThresholdLow->setSizePolicy(sizePolicy1);
        lineEditTrackVolumeThresholdLow->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_19->addWidget(lineEditTrackVolumeThresholdLow);


        gridLayout_3->addLayout(horizontalLayout_19, 0, 1, 1, 1);


        verticalLayout->addLayout(gridLayout_3);

        checkBoxProjectionMap = new QCheckBox(scrollAreaWidgetContents);
        checkBoxProjectionMap->setObjectName(QString::fromUtf8("checkBoxProjectionMap"));

        verticalLayout->addWidget(checkBoxProjectionMap);

        checkBoxShowOutline = new QCheckBox(scrollAreaWidgetContents);
        checkBoxShowOutline->setObjectName(QString::fromUtf8("checkBoxShowOutline"));

        verticalLayout->addWidget(checkBoxShowOutline);

        checkBoxShowContour = new QCheckBox(scrollAreaWidgetContents);
        checkBoxShowContour->setObjectName(QString::fromUtf8("checkBoxShowContour"));

        verticalLayout->addWidget(checkBoxShowContour);

        gridLayoutContour = new QGridLayout();
        gridLayoutContour->setObjectName(QString::fromUtf8("gridLayoutContour"));
        labelContourThresholdLow = new QLabel(scrollAreaWidgetContents);
        labelContourThresholdLow->setObjectName(QString::fromUtf8("labelContourThresholdLow"));

        gridLayoutContour->addWidget(labelContourThresholdLow, 0, 0, 1, 1);

        checkBoxContourExtractAll = new QCheckBox(scrollAreaWidgetContents);
        checkBoxContourExtractAll->setObjectName(QString::fromUtf8("checkBoxContourExtractAll"));

        gridLayoutContour->addWidget(checkBoxContourExtractAll, 2, 1, 1, 1);

        labelContourColor = new QLabel(scrollAreaWidgetContents);
        labelContourColor->setObjectName(QString::fromUtf8("labelContourColor"));
        labelContourColor->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayoutContour->addWidget(labelContourColor, 4, 0, 1, 1);

        colorPickerContour = new QtColorPicker(scrollAreaWidgetContents);
        colorPickerContour->setObjectName(QString::fromUtf8("colorPickerContour"));

        gridLayoutContour->addWidget(colorPickerContour, 4, 1, 1, 1);

        checkBoxUseColorMap = new QCheckBox(scrollAreaWidgetContents);
        checkBoxUseColorMap->setObjectName(QString::fromUtf8("checkBoxUseColorMap"));

        gridLayoutContour->addWidget(checkBoxUseColorMap, 5, 1, 1, 1);

        labelContourThresholdHigh = new QLabel(scrollAreaWidgetContents);
        labelContourThresholdHigh->setObjectName(QString::fromUtf8("labelContourThresholdHigh"));

        gridLayoutContour->addWidget(labelContourThresholdHigh, 1, 0, 1, 1);

        horizontalLayout_13 = new QHBoxLayout();
        horizontalLayout_13->setObjectName(QString::fromUtf8("horizontalLayout_13"));
        sliderContourThresholdHigh = new QSlider(scrollAreaWidgetContents);
        sliderContourThresholdHigh->setObjectName(QString::fromUtf8("sliderContourThresholdHigh"));
        sliderContourThresholdHigh->setMaximum(100);
        sliderContourThresholdHigh->setOrientation(Qt::Horizontal);

        horizontalLayout_13->addWidget(sliderContourThresholdHigh);

        lineEditContourThresholdHigh = new QLineEdit(scrollAreaWidgetContents);
        lineEditContourThresholdHigh->setObjectName(QString::fromUtf8("lineEditContourThresholdHigh"));
        sizePolicy1.setHeightForWidth(lineEditContourThresholdHigh->sizePolicy().hasHeightForWidth());
        lineEditContourThresholdHigh->setSizePolicy(sizePolicy1);
        lineEditContourThresholdHigh->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_13->addWidget(lineEditContourThresholdHigh);


        gridLayoutContour->addLayout(horizontalLayout_13, 1, 1, 1, 1);

        labelSmoothIteration = new QLabel(scrollAreaWidgetContents);
        labelSmoothIteration->setObjectName(QString::fromUtf8("labelSmoothIteration"));

        gridLayoutContour->addWidget(labelSmoothIteration, 3, 0, 1, 1);

        horizontalLayout_15 = new QHBoxLayout();
        horizontalLayout_15->setObjectName(QString::fromUtf8("horizontalLayout_15"));
        sliderContourSmoothIteration = new QSlider(scrollAreaWidgetContents);
        sliderContourSmoothIteration->setObjectName(QString::fromUtf8("sliderContourSmoothIteration"));
        sliderContourSmoothIteration->setMaximum(20);
        sliderContourSmoothIteration->setOrientation(Qt::Horizontal);

        horizontalLayout_15->addWidget(sliderContourSmoothIteration);

        lineEditContourSmoothIteration = new QLineEdit(scrollAreaWidgetContents);
        lineEditContourSmoothIteration->setObjectName(QString::fromUtf8("lineEditContourSmoothIteration"));
        sizePolicy1.setHeightForWidth(lineEditContourSmoothIteration->sizePolicy().hasHeightForWidth());
        lineEditContourSmoothIteration->setSizePolicy(sizePolicy1);
        lineEditContourSmoothIteration->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_15->addWidget(lineEditContourSmoothIteration);


        gridLayoutContour->addLayout(horizontalLayout_15, 3, 1, 1, 1);

        horizontalLayout_16 = new QHBoxLayout();
        horizontalLayout_16->setObjectName(QString::fromUtf8("horizontalLayout_16"));
        pushButtonContourSave = new QPushButton(scrollAreaWidgetContents);
        pushButtonContourSave->setObjectName(QString::fromUtf8("pushButtonContourSave"));

        horizontalLayout_16->addWidget(pushButtonContourSave);

        horizontalSpacer_3 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_16->addItem(horizontalSpacer_3);


        gridLayoutContour->addLayout(horizontalLayout_16, 6, 1, 1, 1);

        horizontalLayout_12 = new QHBoxLayout();
        horizontalLayout_12->setObjectName(QString::fromUtf8("horizontalLayout_12"));
        sliderContourThresholdLow = new QSlider(scrollAreaWidgetContents);
        sliderContourThresholdLow->setObjectName(QString::fromUtf8("sliderContourThresholdLow"));
        sliderContourThresholdLow->setMaximum(100);
        sliderContourThresholdLow->setOrientation(Qt::Horizontal);

        horizontalLayout_12->addWidget(sliderContourThresholdLow);

        lineEditContourThresholdLow = new QLineEdit(scrollAreaWidgetContents);
        lineEditContourThresholdLow->setObjectName(QString::fromUtf8("lineEditContourThresholdLow"));
        sizePolicy1.setHeightForWidth(lineEditContourThresholdLow->sizePolicy().hasHeightForWidth());
        lineEditContourThresholdLow->setSizePolicy(sizePolicy1);
        lineEditContourThresholdLow->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_12->addWidget(lineEditContourThresholdLow);


        gridLayoutContour->addLayout(horizontalLayout_12, 0, 1, 1, 1);


        verticalLayout->addLayout(gridLayoutContour);

        checkBoxShowInfo = new QCheckBox(scrollAreaWidgetContents);
        checkBoxShowInfo->setObjectName(QString::fromUtf8("checkBoxShowInfo"));

        verticalLayout->addWidget(checkBoxShowInfo);

        verticalSpacer = new QSpacerItem(20, 5, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        PanelVolume->setWidget(scrollAreaWidgetContents);

        toolbar->addAction(actionMoveLayerUp);
        toolbar->addAction(actionMoveLayerDown);
        toolbar->addAction(actionLockLayer);
        toolbar2->addAction(actionCopySetting);
        toolbar2->addAction(actionPasteSetting);
        toolbar2->addAction(actionPasteSettingToAll);

        retranslateUi(PanelVolume);
        QObject::connect(checkBoxShowContour, SIGNAL(toggled(bool)), PanelVolume, SLOT(OnCheckShowContour(bool)));
        QObject::connect(comboBoxColorMap, SIGNAL(currentIndexChanged(int)), PanelVolume, SLOT(OnComboColorMap(int)));
        QObject::connect(treeWidgetColorTable, SIGNAL(currentItemChanged(QTreeWidgetItem*,QTreeWidgetItem*)), PanelVolume, SLOT(OnColorTableCurrentItemChanged(QTreeWidgetItem*)));
        QObject::connect(lineEditBrushValue, SIGNAL(textEdited(QString)), PanelVolume, SLOT(OnLineEditBrushValue(QString)));
        QObject::connect(sliderOpacity, SIGNAL(valueChanged(int)), PanelVolume, SLOT(OnSliderOpacity(int)));
        QObject::connect(sliderWindow, SIGNAL(valueChanged(int)), PanelVolume, SLOT(OnSliderWindow(int)));
        QObject::connect(sliderLevel, SIGNAL(valueChanged(int)), PanelVolume, SLOT(OnSliderLevel(int)));
        QObject::connect(sliderMin, SIGNAL(valueChanged(int)), PanelVolume, SLOT(OnSliderMin(int)));
        QObject::connect(sliderMid, SIGNAL(valueChanged(int)), PanelVolume, SLOT(OnSliderMid(int)));
        QObject::connect(sliderMax, SIGNAL(valueChanged(int)), PanelVolume, SLOT(OnSliderMax(int)));
        QObject::connect(sliderOffset, SIGNAL(valueChanged(int)), PanelVolume, SLOT(OnSliderOffset(int)));
        QObject::connect(lineEditWindow, SIGNAL(textEdited(QString)), PanelVolume, SLOT(OnLineEditWindow(QString)));
        QObject::connect(lineEditLevel, SIGNAL(textEdited(QString)), PanelVolume, SLOT(OnLineEditLevel(QString)));
        QObject::connect(lineEditMin, SIGNAL(textEdited(QString)), PanelVolume, SLOT(OnLineEditMin(QString)));
        QObject::connect(lineEditMid, SIGNAL(textEdited(QString)), PanelVolume, SLOT(OnLineEditMid(QString)));
        QObject::connect(lineEditMax, SIGNAL(textEdited(QString)), PanelVolume, SLOT(OnLineEditMax(QString)));
        QObject::connect(lineEditOffset, SIGNAL(textEdited(QString)), PanelVolume, SLOT(OnLineEditOffset(QString)));
        QObject::connect(comboBoxLookUpTable, SIGNAL(currentIndexChanged(int)), PanelVolume, SLOT(OnComboLookupTable(int)));
        QObject::connect(sliderContourThresholdLow, SIGNAL(valueChanged(int)), PanelVolume, SLOT(OnSliderContourMin(int)));
        QObject::connect(sliderContourThresholdHigh, SIGNAL(valueChanged(int)), PanelVolume, SLOT(OnSliderContourMax(int)));
        QObject::connect(sliderContourSmoothIteration, SIGNAL(valueChanged(int)), PanelVolume, SLOT(OnSliderContourSmooth(int)));
        QObject::connect(sliderContourThresholdLow, SIGNAL(sliderReleased()), PanelVolume, SLOT(OnContourValueChanged()));
        QObject::connect(sliderContourThresholdHigh, SIGNAL(sliderReleased()), PanelVolume, SLOT(OnContourValueChanged()));
        QObject::connect(sliderContourSmoothIteration, SIGNAL(sliderReleased()), PanelVolume, SLOT(OnContourValueChanged()));
        QObject::connect(lineEditContourThresholdLow, SIGNAL(returnPressed()), PanelVolume, SLOT(OnContourValueChanged()));
        QObject::connect(lineEditContourThresholdHigh, SIGNAL(returnPressed()), PanelVolume, SLOT(OnContourValueChanged()));
        QObject::connect(lineEditContourSmoothIteration, SIGNAL(returnPressed()), PanelVolume, SLOT(OnContourValueChanged()));
        QObject::connect(pushButtonContourSave, SIGNAL(clicked()), PanelVolume, SLOT(OnContourSave()));
        QObject::connect(actionCopySetting, SIGNAL(triggered()), PanelVolume, SLOT(OnCopySettings()));
        QObject::connect(actionPasteSetting, SIGNAL(triggered()), PanelVolume, SLOT(OnPasteSettings()));
        QObject::connect(actionPasteSettingToAll, SIGNAL(triggered()), PanelVolume, SLOT(OnPasteSettingsToAll()));
        QObject::connect(sliderTrackVolumeThresholdLow, SIGNAL(valueChanged(int)), PanelVolume, SLOT(OnSliderTrackVolumeMin(int)));
        QObject::connect(lineEditTrackVolumeThresholdLow, SIGNAL(returnPressed()), PanelVolume, SLOT(OnTrackVolumeThresholdChanged()));
        QObject::connect(sliderTrackVolumeThresholdLow, SIGNAL(sliderReleased()), PanelVolume, SLOT(OnTrackVolumeThresholdChanged()));
        QObject::connect(treeWidgetColorTable, SIGNAL(itemClicked(QTreeWidgetItem*,int)), PanelVolume, SLOT(UpdateTrackVolumeThreshold()));
        QObject::connect(checkBoxProjectionMap, SIGNAL(toggled(bool)), checkBoxShowOutline, SLOT(setHidden(bool)));
        QObject::connect(checkBoxProjectionMap, SIGNAL(toggled(bool)), checkBoxShowContour, SLOT(setHidden(bool)));
        QObject::connect(checkBoxShowContour, SIGNAL(toggled(bool)), checkBoxProjectionMap, SLOT(setDisabled(bool)));
        QObject::connect(checkBoxShowExistingLabels, SIGNAL(toggled(bool)), PanelVolume, SLOT(OnShowExistingLabelsOnly(bool)));

        QMetaObject::connectSlotsByName(PanelVolume);
    } // setupUi

    void retranslateUi(QScrollArea *PanelVolume)
    {
        PanelVolume->setWindowTitle(QApplication::translate("PanelVolume", "ScrollArea", 0, QApplication::UnicodeUTF8));
        actionMoveLayerUp->setText(QApplication::translate("PanelVolume", "Move Layer Up", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionMoveLayerUp->setToolTip(QApplication::translate("PanelVolume", "Up", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionMoveLayerDown->setText(QApplication::translate("PanelVolume", "Move Layer Down", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionMoveLayerDown->setToolTip(QApplication::translate("PanelVolume", "Down", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionLockLayer->setText(QApplication::translate("PanelVolume", "Lock Layer", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionLockLayer->setToolTip(QApplication::translate("PanelVolume", "Lock", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionCopySetting->setText(QApplication::translate("PanelVolume", "Copy Setting", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionCopySetting->setToolTip(QApplication::translate("PanelVolume", "Copy Setting", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionPasteSetting->setText(QApplication::translate("PanelVolume", "Paste Setting", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionPasteSetting->setToolTip(QApplication::translate("PanelVolume", "Paste Setting", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionPasteSettingToAll->setText(QApplication::translate("PanelVolume", "Paste Setting To All", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionPasteSettingToAll->setToolTip(QApplication::translate("PanelVolume", "Paste Setting To All", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        labelFileName->setText(QApplication::translate("PanelVolume", "File name", 0, QApplication::UnicodeUTF8));
        labelFrame->setText(QApplication::translate("PanelVolume", "Frame", 0, QApplication::UnicodeUTF8));
        checkBoxDisplayVector->setText(QApplication::translate("PanelVolume", "Display as vectors", 0, QApplication::UnicodeUTF8));
        checkBoxDisplayTensor->setText(QApplication::translate("PanelVolume", "Display as tensors", 0, QApplication::UnicodeUTF8));
        labelRenderObject->setText(QApplication::translate("PanelVolume", "Render", 0, QApplication::UnicodeUTF8));
        comboBoxRenderObject->clear();
        comboBoxRenderObject->insertItems(0, QStringList()
         << QApplication::translate("PanelVolume", "Simple line", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "3D bar (slow!)", 0, QApplication::UnicodeUTF8)
        );
        labelInversion->setText(QApplication::translate("PanelVolume", "Inversion", 0, QApplication::UnicodeUTF8));
        comboBoxInversion->clear();
        comboBoxInversion->insertItems(0, QStringList()
         << QApplication::translate("PanelVolume", "None", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "Invert X", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "Invert Y", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "Invert Z", 0, QApplication::UnicodeUTF8)
        );
        labelMask->setText(QApplication::translate("PanelVolume", "Mask", 0, QApplication::UnicodeUTF8));
        comboBoxMask->clear();
        comboBoxMask->insertItems(0, QStringList()
         << QApplication::translate("PanelVolume", "None", 0, QApplication::UnicodeUTF8)
        );
        labelOpacity->setText(QApplication::translate("PanelVolume", "Opacity", 0, QApplication::UnicodeUTF8));
        checkBoxSmooth->setText(QApplication::translate("PanelVolume", "Smooth display", 0, QApplication::UnicodeUTF8));
        checkBoxUpsample->setText(QApplication::translate("PanelVolume", "Upsample display", 0, QApplication::UnicodeUTF8));
        labelColorMap->setText(QApplication::translate("PanelVolume", "Color map", 0, QApplication::UnicodeUTF8));
        comboBoxColorMap->clear();
        comboBoxColorMap->insertItems(0, QStringList()
         << QApplication::translate("PanelVolume", "Grayscale", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "Lookup table", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "Heat", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "Jet", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "GE Color", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "NIH", 0, QApplication::UnicodeUTF8)
        );
        labelLookUpTable->setText(QApplication::translate("PanelVolume", "Lookup table", 0, QApplication::UnicodeUTF8));
        labelDirectionCode->setText(QApplication::translate("PanelVolume", "Color code", 0, QApplication::UnicodeUTF8));
        comboBoxDirectionCode->clear();
        comboBoxDirectionCode->insertItems(0, QStringList()
         << QApplication::translate("PanelVolume", "RAS -> RGB", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "RAS -> RBG", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "RAS -> GRB", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "RAS -> GBR", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "RAS -> BRG", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelVolume", "RAS -> BGR", 0, QApplication::UnicodeUTF8)
        );
        labelWindow->setText(QApplication::translate("PanelVolume", "Window", 0, QApplication::UnicodeUTF8));
        labelLevel->setText(QApplication::translate("PanelVolume", "Level", 0, QApplication::UnicodeUTF8));
        labelMin->setText(QApplication::translate("PanelVolume", "Min", 0, QApplication::UnicodeUTF8));
        labelMid->setText(QApplication::translate("PanelVolume", "Mid", 0, QApplication::UnicodeUTF8));
        labelOffset->setText(QApplication::translate("PanelVolume", "Offset", 0, QApplication::UnicodeUTF8));
        checkBoxTruncate->setText(QApplication::translate("PanelVolume", "Truncate", 0, QApplication::UnicodeUTF8));
        checkBoxClearBackground->setText(QApplication::translate("PanelVolume", "Clear background", 0, QApplication::UnicodeUTF8));
        checkBoxInvert->setText(QApplication::translate("PanelVolume", "Invert", 0, QApplication::UnicodeUTF8));
        checkBoxClearHigher->setText(QApplication::translate("PanelVolume", "Clear higher values", 0, QApplication::UnicodeUTF8));
        labelBrushValue->setText(QApplication::translate("PanelVolume", "Brush value", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        lineEditBrushValue->setToolTip(QApplication::translate("PanelVolume", "Enter partial name to search", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        colorLabelBrushValue->setText(QString());
        labelMax->setText(QApplication::translate("PanelVolume", "Max", 0, QApplication::UnicodeUTF8));
        checkBoxRememberFrame->setText(QApplication::translate("PanelVolume", "Remember each frame's", 0, QApplication::UnicodeUTF8));
        labelRememberFrame->setText(QApplication::translate("PanelVolume", "threshold settings", 0, QApplication::UnicodeUTF8));
        checkBoxShowExistingLabels->setText(QApplication::translate("PanelVolume", "Show existing labels only", 0, QApplication::UnicodeUTF8));
        labelTrackVolumeThreshold->setText(QApplication::translate("PanelVolume", "Threshold", 0, QApplication::UnicodeUTF8));
        checkBoxProjectionMap->setText(QApplication::translate("PanelVolume", "Show intensity projection map", 0, QApplication::UnicodeUTF8));
        checkBoxShowOutline->setText(QApplication::translate("PanelVolume", "Show label outline only (Alt+L)", 0, QApplication::UnicodeUTF8));
        checkBoxShowOutline->setShortcut(QApplication::translate("PanelVolume", "Alt+L", 0, QApplication::UnicodeUTF8));
        checkBoxShowContour->setText(QApplication::translate("PanelVolume", "Show as isosurface in 3D view", 0, QApplication::UnicodeUTF8));
        labelContourThresholdLow->setText(QApplication::translate("PanelVolume", "Low threshold", 0, QApplication::UnicodeUTF8));
        checkBoxContourExtractAll->setText(QApplication::translate("PanelVolume", "Extract all regions", 0, QApplication::UnicodeUTF8));
        labelContourColor->setText(QApplication::translate("PanelVolume", "Surface color", 0, QApplication::UnicodeUTF8));
        checkBoxUseColorMap->setText(QApplication::translate("PanelVolume", "Use color map", 0, QApplication::UnicodeUTF8));
        labelContourThresholdHigh->setText(QApplication::translate("PanelVolume", "High threshold", 0, QApplication::UnicodeUTF8));
        labelSmoothIteration->setText(QApplication::translate("PanelVolume", "Smooth iterations", 0, QApplication::UnicodeUTF8));
        pushButtonContourSave->setText(QApplication::translate("PanelVolume", "Save", 0, QApplication::UnicodeUTF8));
        checkBoxShowInfo->setText(QApplication::translate("PanelVolume", "Show in Info panel", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class PanelVolume: public Ui_PanelVolume {};
} // namespace Ui

QT_END_NAMESPACE

