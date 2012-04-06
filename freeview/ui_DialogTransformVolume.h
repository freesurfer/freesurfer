/********************************************************************************
** Form generated from reading UI file 'DialogTransformVolume.ui'
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
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QScrollBar>
#include <QtGui/QSpacerItem>
#include <QtGui/QTabWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "qtcolorpicker.h"

QT_BEGIN_NAMESPACE

class Ui_DialogTransformVolume
{
public:
    QVBoxLayout *verticalLayout_4;
    QTabWidget *tabWidget;
    QWidget *tabRotate;
    QVBoxLayout *verticalLayout_3;
    QHBoxLayout *horizontalLayout_3;
    QRadioButton *radioButtonRotateManual;
    QRadioButton *radioButtonRotateLandmarks;
    QSpacerItem *horizontalSpacer_4;
    QGroupBox *groupBox;
    QGridLayout *gridLayout;
    QCheckBox *checkBoxRotateX;
    QComboBox *comboBoxRotateX;
    QLineEdit *lineEditRotateX;
    QLabel *label_8;
    QCheckBox *checkBoxRotateY;
    QComboBox *comboBoxRotateY;
    QLineEdit *lineEditRotateY;
    QLabel *label_9;
    QCheckBox *checkBoxRotateZ;
    QComboBox *comboBoxRotateZ;
    QLineEdit *lineEditRotateZ;
    QLabel *label_10;
    QGroupBox *groupBox_2;
    QVBoxLayout *verticalLayout;
    QRadioButton *radioButtonAroundCenter;
    QRadioButton *radioButtonAroundCursor;
    QGroupBox *groupBoxLandmarks;
    QVBoxLayout *verticalLayout_7;
    QLabel *label_14;
    QGridLayout *gridLayout_4;
    QLabel *label_11;
    QPushButton *pushButtonLandmarkPick1;
    QSpacerItem *horizontalSpacer_5;
    QLabel *label_12;
    QtColorPicker *colorPickerLandmark2;
    QPushButton *pushButtonLandmarkPick2;
    QSpacerItem *horizontalSpacer_6;
    QLabel *label_13;
    QPushButton *pushButtonLandmarkPick3;
    QSpacerItem *horizontalSpacer_7;
    QLabel *label_18;
    QtColorPicker *colorPickerLandmark4;
    QPushButton *pushButtonLandmarkPick4;
    QtColorPicker *colorPickerLandmark3;
    QtColorPicker *colorPickerLandmark1;
    QSpacerItem *verticalSpacer_4;
    QLabel *label_15;
    QGridLayout *gridLayout_6;
    QComboBox *comboBoxAxis11;
    QLabel *label_16;
    QComboBox *comboBoxAxis12;
    QLabel *label_17;
    QComboBox *comboBoxAxisTarget1;
    QSpacerItem *horizontalSpacer_11;
    QComboBox *comboBoxAxis21;
    QLabel *label_25;
    QComboBox *comboBoxAxis22;
    QLabel *label_26;
    QComboBox *comboBoxAxisTarget2;
    QSpacerItem *horizontalSpacer_12;
    QGroupBox *groupBox_3;
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout_2;
    QRadioButton *radioButtonNearestNeighbor;
    QRadioButton *radioButtonTrilinear;
    QRadioButton *radioButtonCubic;
    QSpacerItem *horizontalSpacer_2;
    QLabel *label;
    QHBoxLayout *horizontalLayout_4;
    QPushButton *pushButtonApply;
    QSpacerItem *horizontalSpacer_3;
    QSpacerItem *verticalSpacer;
    QWidget *tab_2;
    QVBoxLayout *verticalLayout_5;
    QGridLayout *gridLayout_2;
    QLabel *label_2;
    QScrollBar *scrollBarTranslateX;
    QLineEdit *lineEditTranslateX;
    QLabel *label_3;
    QScrollBar *scrollBarTranslateY;
    QLineEdit *lineEditTranslateY;
    QLabel *label_4;
    QScrollBar *scrollBarTranslateZ;
    QLineEdit *lineEditTranslateZ;
    QSpacerItem *verticalSpacer_2;
    QWidget *tab;
    QVBoxLayout *verticalLayout_6;
    QGridLayout *gridLayout_3;
    QLabel *label_5;
    QScrollBar *scrollBarScaleX;
    QLineEdit *lineEditScaleX;
    QLabel *label_6;
    QScrollBar *scrollBarScaleY;
    QLineEdit *lineEditScaleY;
    QLabel *label_7;
    QScrollBar *scrollBarScaleZ;
    QLineEdit *lineEditScaleZ;
    QSpacerItem *verticalSpacer_3;
    QHBoxLayout *horizontalLayout;
    QPushButton *pushButtonRestore;
    QSpacerItem *horizontalSpacer;
    QPushButton *pushButtonSaveReg;
    QPushButton *pushButtonSaveVolumeAs;

    void setupUi(QDialog *DialogTransformVolume)
    {
        if (DialogTransformVolume->objectName().isEmpty())
            DialogTransformVolume->setObjectName(QString::fromUtf8("DialogTransformVolume"));
        DialogTransformVolume->resize(558, 836);
        verticalLayout_4 = new QVBoxLayout(DialogTransformVolume);
        verticalLayout_4->setSpacing(8);
        verticalLayout_4->setContentsMargins(8, 8, 8, 8);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        verticalLayout_4->setSizeConstraint(QLayout::SetFixedSize);
        tabWidget = new QTabWidget(DialogTransformVolume);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabRotate = new QWidget();
        tabRotate->setObjectName(QString::fromUtf8("tabRotate"));
        verticalLayout_3 = new QVBoxLayout(tabRotate);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        horizontalLayout_3->setContentsMargins(-1, 0, -1, -1);
        radioButtonRotateManual = new QRadioButton(tabRotate);
        radioButtonRotateManual->setObjectName(QString::fromUtf8("radioButtonRotateManual"));
        radioButtonRotateManual->setCheckable(true);
        radioButtonRotateManual->setChecked(true);

        horizontalLayout_3->addWidget(radioButtonRotateManual);

        radioButtonRotateLandmarks = new QRadioButton(tabRotate);
        radioButtonRotateLandmarks->setObjectName(QString::fromUtf8("radioButtonRotateLandmarks"));
        radioButtonRotateLandmarks->setAutoExclusive(true);

        horizontalLayout_3->addWidget(radioButtonRotateLandmarks);

        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_4);


        verticalLayout_3->addLayout(horizontalLayout_3);

        groupBox = new QGroupBox(tabRotate);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        gridLayout = new QGridLayout(groupBox);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setHorizontalSpacing(6);
        checkBoxRotateX = new QCheckBox(groupBox);
        checkBoxRotateX->setObjectName(QString::fromUtf8("checkBoxRotateX"));
        checkBoxRotateX->setChecked(true);

        gridLayout->addWidget(checkBoxRotateX, 0, 0, 1, 1);

        comboBoxRotateX = new QComboBox(groupBox);
        comboBoxRotateX->setObjectName(QString::fromUtf8("comboBoxRotateX"));
        comboBoxRotateX->setMaximumSize(QSize(16777215, 16777215));

        gridLayout->addWidget(comboBoxRotateX, 0, 1, 1, 1);

        lineEditRotateX = new QLineEdit(groupBox);
        lineEditRotateX->setObjectName(QString::fromUtf8("lineEditRotateX"));
        lineEditRotateX->setMaximumSize(QSize(65, 16777215));

        gridLayout->addWidget(lineEditRotateX, 0, 2, 1, 1);

        label_8 = new QLabel(groupBox);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        gridLayout->addWidget(label_8, 0, 3, 1, 1);

        checkBoxRotateY = new QCheckBox(groupBox);
        checkBoxRotateY->setObjectName(QString::fromUtf8("checkBoxRotateY"));

        gridLayout->addWidget(checkBoxRotateY, 1, 0, 1, 1);

        comboBoxRotateY = new QComboBox(groupBox);
        comboBoxRotateY->setObjectName(QString::fromUtf8("comboBoxRotateY"));
        comboBoxRotateY->setEnabled(false);

        gridLayout->addWidget(comboBoxRotateY, 1, 1, 1, 1);

        lineEditRotateY = new QLineEdit(groupBox);
        lineEditRotateY->setObjectName(QString::fromUtf8("lineEditRotateY"));
        lineEditRotateY->setEnabled(false);
        lineEditRotateY->setMaximumSize(QSize(65, 16777215));

        gridLayout->addWidget(lineEditRotateY, 1, 2, 1, 1);

        label_9 = new QLabel(groupBox);
        label_9->setObjectName(QString::fromUtf8("label_9"));
        label_9->setEnabled(false);

        gridLayout->addWidget(label_9, 1, 3, 1, 1);

        checkBoxRotateZ = new QCheckBox(groupBox);
        checkBoxRotateZ->setObjectName(QString::fromUtf8("checkBoxRotateZ"));

        gridLayout->addWidget(checkBoxRotateZ, 2, 0, 1, 1);

        comboBoxRotateZ = new QComboBox(groupBox);
        comboBoxRotateZ->setObjectName(QString::fromUtf8("comboBoxRotateZ"));
        comboBoxRotateZ->setEnabled(false);

        gridLayout->addWidget(comboBoxRotateZ, 2, 1, 1, 1);

        lineEditRotateZ = new QLineEdit(groupBox);
        lineEditRotateZ->setObjectName(QString::fromUtf8("lineEditRotateZ"));
        lineEditRotateZ->setEnabled(false);
        lineEditRotateZ->setMaximumSize(QSize(65, 16777215));

        gridLayout->addWidget(lineEditRotateZ, 2, 2, 1, 1);

        label_10 = new QLabel(groupBox);
        label_10->setObjectName(QString::fromUtf8("label_10"));
        label_10->setEnabled(false);

        gridLayout->addWidget(label_10, 2, 3, 1, 1);

        gridLayout->setColumnStretch(1, 1);

        verticalLayout_3->addWidget(groupBox);

        groupBox_2 = new QGroupBox(tabRotate);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        verticalLayout = new QVBoxLayout(groupBox_2);
        verticalLayout->setSpacing(9);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        radioButtonAroundCenter = new QRadioButton(groupBox_2);
        radioButtonAroundCenter->setObjectName(QString::fromUtf8("radioButtonAroundCenter"));
        radioButtonAroundCenter->setChecked(true);

        verticalLayout->addWidget(radioButtonAroundCenter);

        radioButtonAroundCursor = new QRadioButton(groupBox_2);
        radioButtonAroundCursor->setObjectName(QString::fromUtf8("radioButtonAroundCursor"));

        verticalLayout->addWidget(radioButtonAroundCursor);


        verticalLayout_3->addWidget(groupBox_2);

        groupBoxLandmarks = new QGroupBox(tabRotate);
        groupBoxLandmarks->setObjectName(QString::fromUtf8("groupBoxLandmarks"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(groupBoxLandmarks->sizePolicy().hasHeightForWidth());
        groupBoxLandmarks->setSizePolicy(sizePolicy);
        verticalLayout_7 = new QVBoxLayout(groupBoxLandmarks);
        verticalLayout_7->setObjectName(QString::fromUtf8("verticalLayout_7"));
        label_14 = new QLabel(groupBoxLandmarks);
        label_14->setObjectName(QString::fromUtf8("label_14"));
        label_14->setWordWrap(true);

        verticalLayout_7->addWidget(label_14);

        gridLayout_4 = new QGridLayout();
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        label_11 = new QLabel(groupBoxLandmarks);
        label_11->setObjectName(QString::fromUtf8("label_11"));

        gridLayout_4->addWidget(label_11, 0, 0, 1, 1);

        pushButtonLandmarkPick1 = new QPushButton(groupBoxLandmarks);
        pushButtonLandmarkPick1->setObjectName(QString::fromUtf8("pushButtonLandmarkPick1"));
        pushButtonLandmarkPick1->setCheckable(true);
        pushButtonLandmarkPick1->setAutoExclusive(false);
        pushButtonLandmarkPick1->setAutoDefault(false);

        gridLayout_4->addWidget(pushButtonLandmarkPick1, 0, 2, 1, 1);

        horizontalSpacer_5 = new QSpacerItem(20, 5, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_5, 0, 3, 1, 1);

        label_12 = new QLabel(groupBoxLandmarks);
        label_12->setObjectName(QString::fromUtf8("label_12"));

        gridLayout_4->addWidget(label_12, 0, 4, 1, 1);

        colorPickerLandmark2 = new QtColorPicker(groupBoxLandmarks);
        colorPickerLandmark2->setObjectName(QString::fromUtf8("colorPickerLandmark2"));

        gridLayout_4->addWidget(colorPickerLandmark2, 0, 5, 1, 1);

        pushButtonLandmarkPick2 = new QPushButton(groupBoxLandmarks);
        pushButtonLandmarkPick2->setObjectName(QString::fromUtf8("pushButtonLandmarkPick2"));
        pushButtonLandmarkPick2->setCheckable(true);
        pushButtonLandmarkPick2->setAutoExclusive(false);
        pushButtonLandmarkPick2->setAutoDefault(false);

        gridLayout_4->addWidget(pushButtonLandmarkPick2, 0, 6, 1, 1);

        horizontalSpacer_6 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_6, 0, 7, 1, 1);

        label_13 = new QLabel(groupBoxLandmarks);
        label_13->setObjectName(QString::fromUtf8("label_13"));

        gridLayout_4->addWidget(label_13, 1, 0, 1, 1);

        pushButtonLandmarkPick3 = new QPushButton(groupBoxLandmarks);
        pushButtonLandmarkPick3->setObjectName(QString::fromUtf8("pushButtonLandmarkPick3"));
        pushButtonLandmarkPick3->setCheckable(true);
        pushButtonLandmarkPick3->setAutoExclusive(false);
        pushButtonLandmarkPick3->setAutoDefault(false);

        gridLayout_4->addWidget(pushButtonLandmarkPick3, 1, 2, 1, 1);

        horizontalSpacer_7 = new QSpacerItem(20, 5, QSizePolicy::Fixed, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_7, 1, 3, 1, 1);

        label_18 = new QLabel(groupBoxLandmarks);
        label_18->setObjectName(QString::fromUtf8("label_18"));

        gridLayout_4->addWidget(label_18, 1, 4, 1, 1);

        colorPickerLandmark4 = new QtColorPicker(groupBoxLandmarks);
        colorPickerLandmark4->setObjectName(QString::fromUtf8("colorPickerLandmark4"));

        gridLayout_4->addWidget(colorPickerLandmark4, 1, 5, 1, 1);

        pushButtonLandmarkPick4 = new QPushButton(groupBoxLandmarks);
        pushButtonLandmarkPick4->setObjectName(QString::fromUtf8("pushButtonLandmarkPick4"));
        pushButtonLandmarkPick4->setCheckable(true);
        pushButtonLandmarkPick4->setAutoExclusive(false);
        pushButtonLandmarkPick4->setAutoDefault(false);

        gridLayout_4->addWidget(pushButtonLandmarkPick4, 1, 6, 1, 1);

        colorPickerLandmark3 = new QtColorPicker(groupBoxLandmarks);
        colorPickerLandmark3->setObjectName(QString::fromUtf8("colorPickerLandmark3"));

        gridLayout_4->addWidget(colorPickerLandmark3, 1, 1, 1, 1);

        colorPickerLandmark1 = new QtColorPicker(groupBoxLandmarks);
        colorPickerLandmark1->setObjectName(QString::fromUtf8("colorPickerLandmark1"));

        gridLayout_4->addWidget(colorPickerLandmark1, 0, 1, 1, 1);


        verticalLayout_7->addLayout(gridLayout_4);

        verticalSpacer_4 = new QSpacerItem(20, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout_7->addItem(verticalSpacer_4);

        label_15 = new QLabel(groupBoxLandmarks);
        label_15->setObjectName(QString::fromUtf8("label_15"));
        label_15->setWordWrap(true);

        verticalLayout_7->addWidget(label_15);

        gridLayout_6 = new QGridLayout();
        gridLayout_6->setObjectName(QString::fromUtf8("gridLayout_6"));
        gridLayout_6->setHorizontalSpacing(10);
        comboBoxAxis11 = new QComboBox(groupBoxLandmarks);
        comboBoxAxis11->setObjectName(QString::fromUtf8("comboBoxAxis11"));
        comboBoxAxis11->setMinimumSize(QSize(0, 0));
        comboBoxAxis11->setMaximumSize(QSize(80, 16777215));

        gridLayout_6->addWidget(comboBoxAxis11, 0, 0, 1, 1);

        label_16 = new QLabel(groupBoxLandmarks);
        label_16->setObjectName(QString::fromUtf8("label_16"));

        gridLayout_6->addWidget(label_16, 0, 1, 1, 1);

        comboBoxAxis12 = new QComboBox(groupBoxLandmarks);
        comboBoxAxis12->setObjectName(QString::fromUtf8("comboBoxAxis12"));
        comboBoxAxis12->setMinimumSize(QSize(55, 0));
        comboBoxAxis12->setMaximumSize(QSize(80, 16777215));

        gridLayout_6->addWidget(comboBoxAxis12, 0, 2, 1, 1);

        label_17 = new QLabel(groupBoxLandmarks);
        label_17->setObjectName(QString::fromUtf8("label_17"));

        gridLayout_6->addWidget(label_17, 0, 3, 1, 1);

        comboBoxAxisTarget1 = new QComboBox(groupBoxLandmarks);
        comboBoxAxisTarget1->setObjectName(QString::fromUtf8("comboBoxAxisTarget1"));
        comboBoxAxisTarget1->setMinimumSize(QSize(0, 0));
        comboBoxAxisTarget1->setMaximumSize(QSize(60, 16777215));

        gridLayout_6->addWidget(comboBoxAxisTarget1, 0, 4, 1, 1);

        horizontalSpacer_11 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_6->addItem(horizontalSpacer_11, 0, 5, 1, 1);

        comboBoxAxis21 = new QComboBox(groupBoxLandmarks);
        comboBoxAxis21->setObjectName(QString::fromUtf8("comboBoxAxis21"));
        comboBoxAxis21->setMaximumSize(QSize(80, 16777215));

        gridLayout_6->addWidget(comboBoxAxis21, 1, 0, 1, 1);

        label_25 = new QLabel(groupBoxLandmarks);
        label_25->setObjectName(QString::fromUtf8("label_25"));

        gridLayout_6->addWidget(label_25, 1, 1, 1, 1);

        comboBoxAxis22 = new QComboBox(groupBoxLandmarks);
        comboBoxAxis22->setObjectName(QString::fromUtf8("comboBoxAxis22"));
        comboBoxAxis22->setMaximumSize(QSize(80, 16777215));

        gridLayout_6->addWidget(comboBoxAxis22, 1, 2, 1, 1);

        label_26 = new QLabel(groupBoxLandmarks);
        label_26->setObjectName(QString::fromUtf8("label_26"));

        gridLayout_6->addWidget(label_26, 1, 3, 1, 1);

        comboBoxAxisTarget2 = new QComboBox(groupBoxLandmarks);
        comboBoxAxisTarget2->setObjectName(QString::fromUtf8("comboBoxAxisTarget2"));
        comboBoxAxisTarget2->setMaximumSize(QSize(60, 16777215));

        gridLayout_6->addWidget(comboBoxAxisTarget2, 1, 4, 1, 1);

        horizontalSpacer_12 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_6->addItem(horizontalSpacer_12, 1, 5, 1, 1);


        verticalLayout_7->addLayout(gridLayout_6);


        verticalLayout_3->addWidget(groupBoxLandmarks);

        groupBox_3 = new QGroupBox(tabRotate);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        verticalLayout_2 = new QVBoxLayout(groupBox_3);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        radioButtonNearestNeighbor = new QRadioButton(groupBox_3);
        radioButtonNearestNeighbor->setObjectName(QString::fromUtf8("radioButtonNearestNeighbor"));

        horizontalLayout_2->addWidget(radioButtonNearestNeighbor);

        radioButtonTrilinear = new QRadioButton(groupBox_3);
        radioButtonTrilinear->setObjectName(QString::fromUtf8("radioButtonTrilinear"));
        radioButtonTrilinear->setChecked(true);

        horizontalLayout_2->addWidget(radioButtonTrilinear);

        radioButtonCubic = new QRadioButton(groupBox_3);
        radioButtonCubic->setObjectName(QString::fromUtf8("radioButtonCubic"));
        radioButtonCubic->setChecked(false);

        horizontalLayout_2->addWidget(radioButtonCubic);

        horizontalSpacer_2 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_2);


        verticalLayout_2->addLayout(horizontalLayout_2);

        label = new QLabel(groupBox_3);
        label->setObjectName(QString::fromUtf8("label"));
        label->setScaledContents(false);
        label->setWordWrap(true);

        verticalLayout_2->addWidget(label);


        verticalLayout_3->addWidget(groupBox_3);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(-1, 0, -1, -1);
        pushButtonApply = new QPushButton(tabRotate);
        pushButtonApply->setObjectName(QString::fromUtf8("pushButtonApply"));
        pushButtonApply->setDefault(true);

        horizontalLayout_4->addWidget(pushButtonApply);

        horizontalSpacer_3 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_3);


        verticalLayout_3->addLayout(horizontalLayout_4);

        verticalSpacer = new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_3->addItem(verticalSpacer);

        tabWidget->addTab(tabRotate, QString());
        tab_2 = new QWidget();
        tab_2->setObjectName(QString::fromUtf8("tab_2"));
        verticalLayout_5 = new QVBoxLayout(tab_2);
        verticalLayout_5->setObjectName(QString::fromUtf8("verticalLayout_5"));
        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        gridLayout_2->setHorizontalSpacing(8);
        label_2 = new QLabel(tab_2);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout_2->addWidget(label_2, 0, 0, 1, 1);

        scrollBarTranslateX = new QScrollBar(tab_2);
        scrollBarTranslateX->setObjectName(QString::fromUtf8("scrollBarTranslateX"));
        scrollBarTranslateX->setMaximum(100);
        scrollBarTranslateX->setValue(50);
        scrollBarTranslateX->setOrientation(Qt::Horizontal);

        gridLayout_2->addWidget(scrollBarTranslateX, 0, 1, 1, 1);

        lineEditTranslateX = new QLineEdit(tab_2);
        lineEditTranslateX->setObjectName(QString::fromUtf8("lineEditTranslateX"));
        lineEditTranslateX->setMaximumSize(QSize(65, 16777215));

        gridLayout_2->addWidget(lineEditTranslateX, 0, 2, 1, 1);

        label_3 = new QLabel(tab_2);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout_2->addWidget(label_3, 1, 0, 1, 1);

        scrollBarTranslateY = new QScrollBar(tab_2);
        scrollBarTranslateY->setObjectName(QString::fromUtf8("scrollBarTranslateY"));
        scrollBarTranslateY->setMaximum(100);
        scrollBarTranslateY->setValue(50);
        scrollBarTranslateY->setOrientation(Qt::Horizontal);

        gridLayout_2->addWidget(scrollBarTranslateY, 1, 1, 1, 1);

        lineEditTranslateY = new QLineEdit(tab_2);
        lineEditTranslateY->setObjectName(QString::fromUtf8("lineEditTranslateY"));
        lineEditTranslateY->setMaximumSize(QSize(65, 16777215));

        gridLayout_2->addWidget(lineEditTranslateY, 1, 2, 1, 1);

        label_4 = new QLabel(tab_2);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout_2->addWidget(label_4, 2, 0, 1, 1);

        scrollBarTranslateZ = new QScrollBar(tab_2);
        scrollBarTranslateZ->setObjectName(QString::fromUtf8("scrollBarTranslateZ"));
        scrollBarTranslateZ->setMaximum(100);
        scrollBarTranslateZ->setValue(50);
        scrollBarTranslateZ->setOrientation(Qt::Horizontal);

        gridLayout_2->addWidget(scrollBarTranslateZ, 2, 1, 1, 1);

        lineEditTranslateZ = new QLineEdit(tab_2);
        lineEditTranslateZ->setObjectName(QString::fromUtf8("lineEditTranslateZ"));
        lineEditTranslateZ->setMaximumSize(QSize(65, 16777215));

        gridLayout_2->addWidget(lineEditTranslateZ, 2, 2, 1, 1);

        gridLayout_2->setColumnStretch(1, 1);

        verticalLayout_5->addLayout(gridLayout_2);

        verticalSpacer_2 = new QSpacerItem(13, 283, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_5->addItem(verticalSpacer_2);

        tabWidget->addTab(tab_2, QString());
        tab = new QWidget();
        tab->setObjectName(QString::fromUtf8("tab"));
        verticalLayout_6 = new QVBoxLayout(tab);
        verticalLayout_6->setObjectName(QString::fromUtf8("verticalLayout_6"));
        gridLayout_3 = new QGridLayout();
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        gridLayout_3->setHorizontalSpacing(8);
        gridLayout_3->setVerticalSpacing(6);
        label_5 = new QLabel(tab);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        gridLayout_3->addWidget(label_5, 0, 0, 1, 1);

        scrollBarScaleX = new QScrollBar(tab);
        scrollBarScaleX->setObjectName(QString::fromUtf8("scrollBarScaleX"));
        scrollBarScaleX->setMaximum(100);
        scrollBarScaleX->setValue(50);
        scrollBarScaleX->setOrientation(Qt::Horizontal);

        gridLayout_3->addWidget(scrollBarScaleX, 0, 1, 1, 1);

        lineEditScaleX = new QLineEdit(tab);
        lineEditScaleX->setObjectName(QString::fromUtf8("lineEditScaleX"));
        lineEditScaleX->setMaximumSize(QSize(65, 16777215));

        gridLayout_3->addWidget(lineEditScaleX, 0, 2, 1, 1);

        label_6 = new QLabel(tab);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        gridLayout_3->addWidget(label_6, 1, 0, 1, 1);

        scrollBarScaleY = new QScrollBar(tab);
        scrollBarScaleY->setObjectName(QString::fromUtf8("scrollBarScaleY"));
        scrollBarScaleY->setMaximum(100);
        scrollBarScaleY->setValue(50);
        scrollBarScaleY->setOrientation(Qt::Horizontal);

        gridLayout_3->addWidget(scrollBarScaleY, 1, 1, 1, 1);

        lineEditScaleY = new QLineEdit(tab);
        lineEditScaleY->setObjectName(QString::fromUtf8("lineEditScaleY"));
        lineEditScaleY->setMaximumSize(QSize(65, 16777215));

        gridLayout_3->addWidget(lineEditScaleY, 1, 2, 1, 1);

        label_7 = new QLabel(tab);
        label_7->setObjectName(QString::fromUtf8("label_7"));

        gridLayout_3->addWidget(label_7, 2, 0, 1, 1);

        scrollBarScaleZ = new QScrollBar(tab);
        scrollBarScaleZ->setObjectName(QString::fromUtf8("scrollBarScaleZ"));
        scrollBarScaleZ->setMaximum(100);
        scrollBarScaleZ->setValue(50);
        scrollBarScaleZ->setOrientation(Qt::Horizontal);

        gridLayout_3->addWidget(scrollBarScaleZ, 2, 1, 1, 1);

        lineEditScaleZ = new QLineEdit(tab);
        lineEditScaleZ->setObjectName(QString::fromUtf8("lineEditScaleZ"));
        lineEditScaleZ->setMaximumSize(QSize(65, 16777215));

        gridLayout_3->addWidget(lineEditScaleZ, 2, 2, 1, 1);

        gridLayout_3->setColumnStretch(1, 1);

        verticalLayout_6->addLayout(gridLayout_3);

        verticalSpacer_3 = new QSpacerItem(13, 283, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_6->addItem(verticalSpacer_3);

        tabWidget->addTab(tab, QString());

        verticalLayout_4->addWidget(tabWidget);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        pushButtonRestore = new QPushButton(DialogTransformVolume);
        pushButtonRestore->setObjectName(QString::fromUtf8("pushButtonRestore"));
        pushButtonRestore->setAutoDefault(false);

        horizontalLayout->addWidget(pushButtonRestore);

        horizontalSpacer = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);

        pushButtonSaveReg = new QPushButton(DialogTransformVolume);
        pushButtonSaveReg->setObjectName(QString::fromUtf8("pushButtonSaveReg"));
        pushButtonSaveReg->setAutoDefault(false);

        horizontalLayout->addWidget(pushButtonSaveReg);

        pushButtonSaveVolumeAs = new QPushButton(DialogTransformVolume);
        pushButtonSaveVolumeAs->setObjectName(QString::fromUtf8("pushButtonSaveVolumeAs"));
        pushButtonSaveVolumeAs->setAutoDefault(false);
        pushButtonSaveVolumeAs->setDefault(false);

        horizontalLayout->addWidget(pushButtonSaveVolumeAs);


        verticalLayout_4->addLayout(horizontalLayout);


        retranslateUi(DialogTransformVolume);
        QObject::connect(checkBoxRotateY, SIGNAL(toggled(bool)), comboBoxRotateY, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxRotateY, SIGNAL(toggled(bool)), lineEditRotateY, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxRotateY, SIGNAL(toggled(bool)), label_9, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxRotateX, SIGNAL(toggled(bool)), comboBoxRotateX, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxRotateX, SIGNAL(toggled(bool)), lineEditRotateX, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxRotateX, SIGNAL(toggled(bool)), label_8, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxRotateZ, SIGNAL(toggled(bool)), comboBoxRotateZ, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxRotateZ, SIGNAL(toggled(bool)), lineEditRotateZ, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxRotateZ, SIGNAL(toggled(bool)), label_10, SLOT(setEnabled(bool)));
        QObject::connect(pushButtonApply, SIGNAL(clicked()), DialogTransformVolume, SLOT(OnApply()));
        QObject::connect(pushButtonRestore, SIGNAL(clicked()), DialogTransformVolume, SLOT(OnRestore()));
        QObject::connect(pushButtonSaveReg, SIGNAL(clicked()), DialogTransformVolume, SLOT(OnSaveReg()));
        QObject::connect(scrollBarTranslateX, SIGNAL(valueChanged(int)), DialogTransformVolume, SLOT(OnScrollBarTranslateX(int)));
        QObject::connect(scrollBarTranslateY, SIGNAL(valueChanged(int)), DialogTransformVolume, SLOT(OnScrollBarTranslateY(int)));
        QObject::connect(scrollBarTranslateZ, SIGNAL(valueChanged(int)), DialogTransformVolume, SLOT(OnScrollBarTranslateZ(int)));
        QObject::connect(lineEditTranslateX, SIGNAL(textEdited(QString)), DialogTransformVolume, SLOT(OnLineEditTranslateX(QString)));
        QObject::connect(lineEditTranslateY, SIGNAL(textEdited(QString)), DialogTransformVolume, SLOT(OnLineEditTranslateY(QString)));
        QObject::connect(lineEditTranslateZ, SIGNAL(textEdited(QString)), DialogTransformVolume, SLOT(OnLineEditTranslateZ(QString)));
        QObject::connect(scrollBarScaleX, SIGNAL(valueChanged(int)), DialogTransformVolume, SLOT(OnScrollBarScaleX(int)));
        QObject::connect(scrollBarScaleY, SIGNAL(valueChanged(int)), DialogTransformVolume, SLOT(OnScrollBarScaleY(int)));
        QObject::connect(scrollBarScaleZ, SIGNAL(valueChanged(int)), DialogTransformVolume, SLOT(OnScrollBarScaleZ(int)));
        QObject::connect(lineEditScaleX, SIGNAL(textEdited(QString)), DialogTransformVolume, SLOT(OnLineEditScaleX(QString)));
        QObject::connect(lineEditScaleY, SIGNAL(textEdited(QString)), DialogTransformVolume, SLOT(OnLineEditScaleY(QString)));
        QObject::connect(lineEditScaleZ, SIGNAL(textEdited(QString)), DialogTransformVolume, SLOT(OnLineEditScaleZ(QString)));
        QObject::connect(radioButtonRotateManual, SIGNAL(toggled(bool)), groupBox, SLOT(setShown(bool)));
        QObject::connect(radioButtonRotateManual, SIGNAL(toggled(bool)), groupBox_2, SLOT(setShown(bool)));
        QObject::connect(radioButtonRotateManual, SIGNAL(toggled(bool)), groupBoxLandmarks, SLOT(setHidden(bool)));
        QObject::connect(pushButtonLandmarkPick2, SIGNAL(toggled(bool)), DialogTransformVolume, SLOT(OnButtonLandmarkPick()));
        QObject::connect(pushButtonLandmarkPick3, SIGNAL(toggled(bool)), DialogTransformVolume, SLOT(OnButtonLandmarkPick()));
        QObject::connect(pushButtonLandmarkPick1, SIGNAL(toggled(bool)), DialogTransformVolume, SLOT(OnButtonLandmarkPick()));
        QObject::connect(colorPickerLandmark1, SIGNAL(colorChanged(QColor)), DialogTransformVolume, SLOT(UpdateLandmarkColors()));
        QObject::connect(colorPickerLandmark2, SIGNAL(colorChanged(QColor)), DialogTransformVolume, SLOT(UpdateLandmarkColors()));
        QObject::connect(colorPickerLandmark3, SIGNAL(colorChanged(QColor)), DialogTransformVolume, SLOT(UpdateLandmarkColors()));
        QObject::connect(radioButtonRotateLandmarks, SIGNAL(toggled(bool)), DialogTransformVolume, SLOT(OnRadioButtonLandmark(bool)));
        QObject::connect(pushButtonLandmarkPick4, SIGNAL(toggled(bool)), DialogTransformVolume, SLOT(OnButtonLandmarkPick()));
        QObject::connect(colorPickerLandmark4, SIGNAL(colorChanged(QColor)), DialogTransformVolume, SLOT(UpdateLandmarkColors()));
        QObject::connect(radioButtonNearestNeighbor, SIGNAL(toggled(bool)), DialogTransformVolume, SLOT(OnSampleMethodChanged()));
        QObject::connect(radioButtonTrilinear, SIGNAL(toggled(bool)), DialogTransformVolume, SLOT(OnSampleMethodChanged()));
        QObject::connect(radioButtonCubic, SIGNAL(toggled(bool)), DialogTransformVolume, SLOT(OnSampleMethodChanged()));

        tabWidget->setCurrentIndex(0);
        comboBoxRotateY->setCurrentIndex(1);
        comboBoxRotateZ->setCurrentIndex(2);
        comboBoxAxis11->setCurrentIndex(0);
        comboBoxAxis12->setCurrentIndex(1);
        comboBoxAxis21->setCurrentIndex(2);
        comboBoxAxis22->setCurrentIndex(3);
        comboBoxAxisTarget2->setCurrentIndex(1);


        QMetaObject::connectSlotsByName(DialogTransformVolume);
    } // setupUi

    void retranslateUi(QDialog *DialogTransformVolume)
    {
        DialogTransformVolume->setWindowTitle(QApplication::translate("DialogTransformVolume", "Dialog", 0, QApplication::UnicodeUTF8));
        radioButtonRotateManual->setText(QApplication::translate("DialogTransformVolume", "Manual", 0, QApplication::UnicodeUTF8));
        radioButtonRotateLandmarks->setText(QApplication::translate("DialogTransformVolume", "Use Landmarks", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("DialogTransformVolume", "Rotations", 0, QApplication::UnicodeUTF8));
        checkBoxRotateX->setText(QString());
        comboBoxRotateX->clear();
        comboBoxRotateX->insertItems(0, QStringList()
         << QApplication::translate("DialogTransformVolume", "Rotate X", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "Rotate Y", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "Rotate Z", 0, QApplication::UnicodeUTF8)
        );
        lineEditRotateX->setText(QApplication::translate("DialogTransformVolume", "0", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("DialogTransformVolume", "degrees", 0, QApplication::UnicodeUTF8));
        checkBoxRotateY->setText(QString());
        comboBoxRotateY->clear();
        comboBoxRotateY->insertItems(0, QStringList()
         << QApplication::translate("DialogTransformVolume", "Rotate X", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "Rotate Y", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "Rotate Z", 0, QApplication::UnicodeUTF8)
        );
        lineEditRotateY->setText(QApplication::translate("DialogTransformVolume", "0", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("DialogTransformVolume", "degrees", 0, QApplication::UnicodeUTF8));
        checkBoxRotateZ->setText(QString());
        comboBoxRotateZ->clear();
        comboBoxRotateZ->insertItems(0, QStringList()
         << QApplication::translate("DialogTransformVolume", "Rotate X", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "Rotate Y", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "Rotate Z", 0, QApplication::UnicodeUTF8)
        );
        lineEditRotateZ->setText(QApplication::translate("DialogTransformVolume", "0", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("DialogTransformVolume", "degrees", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("DialogTransformVolume", "Axis Point", 0, QApplication::UnicodeUTF8));
        radioButtonAroundCenter->setText(QApplication::translate("DialogTransformVolume", "Rotate around the center of the volume", 0, QApplication::UnicodeUTF8));
        radioButtonAroundCursor->setText(QApplication::translate("DialogTransformVolume", "Rotate around the cursor", 0, QApplication::UnicodeUTF8));
        groupBoxLandmarks->setTitle(QApplication::translate("DialogTransformVolume", "Landmarks", 0, QApplication::UnicodeUTF8));
        label_14->setText(QApplication::translate("DialogTransformVolume", "First click Pick button to pick landmark positions on the slice view. Need to pick at least 3 landmarks.", 0, QApplication::UnicodeUTF8));
        label_11->setText(QApplication::translate("DialogTransformVolume", "#1", 0, QApplication::UnicodeUTF8));
        pushButtonLandmarkPick1->setText(QApplication::translate("DialogTransformVolume", "Pick", 0, QApplication::UnicodeUTF8));
        label_12->setText(QApplication::translate("DialogTransformVolume", "#2", 0, QApplication::UnicodeUTF8));
        pushButtonLandmarkPick2->setText(QApplication::translate("DialogTransformVolume", "Pick", 0, QApplication::UnicodeUTF8));
        label_13->setText(QApplication::translate("DialogTransformVolume", "#3", 0, QApplication::UnicodeUTF8));
        pushButtonLandmarkPick3->setText(QApplication::translate("DialogTransformVolume", "Pick", 0, QApplication::UnicodeUTF8));
        label_18->setText(QApplication::translate("DialogTransformVolume", "#4", 0, QApplication::UnicodeUTF8));
        pushButtonLandmarkPick4->setText(QApplication::translate("DialogTransformVolume", "Pick", 0, QApplication::UnicodeUTF8));
        label_15->setText(QApplication::translate("DialogTransformVolume", "Then set how to map landmarks to coordinate axis vectors", 0, QApplication::UnicodeUTF8));
        comboBoxAxis11->clear();
        comboBoxAxis11->insertItems(0, QStringList()
         << QApplication::translate("DialogTransformVolume", "#1", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "#2", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "#3", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "#4", 0, QApplication::UnicodeUTF8)
        );
        label_16->setText(QApplication::translate("DialogTransformVolume", "->", 0, QApplication::UnicodeUTF8));
        comboBoxAxis12->clear();
        comboBoxAxis12->insertItems(0, QStringList()
         << QApplication::translate("DialogTransformVolume", "#1", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "#2", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "#3", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "#4", 0, QApplication::UnicodeUTF8)
        );
        label_17->setText(QApplication::translate("DialogTransformVolume", "=", 0, QApplication::UnicodeUTF8));
        comboBoxAxisTarget1->clear();
        comboBoxAxisTarget1->insertItems(0, QStringList()
         << QApplication::translate("DialogTransformVolume", "L->R", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "P->A", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "I->S", 0, QApplication::UnicodeUTF8)
        );
        comboBoxAxis21->clear();
        comboBoxAxis21->insertItems(0, QStringList()
         << QApplication::translate("DialogTransformVolume", "#1", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "#2", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "#3", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "#4", 0, QApplication::UnicodeUTF8)
        );
        label_25->setText(QApplication::translate("DialogTransformVolume", "->", 0, QApplication::UnicodeUTF8));
        comboBoxAxis22->clear();
        comboBoxAxis22->insertItems(0, QStringList()
         << QApplication::translate("DialogTransformVolume", "#1", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "#2", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "#3", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "#4", 0, QApplication::UnicodeUTF8)
        );
        label_26->setText(QApplication::translate("DialogTransformVolume", "=", 0, QApplication::UnicodeUTF8));
        comboBoxAxisTarget2->clear();
        comboBoxAxisTarget2->insertItems(0, QStringList()
         << QApplication::translate("DialogTransformVolume", "L->R", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "P->A", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogTransformVolume", "I->S", 0, QApplication::UnicodeUTF8)
        );
        groupBox_3->setTitle(QApplication::translate("DialogTransformVolume", "Sample Method", 0, QApplication::UnicodeUTF8));
        radioButtonNearestNeighbor->setText(QApplication::translate("DialogTransformVolume", "Nearest Neighbor", 0, QApplication::UnicodeUTF8));
        radioButtonTrilinear->setText(QApplication::translate("DialogTransformVolume", "Trilinear", 0, QApplication::UnicodeUTF8));
        radioButtonCubic->setText(QApplication::translate("DialogTransformVolume", "Cubic", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogTransformVolume", "Interpolation method on label volumes will always be nearest neighbor. Therefore, there is risk to lose voxels after rotation.", 0, QApplication::UnicodeUTF8));
        pushButtonApply->setText(QApplication::translate("DialogTransformVolume", "Apply", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tabRotate), QApplication::translate("DialogTransformVolume", "Rotate", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("DialogTransformVolume", "X (L-R)", 0, QApplication::UnicodeUTF8));
        lineEditTranslateX->setText(QApplication::translate("DialogTransformVolume", "0", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("DialogTransformVolume", "Y (P-A)", 0, QApplication::UnicodeUTF8));
        lineEditTranslateY->setText(QApplication::translate("DialogTransformVolume", "0", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("DialogTransformVolume", "Z (I-S)", 0, QApplication::UnicodeUTF8));
        lineEditTranslateZ->setText(QApplication::translate("DialogTransformVolume", "0", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab_2), QApplication::translate("DialogTransformVolume", "Translate", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("DialogTransformVolume", "X (L-R)", 0, QApplication::UnicodeUTF8));
        lineEditScaleX->setText(QApplication::translate("DialogTransformVolume", "1", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("DialogTransformVolume", "Y (P-A)", 0, QApplication::UnicodeUTF8));
        lineEditScaleY->setText(QApplication::translate("DialogTransformVolume", "1", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("DialogTransformVolume", "Z (I-S)", 0, QApplication::UnicodeUTF8));
        lineEditScaleZ->setText(QApplication::translate("DialogTransformVolume", "1", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tab), QApplication::translate("DialogTransformVolume", "Scale", 0, QApplication::UnicodeUTF8));
        pushButtonRestore->setText(QApplication::translate("DialogTransformVolume", "Restore To Original", 0, QApplication::UnicodeUTF8));
        pushButtonSaveReg->setText(QApplication::translate("DialogTransformVolume", "Save Reg...", 0, QApplication::UnicodeUTF8));
        pushButtonSaveVolumeAs->setText(QApplication::translate("DialogTransformVolume", "Save Volume As...", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogTransformVolume: public Ui_DialogTransformVolume {};
} // namespace Ui

QT_END_NAMESPACE

