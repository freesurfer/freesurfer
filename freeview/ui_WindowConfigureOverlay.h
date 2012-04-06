/********************************************************************************
** Form generated from reading UI file 'WindowConfigureOverlay.ui'
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
#include <QtGui/QDialogButtonBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSlider>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "WidgetHistogram.h"
#include "qtcolorpicker.h"

QT_BEGIN_NAMESPACE

class Ui_WindowConfigureOverlay
{
public:
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout_4;
    QGroupBox *groupBox;
    QHBoxLayout *horizontalLayout;
    QSlider *sliderOpacity;
    QDoubleSpinBox *doubleSpinBoxOpacity;
    QGroupBox *groupBox_5;
    QHBoxLayout *horizontalLayout_5;
    QCheckBox *checkBoxEnableSmooth;
    QSpinBox *spinBoxSmoothSteps;
    QLabel *label_2;
    QSpacerItem *horizontalSpacer;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout;
    QRadioButton *radioButtonHeat;
    QRadioButton *radioButtonColorWheel;
    QRadioButton *radioButtonCustom;
    QSpacerItem *horizontalSpacer_10;
    QPushButton *pushButtonFlip;
    QCheckBox *checkBoxTruncate;
    QCheckBox *checkBoxInverse;
    QCheckBox *checkBoxClearLower;
    QCheckBox *checkBoxClearHigher;
    QGroupBox *groupBox_3;
    QVBoxLayout *verticalLayout;
    WidgetHistogram *widgetHistogram;
    QWidget *widgetHolderThreshold;
    QHBoxLayout *horizontalLayout_2;
    QSpacerItem *horizontalSpacer_2;
    QLabel *labelMin;
    QLineEdit *lineEditMin;
    QSpacerItem *horizontalSpacer_4;
    QLabel *labelMid;
    QLineEdit *lineEditMid;
    QSpacerItem *horizontalSpacer_5;
    QLabel *labelMax;
    QLineEdit *lineEditMax;
    QSpacerItem *horizontalSpacer_3;
    QWidget *widgetHolderAddPoint;
    QHBoxLayout *horizontalLayout_6;
    QSpacerItem *horizontalSpacer_8;
    QLabel *label;
    QLineEdit *lineEditNewPoint;
    QtColorPicker *widgetColorPicker;
    QPushButton *pushButtonAdd;
    QSpacerItem *horizontalSpacer_6;
    QGroupBox *groupBox_4;
    QHBoxLayout *horizontalLayout_3;
    QRadioButton *radioButtonLinear;
    QRadioButton *radioButtonLinearOpaque;
    QRadioButton *radioButtonPiecewise;
    QSpacerItem *horizontalSpacer_9;
    QLabel *label_4;
    QDialogButtonBox *buttonBox;

    void setupUi(QWidget *WindowConfigureOverlay)
    {
        if (WindowConfigureOverlay->objectName().isEmpty())
            WindowConfigureOverlay->setObjectName(QString::fromUtf8("WindowConfigureOverlay"));
        WindowConfigureOverlay->resize(490, 665);
        verticalLayout_2 = new QVBoxLayout(WindowConfigureOverlay);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setSpacing(26);
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        horizontalLayout_4->setContentsMargins(0, 0, -1, -1);
        groupBox = new QGroupBox(WindowConfigureOverlay);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        horizontalLayout = new QHBoxLayout(groupBox);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        sliderOpacity = new QSlider(groupBox);
        sliderOpacity->setObjectName(QString::fromUtf8("sliderOpacity"));
        sliderOpacity->setMaximum(100);
        sliderOpacity->setValue(100);
        sliderOpacity->setOrientation(Qt::Horizontal);

        horizontalLayout->addWidget(sliderOpacity);

        doubleSpinBoxOpacity = new QDoubleSpinBox(groupBox);
        doubleSpinBoxOpacity->setObjectName(QString::fromUtf8("doubleSpinBoxOpacity"));
        doubleSpinBoxOpacity->setMaximum(1);
        doubleSpinBoxOpacity->setSingleStep(0.05);
        doubleSpinBoxOpacity->setValue(1);

        horizontalLayout->addWidget(doubleSpinBoxOpacity);


        horizontalLayout_4->addWidget(groupBox);

        groupBox_5 = new QGroupBox(WindowConfigureOverlay);
        groupBox_5->setObjectName(QString::fromUtf8("groupBox_5"));
        horizontalLayout_5 = new QHBoxLayout(groupBox_5);
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        checkBoxEnableSmooth = new QCheckBox(groupBox_5);
        checkBoxEnableSmooth->setObjectName(QString::fromUtf8("checkBoxEnableSmooth"));

        horizontalLayout_5->addWidget(checkBoxEnableSmooth);

        spinBoxSmoothSteps = new QSpinBox(groupBox_5);
        spinBoxSmoothSteps->setObjectName(QString::fromUtf8("spinBoxSmoothSteps"));
        spinBoxSmoothSteps->setEnabled(false);
        spinBoxSmoothSteps->setMinimum(1);
        spinBoxSmoothSteps->setMaximum(100);

        horizontalLayout_5->addWidget(spinBoxSmoothSteps);

        label_2 = new QLabel(groupBox_5);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setEnabled(false);

        horizontalLayout_5->addWidget(label_2);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_5->addItem(horizontalSpacer);


        horizontalLayout_4->addWidget(groupBox_5);


        verticalLayout_2->addLayout(horizontalLayout_4);

        groupBox_2 = new QGroupBox(WindowConfigureOverlay);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        gridLayout = new QGridLayout(groupBox_2);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        radioButtonHeat = new QRadioButton(groupBox_2);
        radioButtonHeat->setObjectName(QString::fromUtf8("radioButtonHeat"));
        radioButtonHeat->setChecked(true);

        gridLayout->addWidget(radioButtonHeat, 0, 0, 1, 1);

        radioButtonColorWheel = new QRadioButton(groupBox_2);
        radioButtonColorWheel->setObjectName(QString::fromUtf8("radioButtonColorWheel"));

        gridLayout->addWidget(radioButtonColorWheel, 0, 1, 1, 1);

        radioButtonCustom = new QRadioButton(groupBox_2);
        radioButtonCustom->setObjectName(QString::fromUtf8("radioButtonCustom"));

        gridLayout->addWidget(radioButtonCustom, 0, 2, 1, 1);

        horizontalSpacer_10 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_10, 0, 3, 1, 1);

        pushButtonFlip = new QPushButton(groupBox_2);
        pushButtonFlip->setObjectName(QString::fromUtf8("pushButtonFlip"));
        pushButtonFlip->setMaximumSize(QSize(65, 16777215));

        gridLayout->addWidget(pushButtonFlip, 0, 4, 1, 1);

        checkBoxTruncate = new QCheckBox(groupBox_2);
        checkBoxTruncate->setObjectName(QString::fromUtf8("checkBoxTruncate"));

        gridLayout->addWidget(checkBoxTruncate, 1, 0, 1, 1);

        checkBoxInverse = new QCheckBox(groupBox_2);
        checkBoxInverse->setObjectName(QString::fromUtf8("checkBoxInverse"));

        gridLayout->addWidget(checkBoxInverse, 1, 1, 1, 1);

        checkBoxClearLower = new QCheckBox(groupBox_2);
        checkBoxClearLower->setObjectName(QString::fromUtf8("checkBoxClearLower"));

        gridLayout->addWidget(checkBoxClearLower, 1, 2, 1, 1);

        checkBoxClearHigher = new QCheckBox(groupBox_2);
        checkBoxClearHigher->setObjectName(QString::fromUtf8("checkBoxClearHigher"));

        gridLayout->addWidget(checkBoxClearHigher, 1, 3, 1, 2);


        verticalLayout_2->addWidget(groupBox_2);

        groupBox_3 = new QGroupBox(WindowConfigureOverlay);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(1);
        sizePolicy.setHeightForWidth(groupBox_3->sizePolicy().hasHeightForWidth());
        groupBox_3->setSizePolicy(sizePolicy);
        verticalLayout = new QVBoxLayout(groupBox_3);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        widgetHistogram = new WidgetHistogram(groupBox_3);
        widgetHistogram->setObjectName(QString::fromUtf8("widgetHistogram"));
        sizePolicy.setHeightForWidth(widgetHistogram->sizePolicy().hasHeightForWidth());
        widgetHistogram->setSizePolicy(sizePolicy);
        widgetHistogram->setMinimumSize(QSize(0, 200));

        verticalLayout->addWidget(widgetHistogram);

        widgetHolderThreshold = new QWidget(groupBox_3);
        widgetHolderThreshold->setObjectName(QString::fromUtf8("widgetHolderThreshold"));
        horizontalLayout_2 = new QHBoxLayout(widgetHolderThreshold);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalSpacer_2 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_2);

        labelMin = new QLabel(widgetHolderThreshold);
        labelMin->setObjectName(QString::fromUtf8("labelMin"));
        labelMin->setStyleSheet(QString::fromUtf8("color:red;"));

        horizontalLayout_2->addWidget(labelMin);

        lineEditMin = new QLineEdit(widgetHolderThreshold);
        lineEditMin->setObjectName(QString::fromUtf8("lineEditMin"));
        lineEditMin->setMaximumSize(QSize(80, 16777215));
        lineEditMin->setStyleSheet(QString::fromUtf8("color:red;"));

        horizontalLayout_2->addWidget(lineEditMin);

        horizontalSpacer_4 = new QSpacerItem(10, 0, QSizePolicy::Fixed, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_4);

        labelMid = new QLabel(widgetHolderThreshold);
        labelMid->setObjectName(QString::fromUtf8("labelMid"));
        labelMid->setEnabled(false);
        labelMid->setStyleSheet(QString::fromUtf8("color:blue;"));

        horizontalLayout_2->addWidget(labelMid);

        lineEditMid = new QLineEdit(widgetHolderThreshold);
        lineEditMid->setObjectName(QString::fromUtf8("lineEditMid"));
        lineEditMid->setEnabled(false);
        lineEditMid->setMaximumSize(QSize(80, 16777215));
        lineEditMid->setStyleSheet(QString::fromUtf8("color:blue;"));

        horizontalLayout_2->addWidget(lineEditMid);

        horizontalSpacer_5 = new QSpacerItem(10, 0, QSizePolicy::Fixed, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_5);

        labelMax = new QLabel(widgetHolderThreshold);
        labelMax->setObjectName(QString::fromUtf8("labelMax"));
        labelMax->setStyleSheet(QString::fromUtf8("color:green;"));

        horizontalLayout_2->addWidget(labelMax);

        lineEditMax = new QLineEdit(widgetHolderThreshold);
        lineEditMax->setObjectName(QString::fromUtf8("lineEditMax"));
        lineEditMax->setMaximumSize(QSize(80, 16777215));
        lineEditMax->setStyleSheet(QString::fromUtf8("color:green;"));

        horizontalLayout_2->addWidget(lineEditMax);

        horizontalSpacer_3 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_3);


        verticalLayout->addWidget(widgetHolderThreshold);

        widgetHolderAddPoint = new QWidget(groupBox_3);
        widgetHolderAddPoint->setObjectName(QString::fromUtf8("widgetHolderAddPoint"));
        horizontalLayout_6 = new QHBoxLayout(widgetHolderAddPoint);
        horizontalLayout_6->setContentsMargins(0, 0, 0, 0);
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        horizontalSpacer_8 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_6->addItem(horizontalSpacer_8);

        label = new QLabel(widgetHolderAddPoint);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout_6->addWidget(label);

        lineEditNewPoint = new QLineEdit(widgetHolderAddPoint);
        lineEditNewPoint->setObjectName(QString::fromUtf8("lineEditNewPoint"));
        lineEditNewPoint->setMaximumSize(QSize(80, 16777215));

        horizontalLayout_6->addWidget(lineEditNewPoint);

        widgetColorPicker = new QtColorPicker(widgetHolderAddPoint);
        widgetColorPicker->setObjectName(QString::fromUtf8("widgetColorPicker"));

        horizontalLayout_6->addWidget(widgetColorPicker);

        pushButtonAdd = new QPushButton(widgetHolderAddPoint);
        pushButtonAdd->setObjectName(QString::fromUtf8("pushButtonAdd"));

        horizontalLayout_6->addWidget(pushButtonAdd);

        horizontalSpacer_6 = new QSpacerItem(86, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_6->addItem(horizontalSpacer_6);


        verticalLayout->addWidget(widgetHolderAddPoint);

        groupBox_4 = new QGroupBox(groupBox_3);
        groupBox_4->setObjectName(QString::fromUtf8("groupBox_4"));
        horizontalLayout_3 = new QHBoxLayout(groupBox_4);
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        radioButtonLinear = new QRadioButton(groupBox_4);
        radioButtonLinear->setObjectName(QString::fromUtf8("radioButtonLinear"));

        horizontalLayout_3->addWidget(radioButtonLinear);

        radioButtonLinearOpaque = new QRadioButton(groupBox_4);
        radioButtonLinearOpaque->setObjectName(QString::fromUtf8("radioButtonLinearOpaque"));
        radioButtonLinearOpaque->setChecked(true);

        horizontalLayout_3->addWidget(radioButtonLinearOpaque);

        radioButtonPiecewise = new QRadioButton(groupBox_4);
        radioButtonPiecewise->setObjectName(QString::fromUtf8("radioButtonPiecewise"));

        horizontalLayout_3->addWidget(radioButtonPiecewise);

        horizontalSpacer_9 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_9);


        verticalLayout->addWidget(groupBox_4);


        verticalLayout_2->addWidget(groupBox_3);

        label_4 = new QLabel(WindowConfigureOverlay);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setAlignment(Qt::AlignCenter);

        verticalLayout_2->addWidget(label_4);

        buttonBox = new QDialogButtonBox(WindowConfigureOverlay);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setStandardButtons(QDialogButtonBox::Apply|QDialogButtonBox::Close|QDialogButtonBox::Help);

        verticalLayout_2->addWidget(buttonBox);


        retranslateUi(WindowConfigureOverlay);
        QObject::connect(buttonBox, SIGNAL(rejected()), WindowConfigureOverlay, SLOT(close()));
        QObject::connect(sliderOpacity, SIGNAL(valueChanged(int)), WindowConfigureOverlay, SLOT(OnSliderOpacity(int)));
        QObject::connect(doubleSpinBoxOpacity, SIGNAL(valueChanged(double)), WindowConfigureOverlay, SLOT(OnSpinBoxOpacity(double)));
        QObject::connect(checkBoxTruncate, SIGNAL(toggled(bool)), WindowConfigureOverlay, SLOT(UpdateGraph()));
        QObject::connect(checkBoxInverse, SIGNAL(toggled(bool)), WindowConfigureOverlay, SLOT(UpdateGraph()));
        QObject::connect(radioButtonHeat, SIGNAL(toggled(bool)), WindowConfigureOverlay, SLOT(UpdateGraph()));
        QObject::connect(radioButtonColorWheel, SIGNAL(toggled(bool)), WindowConfigureOverlay, SLOT(UpdateGraph()));
        QObject::connect(radioButtonCustom, SIGNAL(toggled(bool)), WindowConfigureOverlay, SLOT(UpdateGraph()));
        QObject::connect(buttonBox, SIGNAL(clicked(QAbstractButton*)), WindowConfigureOverlay, SLOT(OnClicked(QAbstractButton*)));
        QObject::connect(radioButtonColorWheel, SIGNAL(toggled(bool)), checkBoxTruncate, SLOT(setDisabled(bool)));
        QObject::connect(radioButtonCustom, SIGNAL(toggled(bool)), checkBoxTruncate, SLOT(setDisabled(bool)));
        QObject::connect(radioButtonColorWheel, SIGNAL(clicked()), radioButtonLinearOpaque, SLOT(click()));
        QObject::connect(radioButtonLinear, SIGNAL(toggled(bool)), WindowConfigureOverlay, SLOT(UpdateGraph()));
        QObject::connect(radioButtonLinearOpaque, SIGNAL(toggled(bool)), WindowConfigureOverlay, SLOT(UpdateGraph()));
        QObject::connect(radioButtonPiecewise, SIGNAL(toggled(bool)), WindowConfigureOverlay, SLOT(UpdateGraph()));
        QObject::connect(radioButtonColorWheel, SIGNAL(toggled(bool)), groupBox_4, SLOT(setDisabled(bool)));
        QObject::connect(radioButtonCustom, SIGNAL(toggled(bool)), groupBox_4, SLOT(setDisabled(bool)));
        QObject::connect(radioButtonCustom, SIGNAL(clicked()), radioButtonLinearOpaque, SLOT(click()));
        QObject::connect(radioButtonCustom, SIGNAL(toggled(bool)), widgetHolderAddPoint, SLOT(setVisible(bool)));
        QObject::connect(radioButtonCustom, SIGNAL(toggled(bool)), widgetHolderThreshold, SLOT(setHidden(bool)));
        QObject::connect(pushButtonAdd, SIGNAL(clicked()), WindowConfigureOverlay, SLOT(OnButtonAdd()));
        QObject::connect(radioButtonCustom, SIGNAL(toggled(bool)), checkBoxInverse, SLOT(setDisabled(bool)));
        QObject::connect(radioButtonCustom, SIGNAL(toggled(bool)), checkBoxClearLower, SLOT(setVisible(bool)));
        QObject::connect(radioButtonCustom, SIGNAL(toggled(bool)), checkBoxClearHigher, SLOT(setVisible(bool)));
        QObject::connect(checkBoxClearLower, SIGNAL(toggled(bool)), WindowConfigureOverlay, SLOT(UpdateGraph()));
        QObject::connect(checkBoxClearHigher, SIGNAL(toggled(bool)), WindowConfigureOverlay, SLOT(UpdateGraph()));
        QObject::connect(radioButtonCustom, SIGNAL(toggled(bool)), pushButtonFlip, SLOT(setVisible(bool)));
        QObject::connect(pushButtonFlip, SIGNAL(clicked()), widgetHistogram, SLOT(FlipMarkers()));
        QObject::connect(checkBoxEnableSmooth, SIGNAL(toggled(bool)), spinBoxSmoothSteps, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxEnableSmooth, SIGNAL(toggled(bool)), label_2, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxEnableSmooth, SIGNAL(toggled(bool)), WindowConfigureOverlay, SLOT(OnSmoothChanged()));
        QObject::connect(spinBoxSmoothSteps, SIGNAL(valueChanged(int)), WindowConfigureOverlay, SLOT(OnSmoothChanged()));
        QObject::connect(lineEditMin, SIGNAL(textChanged(QString)), WindowConfigureOverlay, SLOT(OnTextThresholdChanged(QString)));
        QObject::connect(lineEditMid, SIGNAL(textChanged(QString)), WindowConfigureOverlay, SLOT(OnTextThresholdChanged(QString)));
        QObject::connect(lineEditMax, SIGNAL(textChanged(QString)), WindowConfigureOverlay, SLOT(OnTextThresholdChanged(QString)));

        QMetaObject::connectSlotsByName(WindowConfigureOverlay);
    } // setupUi

    void retranslateUi(QWidget *WindowConfigureOverlay)
    {
        WindowConfigureOverlay->setWindowTitle(QApplication::translate("WindowConfigureOverlay", "Configure Overlay", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("WindowConfigureOverlay", "Opacity", 0, QApplication::UnicodeUTF8));
        groupBox_5->setTitle(QApplication::translate("WindowConfigureOverlay", "Smooth", 0, QApplication::UnicodeUTF8));
        checkBoxEnableSmooth->setText(QApplication::translate("WindowConfigureOverlay", "Enable", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("WindowConfigureOverlay", "steps", 0, QApplication::UnicodeUTF8));
        groupBox_2->setTitle(QApplication::translate("WindowConfigureOverlay", "Color Scale", 0, QApplication::UnicodeUTF8));
        radioButtonHeat->setText(QApplication::translate("WindowConfigureOverlay", "Heat", 0, QApplication::UnicodeUTF8));
        radioButtonColorWheel->setText(QApplication::translate("WindowConfigureOverlay", "Color Wheel", 0, QApplication::UnicodeUTF8));
        radioButtonCustom->setText(QApplication::translate("WindowConfigureOverlay", "Custom", 0, QApplication::UnicodeUTF8));
        pushButtonFlip->setText(QApplication::translate("WindowConfigureOverlay", "Flip", 0, QApplication::UnicodeUTF8));
        checkBoxTruncate->setText(QApplication::translate("WindowConfigureOverlay", "Truncate", 0, QApplication::UnicodeUTF8));
        checkBoxInverse->setText(QApplication::translate("WindowConfigureOverlay", "Inverse", 0, QApplication::UnicodeUTF8));
        checkBoxClearLower->setText(QApplication::translate("WindowConfigureOverlay", "Clear Lower", 0, QApplication::UnicodeUTF8));
        checkBoxClearHigher->setText(QApplication::translate("WindowConfigureOverlay", "Clear Higher", 0, QApplication::UnicodeUTF8));
        groupBox_3->setTitle(QApplication::translate("WindowConfigureOverlay", "Threshold", 0, QApplication::UnicodeUTF8));
        labelMin->setText(QApplication::translate("WindowConfigureOverlay", "Min", 0, QApplication::UnicodeUTF8));
        labelMid->setText(QApplication::translate("WindowConfigureOverlay", "Mid", 0, QApplication::UnicodeUTF8));
        labelMax->setText(QApplication::translate("WindowConfigureOverlay", "Max", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("WindowConfigureOverlay", "New point", 0, QApplication::UnicodeUTF8));
        pushButtonAdd->setText(QApplication::translate("WindowConfigureOverlay", "Add", 0, QApplication::UnicodeUTF8));
        groupBox_4->setTitle(QApplication::translate("WindowConfigureOverlay", "Method", 0, QApplication::UnicodeUTF8));
        radioButtonLinear->setText(QApplication::translate("WindowConfigureOverlay", "Linear", 0, QApplication::UnicodeUTF8));
        radioButtonLinearOpaque->setText(QApplication::translate("WindowConfigureOverlay", "Linear Opaque", 0, QApplication::UnicodeUTF8));
        radioButtonPiecewise->setText(QApplication::translate("WindowConfigureOverlay", "Piecewise", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("WindowConfigureOverlay", "You must click \"Apply\" to apply the changes to the overlay.", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class WindowConfigureOverlay: public Ui_WindowConfigureOverlay {};
} // namespace Ui

QT_END_NAMESPACE

