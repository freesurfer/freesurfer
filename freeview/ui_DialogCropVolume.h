/********************************************************************************
** Form generated from reading UI file 'DialogCropVolume.ui'
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
#include <QtGui/QDialog>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogCropVolume
{
public:
    QVBoxLayout *verticalLayout;
    QGroupBox *groupBox;
    QGridLayout *gridLayout;
    QLabel *label_2;
    QSpinBox *spinBoxMinX;
    QSpinBox *spinBoxMaxX;
    QLabel *label_3;
    QSpinBox *spinBoxMinY;
    QSpinBox *spinBoxMaxY;
    QLabel *label_4;
    QSpinBox *spinBoxMinZ;
    QSpinBox *spinBoxMaxZ;
    QSpacerItem *verticalSpacer;
    QGridLayout *gridLayout_2;
    QPushButton *pushButtonApply;
    QPushButton *pushButtonReset;
    QPushButton *pushButtonSaveAs;

    void setupUi(QDialog *DialogCropVolume)
    {
        if (DialogCropVolume->objectName().isEmpty())
            DialogCropVolume->setObjectName(QString::fromUtf8("DialogCropVolume"));
        DialogCropVolume->resize(210, 247);
        verticalLayout = new QVBoxLayout(DialogCropVolume);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        groupBox = new QGroupBox(DialogCropVolume);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        gridLayout = new QGridLayout(groupBox);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        label_2 = new QLabel(groupBox);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 0, 0, 1, 1);

        spinBoxMinX = new QSpinBox(groupBox);
        spinBoxMinX->setObjectName(QString::fromUtf8("spinBoxMinX"));

        gridLayout->addWidget(spinBoxMinX, 0, 1, 1, 1);

        spinBoxMaxX = new QSpinBox(groupBox);
        spinBoxMaxX->setObjectName(QString::fromUtf8("spinBoxMaxX"));

        gridLayout->addWidget(spinBoxMaxX, 0, 2, 1, 1);

        label_3 = new QLabel(groupBox);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout->addWidget(label_3, 1, 0, 1, 1);

        spinBoxMinY = new QSpinBox(groupBox);
        spinBoxMinY->setObjectName(QString::fromUtf8("spinBoxMinY"));

        gridLayout->addWidget(spinBoxMinY, 1, 1, 1, 1);

        spinBoxMaxY = new QSpinBox(groupBox);
        spinBoxMaxY->setObjectName(QString::fromUtf8("spinBoxMaxY"));

        gridLayout->addWidget(spinBoxMaxY, 1, 2, 1, 1);

        label_4 = new QLabel(groupBox);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout->addWidget(label_4, 2, 0, 1, 1);

        spinBoxMinZ = new QSpinBox(groupBox);
        spinBoxMinZ->setObjectName(QString::fromUtf8("spinBoxMinZ"));

        gridLayout->addWidget(spinBoxMinZ, 2, 1, 1, 1);

        spinBoxMaxZ = new QSpinBox(groupBox);
        spinBoxMaxZ->setObjectName(QString::fromUtf8("spinBoxMaxZ"));

        gridLayout->addWidget(spinBoxMaxZ, 2, 2, 1, 1);


        verticalLayout->addWidget(groupBox);

        verticalSpacer = new QSpacerItem(0, 10, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        pushButtonApply = new QPushButton(DialogCropVolume);
        pushButtonApply->setObjectName(QString::fromUtf8("pushButtonApply"));

        gridLayout_2->addWidget(pushButtonApply, 0, 0, 1, 1);

        pushButtonReset = new QPushButton(DialogCropVolume);
        pushButtonReset->setObjectName(QString::fromUtf8("pushButtonReset"));

        gridLayout_2->addWidget(pushButtonReset, 0, 1, 1, 1);

        pushButtonSaveAs = new QPushButton(DialogCropVolume);
        pushButtonSaveAs->setObjectName(QString::fromUtf8("pushButtonSaveAs"));

        gridLayout_2->addWidget(pushButtonSaveAs, 1, 0, 1, 1);


        verticalLayout->addLayout(gridLayout_2);


        retranslateUi(DialogCropVolume);
        QObject::connect(spinBoxMinX, SIGNAL(valueChanged(int)), DialogCropVolume, SLOT(OnSpinRange(int)));
        QObject::connect(spinBoxMaxX, SIGNAL(valueChanged(int)), DialogCropVolume, SLOT(OnSpinRange(int)));
        QObject::connect(spinBoxMinY, SIGNAL(valueChanged(int)), DialogCropVolume, SLOT(OnSpinRange(int)));
        QObject::connect(spinBoxMaxY, SIGNAL(valueChanged(int)), DialogCropVolume, SLOT(OnSpinRange(int)));
        QObject::connect(spinBoxMinZ, SIGNAL(valueChanged(int)), DialogCropVolume, SLOT(OnSpinRange(int)));
        QObject::connect(spinBoxMaxZ, SIGNAL(valueChanged(int)), DialogCropVolume, SLOT(OnSpinRange(int)));

        QMetaObject::connectSlotsByName(DialogCropVolume);
    } // setupUi

    void retranslateUi(QDialog *DialogCropVolume)
    {
        DialogCropVolume->setWindowTitle(QApplication::translate("DialogCropVolume", "Crop Volume", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("DialogCropVolume", "Cropping range", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("DialogCropVolume", "L-R", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("DialogCropVolume", "P-A", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("DialogCropVolume", "I-S", 0, QApplication::UnicodeUTF8));
        pushButtonApply->setText(QApplication::translate("DialogCropVolume", "Apply", 0, QApplication::UnicodeUTF8));
        pushButtonReset->setText(QApplication::translate("DialogCropVolume", "Reset", 0, QApplication::UnicodeUTF8));
        pushButtonSaveAs->setText(QApplication::translate("DialogCropVolume", "Save As", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogCropVolume: public Ui_DialogCropVolume {};
} // namespace Ui

QT_END_NAMESPACE

