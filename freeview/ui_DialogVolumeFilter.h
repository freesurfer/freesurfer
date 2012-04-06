/********************************************************************************
** Form generated from reading UI file 'DialogVolumeFilter.ui'
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
#include <QtGui/QDialogButtonBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogVolumeFilter
{
public:
    QVBoxLayout *verticalLayout;
    QGridLayout *gridLayout;
    QLabel *labelKernelSize;
    QSpinBox *spinBoxKernelSize;
    QLabel *labelSigma;
    QLineEdit *lineEditSigma;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogVolumeFilter)
    {
        if (DialogVolumeFilter->objectName().isEmpty())
            DialogVolumeFilter->setObjectName(QString::fromUtf8("DialogVolumeFilter"));
        DialogVolumeFilter->resize(176, 139);
        verticalLayout = new QVBoxLayout(DialogVolumeFilter);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        labelKernelSize = new QLabel(DialogVolumeFilter);
        labelKernelSize->setObjectName(QString::fromUtf8("labelKernelSize"));

        gridLayout->addWidget(labelKernelSize, 0, 0, 1, 1);

        spinBoxKernelSize = new QSpinBox(DialogVolumeFilter);
        spinBoxKernelSize->setObjectName(QString::fromUtf8("spinBoxKernelSize"));
        spinBoxKernelSize->setMaximumSize(QSize(55, 16777215));
        spinBoxKernelSize->setMinimum(3);
        spinBoxKernelSize->setSingleStep(1);

        gridLayout->addWidget(spinBoxKernelSize, 0, 1, 1, 1);

        labelSigma = new QLabel(DialogVolumeFilter);
        labelSigma->setObjectName(QString::fromUtf8("labelSigma"));

        gridLayout->addWidget(labelSigma, 1, 0, 1, 1);

        lineEditSigma = new QLineEdit(DialogVolumeFilter);
        lineEditSigma->setObjectName(QString::fromUtf8("lineEditSigma"));
        lineEditSigma->setMaximumSize(QSize(55, 16777215));

        gridLayout->addWidget(lineEditSigma, 1, 1, 1, 1);


        verticalLayout->addLayout(gridLayout);

        buttonBox = new QDialogButtonBox(DialogVolumeFilter);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogVolumeFilter);
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogVolumeFilter, SLOT(OnOK()));
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogVolumeFilter, SLOT(reject()));

        QMetaObject::connectSlotsByName(DialogVolumeFilter);
    } // setupUi

    void retranslateUi(QDialog *DialogVolumeFilter)
    {
        DialogVolumeFilter->setWindowTitle(QApplication::translate("DialogVolumeFilter", "Volume Filter", 0, QApplication::UnicodeUTF8));
        labelKernelSize->setText(QApplication::translate("DialogVolumeFilter", "Kernel size", 0, QApplication::UnicodeUTF8));
        labelSigma->setText(QApplication::translate("DialogVolumeFilter", "Gaussian sigma", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogVolumeFilter: public Ui_DialogVolumeFilter {};
} // namespace Ui

QT_END_NAMESPACE

