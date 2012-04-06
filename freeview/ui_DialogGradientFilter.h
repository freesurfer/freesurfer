/********************************************************************************
** Form generated from reading UI file 'DialogGradientFilter.ui'
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
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogGradientFilter
{
public:
    QVBoxLayout *verticalLayout;
    QCheckBox *checkBoxSmooth;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QDoubleSpinBox *doubleSpinBoxSD;
    QLabel *label_2;
    QSpacerItem *horizontalSpacer;
    QSpacerItem *verticalSpacer;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogGradientFilter)
    {
        if (DialogGradientFilter->objectName().isEmpty())
            DialogGradientFilter->setObjectName(QString::fromUtf8("DialogGradientFilter"));
        DialogGradientFilter->resize(334, 143);
        verticalLayout = new QVBoxLayout(DialogGradientFilter);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        checkBoxSmooth = new QCheckBox(DialogGradientFilter);
        checkBoxSmooth->setObjectName(QString::fromUtf8("checkBoxSmooth"));

        verticalLayout->addWidget(checkBoxSmooth);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        label = new QLabel(DialogGradientFilter);
        label->setObjectName(QString::fromUtf8("label"));
        label->setEnabled(false);

        horizontalLayout->addWidget(label);

        doubleSpinBoxSD = new QDoubleSpinBox(DialogGradientFilter);
        doubleSpinBoxSD->setObjectName(QString::fromUtf8("doubleSpinBoxSD"));
        doubleSpinBoxSD->setEnabled(false);
        doubleSpinBoxSD->setMaximum(10);
        doubleSpinBoxSD->setSingleStep(0.05);
        doubleSpinBoxSD->setValue(1);

        horizontalLayout->addWidget(doubleSpinBoxSD);

        label_2 = new QLabel(DialogGradientFilter);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setEnabled(false);

        horizontalLayout->addWidget(label_2);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);


        verticalLayout->addLayout(horizontalLayout);

        verticalSpacer = new QSpacerItem(20, 26, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        buttonBox = new QDialogButtonBox(DialogGradientFilter);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogGradientFilter);
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogGradientFilter, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogGradientFilter, SLOT(reject()));
        QObject::connect(checkBoxSmooth, SIGNAL(toggled(bool)), label, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxSmooth, SIGNAL(toggled(bool)), doubleSpinBoxSD, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxSmooth, SIGNAL(toggled(bool)), label_2, SLOT(setEnabled(bool)));

        QMetaObject::connectSlotsByName(DialogGradientFilter);
    } // setupUi

    void retranslateUi(QDialog *DialogGradientFilter)
    {
        DialogGradientFilter->setWindowTitle(QApplication::translate("DialogGradientFilter", "Gradient Filter", 0, QApplication::UnicodeUTF8));
        checkBoxSmooth->setText(QApplication::translate("DialogGradientFilter", "Apply Gaussian smoothing before filtering", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogGradientFilter", "Standard deviation", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("DialogGradientFilter", "voxel(s)", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogGradientFilter: public Ui_DialogGradientFilter {};
} // namespace Ui

QT_END_NAMESPACE

