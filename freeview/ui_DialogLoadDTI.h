/********************************************************************************
** Form generated from reading UI file 'DialogLoadDTI.ui'
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
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QSpacerItem>
#include <QtGui/QToolButton>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogLoadDTI
{
public:
    QVBoxLayout *verticalLayout;
    QLabel *label;
    QHBoxLayout *horizontalLayout;
    QLineEdit *lineEditFA;
    QToolButton *toolButtonOpenFA;
    QSpacerItem *verticalSpacer_2;
    QLabel *label_2;
    QHBoxLayout *horizontalLayout_2;
    QLineEdit *lineEditVector;
    QToolButton *toolButtonOpenVector;
    QSpacerItem *verticalSpacer_3;
    QCheckBox *checkBoxResample;
    QCheckBox *checkBoxRegistration;
    QHBoxLayout *horizontalLayout_3;
    QLineEdit *lineEditRegistration;
    QToolButton *toolButtonOpenRegistration;
    QSpacerItem *verticalSpacer;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogLoadDTI)
    {
        if (DialogLoadDTI->objectName().isEmpty())
            DialogLoadDTI->setObjectName(QString::fromUtf8("DialogLoadDTI"));
        DialogLoadDTI->resize(376, 305);
        verticalLayout = new QVBoxLayout(DialogLoadDTI);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        label = new QLabel(DialogLoadDTI);
        label->setObjectName(QString::fromUtf8("label"));

        verticalLayout->addWidget(label);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        lineEditFA = new QLineEdit(DialogLoadDTI);
        lineEditFA->setObjectName(QString::fromUtf8("lineEditFA"));
        lineEditFA->setMinimumSize(QSize(320, 0));

        horizontalLayout->addWidget(lineEditFA);

        toolButtonOpenFA = new QToolButton(DialogLoadDTI);
        toolButtonOpenFA->setObjectName(QString::fromUtf8("toolButtonOpenFA"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/file_open_16.png"), QSize(), QIcon::Normal, QIcon::Off);
        toolButtonOpenFA->setIcon(icon);

        horizontalLayout->addWidget(toolButtonOpenFA);


        verticalLayout->addLayout(horizontalLayout);

        verticalSpacer_2 = new QSpacerItem(0, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer_2);

        label_2 = new QLabel(DialogLoadDTI);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        verticalLayout->addWidget(label_2);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        lineEditVector = new QLineEdit(DialogLoadDTI);
        lineEditVector->setObjectName(QString::fromUtf8("lineEditVector"));

        horizontalLayout_2->addWidget(lineEditVector);

        toolButtonOpenVector = new QToolButton(DialogLoadDTI);
        toolButtonOpenVector->setObjectName(QString::fromUtf8("toolButtonOpenVector"));
        toolButtonOpenVector->setIcon(icon);

        horizontalLayout_2->addWidget(toolButtonOpenVector);


        verticalLayout->addLayout(horizontalLayout_2);

        verticalSpacer_3 = new QSpacerItem(0, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer_3);

        checkBoxResample = new QCheckBox(DialogLoadDTI);
        checkBoxResample->setObjectName(QString::fromUtf8("checkBoxResample"));

        verticalLayout->addWidget(checkBoxResample);

        checkBoxRegistration = new QCheckBox(DialogLoadDTI);
        checkBoxRegistration->setObjectName(QString::fromUtf8("checkBoxRegistration"));

        verticalLayout->addWidget(checkBoxRegistration);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        lineEditRegistration = new QLineEdit(DialogLoadDTI);
        lineEditRegistration->setObjectName(QString::fromUtf8("lineEditRegistration"));
        lineEditRegistration->setEnabled(false);

        horizontalLayout_3->addWidget(lineEditRegistration);

        toolButtonOpenRegistration = new QToolButton(DialogLoadDTI);
        toolButtonOpenRegistration->setObjectName(QString::fromUtf8("toolButtonOpenRegistration"));
        toolButtonOpenRegistration->setEnabled(false);
        toolButtonOpenRegistration->setIcon(icon);

        horizontalLayout_3->addWidget(toolButtonOpenRegistration);


        verticalLayout->addLayout(horizontalLayout_3);

        verticalSpacer = new QSpacerItem(0, 10, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        buttonBox = new QDialogButtonBox(DialogLoadDTI);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogLoadDTI);
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogLoadDTI, SLOT(reject()));
        QObject::connect(checkBoxRegistration, SIGNAL(toggled(bool)), lineEditRegistration, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxRegistration, SIGNAL(toggled(bool)), toolButtonOpenRegistration, SLOT(setEnabled(bool)));
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogLoadDTI, SLOT(OnOK()));
        QObject::connect(toolButtonOpenRegistration, SIGNAL(clicked()), DialogLoadDTI, SLOT(OnButtonRegistration()));
        QObject::connect(toolButtonOpenFA, SIGNAL(clicked()), DialogLoadDTI, SLOT(OnButtonFA()));
        QObject::connect(toolButtonOpenVector, SIGNAL(clicked()), DialogLoadDTI, SLOT(OnButtonVector()));

        QMetaObject::connectSlotsByName(DialogLoadDTI);
    } // setupUi

    void retranslateUi(QDialog *DialogLoadDTI)
    {
        DialogLoadDTI->setWindowTitle(QApplication::translate("DialogLoadDTI", "Load DTI Volume", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogLoadDTI", "Select FA data file", 0, QApplication::UnicodeUTF8));
        toolButtonOpenFA->setText(QString());
        label_2->setText(QApplication::translate("DialogLoadDTI", "Select DTI vector file", 0, QApplication::UnicodeUTF8));
        toolButtonOpenVector->setText(QString());
        checkBoxResample->setText(QApplication::translate("DialogLoadDTI", "Resample to standard RAS space", 0, QApplication::UnicodeUTF8));
        checkBoxRegistration->setText(QApplication::translate("DialogLoadDTI", "Apply registration file", 0, QApplication::UnicodeUTF8));
        toolButtonOpenRegistration->setText(QString());
    } // retranslateUi

};

namespace Ui {
    class DialogLoadDTI: public Ui_DialogLoadDTI {};
} // namespace Ui

QT_END_NAMESPACE

