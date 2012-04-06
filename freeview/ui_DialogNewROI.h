/********************************************************************************
** Form generated from reading UI file 'DialogNewROI.ui'
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
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogNewROI
{
public:
    QVBoxLayout *verticalLayout;
    QLabel *label;
    QLineEdit *lineEditName;
    QSpacerItem *verticalSpacer_2;
    QLabel *label_2;
    QComboBox *comboBoxTemplate;
    QSpacerItem *verticalSpacer;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogNewROI)
    {
        if (DialogNewROI->objectName().isEmpty())
            DialogNewROI->setObjectName(QString::fromUtf8("DialogNewROI"));
        DialogNewROI->resize(260, 198);
        verticalLayout = new QVBoxLayout(DialogNewROI);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        label = new QLabel(DialogNewROI);
        label->setObjectName(QString::fromUtf8("label"));

        verticalLayout->addWidget(label);

        lineEditName = new QLineEdit(DialogNewROI);
        lineEditName->setObjectName(QString::fromUtf8("lineEditName"));

        verticalLayout->addWidget(lineEditName);

        verticalSpacer_2 = new QSpacerItem(0, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer_2);

        label_2 = new QLabel(DialogNewROI);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        verticalLayout->addWidget(label_2);

        comboBoxTemplate = new QComboBox(DialogNewROI);
        comboBoxTemplate->setObjectName(QString::fromUtf8("comboBoxTemplate"));

        verticalLayout->addWidget(comboBoxTemplate);

        verticalSpacer = new QSpacerItem(0, 10, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        buttonBox = new QDialogButtonBox(DialogNewROI);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogNewROI);
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogNewROI, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogNewROI, SLOT(reject()));
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogNewROI, SLOT(OnOK()));

        QMetaObject::connectSlotsByName(DialogNewROI);
    } // setupUi

    void retranslateUi(QDialog *DialogNewROI)
    {
        DialogNewROI->setWindowTitle(QApplication::translate("DialogNewROI", "New ROI", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogNewROI", "Enter the name of the new ROI", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("DialogNewROI", "Select the template volume", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogNewROI: public Ui_DialogNewROI {};
} // namespace Ui

QT_END_NAMESPACE

