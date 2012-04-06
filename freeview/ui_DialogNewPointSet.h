/********************************************************************************
** Form generated from reading UI file 'DialogNewPointSet.ui'
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
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QRadioButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogNewPointSet
{
public:
    QVBoxLayout *verticalLayout;
    QLabel *label;
    QLineEdit *lineEditName;
    QHBoxLayout *horizontalLayout;
    QRadioButton *radioButtonControlPoint;
    QRadioButton *radioButtonWayPoint;
    QSpacerItem *horizontalSpacer;
    QSpacerItem *verticalSpacer_2;
    QLabel *label_2;
    QComboBox *comboBoxTemplate;
    QSpacerItem *verticalSpacer;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogNewPointSet)
    {
        if (DialogNewPointSet->objectName().isEmpty())
            DialogNewPointSet->setObjectName(QString::fromUtf8("DialogNewPointSet"));
        DialogNewPointSet->resize(323, 226);
        verticalLayout = new QVBoxLayout(DialogNewPointSet);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        label = new QLabel(DialogNewPointSet);
        label->setObjectName(QString::fromUtf8("label"));

        verticalLayout->addWidget(label);

        lineEditName = new QLineEdit(DialogNewPointSet);
        lineEditName->setObjectName(QString::fromUtf8("lineEditName"));

        verticalLayout->addWidget(lineEditName);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        radioButtonControlPoint = new QRadioButton(DialogNewPointSet);
        radioButtonControlPoint->setObjectName(QString::fromUtf8("radioButtonControlPoint"));
        radioButtonControlPoint->setChecked(true);

        horizontalLayout->addWidget(radioButtonControlPoint);

        radioButtonWayPoint = new QRadioButton(DialogNewPointSet);
        radioButtonWayPoint->setObjectName(QString::fromUtf8("radioButtonWayPoint"));
        radioButtonWayPoint->setChecked(false);

        horizontalLayout->addWidget(radioButtonWayPoint);

        horizontalSpacer = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);


        verticalLayout->addLayout(horizontalLayout);

        verticalSpacer_2 = new QSpacerItem(0, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer_2);

        label_2 = new QLabel(DialogNewPointSet);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        verticalLayout->addWidget(label_2);

        comboBoxTemplate = new QComboBox(DialogNewPointSet);
        comboBoxTemplate->setObjectName(QString::fromUtf8("comboBoxTemplate"));

        verticalLayout->addWidget(comboBoxTemplate);

        verticalSpacer = new QSpacerItem(0, 10, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        buttonBox = new QDialogButtonBox(DialogNewPointSet);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogNewPointSet);
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogNewPointSet, SLOT(reject()));
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogNewPointSet, SLOT(OnOK()));

        QMetaObject::connectSlotsByName(DialogNewPointSet);
    } // setupUi

    void retranslateUi(QDialog *DialogNewPointSet)
    {
        DialogNewPointSet->setWindowTitle(QApplication::translate("DialogNewPointSet", "New Point Set", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogNewPointSet", "Enter the name of the new point set", 0, QApplication::UnicodeUTF8));
        radioButtonControlPoint->setText(QApplication::translate("DialogNewPointSet", "Control points", 0, QApplication::UnicodeUTF8));
        radioButtonWayPoint->setText(QApplication::translate("DialogNewPointSet", "Way points", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("DialogNewPointSet", "Select template volume", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogNewPointSet: public Ui_DialogNewPointSet {};
} // namespace Ui

QT_END_NAMESPACE

