/********************************************************************************
** Form generated from reading UI file 'DialogLoadPointSet.ui'
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
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QRadioButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QToolButton>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogLoadPointSet
{
public:
    QVBoxLayout *verticalLayout;
    QLabel *label;
    QHBoxLayout *horizontalLayout;
    QLineEdit *lineEditFileName;
    QToolButton *toolButtonOpen;
    QSpacerItem *verticalSpacer_2;
    QLabel *label_2;
    QHBoxLayout *horizontalLayout_2;
    QRadioButton *radioButtonAuto;
    QRadioButton *radioButtonWayPoint;
    QRadioButton *radioButtonControlPoint;
    QSpacerItem *horizontalSpacer;
    QSpacerItem *verticalSpacer;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogLoadPointSet)
    {
        if (DialogLoadPointSet->objectName().isEmpty())
            DialogLoadPointSet->setObjectName(QString::fromUtf8("DialogLoadPointSet"));
        DialogLoadPointSet->resize(352, 194);
        verticalLayout = new QVBoxLayout(DialogLoadPointSet);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        label = new QLabel(DialogLoadPointSet);
        label->setObjectName(QString::fromUtf8("label"));

        verticalLayout->addWidget(label);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        lineEditFileName = new QLineEdit(DialogLoadPointSet);
        lineEditFileName->setObjectName(QString::fromUtf8("lineEditFileName"));
        lineEditFileName->setMinimumSize(QSize(300, 0));

        horizontalLayout->addWidget(lineEditFileName);

        toolButtonOpen = new QToolButton(DialogLoadPointSet);
        toolButtonOpen->setObjectName(QString::fromUtf8("toolButtonOpen"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/file_open_16.png"), QSize(), QIcon::Normal, QIcon::Off);
        toolButtonOpen->setIcon(icon);

        horizontalLayout->addWidget(toolButtonOpen);


        verticalLayout->addLayout(horizontalLayout);

        verticalSpacer_2 = new QSpacerItem(20, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer_2);

        label_2 = new QLabel(DialogLoadPointSet);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        verticalLayout->addWidget(label_2);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        radioButtonAuto = new QRadioButton(DialogLoadPointSet);
        radioButtonAuto->setObjectName(QString::fromUtf8("radioButtonAuto"));
        radioButtonAuto->setChecked(true);

        horizontalLayout_2->addWidget(radioButtonAuto);

        radioButtonWayPoint = new QRadioButton(DialogLoadPointSet);
        radioButtonWayPoint->setObjectName(QString::fromUtf8("radioButtonWayPoint"));

        horizontalLayout_2->addWidget(radioButtonWayPoint);

        radioButtonControlPoint = new QRadioButton(DialogLoadPointSet);
        radioButtonControlPoint->setObjectName(QString::fromUtf8("radioButtonControlPoint"));

        horizontalLayout_2->addWidget(radioButtonControlPoint);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer);


        verticalLayout->addLayout(horizontalLayout_2);

        verticalSpacer = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        buttonBox = new QDialogButtonBox(DialogLoadPointSet);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogLoadPointSet);
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogLoadPointSet, SLOT(reject()));
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogLoadPointSet, SLOT(OnOK()));
        QObject::connect(toolButtonOpen, SIGNAL(clicked()), DialogLoadPointSet, SLOT(OnButtonOpen()));

        QMetaObject::connectSlotsByName(DialogLoadPointSet);
    } // setupUi

    void retranslateUi(QDialog *DialogLoadPointSet)
    {
        DialogLoadPointSet->setWindowTitle(QApplication::translate("DialogLoadPointSet", "Load Point Set", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogLoadPointSet", "Enter file name", 0, QApplication::UnicodeUTF8));
        toolButtonOpen->setText(QString());
        label_2->setText(QApplication::translate("DialogLoadPointSet", "Load as", 0, QApplication::UnicodeUTF8));
        radioButtonAuto->setText(QApplication::translate("DialogLoadPointSet", "Auto", 0, QApplication::UnicodeUTF8));
        radioButtonWayPoint->setText(QApplication::translate("DialogLoadPointSet", "Way Points", 0, QApplication::UnicodeUTF8));
        radioButtonControlPoint->setText(QApplication::translate("DialogLoadPointSet", "Control Points", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogLoadPointSet: public Ui_DialogLoadPointSet {};
} // namespace Ui

QT_END_NAMESPACE

