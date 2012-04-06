/********************************************************************************
** Form generated from reading UI file 'DialogSavePointSet.ui'
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

class Ui_DialogSavePointSet
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
    QRadioButton *radioButtonWayPoint;
    QRadioButton *radioButtonControlPoint;
    QSpacerItem *horizontalSpacer;
    QSpacerItem *verticalSpacer;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogSavePointSet)
    {
        if (DialogSavePointSet->objectName().isEmpty())
            DialogSavePointSet->setObjectName(QString::fromUtf8("DialogSavePointSet"));
        DialogSavePointSet->resize(400, 217);
        verticalLayout = new QVBoxLayout(DialogSavePointSet);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        label = new QLabel(DialogSavePointSet);
        label->setObjectName(QString::fromUtf8("label"));

        verticalLayout->addWidget(label);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        lineEditFileName = new QLineEdit(DialogSavePointSet);
        lineEditFileName->setObjectName(QString::fromUtf8("lineEditFileName"));
        lineEditFileName->setMinimumSize(QSize(300, 0));

        horizontalLayout->addWidget(lineEditFileName);

        toolButtonOpen = new QToolButton(DialogSavePointSet);
        toolButtonOpen->setObjectName(QString::fromUtf8("toolButtonOpen"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/file_open_16.png"), QSize(), QIcon::Normal, QIcon::Off);
        toolButtonOpen->setIcon(icon);

        horizontalLayout->addWidget(toolButtonOpen);


        verticalLayout->addLayout(horizontalLayout);

        verticalSpacer_2 = new QSpacerItem(20, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer_2);

        label_2 = new QLabel(DialogSavePointSet);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        verticalLayout->addWidget(label_2);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        radioButtonWayPoint = new QRadioButton(DialogSavePointSet);
        radioButtonWayPoint->setObjectName(QString::fromUtf8("radioButtonWayPoint"));
        radioButtonWayPoint->setChecked(true);

        horizontalLayout_2->addWidget(radioButtonWayPoint);

        radioButtonControlPoint = new QRadioButton(DialogSavePointSet);
        radioButtonControlPoint->setObjectName(QString::fromUtf8("radioButtonControlPoint"));

        horizontalLayout_2->addWidget(radioButtonControlPoint);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer);


        verticalLayout->addLayout(horizontalLayout_2);

        verticalSpacer = new QSpacerItem(20, 20, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        buttonBox = new QDialogButtonBox(DialogSavePointSet);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogSavePointSet);
        QObject::connect(toolButtonOpen, SIGNAL(clicked()), DialogSavePointSet, SLOT(OnOpen()));
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogSavePointSet, SLOT(OnOK()));
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogSavePointSet, SLOT(reject()));

        QMetaObject::connectSlotsByName(DialogSavePointSet);
    } // setupUi

    void retranslateUi(QDialog *DialogSavePointSet)
    {
        DialogSavePointSet->setWindowTitle(QApplication::translate("DialogSavePointSet", "Save Point Set", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogSavePointSet", "Enter file name", 0, QApplication::UnicodeUTF8));
        toolButtonOpen->setText(QString());
        label_2->setText(QApplication::translate("DialogSavePointSet", "Save as", 0, QApplication::UnicodeUTF8));
        radioButtonWayPoint->setText(QApplication::translate("DialogSavePointSet", "Way Points", 0, QApplication::UnicodeUTF8));
        radioButtonControlPoint->setText(QApplication::translate("DialogSavePointSet", "Control Points", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogSavePointSet: public Ui_DialogSavePointSet {};
} // namespace Ui

QT_END_NAMESPACE

