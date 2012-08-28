/********************************************************************************
** Form generated from reading UI file 'DialogLoadSurfaceOverlay.ui'
**
** Created: Tue Aug 28 14:49:01 2012
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

class Ui_DialogLoadSurfaceOverlay
{
public:
    QVBoxLayout *verticalLayout;
    QLabel *label;
    QHBoxLayout *horizontalLayout;
    QLineEdit *lineEditFile;
    QToolButton *toolButton;
    QCheckBox *checkBoxRegistration;
    QHBoxLayout *horizontalLayout_2;
    QLineEdit *lineEditRegistration;
    QToolButton *toolButtonRegistration;
    QSpacerItem *verticalSpacer;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogLoadSurfaceOverlay)
    {
        if (DialogLoadSurfaceOverlay->objectName().isEmpty())
            DialogLoadSurfaceOverlay->setObjectName(QString::fromUtf8("DialogLoadSurfaceOverlay"));
        DialogLoadSurfaceOverlay->resize(357, 193);
        verticalLayout = new QVBoxLayout(DialogLoadSurfaceOverlay);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        label = new QLabel(DialogLoadSurfaceOverlay);
        label->setObjectName(QString::fromUtf8("label"));

        verticalLayout->addWidget(label);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        lineEditFile = new QLineEdit(DialogLoadSurfaceOverlay);
        lineEditFile->setObjectName(QString::fromUtf8("lineEditFile"));
        lineEditFile->setMinimumSize(QSize(300, 0));

        horizontalLayout->addWidget(lineEditFile);

        toolButton = new QToolButton(DialogLoadSurfaceOverlay);
        toolButton->setObjectName(QString::fromUtf8("toolButton"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/file_open_16.png"), QSize(), QIcon::Normal, QIcon::Off);
        toolButton->setIcon(icon);

        horizontalLayout->addWidget(toolButton);


        verticalLayout->addLayout(horizontalLayout);

        checkBoxRegistration = new QCheckBox(DialogLoadSurfaceOverlay);
        checkBoxRegistration->setObjectName(QString::fromUtf8("checkBoxRegistration"));

        verticalLayout->addWidget(checkBoxRegistration);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        lineEditRegistration = new QLineEdit(DialogLoadSurfaceOverlay);
        lineEditRegistration->setObjectName(QString::fromUtf8("lineEditRegistration"));
        lineEditRegistration->setEnabled(false);

        horizontalLayout_2->addWidget(lineEditRegistration);

        toolButtonRegistration = new QToolButton(DialogLoadSurfaceOverlay);
        toolButtonRegistration->setObjectName(QString::fromUtf8("toolButtonRegistration"));
        toolButtonRegistration->setEnabled(false);
        toolButtonRegistration->setIcon(icon);

        horizontalLayout_2->addWidget(toolButtonRegistration);


        verticalLayout->addLayout(horizontalLayout_2);

        verticalSpacer = new QSpacerItem(20, 15, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        buttonBox = new QDialogButtonBox(DialogLoadSurfaceOverlay);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);
        buttonBox->setCenterButtons(true);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogLoadSurfaceOverlay);
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogLoadSurfaceOverlay, SLOT(OnOK()));
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogLoadSurfaceOverlay, SLOT(reject()));
        QObject::connect(checkBoxRegistration, SIGNAL(toggled(bool)), lineEditRegistration, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxRegistration, SIGNAL(toggled(bool)), toolButtonRegistration, SLOT(setEnabled(bool)));
        QObject::connect(toolButton, SIGNAL(clicked()), DialogLoadSurfaceOverlay, SLOT(OnButtonOpen()));
        QObject::connect(toolButtonRegistration, SIGNAL(clicked()), DialogLoadSurfaceOverlay, SLOT(OnButtonRegistration()));

        QMetaObject::connectSlotsByName(DialogLoadSurfaceOverlay);
    } // setupUi

    void retranslateUi(QDialog *DialogLoadSurfaceOverlay)
    {
        DialogLoadSurfaceOverlay->setWindowTitle(QApplication::translate("DialogLoadSurfaceOverlay", "Dialog", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogLoadSurfaceOverlay", "Select volume file", 0, QApplication::UnicodeUTF8));
        toolButton->setText(QApplication::translate("DialogLoadSurfaceOverlay", "...", 0, QApplication::UnicodeUTF8));
        checkBoxRegistration->setText(QApplication::translate("DialogLoadSurfaceOverlay", "Apply registration", 0, QApplication::UnicodeUTF8));
        toolButtonRegistration->setText(QApplication::translate("DialogLoadSurfaceOverlay", "...", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogLoadSurfaceOverlay: public Ui_DialogLoadSurfaceOverlay {};
} // namespace Ui

QT_END_NAMESPACE

