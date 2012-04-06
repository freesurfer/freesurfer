/********************************************************************************
** Form generated from reading UI file 'DialogSaveVolume.ui'
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
#include <QtGui/QLineEdit>
#include <QtGui/QSpacerItem>
#include <QtGui/QToolButton>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogSaveVolume
{
public:
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QLineEdit *lineEditFileName;
    QToolButton *toolButtonOpen;
    QCheckBox *checkBoxNoResample;
    QCheckBox *checkBoxCrop;
    QSpacerItem *verticalSpacer;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogSaveVolume)
    {
        if (DialogSaveVolume->objectName().isEmpty())
            DialogSaveVolume->setObjectName(QString::fromUtf8("DialogSaveVolume"));
        DialogSaveVolume->resize(401, 163);
        verticalLayout = new QVBoxLayout(DialogSaveVolume);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        lineEditFileName = new QLineEdit(DialogSaveVolume);
        lineEditFileName->setObjectName(QString::fromUtf8("lineEditFileName"));
        lineEditFileName->setMinimumSize(QSize(300, 0));

        horizontalLayout->addWidget(lineEditFileName);

        toolButtonOpen = new QToolButton(DialogSaveVolume);
        toolButtonOpen->setObjectName(QString::fromUtf8("toolButtonOpen"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/file_open_16.png"), QSize(), QIcon::Normal, QIcon::Off);
        toolButtonOpen->setIcon(icon);

        horizontalLayout->addWidget(toolButtonOpen);


        verticalLayout->addLayout(horizontalLayout);

        checkBoxNoResample = new QCheckBox(DialogSaveVolume);
        checkBoxNoResample->setObjectName(QString::fromUtf8("checkBoxNoResample"));

        verticalLayout->addWidget(checkBoxNoResample);

        checkBoxCrop = new QCheckBox(DialogSaveVolume);
        checkBoxCrop->setObjectName(QString::fromUtf8("checkBoxCrop"));

        verticalLayout->addWidget(checkBoxCrop);

        verticalSpacer = new QSpacerItem(20, 16, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer);

        buttonBox = new QDialogButtonBox(DialogSaveVolume);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogSaveVolume);
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogSaveVolume, SLOT(reject()));
        QObject::connect(toolButtonOpen, SIGNAL(clicked()), DialogSaveVolume, SLOT(OnOpen()));
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogSaveVolume, SLOT(OnOK()));
        QObject::connect(checkBoxNoResample, SIGNAL(toggled(bool)), checkBoxCrop, SLOT(setHidden(bool)));

        QMetaObject::connectSlotsByName(DialogSaveVolume);
    } // setupUi

    void retranslateUi(QDialog *DialogSaveVolume)
    {
        DialogSaveVolume->setWindowTitle(QApplication::translate("DialogSaveVolume", "Save Volume As", 0, QApplication::UnicodeUTF8));
        toolButtonOpen->setText(QApplication::translate("DialogSaveVolume", "...", 0, QApplication::UnicodeUTF8));
        checkBoxNoResample->setText(QApplication::translate("DialogSaveVolume", "Do not resample when saving", 0, QApplication::UnicodeUTF8));
        checkBoxCrop->setText(QApplication::translate("DialogSaveVolume", "Crop to original volume size", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogSaveVolume: public Ui_DialogSaveVolume {};
} // namespace Ui

QT_END_NAMESPACE

