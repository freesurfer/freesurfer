/********************************************************************************
** Form generated from reading UI file 'DialogReloadLayer.ui'
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
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogReloadLayer
{
public:
    QVBoxLayout *verticalLayout;
    QLabel *labelMessage;
    QCheckBox *checkBoxCloseFirst;
    QSpacerItem *verticalSpacer;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogReloadLayer)
    {
        if (DialogReloadLayer->objectName().isEmpty())
            DialogReloadLayer->setObjectName(QString::fromUtf8("DialogReloadLayer"));
        DialogReloadLayer->resize(400, 141);
        verticalLayout = new QVBoxLayout(DialogReloadLayer);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        labelMessage = new QLabel(DialogReloadLayer);
        labelMessage->setObjectName(QString::fromUtf8("labelMessage"));

        verticalLayout->addWidget(labelMessage);

        checkBoxCloseFirst = new QCheckBox(DialogReloadLayer);
        checkBoxCloseFirst->setObjectName(QString::fromUtf8("checkBoxCloseFirst"));

        verticalLayout->addWidget(checkBoxCloseFirst);

        verticalSpacer = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer);

        buttonBox = new QDialogButtonBox(DialogReloadLayer);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogReloadLayer);
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogReloadLayer, SLOT(accept()));
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogReloadLayer, SLOT(reject()));

        QMetaObject::connectSlotsByName(DialogReloadLayer);
    } // setupUi

    void retranslateUi(QDialog *DialogReloadLayer)
    {
        DialogReloadLayer->setWindowTitle(QApplication::translate("DialogReloadLayer", "Dialog", 0, QApplication::UnicodeUTF8));
        labelMessage->setText(QApplication::translate("DialogReloadLayer", "Reload text", 0, QApplication::UnicodeUTF8));
        checkBoxCloseFirst->setText(QApplication::translate("DialogReloadLayer", "Close the current layer before reloading", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogReloadLayer: public Ui_DialogReloadLayer {};
} // namespace Ui

QT_END_NAMESPACE

