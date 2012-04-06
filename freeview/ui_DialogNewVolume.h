/********************************************************************************
** Form generated from reading UI file 'DialogNewVolume.ui'
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
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogNewVolume
{
public:
    QVBoxLayout *verticalLayout;
    QLabel *label;
    QLineEdit *lineEditName;
    QSpacerItem *verticalSpacer_2;
    QLabel *label_2;
    QComboBox *comboBoxDataType;
    QSpacerItem *verticalSpacer_3;
    QLabel *label_3;
    QComboBox *comboBoxTemplate;
    QGridLayout *gridLayout;
    QCheckBox *checkBoxCopyData;
    QCheckBox *checkBoxDummySphere;
    QCheckBox *checkBoxDummyCube;
    QCheckBox *checkBoxMask;
    QSpacerItem *verticalSpacer;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogNewVolume)
    {
        if (DialogNewVolume->objectName().isEmpty())
            DialogNewVolume->setObjectName(QString::fromUtf8("DialogNewVolume"));
        DialogNewVolume->resize(295, 325);
        verticalLayout = new QVBoxLayout(DialogNewVolume);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        verticalLayout->setContentsMargins(15, -1, 15, -1);
        label = new QLabel(DialogNewVolume);
        label->setObjectName(QString::fromUtf8("label"));

        verticalLayout->addWidget(label);

        lineEditName = new QLineEdit(DialogNewVolume);
        lineEditName->setObjectName(QString::fromUtf8("lineEditName"));

        verticalLayout->addWidget(lineEditName);

        verticalSpacer_2 = new QSpacerItem(0, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer_2);

        label_2 = new QLabel(DialogNewVolume);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        verticalLayout->addWidget(label_2);

        comboBoxDataType = new QComboBox(DialogNewVolume);
        comboBoxDataType->setObjectName(QString::fromUtf8("comboBoxDataType"));

        verticalLayout->addWidget(comboBoxDataType);

        verticalSpacer_3 = new QSpacerItem(0, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer_3);

        label_3 = new QLabel(DialogNewVolume);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        verticalLayout->addWidget(label_3);

        comboBoxTemplate = new QComboBox(DialogNewVolume);
        comboBoxTemplate->setObjectName(QString::fromUtf8("comboBoxTemplate"));

        verticalLayout->addWidget(comboBoxTemplate);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setContentsMargins(-1, 0, -1, -1);
        checkBoxCopyData = new QCheckBox(DialogNewVolume);
        checkBoxCopyData->setObjectName(QString::fromUtf8("checkBoxCopyData"));

        gridLayout->addWidget(checkBoxCopyData, 0, 0, 1, 2);

        checkBoxDummySphere = new QCheckBox(DialogNewVolume);
        checkBoxDummySphere->setObjectName(QString::fromUtf8("checkBoxDummySphere"));

        gridLayout->addWidget(checkBoxDummySphere, 2, 0, 1, 1);

        checkBoxDummyCube = new QCheckBox(DialogNewVolume);
        checkBoxDummyCube->setObjectName(QString::fromUtf8("checkBoxDummyCube"));

        gridLayout->addWidget(checkBoxDummyCube, 2, 1, 1, 1);

        checkBoxMask = new QCheckBox(DialogNewVolume);
        checkBoxMask->setObjectName(QString::fromUtf8("checkBoxMask"));

        gridLayout->addWidget(checkBoxMask, 1, 0, 1, 2);


        verticalLayout->addLayout(gridLayout);

        verticalSpacer = new QSpacerItem(0, 10, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        buttonBox = new QDialogButtonBox(DialogNewVolume);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogNewVolume);
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogNewVolume, SLOT(reject()));
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogNewVolume, SLOT(OnOK()));
        QObject::connect(checkBoxCopyData, SIGNAL(toggled(bool)), DialogNewVolume, SLOT(OnToggleCopyVoxelData(bool)));
        QObject::connect(checkBoxCopyData, SIGNAL(toggled(bool)), DialogNewVolume, SLOT(OnToggleVoxelDataOption(bool)));
        QObject::connect(checkBoxDummySphere, SIGNAL(toggled(bool)), DialogNewVolume, SLOT(OnToggleVoxelDataOption(bool)));
        QObject::connect(checkBoxDummyCube, SIGNAL(toggled(bool)), DialogNewVolume, SLOT(OnToggleVoxelDataOption(bool)));
        QObject::connect(checkBoxMask, SIGNAL(toggled(bool)), DialogNewVolume, SLOT(OnToggleVoxelDataOption(bool)));
        QObject::connect(checkBoxMask, SIGNAL(toggled(bool)), DialogNewVolume, SLOT(OnToggleMask(bool)));

        comboBoxDataType->setCurrentIndex(1);


        QMetaObject::connectSlotsByName(DialogNewVolume);
    } // setupUi

    void retranslateUi(QDialog *DialogNewVolume)
    {
        DialogNewVolume->setWindowTitle(QApplication::translate("DialogNewVolume", "New Volume", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogNewVolume", "Enter the name of the new volume", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("DialogNewVolume", "Select data type", 0, QApplication::UnicodeUTF8));
        comboBoxDataType->clear();
        comboBoxDataType->insertItems(0, QStringList()
         << QApplication::translate("DialogNewVolume", "UCHAR", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogNewVolume", "INT (for new label volume)", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogNewVolume", "LONG", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogNewVolume", "FLOAT", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogNewVolume", "SHORT", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogNewVolume", "Same as template volume", 0, QApplication::UnicodeUTF8)
        );
        label_3->setText(QApplication::translate("DialogNewVolume", "Select template volume", 0, QApplication::UnicodeUTF8));
        checkBoxCopyData->setText(QApplication::translate("DialogNewVolume", "Copy voxel data from template", 0, QApplication::UnicodeUTF8));
        checkBoxDummySphere->setText(QApplication::translate("DialogNewVolume", "Dummy Sphere", 0, QApplication::UnicodeUTF8));
        checkBoxDummyCube->setText(QApplication::translate("DialogNewVolume", "Dummy Cube", 0, QApplication::UnicodeUTF8));
        checkBoxMask->setText(QApplication::translate("DialogNewVolume", "Create as mask", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogNewVolume: public Ui_DialogNewVolume {};
} // namespace Ui

QT_END_NAMESPACE

