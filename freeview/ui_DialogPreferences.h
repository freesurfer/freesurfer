/********************************************************************************
** Form generated from reading UI file 'DialogPreferences.ui'
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
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QRadioButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>
#include "qtcolorpicker.h"

QT_BEGIN_NAMESPACE

class Ui_DialogPreferences
{
public:
    QVBoxLayout *verticalLayout;
    QGroupBox *groupBox;
    QGridLayout *gridLayout;
    QLabel *label;
    QtColorPicker *colorPickerBackground;
    QSpacerItem *horizontalSpacer;
    QLabel *label_2;
    QtColorPicker *colorPickerCursor;
    QSpacerItem *horizontalSpacer_2;
    QLabel *label_3;
    QComboBox *comboBoxCursorStyle;
    QSpacerItem *horizontalSpacer_3;
    QGroupBox *groupBox_2;
    QVBoxLayout *verticalLayout_2;
    QCheckBox *checkBoxSaveCopy;
    QGroupBox *groupBox_3;
    QVBoxLayout *verticalLayout_3;
    QCheckBox *checkBoxSyncZoom;
    QGroupBox *groupBox_4;
    QHBoxLayout *horizontalLayout;
    QRadioButton *radioButtonThemeDark;
    QRadioButton *radioButtonThemeLight;
    QSpacerItem *horizontalSpacer_4;
    QGroupBox *groupBoxMac;
    QVBoxLayout *verticalLayout_4;
    QCheckBox *checkBoxCommandKey;
    QCheckBox *checkBoxMacUnified;
    QSpacerItem *verticalSpacer;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogPreferences)
    {
        if (DialogPreferences->objectName().isEmpty())
            DialogPreferences->setObjectName(QString::fromUtf8("DialogPreferences"));
        DialogPreferences->setWindowModality(Qt::NonModal);
        DialogPreferences->resize(386, 464);
        DialogPreferences->setModal(false);
        verticalLayout = new QVBoxLayout(DialogPreferences);
        verticalLayout->setContentsMargins(12, 12, 12, 12);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        groupBox = new QGroupBox(DialogPreferences);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        gridLayout = new QGridLayout(groupBox);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        label = new QLabel(groupBox);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        colorPickerBackground = new QtColorPicker(groupBox);
        colorPickerBackground->setObjectName(QString::fromUtf8("colorPickerBackground"));

        gridLayout->addWidget(colorPickerBackground, 0, 1, 1, 1);

        horizontalSpacer = new QSpacerItem(118, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer, 0, 2, 1, 1);

        label_2 = new QLabel(groupBox);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 1, 0, 1, 1);

        colorPickerCursor = new QtColorPicker(groupBox);
        colorPickerCursor->setObjectName(QString::fromUtf8("colorPickerCursor"));

        gridLayout->addWidget(colorPickerCursor, 1, 1, 1, 1);

        horizontalSpacer_2 = new QSpacerItem(118, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_2, 1, 2, 1, 1);

        label_3 = new QLabel(groupBox);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout->addWidget(label_3, 2, 0, 1, 1);

        comboBoxCursorStyle = new QComboBox(groupBox);
        comboBoxCursorStyle->setObjectName(QString::fromUtf8("comboBoxCursorStyle"));

        gridLayout->addWidget(comboBoxCursorStyle, 2, 1, 1, 1);

        horizontalSpacer_3 = new QSpacerItem(118, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_3, 2, 2, 1, 1);


        verticalLayout->addWidget(groupBox);

        groupBox_2 = new QGroupBox(DialogPreferences);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        verticalLayout_2 = new QVBoxLayout(groupBox_2);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        checkBoxSaveCopy = new QCheckBox(groupBox_2);
        checkBoxSaveCopy->setObjectName(QString::fromUtf8("checkBoxSaveCopy"));
        checkBoxSaveCopy->setChecked(true);

        verticalLayout_2->addWidget(checkBoxSaveCopy);


        verticalLayout->addWidget(groupBox_2);

        groupBox_3 = new QGroupBox(DialogPreferences);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        groupBox_3->setFlat(false);
        verticalLayout_3 = new QVBoxLayout(groupBox_3);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        checkBoxSyncZoom = new QCheckBox(groupBox_3);
        checkBoxSyncZoom->setObjectName(QString::fromUtf8("checkBoxSyncZoom"));
        checkBoxSyncZoom->setChecked(true);

        verticalLayout_3->addWidget(checkBoxSyncZoom);


        verticalLayout->addWidget(groupBox_3);

        groupBox_4 = new QGroupBox(DialogPreferences);
        groupBox_4->setObjectName(QString::fromUtf8("groupBox_4"));
        horizontalLayout = new QHBoxLayout(groupBox_4);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        radioButtonThemeDark = new QRadioButton(groupBox_4);
        radioButtonThemeDark->setObjectName(QString::fromUtf8("radioButtonThemeDark"));
        radioButtonThemeDark->setChecked(true);

        horizontalLayout->addWidget(radioButtonThemeDark);

        radioButtonThemeLight = new QRadioButton(groupBox_4);
        radioButtonThemeLight->setObjectName(QString::fromUtf8("radioButtonThemeLight"));

        horizontalLayout->addWidget(radioButtonThemeLight);

        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_4);


        verticalLayout->addWidget(groupBox_4);

        groupBoxMac = new QGroupBox(DialogPreferences);
        groupBoxMac->setObjectName(QString::fromUtf8("groupBoxMac"));
        verticalLayout_4 = new QVBoxLayout(groupBoxMac);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        checkBoxCommandKey = new QCheckBox(groupBoxMac);
        checkBoxCommandKey->setObjectName(QString::fromUtf8("checkBoxCommandKey"));

        verticalLayout_4->addWidget(checkBoxCommandKey);

        checkBoxMacUnified = new QCheckBox(groupBoxMac);
        checkBoxMacUnified->setObjectName(QString::fromUtf8("checkBoxMacUnified"));

        verticalLayout_4->addWidget(checkBoxMacUnified);


        verticalLayout->addWidget(groupBoxMac);

        verticalSpacer = new QSpacerItem(20, 5, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        buttonBox = new QDialogButtonBox(DialogPreferences);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Close|QDialogButtonBox::Reset);
        buttonBox->setCenterButtons(false);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogPreferences);
        QObject::connect(buttonBox, SIGNAL(clicked(QAbstractButton*)), DialogPreferences, SLOT(OnClicked(QAbstractButton*)));
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogPreferences, SLOT(reject()));

        QMetaObject::connectSlotsByName(DialogPreferences);
    } // setupUi

    void retranslateUi(QDialog *DialogPreferences)
    {
        DialogPreferences->setWindowTitle(QApplication::translate("DialogPreferences", "Preferences", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("DialogPreferences", "View", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogPreferences", "Background color", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("DialogPreferences", "Cursor color", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("DialogPreferences", "Cursor style", 0, QApplication::UnicodeUTF8));
        comboBoxCursorStyle->clear();
        comboBoxCursorStyle->insertItems(0, QStringList()
         << QApplication::translate("DialogPreferences", "Short", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogPreferences", "Long", 0, QApplication::UnicodeUTF8)
        );
        groupBox_2->setTitle(QApplication::translate("DialogPreferences", "Edit", 0, QApplication::UnicodeUTF8));
        checkBoxSaveCopy->setText(QApplication::translate("DialogPreferences", "Save a copy when overwriting an existing file", 0, QApplication::UnicodeUTF8));
        groupBox_3->setTitle(QApplication::translate("DialogPreferences", "2D Display", 0, QApplication::UnicodeUTF8));
        checkBoxSyncZoom->setText(QApplication::translate("DialogPreferences", "Sync zoom factor of 2D views", 0, QApplication::UnicodeUTF8));
        groupBox_4->setTitle(QApplication::translate("DialogPreferences", "Command Console", 0, QApplication::UnicodeUTF8));
        radioButtonThemeDark->setText(QApplication::translate("DialogPreferences", "Dark", 0, QApplication::UnicodeUTF8));
        radioButtonThemeLight->setText(QApplication::translate("DialogPreferences", "Light", 0, QApplication::UnicodeUTF8));
        groupBoxMac->setTitle(QApplication::translate("DialogPreferences", "Mac Settings", 0, QApplication::UnicodeUTF8));
        checkBoxCommandKey->setText(QApplication::translate("DialogPreferences", "Use Command key as equivalent Ctrl key on Linux", 0, QApplication::UnicodeUTF8));
        checkBoxMacUnified->setText(QApplication::translate("DialogPreferences", "Unified title and toolbar on main window", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogPreferences: public Ui_DialogPreferences {};
} // namespace Ui

QT_END_NAMESPACE

