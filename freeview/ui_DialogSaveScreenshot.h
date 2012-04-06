/********************************************************************************
** Form generated from reading UI file 'DialogSaveScreenshot.ui'
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
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QToolButton>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogSaveScreenshot
{
public:
    QVBoxLayout *verticalLayout_2;
    QLabel *label_3;
    QHBoxLayout *horizontalLayout_2;
    QLineEdit *lineEditFileName;
    QToolButton *toolButtonOpen;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout;
    QHBoxLayout *horizontalLayout;
    QLabel *label_2;
    QSpinBox *spinBoxMagnification;
    QSpacerItem *horizontalSpacer;
    QCheckBox *checkBoxHideAnnotation;
    QCheckBox *checkBoxHideCursor;
    QCheckBox *checkBoxAntiAliasing;
    QCheckBox *checkBoxKeepWindow;
    QLabel *label;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogSaveScreenshot)
    {
        if (DialogSaveScreenshot->objectName().isEmpty())
            DialogSaveScreenshot->setObjectName(QString::fromUtf8("DialogSaveScreenshot"));
        DialogSaveScreenshot->resize(366, 322);
        DialogSaveScreenshot->setLayoutDirection(Qt::LeftToRight);
        verticalLayout_2 = new QVBoxLayout(DialogSaveScreenshot);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        verticalLayout_2->setSizeConstraint(QLayout::SetFixedSize);
        label_3 = new QLabel(DialogSaveScreenshot);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        verticalLayout_2->addWidget(label_3);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        lineEditFileName = new QLineEdit(DialogSaveScreenshot);
        lineEditFileName->setObjectName(QString::fromUtf8("lineEditFileName"));

        horizontalLayout_2->addWidget(lineEditFileName);

        toolButtonOpen = new QToolButton(DialogSaveScreenshot);
        toolButtonOpen->setObjectName(QString::fromUtf8("toolButtonOpen"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/file_open_16.png"), QSize(), QIcon::Active, QIcon::On);
        toolButtonOpen->setIcon(icon);

        horizontalLayout_2->addWidget(toolButtonOpen);


        verticalLayout_2->addLayout(horizontalLayout_2);

        groupBox = new QGroupBox(DialogSaveScreenshot);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        verticalLayout = new QVBoxLayout(groupBox);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        label_2 = new QLabel(groupBox);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        horizontalLayout->addWidget(label_2);

        spinBoxMagnification = new QSpinBox(groupBox);
        spinBoxMagnification->setObjectName(QString::fromUtf8("spinBoxMagnification"));
        spinBoxMagnification->setMinimum(1);
        spinBoxMagnification->setMaximum(10);

        horizontalLayout->addWidget(spinBoxMagnification);

        horizontalSpacer = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);


        verticalLayout->addLayout(horizontalLayout);

        checkBoxHideAnnotation = new QCheckBox(groupBox);
        checkBoxHideAnnotation->setObjectName(QString::fromUtf8("checkBoxHideAnnotation"));

        verticalLayout->addWidget(checkBoxHideAnnotation);

        checkBoxHideCursor = new QCheckBox(groupBox);
        checkBoxHideCursor->setObjectName(QString::fromUtf8("checkBoxHideCursor"));

        verticalLayout->addWidget(checkBoxHideCursor);

        checkBoxAntiAliasing = new QCheckBox(groupBox);
        checkBoxAntiAliasing->setObjectName(QString::fromUtf8("checkBoxAntiAliasing"));

        verticalLayout->addWidget(checkBoxAntiAliasing);

        checkBoxKeepWindow = new QCheckBox(groupBox);
        checkBoxKeepWindow->setObjectName(QString::fromUtf8("checkBoxKeepWindow"));

        verticalLayout->addWidget(checkBoxKeepWindow);


        verticalLayout_2->addWidget(groupBox);

        label = new QLabel(DialogSaveScreenshot);
        label->setObjectName(QString::fromUtf8("label"));
        label->setWordWrap(true);

        verticalLayout_2->addWidget(label);

        buttonBox = new QDialogButtonBox(DialogSaveScreenshot);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Save);

        verticalLayout_2->addWidget(buttonBox);


        retranslateUi(DialogSaveScreenshot);
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogSaveScreenshot, SLOT(OnSave()));
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogSaveScreenshot, SLOT(reject()));
        QObject::connect(toolButtonOpen, SIGNAL(clicked()), DialogSaveScreenshot, SLOT(OnOpen()));

        QMetaObject::connectSlotsByName(DialogSaveScreenshot);
    } // setupUi

    void retranslateUi(QDialog *DialogSaveScreenshot)
    {
        DialogSaveScreenshot->setWindowTitle(QApplication::translate("DialogSaveScreenshot", "Save Screenshots", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("DialogSaveScreenshot", "Enter file name", 0, QApplication::UnicodeUTF8));
        toolButtonOpen->setText(QString());
        groupBox->setTitle(QApplication::translate("DialogSaveScreenshot", "Options", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("DialogSaveScreenshot", "Magnification factor", 0, QApplication::UnicodeUTF8));
        checkBoxHideAnnotation->setText(QApplication::translate("DialogSaveScreenshot", "Hide annotations", 0, QApplication::UnicodeUTF8));
        checkBoxHideCursor->setText(QApplication::translate("DialogSaveScreenshot", "Hide cursor", 0, QApplication::UnicodeUTF8));
        checkBoxAntiAliasing->setText(QApplication::translate("DialogSaveScreenshot", "Enable anti-aliasing (3D view only)", 0, QApplication::UnicodeUTF8));
        checkBoxKeepWindow->setText(QApplication::translate("DialogSaveScreenshot", "Do not close this window after save", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogSaveScreenshot", "For best result, move this window out of the view that you are taking the snapshot of.", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogSaveScreenshot: public Ui_DialogSaveScreenshot {};
} // namespace Ui

QT_END_NAMESPACE

