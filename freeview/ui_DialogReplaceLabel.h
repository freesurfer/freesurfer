/********************************************************************************
** Form generated from reading UI file 'DialogReplaceLabel.ui'
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
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogReplaceLabel
{
public:
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label;
    QLineEdit *lineEditOriginalValue;
    QSpacerItem *horizontalSpacer_3;
    QLabel *label_2;
    QLineEdit *lineEditNewValue;
    QGroupBox *groupBox;
    QVBoxLayout *verticalLayout;
    QRadioButton *radioButtonAllSlices;
    QRadioButton *radioButtonCurrentSlice;
    QSpacerItem *verticalSpacer;
    QHBoxLayout *horizontalLayout;
    QSpacerItem *horizontalSpacer;
    QPushButton *pushButtonCancel;
    QPushButton *pushButtonReplace;
    QSpacerItem *horizontalSpacer_2;

    void setupUi(QDialog *DialogReplaceLabel)
    {
        if (DialogReplaceLabel->objectName().isEmpty())
            DialogReplaceLabel->setObjectName(QString::fromUtf8("DialogReplaceLabel"));
        DialogReplaceLabel->resize(346, 191);
        verticalLayout_2 = new QVBoxLayout(DialogReplaceLabel);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        verticalLayout_2->setSizeConstraint(QLayout::SetFixedSize);
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        label = new QLabel(DialogReplaceLabel);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout_2->addWidget(label);

        lineEditOriginalValue = new QLineEdit(DialogReplaceLabel);
        lineEditOriginalValue->setObjectName(QString::fromUtf8("lineEditOriginalValue"));
        lineEditOriginalValue->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_2->addWidget(lineEditOriginalValue);

        horizontalSpacer_3 = new QSpacerItem(15, 5, QSizePolicy::Fixed, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_3);

        label_2 = new QLabel(DialogReplaceLabel);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        horizontalLayout_2->addWidget(label_2);

        lineEditNewValue = new QLineEdit(DialogReplaceLabel);
        lineEditNewValue->setObjectName(QString::fromUtf8("lineEditNewValue"));
        lineEditNewValue->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_2->addWidget(lineEditNewValue);


        verticalLayout_2->addLayout(horizontalLayout_2);

        groupBox = new QGroupBox(DialogReplaceLabel);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        verticalLayout = new QVBoxLayout(groupBox);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        radioButtonAllSlices = new QRadioButton(groupBox);
        radioButtonAllSlices->setObjectName(QString::fromUtf8("radioButtonAllSlices"));
        radioButtonAllSlices->setChecked(true);

        verticalLayout->addWidget(radioButtonAllSlices);

        radioButtonCurrentSlice = new QRadioButton(groupBox);
        radioButtonCurrentSlice->setObjectName(QString::fromUtf8("radioButtonCurrentSlice"));

        verticalLayout->addWidget(radioButtonCurrentSlice);


        verticalLayout_2->addWidget(groupBox);

        verticalSpacer = new QSpacerItem(20, 8, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_2->addItem(verticalSpacer);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);

        pushButtonCancel = new QPushButton(DialogReplaceLabel);
        pushButtonCancel->setObjectName(QString::fromUtf8("pushButtonCancel"));
        pushButtonCancel->setDefault(true);

        horizontalLayout->addWidget(pushButtonCancel);

        pushButtonReplace = new QPushButton(DialogReplaceLabel);
        pushButtonReplace->setObjectName(QString::fromUtf8("pushButtonReplace"));

        horizontalLayout->addWidget(pushButtonReplace);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_2);


        verticalLayout_2->addLayout(horizontalLayout);


        retranslateUi(DialogReplaceLabel);
        QObject::connect(pushButtonCancel, SIGNAL(clicked()), DialogReplaceLabel, SLOT(close()));
        QObject::connect(pushButtonReplace, SIGNAL(clicked()), DialogReplaceLabel, SLOT(OnReplace()));

        QMetaObject::connectSlotsByName(DialogReplaceLabel);
    } // setupUi

    void retranslateUi(QDialog *DialogReplaceLabel)
    {
        DialogReplaceLabel->setWindowTitle(QApplication::translate("DialogReplaceLabel", "Replace Label", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogReplaceLabel", "Original value", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("DialogReplaceLabel", "New value", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QApplication::translate("DialogReplaceLabel", "Options", 0, QApplication::UnicodeUTF8));
        radioButtonAllSlices->setText(QApplication::translate("DialogReplaceLabel", "Replace in all slices", 0, QApplication::UnicodeUTF8));
        radioButtonCurrentSlice->setText(QApplication::translate("DialogReplaceLabel", "Replace in current slice only", 0, QApplication::UnicodeUTF8));
        pushButtonCancel->setText(QApplication::translate("DialogReplaceLabel", "Cancel", 0, QApplication::UnicodeUTF8));
        pushButtonReplace->setText(QApplication::translate("DialogReplaceLabel", "Replace", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogReplaceLabel: public Ui_DialogReplaceLabel {};
} // namespace Ui

QT_END_NAMESPACE

