/********************************************************************************
** Form generated from reading UI file 'DialogSmoothSurface.ui'
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
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogSmoothSurface
{
public:
    QVBoxLayout *verticalLayout;
    QGridLayout *gridLayout;
    QLabel *label;
    QLineEdit *lineEditIterations;
    QLabel *labelLambda;
    QLineEdit *lineEditLambda;
    QLabel *labelFrequencyCutoff;
    QLineEdit *lineEditFrequencyCutoff;
    QSpacerItem *horizontalSpacer_4;
    QLabel *labelHelper;
    QComboBox *comboBoxMethod;
    QLabel *label_5;
    QSpacerItem *horizontalSpacer_3;
    QSpacerItem *verticalSpacer;
    QHBoxLayout *horizontalLayout;
    QSpacerItem *horizontalSpacer;
    QPushButton *pushButton;
    QPushButton *pushButton_2;
    QSpacerItem *horizontalSpacer_2;

    void setupUi(QDialog *DialogSmoothSurface)
    {
        if (DialogSmoothSurface->objectName().isEmpty())
            DialogSmoothSurface->setObjectName(QString::fromUtf8("DialogSmoothSurface"));
        DialogSmoothSurface->resize(314, 206);
        verticalLayout = new QVBoxLayout(DialogSmoothSurface);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        label = new QLabel(DialogSmoothSurface);
        label->setObjectName(QString::fromUtf8("label"));
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label, 1, 1, 1, 1);

        lineEditIterations = new QLineEdit(DialogSmoothSurface);
        lineEditIterations->setObjectName(QString::fromUtf8("lineEditIterations"));
        lineEditIterations->setMaximumSize(QSize(167777, 16777215));

        gridLayout->addWidget(lineEditIterations, 1, 2, 1, 1);

        labelLambda = new QLabel(DialogSmoothSurface);
        labelLambda->setObjectName(QString::fromUtf8("labelLambda"));
        labelLambda->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelLambda, 2, 1, 1, 1);

        lineEditLambda = new QLineEdit(DialogSmoothSurface);
        lineEditLambda->setObjectName(QString::fromUtf8("lineEditLambda"));
        lineEditLambda->setMaximumSize(QSize(167777, 16777215));

        gridLayout->addWidget(lineEditLambda, 2, 2, 1, 1);

        labelFrequencyCutoff = new QLabel(DialogSmoothSurface);
        labelFrequencyCutoff->setObjectName(QString::fromUtf8("labelFrequencyCutoff"));
        labelFrequencyCutoff->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelFrequencyCutoff, 3, 1, 1, 1);

        lineEditFrequencyCutoff = new QLineEdit(DialogSmoothSurface);
        lineEditFrequencyCutoff->setObjectName(QString::fromUtf8("lineEditFrequencyCutoff"));
        lineEditFrequencyCutoff->setMaximumSize(QSize(167777, 16777215));

        gridLayout->addWidget(lineEditFrequencyCutoff, 3, 2, 1, 1);

        horizontalSpacer_4 = new QSpacerItem(5, 5, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_4, 1, 0, 1, 1);

        labelHelper = new QLabel(DialogSmoothSurface);
        labelHelper->setObjectName(QString::fromUtf8("labelHelper"));
        labelHelper->setStyleSheet(QString::fromUtf8("color: gray;"));

        gridLayout->addWidget(labelHelper, 2, 3, 1, 1);

        comboBoxMethod = new QComboBox(DialogSmoothSurface);
        comboBoxMethod->setObjectName(QString::fromUtf8("comboBoxMethod"));

        gridLayout->addWidget(comboBoxMethod, 0, 2, 1, 1);

        label_5 = new QLabel(DialogSmoothSurface);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_5, 0, 1, 1, 1);

        horizontalSpacer_3 = new QSpacerItem(5, 5, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_3, 2, 4, 1, 1);


        verticalLayout->addLayout(gridLayout);

        verticalSpacer = new QSpacerItem(20, 15, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalSpacer = new QSpacerItem(5, 5, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);

        pushButton = new QPushButton(DialogSmoothSurface);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));

        horizontalLayout->addWidget(pushButton);

        pushButton_2 = new QPushButton(DialogSmoothSurface);
        pushButton_2->setObjectName(QString::fromUtf8("pushButton_2"));

        horizontalLayout->addWidget(pushButton_2);

        horizontalSpacer_2 = new QSpacerItem(5, 5, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_2);


        verticalLayout->addLayout(horizontalLayout);


        retranslateUi(DialogSmoothSurface);
        QObject::connect(pushButton_2, SIGNAL(clicked()), DialogSmoothSurface, SLOT(close()));
        QObject::connect(pushButton, SIGNAL(clicked()), DialogSmoothSurface, SLOT(OnApply()));
        QObject::connect(comboBoxMethod, SIGNAL(currentIndexChanged(int)), DialogSmoothSurface, SLOT(OnMethod(int)));

        QMetaObject::connectSlotsByName(DialogSmoothSurface);
    } // setupUi

    void retranslateUi(QDialog *DialogSmoothSurface)
    {
        DialogSmoothSurface->setWindowTitle(QApplication::translate("DialogSmoothSurface", "Smooth Surface", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogSmoothSurface", "Number of iterations:", 0, QApplication::UnicodeUTF8));
        labelLambda->setText(QApplication::translate("DialogSmoothSurface", "Lambda:", 0, QApplication::UnicodeUTF8));
        lineEditLambda->setText(QApplication::translate("DialogSmoothSurface", "0.3", 0, QApplication::UnicodeUTF8));
        labelFrequencyCutoff->setText(QApplication::translate("DialogSmoothSurface", "Frequency cutoff:", 0, QApplication::UnicodeUTF8));
        lineEditFrequencyCutoff->setText(QApplication::translate("DialogSmoothSurface", "0.5", 0, QApplication::UnicodeUTF8));
        labelHelper->setText(QApplication::translate("DialogSmoothSurface", "0.0~1.0", 0, QApplication::UnicodeUTF8));
        comboBoxMethod->clear();
        comboBoxMethod->insertItems(0, QStringList()
         << QApplication::translate("DialogSmoothSurface", "Taubin", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogSmoothSurface", "Standard", 0, QApplication::UnicodeUTF8)
        );
        label_5->setText(QApplication::translate("DialogSmoothSurface", "Method:", 0, QApplication::UnicodeUTF8));
        pushButton->setText(QApplication::translate("DialogSmoothSurface", "Apply", 0, QApplication::UnicodeUTF8));
        pushButton_2->setText(QApplication::translate("DialogSmoothSurface", "Close", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogSmoothSurface: public Ui_DialogSmoothSurface {};
} // namespace Ui

QT_END_NAMESPACE

