/********************************************************************************
** Form generated from reading UI file 'DialogLabelStats.ui'
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
#include <QtGui/QGridLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_DialogLabelStats
{
public:
    QVBoxLayout *verticalLayout;
    QGridLayout *gridLayout;
    QLabel *label;
    QLabel *labelLabel;
    QLabel *label_3;
    QLabel *labelCount;
    QLabel *label_5;
    QLabel *labelArea;
    QLabel *label_2;
    QLabel *labelMean;
    QSpacerItem *verticalSpacer;

    void setupUi(QWidget *DialogLabelStats)
    {
        if (DialogLabelStats->objectName().isEmpty())
            DialogLabelStats->setObjectName(QString::fromUtf8("DialogLabelStats"));
        DialogLabelStats->resize(200, 143);
        verticalLayout = new QVBoxLayout(DialogLabelStats);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        verticalLayout->setContentsMargins(20, 15, 20, 20);
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setHorizontalSpacing(8);
        label = new QLabel(DialogLabelStats);
        label->setObjectName(QString::fromUtf8("label"));
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label, 0, 0, 1, 1);

        labelLabel = new QLabel(DialogLabelStats);
        labelLabel->setObjectName(QString::fromUtf8("labelLabel"));
        labelLabel->setTextInteractionFlags(Qt::LinksAccessibleByMouse|Qt::TextSelectableByMouse);

        gridLayout->addWidget(labelLabel, 0, 1, 1, 1);

        label_3 = new QLabel(DialogLabelStats);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_3, 1, 0, 1, 1);

        labelCount = new QLabel(DialogLabelStats);
        labelCount->setObjectName(QString::fromUtf8("labelCount"));
        labelCount->setTextInteractionFlags(Qt::LinksAccessibleByMouse|Qt::TextSelectableByMouse);

        gridLayout->addWidget(labelCount, 1, 1, 1, 1);

        label_5 = new QLabel(DialogLabelStats);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_5, 2, 0, 1, 1);

        labelArea = new QLabel(DialogLabelStats);
        labelArea->setObjectName(QString::fromUtf8("labelArea"));
        labelArea->setTextInteractionFlags(Qt::LinksAccessibleByMouse|Qt::TextSelectableByMouse);

        gridLayout->addWidget(labelArea, 2, 1, 1, 1);

        label_2 = new QLabel(DialogLabelStats);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_2, 3, 0, 1, 1);

        labelMean = new QLabel(DialogLabelStats);
        labelMean->setObjectName(QString::fromUtf8("labelMean"));
        labelMean->setTextInteractionFlags(Qt::LinksAccessibleByMouse|Qt::TextSelectableByMouse);

        gridLayout->addWidget(labelMean, 3, 1, 1, 1);

        gridLayout->setColumnStretch(1, 1);

        verticalLayout->addLayout(gridLayout);

        verticalSpacer = new QSpacerItem(5, 5, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);


        retranslateUi(DialogLabelStats);

        QMetaObject::connectSlotsByName(DialogLabelStats);
    } // setupUi

    void retranslateUi(QWidget *DialogLabelStats)
    {
        DialogLabelStats->setWindowTitle(QApplication::translate("DialogLabelStats", "Label Stats", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogLabelStats", "Label:", 0, QApplication::UnicodeUTF8));
        labelLabel->setText(QApplication::translate("DialogLabelStats", "TextLabel", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("DialogLabelStats", "Count:", 0, QApplication::UnicodeUTF8));
        labelCount->setText(QApplication::translate("DialogLabelStats", "TextLabel", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("DialogLabelStats", "Area:", 0, QApplication::UnicodeUTF8));
        labelArea->setText(QApplication::translate("DialogLabelStats", "TextLabel", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("DialogLabelStats", "Mean:", 0, QApplication::UnicodeUTF8));
        labelMean->setText(QApplication::translate("DialogLabelStats", "TextLabel", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogLabelStats: public Ui_DialogLabelStats {};
} // namespace Ui

QT_END_NAMESPACE

