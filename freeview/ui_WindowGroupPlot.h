/********************************************************************************
** Form generated from reading UI file 'WindowGroupPlot.ui'
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
#include <QtGui/QComboBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QSpacerItem>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "WidgetGroupPlot.h"

QT_BEGIN_NAMESPACE

class Ui_WindowGroupPlot
{
public:
    QVBoxLayout *verticalLayout;
    WidgetGroupPlot *widgetPlot;
    QHBoxLayout *horizontalLayout;
    QSpacerItem *horizontalSpacer;
    QLabel *label;
    QComboBox *comboBoxVariable;
    QSpacerItem *horizontalSpacer_3;
    QLabel *label_2;
    QComboBox *comboBoxPlotStyle;
    QSpacerItem *horizontalSpacer_2;

    void setupUi(QWidget *WindowGroupPlot)
    {
        if (WindowGroupPlot->objectName().isEmpty())
            WindowGroupPlot->setObjectName(QString::fromUtf8("WindowGroupPlot"));
        WindowGroupPlot->resize(505, 367);
        verticalLayout = new QVBoxLayout(WindowGroupPlot);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        widgetPlot = new WidgetGroupPlot(WindowGroupPlot);
        widgetPlot->setObjectName(QString::fromUtf8("widgetPlot"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(1);
        sizePolicy.setHeightForWidth(widgetPlot->sizePolicy().hasHeightForWidth());
        widgetPlot->setSizePolicy(sizePolicy);
        widgetPlot->setMinimumSize(QSize(400, 300));

        verticalLayout->addWidget(widgetPlot);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalSpacer = new QSpacerItem(5, 5, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);

        label = new QLabel(WindowGroupPlot);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout->addWidget(label);

        comboBoxVariable = new QComboBox(WindowGroupPlot);
        comboBoxVariable->setObjectName(QString::fromUtf8("comboBoxVariable"));

        horizontalLayout->addWidget(comboBoxVariable);

        horizontalSpacer_3 = new QSpacerItem(10, 5, QSizePolicy::Fixed, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_3);

        label_2 = new QLabel(WindowGroupPlot);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        horizontalLayout->addWidget(label_2);

        comboBoxPlotStyle = new QComboBox(WindowGroupPlot);
        comboBoxPlotStyle->setObjectName(QString::fromUtf8("comboBoxPlotStyle"));

        horizontalLayout->addWidget(comboBoxPlotStyle);

        horizontalSpacer_2 = new QSpacerItem(5, 5, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_2);


        verticalLayout->addLayout(horizontalLayout);


        retranslateUi(WindowGroupPlot);
        QObject::connect(comboBoxVariable, SIGNAL(currentIndexChanged(int)), widgetPlot, SLOT(SetCurrentVariableIndex(int)));
        QObject::connect(comboBoxPlotStyle, SIGNAL(currentIndexChanged(int)), widgetPlot, SLOT(SetPlotType(int)));

        QMetaObject::connectSlotsByName(WindowGroupPlot);
    } // setupUi

    void retranslateUi(QWidget *WindowGroupPlot)
    {
        WindowGroupPlot->setWindowTitle(QApplication::translate("WindowGroupPlot", "Group Plot", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("WindowGroupPlot", "Variable", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("WindowGroupPlot", "Style", 0, QApplication::UnicodeUTF8));
        comboBoxPlotStyle->clear();
        comboBoxPlotStyle->insertItems(0, QStringList()
         << QApplication::translate("WindowGroupPlot", "Point", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("WindowGroupPlot", "Histogram", 0, QApplication::UnicodeUTF8)
        );
    } // retranslateUi

};

namespace Ui {
    class WindowGroupPlot: public Ui_WindowGroupPlot {};
} // namespace Ui

QT_END_NAMESPACE

