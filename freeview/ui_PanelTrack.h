/********************************************************************************
** Form generated from reading UI file 'PanelTrack.ui'
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
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QScrollArea>
#include <QtGui/QSpacerItem>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "LayerTreeWidget.h"
#include "qtcolorpicker.h"

QT_BEGIN_NAMESPACE

class Ui_PanelTrack
{
public:
    QWidget *scrollAreaWidgetContents;
    QVBoxLayout *verticalLayout;
    QVBoxLayout *verticalLayoutToolbar;
    QToolBar *toolbar;
    LayerTreeWidget *treeWidgetLayers;
    QGridLayout *gridLayout;
    QLabel *labelFileName;
    QLineEdit *lineEditFileName;
    QLabel *label;
    QComboBox *comboBoxColorCode;
    QLabel *labelSolidColor;
    QHBoxLayout *horizontalLayout;
    QtColorPicker *colorPickerSolidColor;
    QSpacerItem *horizontalSpacer;
    QLabel *labelDirectionScheme;
    QComboBox *comboBoxDirectionScheme;
    QLabel *labelDirectionMapping;
    QComboBox *comboBoxDirectionMapping;
    QLabel *label_2;
    QComboBox *comboBoxRenderRep;
    QSpacerItem *verticalSpacer;

    void setupUi(QScrollArea *PanelTrack)
    {
        if (PanelTrack->objectName().isEmpty())
            PanelTrack->setObjectName(QString::fromUtf8("PanelTrack"));
        PanelTrack->resize(374, 449);
        PanelTrack->setFrameShape(QFrame::NoFrame);
        PanelTrack->setWidgetResizable(true);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QString::fromUtf8("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 374, 449));
        verticalLayout = new QVBoxLayout(scrollAreaWidgetContents);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(5, 0, 5, 5);
        verticalLayoutToolbar = new QVBoxLayout();
        verticalLayoutToolbar->setSpacing(0);
        verticalLayoutToolbar->setObjectName(QString::fromUtf8("verticalLayoutToolbar"));
        verticalLayoutToolbar->setContentsMargins(0, -1, -1, 0);
        toolbar = new QToolBar(scrollAreaWidgetContents);
        toolbar->setObjectName(QString::fromUtf8("toolbar"));
        toolbar->setIconSize(QSize(24, 24));
        toolbar->setFloatable(false);

        verticalLayoutToolbar->addWidget(toolbar);

        treeWidgetLayers = new LayerTreeWidget(scrollAreaWidgetContents);
        QTreeWidgetItem *__qtreewidgetitem = new QTreeWidgetItem();
        __qtreewidgetitem->setText(0, QString::fromUtf8("1"));
        treeWidgetLayers->setHeaderItem(__qtreewidgetitem);
        treeWidgetLayers->setObjectName(QString::fromUtf8("treeWidgetLayers"));
        QSizePolicy sizePolicy(QSizePolicy::Expanding, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(treeWidgetLayers->sizePolicy().hasHeightForWidth());
        treeWidgetLayers->setSizePolicy(sizePolicy);
        treeWidgetLayers->setMaximumSize(QSize(16777215, 145));
        treeWidgetLayers->setRootIsDecorated(false);
        treeWidgetLayers->setUniformRowHeights(true);
        treeWidgetLayers->header()->setVisible(false);

        verticalLayoutToolbar->addWidget(treeWidgetLayers);


        verticalLayout->addLayout(verticalLayoutToolbar);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        labelFileName = new QLabel(scrollAreaWidgetContents);
        labelFileName->setObjectName(QString::fromUtf8("labelFileName"));

        gridLayout->addWidget(labelFileName, 1, 0, 1, 1);

        lineEditFileName = new QLineEdit(scrollAreaWidgetContents);
        lineEditFileName->setObjectName(QString::fromUtf8("lineEditFileName"));
        lineEditFileName->setReadOnly(true);

        gridLayout->addWidget(lineEditFileName, 1, 1, 1, 1);

        label = new QLabel(scrollAreaWidgetContents);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 3, 0, 1, 1);

        comboBoxColorCode = new QComboBox(scrollAreaWidgetContents);
        comboBoxColorCode->setObjectName(QString::fromUtf8("comboBoxColorCode"));

        gridLayout->addWidget(comboBoxColorCode, 3, 1, 1, 1);

        labelSolidColor = new QLabel(scrollAreaWidgetContents);
        labelSolidColor->setObjectName(QString::fromUtf8("labelSolidColor"));
        labelSolidColor->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelSolidColor, 4, 0, 1, 1);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        colorPickerSolidColor = new QtColorPicker(scrollAreaWidgetContents);
        colorPickerSolidColor->setObjectName(QString::fromUtf8("colorPickerSolidColor"));

        horizontalLayout->addWidget(colorPickerSolidColor);

        horizontalSpacer = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);


        gridLayout->addLayout(horizontalLayout, 4, 1, 1, 1);

        labelDirectionScheme = new QLabel(scrollAreaWidgetContents);
        labelDirectionScheme->setObjectName(QString::fromUtf8("labelDirectionScheme"));
        labelDirectionScheme->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelDirectionScheme, 5, 0, 1, 1);

        comboBoxDirectionScheme = new QComboBox(scrollAreaWidgetContents);
        comboBoxDirectionScheme->setObjectName(QString::fromUtf8("comboBoxDirectionScheme"));

        gridLayout->addWidget(comboBoxDirectionScheme, 5, 1, 1, 1);

        labelDirectionMapping = new QLabel(scrollAreaWidgetContents);
        labelDirectionMapping->setObjectName(QString::fromUtf8("labelDirectionMapping"));
        labelDirectionMapping->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(labelDirectionMapping, 6, 0, 1, 1);

        comboBoxDirectionMapping = new QComboBox(scrollAreaWidgetContents);
        comboBoxDirectionMapping->setObjectName(QString::fromUtf8("comboBoxDirectionMapping"));

        gridLayout->addWidget(comboBoxDirectionMapping, 6, 1, 1, 1);

        label_2 = new QLabel(scrollAreaWidgetContents);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 2, 0, 1, 1);

        comboBoxRenderRep = new QComboBox(scrollAreaWidgetContents);
        comboBoxRenderRep->setObjectName(QString::fromUtf8("comboBoxRenderRep"));

        gridLayout->addWidget(comboBoxRenderRep, 2, 1, 1, 1);


        verticalLayout->addLayout(gridLayout);

        verticalSpacer = new QSpacerItem(0, 240, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        PanelTrack->setWidget(scrollAreaWidgetContents);

        retranslateUi(PanelTrack);

        QMetaObject::connectSlotsByName(PanelTrack);
    } // setupUi

    void retranslateUi(QScrollArea *PanelTrack)
    {
        PanelTrack->setWindowTitle(QApplication::translate("PanelTrack", "ScrollArea", 0, QApplication::UnicodeUTF8));
        labelFileName->setText(QApplication::translate("PanelTrack", "File name", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("PanelTrack", "Color code", 0, QApplication::UnicodeUTF8));
        comboBoxColorCode->clear();
        comboBoxColorCode->insertItems(0, QStringList()
         << QApplication::translate("PanelTrack", "Directional", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelTrack", "Solid Color", 0, QApplication::UnicodeUTF8)
        );
        labelSolidColor->setText(QApplication::translate("PanelTrack", "Color", 0, QApplication::UnicodeUTF8));
        labelDirectionScheme->setText(QApplication::translate("PanelTrack", "Scheme", 0, QApplication::UnicodeUTF8));
        comboBoxDirectionScheme->clear();
        comboBoxDirectionScheme->insertItems(0, QStringList()
         << QApplication::translate("PanelTrack", "End Points", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelTrack", "Mid Segment", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelTrack", "Every Segment", 0, QApplication::UnicodeUTF8)
        );
        labelDirectionMapping->setText(QApplication::translate("PanelTrack", "Mapping", 0, QApplication::UnicodeUTF8));
        comboBoxDirectionMapping->clear();
        comboBoxDirectionMapping->insertItems(0, QStringList()
         << QApplication::translate("PanelTrack", "RAS -> RGB", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelTrack", "RAS -> RBG", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelTrack", "RAS -> GRB", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelTrack", "RAS -> GBR", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelTrack", "RAS -> BRG", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelTrack", "RAS -> BGR", 0, QApplication::UnicodeUTF8)
        );
        label_2->setText(QApplication::translate("PanelTrack", "Render", 0, QApplication::UnicodeUTF8));
        comboBoxRenderRep->clear();
        comboBoxRenderRep->insertItems(0, QStringList()
         << QApplication::translate("PanelTrack", "Lines", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("PanelTrack", "Tubes (Slow!)", 0, QApplication::UnicodeUTF8)
        );
    } // retranslateUi

};

namespace Ui {
    class PanelTrack: public Ui_PanelTrack {};
} // namespace Ui

QT_END_NAMESPACE

