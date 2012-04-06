/********************************************************************************
** Form generated from reading UI file 'PanelROI.ui'
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
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QScrollArea>
#include <QtGui/QSlider>
#include <QtGui/QSpacerItem>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "LayerTreeWidget.h"
#include "qtcolorpicker.h"

QT_BEGIN_NAMESPACE

class Ui_PanelROI
{
public:
    QAction *actionMoveLayerUp;
    QAction *actionMoveLayerDown;
    QWidget *scrollAreaWidgetContents;
    QVBoxLayout *verticalLayout_2;
    QVBoxLayout *verticalLayoutToolbar;
    QToolBar *toolbar;
    LayerTreeWidget *treeWidgetLayers;
    QGridLayout *gridLayout;
    QLabel *labelFileName;
    QLineEdit *lineEditFileName;
    QLabel *labelOpacity;
    QHBoxLayout *horizontalLayout_2;
    QSlider *sliderOpacity;
    QDoubleSpinBox *doubleSpinBoxOpacity;
    QHBoxLayout *horizontalLayout;
    QtColorPicker *colorPickerColor;
    QSpacerItem *horizontalSpacer;
    QLabel *labelColor;
    QSpacerItem *verticalSpacer;

    void setupUi(QScrollArea *PanelROI)
    {
        if (PanelROI->objectName().isEmpty())
            PanelROI->setObjectName(QString::fromUtf8("PanelROI"));
        PanelROI->resize(317, 901);
        PanelROI->setFrameShape(QFrame::NoFrame);
        PanelROI->setWidgetResizable(true);
        actionMoveLayerUp = new QAction(PanelROI);
        actionMoveLayerUp->setObjectName(QString::fromUtf8("actionMoveLayerUp"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/move_up.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionMoveLayerUp->setIcon(icon);
        actionMoveLayerDown = new QAction(PanelROI);
        actionMoveLayerDown->setObjectName(QString::fromUtf8("actionMoveLayerDown"));
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/resource/icons/move_down.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionMoveLayerDown->setIcon(icon1);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QString::fromUtf8("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 317, 901));
        verticalLayout_2 = new QVBoxLayout(scrollAreaWidgetContents);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(5, 0, 5, 5);
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


        verticalLayout_2->addLayout(verticalLayoutToolbar);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        labelFileName = new QLabel(scrollAreaWidgetContents);
        labelFileName->setObjectName(QString::fromUtf8("labelFileName"));

        gridLayout->addWidget(labelFileName, 1, 0, 1, 1);

        lineEditFileName = new QLineEdit(scrollAreaWidgetContents);
        lineEditFileName->setObjectName(QString::fromUtf8("lineEditFileName"));

        gridLayout->addWidget(lineEditFileName, 1, 1, 1, 1);

        labelOpacity = new QLabel(scrollAreaWidgetContents);
        labelOpacity->setObjectName(QString::fromUtf8("labelOpacity"));

        gridLayout->addWidget(labelOpacity, 2, 0, 1, 1);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        sliderOpacity = new QSlider(scrollAreaWidgetContents);
        sliderOpacity->setObjectName(QString::fromUtf8("sliderOpacity"));
        sliderOpacity->setMaximum(100);
        sliderOpacity->setValue(100);
        sliderOpacity->setOrientation(Qt::Horizontal);

        horizontalLayout_2->addWidget(sliderOpacity);

        doubleSpinBoxOpacity = new QDoubleSpinBox(scrollAreaWidgetContents);
        doubleSpinBoxOpacity->setObjectName(QString::fromUtf8("doubleSpinBoxOpacity"));
        doubleSpinBoxOpacity->setMaximum(1);
        doubleSpinBoxOpacity->setSingleStep(0.1);
        doubleSpinBoxOpacity->setValue(1);

        horizontalLayout_2->addWidget(doubleSpinBoxOpacity);


        gridLayout->addLayout(horizontalLayout_2, 2, 1, 1, 1);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        colorPickerColor = new QtColorPicker(scrollAreaWidgetContents);
        colorPickerColor->setObjectName(QString::fromUtf8("colorPickerColor"));

        horizontalLayout->addWidget(colorPickerColor);

        horizontalSpacer = new QSpacerItem(1, 1, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);


        gridLayout->addLayout(horizontalLayout, 3, 1, 1, 1);

        labelColor = new QLabel(scrollAreaWidgetContents);
        labelColor->setObjectName(QString::fromUtf8("labelColor"));

        gridLayout->addWidget(labelColor, 3, 0, 1, 1);


        verticalLayout_2->addLayout(gridLayout);

        verticalSpacer = new QSpacerItem(0, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_2->addItem(verticalSpacer);

        PanelROI->setWidget(scrollAreaWidgetContents);

        toolbar->addAction(actionMoveLayerUp);
        toolbar->addAction(actionMoveLayerDown);

        retranslateUi(PanelROI);
        QObject::connect(sliderOpacity, SIGNAL(valueChanged(int)), PanelROI, SLOT(OnSliderOpacity(int)));

        QMetaObject::connectSlotsByName(PanelROI);
    } // setupUi

    void retranslateUi(QScrollArea *PanelROI)
    {
        PanelROI->setWindowTitle(QApplication::translate("PanelROI", "ScrollArea", 0, QApplication::UnicodeUTF8));
        actionMoveLayerUp->setText(QApplication::translate("PanelROI", "Move Layer Up", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionMoveLayerUp->setToolTip(QApplication::translate("PanelROI", "Up", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionMoveLayerDown->setText(QApplication::translate("PanelROI", "Move Layer Down", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionMoveLayerDown->setToolTip(QApplication::translate("PanelROI", "Down", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        labelFileName->setText(QApplication::translate("PanelROI", "File name", 0, QApplication::UnicodeUTF8));
        labelOpacity->setText(QApplication::translate("PanelROI", "Opacity", 0, QApplication::UnicodeUTF8));
        labelColor->setText(QApplication::translate("PanelROI", "Color", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class PanelROI: public Ui_PanelROI {};
} // namespace Ui

QT_END_NAMESPACE

