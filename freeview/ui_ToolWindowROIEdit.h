/********************************************************************************
** Form generated from reading UI file 'ToolWindowROIEdit.ui'
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
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QToolBar>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ToolWindowROIEdit
{
public:
    QAction *actionFreeHand;
    QAction *actionPolyLine;
    QAction *actionLiveWire;
    QAction *actionFill;
    QVBoxLayout *verticalLayout;
    QToolBar *toolbar;
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QSpinBox *spinBoxBrushSize;
    QSpacerItem *horizontalSpacer;

    void setupUi(QWidget *ToolWindowROIEdit)
    {
        if (ToolWindowROIEdit->objectName().isEmpty())
            ToolWindowROIEdit->setObjectName(QString::fromUtf8("ToolWindowROIEdit"));
        ToolWindowROIEdit->resize(190, 64);
        actionFreeHand = new QAction(ToolWindowROIEdit);
        actionFreeHand->setObjectName(QString::fromUtf8("actionFreeHand"));
        actionFreeHand->setCheckable(true);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/voxel_draw_freehand.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionFreeHand->setIcon(icon);
        actionPolyLine = new QAction(ToolWindowROIEdit);
        actionPolyLine->setObjectName(QString::fromUtf8("actionPolyLine"));
        actionPolyLine->setCheckable(true);
        QIcon icon1;
        icon1.addFile(QString::fromUtf8(":/resource/icons/voxel_draw_polyline.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionPolyLine->setIcon(icon1);
        actionLiveWire = new QAction(ToolWindowROIEdit);
        actionLiveWire->setObjectName(QString::fromUtf8("actionLiveWire"));
        actionLiveWire->setCheckable(true);
        QIcon icon2;
        icon2.addFile(QString::fromUtf8(":/resource/icons/voxel_draw_livewire.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionLiveWire->setIcon(icon2);
        actionFill = new QAction(ToolWindowROIEdit);
        actionFill->setObjectName(QString::fromUtf8("actionFill"));
        actionFill->setCheckable(true);
        QIcon icon3;
        icon3.addFile(QString::fromUtf8(":/resource/icons/voxel_draw_fill.png"), QSize(), QIcon::Normal, QIcon::Off);
        actionFill->setIcon(icon3);
        verticalLayout = new QVBoxLayout(ToolWindowROIEdit);
        verticalLayout->setSpacing(0);
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        toolbar = new QToolBar(ToolWindowROIEdit);
        toolbar->setObjectName(QString::fromUtf8("toolbar"));
        toolbar->setIconSize(QSize(24, 24));
        toolbar->setFloatable(false);

        verticalLayout->addWidget(toolbar);

        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        label = new QLabel(ToolWindowROIEdit);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout->addWidget(label);

        spinBoxBrushSize = new QSpinBox(ToolWindowROIEdit);
        spinBoxBrushSize->setObjectName(QString::fromUtf8("spinBoxBrushSize"));
        spinBoxBrushSize->setMinimum(1);
        spinBoxBrushSize->setMaximum(100);

        horizontalLayout->addWidget(spinBoxBrushSize);

        horizontalSpacer = new QSpacerItem(1, 1, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);


        verticalLayout_2->addLayout(horizontalLayout);


        verticalLayout->addLayout(verticalLayout_2);


        toolbar->addAction(actionFreeHand);
        toolbar->addAction(actionPolyLine);
        toolbar->addAction(actionLiveWire);
        toolbar->addAction(actionFill);

        retranslateUi(ToolWindowROIEdit);

        QMetaObject::connectSlotsByName(ToolWindowROIEdit);
    } // setupUi

    void retranslateUi(QWidget *ToolWindowROIEdit)
    {
        ToolWindowROIEdit->setWindowTitle(QApplication::translate("ToolWindowROIEdit", "ROI Edit", 0, QApplication::UnicodeUTF8));
        actionFreeHand->setText(QApplication::translate("ToolWindowROIEdit", "FreeHand", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionFreeHand->setToolTip(QApplication::translate("ToolWindowROIEdit", "Freehand", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionPolyLine->setText(QApplication::translate("ToolWindowROIEdit", "PolyLine", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionPolyLine->setToolTip(QApplication::translate("ToolWindowROIEdit", "Polyline", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionLiveWire->setText(QApplication::translate("ToolWindowROIEdit", "LiveWire", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionLiveWire->setToolTip(QApplication::translate("ToolWindowROIEdit", "Livewire", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        actionFill->setText(QApplication::translate("ToolWindowROIEdit", "Fill", 0, QApplication::UnicodeUTF8));
#ifndef QT_NO_TOOLTIP
        actionFill->setToolTip(QApplication::translate("ToolWindowROIEdit", "Fill", 0, QApplication::UnicodeUTF8));
#endif // QT_NO_TOOLTIP
        label->setText(QApplication::translate("ToolWindowROIEdit", "Brush Size:", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class ToolWindowROIEdit: public Ui_ToolWindowROIEdit {};
} // namespace Ui

QT_END_NAMESPACE

