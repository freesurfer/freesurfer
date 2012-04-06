/********************************************************************************
** Form generated from reading UI file 'DialogRepositionSurface.ui'
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
#include <QtGui/QDialog>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QRadioButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QTabWidget>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_DialogRepositionSurface
{
public:
    QVBoxLayout *verticalLayout;
    QTabWidget *tabWidget;
    QWidget *tabVolume;
    QVBoxLayout *verticalLayout_2;
    QGridLayout *gridLayout;
    QLabel *label;
    QLineEdit *lineEditVertex;
    QComboBox *comboBoxTarget;
    QLineEdit *lineEditTarget;
    QLabel *label_2;
    QLineEdit *lineEditSize;
    QLabel *label_3;
    QLineEdit *lineEditSigma;
    QLabel *label_6;
    QLabel *label_7;
    QSpacerItem *verticalSpacer;
    QWidget *tabCoordinate;
    QGridLayout *gridLayout_2;
    QLabel *label_4;
    QHBoxLayout *horizontalLayout_5;
    QLineEdit *lineEditVertex2;
    QSpacerItem *horizontalSpacer_6;
    QLabel *label_5;
    QHBoxLayout *horizontalLayout_3;
    QLineEdit *lineEditCoordX;
    QLineEdit *lineEditCoordY;
    QLineEdit *lineEditCoordZ;
    QSpacerItem *horizontalSpacer_8;
    QHBoxLayout *horizontalLayout_4;
    QRadioButton *radioButtonCoordSurface;
    QRadioButton *radioButtonCoordRAS;
    QSpacerItem *horizontalSpacer_7;
    QSpacerItem *verticalSpacer_2;
    QHBoxLayout *horizontalLayout;
    QSpacerItem *horizontalSpacer_2;
    QPushButton *pushButtonApply;
    QPushButton *pushButtonSave;
    QPushButton *pushButtonSaveAs;
    QSpacerItem *horizontalSpacer_3;
    QHBoxLayout *horizontalLayout_2;
    QSpacerItem *horizontalSpacer_4;
    QPushButton *pushButtonUndo;
    QPushButton *pushButtonClose;
    QSpacerItem *horizontalSpacer_5;

    void setupUi(QDialog *DialogRepositionSurface)
    {
        if (DialogRepositionSurface->objectName().isEmpty())
            DialogRepositionSurface->setObjectName(QString::fromUtf8("DialogRepositionSurface"));
        DialogRepositionSurface->resize(328, 278);
        verticalLayout = new QVBoxLayout(DialogRepositionSurface);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        tabWidget = new QTabWidget(DialogRepositionSurface);
        tabWidget->setObjectName(QString::fromUtf8("tabWidget"));
        tabVolume = new QWidget();
        tabVolume->setObjectName(QString::fromUtf8("tabVolume"));
        verticalLayout_2 = new QVBoxLayout(tabVolume);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        label = new QLabel(tabVolume);
        label->setObjectName(QString::fromUtf8("label"));
        label->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label, 0, 0, 1, 1);

        lineEditVertex = new QLineEdit(tabVolume);
        lineEditVertex->setObjectName(QString::fromUtf8("lineEditVertex"));
        lineEditVertex->setMaximumSize(QSize(80, 16777215));

        gridLayout->addWidget(lineEditVertex, 0, 1, 1, 1);

        comboBoxTarget = new QComboBox(tabVolume);
        comboBoxTarget->setObjectName(QString::fromUtf8("comboBoxTarget"));

        gridLayout->addWidget(comboBoxTarget, 1, 0, 1, 1);

        lineEditTarget = new QLineEdit(tabVolume);
        lineEditTarget->setObjectName(QString::fromUtf8("lineEditTarget"));
        lineEditTarget->setMaximumSize(QSize(80, 16777215));

        gridLayout->addWidget(lineEditTarget, 1, 1, 1, 1);

        label_2 = new QLabel(tabVolume);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_2, 2, 0, 1, 1);

        lineEditSize = new QLineEdit(tabVolume);
        lineEditSize->setObjectName(QString::fromUtf8("lineEditSize"));
        lineEditSize->setMaximumSize(QSize(80, 16777215));

        gridLayout->addWidget(lineEditSize, 2, 1, 1, 1);

        label_3 = new QLabel(tabVolume);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        label_3->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout->addWidget(label_3, 3, 0, 1, 1);

        lineEditSigma = new QLineEdit(tabVolume);
        lineEditSigma->setObjectName(QString::fromUtf8("lineEditSigma"));
        lineEditSigma->setMaximumSize(QSize(80, 16777215));

        gridLayout->addWidget(lineEditSigma, 3, 1, 1, 1);

        label_6 = new QLabel(tabVolume);
        label_6->setObjectName(QString::fromUtf8("label_6"));
        label_6->setStyleSheet(QString::fromUtf8("color:gray;"));

        gridLayout->addWidget(label_6, 0, 2, 1, 1);

        label_7 = new QLabel(tabVolume);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        label_7->setStyleSheet(QString::fromUtf8("color: gray;"));

        gridLayout->addWidget(label_7, 1, 2, 1, 1);


        verticalLayout_2->addLayout(gridLayout);

        verticalSpacer = new QSpacerItem(5, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_2->addItem(verticalSpacer);

        tabWidget->addTab(tabVolume, QString());
        tabCoordinate = new QWidget();
        tabCoordinate->setObjectName(QString::fromUtf8("tabCoordinate"));
        gridLayout_2 = new QGridLayout(tabCoordinate);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        label_4 = new QLabel(tabCoordinate);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_4, 0, 0, 1, 1);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        lineEditVertex2 = new QLineEdit(tabCoordinate);
        lineEditVertex2->setObjectName(QString::fromUtf8("lineEditVertex2"));
        lineEditVertex2->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_5->addWidget(lineEditVertex2);

        horizontalSpacer_6 = new QSpacerItem(1, 1, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_5->addItem(horizontalSpacer_6);


        gridLayout_2->addLayout(horizontalLayout_5, 0, 1, 1, 1);

        label_5 = new QLabel(tabCoordinate);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setAlignment(Qt::AlignRight|Qt::AlignTrailing|Qt::AlignVCenter);

        gridLayout_2->addWidget(label_5, 1, 0, 1, 1);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        lineEditCoordX = new QLineEdit(tabCoordinate);
        lineEditCoordX->setObjectName(QString::fromUtf8("lineEditCoordX"));
        lineEditCoordX->setMinimumSize(QSize(60, 0));
        lineEditCoordX->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_3->addWidget(lineEditCoordX);

        lineEditCoordY = new QLineEdit(tabCoordinate);
        lineEditCoordY->setObjectName(QString::fromUtf8("lineEditCoordY"));
        lineEditCoordY->setMinimumSize(QSize(60, 0));
        lineEditCoordY->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_3->addWidget(lineEditCoordY);

        lineEditCoordZ = new QLineEdit(tabCoordinate);
        lineEditCoordZ->setObjectName(QString::fromUtf8("lineEditCoordZ"));
        lineEditCoordZ->setMinimumSize(QSize(60, 0));
        lineEditCoordZ->setMaximumSize(QSize(60, 16777215));

        horizontalLayout_3->addWidget(lineEditCoordZ);

        horizontalSpacer_8 = new QSpacerItem(1, 1, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_3->addItem(horizontalSpacer_8);


        gridLayout_2->addLayout(horizontalLayout_3, 1, 1, 1, 1);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        radioButtonCoordSurface = new QRadioButton(tabCoordinate);
        radioButtonCoordSurface->setObjectName(QString::fromUtf8("radioButtonCoordSurface"));
        radioButtonCoordSurface->setChecked(true);

        horizontalLayout_4->addWidget(radioButtonCoordSurface);

        radioButtonCoordRAS = new QRadioButton(tabCoordinate);
        radioButtonCoordRAS->setObjectName(QString::fromUtf8("radioButtonCoordRAS"));
        radioButtonCoordRAS->setChecked(false);

        horizontalLayout_4->addWidget(radioButtonCoordRAS);

        horizontalSpacer_7 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_4->addItem(horizontalSpacer_7);


        gridLayout_2->addLayout(horizontalLayout_4, 2, 1, 1, 1);

        verticalSpacer_2 = new QSpacerItem(5, 5, QSizePolicy::Minimum, QSizePolicy::Expanding);

        gridLayout_2->addItem(verticalSpacer_2, 3, 1, 1, 1);

        tabWidget->addTab(tabCoordinate, QString());

        verticalLayout->addWidget(tabWidget);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        horizontalSpacer_2 = new QSpacerItem(5, 13, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_2);

        pushButtonApply = new QPushButton(DialogRepositionSurface);
        pushButtonApply->setObjectName(QString::fromUtf8("pushButtonApply"));

        horizontalLayout->addWidget(pushButtonApply);

        pushButtonSave = new QPushButton(DialogRepositionSurface);
        pushButtonSave->setObjectName(QString::fromUtf8("pushButtonSave"));
        pushButtonSave->setEnabled(false);

        horizontalLayout->addWidget(pushButtonSave);

        pushButtonSaveAs = new QPushButton(DialogRepositionSurface);
        pushButtonSaveAs->setObjectName(QString::fromUtf8("pushButtonSaveAs"));

        horizontalLayout->addWidget(pushButtonSaveAs);

        horizontalSpacer_3 = new QSpacerItem(5, 5, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_3);


        verticalLayout->addLayout(horizontalLayout);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalSpacer_4 = new QSpacerItem(5, 5, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_4);

        pushButtonUndo = new QPushButton(DialogRepositionSurface);
        pushButtonUndo->setObjectName(QString::fromUtf8("pushButtonUndo"));
        pushButtonUndo->setEnabled(false);

        horizontalLayout_2->addWidget(pushButtonUndo);

        pushButtonClose = new QPushButton(DialogRepositionSurface);
        pushButtonClose->setObjectName(QString::fromUtf8("pushButtonClose"));

        horizontalLayout_2->addWidget(pushButtonClose);

        horizontalSpacer_5 = new QSpacerItem(5, 5, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_5);


        verticalLayout->addLayout(horizontalLayout_2);


        retranslateUi(DialogRepositionSurface);
        QObject::connect(pushButtonClose, SIGNAL(clicked()), DialogRepositionSurface, SLOT(hide()));
        QObject::connect(comboBoxTarget, SIGNAL(currentIndexChanged(int)), DialogRepositionSurface, SLOT(OnComboTarget(int)));
        QObject::connect(pushButtonApply, SIGNAL(clicked()), DialogRepositionSurface, SLOT(OnApply()));
        QObject::connect(pushButtonSave, SIGNAL(clicked()), DialogRepositionSurface, SLOT(OnSave()));
        QObject::connect(pushButtonSaveAs, SIGNAL(clicked()), DialogRepositionSurface, SLOT(OnSaveAs()));
        QObject::connect(pushButtonUndo, SIGNAL(clicked()), DialogRepositionSurface, SLOT(OnUndo()));
        QObject::connect(radioButtonCoordRAS, SIGNAL(toggled(bool)), DialogRepositionSurface, SLOT(OnCoordinateTypeChanged()));
        QObject::connect(radioButtonCoordSurface, SIGNAL(toggled(bool)), DialogRepositionSurface, SLOT(OnCoordinateTypeChanged()));

        tabWidget->setCurrentIndex(0);


        QMetaObject::connectSlotsByName(DialogRepositionSurface);
    } // setupUi

    void retranslateUi(QDialog *DialogRepositionSurface)
    {
        DialogRepositionSurface->setWindowTitle(QApplication::translate("DialogRepositionSurface", "Reposition Surface", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogRepositionSurface", "Vertex", 0, QApplication::UnicodeUTF8));
        comboBoxTarget->clear();
        comboBoxTarget->insertItems(0, QStringList()
         << QApplication::translate("DialogRepositionSurface", "Intensity", 0, QApplication::UnicodeUTF8)
        );
        label_2->setText(QApplication::translate("DialogRepositionSurface", "Size", 0, QApplication::UnicodeUTF8));
        lineEditSize->setText(QApplication::translate("DialogRepositionSurface", "1", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("DialogRepositionSurface", "Sigma", 0, QApplication::UnicodeUTF8));
        lineEditSigma->setText(QApplication::translate("DialogRepositionSurface", "2.0", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("DialogRepositionSurface", "Shift+Left Btn", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("DialogRepositionSurface", "Shift+Ctrl+Left Btn", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tabVolume), QApplication::translate("DialogRepositionSurface", "Volume Based", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("DialogRepositionSurface", "Vertex", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("DialogRepositionSurface", "Coordinate", 0, QApplication::UnicodeUTF8));
        radioButtonCoordSurface->setText(QApplication::translate("DialogRepositionSurface", "Surface", 0, QApplication::UnicodeUTF8));
        radioButtonCoordRAS->setText(QApplication::translate("DialogRepositionSurface", "RAS", 0, QApplication::UnicodeUTF8));
        tabWidget->setTabText(tabWidget->indexOf(tabCoordinate), QApplication::translate("DialogRepositionSurface", "Edit Coordinate", 0, QApplication::UnicodeUTF8));
        pushButtonApply->setText(QApplication::translate("DialogRepositionSurface", "Apply", 0, QApplication::UnicodeUTF8));
        pushButtonSave->setText(QApplication::translate("DialogRepositionSurface", "Save", 0, QApplication::UnicodeUTF8));
        pushButtonSaveAs->setText(QApplication::translate("DialogRepositionSurface", "Save As", 0, QApplication::UnicodeUTF8));
        pushButtonUndo->setText(QApplication::translate("DialogRepositionSurface", "Undo", 0, QApplication::UnicodeUTF8));
        pushButtonClose->setText(QApplication::translate("DialogRepositionSurface", "Close", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogRepositionSurface: public Ui_DialogRepositionSurface {};
} // namespace Ui

QT_END_NAMESPACE

