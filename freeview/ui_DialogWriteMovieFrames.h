/********************************************************************************
** Form generated from reading UI file 'DialogWriteMovieFrames.ui'
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
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QToolButton>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogWriteMovieFrames
{
public:
    QVBoxLayout *verticalLayout;
    QLabel *label_5;
    QGridLayout *gridLayout;
    QLabel *label;
    QLineEdit *lineEditOutputLocation;
    QToolButton *toolButtonOpen;
    QLabel *label_2;
    QComboBox *comboBoxExtension;
    QSpacerItem *horizontalSpacer_2;
    QLabel *labelAngleStep;
    QDoubleSpinBox *doubleSpinBoxAngleStep;
    QLabel *labelSliceStartNumber;
    QSpinBox *spinBoxSliceStart;
    QLabel *labelSliceStep;
    QSpinBox *spinBoxSliceStep;
    QSpacerItem *verticalSpacer_2;
    QLabel *label_4;
    QSpacerItem *verticalSpacer;
    QHBoxLayout *horizontalLayout;
    QPushButton *pushButtonAbort;
    QSpacerItem *horizontalSpacer;
    QPushButton *pushButtonWrite;
    QPushButton *pushButtonClose;

    void setupUi(QDialog *DialogWriteMovieFrames)
    {
        if (DialogWriteMovieFrames->objectName().isEmpty())
            DialogWriteMovieFrames->setObjectName(QString::fromUtf8("DialogWriteMovieFrames"));
        DialogWriteMovieFrames->resize(426, 328);
        verticalLayout = new QVBoxLayout(DialogWriteMovieFrames);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        label_5 = new QLabel(DialogWriteMovieFrames);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        label_5->setWordWrap(true);

        verticalLayout->addWidget(label_5);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        label = new QLabel(DialogWriteMovieFrames);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        lineEditOutputLocation = new QLineEdit(DialogWriteMovieFrames);
        lineEditOutputLocation->setObjectName(QString::fromUtf8("lineEditOutputLocation"));
        lineEditOutputLocation->setMinimumSize(QSize(215, 0));

        gridLayout->addWidget(lineEditOutputLocation, 0, 1, 1, 2);

        toolButtonOpen = new QToolButton(DialogWriteMovieFrames);
        toolButtonOpen->setObjectName(QString::fromUtf8("toolButtonOpen"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/file_open_16.png"), QSize(), QIcon::Normal, QIcon::Off);
        toolButtonOpen->setIcon(icon);

        gridLayout->addWidget(toolButtonOpen, 0, 3, 1, 1);

        label_2 = new QLabel(DialogWriteMovieFrames);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 1, 0, 1, 1);

        comboBoxExtension = new QComboBox(DialogWriteMovieFrames);
        comboBoxExtension->setObjectName(QString::fromUtf8("comboBoxExtension"));
        comboBoxExtension->setMaximumSize(QSize(80, 16777215));

        gridLayout->addWidget(comboBoxExtension, 1, 1, 1, 1);

        horizontalSpacer_2 = new QSpacerItem(0, 0, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_2, 1, 2, 1, 1);

        labelAngleStep = new QLabel(DialogWriteMovieFrames);
        labelAngleStep->setObjectName(QString::fromUtf8("labelAngleStep"));

        gridLayout->addWidget(labelAngleStep, 2, 0, 1, 1);

        doubleSpinBoxAngleStep = new QDoubleSpinBox(DialogWriteMovieFrames);
        doubleSpinBoxAngleStep->setObjectName(QString::fromUtf8("doubleSpinBoxAngleStep"));
        doubleSpinBoxAngleStep->setDecimals(2);
        doubleSpinBoxAngleStep->setMinimum(-100);
        doubleSpinBoxAngleStep->setMaximum(100);
        doubleSpinBoxAngleStep->setValue(1);

        gridLayout->addWidget(doubleSpinBoxAngleStep, 2, 1, 1, 1);

        labelSliceStartNumber = new QLabel(DialogWriteMovieFrames);
        labelSliceStartNumber->setObjectName(QString::fromUtf8("labelSliceStartNumber"));

        gridLayout->addWidget(labelSliceStartNumber, 3, 0, 1, 1);

        spinBoxSliceStart = new QSpinBox(DialogWriteMovieFrames);
        spinBoxSliceStart->setObjectName(QString::fromUtf8("spinBoxSliceStart"));
        spinBoxSliceStart->setMinimum(0);
        spinBoxSliceStart->setMaximum(10000);
        spinBoxSliceStart->setValue(0);

        gridLayout->addWidget(spinBoxSliceStart, 3, 1, 1, 1);

        labelSliceStep = new QLabel(DialogWriteMovieFrames);
        labelSliceStep->setObjectName(QString::fromUtf8("labelSliceStep"));

        gridLayout->addWidget(labelSliceStep, 4, 0, 1, 1);

        spinBoxSliceStep = new QSpinBox(DialogWriteMovieFrames);
        spinBoxSliceStep->setObjectName(QString::fromUtf8("spinBoxSliceStep"));
        spinBoxSliceStep->setMinimum(1);
        spinBoxSliceStep->setMaximum(100);

        gridLayout->addWidget(spinBoxSliceStep, 4, 1, 1, 1);


        verticalLayout->addLayout(gridLayout);

        verticalSpacer_2 = new QSpacerItem(20, 5, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer_2);

        label_4 = new QLabel(DialogWriteMovieFrames);
        label_4->setObjectName(QString::fromUtf8("label_4"));
        label_4->setWordWrap(true);

        verticalLayout->addWidget(label_4);

        verticalSpacer = new QSpacerItem(10, 11, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        pushButtonAbort = new QPushButton(DialogWriteMovieFrames);
        pushButtonAbort->setObjectName(QString::fromUtf8("pushButtonAbort"));
        pushButtonAbort->setEnabled(false);

        horizontalLayout->addWidget(pushButtonAbort);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);

        pushButtonWrite = new QPushButton(DialogWriteMovieFrames);
        pushButtonWrite->setObjectName(QString::fromUtf8("pushButtonWrite"));
        pushButtonWrite->setEnabled(true);

        horizontalLayout->addWidget(pushButtonWrite);

        pushButtonClose = new QPushButton(DialogWriteMovieFrames);
        pushButtonClose->setObjectName(QString::fromUtf8("pushButtonClose"));

        horizontalLayout->addWidget(pushButtonClose);


        verticalLayout->addLayout(horizontalLayout);


        retranslateUi(DialogWriteMovieFrames);
        QObject::connect(pushButtonClose, SIGNAL(clicked()), DialogWriteMovieFrames, SLOT(close()));
        QObject::connect(toolButtonOpen, SIGNAL(clicked()), DialogWriteMovieFrames, SLOT(OnOpen()));
        QObject::connect(pushButtonAbort, SIGNAL(clicked()), DialogWriteMovieFrames, SLOT(OnAbort()));
        QObject::connect(pushButtonWrite, SIGNAL(clicked()), DialogWriteMovieFrames, SLOT(OnWrite()));
        QObject::connect(pushButtonClose, SIGNAL(clicked()), DialogWriteMovieFrames, SLOT(close()));

        QMetaObject::connectSlotsByName(DialogWriteMovieFrames);
    } // setupUi

    void retranslateUi(QDialog *DialogWriteMovieFrames)
    {
        DialogWriteMovieFrames->setWindowTitle(QApplication::translate("DialogWriteMovieFrames", "Write Movie Frames", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("DialogWriteMovieFrames", "Move this window out of the main viewport before starting to write.", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogWriteMovieFrames", "Output location", 0, QApplication::UnicodeUTF8));
        toolButtonOpen->setText(QString());
        label_2->setText(QApplication::translate("DialogWriteMovieFrames", "Filename extentsion", 0, QApplication::UnicodeUTF8));
        comboBoxExtension->clear();
        comboBoxExtension->insertItems(0, QStringList()
         << QApplication::translate("DialogWriteMovieFrames", "png", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogWriteMovieFrames", "jpeg", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogWriteMovieFrames", "tiff", 0, QApplication::UnicodeUTF8)
        );
        labelAngleStep->setText(QApplication::translate("DialogWriteMovieFrames", "Angle step size", 0, QApplication::UnicodeUTF8));
        labelSliceStartNumber->setText(QApplication::translate("DialogWriteMovieFrames", "Slice start nubmer", 0, QApplication::UnicodeUTF8));
        labelSliceStep->setText(QApplication::translate("DialogWriteMovieFrames", "Slice step size", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("DialogWriteMovieFrames", "NOTE: Angle step is for spinning in 3D view only. If the main viewport is 2D slice view, it will be ignored. A slice fly-through will be saved instead.", 0, QApplication::UnicodeUTF8));
        pushButtonAbort->setText(QApplication::translate("DialogWriteMovieFrames", "Abort", 0, QApplication::UnicodeUTF8));
        pushButtonWrite->setText(QApplication::translate("DialogWriteMovieFrames", "Write", 0, QApplication::UnicodeUTF8));
        pushButtonClose->setText(QApplication::translate("DialogWriteMovieFrames", "Close", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogWriteMovieFrames: public Ui_DialogWriteMovieFrames {};
} // namespace Ui

QT_END_NAMESPACE

