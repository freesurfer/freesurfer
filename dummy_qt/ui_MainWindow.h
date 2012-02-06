/********************************************************************************
** Form generated from reading UI file 'MainWindow.ui'
**
** Created: Thu Jan 26 20:00:21 2012
**      by: Qt User Interface Compiler version 4.7.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QTextEdit>
#include <QtGui/QToolButton>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralWidget;
    QVBoxLayout *verticalLayout;
    QGridLayout *gridLayout;
    QLabel *label;
    QLineEdit *lineEditInput;
    QToolButton *toolButtonInput;
    QLabel *label_2;
    QLineEdit *lineEditOutput;
    QToolButton *toolButtonOutput;
    QSpacerItem *verticalSpacer;
    QHBoxLayout *horizontalLayout;
    QPushButton *pushButtonRun;
    QSpacerItem *horizontalSpacer;
    QPushButton *pushButtonAbort;
    QTextEdit *textEditLog;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(347, 312);
        centralWidget = new QWidget(MainWindow);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        verticalLayout = new QVBoxLayout(centralWidget);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        gridLayout = new QGridLayout();
        gridLayout->setSpacing(6);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        label = new QLabel(centralWidget);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 0, 0, 1, 1);

        lineEditInput = new QLineEdit(centralWidget);
        lineEditInput->setObjectName(QString::fromUtf8("lineEditInput"));

        gridLayout->addWidget(lineEditInput, 0, 1, 1, 1);

        toolButtonInput = new QToolButton(centralWidget);
        toolButtonInput->setObjectName(QString::fromUtf8("toolButtonInput"));
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/images/file_open_16.png"), QSize(), QIcon::Normal, QIcon::Off);
        toolButtonInput->setIcon(icon);

        gridLayout->addWidget(toolButtonInput, 0, 2, 1, 1);

        label_2 = new QLabel(centralWidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 1, 0, 1, 1);

        lineEditOutput = new QLineEdit(centralWidget);
        lineEditOutput->setObjectName(QString::fromUtf8("lineEditOutput"));

        gridLayout->addWidget(lineEditOutput, 1, 1, 1, 1);

        toolButtonOutput = new QToolButton(centralWidget);
        toolButtonOutput->setObjectName(QString::fromUtf8("toolButtonOutput"));
        toolButtonOutput->setIcon(icon);

        gridLayout->addWidget(toolButtonOutput, 1, 2, 1, 1);


        verticalLayout->addLayout(gridLayout);

        verticalSpacer = new QSpacerItem(20, 10, QSizePolicy::Minimum, QSizePolicy::Fixed);

        verticalLayout->addItem(verticalSpacer);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        pushButtonRun = new QPushButton(centralWidget);
        pushButtonRun->setObjectName(QString::fromUtf8("pushButtonRun"));
        pushButtonRun->setEnabled(false);

        horizontalLayout->addWidget(pushButtonRun);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout->addItem(horizontalSpacer);

        pushButtonAbort = new QPushButton(centralWidget);
        pushButtonAbort->setObjectName(QString::fromUtf8("pushButtonAbort"));
        pushButtonAbort->setEnabled(false);

        horizontalLayout->addWidget(pushButtonAbort);


        verticalLayout->addLayout(horizontalLayout);

        textEditLog = new QTextEdit(centralWidget);
        textEditLog->setObjectName(QString::fromUtf8("textEditLog"));
        textEditLog->setReadOnly(true);

        verticalLayout->addWidget(textEditLog);

        MainWindow->setCentralWidget(centralWidget);

        retranslateUi(MainWindow);
        QObject::connect(toolButtonInput, SIGNAL(clicked()), MainWindow, SLOT(OnButtonOpenInput()));
        QObject::connect(toolButtonOutput, SIGNAL(clicked()), MainWindow, SLOT(OnButtonOpenOutput()));
        QObject::connect(pushButtonRun, SIGNAL(clicked()), MainWindow, SLOT(OnButtonRun()));
        QObject::connect(pushButtonAbort, SIGNAL(clicked()), MainWindow, SLOT(OnButtonAbort()));
        QObject::connect(lineEditInput, SIGNAL(textChanged(QString)), MainWindow, SLOT(UpdateStatus()));
        QObject::connect(lineEditOutput, SIGNAL(textChanged(QString)), MainWindow, SLOT(UpdateStatus()));

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "mri_convert", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "Input", 0, QApplication::UnicodeUTF8));
        toolButtonInput->setText(QApplication::translate("MainWindow", "Open", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("MainWindow", "Output", 0, QApplication::UnicodeUTF8));
        toolButtonOutput->setText(QApplication::translate("MainWindow", "Open", 0, QApplication::UnicodeUTF8));
        pushButtonRun->setText(QApplication::translate("MainWindow", "Run", 0, QApplication::UnicodeUTF8));
        pushButtonAbort->setText(QApplication::translate("MainWindow", "Abort", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

