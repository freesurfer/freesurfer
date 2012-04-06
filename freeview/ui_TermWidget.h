/********************************************************************************
** Form generated from reading UI file 'TermWidget.ui'
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
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QTextEdit>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "CommandEdit.h"

QT_BEGIN_NAMESPACE

class Ui_TermWidget
{
public:
    QVBoxLayout *verticalLayout;
    QTextEdit *textLog;
    CommandEdit *textCommand;
    QLabel *label;

    void setupUi(QWidget *TermWidget)
    {
        if (TermWidget->objectName().isEmpty())
            TermWidget->setObjectName(QString::fromUtf8("TermWidget"));
        TermWidget->resize(615, 254);
        TermWidget->setWindowOpacity(1);
        verticalLayout = new QVBoxLayout(TermWidget);
        verticalLayout->setSpacing(0);
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        textLog = new QTextEdit(TermWidget);
        textLog->setObjectName(QString::fromUtf8("textLog"));
        QFont font;
        font.setFamily(QString::fromUtf8("Menlo"));
        font.setStyleStrategy(QFont::PreferDefault);
        textLog->setFont(font);
        textLog->setStyleSheet(QString::fromUtf8(""));
        textLog->setFrameShape(QFrame::NoFrame);
        textLog->setReadOnly(true);

        verticalLayout->addWidget(textLog);

        textCommand = new CommandEdit(TermWidget);
        textCommand->setObjectName(QString::fromUtf8("textCommand"));
        textCommand->setMaximumSize(QSize(16777215, 70));
        textCommand->setFont(font);
        textCommand->setFocusPolicy(Qt::StrongFocus);
        textCommand->setStyleSheet(QString::fromUtf8(""));
        textCommand->setFrameShape(QFrame::StyledPanel);
        textCommand->setFrameShadow(QFrame::Raised);
        textCommand->setReadOnly(false);

        verticalLayout->addWidget(textCommand);

        label = new QLabel(TermWidget);
        label->setObjectName(QString::fromUtf8("label"));
        label->setStyleSheet(QString::fromUtf8("color:#444;font-size:11px;"));
        label->setWordWrap(true);

        verticalLayout->addWidget(label);


        retranslateUi(TermWidget);

        QMetaObject::connectSlotsByName(TermWidget);
    } // setupUi

    void retranslateUi(QWidget *TermWidget)
    {
        TermWidget->setWindowTitle(QApplication::translate("TermWidget", "Command Console", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("TermWidget", "This is NOT a shell. Only FreeView commands are supported. Type \"-h\" for all the available commands.", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class TermWidget: public Ui_TermWidget {};
} // namespace Ui

QT_END_NAMESPACE

