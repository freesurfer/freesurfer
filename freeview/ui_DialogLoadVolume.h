/********************************************************************************
** Form generated from reading UI file 'DialogLoadVolume.ui'
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
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QDialog>
#include <QtGui/QDialogButtonBox>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QRadioButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QToolButton>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_DialogLoadVolume
{
public:
    QVBoxLayout *verticalLayout;
    QLabel *label;
    QHBoxLayout *horizontalLayout;
    QComboBox *comboBoxFilenames;
    QToolButton *toolButtonOpen;
    QCheckBox *checkBoxResampleToRAS;
    QCheckBox *checkBoxRegistration;
    QHBoxLayout *horizontalLayout_2;
    QLineEdit *lineEditRegistration;
    QToolButton *toolButtonOpenRegistration;
    QGroupBox *groupBox;
    QGridLayout *gridLayout;
    QRadioButton *radioNearest;
    QRadioButton *radioTrilinear;
    QSpacerItem *horizontalSpacer;
    QLabel *label_2;
    QRadioButton *radioCubic;
    QGridLayout *gridLayout_2;
    QLabel *label_3;
    QComboBox *comboBoxColorMap;
    QLabel *labelLUT;
    QComboBox *comboBoxLUT;
    QSpacerItem *verticalSpacer;
    QDialogButtonBox *buttonBox;

    void setupUi(QDialog *DialogLoadVolume)
    {
        if (DialogLoadVolume->objectName().isEmpty())
            DialogLoadVolume->setObjectName(QString::fromUtf8("DialogLoadVolume"));
        DialogLoadVolume->resize(325, 395);
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(DialogLoadVolume->sizePolicy().hasHeightForWidth());
        DialogLoadVolume->setSizePolicy(sizePolicy);
        DialogLoadVolume->setMinimumSize(QSize(0, 0));
        DialogLoadVolume->setMaximumSize(QSize(1000000, 1000000));
        DialogLoadVolume->setModal(true);
        verticalLayout = new QVBoxLayout(DialogLoadVolume);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setSizeConstraint(QLayout::SetFixedSize);
        label = new QLabel(DialogLoadVolume);
        label->setObjectName(QString::fromUtf8("label"));

        verticalLayout->addWidget(label);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        comboBoxFilenames = new QComboBox(DialogLoadVolume);
        comboBoxFilenames->setObjectName(QString::fromUtf8("comboBoxFilenames"));
        comboBoxFilenames->setMaximumSize(QSize(400, 16777215));
        comboBoxFilenames->setEditable(true);

        horizontalLayout->addWidget(comboBoxFilenames);

        toolButtonOpen = new QToolButton(DialogLoadVolume);
        toolButtonOpen->setObjectName(QString::fromUtf8("toolButtonOpen"));
        QSizePolicy sizePolicy1(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(toolButtonOpen->sizePolicy().hasHeightForWidth());
        toolButtonOpen->setSizePolicy(sizePolicy1);
        QIcon icon;
        icon.addFile(QString::fromUtf8(":/resource/icons/file_open_16.png"), QSize(), QIcon::Normal, QIcon::Off);
        toolButtonOpen->setIcon(icon);

        horizontalLayout->addWidget(toolButtonOpen);


        verticalLayout->addLayout(horizontalLayout);

        checkBoxResampleToRAS = new QCheckBox(DialogLoadVolume);
        checkBoxResampleToRAS->setObjectName(QString::fromUtf8("checkBoxResampleToRAS"));
        QSizePolicy sizePolicy2(QSizePolicy::Preferred, QSizePolicy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(checkBoxResampleToRAS->sizePolicy().hasHeightForWidth());
        checkBoxResampleToRAS->setSizePolicy(sizePolicy2);

        verticalLayout->addWidget(checkBoxResampleToRAS);

        checkBoxRegistration = new QCheckBox(DialogLoadVolume);
        checkBoxRegistration->setObjectName(QString::fromUtf8("checkBoxRegistration"));

        verticalLayout->addWidget(checkBoxRegistration);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        lineEditRegistration = new QLineEdit(DialogLoadVolume);
        lineEditRegistration->setObjectName(QString::fromUtf8("lineEditRegistration"));
        lineEditRegistration->setEnabled(false);
        lineEditRegistration->setMinimumSize(QSize(0, 0));

        horizontalLayout_2->addWidget(lineEditRegistration);

        toolButtonOpenRegistration = new QToolButton(DialogLoadVolume);
        toolButtonOpenRegistration->setObjectName(QString::fromUtf8("toolButtonOpenRegistration"));
        toolButtonOpenRegistration->setEnabled(false);
        toolButtonOpenRegistration->setIcon(icon);

        horizontalLayout_2->addWidget(toolButtonOpenRegistration);


        verticalLayout->addLayout(horizontalLayout_2);

        groupBox = new QGroupBox(DialogLoadVolume);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        groupBox->setEnabled(true);
        QFont font;
        font.setBold(false);
        font.setWeight(50);
        font.setKerning(false);
        groupBox->setFont(font);
        groupBox->setTitle(QString::fromUtf8("Sample method"));
        groupBox->setFlat(false);
        groupBox->setCheckable(false);
        gridLayout = new QGridLayout(groupBox);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        radioNearest = new QRadioButton(groupBox);
        radioNearest->setObjectName(QString::fromUtf8("radioNearest"));
        radioNearest->setChecked(true);

        gridLayout->addWidget(radioNearest, 0, 0, 1, 1);

        radioTrilinear = new QRadioButton(groupBox);
        radioTrilinear->setObjectName(QString::fromUtf8("radioTrilinear"));

        gridLayout->addWidget(radioTrilinear, 0, 1, 1, 1);

        horizontalSpacer = new QSpacerItem(115, 17, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer, 0, 3, 1, 1);

        label_2 = new QLabel(groupBox);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setTextFormat(Qt::PlainText);
        label_2->setAlignment(Qt::AlignLeading|Qt::AlignLeft|Qt::AlignTop);
        label_2->setWordWrap(true);

        gridLayout->addWidget(label_2, 1, 0, 1, 4);

        radioCubic = new QRadioButton(groupBox);
        radioCubic->setObjectName(QString::fromUtf8("radioCubic"));

        gridLayout->addWidget(radioCubic, 0, 2, 1, 1);


        verticalLayout->addWidget(groupBox);

        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        label_3 = new QLabel(DialogLoadVolume);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout_2->addWidget(label_3, 0, 0, 1, 1);

        comboBoxColorMap = new QComboBox(DialogLoadVolume);
        comboBoxColorMap->setObjectName(QString::fromUtf8("comboBoxColorMap"));

        gridLayout_2->addWidget(comboBoxColorMap, 0, 1, 1, 1);

        labelLUT = new QLabel(DialogLoadVolume);
        labelLUT->setObjectName(QString::fromUtf8("labelLUT"));
        labelLUT->setEnabled(false);

        gridLayout_2->addWidget(labelLUT, 1, 0, 1, 1);

        comboBoxLUT = new QComboBox(DialogLoadVolume);
        comboBoxLUT->setObjectName(QString::fromUtf8("comboBoxLUT"));
        comboBoxLUT->setEnabled(false);
        comboBoxLUT->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_2->addWidget(comboBoxLUT, 1, 1, 1, 1);

        gridLayout_2->setColumnStretch(1, 1);

        verticalLayout->addLayout(gridLayout_2);

        verticalSpacer = new QSpacerItem(0, 5, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);

        buttonBox = new QDialogButtonBox(DialogLoadVolume);
        buttonBox->setObjectName(QString::fromUtf8("buttonBox"));
        buttonBox->setOrientation(Qt::Horizontal);
        buttonBox->setStandardButtons(QDialogButtonBox::Cancel|QDialogButtonBox::Ok);

        verticalLayout->addWidget(buttonBox);


        retranslateUi(DialogLoadVolume);
        QObject::connect(buttonBox, SIGNAL(rejected()), DialogLoadVolume, SLOT(reject()));
        QObject::connect(checkBoxRegistration, SIGNAL(toggled(bool)), lineEditRegistration, SLOT(setEnabled(bool)));
        QObject::connect(checkBoxRegistration, SIGNAL(toggled(bool)), toolButtonOpenRegistration, SLOT(setEnabled(bool)));
        QObject::connect(toolButtonOpen, SIGNAL(clicked()), DialogLoadVolume, SLOT(OnOpen()));
        QObject::connect(toolButtonOpenRegistration, SIGNAL(clicked()), DialogLoadVolume, SLOT(OnOpenRegistration()));
        QObject::connect(comboBoxColorMap, SIGNAL(currentIndexChanged(int)), DialogLoadVolume, SLOT(OnColorMap(int)));
        QObject::connect(comboBoxLUT, SIGNAL(currentIndexChanged(int)), DialogLoadVolume, SLOT(OnLUT(int)));
        QObject::connect(buttonBox, SIGNAL(accepted()), DialogLoadVolume, SLOT(OnOK()));

        QMetaObject::connectSlotsByName(DialogLoadVolume);
    } // setupUi

    void retranslateUi(QDialog *DialogLoadVolume)
    {
        DialogLoadVolume->setWindowTitle(QApplication::translate("DialogLoadVolume", "Load Volume", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("DialogLoadVolume", "Select volume file", 0, QApplication::UnicodeUTF8));
        toolButtonOpen->setText(QApplication::translate("DialogLoadVolume", "...", 0, QApplication::UnicodeUTF8));
        checkBoxResampleToRAS->setText(QApplication::translate("DialogLoadVolume", "Resample to standard RAS space", 0, QApplication::UnicodeUTF8));
        checkBoxRegistration->setText(QApplication::translate("DialogLoadVolume", "Apply registration file", 0, QApplication::UnicodeUTF8));
        toolButtonOpenRegistration->setText(QApplication::translate("DialogLoadVolume", "...", 0, QApplication::UnicodeUTF8));
        radioNearest->setText(QApplication::translate("DialogLoadVolume", "Nearest neighbor", 0, QApplication::UnicodeUTF8));
        radioTrilinear->setText(QApplication::translate("DialogLoadVolume", "Trilinear", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("DialogLoadVolume", "ALWAYS choose 'Nearest Neighbor' if you are loading a label volume.", 0, QApplication::UnicodeUTF8));
        radioCubic->setText(QApplication::translate("DialogLoadVolume", "Cubic", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("DialogLoadVolume", "Color map", 0, QApplication::UnicodeUTF8));
        comboBoxColorMap->clear();
        comboBoxColorMap->insertItems(0, QStringList()
         << QApplication::translate("DialogLoadVolume", "Grayscale", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogLoadVolume", "Lookup Table", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogLoadVolume", "Heat", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogLoadVolume", "Jet", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogLoadVolume", "GE Color", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("DialogLoadVolume", "NIH", 0, QApplication::UnicodeUTF8)
        );
        labelLUT->setText(QApplication::translate("DialogLoadVolume", "Look up table", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class DialogLoadVolume: public Ui_DialogLoadVolume {};
} // namespace Ui

QT_END_NAMESPACE

