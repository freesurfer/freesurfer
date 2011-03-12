#include "WindowQuickReference.h"
#include "ui_WindowQuickReference.h"
#include <QFile>
#include <QSettings>

WindowQuickReference::WindowQuickReference(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::WindowQuickReference)
{
    ui->setupUi(this);
    setWindowFlags( Qt::Tool );
    QFile file(":/resource/QuickRef.html");
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    ui->textBrowser->setHtml(file.readAll());

    QSettings settings;
    restoreGeometry(settings.value("WindowQuickRef/Geometry").toByteArray());
}

WindowQuickReference::~WindowQuickReference()
{
    QSettings settings;
    settings.setValue("WindowQuickRef/Geometry", this->saveGeometry());
    delete ui;
}
