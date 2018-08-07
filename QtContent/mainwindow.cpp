#include <QtWidgets>
#include <QTextStream>
#include "mainwindow.h"
#include "configurationfileform.h"

MainWindow::MainWindow()
    : textEdit(new QPlainTextEdit)
{
    setCentralWidget(textEdit);

    createActions();
    createStatusBar();
    table = NULL;

}

MainWindow::~MainWindow()
{
    if (table != NULL)
    {
        delete table;
        table = NULL;
    }
}

void MainWindow::createStatusBar()
{
    statusBar()->showMessage(tr("Ready"));
}

void MainWindow::createActions()
{
    QMenu *fileMenu = menuBar()->addMenu(tr("&File"));
    QToolBar *fileToolBar = addToolBar(tr("File"));
    const QIcon newIcon = QIcon::fromTheme("document-new", QIcon(":/images/new.png"));
    QAction *newAct = new QAction(newIcon, tr("&New"), this);
    newAct->setShortcuts(QKeySequence::New);
    newAct->setStatusTip(tr("Create a new file"));
    connect(newAct, &QAction::triggered, this, &MainWindow::newFile);
    fileMenu->addAction(newAct);
    fileToolBar->addAction(newAct);

    const QIcon openIcon = QIcon::fromTheme("document-open", QIcon(":/images/open.png"));
    QAction *openAct = new QAction(openIcon, tr("&Open..."), this);
    openAct->setShortcuts(QKeySequence::Open);
    openAct->setStatusTip(tr("Open an existing file"));
    connect(openAct, &QAction::triggered, this, &MainWindow::open);
    fileMenu->addAction(openAct);
    fileToolBar->addAction(openAct);

}

void MainWindow::newFile()
{

}

void MainWindow::open()
{
    QString fileName = QFileDialog::getOpenFileName(this);
    // QTextStream(stdout) << fileName << endl;
    int length = fileName.length();

    if (length > 7)
    {
        QString end_of_file = fileName;
        end_of_file.remove(0, length - 7);

        // checks if we're looking at a "dcm" file... if so we call the function that will write the configuration file
        if (end_of_file == ".nii.gz")
        {
            table = new ConfigurationFileForm();
            table->loadDefaults(fileName);
            this->hide();
            table->show();
        }

        // if its not a "dcm" file... we just open it as we normally would
        else
        {
            loadFile(fileName);
        }
    }
}

void MainWindow::loadFile(const QString &fileName)
{
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(QDir::toNativeSeparators(fileName), file.errorString()));
        return;
    }

    QTextStream in(&file);

#ifndef QT_NO_CURSOR

    QApplication::setOverrideCursor(Qt::WaitCursor);

#endif

    textEdit->setPlainText(in.readAll());

#ifndef QT_NO_CURSOR

    QApplication::restoreOverrideCursor();
#endif

    setCurrentFile(fileName);
    statusBar()->showMessage(tr("File loaded"), 2000);
}

void MainWindow::setCurrentFile(const QString &fileName)
{
    curFile = fileName;
    textEdit->document()->setModified(false);
    setWindowModified(false);

    QString shownName = curFile;
    if (curFile.isEmpty())
        shownName = "untitled.txt";
    setWindowFilePath(shownName);
}

// this function checks if there is any .bvec file in the parent directory, if it is... they might be used for the subjects (if said subjects dont have any gradient tables or coordinates of their own)
QString bvecSeeker(QFileInfoList list, QFileInfo current)
{
    for (int i = 0; i < list.size(); ++i)
    {
        current = list[i];
        if (current.isFile())
        {
            if (current.fileName() == "gradients.txt")
            {
                return current.filePath();
            }
        }
    }

    return "";
}

// this function checks if there is any .bval file in the parent directory, if it is... they might be used for the subjects (if said subjects dont have any gradient tables or coordinates of their own)
QString bvalSeeker(QFileInfoList list, QFileInfo current)
{
    for (int i = 0; i < list.size(); ++i)
    {
        current = list[i];
        if (current.isFile())
        {
            if (current.fileName() == "bvalues.txt")
            {
                return current.filePath();
            }
        }
    }

    return "";
}
