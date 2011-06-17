#ifndef MAINAPPLICATION_H
#define MAINAPPLICATION_H

#include <QApplication>

class MainApplication : public QApplication
{
    Q_OBJECT
public:
    explicit MainApplication( int & argc, char ** argv );

signals:
    void GlobalProgress(int n);

public slots:
    void EmitProgress(int n)
    {
      emit GlobalProgress(n);
    }
};

#endif // MAINAPPLICATION_H
