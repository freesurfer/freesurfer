#ifndef CONFIGURATIONFILEFORM_H
#define CONFIGURATIONFILEFORM_H

#include <QFileInfo>
#include <QVector>
#include <QDir>
#include "QDialog"
#include "ui_configurationfileform.h"


bool validBval(QFileInfo file);
bool validBvec(QFileInfo file);
bool validTrk(QFileInfo file);
bool subjectSeeker(QString file_path, QString subject, QString dir_path);

class ConfigurationFileForm : public QDialog
{
    Q_OBJECT

public:
    explicit ConfigurationFileForm(QDialog *parent = 0);

    void loadDefaults(QString fileName);
    void loadDefaults2(QString fileName);

private slots:
    void on_pushButton_clicked(bool checked);
    void on_pushButton_3_clicked(bool checked);
    void on_pushButton_2_clicked(bool checked);
    void on_pushButton_4_clicked(bool checked);
    void on_pushButton_5_clicked(bool checked);
    void on_buttonBox_2_rejected();
    bool validDcm(QFileInfo file);
    void on_checkBox_26_stateChanged(int state);
    bool prelimWriteFile();
    void writeFile();
    QString bvecSeeker(QFileInfoList list);
    QString bvalSeeker(QFileInfoList list);
    bool replace_subject_dcm();
    bool replace_subject_bvec();
    bool replace_subject_bval();
    void matcher(QDir directory, QString pattern, QVector<QString>& vector);

private:
    Ui::Dialog ui;
    QString original_directory_path;
    QString original_parent_path;
    QString file_root;
    QString original_dcm;
    QString original_bvec;
    QString original_bval;
    QVector<QString> subjects;
    QVector<QString> dcms;
    QVector<QString> bvecs;
    QVector<QString> bvals;
    int length_of_root_directory;
};

#endif // CONFIGURATIONFILEFORM_H
