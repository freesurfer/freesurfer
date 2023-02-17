#ifndef DIALOGSAVEALLVOLUMES_H
#define DIALOGSAVEALLVOLUMES_H

#include <QDialog>
#include <QThread>

class LayerMRI;

namespace Ui {
class DialogSaveAllVolumes;
}

class ThreadSaveAll : public QThread
{
  Q_OBJECT
public:
  explicit ThreadSaveAll(QObject *parent = 0);

  void SaveAllVolumes(const QList<LayerMRI*>& layers);

signals:
  void Progress( int n );
  void Saved(LayerMRI* mri, const QString& error_msg = "");

protected:
  void run();

private:
  QList<LayerMRI*> m_listMRIs;
};

class DialogSaveAllVolumes : public QDialog
{
  Q_OBJECT

public:
  explicit DialogSaveAllVolumes(QWidget *parent = nullptr);
  ~DialogSaveAllVolumes();

public slots:
  void OnButtonOpen();
  void OnButtonSave();
  void OnVolumeSaved(LayerMRI* mri, const QString& error_msg = "");
  void OnFinished();

private:
  Ui::DialogSaveAllVolumes *ui;
  QString m_strLastDir;
  ThreadSaveAll*  m_thread;
};

#endif // DIALOGSAVEALLVOLUMES_H
