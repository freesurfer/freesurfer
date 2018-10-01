#ifndef DialogCustomFill_H
#define DialogCustomFill_H

#include <QDialog>
#include <QVariantMap>

namespace Ui {
class DialogCustomFill;
}

class DialogCustomFill : public QDialog
{
  Q_OBJECT

public:
  explicit DialogCustomFill(QWidget *parent = 0);
  ~DialogCustomFill();

signals:
  void CustomFillTriggered(const QVariantMap& options);

public slots:
  void OnButtonFill();

private:
  Ui::DialogCustomFill *ui;
};

#endif // DialogCustomFill_H
