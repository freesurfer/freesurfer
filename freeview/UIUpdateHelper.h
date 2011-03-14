#ifndef UIUPDATEHELPER_H
#define UIUPDATEHELPER_H

#include <QList>

class QWidget;
class QSpinBox;
class QDoubleSpinBox;
class QLineEdit;

class UIUpdateHelper
{
public:
    UIUpdateHelper();

    void ShowWidgets( const QList<QWidget*>& widgets, bool bShow );
    void EnableWidgets( const QList<QWidget*>& widgets, bool bEnable );

    void ChangeLineEditText( QLineEdit* w, const QString& strg );
    void ChangeLineEditNumber( QLineEdit* w, double val );
    void ChangeSpinBoxValue( QSpinBox* w, int nVal );
    void ChangeDoubleSpinBoxValue( QDoubleSpinBox* w, double dVal );
    void BlockAllSignals(QWidget* w, bool bBlock);
};

#endif // UIUPDATEHELPER_H
