#ifndef JSON_H
#define JSON_H

#include <QString>
#include <QVariantMap>
#include <QtScript/QScriptEngine>

class Json
{
public:
  Json();

  QString encode(const QVariantMap &map);

  QVariantMap decode(const QString &jsonStr);

  QScriptValue encodeInner(const QVariantMap &map, QScriptEngine* engine);

  QVariantMap decodeInner(QScriptValue object);

  QList<QVariant> decodeInnerToList(QScriptValue arrayValue);
};

#endif // JSON_H
