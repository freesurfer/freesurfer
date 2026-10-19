#ifndef COMMONDEF_H
#define COMMONDEF_H

#define PY_COMMAND "fspython"

// QFont::setFamilies() needs Qt 5.13; FreeSurfer still builds on Qt 5.12 for
// several platforms. On 5.13+ this calls setFamilies(); before that it picks
// the first installed family. Callers also set a style hint as a last resort.
#include <QFont>
#include <QStringList>
#if QT_VERSION < QT_VERSION_CHECK(5, 13, 0)
#include <QFontDatabase>
#endif

inline void FsSetFontFamilies(QFont& font, const QStringList& families)
{
#if QT_VERSION >= QT_VERSION_CHECK(5, 13, 0)
  font.setFamilies(families);
#else
  QFontDatabase db;
  for (int i = 0; i < families.size(); i++)
  {
    if (db.hasFamily(families[i]))
    {
      font.setFamily(families[i]);
      return;
    }
  }
#endif
}

#endif // COMMONDEF_H
