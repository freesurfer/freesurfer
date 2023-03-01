#ifndef MIGRATION_DEFS_H
#define MIGRATION_DEFS_H

#if (QT_VERSION < QT_VERSION_CHECK(5,15,0))
#define MD_SkipEmptyParts QString::SkipEmptyParts
#else
#define MD_SkipEmptyParts Qt::SkipEmptyParts
#endif

#endif // MIGRATIONDEFS_H
