#include "ProgressCallback.h"
#include "MainApplication.h"
#include <QDebug>

void ProgressCallback(int nProgress)
{
  static int n = 0;
  if (n != nProgress)
  {
    n = nProgress;
    qobject_cast<MainApplication*>(qApp)->EmitProgress(nProgress);
  }
}

