#ifndef VOLUMEFILTEROPTIMAL_H
#define VOLUMEFILTEROPTIMAL_H

#include "VolumeFilter.h"
#include <QList>

class LayerMRI;
class LayerROI;

class VolumeFilterOptimal : public VolumeFilter
{
public:
  VolumeFilterOptimal( QList<LayerMRI*> inputs, QList<LayerROI*> input_labels, LayerMRI* output, QObject* parent = 0 );

  QString GetName()
  {
    return "Optimal";
  }

protected:
  bool Execute();

  QList<LayerMRI*>  m_inputMRIs;
  QList<LayerROI*>  m_inputROIs;
};

#endif // VOLUMEFILTEROPTIMAL_H
