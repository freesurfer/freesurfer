#ifndef GEODESICMATTING_H
#define GEODESICMATTING_H

#include <math.h>
#include <vector>
#include <QObject>
#include <QMutex>

class GeodesicMatting : public QObject
{
  Q_OBJECT
public:
  GeodesicMatting(QObject* parent = NULL );

//  bool Compute(int* dim, double* mri_in, double* mri_range_in,
//               unsigned char* seeds_in, std::vector<unsigned char>& label_list, unsigned char* seeds_out);
  bool ComputeWithBinning(int* dim, double* voxel_size, double* mri_in, double* mri_range_in,
                          unsigned char* seeds_in, std::vector<unsigned char>& label_list, unsigned char* seeds_out);
  void Abort();

  QString GetErrorMessage()
  {
    return m_strErrorMessage;
  }

signals:
  void Progress(double percentage);

private:
  double Interpolate(const std::vector<double>& v, const std::vector<double>& hf, double val);
  void Dilate(int* dim, unsigned char* in, unsigned char* out);
  double ComputeNeighDist(double** lHood, int nLabels, unsigned char* KNOWN, double* D, int* dim, int it, int jt, int kt, double* scale);
  int GetMinValIndex(const std::vector<double>& vals);
  int GetMinValIndexInSorted(const std::vector<double>& vals, const std::vector<std::size_t>& sorted_idx);

  bool m_bAbort;
  QString m_strErrorMessage;
  QMutex  m_mutex;
};

#endif // GEODESICMATTING_H
