#ifndef GEODESICMATTING_H
#define GEODESICMATTING_H

#include <math.h>
#include <vector>

#include <QMutex>

class GeodesicMatting
{
public:
  GeodesicMatting();

  bool Compute(int* dim, double* mri_in, double* mri_range_in,
               unsigned char* seeds_in, std::vector<unsigned char>& label_list, unsigned char* seeds_out);
  bool ComputeWithBinning(int* dim, double* mri_in, double* mri_range_in,
                          unsigned char* seeds_in, std::vector<unsigned char>& label_list, unsigned char* seeds_out);
  void Abort();

private:
  double Interpolate(const std::vector<double>& v, const std::vector<double>& hf, double val);
  void Dilate(int* dim, unsigned char* in, unsigned char* out);
  double ComputeNeighDist(double** lHood, int nLabels, unsigned char* KNOWN, double* D, int* dim, int it, int jt, int kt);
  int GetMinValIndex(const std::vector<double>& vals);
  int GetMinValIndexInSorted(const std::vector<double>& vals, const std::vector<std::size_t>& sorted_idx);

  bool m_bAbort;
  QMutex  m_mutex;
};

#endif // GEODESICMATTING_H
