#include "GeodesicMatting.h"
#include "stdlib.h"
#include <cstring>
#include "kde.h"
#include <QtDebug>
#include <QElapsedTimer>

#define LARGENUMBER 1e12

using namespace std;

GeodesicMatting::GeodesicMatting()
{

}

double GeodesicMatting::Interpolate(const std::vector<double> &v, const std::vector<double> &hf, double val)
{
  if (val <= v[0])
    return hf[0];
  else if (val >= v[v.size()-1])
    return hf[hf.size()-1];
  else
  {
    double interval = v[1] - v[0];
    int n = (int)((val-v[0])/interval);
    return (val-v[n])/interval*hf[n+1] + (v[n+1]-val)/interval*hf[n];
  }
}

// dilate and return volume with only new voxels
void GeodesicMatting::Dilate(int *dim, unsigned char *p_in, unsigned char *p_out)
{
  size_t vol_size = ((size_t)dim[0])*dim[1]*dim[2];
  memset(p_out, 0, vol_size);
  for (size_t i = 0; i < dim[0]; i++)
  {
    for (size_t j = 0; j < dim[1]; j++)
    {
      for (size_t k = 0; k < dim[2]; k++)
      {
        size_t idx0 = k*dim[0]*dim[1] + j*dim[0] + i;
        if (p_in[idx0])
        {
          if (i > 0 && !p_in[idx0-1])
          {
            p_out[idx0-1] = 1;
          }
          if (i < dim[0]-1 && !p_in[idx0+1])
          {
            p_out[idx0+1] = 1;
          }
          if (j > 0 && !p_in[idx0-dim[0]])
          {
            p_out[idx0-dim[0]] = 1;
          }
          if (j < dim[1]-1 && !p_in[idx0+dim[0]])
          {
            p_out[idx0+dim[0]] = 1;
          }
          if (k > 0 && !p_in[idx0-dim[0]*dim[1]])
          {
            p_out[idx0-dim[0]*dim[1]] = 1;
          }
          if (k < dim[2]-1 && !p_in[idx0+dim[0]*dim[1]])
          {
            p_out[idx0+dim[0]*dim[1]] = 1;
          }
        }
      }
    }
  }
}

double GeodesicMatting::ComputeNeighDist(double **lHood, int nLabels, unsigned char *KNOWN, double *D, int *dim, int it, int jt, int kt)
{
  size_t idx = ((size_t)it) + jt*dim[0] + kt*dim[0]*dim[1];
  double probs1[10], d = 1e12;
  for (int i = 0; i < nLabels; i++)
    probs1[i] = lHood[i][idx];
  if (it > 0 && KNOWN[idx-1])
  {
    double sum = 0;
    for (int i = 0; i < nLabels; i++)
      sum += qAbs(lHood[i][idx-1] - probs1[i]);
    d = qMin(d, D[idx-1]+sum);
  }
  if (it < dim[0]-1 && KNOWN[idx+1])
  {
    double sum = 0;
    for (int i = 0; i < nLabels; i++)
      sum += qAbs(lHood[i][idx+1] - probs1[i]);
    d = qMin(d, D[idx+1]+sum);
  }
  if (jt > 0 && KNOWN[idx-dim[0]])
  {
    double sum = 0;
    for (int i = 0; i < nLabels; i++)
      sum += qAbs(lHood[i][idx-dim[0]] - probs1[i]);
    d = qMin(d, D[idx-dim[0]]+sum);
  }
  if (jt < dim[1]-1 && KNOWN[idx+dim[0]])
  {
    double sum = 0;
    for (int i = 0; i < nLabels; i++)
      sum += qAbs(lHood[i][idx+dim[0]] - probs1[i]);
    d = qMin(d, D[idx+dim[0]]+sum);
  }
  if (kt > 0 && KNOWN[idx-dim[0]*dim[1]])
  {
    double sum = 0;
    for (int i = 0; i < nLabels; i++)
      sum += qAbs(lHood[i][idx-dim[0]*dim[1]] - probs1[i]);
    d = qMin(d, D[idx-dim[0]*dim[1]]+sum);
  }
  if (kt < dim[2]-1 && KNOWN[idx+dim[0]*dim[1]])
  {
    double sum = 0;
    for (int i = 0; i < nLabels; i++)
      sum += qAbs(lHood[i][idx+dim[0]*dim[1]] - probs1[i]);
    d = qMin(d, D[idx+dim[0]*dim[1]]+sum);
  }
  return d;
}

int FindMinFromSorted(const std::vector<double>& vals, const std::vector<size_t>& sorted_idx, int nStart, int nEnd)
{
  if (nEnd < nStart) return 0;
  if (nEnd == nStart) return nStart;

  int nMid = nStart + (nEnd-nStart)/2;

  if (nMid < nEnd && vals[sorted_idx[nMid+1]] < vals[sorted_idx[nMid]])
    return nMid+1;

  if (nMid > nStart && vals[sorted_idx[nMid]] < vals[sorted_idx[nMid-1]])
    return nMid;

  if (vals[sorted_idx[nEnd]] > vals[sorted_idx[nMid]])
    return FindMinFromSorted(vals, sorted_idx, nStart, nMid-1);

  return FindMinFromSorted(vals, sorted_idx, nMid+1, nEnd);
}

int BinarySearch(double val, std::vector<double> &vals, int low, int high)
{
  if (high <= low)
    return (val > vals[low])?(low+1):low;

  int mid = (low + high)/2;

  if (val == vals[mid])
    return mid+1;

  if (val > vals[mid])
    return BinarySearch(val, vals, mid+1, high);
  return BinarySearch(val, vals, low, mid-1);
}

int BinarySearch(double val, std::vector<double> &vals)
{
  return BinarySearch(val, vals, 0, vals.size()-1);
}

int GeodesicMatting::GetMinValIndexInSorted(const std::vector<double>& vals, const std::vector<size_t>& sorted_idx)
{
  //  double d_min = 1e10;
  //  int idx = 0;
  //  for (size_t i = 0; i < vals.size(); i++)
  //  {
  //    if (d_min > vals[i])
  //    {
  //      d_min = vals[i];
  //      idx = i;
  //    }
  //  }
  //  return idx;
  return FindMinFromSorted(vals, sorted_idx, 0, vals.size()-1);
}

int GeodesicMatting::GetMinValIndex(const std::vector<double>& vals)
{
  double d_min = 1e10;
  int idx = 0;
  for (size_t i = 0; i < vals.size(); i++)
  {
    if (d_min > vals[i])
    {
      d_min = vals[i];
      idx = i;
    }
  }
  return idx;
}

template <typename T> vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

bool GeodesicMatting::Compute(int *dim, double *mri_in, double* mri_range_in, unsigned char *seeds_in, std::vector<unsigned char>& label_list,
                              unsigned char *seeds_out)
{
  double** lHood = new double*[label_list.size()];
  double** DISTS = new double*[label_list.size()];
  size_t vol_size = ((size_t)dim[0])*dim[1]*dim[2];
  for (size_t i = 0; i < label_list.size(); i++)
  {
    lHood[i] = new double[vol_size];
    memset(lHood[i], 0, vol_size*sizeof(double));
    DISTS[i] = new double[vol_size];
    memset(DISTS[i], 0, vol_size*sizeof(double));
  }

  int n_steps = 15;
  std::vector<double> v;
  double step_size = (mri_range_in[1]-mri_range_in[0])/n_steps;
  for (int i = 0; i <= n_steps; i++)
    v.push_back(mri_range_in[0] + i*step_size);

  // compute probabilities
  for (size_t l = 0; l < label_list.size(); l++)
  {
    KDE kde;
    for (size_t i = 0; i < vol_size; i++)
    {
      if (seeds_in[i] == label_list[l])
      {
        kde.add_data(mri_in[i]);
        //  qDebug() << mri_in[i];
      }
    }

    std::vector<double> hf;
    for (size_t i = 0; i < v.size(); i++)
    {
      hf.push_back(kde.pdf(v[i]));
      //  qDebug() << i << hf[i];
    }

    for (size_t i = 0; i < vol_size; i++)
    {
      lHood[l][i] = Interpolate(v, hf, mri_in[i]);
    }

    //    if (l == 0)
    //    {
    //      for (size_t i = 0; i < hf.size(); i++)
    //        qDebug() << hf[i];
    //    }
  }

  // normalize probabilities
  for (size_t i = 0; i < vol_size; i++)
  {
    double sum = 0;
    for (size_t j = 0; j < label_list.size(); j++)
      sum += lHood[j][i];
    if (sum > 0)
    {
      for (size_t j = 0; j < label_list.size(); j++)
        lHood[j][i] /= sum;
    }
  }

  double* D = new double[vol_size];
  double* TRIALVALS = new double[vol_size];
  unsigned char* KNOWN = new unsigned char[vol_size];
  unsigned char* TRIAL = new unsigned char[vol_size];
  QElapsedTimer timer;
  timer.start();
  for (size_t l = 0; l < label_list.size(); l++)
  {
    // initialize
    memset(D, 0, sizeof(double)*vol_size);
    memset(KNOWN, 0, vol_size);
    for (size_t i = 0; i < vol_size; i++)
    {
      if (seeds_in[i] == label_list[l])
      {
        KNOWN[i] = 1;
      }
    }
    Dilate(dim, KNOWN, TRIAL);
    for (size_t i = 0; i < vol_size; i++)
    {
      if (TRIAL[i])
      {
        TRIALVALS[i] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i%dim[0], (i/dim[0])%dim[1], i/(dim[0]*dim[1]));;
      }
      else
      {
        TRIALVALS[i]=LARGENUMBER;
      }
    }

    // fast marching
    bool ready = false;
    int npix = 0;
    while (!ready)
    {
      npix++;
      if (npix%10000 == 0)
        qDebug() << npix;

      double mini = LARGENUMBER;
      int idx = -1;
      for (size_t i = 0; i < vol_size; i++)
      {
        if (TRIALVALS[i]<mini)
        {
          mini=TRIALVALS[i];
          idx=i;
        }
      }
      if (idx < 0 || mini >= LARGENUMBER)
        break;

      KNOWN[idx]=1;
      TRIAL[idx]=0;
      TRIALVALS[idx]=LARGENUMBER;
      D[idx]=mini;

      size_t i = idx%dim[0];
      size_t j = (idx/dim[0])%dim[1];
      size_t k =  idx/(dim[0]*dim[1]);

      if (i > 0 && KNOWN[idx-1] == 0)
      {
        //        TRIAL[idx-1] = 1;
        //        TRIALVALS[idx-1] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i-1, j, k);

        if (TRIAL[idx-1] == 0)
        {
          TRIAL[idx-1] = 1;
          TRIALVALS[idx-1] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i-1, j, k);
        }
        else
        {
          double step = 0;
          for (int l = 0; l < label_list.size(); l++)
            step += qAbs(lHood[l][idx-1] - lHood[l][idx]);

          TRIALVALS[idx-1] = qMin(TRIALVALS[idx-1],mini+step);
        }
      }

      if (i < dim[0]-1 && KNOWN[idx+1] == 0)
      {
        //        TRIAL[idx+1] = 1;
        //        TRIALVALS[idx+1] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i+1, j, k);

        if (TRIAL[idx+1] == 0)
        {
          TRIAL[idx+1] = 1;
          TRIALVALS[idx+1] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i+1, j, k);
        }
        else
        {
          double step = 0;
          for (int l = 0; l < label_list.size(); l++)
            step += qAbs(lHood[l][idx+1] - lHood[l][idx]);

          TRIALVALS[idx+1] = qMin(TRIALVALS[idx+1],mini+step);
        }
      }

      if (j > 0 && KNOWN[idx-dim[0]] == 0)
      {
        //        TRIAL[idx-dim[0]] = 1;
        //        TRIALVALS[idx-dim[0]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j-1, k);

        if (TRIAL[idx-dim[0]] == 0)
        {
          TRIAL[idx-dim[0]] = 1;
          TRIALVALS[idx-dim[0]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j-1, k);
        }
        else
        {
          double step = 0;
          for (int l = 0; l < label_list.size(); l++)
            step += qAbs(lHood[l][idx-dim[0]] - lHood[l][idx]);

          TRIALVALS[idx-dim[0]] = qMin(TRIALVALS[idx-dim[0]],mini+step);
        }
      }

      if (j < dim[1]-1 && KNOWN[idx+dim[0]] == 0)
      {
        //        TRIAL[idx+dim[0]] = 1;
        //        TRIALVALS[idx+dim[0]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j+1, k);

        if (TRIAL[idx+dim[0]] == 0)
        {
          TRIAL[idx+dim[0]] = 1;
          TRIALVALS[idx+dim[0]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j+1, k);
        }
        else
        {
          double step = 0;
          for (int l = 0; l < label_list.size(); l++)
            step += qAbs(lHood[l][idx+dim[0]] - lHood[l][idx]);

          TRIALVALS[idx+dim[0]] = qMin(TRIALVALS[idx+dim[0]],mini+step);
        }
      }

      if (k > 0 && KNOWN[idx-dim[0]*dim[1]] == 0)
      {
        //        TRIAL[idx-dim[0]*dim[1]] = 1;
        //        TRIALVALS[idx-dim[0]*dim[1]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j, k-1);

        if (TRIAL[idx-dim[0]*dim[1]] == 0)
        {
          TRIAL[idx-dim[0]*dim[1]] = 1;
          TRIALVALS[idx-dim[0]*dim[1]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j, k-1);
        }
        else
        {
          double step = 0;
          for (int l = 0; l < label_list.size(); l++)
            step += qAbs(lHood[l][idx-dim[0]*dim[1]] - lHood[l][idx]);

          TRIALVALS[idx-dim[0]*dim[1]] = qMin(TRIALVALS[idx-dim[0]*dim[1]],mini+step);
        }
      }

      if (k < dim[2]-1  && KNOWN[idx+dim[0]*dim[1]] == 0)
      {
        //        TRIAL[idx+dim[0]*dim[1]] = 1;
        //        TRIALVALS[idx+dim[0]*dim[1]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j, k+1);

        if (TRIAL[idx+dim[0]*dim[1]] == 0)
        {
          TRIAL[idx+dim[0]*dim[1]] = 1;
          TRIALVALS[idx+dim[0]*dim[1]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j, k+1);
        }
        else
        {
          double step = 0;
          for (int l = 0; l < label_list.size(); l++)
            step += qAbs(lHood[l][idx+dim[0]*dim[1]] - lHood[l][idx]);

          TRIALVALS[idx+dim[0]*dim[1]] = qMin(TRIALVALS[idx+dim[0]*dim[1]],mini+step);
        }
      }

      if (npix>=vol_size)
        ready = true;
    }
    memcpy(DISTS[l], D, sizeof(double)*vol_size);
  }
  qDebug() << timer.elapsed()/1000.;

  memset(seeds_out, 0, vol_size);
  for (size_t i = 0; i < vol_size; i++)
  {
    double dmin = 1e10;
    int n = 0;
    for (size_t l = 0; l < label_list.size(); l++)
    {
      if (DISTS[l][i] < dmin)
      {
        dmin = DISTS[l][i];
        n = l;
      }
    }
    if (n != label_list.size()-1)
      seeds_out[i] = label_list[n];
  }
  delete[] D;
  delete[] TRIALVALS;
  delete[] KNOWN;
  delete[] TRIAL;

  for (size_t i = 0; i < label_list.size(); i++)
  {
    delete[] lHood[i];
    delete[] DISTS[i];
  }
  delete[] lHood;
  delete[] DISTS;

  return true;
}

bool GeodesicMatting::ComputeWithSorting(int *dim, double *mri_in, double* mri_range_in, unsigned char *seeds_in, std::vector<unsigned char>& label_list,
                                         unsigned char *seeds_out)
{
  double** lHood = new double*[label_list.size()];
  double** DISTS = new double*[label_list.size()];
  size_t vol_size = ((size_t)dim[0])*dim[1]*dim[2];
  for (size_t i = 0; i < label_list.size(); i++)
  {
    lHood[i] = new double[vol_size];
    memset(lHood[i], 0, vol_size*sizeof(double));
    DISTS[i] = new double[vol_size];
    memset(DISTS[i], 0, vol_size*sizeof(double));
  }

  int n_steps = 15;
  std::vector<double> v;
  double step_size = (mri_range_in[1]-mri_range_in[0])/n_steps;
  for (int i = 0; i <= n_steps; i++)
    v.push_back(mri_range_in[0] + i*step_size);

  // compute probabilities
  for (size_t l = 0; l < label_list.size(); l++)
  {
    KDE kde;
    for (size_t i = 0; i < vol_size; i++)
    {
      if (seeds_in[i] == label_list[l])
      {
        kde.add_data(mri_in[i]);
        //  qDebug() << mri_in[i];
      }
    }

    std::vector<double> hf;
    for (size_t i = 0; i < v.size(); i++)
    {
      hf.push_back(kde.pdf(v[i]));
      //  qDebug() << i << hf[i];
    }

    for (size_t i = 0; i < vol_size; i++)
    {
      lHood[l][i] = Interpolate(v, hf, mri_in[i]);
    }

    //    if (l == 0)
    //    {
    //      for (size_t i = 0; i < hf.size(); i++)
    //        qDebug() << hf[i];
    //    }
  }

  // normalize probabilities
  for (size_t i = 0; i < vol_size; i++)
  {
    double sum = 0;
    for (size_t j = 0; j < label_list.size(); j++)
      sum += lHood[j][i];
    if (sum > 0)
    {
      for (size_t j = 0; j < label_list.size(); j++)
        lHood[j][i] /= sum;
    }
  }

  double* D = new double[vol_size];
  unsigned char* KNOWN = new unsigned char[vol_size];
  unsigned char* TRIAL = new unsigned char[vol_size];
  double* TRIALVALS = new double[vol_size];
  QElapsedTimer timer;
  timer.start();
  for (size_t l = 0; l < label_list.size(); l++)
  {
    // initialize
    memset(D, 0, sizeof(double)*vol_size);
    memset(KNOWN, 0, vol_size);
    for (size_t i = 0; i < vol_size; i++)
    {
      if (seeds_in[i] == label_list[l])
      {
        KNOWN[i] = 1;
      }
    }
    Dilate(dim, KNOWN, TRIAL);
    std::vector<size_t> idxs;
    for (size_t i = 0; i < vol_size; i++)
    {
      if (TRIAL[i])
      {
        TRIALVALS[i] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i%dim[0], (i/dim[0])%dim[1], i/(dim[0]*dim[1]));
        idxs.push_back(i);
      }
      else
      {
        TRIALVALS[i] = LARGENUMBER;
      }
    }

    //    std::vector<size_t> sorted_idx = sort_indexes(vals);

    // fast marching
    bool ready = false;
    int npix = 0;
    while (!ready)
    {
      npix++;
      if (npix%10000 == 0)
        qDebug() << npix;

      double mini = LARGENUMBER;
      size_t idx = -1;
      int n_idxs = -1;
      //      for (size_t i = 0; i < idxs.size(); i++)
      //      {
      //        if (TRIALVALS[idxs[i]] < mini)
      //        {
      //          mini = TRIALVALS[idxs[i]];
      //          idx = idxs[i];
      //          n_idxs = i;
      //        }
      //      }

      for (size_t i = 0; i < vol_size; i++)
      {
        if (TRIALVALS[i]<mini)
        {
          mini=TRIALVALS[i];
          idx=i;
        }
      }

      if (idx < 0 || mini >= LARGENUMBER)
        break;

      KNOWN[idx] = 1;
      TRIAL[idx] = 0;
      TRIALVALS[idx] = LARGENUMBER;
      D[idx] = mini;
      //      idxs.erase(idxs.begin()+n_idxs);

      size_t i = idx%dim[0];
      size_t j = (idx/dim[0])%dim[1];
      size_t k =  idx/(dim[0]*dim[1]);

      if (i > 0 && KNOWN[idx-1] == 0)
      {
        //        TRIAL[idx-1] = 1;
        //        TRIALVALS[idx-1] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i-1, j, k);

        if (TRIAL[idx-1] == 0)
        {
          TRIAL[idx-1] = 1;
          TRIALVALS[idx-1] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i-1, j, k);
        }
        else
        {
          double step = 0;
          for (int l = 0; l < label_list.size(); l++)
            step += qAbs(lHood[l][idx-1] - lHood[l][idx]);

          TRIALVALS[idx-1] = qMin(TRIALVALS[idx-1],mini+step);
        }
        //  idxs.push_back(idx-1);
      }

      if (i < dim[0]-1 && KNOWN[idx+1] == 0)
      {
        //        TRIAL[idx+1] = 1;
        //        TRIALVALS[idx+1] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i+1, j, k);

        if (TRIAL[idx+1] == 0)
        {
          TRIAL[idx+1] = 1;
          TRIALVALS[idx+1] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i+1, j, k);
        }
        else
        {
          double step = 0;
          for (int l = 0; l < label_list.size(); l++)
            step += qAbs(lHood[l][idx+1] - lHood[l][idx]);

          TRIALVALS[idx+1] = qMin(TRIALVALS[idx+1],mini+step);
        }
        //   idxs.push_back(idx+1);
      }

      if (j > 0 && KNOWN[idx-dim[0]] == 0)
      {
        //        TRIAL[idx-dim[0]] = 1;
        //        TRIALVALS[idx-dim[0]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j-1, k);

        if (TRIAL[idx-dim[0]] == 0)
        {
          TRIAL[idx-dim[0]] = 1;
          TRIALVALS[idx-dim[0]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j-1, k);
        }
        else
        {
          double step = 0;
          for (int l = 0; l < label_list.size(); l++)
            step += qAbs(lHood[l][idx-dim[0]] - lHood[l][idx]);

          TRIALVALS[idx-dim[0]] = qMin(TRIALVALS[idx-dim[0]],mini+step);
        }
        //  idxs.push_back(idx-dim[0]);
      }

      if (j < dim[1]-1 && KNOWN[idx+dim[0]] == 0)
      {
        //        TRIAL[idx+dim[0]] = 1;
        //        TRIALVALS[idx+dim[0]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j+1, k);

        if (TRIAL[idx+dim[0]] == 0)
        {
          TRIAL[idx+dim[0]] = 1;
          TRIALVALS[idx+dim[0]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j+1, k);
        }
        else
        {
          double step = 0;
          for (int l = 0; l < label_list.size(); l++)
            step += qAbs(lHood[l][idx+dim[0]] - lHood[l][idx]);

          TRIALVALS[idx+dim[0]] = qMin(TRIALVALS[idx+dim[0]],mini+step);
        }
        //  idxs.push_back(idx+dim[0]);
      }

      if (k > 0 && KNOWN[idx-dim[0]*dim[1]] == 0)
      {
        //        TRIAL[idx-dim[0]*dim[1]] = 1;
        //        TRIALVALS[idx-dim[0]*dim[1]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j, k-1);

        if (TRIAL[idx-dim[0]*dim[1]] == 0)
        {
          TRIAL[idx-dim[0]*dim[1]] = 1;
          TRIALVALS[idx-dim[0]*dim[1]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j, k-1);
        }
        else
        {
          double step = 0;
          for (int l = 0; l < label_list.size(); l++)
            step += qAbs(lHood[l][idx-dim[0]*dim[1]] - lHood[l][idx]);

          TRIALVALS[idx-dim[0]*dim[1]] = qMin(TRIALVALS[idx-dim[0]*dim[1]],mini+step);
        }
        //  idxs.push_back(idx-dim[0]*dim[1]);
      }

      if (k < dim[2]-1  && KNOWN[idx+dim[0]*dim[1]] == 0)
      {
        //        TRIAL[idx+dim[0]*dim[1]] = 1;
        //        TRIALVALS[idx+dim[0]*dim[1]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j, k+1);

        if (TRIAL[idx+dim[0]*dim[1]] == 0)
        {
          TRIAL[idx+dim[0]*dim[1]] = 1;
          TRIALVALS[idx+dim[0]*dim[1]] = ComputeNeighDist(lHood, label_list.size(), KNOWN, D, dim, i, j, k+1);
        }
        else
        {
          double step = 0;
          for (int l = 0; l < label_list.size(); l++)
            step += qAbs(lHood[l][idx+dim[0]*dim[1]] - lHood[l][idx]);

          TRIALVALS[idx+dim[0]*dim[1]] = qMin(TRIALVALS[idx+dim[0]*dim[1]],mini+step);
        }
        //  idxs.push_back(idx+dim[0]*dim[1]);
      }

      if (npix>=vol_size)
        ready = true;
    }
    memcpy(DISTS[l], D, sizeof(double)*vol_size);
  }
  qDebug() << timer.elapsed()/1000.;

  memset(seeds_out, 0, vol_size);
  for (size_t i = 0; i < vol_size; i++)
  {
    double dmin = 1e10;
    int n = 0;
    for (size_t l = 0; l < label_list.size(); l++)
    {
      if (DISTS[l][i] < dmin)
      {
        dmin = DISTS[l][i];
        n = l;
      }
    }
    if (n != label_list.size()-1)
      seeds_out[i] = label_list[n];
  }
  delete[] D;
  delete[] KNOWN;
  delete[] TRIAL;
  delete[] TRIALVALS;

  for (size_t i = 0; i < label_list.size(); i++)
  {
    delete[] lHood[i];
    delete[] DISTS[i];
  }
  delete[] lHood;
  delete[] DISTS;
  return true;
}

