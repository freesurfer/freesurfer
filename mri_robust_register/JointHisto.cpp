/**
 * @brief A class for a joint histogram of two images
 *
 */

/*
 * Original Author: Martin Reuter
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include "JointHisto.h"
#if ITK_VERSION_MAJOR >= 5
#include <iostream.h>
#include <fstream.h>
#include <vcl_compiler.h>
#else
#include <vcl_iostream.h>
#include <vcl_fstream.h>
#endif
#include <vnl/vnl_matlab_print.h>
#include <iomanip>

using namespace std;

void JointHisto::create(MRI *mri1, MRI * mri2, int d1, int d2, int d3)
// images are expected to be uchar 0..255 and have uniform voxels and same dimensions
{
  //cout << " JointHisto::create " << endl;

  int width = mri1->width;
  int height = mri1->height;
  int depth = mri1->depth;
  assert(width == mri2->width);
  assert(height == mri2->height);
  assert(depth == mri2->depth);
  assert(mri1->type == MRI_UCHAR);
  assert(mri2->type == MRI_UCHAR);
  assert(mri1->xsize == mri1->ysize);
  assert(mri1->xsize == mri1->zsize);
  assert(mri1->xsize == mri2->xsize);
  assert(mri1->xsize == mri2->ysize);
  assert(mri1->xsize == mri2->zsize);

  histo.fill(0.0);
  haseps = false;

  int x, y, z;
  int isum = 0;
  for (z = 0; z <= depth - d3; z += d3)
    for (y = 0; y <= height - d2; y += d2)
      for (x = 0; x <= width - d1; x += d1)
      {
        histo[(int) MRIvox(mri1,x,y,z) ][(int)MRIvox(mri2,x,y,z)] += 1;
        isum++;
      }
  sum = isum;
  //cout <<"Sum: " << sum << endl;

}

void JointHisto::save(const std::string & fname, const std::string & n)
{
  vcl_ostream* s = new vcl_ofstream(fname.c_str());
  if (!s || !(*s))
  {
    vcl_cerr << "Cannot open " << fname << " for writing\n";
  }
  else
  {
    vnl_matlab_print(*s,histo,n.c_str(),vnl_matlab_print_format_long);
    *s << "\n";
  }
  delete s; s = 0;
}

double JointHisto::clip(double thres)
{

  int i, j;
  double csum = 0.0;
  thres *= sum;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      if (histo[i][j] < thres)
      {
        csum += histo[i][j];
        histo[i][j] = 0.0;
      }
    }

  double cratio = csum / sum;

  sum -= csum;

  return cratio;

}

pair<double, double> JointHisto::getMinMax()
{
  int i, j;
  pair<double, double> mm(-1, -1);
  double tmp;
  if (n == 0)
    return mm;
  mm.first = histo[0][0];
  mm.second = mm.first;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      tmp = histo[i][j];
      if (mm.first > tmp)
        mm.first = tmp;
      if (mm.second < tmp)
        mm.second = tmp;
    }

  return mm;
}

void JointHisto::create(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3)
{
//  cout << " JointHisto::create " << endl;
//   vnl_matlab_print(vcl_cerr,Msi,"Msi",vnl_matlab_print_format_long);std::cerr << std::endl;
//   vnl_matlab_print(vcl_cerr,Mti,"Mti",vnl_matlab_print_format_long);std::cerr << std::endl;
//  cout << " d1: " << d1 << " d2: " << d2 << " d3: " << d3 << endl;
  //   assert (mriS->type == MRI_UCHAR);
//   assert (mriT->type == MRI_UCHAR);
  //int dim = 4;
  int dim = 1;
  n = 256 / dim;
  int nm1 = n - 1;
  histo.set_size(n, n);
  histo.fill(0.0);
  haseps = false;
  sum = 0.0;

//  MRIwrite(mriS,"histoS.mgz");
//  MRIwrite(mriT,"histoT.mgz");
//  vnl_matlab_print(vcl_cerr,Msi,"Msi",vnl_matlab_print_format_long);
//  vnl_matlab_print(vcl_cerr,Mti,"Mti",vnl_matlab_print_format_long);

//   static double uniform[101] = { 0.825960894329359, 0.391799656122844, 0.606822280998442, 0.923192169644637,
//                                  0.772662867376799, 0.394913419972908, 0.562648884384293, 0.821213453107628,
//                                  0.077777173551299, 0.536306917657496, 0.514107239245956, 0.761227511728186,
//                                  0.978425890910315, 0.157030594951506, 0.147496203406248, 0.783403101626463,
//                                  0.637363106003275, 0.257349279177570, 0.109403726553149, 0.136954287985675,
//                                  0.654171911591805, 0.204476823923606, 0.164700938803686, 0.566376934840801,
//                                  0.854095286894540, 0.430581644078462, 0.397193805113614, 0.045536084142315,
//                                  0.140620297119782, 0.724179303888238, 0.210184975302471, 0.346103835004574,
//                                  0.226683913788046, 0.308131967172401, 0.451070768412025, 0.831843821825446,
//                                  0.086760750364257, 0.854135129404520, 0.742590231477529, 0.053858310273945,
//                                  0.122791324135168, 0.526967625411493, 0.320307444448230, 0.520062337463421,
//                                  0.878893684505441, 0.443226585521493, 0.320014638636649, 0.868171615160687,
//                                  0.764138797526944, 0.864611801159897, 0.083244000805187, 0.340050247951733,
//                                  0.581607039757426, 0.291315250124114, 0.467637373824935, 0.596640965216293,
//                                  0.897387480194169, 0.121159463814373, 0.018178277226624, 0.747666956017887,
//                                  0.670108792666926, 0.092391423049263, 0.917528192569663, 0.347581829497044,
//                                  0.857229173664249, 0.632564501155682, 0.363593396666436, 0.385827512659221,
//                                  0.003355759148503, 0.969283884272977, 0.459201707491767, 0.800715562489883,
//                                  0.296945204007638, 0.619917383791628, 0.748749776544515, 0.394902341568457,
//                                  0.354479046607957, 0.867388084434840, 0.627567204750354, 0.984294463747630,
//                                  0.824472893444168, 0.295272747204669, 0.927251152590456, 0.119266587069099,
//                                  0.241265513320426, 0.527775796652046, 0.060203502196487, 0.835363100840541,
//                                  0.148398316485924, 0.121273963075904, 0.683207184266160, 0.500002874704566,
//                                  0.561265939626174, 0.234256657698203, 0.486854859925734, 0.071141206079346,
//                                  0.442630693187859, 0.327200604299592, 0.226827609433358, 0.971944076189183,
//                                  0.302612670611030 };
//   int opos = 0;
// 
//    int    upos=0;
//   bool   ssamp = (d1>1 || d2 > 1 || d3 >1);

  int x, y, z;
  double ux, uy, uz;
  double xs, ys, zs;
  double xt, yt, zt;
  double vs, vt;
  int ivs, ivt;
  double sdiff, tdiff;
  double dt[3] =
    { (double)mriT->width, (double)mriT->height, (double)mriT->depth };
  // double ds[3] = { mriS->width, mriS->height, mriS->depth};
  mriS->outside_val = -1;
  mriT->outside_val = -1;
  //int count = 0;
//  cout <<" dt " << dt[0] << " " << dt[1] << " " << dt[2] << endl;
  for (z = 0; z < dt[2] - d3 + 1; z += d3)
  {
    for (y = 0; y < dt[1] - d2 + 1; y += d2)
    {
      for (x = 0; x < dt[0] - d1 + 1; x += d1)
      {
//         if (ssamp)
//         {
//           ux  = x + uniform[upos]*d1;
//           upos = (upos+1)%101;
//           uy  = y + uniform[upos]*d2;
//           upos = (upos+1)%101;
//           uz  = z + uniform[upos]*d3;
//           upos = (upos+1)%101;
//         }
//         else
        {
          ux = x;
          uy = y;
          uz = z;
        }
        xt = Mti[0][0] * ux + Mti[0][1] * uy + Mti[0][2] * uz + Mti[0][3];
        yt = Mti[1][0] * ux + Mti[1][1] * uy + Mti[1][2] * uz + Mti[1][3];
        zt = Mti[2][0] * ux + Mti[2][1] * uy + Mti[2][2] * uz + Mti[2][3];

        xs = Msi[0][0] * ux + Msi[0][1] * uy + Msi[0][2] * uz + Msi[0][3];
        ys = Msi[1][0] * ux + Msi[1][1] * uy + Msi[1][2] * uz + Msi[1][3];
        zs = Msi[2][0] * ux + Msi[2][1] * uy + Msi[2][2] * uz + Msi[2][3];

        //if (zs>=0.0 && zs<ds[2]-1 && ys>=0.0 && ys<ds[1]-1 && xs>=0.0 && xs<ds[0]-1 &&
        //    zt>=0.0 && zt<dt[2]-1 && yt>=0.0 && yt<dt[1]-1 && xt>=0.0 && xt<dt[0]-1 )
        {
          MRIsampleVolumeFrame(mriS, xs, ys, zs, 0, &vs);
          if (vs == -1)
            continue;
          if (dim != 1)
            vs /= dim;
          //if (mask && vs == 0.0) continue;  // bad idea, try registering skull stripped to full
          ivs = (int) floor(vs);
          MRIsampleVolumeFrame(mriT, xt, yt, zt, 0, &vt);
          if (vt == -1)
            continue;
          if (dim != 1)
            vt /= dim;
          ivt = (int) floor(vt);
          //cout << "ivs " << ivs << " ivt " << ivt << endl;
          if (ivs < 0 || ivt < 0 || ivs >= n || ivt >= n)
          {
            cout << "ERROR at ( " << x << " " << y << " " << z << " ) " << endl;
            cout << " PS: ( " << xs << " " << ys << " " << zs << " )  PT: ( "
                << xt << " " << yt << " " << zt << " )" << endl;
            cout << " ivs: " << ivs << "  ivt: " << ivt << "  n: " << n << endl;
            exit(1);
          }
          sdiff = vs - ivs;
          tdiff = vt - ivt;
          //cout << "ivs " << ivs << " ivt " << ivt << endl;
          // distribute peak among bins symetrically:
          histo[ivs][ivt] += (1 - sdiff) * (1 - tdiff);
          if (ivs < nm1)
            histo[ivs + 1][ivt] += sdiff * (1 - tdiff);
          if (ivt < nm1)
            histo[ivs][ivt + 1] += (1 - sdiff) * tdiff;
          if (ivs < nm1 && ivt < nm1)
            histo[ivs + 1][ivt + 1] += sdiff * tdiff;
          //count++;
          sum = sum + 1.0;
        }
      }
    }
  }
  //cout << " count: " << count << endl;
  if (sum == 0.0)
  {
    cout << " WARNING JointHisto::create : histogram is empty (no overlap?) "
        << endl;
    //cout << " ERROR JointHisto::create : histogram is empty? Writing debug images. " << endl;
    //MRIwrite(mriS,"mriS_error.mgz");
    //MRIwrite(mriT,"mriT_error.mgz");
    //exit(1);
  }
}

MRI * JointHisto::locate(MRI *mriS, MRI* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3)
{
  //int dim = 4;
  int dim = 1;
  n = 256 / dim;

  normalize();
  pair<double, double> minmax = getMinMax();
  double min = minmax.first;
  double scale = 1.0 / (minmax.second - minmax.first);

//  MRIwrite(mriS,"histoS.mgz");
//  MRIwrite(mriT,"histoT.mgz");
//  vnl_matlab_print(vcl_cerr,Msi,"Msi",vnl_matlab_print_format_long);
//  vnl_matlab_print(vcl_cerr,Mti,"Mti",vnl_matlab_print_format_long);

//   static double uniform[101] = { 0.825960894329359, 0.391799656122844, 0.606822280998442, 0.923192169644637,
//                                  0.772662867376799, 0.394913419972908, 0.562648884384293, 0.821213453107628,
//                                  0.077777173551299, 0.536306917657496, 0.514107239245956, 0.761227511728186,
//                                  0.978425890910315, 0.157030594951506, 0.147496203406248, 0.783403101626463,
//                                  0.637363106003275, 0.257349279177570, 0.109403726553149, 0.136954287985675,
//                                  0.654171911591805, 0.204476823923606, 0.164700938803686, 0.566376934840801,
//                                  0.854095286894540, 0.430581644078462, 0.397193805113614, 0.045536084142315,
//                                  0.140620297119782, 0.724179303888238, 0.210184975302471, 0.346103835004574,
//                                  0.226683913788046, 0.308131967172401, 0.451070768412025, 0.831843821825446,
//                                  0.086760750364257, 0.854135129404520, 0.742590231477529, 0.053858310273945,
//                                  0.122791324135168, 0.526967625411493, 0.320307444448230, 0.520062337463421,
//                                  0.878893684505441, 0.443226585521493, 0.320014638636649, 0.868171615160687,
//                                  0.764138797526944, 0.864611801159897, 0.083244000805187, 0.340050247951733,
//                                  0.581607039757426, 0.291315250124114, 0.467637373824935, 0.596640965216293,
//                                  0.897387480194169, 0.121159463814373, 0.018178277226624, 0.747666956017887,
//                                  0.670108792666926, 0.092391423049263, 0.917528192569663, 0.347581829497044,
//                                  0.857229173664249, 0.632564501155682, 0.363593396666436, 0.385827512659221,
//                                  0.003355759148503, 0.969283884272977, 0.459201707491767, 0.800715562489883,
//                                  0.296945204007638, 0.619917383791628, 0.748749776544515, 0.394902341568457,
//                                  0.354479046607957, 0.867388084434840, 0.627567204750354, 0.984294463747630,
//                                  0.824472893444168, 0.295272747204669, 0.927251152590456, 0.119266587069099,
//                                  0.241265513320426, 0.527775796652046, 0.060203502196487, 0.835363100840541,
//                                  0.148398316485924, 0.121273963075904, 0.683207184266160, 0.500002874704566,
//                                  0.561265939626174, 0.234256657698203, 0.486854859925734, 0.071141206079346,
//                                  0.442630693187859, 0.327200604299592, 0.226827609433358, 0.971944076189183,
//                                  0.302612670611030 };
//   int opos = 0;
// 
//    int    upos=0;
//   bool   ssamp = (d1>1 || d2 > 1 || d3 >1);

  int x, y, z, i, maxi;
  double ux, uy, uz;
  double xs, ys, zs;
  double xt, yt, zt;
  double vs, vt;
  int ivs, ivt;
  double sdiff, tdiff;
  double dt[3] =
  { (double)mriT->width, (double)mriT->height, (double)mriT->depth };
  double four[4];
  double maxfour;
  MRI * weights = MRIcloneDifferentType(mriT, MRI_FLOAT);

  // double ds[3] = { mriS->width, mriS->height, mriS->depth};
  mriS->outside_val = -1;
  mriT->outside_val = -1;
  //int count = 0;
  //cout <<" df " << df[0] << " " << df[1] << " " << df[2] << endl;
  for (z = 0; z < dt[2] - d3 + 1; z += d3)
  {
    for (y = 0; y < dt[1] - d2 + 1; y += d2)
    {
      for (x = 0; x < dt[0] - d1 + 1; x += d1)
      {
//         if (ssamp)
//         {
//           ux  = x + uniform[upos]*d1;
//           upos = (upos+1)%101;
//           uy  = y + uniform[upos]*d2;
//           upos = (upos+1)%101;
//           uz  = z + uniform[upos]*d3;
//           upos = (upos+1)%101;
//         }
//         else
        {
          ux = x;
          uy = y;
          uz = z;
        }
        xt = Mti[0][0] * ux + Mti[0][1] * uy + Mti[0][2] * uz + Mti[0][3];
        yt = Mti[1][0] * ux + Mti[1][1] * uy + Mti[1][2] * uz + Mti[1][3];
        zt = Mti[2][0] * ux + Mti[2][1] * uy + Mti[2][2] * uz + Mti[2][3];

        xs = Msi[0][0] * ux + Msi[0][1] * uy + Msi[0][2] * uz + Msi[0][3];
        ys = Msi[1][0] * ux + Msi[1][1] * uy + Msi[1][2] * uz + Msi[1][3];
        zs = Msi[2][0] * ux + Msi[2][1] * uy + Msi[2][2] * uz + Msi[2][3];

        //if (zs>=0.0 && zs<ds[2]-1 && ys>=0.0 && ys<ds[1]-1 && xs>=0.0 && xs<ds[0]-1 &&
        //    zt>=0.0 && zt<dt[2]-1 && yt>=0.0 && yt<dt[1]-1 && xt>=0.0 && xt<dt[0]-1 )
        {
          MRIsampleVolumeFrame(mriS, xs, ys, zs, 0, &vs);
          if (vs == -1)
            continue;
          if (dim != 1)
            vs /= dim;
          //if (mask && vs == 0.0) continue;  // bad idea, try registering skull stripped to full
          ivs = (int) floor(vs);
          MRIsampleVolumeFrame(mriT, xt, yt, zt, 0, &vt);
          if (vt == -1)
            continue;
          if (dim != 1)
            vt /= dim;
          ivt = (int) floor(vt);
          if (ivs < 0 || ivt < 0 || ivs >= n || ivt >= n)
          {
            cout << "ERROR at ( " << x << " " << y << " " << z << " ) " << endl;
            cout << " PS: ( " << xs << " " << ys << " " << zs << " )  PT: ( "
                << xt << " " << yt << " " << zt << " )" << endl;
            cout << " ivs: " << ivs << "  ivt: " << ivt << "  n: " << n << endl;
            exit(1);
          }
          sdiff = vs - ivs;
          tdiff = vt - ivt;
          // distribute peak among bins symetrically:
          four[0] = (1 - sdiff) * (1 - tdiff);
          four[1] = sdiff * (1 - tdiff);
          four[2] = (1 - sdiff) * tdiff;
          four[3] = sdiff * tdiff;
          maxfour = four[3];
          maxi = 3;
          for (i = 0; i < 3; i++)
            if (four[i] > maxfour)
            {
              maxfour = four[i];
              maxi = i;
            }

          switch (maxi)
          {
          case 0:
            MRIsetVoxVal(weights, x, y, z, 0,
                1 - (histo[ivs][ivt] - min) * scale);
            break;
          case 1:
            MRIsetVoxVal(weights, x, y, z, 0,
                1 - (histo[ivs + 1][ivt] - min) * scale);
            break;
          case 2:
            MRIsetVoxVal(weights, x, y, z, 0,
                1 - (histo[ivs][ivt + 1] - min) * scale);
            break;
          case 3:
            MRIsetVoxVal(weights, x, y, z, 0,
                1 - (histo[ivs + 1][ivt + 1] - min) * scale);
            break;

          }
        }
      }
    }
  }
  //cout << " count: " << count << endl;
  if (sum == 0.0)
  {
    cout << " ERROR JointHisto::create : histogram is empty? " << endl;
    exit(1);
  }

  return weights;
}

void JointHisto::create(MRI_BSPLINE *mriS, MRI_BSPLINE* mriT,
    const vnl_matrix_fixed<double, 4, 4>& Msi,
    const vnl_matrix_fixed<double, 4, 4>& Mti, int d1, int d2, int d3)
{
  //int dim = 4;
  int dim = 1;
  n = 256 / dim;
  int nm1 = n - 1;
  histo.set_size(n, n);
  histo.fill(0.0);
  haseps = false;
  sum = 0.0;

//  MRIwrite(mriS,"histoS.mgz");
//  MRIwrite(mriT,"histoT.mgz");
//  vnl_matlab_print(vcl_cerr,Msi,"Msi",vnl_matlab_print_format_long);
//  vnl_matlab_print(vcl_cerr,Mti,"Mti",vnl_matlab_print_format_long);

//   static double uniform[101] = { 0.825960894329359, 0.391799656122844, 0.606822280998442, 0.923192169644637,
//                                  0.772662867376799, 0.394913419972908, 0.562648884384293, 0.821213453107628,
//                                  0.077777173551299, 0.536306917657496, 0.514107239245956, 0.761227511728186,
//                                  0.978425890910315, 0.157030594951506, 0.147496203406248, 0.783403101626463,
//                                  0.637363106003275, 0.257349279177570, 0.109403726553149, 0.136954287985675,
//                                  0.654171911591805, 0.204476823923606, 0.164700938803686, 0.566376934840801,
//                                  0.854095286894540, 0.430581644078462, 0.397193805113614, 0.045536084142315,
//                                  0.140620297119782, 0.724179303888238, 0.210184975302471, 0.346103835004574,
//                                  0.226683913788046, 0.308131967172401, 0.451070768412025, 0.831843821825446,
//                                  0.086760750364257, 0.854135129404520, 0.742590231477529, 0.053858310273945,
//                                  0.122791324135168, 0.526967625411493, 0.320307444448230, 0.520062337463421,
//                                  0.878893684505441, 0.443226585521493, 0.320014638636649, 0.868171615160687,
//                                  0.764138797526944, 0.864611801159897, 0.083244000805187, 0.340050247951733,
//                                  0.581607039757426, 0.291315250124114, 0.467637373824935, 0.596640965216293,
//                                  0.897387480194169, 0.121159463814373, 0.018178277226624, 0.747666956017887,
//                                  0.670108792666926, 0.092391423049263, 0.917528192569663, 0.347581829497044,
//                                  0.857229173664249, 0.632564501155682, 0.363593396666436, 0.385827512659221,
//                                  0.003355759148503, 0.969283884272977, 0.459201707491767, 0.800715562489883,
//                                  0.296945204007638, 0.619917383791628, 0.748749776544515, 0.394902341568457,
//                                  0.354479046607957, 0.867388084434840, 0.627567204750354, 0.984294463747630,
//                                  0.824472893444168, 0.295272747204669, 0.927251152590456, 0.119266587069099,
//                                  0.241265513320426, 0.527775796652046, 0.060203502196487, 0.835363100840541,
//                                  0.148398316485924, 0.121273963075904, 0.683207184266160, 0.500002874704566,
//                                  0.561265939626174, 0.234256657698203, 0.486854859925734, 0.071141206079346,
//                                  0.442630693187859, 0.327200604299592, 0.226827609433358, 0.971944076189183,
//                                  0.302612670611030 };
//   int opos = 0;
// 
//    int    upos=0;
//   bool   ssamp = (d1>1 || d2 > 1 || d3 >1);

  int x, y, z;
  double ux, uy, uz;
  double xs, ys, zs;
  double xt, yt, zt;
  double vs, vt;
  int ivs, ivt;
  double sdiff, tdiff;
  double dt[3] =
  { (double)mriT->coeff->width, (double)mriT->coeff->height, (double)mriT->coeff->depth };
  // double ds[3] = { mriS->width, mriS->height, mriS->depth};
  mriS->coeff->outside_val = -1;
  mriT->coeff->outside_val = -1;
  //int count = 0;
  //cout <<" df " << df[0] << " " << df[1] << " " << df[2] << endl;
  for (z = 0; z < dt[2] - d3 + 1; z += d3)
  {
    for (y = 0; y < dt[1] - d2 + 1; y += d2)
    {
      for (x = 0; x < dt[0] - d1 + 1; x += d1)
      {
//         if (ssamp)
//         {
//           ux  = x + uniform[upos]*d1;
//           upos = (upos+1)%101;
//           uy  = y + uniform[upos]*d2;
//           upos = (upos+1)%101;
//           uz  = z + uniform[upos]*d3;
//           upos = (upos+1)%101;
//         }
//         else
        {
          ux = x;
          uy = y;
          uz = z;
        }
        xt = Mti[0][0] * ux + Mti[0][1] * uy + Mti[0][2] * uz + Mti[0][3];
        yt = Mti[1][0] * ux + Mti[1][1] * uy + Mti[1][2] * uz + Mti[1][3];
        zt = Mti[2][0] * ux + Mti[2][1] * uy + Mti[2][2] * uz + Mti[2][3];

        xs = Msi[0][0] * ux + Msi[0][1] * uy + Msi[0][2] * uz + Msi[0][3];
        ys = Msi[1][0] * ux + Msi[1][1] * uy + Msi[1][2] * uz + Msi[1][3];
        zs = Msi[2][0] * ux + Msi[2][1] * uy + Msi[2][2] * uz + Msi[2][3];

        //if (zs>=0.0 && zs<ds[2]-1 && ys>=0.0 && ys<ds[1]-1 && xs>=0.0 && xs<ds[0]-1 &&
        //    zt>=0.0 && zt<dt[2]-1 && yt>=0.0 && yt<dt[1]-1 && xt>=0.0 && xt<dt[0]-1 )
        {
          MRIsampleBSpline(mriS, xs, ys, zs, 0, &vs);
          if (vs == -1)
            continue;
          if (dim != 1)
            vs /= dim;
          //if (mask && vs == 0.0) continue;  // bad idea, try registering skull stripped to full
          ivs = (int) floor(vs);
          MRIsampleBSpline(mriT, xt, yt, zt, 0, &vt);
          if (vt == -1)
            continue;
          if (dim != 1)
            vt /= dim;
          ivt = (int) floor(vt);
          if (ivs < 0 || ivt < 0 || ivs >= n || ivt >= n)
          {
            cout << "ERROR at ( " << x << " " << y << " " << z << " ) " << endl;
            cout << " PS: ( " << xs << " " << ys << " " << zs << " )  PT: ( "
                << xt << " " << yt << " " << zt << " )" << endl;
            cout << " ivs: " << ivs << "  ivt: " << ivt << "  n: " << n << endl;
            exit(1);
          }
          sdiff = vs - ivs;
          tdiff = vt - ivt;
          // distribute peak among bins symetrically:
          histo[ivs][ivt] += (1 - sdiff) * (1 - tdiff);
          if (ivs < nm1)
            histo[ivs + 1][ivt] += sdiff * (1 - tdiff);
          if (ivt < nm1)
            histo[ivs][ivt + 1] += (1 - sdiff) * tdiff;
          if (ivs < nm1 && ivt < nm1)
            histo[ivs + 1][ivt + 1] += sdiff * tdiff;
          //count++;
          sum = sum + 1.0;
        }
      }
    }
  }
  //cout << " count: " << count << endl;
  if (sum == 0.0)
  {
    cout << " ERROR JointHisto::create : histogram is empty? " << endl;
    exit(1);
  }
}

void JointHisto::set(const vnl_matrix<double> & h)
// set histo matrix directly (for debugging)
{
  assert(h.rows() ==h.cols());
  histo = h;
  n = h.rows();
  haseps = false;
  sum = 0.0;
  double sum1 = 0.0;
  int i, j;
  for (i = 0; i < n; i++)
  {
    sum1 = 0.0;
    for (j = 0; j < n; j++)
    {
      sum1 += histo[i][j];
    }
    sum += sum1;
  }
}

void JointHisto::smooth(double fwhm)
{
  // compute sigma from fwhm:
  double sm = fwhm / sqrt(8 * log(2));
  if (sm < 0.001)
    sm = 0.001;
  int t = (int) floor(3 * sm + 0.5); // rounding 3*sm
  // create kernel
  double filter[2 * t + 1];
  double fsum = 0.0;
  double sm2 = -sm * sm;
  //cout << " sm = " << sm << endl;
  for (int i = -t; i <= t; i++)
  {
    filter[i + t] = exp((i * i) / sm2);
    fsum += filter[i + t];
  }
  // normalize
  for (int i = 0; i <= 2 * t; i++)
  {
    filter[i] /= fsum;
    //cout << setprecision(16)<< " f[" << i << "]= " << filter[i] << endl;
  }

  // convolve 2D
  vnl_matrix<double> htmp(n, n);
  double dtmp;
  int i, j, k, count;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      dtmp = 0;
      count = -1;
      for (k = j - t; k <= j + t; k++)
      {
        count++;
        if (k < 0 || k >= n)
          continue;
        dtmp += histo[i][k] * filter[count];
      }
      htmp[i][j] = dtmp;
    }
  sum = 0.0;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      dtmp = 0;
      count = -1;
      for (k = i - t; k <= i + t; k++)
      {
        count++;
        if (k < 0 || k >= n)
          continue;
        dtmp += htmp[k][j] * filter[count];
      }
      histo[i][j] = dtmp;
      sum += dtmp;
    }
}

double JointHisto::computeMI()
// compute Mutual Information
// Collignon, Maes, Delaere, Vandermeulen, Suetens & Marchal (1995).
// "Automated multi-modality image registration based on information theory".
// In Bizais, Barillot & Di Paola, editors, Proc. Information Processing
// in Medical Imaging, pages 263--274, Dordrecht, The Netherlands, 1995.
// Kluwer Academic Publishers.
//
// Wells III, Viola, Atsumi, Nakajima & Kikinis (1996).
// "Multi-modal volume registration by maximisation of mutual information".
// Medical Image Analysis, 1(1):35-51, 1996. 
// mi =  sum_{i,j} [ P(i,j) log2 (P(i,j) /(Pr(i) Pc(j)) ) ]
// mi = H(M) + H(N) - H(M,N)
{
  int i, j;
  double mi = 0.0;
  double d;
  //addeps(2.2204E-16);
  normalize();
  computeRCsums();
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      d = histo[i][j];
      if (d != 0.0)
        mi += d * log2(d / (rowsum[i] * colsum[j]));
    }
  return mi;
}

double JointHisto::computeECC()
// compute Entropy Correlation Coefficient
// F Maes, A Collignon, D Vandermeulen, G Marchal & P Suetens (1997).
// "Multimodality image registration by maximisation of mutual
// information". IEEE Transactions on Medical Imaging 16(2):187-198
//    H   = H.*log2(H./(s2*s1));
//    mi  = sum(H(:));
//    ecc = -2*mi/(sum(s1.*log2(s1))+sum(s2.*log2(s2)));
// mi =  sum_{i,j} [ P(i,j) log2 (P(i,j) / Pr(i) Pc(j)) ]
// mi = H(M) + H(N) - H(M,N)
// ecc = -2*mi / [sum_i(Pr(i) log2(Pr(i))) + sum_i(Pc(i) log2(Pc(i))) ]
// ecc = -2*mi / (H(M) + H(N))
// ecc = -2 [ 1 - H(M,N) / (H(M) + H(N))  ] 
// ecc = -2 + 2 H(M,N) / (H(M) + H(N))  
{
  //double mi = computeMI();
  double msum = 0.0;
  double mi = 0.0;
  double d;
  int i, j;
  //addeps(2.2204E-16);
  normalize();
  computeRCsums();
  for (i = 0; i < n; i++)
  {
    if (rowsum[i] != 0.0)
      msum += rowsum[i] * log2(rowsum[i]);
    if (colsum[i] != 0.0)
      msum += colsum[i] * log2(colsum[i]);
    for (j = 0; j < n; j++)
    {
      d = histo[i][j];
      if (d != 0.0)
        mi += d * log2(d / (rowsum[i] * colsum[j]));
    }
  }
  return -2 * mi / msum;
//  return (-2.0 + 2.0/computeNMI());
}

double JointHisto::computeNMI()
// compute Normalised Mutual Information
// Studholme,  Hill & Hawkes (1998).
// "A normalized entropy measure of 3-D medical image alignment".
// in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.             
// nmi = [sum_i(Pr(i) log2(Pr(i))) + sum_i(Pc(i) log2(Pc(i))) ] / sum_ij (P(i,j) log2(P(i,j)))
// (H(M) + H(N)) / H(M,N)
{
  double s1 = 0;
  double s2 = 0;
  int i, j;
  double d;
  if (sum == 0.0)
    return -100;
  //addeps(2.2204E-16);
  normalize();
  computeRCsums();

  for (i = 0; i < n; i++)
  {
    if (rowsum[i] != 0.0)
      s1 += rowsum[i] * log2(rowsum[i]);
    if (colsum[i] != 0.0)
      s1 += colsum[i] * log2(colsum[i]);
    for (j = 0; j < n; j++)
    {
      d = histo[i][j];
      if (d != 0.0)
        s2 += d * log2(d);
    }
  }

  return s1 / s2;
}

double JointHisto::computeNCC()
// compute Normalised Cross Correlation
// mean1 = sum_i Pr(i) * i
// mean2 = sum_i Pc(i) * i
// sig1  = sqrt(sum_i Pr(i) * (i-mean1)^2)
// sig2  = sqrt(sum_i Pc(i) * (i-mean2)^2)
// ncc   = (sum_ij P(i,j) (i-mean1) (j-mean2)) / (sig1 * sig2)
{
  double m1 = 0.0;
  double m2 = 0.0;
  int i, j;
  //addeps(2.2204E-16);
  //addeps(1);
  normalize();
  computeRCsums();
  for (i = 0; i < n; i++)
  {
    m1 += rowsum[i] * i;
    m2 += colsum[i] * i;
  }
  //cout << " m1 : " << m1 << "  m2: " << m2 << endl;

  double sig1 = 0.0;
  double sig2 = 0.0;
  for (i = 0; i < n; i++)
  {
    sig1 += rowsum[i] * (i - m1) * (i - m1);
    sig2 += colsum[i] * (i - m2) * (i - m2);
  }
  sig1 = sqrt(sig1);
  sig2 = sqrt(sig2);
  //cout << " sig1 : " << sig1 << "  sig2: " << sig2 << endl;

  double ncc = 0.0;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      ncc += histo[i][j] * (i - m1) * (j - m2);

  return ncc / (sig1 * sig2);

}

double JointHisto::computeSCR()
// symmetric correlation ration (CR(I1, I2) + CR(I2,I1))
// Roche, Malandain, Pennec, Ayache (1998)
// "The Correlation Ratio as a New Similarity Measure for
// Multimodal Image Registration"
// in Proc. MICCAI'98, vol 1496 of LNCS, 1115-1124
//
// CR(Ic,Ir) := 1 - (1/sig) sum_i (sig_i Pc(i))
// where sig = sum_i(i^2 Pr(i)) - [sum_i (i Pr(i))]^2
// and  sig_i = (1/Pc(i)) sum_j j^2 P(j,i) - [ 1/Pc(i) sum_j j P(j,i) ]^2
// or equiv:
// sig = sum_i ( Pr(i) (i-m)^2 )   with  m = sum_i Pr(i) i
// sig_i =  (1/Pc(i)) sum_j P(j,i) (j-m_i)^2  with m_i=sum_j P(j,i) j
{

  normalize();
  computeRCsums();
  // compute means
  double m1 = 0.0;
  double m2 = 0.0;
  int i, j;
  for (i = 0; i < n; i++)
  {
    m1 += rowsum[i] * i;
    m2 += colsum[i] * i;
  }
  // compute sig^2
  double sig1 = 0.0;
  double sig2 = 0.0;
  for (i = 0; i < n; i++)
  {
    sig1 += rowsum[i] * (i - m1) * (i - m1);
    sig2 += colsum[i] * (i - m2) * (i - m2);
  }

  double mi1, sigi1;
  double sum1 = 0.0;
  for (i = 0; i < n; i++)
  {
    if (colsum[i] == 0.0)
      continue;
    mi1 = 0.0;
    for (j = 0; j < n; j++)
      mi1 += j * histo[j][i];
    mi1 /= colsum[i];

    sigi1 = 0.0;
    for (j = 0; j < n; j++)
      sigi1 += histo[j][i] * (j - mi1) * (j - mi1);

    sum1 += sigi1;
  }
  sum1 /= sig1;
  sum1 = 1.0 - sum1;

  double mi2, sigi2;
  double sum2 = 0.0;
  for (i = 0; i < n; i++)
  {
    if (rowsum[i] == 0.0)
      continue;
    mi2 = 0.0;
    for (j = 0; j < n; j++)
      mi2 += j * histo[i][j];
    mi2 /= rowsum[i];

    sigi2 = 0.0;
    for (j = 0; j < n; j++)
      sigi2 += histo[i][j] * (j - mi2) * (j - mi2);

    sum2 += sigi2;
  }
  sum2 /= sig2;
  sum2 = 1.0 - sum2;

  return sum1 + sum2;

// // same but slower:
//   int i,j;
//   //addeps(1);
//   normalize();
//   computeRCsums();
//   
//   // do everything twice to make stuff symmetric
//   double s1t1 = 0.0; // sum_i (i^2 P(i))
//   double s1t2 = 0.0; // sum_i (i P(i))
//   double s2t1 = 0.0;
//   double s2t2 = 0.0;
//   double t;
//   for (i = 0; i<n; i++)
//   {
//      t = rowsum[i] * i;
//      s1t1 += t*i;
//      s1t2 += t;
//      
//      t = colsum[i] * i;
//      s2t1 += t*i;
//      s2t2 += t;
//   }
//   double sig1 = s1t1 - s1t2 * s1t2;
//   double sig2 = s2t1 - s2t2 * s2t2;
//   //cout << " sig1: " << sig1 << endl;
//   //cout << " sig2: " << sig2 << endl;
// 
// 
//   double siga1t1;
//   double siga1t2;
//   double sum1 = 0.0;
//   for (i = 0; i<n; i++)
//   {    
//     if (colsum[i]==0.0) continue;
//     siga1t1 = 0.0; siga1t2 = 0.0;
//     for (j = 0; j<n; j++)
//     {
//       t = j * histo[j][i];
//       siga1t1 += j * t;
//       siga1t2 += t;
//     }
//     sum1 += siga1t1 - siga1t2 * siga1t2 / colsum[i];
//   }
// 
//   double siga2t1;
//   double siga2t2;
//   double sum2 = 0.0;
//   for (i = 0; i<n; i++)
//   {    
//     if (rowsum[i] == 0.0) continue;
//     siga2t1 = 0.0; siga2t2 = 0.0;
//     for (j = 0; j<n; j++)
//     {
//       t = j * histo[i][j];
//       siga2t1 += j * t;
//       siga2t2 += t;
//     }
//     sum2 += siga2t1 - siga2t2 * siga2t2 / rowsum[i];
//   }
// //cout << " sum1: " << sum1 << endl;
// //cout << " sum2: " << sum2 << endl;
// 
//   return 2.0 - (sum1/sig1) - (sum2/sig2);
//   
}

double JointHisto::computeLS()
{
  double ssd = 0.0;
  int i, j, k;
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
    {
      if (i == j)
        continue;
      k = i - j;
      ssd += k * k * histo[i][j];
    }
  return ssd;
}
