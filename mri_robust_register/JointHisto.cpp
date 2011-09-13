/**
 * @file JointHisto.cpp
 * @brief A class for a joint histogram of two images
 *
 */

/*
 * Original Author: Martin Reuter
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2011/09/13 03:08:25 $
 *    $Revision: 1.1 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <vcl_iostream.h>
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
  assert (width == mri2->width);
  assert (height == mri2->height);
  assert (depth == mri2->depth);
  assert (mri1->type == MRI_UCHAR);
  assert (mri2->type == MRI_UCHAR);
  assert (mri1->xsize == mri1->ysize);
  assert (mri1->xsize == mri1->zsize);
  assert (mri1->xsize == mri2->xsize);
  assert (mri1->xsize == mri2->ysize);
  assert (mri1->xsize == mri2->zsize);

  histo.fill(0.0);
  haseps = false;

  int x,y,z;
  int isum = 0;
  for (z=0;z<=depth - d3; z+=d3)
  for (y=0;y<=height- d2; y+=d2)
  for (x=0;x<=width - d1; x+=d1)
  {
    histo[(int)MRIvox(mri1,x,y,z)][(int)MRIvox(mri2,x,y,z)] += 1;
    isum++; 
  }
  sum = isum;
  //cout <<"Sum: " << sum << endl;
  
}

float samp(MRI * mri, float x, float y, float z)
{
//cout << " samp: " << x << " " << y <<" " << z << endl;
	int ix, iy, iz;
	float dx1, dy1, dz1, dx2, dy2, dz2;
	float k111,k112,k121,k122,k211,k212,k221,k222;
	float vf;
	//unsigned char *ff;

	ix = (int)floor(x); dx1=x-ix; dx2=1.0-dx1;
	iy = (int)floor(y); dy1=y-iy; dy2=1.0-dy1;
	iz = (int)floor(z); dz1=z-iz; dz2=1.0-dz1;

  k222 = MRIgetVoxVal(mri,ix,iy,iz,0);
  k122 = MRIgetVoxVal(mri,ix+1,iy,iz,0);
  k212 = MRIgetVoxVal(mri,ix,iy+1,iz,0);
  k112 = MRIgetVoxVal(mri,ix+1,iy+1,iz,0);
  k221 = MRIgetVoxVal(mri,ix,iy,iz+1,0);
  k121 = MRIgetVoxVal(mri,ix+1,iy,iz+1,0);
  k211 = MRIgetVoxVal(mri,ix,iy+1,iz+1,0);
  k111 = MRIgetVoxVal(mri,ix+1,iy+1,iz+1,0);

	vf = (((k222*dx2+k122*dx1)*dy2       +
	       (k212*dx2+k112*dx1)*dy1))*dz2 +
	     (((k221*dx2+k121*dx1)*dy2       +
	       (k211*dx2+k111*dx1)*dy1))*dz1;
	return(vf);
}

void JointHisto::create(MRI *mriS, MRI* mriT,
           const vnl_matrix_fixed < double, 4,4 >& Msi,
           const vnl_matrix_fixed < double, 4,4 >& Mti,
           int d1, int d2, int d3)
{
  assert (mriS->type == MRI_UCHAR);
  assert (mriT->type == MRI_UCHAR);

  //double eps = 2.2204E-16;
  //histo.fill(eps);
  histo.fill(0);
  haseps = false;
    
//   static float ran[97] = {0.656619,0.891183,0.488144,0.992646,0.373326,0.531378,0.181316,0.501944,0.422195,
// 	                        0.660427,0.673653,0.95733,0.191866,0.111216,0.565054,0.969166,0.0237439,0.870216,
// 	                        0.0268766,0.519529,0.192291,0.715689,0.250673,0.933865,0.137189,0.521622,0.895202,
// 	                        0.942387,0.335083,0.437364,0.471156,0.14931,0.135864,0.532498,0.725789,0.398703,
// 	                        0.358419,0.285279,0.868635,0.626413,0.241172,0.978082,0.640501,0.229849,0.681335,
// 	                        0.665823,0.134718,0.0224933,0.262199,0.116515,0.0693182,0.85293,0.180331,0.0324186,
// 	                        0.733926,0.536517,0.27603,0.368458,0.0128863,0.889206,0.866021,0.254247,0.569481,
// 	                        0.159265,0.594364,0.3311,0.658613,0.863634,0.567623,0.980481,0.791832,0.152594,
// 	                        0.833027,0.191863,0.638987,0.669,0.772088,0.379818,0.441585,0.48306,0.608106,
// 	                        0.175996,0.00202556,0.790224,0.513609,0.213229,0.10345,0.157337,0.407515,0.407757,
// 	                        0.0526927,0.941815,0.149972,0.384374,0.311059,0.168534,0.896648};
// 	int    iran=0;
	double x,y,z;
	double rx, ry, rz;
  double xs,ys,zs;
  double xt,yt,zt;
	double vs,vt;
	int    ivs, ivt;
  double sdiff, tdiff;
  double dt[3] = { mriT->width, mriT->height, mriT->depth};
  double ds[3] = { mriS->width, mriS->height, mriS->depth};
  //int count = 0;
  //cout <<" df " << df[0] << " " << df[1] << " " << df[2] << endl;
	for(z=0.0; z<dt[2]-d3-1; z+=d3)
	{
		for(y=0.0; y<dt[1]-d2-1; y+=d2)
		{
			for(x=0.0; x<dt[0]-d1-1; x+=d1)
			{

				//rx  = x + ran[iran]*d1;
        //iran = (iran+1)%97;
				//ry  = y + ran[iran]*d2;
        //iran = (iran+1)%97;
				//rz  = z + ran[iran]*d3;
        //iran = (iran+1)%97;
        
        rx = x; ry = y; rz = z;

				xt  = Mti[0][0]*rx + Mti[0][1]*ry + Mti[0][2]*rz + Mti[0][3];
				yt  = Mti[1][0]*rx + Mti[1][1]*ry + Mti[1][2]*rz + Mti[1][3];
				zt  = Mti[2][0]*rx + Mti[2][1]*ry + Mti[2][2]*rz + Mti[2][3];

				xs  = Msi[0][0]*rx + Msi[0][1]*ry + Msi[0][2]*rz + Msi[0][3];
				ys  = Msi[1][0]*rx + Msi[1][1]*ry + Msi[1][2]*rz + Msi[1][3];
				zs  = Msi[2][0]*rx + Msi[2][1]*ry + Msi[2][2]*rz + Msi[2][3];
        //cout << "( " << x << " " << y << " " << z << " )  ( " << xp << " " << yp << " " << zp << " )" << endl;
				if (zs>=0.0 && zs<ds[2]-1 && ys>=0.0 && ys<ds[1]-1 && xs>=0.0 && xs<ds[0]-1 &&
            zt>=0.0 && zt<dt[2]-1 && yt>=0.0 && yt<dt[1]-1 && xt>=0.0 && xt<dt[0]-1 )
				{
					//vf  = samp(f, xp,yp,zp);
          MRIsampleVolumeFrame(mriS,xs,ys,zs,0,&vs);
					ivs = (int)floor(vs);
          //vg = samp(g,rx,ry,rz);
          MRIsampleVolumeFrame(mriT,xt,yt,zt,0,&vt);
					//ivg = (int)floor(vg+0.5);
					ivt = (int)floor(vt);
          //cout <<" ivf : " << ivf << " ivg: " << ivg << endl;
          assert (ivs >=0);
          assert (ivt >=0);
          assert (ivs <256);
          assert (ivt <256);
          sdiff  = vs-ivs;
          tdiff  = vt-ivt;
          // distribute peak among bins symetrically:
					histo[ivs][ivt] += (1-sdiff)*(1-tdiff);
					if (ivs<255)
						histo[ivs+1][ivt] += sdiff*(1-tdiff);
					if (ivt<255)
						histo[ivs][ivt+1] += (1-sdiff)*tdiff;
					if (ivs<255 && ivt<255)
						histo[ivs+1][ivt+1] += sdiff*tdiff;
          //count++;
				}
			}
		}
	}  
  //cout << " count: " << count << endl;
}

void JointHisto::set(const  vnl_matrix < double > & h)
// set histo matrix directly (for debugging)
{
  assert(h.rows() ==h.cols());
  histo = h;
  n=h.rows();
  haseps = false;
  sum = 0.0;
  double sum1 = 0.0;
  int i,j;
  for (i=0;i<n;i++)
  {
    sum1 = 0.0;
    for (j=0;j<n;j++)
    {
      sum1 += histo[i][j];
    }
    sum += sum1;
  }
}

void JointHisto::smooth(double fwhm)
{
  // compute sigma from fwhm:
  double sm = fwhm/sqrt(8*log(2));
  if (sm < 0.001) sm = 0.001;
  int t = (int)floor(3*sm+0.5); // rounding 3*sm
  // create kernel
  double  filter[2*t+1];
  double fsum = 0.0;
  double sm2 = -sm*sm;
  //cout << " sm = " << sm << endl;
  for (int i = -t; i<=t; i++)
  {
    filter[i+t] = exp((i*i)/sm2);
    fsum += filter[i+t];
  }
  // normalize
  for (int i=0;i<=2*t;i++)
  {
    filter[i] /= fsum;
    //cout << setprecision(16)<< " f[" << i << "]= " << filter[i] << endl;
  }

  // convolve 2D
  vnl_matrix < double > htmp(n,n);
  double dtmp;
  int i,j,k,count;
  for (i = 0;i<n;i++)
  for (j = 0;j<n;j++)
  {
    dtmp = 0;
    count = -1;
    for (k=j-t;k<=j+t;k++)
    {
      count++;
      if (k<0 || k>=n) continue;
      dtmp += histo[i][k] * filter[count];
    }
    htmp[i][j] = dtmp;
  }
  sum = 0.0;
  for (i = 0;i<n;i++)
  for (j = 0;j<n;j++)
  {
    dtmp = 0;
    count = -1;
    for (k=i-t;k<=i+t;k++)
    {
      count++;
      if (k<0 || k>=n) continue;
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
//		H   = H.*log2(H./(s2*s1));
//		mi  = sum(H(:));
{
  int i,j;
  double mi =0.0;
  double d;
  addeps(2.2204E-16);
  normalize();
  computeRCsums();
  for (i=0;i<n;i++)
  for (j=0;j<n;j++)
  {
    d = histo[i][j];
    mi += d * log2(d / (rowsum[i]*colsum[j]));
  }
  return mi;
}

double JointHisto::computeECC()
// compute Entropy Correlation Coefficient
// F Maes, A Collignon, D Vandermeulen, G Marchal & P Suetens (1997).
// "Multimodality image registration by maximisation of mutual
// information". IEEE Transactions on Medical Imaging 16(2):187-198
//		H   = H.*log2(H./(s2*s1));
//		mi  = sum(H(:));
//		ecc = -2*mi/(sum(s1.*log2(s1))+sum(s2.*log2(s2)));
{
  double mi = computeMI();
  double ecc = 0;
  addeps(2.2204E-16);
  normalize();
  computeRCsums();
  for (int i=0;i<n;i++)
    ecc += rowsum[i]*log2(rowsum[i])+colsum[i]*log2(colsum[i]);
 
  return -2*mi/ecc;
}

double JointHisto::computeNMI()
// compute Normalised Mutual Information
// Studholme,  Hill & Hawkes (1998).
// "A normalized entropy measure of 3-D medical image alignment".
// in Proc. Medical Imaging 1998, vol. 3338, San Diego, CA, pp. 132-143.             
//		nmi = (sum(s1.*log2(s1))+sum(s2.*log2(s2)))/sum(sum(H.*log2(H)));
{
  double s1 = 0;
  double s2 = 0;
  int i,j;
  double d;
  addeps(2.2204E-16);
  normalize();
  computeRCsums();
    
  for (i = 0;i<n;i++)
  {
    s1 += rowsum[i]*log2(rowsum[i])+colsum[i]*log2(colsum[i]);
    for (j=0;j<n;j++)
    {
      d = histo[i][j];
      s2 += d * log2(d);
    }
  }

  return  s1 / s2;
}

double JointHisto::computeNCC()
// compute Normalised Cross Correlation
//		i     = 1:size(H,1);
//		j     = 1:size(H,2);
//		m1    = sum(s2.*i');
//		m2    = sum(s1.*j);
//		sig1  = sqrt(sum(s2.*(i'-m1).^2));
//		sig2  = sqrt(sum(s1.*(j -m2).^2));
//		[i,j] = ndgrid(i-m1,j-m2);
//		ncc   = sum(sum(H.*i.*j))/(sig1*sig2);
{
cout << " UNTESTED " << endl;
  double m1 = 0;
  double m2 = 0;
  int i,j;
  addeps(2.2204E-16);
  normalize();
  computeRCsums();
  for (i = 0;i<n;i++)
  {
     m1 += rowsum[i] * i;
     m2 += colsum[i] * i;
  }
  //cout << " m1 : " << m1 << "  m2: " << m2 << endl;
  
  double sig1 = 0;
  double sig2 = 0;
  for (i = 0;i<n;i++)
  {
     sig1 += rowsum[i] * (i-m1) * (i-m1);
     sig2 += colsum[i] * (i-m2) * (i-m2);
  }
  sig1 = sqrt(sig1);
  sig2 = sqrt(sig2);
  //cout << " sig1 : " << sig1 << "  sig2: " << sig2 << endl;
  
  double ncc = 0;
  for (i = 0;i<n;i++)
  for (j = 0;j<n;j++)
     ncc += histo[i][j] * (i-m1) * (j-m2);
  
  return ncc/(sig1*sig2);

}

double JointHisto::computeLS()
{
  double ssd = 0.0;
  int i,j,k;
  for (i = 0;i<n;i++)
  for (j = 0;j<n;j++)
  {
    if (i==j) continue;
    k=i-j;
    ssd += k*histo[i][j];
  }
  return ssd;
}
