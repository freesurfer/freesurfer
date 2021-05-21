/**
 * @brief A programm to compute gradient information
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
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cassert>

#include "error.h"
#include "macros.h"
#include "mri.h"
#include "histo.h"
#include "matrix.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "version.h"
#include "transform.h"

using namespace std;
const char *Progname = NULL;

double scanX(MRI* mri_grad, MRI* mri_mask)
{
  if (mri_mask)
    assert(mri_grad->width == mri_mask->width);
  if (mri_mask)
    assert(mri_grad->height == mri_mask->height);
  if (mri_mask)
    assert(mri_grad->depth == mri_mask->depth);

  // find vert edges in gradx
  double thres = 5.0;
  int w, h, d;

  double lastg = -1;
  double currg = -1;
  double nextg;
  double edgewsum = 0.0;
  int edgecount = 0;
  int left = -1;
  int right = -1;
  double max = 0;
  for (d = 0; d < mri_grad->depth; d++)
    for (h = 0; h < mri_grad->height; h++)
//d=128;
//h=128;
    {
      //cout << " d: " << d << "  h: "<< h << endl;

      lastg = -1;
      currg = -1;
      left = -1;
      right = -1;
      max = 0;
      for (w = 0; w < mri_grad->width; w++)
      {
        if (mri_mask != NULL)
        {
          if (MRIgetVoxVal(mri_mask, w, h, d, 0) <= 0)
          {
            lastg = -1;
            currg = -1;
            left = -1;
            right = -1;
            max = 0;
            continue;
          }
        }
        if (lastg < 0)
        {
          lastg = fabs(MRIgetVoxVal(mri_grad, w, h, d, 0));
          continue;
        }

        if (currg < 0)
        {
          currg = fabs(MRIgetVoxVal(mri_grad, w, h, d, 0));
          if (currg > lastg)
            left = w - 1;
          continue;
        }

        nextg = fabs(MRIgetVoxVal(mri_grad, w, h, d, 0));
        //cout << " .. " << lastg << " " << currg << " " << nextg << "  left: " << left << " max: " << max << " right: " << right <<endl;
        // check direction change
        if (left < 0) // not on an edge
        {
          if (nextg > currg) // now going up, set left
          {
            left = w - 1;
            max = nextg; // reset max starting here
          }
          // continue in other cases (down or sideways)
        }
        else //we have a left minimum
        {

          if (nextg == currg && currg < max && right < 0) // going sideways after max passing max
            right = w - 1; // possible right boundary

          if (nextg > currg && currg < max) // going up after passing max
          {
            if (right < 0)
              right = w - 1; // else right was set when going sideways

            if (max >= thres)
            {
              // measure edge
              edgecount++;
              edgewsum += right - left;
              //cout << " measure edge max: " << max << " L: " << left << " R: " << right << " w: " << right-left<< endl;
            }
            left = w - 1;
            right = -1;
            max = nextg;
          }
          // other cases are 1. going up/sideways before reaching max or going down after max

          // finally set new max
          if (nextg > max)
            max = nextg; // store higher max
        }

        lastg = currg;
        currg = nextg;

      }
    }
  edgewsum /= edgecount;
  cout << " blur: " << edgewsum << " edges: " << edgecount << endl;
  return edgewsum;
}

double scanY(MRI* mri_grad, MRI* mri_mask)
{
  if (mri_mask)
    assert(mri_grad->width == mri_mask->width);
  if (mri_mask)
    assert(mri_grad->height == mri_mask->height);
  if (mri_mask)
    assert(mri_grad->depth == mri_mask->depth);

  // find vert edges in gradx
  double thres = 5.0;
  int w, h, d;

  double lastg = -1;
  double currg = -1;
  double nextg;
  double edgewsum = 0.0;
  int edgecount = 0;
  int left = -1;
  int right = -1;
  double max = 0;
  for (d = 0; d < mri_grad->depth; d++)
    for (w = 0; w < mri_grad->width; w++)
//d=128;
//h=128;
    {
      //cout << " d: " << d << "  h: "<< h << endl;

      lastg = -1;
      currg = -1;
      left = -1;
      right = -1;
      max = 0;
      for (h = 0; h < mri_grad->height; h++)
      {
        if (mri_mask != NULL)
        {
          if (MRIgetVoxVal(mri_mask, w, h, d, 0) <= 0)
          {
            lastg = -1;
            currg = -1;
            left = -1;
            right = -1;
            max = 0;
            continue;
          }
        }
        if (lastg < 0)
        {
          lastg = fabs(MRIgetVoxVal(mri_grad, w, h, d, 1));
          continue;
        }

        if (currg < 0)
        {
          currg = fabs(MRIgetVoxVal(mri_grad, w, h, d, 1));
          if (currg > lastg)
            left = h - 1;
          continue;
        }

        nextg = fabs(MRIgetVoxVal(mri_grad, w, h, d, 1));
        //cout << " .. " << lastg << " " << currg << " " << nextg << "  left: " << left << " max: " << max << " right: " << right <<endl;
        // check direction change
        if (left < 0) // not on an edge
        {
          if (nextg > currg) // now going up, set left
          {
            left = h - 1;
            max = nextg; // reset max starting here
          }
          // continue in other cases (down or sideways)
        }
        else //we have a left minimum
        {

          if (nextg == currg && currg < max && right < 0) // going sideways after max passing max
            right = h - 1; // possible right boundary

          if (nextg > currg && currg < max) // going up after passing max
          {
            if (right < 0)
              right = h - 1; // else right was set when going sideways

            if (max >= thres)
            {
              // measure edge
              edgecount++;
              edgewsum += right - left;
              //cout << " measure edge max: " << max << " L: " << left << " R: " << right << " w: " << right-left<< endl;
            }
            left = h - 1;
            right = -1;
            max = nextg;
          }
          // other cases are 1. going up/sideways before reaching max or going down after max

          // finally set new max
          if (nextg > max)
            max = nextg; // store higher max
        }

        lastg = currg;
        currg = nextg;

      }
    }
  edgewsum /= edgecount;
  cout << " blur: " << edgewsum << " edges: " << edgecount << endl;
  return edgecount;
}

double scanZ(MRI* mri_grad, MRI* mri_mask)
{
  if (mri_mask)
    assert(mri_grad->width == mri_mask->width);
  if (mri_mask)
    assert(mri_grad->height == mri_mask->height);
  if (mri_mask)
    assert(mri_grad->depth == mri_mask->depth);

  // find vert edges in gradx
  double thres = 5.0;
  int w, h, d;

  double lastg = -1;
  double currg = -1;
  double nextg;
  double edgewsum = 0.0;
  int edgecount = 0;
  int left = -1;
  int right = -1;
  double max = 0;
  for (w = 0; w < mri_grad->width; w++)
    for (h = 0; h < mri_grad->height; h++)
//d=128;
//h=128;
    {
      //cout << " d: " << d << "  h: "<< h << endl;

      lastg = -1;
      currg = -1;
      left = -1;
      right = -1;
      max = 0;
      for (d = 0; d < mri_grad->depth; d++)
      {
        if (mri_mask != NULL)
        {
          if (MRIgetVoxVal(mri_mask, w, h, d, 0) <= 0)
          {
            lastg = -1;
            currg = -1;
            left = -1;
            right = -1;
            max = 0;
            continue;
          }
        }
        if (lastg < 0)
        {
          lastg = fabs(MRIgetVoxVal(mri_grad, w, h, d, 2));
          continue;
        }

        if (currg < 0)
        {
          currg = fabs(MRIgetVoxVal(mri_grad, w, h, d, 2));
          if (currg > lastg)
            left = d - 1;
          continue;
        }

        nextg = fabs(MRIgetVoxVal(mri_grad, w, h, d, 2));
        //cout << " .. " << lastg << " " << currg << " " << nextg << "  left: " << left << " max: " << max << " right: " << right <<endl;
        // check direction change
        if (left < 0) // not on an edge
        {
          if (nextg > currg) // now going up, set left
          {
            left = d - 1;
            max = nextg; // reset max starting here
          }
          // continue in other cases (down or sideways)
        }
        else //we have a left minimum
        {

          if (nextg == currg && currg < max && right < 0) // going sideways after max passing max
            right = d - 1; // possible right boundary

          if (nextg > currg && currg < max) // going up after passing max
          {
            if (right < 0)
              right = d - 1; // else right was set when going sideways

            if (max >= thres)
            {
              // measure edge
              edgecount++;
              edgewsum += right - left;
              //cout << " measure edge max: " << max << " L: " << left << " R: " << right << " w: " << right-left<< endl;
            }
            left = d - 1;
            right = -1;
            max = nextg;
          }
          // other cases are 1. going up/sideways before reaching max or going down after max

          // finally set new max
          if (nextg > max)
            max = nextg; // store higher max
        }

        lastg = currg;
        currg = nextg;

      }
    }
  edgewsum /= edgecount;
  cout << " blur: " << edgewsum << " edges: " << edgecount << endl;
  return edgecount;
}

double scanXold(MRI * mri_in, MRI* mri_grad, MRI* mri_mask)
{
  assert(mri_in->width == mri_grad->width);
  if (mri_mask)
    assert(mri_in->width == mri_mask->width);
  assert(mri_in->height == mri_grad->height);
  if (mri_mask)
    assert(mri_in->height == mri_mask->height);
  assert(mri_in->depth == mri_grad->depth);
  if (mri_mask)
    assert(mri_in->depth == mri_mask->depth);

  // find vert edges in gradx
  double thres = 5.0;
  int w, h, d;
  int extremum = -1;
  bool onedge = false;
  double lasti = -1;
  double lastg = -1;
  double curri = -1;
  double currg = -1;
  double nexti, nextg;
  double edgewsum = 0.0;
  int edgecount = 0;
  for (d = 0; d < mri_in->depth; d++)
    for (h = 0; h < mri_in->height; h++)
//d=128;
//h=128;
    {
      extremum = -1;
      onedge = false;
      lasti = -1;
      curri = -1;
      for (w = 0; w < mri_in->width; w++)
      {
        if (mri_mask != NULL)
        {
          if (MRIgetVoxVal(mri_mask, w, h, d, 0) <= 0)
          {
            extremum = -1;
            onedge = false;
            lasti = -1;
            curri = -1;
            continue;
          }
        }
        if (lasti < 0)
        {
          lasti = MRIgetVoxVal(mri_in, w, h, d, 0);
          lastg = fabs(MRIgetVoxVal(mri_grad, w, h, d, 0));
          extremum = w;
          continue;
        }

        if (curri < 0)
        {
          curri = MRIgetVoxVal(mri_in, w, h, d, 0);
          currg = fabs(MRIgetVoxVal(mri_grad, w, h, d, 0));
          continue;
        }

        nexti = MRIgetVoxVal(mri_in, w, h, d, 0);
        nextg = fabs(MRIgetVoxVal(mri_grad, w, h, d, 0));

        if (onedge)
        {
          assert(extremum >= 0);
          if (currg > lastg && currg > nextg)
          {
            cout << " onedge and found edge again???" << endl;
            assert(1 == 2);
          }

          // found closing extremum at w-1
          if ((curri > lasti && curri > nexti)
              || (curri < lasti && curri < nexti))
          {
            edgecount++;
            edgewsum += w - 1 - extremum;
            extremum = w - 1;
            onedge = false;
          }
        }
        else // not on edge:
        {
          if (curri > lasti && curri > nexti)
            extremum = w - 1;
          if (curri < lasti && curri < nexti)
            extremum = w - 1;

          if (currg > lastg && currg > nextg && currg > thres)
            onedge = true;
        }

        lasti = curri;
        curri = nexti;
        lastg = currg;
        currg = nextg;

      }
    }
  edgewsum /= edgecount;
  //cout << " blur: " << edgewsum << " edges: " << edgecount << endl;
  return edgecount;
}

int main(int argc, char *argv[])
{

  Progname = argv[0];

  if (argc < 2)
  {
    cout << endl;
    cout << argv[0] << " image.mgz" << endl;
    cout << endl;
//    cout << "    norm-div  (=1)  divide final distance by this (e.g. step adjustment)" << endl;
//    cout << "    dist-type " << endl;
//    cout << "       1  (default) Rigid Transform Distance (||log(R)|| + ||T||)" << endl;
//    cout << "       2            Affine Transform Distance (RMS) " << endl;
//    cout << "       3            8-corners mean distance after transform " << endl;
//    cout << "    invert1         1 true, 0 false (default)" << endl;
//    cout << endl;
    exit(1);
  }
  string mrif = argv[1];

  MRI* mri_in = MRIread(mrif.c_str());


  MRI* mri_mag = MRIalloc(mri_in->width, mri_in->height, mri_in->depth,
      MRI_FLOAT);
  MRI* mri_grad = MRIsobel(mri_in, NULL, mri_mag);

//  MRIwriteFrame(mri_grad,"sobel_grad1.mgz",0);
//  MRIwriteFrame(mri_grad,"sobel_grad2.mgz",1);
//  MRIwriteFrame(mri_grad,"sobel_grad3.mgz",2);
//  MRIwrite(mri_mag,"sobel_mag.mgz");

  int dd, hh, ww;
  double avg = 0.0;
  float val;
  int count = 0;
  int all = mri_mag->depth * mri_mag->height * mri_mag->width;
  for (dd = 0; dd < mri_mag->depth; dd++)
    for (hh = 0; hh < mri_mag->height; hh++)
      for (ww = 0; ww < mri_mag->width; ww++)
      {
        val = MRIgetVoxVal(mri_mag, ww, hh, dd, 0);
        if (val > 30)
        {
          count++;
          avg += val;
        }

        //   avg += tanh(val/10.0);
      }
  avg /= count;
  cout << "avg mag : " << avg << endl;
  cout << "count   : " << count << endl;

//  double b1 =  
//scanX(mri_grad,mri_mask);

//  double b2 =  scanY(mri_grad,mri_mask);
//  double b3 = scanZ(mri_grad,mri_mask);
//scanX(mri_mag,mri_mask);

//  double bmax = b1;
//  if (b2 > bmax) b2=bmax;
//  if (b3 > bmax) b3=bmax;

//  cout << "avg: " << (b1+b2+b3)/3.0 <<" max: " << bmax << endl;
//  scanX(mri_in,mri_grad,NULL);

//  exit(0);

  MRIfree(&mri_grad);

  MRIsetVoxVal(mri_mag, 0, 0, 0, 0, 77);
  int n = 2;
  HISTOGRAM *h = MRIhistogram(mri_mag, n);

  cout << " h = [ ";
  for (int ii = 0; ii < n; ii++)
  {
//    cout << "i " << ii << " : " << h->counts[ii] <<  " uval: " << h->bins[ii] << endl;
    cout << h->bins[ii] << " , " << h->counts[ii];
    if (ii < n - 1)
      cout << " ; " << endl;
    else
      cout << "] ; " << endl;
  }

  cout << "pct: " << 100 * h->counts[n - 1] / all << endl;

  cout << " entropy: " << HISTOgetEntropy(h) << endl;

}
