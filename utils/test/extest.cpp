/**
 * @brief test fread*Ex() routines
 *
 */
/*
 * Original Author: Y. Tosa
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

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cstdlib>
#include <typeinfo>

extern "C"
{
#include "fio.h"
}

const int MAX = 5;
const double MAXD=1.7976931348623157e+308;
const double MIND=2.2250738585072014e-308;
const char* Progname = "extest";

using namespace std;

template <typename T> bool checkOK(const T &p, const T &q)
{
  if (p != q)
  {
    cerr << "type is " << typeid(p).name() << ":";
    cerr << "value differed " << p << " vs. " << q << endl;
    return false;
  }
  else
  {
    cout << "type is " << typeid(p).name() << ":";
    cout << "read " << p << "  orig: " << q << endl;
    return true;
  }
}

template <> bool checkOK(const double &p, const double &q)
{
  if ((p-q) > 0.00000000001)
  {
    cerr << "type is " << typeid(p).name() << ":";
    cerr << "value differed " << p << " vs. " << q << endl;
    return false;
  }
  else
  {
    cout << "type is " << typeid(p).name() << ":";
    cout << "read " << p << "  orig: " << q << endl;
    return true;
  }
}

template <> bool checkOK(const float &p, const float &q)
{
  if ((p-q) > 0.0000001)
  {
    cerr << "type is " << typeid(p).name() << ":";
    cerr << "value differed " << p << " vs. " << q << endl;
    return false;
  }
  else
  {
    cout << "type is " << typeid(p).name() << ":";
    cout << "read " << p << "  orig: " << q << endl;
    return true;
  }
}

int main(int argc, char *argv[])
{
  FILE *in, *out;
  int i[MAX];
  short s[MAX];
  float f[MAX];
  double d[MAX];

  int ii[MAX];
  short ss[MAX];
  float ff[MAX];
  double dd[MAX];

  const char *file = "./extest.out";

  out = fopen(file, "w");
  if (!out)
  {
    cerr << "could not open " << file << endl;
    exit (-1);
  }
  i[0] = 1;
  i[1] = 256;
  i[2] = 256*256*256;
  s[0] = 1;
  s[1] = (short)(256*256);
  s[2] = (short)(i[2]-1);
  f[0] = 1.;
  f[1] = 1.-0.001;
  f[2] = 1./3.;
  d[0] = 1.;
  d[1] = 1.-0.001;
  d[2] = 1./3.;
  d[3] = MAXD;
  d[4] = MIND;

  fwriteInt(i[0], out);
  fwriteShort(s[0], out);
  fwriteFloat(f[0], out);
  fwriteDouble(d[0], out);
  fwriteDouble(d[1], out);
  fwriteFloat(f[1], out);
  fwriteShort(s[1], out);
  fwriteInt(i[1], out);
  fwriteDouble(d[2], out);
  fwriteFloat(f[2], out);
  fwriteShort(s[2], out);
  fwriteInt(i[2], out);
  fwriteDouble(d[3], out);
  fwriteDouble(d[4], out);

  fclose(out);

  in = fopen(file, "rb");
  if (freadIntEx(&ii[0], in))
  {
    if (!checkOK(i[0], ii[0]))
      goto finalize;
  }
  if (freadShortEx(&ss[0], in))
  {
    if (!checkOK(s[0], ss[0]))
      goto finalize;
  }
  if (freadFloatEx(&ff[0], in))
  {
    if (!checkOK(f[0], ff[0]))
      goto finalize;
  }
  if (freadDoubleEx(&dd[0], in))
  {
    if (!checkOK(d[0], dd[0]))
      goto finalize;
  }
  if (freadDoubleEx(&dd[1], in))
  {
    if (!checkOK(d[1], dd[1]))
      goto finalize;
  }
  if (freadFloatEx(&ff[1], in))
  {
    if (!checkOK(f[1], ff[1]))
      goto finalize;
  }
  if (freadShortEx(&ss[1], in))
  {
    if (!checkOK(s[1], ss[1]))
      goto finalize;
  }
  if (freadIntEx(&ii[1], in))
  {
    if (!checkOK(i[1], ii[1]))
      goto finalize;
  }
  if (freadDoubleEx(&dd[2], in))
  {
    if (!checkOK(d[2], dd[2]))
      goto finalize;
  }
  if (freadFloatEx(&ff[2], in))
  {
    if (!checkOK(f[2], ff[2]))
      goto finalize;
  }
  if (freadShortEx(&ss[2], in))
  {
    if (!checkOK(s[2], ss[2]))
      goto finalize;
  }
  if (freadIntEx(&ii[2], in))
  {
    if (!checkOK(i[2], ii[2]))
      goto finalize;
  }
  if (freadDoubleEx(&dd[3], in))
  {
    if (!checkOK(d[3], dd[3]))
      goto finalize;
  }
  if (freadDoubleEx(&dd[4], in))
  {
    if (!checkOK(d[4], dd[4]))
      goto finalize;
  }

  fclose(in);
  cout << "No error" << endl;
  return 0;

finalize:
  cerr << "Error" << endl;
  exit(-1);
}
