//
// mnitest.cpp
//
// testing three kinds of mni xfm reading
//
//
#include <iostream>
#include <iomanip>

extern "C" {

#include "transform.h"

  char *Progname="mnitest";
}

char *file1="./bruce.xfm";
char *file2="./bruce2.xfm";
char *file3="./david.xfm";
char *file4="./tosa.xfm";

using namespace std;

int main(int argc, char *argv[])
{
  LTA *lta1, *lta2, *lta3, *lta4;
  cout << "////////////////////////////////////" << endl;
  cout << "Old way of reading -----------------" << endl;
  cout << "////////////////////////////////////" << endl;
  cout << "bruce.xfm" << endl;
  lta1 = LTAread(file1);
  LTAprint(stdout, lta1);
  cout << "bruce2.xfm" << endl;
  lta2 = LTAread(file2);
  LTAprint(stdout, lta2);
  cout << "david.xfm" << endl;
  lta3 = LTAread(file3);
  LTAprint(stdout, lta3);
  cout << "tosa.xfm" << endl;
  lta4 = LTAread(file4);
  LTAprint(stdout, lta4);
  LTAfree(&lta1);
  LTAfree(&lta2);
  LTAfree(&lta3);
  LTAfree(&lta4);

  cout << "////////////////////////////////////" << endl;
  cout << "New way of reading -----------------" << endl;
  cout << "////////////////////////////////////" << endl;
  cout << "bruce.xfm" << endl;
  lta1 = LTAreadEx(file1);
  LTAprint(stdout, lta1);
  cout << "bruce2.xfm" << endl;
  lta2 = LTAreadEx(file2);
  LTAprint(stdout, lta2);
  cout << "david.xfm" << endl;
  lta3 = LTAreadEx(file3);
  LTAprint(stdout, lta3);
  cout << "tosa.xfm" << endl;
  lta4 = LTAreadEx(file4);
  LTAprint(stdout, lta4);
  LTAfree(&lta1);
  LTAfree(&lta2);
  LTAfree(&lta3);
  LTAfree(&lta4);
}
