//
// gcaread
//

#include <iostream>

extern "C" {
#include "mri.h"
#include "gca.h"

  char *Progname = "gcaread";
}

using namespace std;

int main(int argc, char *argv[])
{
  if (argc < 2)
  {
    cout << "Usage: gcaread <gcafile>" << endl;
    return -1;
  }
  GCA *gca = GCAread(argv[1]);
  if (gca == 0)
  {
    cout << "could not open file " << argv[1] << endl;
  } 
  else
    GCAfree(&gca);
}
