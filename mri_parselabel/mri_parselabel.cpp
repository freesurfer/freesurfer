//
// mri_parselabel.cpp
//
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>


extern "C" {
#include "mri.h"
#include "version.h"
  char *Progname;
}

using namespace std;

int useRealRAS = 0;

class Vertex
{
public:
  Vertex(double x, double y, double z, double value) 
    : x_(x), y_(y), z_(z), value_(value) {}

  double x_; 
  double y_; 
  double z_;
  double value_;
};

static int get_option(int argc, char *argv[]) ;

int main(int argc, char *argv[])
{
  int ac, nargs;
  char **av;
  Progname=argv[0];

  nargs = handle_version_option (argc, argv, "$Id: mri_parselabel.cpp,v 1.1 2004/06/01 20:14:05 tosa Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  if (argc < 5)
  {
    cerr << "Usage: mri_parselabel [option] <labelfile> <positionscaling> <volfileforlabel> <outputvol> <greyforlabel>" << endl;
    cerr << "option: -cras : use scanner ras value for label position" << endl;
    cerr << "                default is to use the surface ras value."<< endl;
    return -1;
  }
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  cout << "---------------------------------------------------------------" << endl;
  cout << "Inputs: " << endl;
  cout << "label file : " << argv[1] << endl;
  cout << "scale      : " << argv[2] << endl;
  cout << "label vol  : " << argv[3] << endl;
  cout << "output vol : " << argv[4] << endl;
  cout << "greyvalue  : " << argv[5] << endl;
  cout << "---------------------------------------------------------------" << endl;

  int val = (int) atof(argv[5]);
  if (val < 0 || val > 255)
  {
    cout << "grey value must be between 0 and 255.  you gave " << argv[5] << endl;
    return -1;
  }
  ifstream flabel(argv[1], ios::in);
  if (!flabel.good())
  {
    cerr << "Could not open file for label: " << argv[1] << endl;
    return -1;
  }  
  // get the scaling
  double scale = atof(argv[2]);

  // read all vertex positions
  vector<Vertex> vertices;
  char buf[1024];
  flabel.getline(buf, sizeof(buf));
  cout << "FileInfo: " << buf << endl;
  int numvertices;
  flabel >> numvertices;
  cout << "number of vertices are: " << numvertices << endl;
  while (flabel.good())
  {
    int num;
    double x, y, z;
    double value;
    flabel >> num >> x >> y >> z >> value;
    // cout << "(" << x << ", " << y << ", " << z << ")" << endl;
    Vertex v(scale*x, scale*y, scale*z, value);
    vertices.push_back(v);
  }
  flabel.close();

  MRI *mriIn=MRIread(argv[3]);
  if (!mriIn)
  {
    cerr << "Could not open input volume : " << argv[3] << endl;
    return -1;
  }
  MRI *mriOut = MRIalloc(mriIn->width, mriIn->height, mriIn->depth, MRI_UCHAR);
  MRIcopyHeader(mriIn, mriOut);
  Real xv, yv, zv;
  double value;
  for (size_t i=0; i < vertices.size(); ++i)
  {
    switch(useRealRAS)
    {
    case 0:
      MRIsurfaceRASToVoxel(mriIn, vertices[i].x_, vertices[i].y_, vertices[i].z_,
			   &xv, &yv, &zv);
      break;
    case 1:
      MRIworldToVoxel(mriIn, vertices[i].x_, vertices[i].y_, vertices[i].z_,
			   &xv, &yv, &zv);
      break;
    default:
      cerr << "this should not happen.  useRealRAS is " << useRealRAS << endl;
      return -1;
    }
    // verify the voxel position
    if (xv < 0 || yv < 0 || zv < 0 
	|| xv > mriIn->width-1 || yv > mriIn->height-1 || zv > mriIn->depth -1)
    {
      cerr << "Vertex " << i 
	   << " has invalid voxel position (" << xv << ", " << yv << ", " << zv << ") from "
	   << " (" << vertices[i].x_ << ", " << vertices[i].y_ << ", " << vertices[i].z_ 
	   << ") " << endl;

      MRIfree(&mriIn);
      MRIfree(&mriOut);
  
      return -1;
    }

    switch(mriIn->type)
    {
    case MRI_UCHAR:
      value = MRIvox(mriIn, nint(xv), nint(yv), nint(zv)); break;
    case MRI_INT:
      value = MRIIvox(mriIn, nint(xv), nint(yv), nint(zv)); break;
    case MRI_FLOAT:
      value = MRIFvox(mriIn, nint(xv), nint(yv), nint(zv)); break;
    case MRI_SHORT:
      value = MRISvox(mriIn, nint(xv), nint(yv), nint(zv)); break;
    case MRI_LONG:
      value = MRILvox(mriIn, nint(xv), nint(yv), nint(zv)); break;
    }
    if (!FZERO(val - vertices[i].value_))
    {
      cout << "Grey scale does not match: " << endl;
      cout << "label position : (" << vertices[i].x_ << ", "
	   << vertices[i].y_ <<", " << vertices[i].z_ << ") with value " << vertices[i].value_ << endl;
      cout << "voxel position : (" << nint(xv) << ", " << nint(yv) << ", " << nint(zv)
	   << ")  has the value " << value << endl;
    }
    MRIvox(mriOut, nint(xv), nint(yv), nint(zv)) = (unsigned char) val;
  }
  MRIwrite(mriOut, argv[4]);

  MRIfree(&mriIn);
  MRIfree(&mriOut);
  
  return 0;
}

int get_option(int argc, char *argv[])
{
  char *option ;
  int nargs;

  option=argv[1]+1; // pass -
  if (strcmp(option, "cras")==0)
  {
    cout << "use the scanner RAS coordinates for label position" << endl;
    useRealRAS=1;
  }
  nargs=0;
  return (nargs);
}
    
  
