//
// mri_parselabel.cpp
//
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>

extern "C" {
#include "mri.h"
#include "version.h"
#include "macros.h"
  char *Progname;
}

using namespace std;

int useRealRAS = 0;
int fillup=0;
double scale = 1;

// voxel size
double vx;
double vy;
double vz;

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


bool operator==(const Vertex &a, const Vertex &b)
{
  // don't compare values
  if ((fabs(a.x_ - b.x_) < vx/2.)
      && (fabs(a.y_ - b.y_) < vy/2.)
      && (fabs(a.z_ - b.z_) < vz/2.))
    return true;
  else
    return false;
}

// used for sort 
class LessX
{
public:
  bool operator() (const Vertex &a, const Vertex &b)
  {
    return ((a.x_ - b.x_) < 0);
  }
};

static int get_option(int argc, char *argv[]) ;

int (*voxToRAS)(MRI *mri, Real x, Real y, Real z, Real *xs, Real *ys, Real *zs);
int (*rasToVox)(MRI *mri, Real x, Real y, Real z, Real *xs, Real *ys, Real *zs);

int main(int argc, char *argv[])
{
  int nargs;
  Progname=argv[0];

  nargs = handle_version_option (argc, argv, "$Id: mri_parselabel.cpp,v 1.4 2004/06/02 21:30:51 tosa Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  // option parse
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  // check 
  if (argc < 5)
  {
    cerr << "Usage: mri_parselabel [option] <labelfile> <positionscaling> <volfileforlabel> <outputvol> <greyforlabel>" << endl;
    cerr << "option: -cras        : use scanner ras value for label position" << endl;
    cerr << "                       default is to use the surface ras value."<< endl;
    cerr << "        -scale <val> : uses <val> to scale label position" << endl;
    cerr << "        -fillup      : try to verify none of the label positions are missed." << endl;
    cerr << "                       takes a long time.....                               " << endl;
    return -1;
  }
  cout << "---------------------------------------------------------------" << endl;
  cout << "Inputs: " << endl;
  cout << "label file : " << argv[1] << endl;
  cout << "label vol  : " << argv[2] << endl;
  cout << "output vol : " << argv[3] << endl;
  cout << "greyvalue  : " << argv[4] << endl;
  cout << "scale      : " << scale << endl;
  cout << "---------------------------------------------------------------" << endl;

  int val = (int) atof(argv[4]);
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

  MRI *mriIn=MRIread(argv[2]);
  if (!mriIn)
  {
    cerr << "Could not open input volume : " << argv[2] << endl;
    return -1;
  }
  vx = mriIn->xsize;
  vy = mriIn->ysize;
  vz = mriIn->zsize;

  // read all vertex positions
  vector<Vertex> vertices;
  char buf[1024];
  flabel.getline(buf, sizeof(buf));
  cout << "FileInfo: " << buf << endl;
  int numvertices;
  flabel >> numvertices;
  cout << "number of vertices : " << numvertices << endl;
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

  // sort
  cout << "sorting ..." << endl;
  sort(vertices.begin(), vertices.end(), LessX());
  cout << "sorting done" << endl;
  
  MRI *mriOut = MRIalloc(mriIn->width, mriIn->height, mriIn->depth, MRI_UCHAR);
  MRIcopyHeader(mriIn, mriOut);

  // set voxToRAS and rasToVox transform function
  Real xv, yv, zv;
  if (useRealRAS)
  {
    voxToRAS = MRIvoxelToWorld;
    rasToVox = MRIworldToVoxel;
  }
  else
  {
    voxToRAS = MRIvoxelToSurfaceRAS;
    rasToVox = MRIsurfaceRASToVoxel;
  }

  cout << "naive filling...." << endl;
  // first use the naive approach
  for (size_t i=0; i < vertices.size(); ++i)
  {
    rasToVox(mriIn, vertices[i].x_, vertices[i].y_, vertices[i].z_, &xv, &yv, &zv);
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
    MRIvox(mriOut, nint(xv), nint(yv), nint(zv)) = (unsigned char) val;
  }
  cout << "naive filling done" << endl;
  MRIfree(&mriIn);
  
  if (!fillup)
    goto here;

  cout << "check missed positions...." << endl;
  // now check whether missing points exist or not
  {
    unsigned char curval =0;
    int count =0;
    int addedcount=0;
    int xi, yi, zi;
    for (int z = 0; z < mriOut->depth; z++)
    {
      for (int y = 0; y < mriOut->height; y++)
      {
	for (int x = 0; x < mriOut->width; x++)
	{
	  curval = MRIvox(mriOut, x, y, z);
	  // this point is marked, then looking around to see other points
	  // may belong to the label list.
	  if (curval == val)
	  {
	    // look around 3 x 3 x 3 neighbors
	    /////////////////////////////////////////////
	    for (int z0 = -1 ; z0 <= 1 ; z0++)
	    {
	      zi = mriOut->zi[z+z0] ;
	      for (int y0 = -1 ; y0 <= 1 ; y0++)
	      {
		yi = mriOut->yi[y+y0] ;
		for (int x0 = -1 ; x0 <= 1 ; x0++)
		{
		  xi = mriOut->xi[x+x0] ;
		  Real xl, yl, zl;
		  voxToRAS(mriOut, xi, yi, zi, &xl, &yl, &zl);
		  // search
		  for (size_t i=0; i < vertices.size(); ++i)
		  {
		    // since I sorted by X, the vertex position > current position + 2x vx
		    // don't do anymore
		    if ((vertices[i].x_) > (xl + 2*vx))
		      goto nextpoint;

		    if (vertices[i] == Vertex(xl, yl, zl, 0))
		    {
		      if (xi != x || yi != y || zi != z)
		      {
			// if it is not set
			if (MRIvox(mriOut, xi, yi, zi) != val)
			{
			  MRIvox(mriOut, xi, yi, zi) = val;
			  cout << "\nadded another point: ( " 
			       << xi << ", " << yi << ", " << zi <<")" << endl; 
			  addedcount++;
			}
		      }
		      // same point
		      else if (xi == x && yi == y && zi == z)
		      {
			count++;
			if (count%1000==0)
			{
			  cout << ".";
			  cout.flush();
			}
		      }
		    }		    
		  }// for i
		nextpoint:
		  continue;
		}
	      }
	    } /////////////////////////////////////////////////////////////////
	  } // x
	} // y
      } // z
    }
    cout << endl;
    cout << "total voxels = " << count + addedcount << endl;
    cout << "     initial = " << count << endl;
    cout << "       added = " << addedcount << endl;
  }
 here:

  MRIwrite(mriOut, argv[3]);

  MRIfree(&mriOut);
  
  return 0;
}

int get_option(int argc, char *argv[])
{
  char *option ;
  int nargs = 0;

  option=argv[1]+1; // pass -
  if (strcmp(option, "cras")==0)
  {
    cout << "use the scanner RAS coordinates for label position." << endl;
    useRealRAS=1;
  }
  else if (strcmp(option,"fillup") == 0)
  {
    cout << "check to make sure all points are filled." << endl;
    fillup = 1;
  }
  else if (strcmp(option,"scale") == 0)
  {
    cout << "we scale the label positions by " << argv[2] << endl;
    scale = atof(argv[2]);
    nargs=1;
  }
  return (nargs);
}
    
  
