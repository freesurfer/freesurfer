//
// mri_parselabel.cpp
//
//
#include <iostream>
#include <iomanip>

// workaround for gcc2.96
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif

#include <string>

#include <fstream>
#include <vector>
#include <algorithm>

extern "C" {
#include "mri.h"
#include "version.h"
#include "macros.h"
#include "transform.h"
#include "talairachex.h"
  char *Progname;
}

using namespace std;

int useRealRAS = 0;
int fillup=0;
double scale = 1;
string xfname;
int invert=0;

// voxel size
double vx;
double vy;
double vz;

class Vertex
{
public:
  Vertex(double x=0, double y=0, double z=0, double value=0) 
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
MRI *erodeRegion(MRI *mri, MRI *cur, unsigned char val);
int selectAddedRegion(MRI *mri, MRI *orig, unsigned char val);

int (*voxToRAS)(MRI *mri, Real x, Real y, Real z, Real *xs, Real *ys, Real *zs);
int (*rasToVox)(MRI *mri, Real x, Real y, Real z, Real *xs, Real *ys, Real *zs);

#define V4_LOAD(v, x, y, z, r)  (VECTOR_ELT(v,1)=x, VECTOR_ELT(v,2)=y, \
                                  VECTOR_ELT(v,3)=z, VECTOR_ELT(v,4)=r) ;

int main(int argc, char *argv[])
{
  int nargs;
  Progname=argv[0];

  nargs = handle_version_option (argc, argv, "$Id: mri_parselabel.cpp,v 1.10 2004/06/10 16:53:57 tosa Exp $", "$Name:  $");
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
    cerr << "        -xfm <xfm>   : use the xfm to transform the vertices " << endl;
    cerr << "                     : xfm must be from highres to lowres." << endl;
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
  if (xfname.size() != 0)
  {
    cout << "xfm        : " << xfname.c_str() << endl;
    cout << "invert     : " << invert << endl;
  }
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

  //////////////////////////////////////////////////////////////////
  LTA *lta = 0;
  if (xfname.size())
  {
    lta = LTAreadEx(const_cast<char *>(xfname.c_str()));
    if (!lta->xforms[0].src.valid)
    {
      cerr << "could not find the src volume" << lta->xforms[0].src.fname << endl;
      cerr << " need to calculate the correct transform" << endl;
      return -1;
    }
    if (!lta->xforms[0].dst.valid)
    {
      cerr << "could not find the dst volume" << lta->xforms[0].dst.fname << endl;
      cerr << " need to calculate the correct transform" << endl;
      return -1;
    }
  }

  MRI *mriIn=MRIread(argv[2]);
  if (!mriIn)
  {
    cerr << "Could not open input volume : " << argv[2] << endl;
    return -1;
  }
  // setting the voxel size
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
  cout << "reading vertices ..." << endl;

  //////////////////////////////////////////////////////////////////////
  // reading vertices 
  //////////////////////////////////////////////////////////////////////
  VECTOR *vIn=0;
  VECTOR *vOut=0;
  if (xfname.size())
  {
    vIn = VectorAlloc(4, MATRIX_REAL);
    vOut = VectorAlloc(4, MATRIX_REAL);
  }
  MATRIX *hSRASTolSRAS = 0;
  if (xfname.size())
  {
    //  conformed -------> surfaceRAS
    //      |                  |
    //      V                  V  RASFromSRAS
    //    hres    ------->   RAS
    //      |                  |   xfm
    //      V                  V
    //    lres    ------->   RAS
    //      |                  |  sRASFromRAS
    //      V                  V
    //  conformed --------> surfaceRAS

    cout << "calculate src surfaceRAS to dst surfaceRAS matrix" << endl;
    // must convert to RAS first
    MATRIX *RASFromSRAS = MatrixAlloc(4, 4, MATRIX_REAL);
    MatrixIdentity(4, RASFromSRAS);
    *MATRIX_RELT(RASFromSRAS, 1,4) = lta->xforms[0].src.c_r;
    *MATRIX_RELT(RASFromSRAS, 2,4) = lta->xforms[0].src.c_a;
    *MATRIX_RELT(RASFromSRAS, 3,4) = lta->xforms[0].src.c_s;
    
    MATRIX *sRASFromRAS;
    sRASFromRAS = MatrixAlloc(4, 4, MATRIX_REAL);
    MatrixIdentity(4, sRASFromRAS);
    *MATRIX_RELT(sRASFromRAS, 1,4) = - lta->xforms[0].dst.c_r;
    *MATRIX_RELT(sRASFromRAS, 2,4) = - lta->xforms[0].dst.c_a;
    *MATRIX_RELT(sRASFromRAS, 3,4) = - lta->xforms[0].dst.c_s;
    
    MATRIX *tmpM = MatrixMultiply(lta->xforms[0].m_L, RASFromSRAS, NULL);
    hSRASTolSRAS = MatrixMultiply(sRASFromRAS, tmpM, NULL); 

    MatrixFree(&RASFromSRAS);
    MatrixFree(&sRASFromRAS);
    MatrixFree(&tmpM);
    cout << "calculation done." << endl;
  }
  // read label
  int count =0;
  while (flabel.good())
  {
    Vertex v;
    int num;
    double x, y, z;
    double value;
    flabel >> num >> x >> y >> z >> value;
    // cout << "(" << x << ", " << y << ", " << z << ")" << endl;
    if (xfname.size())
    {
      double xt, yt, zt;
      TransformWithMatrix(hSRASTolSRAS, scale*x, scale*y, scale*z, &xt, &yt, &zt);
      // create a vector
      v = Vertex(xt, yt, zt, value);
      if (count < 10)
      {
	cout << "transformed: (" << scale*x << ", " << scale*y << ", " << scale*z << ") " 
	     << " to  (" << xt << ", " << yt << ", " << zt << ") " << endl;
      }
      count++;
    }
    else
    {
      V4_LOAD(vIn, scale*x, scale*y, scale*z, 1.);
      v = Vertex(scale*x, scale*y, scale*z, value);
    }
    vertices.push_back(v);
  }
  count = 0;
  flabel.close();
  // free vector
  if (vIn)
    VectorFree(&vIn);
  if (vOut)
    VectorFree(&vOut);
  if (lta)
    LTAfree(&lta);
  if (hSRASTolSRAS)
    MatrixFree(&hSRASTolSRAS);

#if 0
  // sort
  cout << "sorting ..." << endl;
  sort(vertices.begin(), vertices.end(), LessX());
  cout << "sorting done" << endl;
#endif

  /////////////////////////////////////////////////////////////////////////////////
  // create a volume
  /////////////////////////////////////////////////////////////////////////////////
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

  cout << "filling...." << endl;
  int added = 0;
  // first use the naive approach
  for (size_t i=0; i < vertices.size(); ++i)
  {
    rasToVox(mriIn, vertices[i].x_, vertices[i].y_, vertices[i].z_, &xv, &yv, &zv);
    if (count < 10)
    {
      cout << "get voxel at (" << xv << ", " << yv << ", " << zv << ")" 
	   << "  from vertex (" << vertices[i].x_ << ", " << vertices[i].y_ << ", " << vertices[i].z_ << ") " << endl;
    }
    count++;

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
    //////////////////////////////////////////////////////////////////////////
    // check possible others.  look around 3x3x3 neighbors
    for (int z0=-1; z0 < 2; ++z0)
      for (int y0=-1; y0 < 2; ++y0)
	for (int x0=-1; x0 < 2; ++x0)
	{
	  int x = nint(xv) + x0;
	  int y = nint(yv) + y0;
	  int z = nint(zv) + z0;
	  Real xr, yr, zr;
	  // go back to RAS to see whether the point is close enough
	  voxToRAS(mriIn, x, y, z, &xr, &yr, &zr);
	  // if the difference is 1/2 voxel size mm then equal
	  if (Vertex(xr, yr, zr) == vertices[i] && (x0+y0+z0) != 0)
	  {
	    MRIvox(mriOut, x, y, z ) = (unsigned char) val;
	    cout << "\nadded another voxel point: ( " 
		 << x << ", " << y << ", " << z <<")" << "for ras point: (" 
		 << xr << ", " << yr << ", " << zr << ")"  
	         << " for vertex (" << vertices[i].x_ << ", " << vertices[i].y_ << ", " << vertices[i].z_ << ") "
	         << endl; 
	    added++;
	  }
	}
    //////////////////////////////////////////////////////////////////////////
  }
  cout << "filling done" << endl;
  if (added != 0)
    cout << "added points " << added << endl;
  MRIfree(&mriIn);

  MRI *mask =0;

  // the following is unnecessary, since I check the neightbor above
  if (fillup)
  {

    // create a mask with eroded region
    mask = MRIcopy(mriOut, NULL);
    cout << "eroding region..." << endl;
    // add 3x3x3 neighbors
    mask = erodeRegion(mask, mriOut, val);

    // MRIwrite(mask, "./eroded.mgh");

    // get added points only
    cout << "select only added points" << endl;
    int addedPoints = selectAddedRegion(mask, mriOut, val);

    int fivepercent = addedPoints/20;

    cout << "check missed positions...." << endl;

    // now check whether missing points exist or not
    unsigned char curval =0;
    int addedcount=0;
    int totcount = 0;
    int percentage = 0;
    Real xl, yl, zl;
    for (int z = 0; z < mriOut->depth; z++)
    {
      for (int y = 0; y < mriOut->height; y++)
      {
	for (int x = 0; x < mriOut->width; x++)
	{
	  // mask has only added points marked with val
	  curval = MRIvox(mask, x, y, z);
	  // get the RAS point
	  voxToRAS(mriOut, x, y, z, &xl, &yl, &zl);

	  // this point is marked, then looking around to see other points
	  // may belong to the label list.
	  if (curval == val)
	  {
	    totcount++;
	    if (totcount % fivepercent == 0)
	    {
	      percentage += 5;
	      cout << percentage << "% done" << endl;
	    }
	    // search
	    for (size_t i=0; i < vertices.size(); ++i)
	    {
	      // since I sorted by X, the vertex position > current position + 2x vx
	      // don't do anymore
	      if ((vertices[i].x_) > (xl + 2*vx))
		goto nextpoint;
	      
	      if (vertices[i] == Vertex(xl, yl, zl, 0))
	      {
		// if it is not set
		if (MRIvox(mriOut, x, y, z) != val)
		{
		  MRIvox(mriOut, x, y, z) = val;
		  cout << "\nadded another point: ( " 
		       << x << ", " << y << ", " << z <<")" << endl; 
		  addedcount++;
		}
	      } // for i		    
	    nextpoint:
	      continue;
	    }
	  }
	} // x
      } // y
    } // z
    cout << endl;
    cout << "       added = " << addedcount << endl;
  }
  
  MRIwrite(mriOut, argv[3]);

  MRIfree(&mriOut);

  if (mask)
    MRIfree(&mask);

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
    cout << "check to make sure all points are filled.... no longer necessary" << endl;
    fillup = 1;
  }
  else if (strcmp(option,"scale") == 0)
  {
    cout << "we scale the label positions by " << argv[2] << endl;
    scale = atof(argv[2]);
    nargs=1;
  }
  else if (strcmp(option, "xfm") == 0)
  {
    cout << "use xfm " << argv[2] << endl;
    cout << "make sure that this is the xfm from high res to low res" << endl;
    cout << "if the other way, you must use -invert option also" << endl;
    xfname = argv[2];
    nargs=1;
  }
  else if (strcmp(option, "invert")==0)
  {
    cout << "xfm will be inverted" << endl;
    invert=1;
    nargs=0;
  }
  return (nargs);
}
    
MRI *erodeRegion(MRI *mask, MRI *current, unsigned char val)
{  
  int xi, yi, zi;
  unsigned char curval;
  for (int z = 0; z < mask->depth; z++)
  {
    for (int y = 0; y < mask->height; y++)
    {
      for (int x = 0; x < mask->width; x++)
      {
	curval = MRIvox(current, x, y, z);
	// this point is marked, then looking around to see other points
	// may belong to the label list.
	if (curval == val)
	{
	  // look around 3 x 3 x 3 neighbors
	  /////////////////////////////////////////////
	  for (int z0 = -1 ; z0 <= 1 ; z0++)
	  {
	    zi = mask->zi[z+z0] ;
	    for (int y0 = -1 ; y0 <= 1 ; y0++)
	    {
	      yi = mask->yi[y+y0] ;
	      for (int x0 = -1 ; x0 <= 1 ; x0++)
	      {
		xi = mask->xi[x+x0] ;
		MRIvox(mask, xi, yi, zi) = val;
	      }
	    }
	  } 
	} 
      } // x
    } // y
  } // z
  return mask;
}

int selectAddedRegion(MRI *mask, MRI *orig, unsigned char val)
{
  int totcount = 0;
  int origcount = 0;
  unsigned char mval;
  unsigned char curval;
  for (int z = 0; z < mask->depth; z++)
  {
    for (int y = 0; y < mask->height; y++)
    {
      for (int x = 0; x < mask->width; x++)
      {
	mval = MRIvox(mask, x, y, z);
	curval = MRIvox(orig, x, y, z);
	
	// this point is marked, then looking around to see other points
	// may belong to the label list.
	if (mval == val)
	{
	  totcount++;
	  if (curval == val)
	  {
	    MRIvox(mask, x, y, z) = 0;
	    origcount++;
	  }
	} 
      } // x
    } // y
  } // z
  printf("added points = %d \n", totcount-origcount);
  printf("total count = %d, naive count = %d\n", totcount, origcount);
  printf("total voxel count = %d\n", mask->depth*mask->height*mask->width);
  return (totcount-origcount);
}
