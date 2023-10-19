/**
 * @brief Uses the 4 surfaces of a scan to construct a mask volume
 *
 * Uses the 4 surfaces of a scan to construct a mask volume showing the
 * position of each voxel with respect to the surfaces - GM, WM, LH or RH.
 */
/*
 * Original Author: Gheorghe Postelnicu
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

// STL
#include <queue>
#include <string>
#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>

// option to not use VTK libs to perform surface/volume intersection detection
// this is useful to build a version w/o the VTK/GL dependency, which is
// problematic for cluster machines which may not have those libs
#ifndef NO_VTK
// VTK
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkStructuredPoints.h>
#include <vtkImplicitModeller.h>
#include <vtkOBBTree.h>
#endif
#include "cmd_line_interface.h"

// FS


#include "fsenv.h"
#include "mrisurf.h"
#include "mri.h"
#include "mri2.h"
#include "error.h"
#include "cma.h"
#include "diag.h"
#include "macros.h"
#include "gca.h"
;
const char *Progname;

// static function declarations
// forward declaration
struct IoParams;
using namespace std;

class IoError : public std::exception
{
public:
  IoError(const std::string& what)
      : m_what(what)
  {}
  virtual ~IoError() throw()
  {}
  using std::exception::what;
  virtual const char* what()
  {
    return m_what.c_str();
  }
private:
  std::string m_what;
};

/*
  loads all the input files
  the return value is the output path for the volume files
*/
std::string
operator/(const std::string& a,
          const std::string& b)
{
  std::string s(a);
  s += "/" + b;
  return s;
}

std::string
LoadInputFiles(const IoParams& params,
               MRI*& mriInput,
               MRIS*& surfLeftWhite,
               MRIS*& surfLeftPial,
               MRIS*& surfRightWhite,
               MRIS*& surfRightPial);

void ComputeSurfaceDistanceFunction
(MRIS* mris, //input surface
 float thickness, // cap value for the precise distance value
 MRI* mriInOut); // output mri structure

MRI* CreateHemiMask(MRI* dpial, MRI* dwhite,
                    const unsigned char lblWhite,
                    const unsigned char lblRibbon,
                    const unsigned char lblBackground);

MRI* CombineMasks(MRI* maskOne,
                  MRI* maskTwo,
                  const unsigned char lblBackground);

//  return a binary mask
MRI* FilterLabel(MRI* initialMask,
                 const unsigned char lbl);

/*
  IO structure
  
  The user can either specify an input subject (thus implicitly
  setting the names of the in-out files)
  or use an advanced mode and manually specify the files
*/
struct IoParams
{
  typedef std::string StringType;

  StringType  subject;
  StringType  subjectsDir;

  StringType  templateMri;
  StringType  surfWhiteRoot;
  StringType  surfPialRoot;

  StringType  outRoot;

  unsigned char labelLeftWhite;
  unsigned char labelLeftRibbon;
  unsigned char labelRightWhite;
  unsigned char labelRightRibbon;
  unsigned char labelBackground;

  float capValue;

  bool bSaveDistance;
  bool bSaveRibbon;
  bool bEditAseg;

  IoParams();
  void parse(int ac, char* av[]);
};

int
main(int ac, char* av[])
{
  // parse command-line
  IoParams params;
  try
  {
    params.parse(ac,av);
  }
  catch (std::exception& excp)
  {
    std::cerr << " Exception caught while parsing the command-line\n"
    << excp.what() << std::endl;
    exit(1);
  }

  // process input files
  // will also resolve the paths depending on the mode of the application
  // (namely if the subject option has been specified or not)
  MRI* mriTemplate;
  MRIS* surfLeftWhite = NULL;
  MRIS* surfLeftPial = NULL;
  MRIS* surfRightPial = NULL;
  MRIS* surfRightWhite = NULL;

  std::string outputPath;

  try
  {
    outputPath = LoadInputFiles(params,
                                mriTemplate,
                                surfLeftWhite,
                                surfLeftPial,
                                surfRightWhite,
                                surfRightPial);
  }
  catch (std::exception& e)
  {
    std::cerr << " Exception caught while processing input files \n"
    << e.what() << std::endl;
    exit(1);
  }

  if (params.bEditAseg)
  {
    std::string subjDir = params.subjectsDir / params.subject;
    std::string pathSurf(subjDir / "surf"),
    pathMriOutput = outputPath / "aseg.ribbon.mgz";

    printf("inserting LH into aseg...\n") ;
    insert_ribbon_into_aseg(mriTemplate, mriTemplate,
                            surfLeftWhite, surfLeftPial, LEFT_HEMISPHERE) ;
    printf("inserting RH into aseg...\n") ;
    insert_ribbon_into_aseg(mriTemplate, mriTemplate,
                            surfRightWhite, surfRightPial, RIGHT_HEMISPHERE);
    printf("writing output to %s\n",
           (const_cast<char*>( pathMriOutput.c_str() )));
    MRIwrite(mriTemplate,  (const_cast<char*>( pathMriOutput.c_str() ))) ;
    exit(0) ;
  }

  /*
    Process LEFT hemisphere
  */

  //---------------------
  // proces white surface - convert to voxel-space
#ifndef NO_VTK
  MRIConvertSurfaceVertexCoordinates(surfLeftWhite, mriTemplate);
#endif
  // allocate distance
  MRI* dLeftWhite = MRIalloc( mriTemplate->width,
                              mriTemplate->height,
                              mriTemplate->depth,
                              MRI_FLOAT );
  MRIcopyHeader(mriTemplate, dLeftWhite);

  // Computes the signed distance to given surface. Sign indicates
  // whether it is on the inside or outside. params.capValue -
  // saturation/clip value for distance.
#ifndef NO_VTK
  ComputeSurfaceDistanceFunction(surfLeftWhite,
                                 params.capValue,
                                 dLeftWhite);
#else
  std::cout << "computing distance to left white surface \n" ;
  MRIScomputeDistanceToSurface(surfLeftWhite, dLeftWhite, mriTemplate->xsize) ;
#endif
  // if the option is there, output distance
  if ( params.bSaveDistance )
    MRIwrite
    ( dLeftWhite,
      const_cast<char*>( (outputPath / "lh.dwhite." +
                          params.outRoot + ".mgz").c_str() )
    );

  //-----------------------
  // process pial surface
#ifndef NO_VTK
  MRIConvertSurfaceVertexCoordinates(surfLeftPial,  mriTemplate);
#endif
  MRI* dLeftPial = MRIalloc( mriTemplate->width,
                             mriTemplate->height,
                             mriTemplate->depth,
                             MRI_FLOAT);
  MRIcopyHeader(mriTemplate,dLeftPial);
#ifndef NO_VTK
  ComputeSurfaceDistanceFunction(surfLeftPial,
                                 params.capValue,
                                 dLeftPial);
#else
  std::cout << "computing distance to left pial surface \n" ;
  MRIScomputeDistanceToSurface(surfLeftPial, dLeftPial, mriTemplate->xsize) ;
#endif
  if ( params.bSaveDistance )
    MRIwrite
    ( dLeftPial,
      const_cast<char*>( (outputPath / "lh.dpial." +
                          params.outRoot + ".mgz").c_str() )
    );

  // combine them and create a mask for the left hemi. Must be
  // outside of white and inside pial. Creates labels for WM and Ribbon.
  MRI* maskLeftHemi = CreateHemiMask(dLeftPial,dLeftWhite,
                                     params.labelLeftWhite,
                                     params.labelLeftRibbon,
                                     params.labelBackground);
  // no need for the hemi distances anymore
  MRIfree(&dLeftWhite);
  MRIfree(&dLeftPial);

  /*
    Process RIGHT hemi
  */

  //-------------------
  // process white
#ifndef NO_VTK
  MRIConvertSurfaceVertexCoordinates(surfRightWhite,mriTemplate);
#endif
  MRI* dRightWhite = MRIalloc( mriTemplate->width,
                               mriTemplate->height,
                               mriTemplate->depth,
                               MRI_FLOAT);
  MRIcopyHeader(mriTemplate, dRightWhite);
#ifndef NO_VTK
  ComputeSurfaceDistanceFunction( surfRightWhite,
                                  params.capValue,
                                  dRightWhite);
#else
  std::cout << "computing distance to right white surface \n" ;
  MRIScomputeDistanceToSurface(surfRightWhite, 
                               dRightWhite, 
                               mriTemplate->xsize) ;
#endif
  if ( params.bSaveDistance )
    MRIwrite
    ( dRightWhite,
      const_cast<char*>( (outputPath / "rh.dwhite." +
                          params.outRoot + ".mgz").c_str() )
    );

  //--------------------
  // process pial
#ifndef NO_VTK
  MRIConvertSurfaceVertexCoordinates(surfRightPial, mriTemplate);
#endif
  MRI* dRightPial = MRIalloc( mriTemplate->width,
                              mriTemplate->height,
                              mriTemplate->depth,
                              MRI_FLOAT);
  MRIcopyHeader(mriTemplate, dRightPial);
#ifndef NO_VTK
  ComputeSurfaceDistanceFunction(surfRightPial,
                                 params.capValue,
                                 dRightPial);
#else
  std::cout << "computing distance to right pial surface \n" ;
  MRIScomputeDistanceToSurface(surfRightPial, dRightPial, mriTemplate->xsize) ;
#endif
  if ( params.bSaveDistance )
    MRIwrite
    ( dRightPial,
      const_cast<char*>( (outputPath / "rh.dpial." +
                          params.outRoot + ".mgz").c_str() )
    );
  // compute hemi mask
  MRI* maskRightHemi = CreateHemiMask(dRightPial, dRightWhite,
                                      params.labelRightWhite,
                                      params.labelRightRibbon,
                                      params.labelBackground);
  // no need for the hemi distances anymore
  MRIfree(&dRightWhite);
  MRIfree(&dRightPial);

  /*

  finally combine the two created masks -- need to resolve overlap

  */
  MRI* finalMask = CombineMasks(maskLeftHemi, maskRightHemi,
                                params.labelBackground);
  MRIcopyHeader( mriTemplate, finalMask);
  // write final mask

  std::cout << "writing volume " <<  
    const_cast<char*>( (outputPath / (params.outRoot +
                                      ".mgz")).c_str() ) << endl;

  MRIwrite( finalMask,
            const_cast<char*>( (outputPath / (params.outRoot +
                                              ".mgz")).c_str() )
          );
  // sanity-check: make sure location 0,0,0 is background (not brain)
  if ( MRIgetVoxVal(finalMask,0,0,0,0) != 0 )
  {
    cerr << "ERROR: ribbon has non-zero value at location 0,0,0"
         << endl;
    exit(1);      
  }

  /*
    if present, also write the ribbon masks by themselves
  */
  if ( params.bSaveRibbon )
  {
    std::cout << " writing ribbon files\n";
    // filter the mask of the left hemi
    MRI* ribbon = FilterLabel(maskLeftHemi,
                              params.labelLeftRibbon);
    MRIcopyHeader( mriTemplate, ribbon);
    MRIwrite( ribbon,
              const_cast<char*>( (outputPath / "lh." +
                                  params.outRoot + ".mgz").c_str()  )
            );
    // sanity-check: make sure location 0,0,0 is background (not brain)
    if ( MRIgetVoxVal(ribbon,0,0,0,0) != 0 )
    {
      cerr << "ERROR: lh ribbon has non-zero value at location 0,0,0"
           << endl;
      exit(1);      
    }
    MRIfree(&ribbon);

    ribbon = FilterLabel(maskRightHemi,
                         params.labelRightRibbon);
    MRIcopyHeader( mriTemplate, ribbon);
    MRIwrite( ribbon,
              const_cast<char*>( (outputPath / "rh." +
                                  params.outRoot + ".mgz").c_str() )
            );
    // sanity-check: make sure location 0,0,0 is background (not brain)
    if ( MRIgetVoxVal(ribbon,0,0,0,0) != 0 )
    {
      cerr << "ERROR: rh ribbon has non-zero value at location 0,0,0"
           << endl;
      exit(1);      
    }
    MRIfree(&ribbon);
  }

  return 0;
}


IoParams::IoParams()
{
  labelBackground = 0;

  labelLeftWhite = 20;
  labelRightWhite = 120;

  labelLeftRibbon = 10;
  labelRightRibbon = 110;

  capValue = 3;
  bSaveDistance = false;
  bEditAseg = false ;
  bSaveRibbon = false;

  outRoot = "ribbon";
  surfWhiteRoot = "white";
  surfPialRoot = "pial";

  char *sd = FSENVgetSUBJECTS_DIR();
  if (NULL == sd) subjectsDir = "";
  else subjectsDir = sd;
}

void
IoParams::parse(int ac, char* av[])
{
  std::string sl = "left_";
  std::string sr = "right_";
  std::string srib = "ribbon";
  std::string sw = "white";
  std::string ssurf = "surf_";
  std::string slbl = "label_";

  int iLeftWhite(labelLeftWhite),
  iLeftRibbon(labelLeftRibbon),
  iRightWhite(labelRightWhite),
  iRightRibbon(labelRightRibbon),
  iBackground(labelBackground);

  std::string strUse =  "surface root name (i.e. <subject>/surf/$hemi.<NAME>";
  CCmdLineInterface interface(av[0]);
  bool showHelp(false);

  interface.AddOptionBool
  ( "help", &showHelp, "display help message");
  interface.AddOptionString
  ( (ssurf+sw).c_str(), &surfWhiteRoot,
    (strUse + " - default value is white").c_str()
  );
  interface.AddOptionString
  ( (ssurf+"pial").c_str(), &surfPialRoot,
    (strUse + " - default value is pial").c_str()
  );
  interface.AddOptionString
  ( "out_root", &outRoot,
    "default value is ribbon - output will then be mri/ribbon.mgz and "
    "mri/lh.ribbon.mgz and mri/rh.ribbon.mgz "
    "(last 2 if -save_ribbon is used)"
  );
  interface.AddOptionInt
  ( (slbl+"background").c_str(), &iBackground,
    "override default value for background label value (0)"
  );
  interface.AddOptionInt
  ( (slbl+sl+sw).c_str(), &iLeftWhite,
    "override default value for left white matter label - 20"
  );
  interface.AddOptionInt
  ( (slbl+sl+srib).c_str(), &iLeftRibbon,
    "override default value for left ribbon label - 10"
  );
  interface.AddOptionInt
  ( (slbl+sr+sw).c_str(), &iRightWhite,
    "override default value for right white matter label - 120"
  );
  interface.AddOptionInt
  ( (slbl+sr+srib).c_str(), &iRightRibbon,
    "override default value for right ribbon label - 110"
  );
  interface.AddOptionFloat
  ( "cap_distance", &capValue,
    (char*)"maximum distance up to which the signed distance function "
    "computation is accurate"
  );
  interface.AddOptionBool
  ( "save_distance", &bSaveDistance,
    "option to save the signed distance function as ?h.dwhite.mgz "
    "?h.dpial.mgz in the mri directory"
  );
  interface.AddOptionBool
  ( "edit_aseg", &bEditAseg,
    "option to edit the aseg using the ribbons and save to "
    "aseg.ribbon.mgz in the mri directory"
  );
  interface.AddOptionBool
  ( "save_ribbon", &bSaveRibbon,
    "option to save just the ribbon for the hemispheres - "
    "in the format ?h.ribbon.mgz"
  );
  interface.AddOptionString
  ( "sd", &subjectsDir,
    "option to specify SUBJECTS_DIR, default is read from enviro"
  );
  interface.AddIoItem
  (&subject, " subject (required param!)");

  // if ac == 0, then print complete help
  if ( ac == 1 )
  {
    std::cout <<
    "\n"
    " Computes a volume mask, at the same resolution as the\n"
    " <subject>/mri/brain.mgz.  The volume mask contains 4 values:\n"
    "   LH_WM (default 10)\n"
    "   LH_GM (default 100)\n"
    "   RH_WM (default 20)\n"
    "   RH_GM (default 200)\n"
    " The algorithm uses the 4 surfaces situated in\n"
    " <subject>/surf/[lh|rh].[white|pial].surf and labels voxels\n"
    " based on the signed-distance function from the surface.\n";
    interface.PrintHelp();
    exit(0);
  }

  interface.Parse(ac,av);
  if ( showHelp )
  {
    exit(0);
  }

  labelLeftWhite = (unsigned char)(iLeftWhite);
  labelLeftRibbon = (unsigned char)(iLeftRibbon);
  labelRightWhite = (unsigned char)(iRightWhite);
  labelRightRibbon = (unsigned char)(iRightRibbon);
  labelBackground = (unsigned char)(iBackground);

}


std::string
LoadInputFiles(const IoParams& params,
               MRI*& mriTemplate,
               MRIS*& surfLeftWhite,
               MRIS*& surfLeftPial,
               MRIS*& surfRightWhite,
               MRIS*& surfRightPial)
{
  // determine the mode of the application and infer the input file names

  // declare path objects
  std::string pathSurfLeftWhite,
  pathSurfLeftPial,
  pathSurfRightWhite,
  pathSurfRightPial,
  pathMriInput,
  pathOutput;

  if ( params.subjectsDir.empty() )
  {
    cerr << "SUBJECTS_DIR not found. Use --sd <dir>, or set SUBJECTS_DIR"
    << endl;
    exit(1);
  }
  else
  {
    cout << "SUBJECTS_DIR is " << params.subjectsDir << endl;
  }

  if ( !params.subject.empty() ) // application is in subject-mode
  {
    std::string subjDir = params.subjectsDir / params.subject;
    std::string pathSurf(subjDir / "surf");
    pathSurfLeftWhite = pathSurf / "lh.white";
    pathSurfLeftPial = pathSurf / "lh.pial";
    pathSurfRightWhite = pathSurf / "rh.white";
    pathSurfRightPial = pathSurf / "rh.pial";
    pathMriInput = subjDir / "mri" / "aseg.mgz";

    pathOutput = subjDir / "mri";
  }

  printf("loading input data...\n") ;

  // load actual files now
  surfLeftWhite = MRISread
                  ( const_cast<char*>( pathSurfLeftWhite.c_str() )
                  );
  if ( !surfLeftWhite )
    throw IoError( std::string("failed to read left white surface ")
                   + pathSurfLeftWhite );

  surfLeftPial  = MRISread
                  ( const_cast<char*>( pathSurfLeftPial.c_str() )
                  );
  if ( !surfLeftPial )
    throw IoError( std::string("failed to read left pial surface ")
                   + pathSurfLeftPial );

  surfRightWhite = MRISread
                   ( const_cast<char*>( pathSurfRightWhite.c_str() )
                   );
  if ( !surfRightWhite )
    throw IoError( std::string("failed to read right white surface ")
                   + pathSurfRightWhite );

  surfRightPial = MRISread
                  ( const_cast<char*>( pathSurfRightPial.c_str() )
                  );
  if ( !surfRightPial )
    throw IoError( std::string("failed to read right pial surface ")
                   + pathSurfRightPial );

  mriTemplate = MRIread
                ( const_cast<char*>( pathMriInput.c_str() )
                );
  if ( !mriTemplate )
    throw IoError( std::string("failed to read template mri ")
                   + pathMriInput );

  return pathOutput;
}

#ifndef NO_VTK
void
ComputeSurfaceDistanceFunction(MRIS* mris,
                               float thickness,
                               MRI* vol)
{
  vtkPoints* points;
  vtkCellArray* faces;
  vtkPolyData* mesh;

  points = vtkPoints::New();
  points->SetNumberOfPoints( mris->nvertices );
  VERTEX* pvtx = &( mris->vertices[0] );
  for (unsigned int ui(0), nvertices(mris->nvertices);
       ui < nvertices; ++ui, ++pvtx )
  {
    points->SetPoint(ui, pvtx->x, pvtx->y, pvtx->z);
  } // next point (vertex)

  faces = vtkCellArray::New();
  faces->Allocate( mris->nfaces );
  FACE* pface = &( mris->faces[0] );
  for ( unsigned int ui(0), nfaces(mris->nfaces);
        ui < nfaces; ++ui, ++pface )
  {
    faces->InsertNextCell((long long int)VERTICES_PER_FACE);
    for ( int i = 0; i < VERTICES_PER_FACE; i++ )
    {
      faces->InsertCellPoint(pface->v[i]);
    }
  } // next cell (face)

  mesh = vtkPolyData::New();
  mesh->SetPoints(points);
  mesh->SetPolys(faces);
  mesh->Modified();

  vtkStructuredPoints* sp = vtkStructuredPoints::New();
  int dims[3] = { vol->width,
                  vol->height,
                  vol->depth };
  sp->SetDimensions(dims);

  // compute the unsigned distance
  vtkImplicitModeller* implicit = vtkImplicitModeller::New();
#if VTK_MAJOR_VERSION > 5
  implicit->SetInputData(mesh);
#else
  implicit->SetInput(mesh);
#endif
  double bounds[6] =
    { 0, (double)dims[0]-1.0, 0.0, (double)dims[1]-1.0, 0.0, (double)dims[2]-1.0 };
  implicit->SetModelBounds(bounds);
  implicit->SetMaximumDistance( thickness );
  implicit->SetSampleDimensions(dims);
  implicit->AdjustBoundsOff();
  implicit->CappingOff();
  implicit->SetCapValue(thickness);
  implicit->SetProcessModeToPerVoxel();
  implicit->SetNumberOfThreads( 8 );
  implicit->Update();

  sp->ShallowCopy( implicit->GetOutput() );

  // determine if voxels are inside or outside
  // minimize no calls to  vtkOBBTree::InsideOrOutside - costly

  vtkOBBTree* obb = vtkOBBTree::New();
  obb->CacheCellBoundsOff();
  obb->SetDataSet(mesh);
  obb->AutomaticOn();
  obb->SetTolerance(0);
  obb->BuildLocator();

  // array for voxels inside, outside or yet to determin
  int npoints = sp->GetNumberOfPoints();
  int *inout = new int[npoints];
  memset(inout, 0, npoints*sizeof(int));
  float *tab = (float*)sp->GetScalarPointer();

  std::queue<vtkIdType> filo;

  vtkIdType incx,incy,incz;
  sp->GetIncrements(incx,incy,incz);

  // iterate thru each volume point, apply sign to tab
  for (vtkIdType i=0; i<npoints; ++i)
  {
    if (inout[i]) continue; // skip if already visited

    double pt[3];
    sp->GetPoint(i,pt);
    const int _inout = obb->InsideOrOutside(pt); // sign
    if ( _inout == 0 ) std::cerr << "Warning: voxel could not be classified\n";

    filo.push(i);
    while (!filo.empty())
    {
      const vtkIdType j = filo.front();
      filo.pop();
      if (inout[j]) continue;

      inout[j] = _inout;
      const float dist = tab[j]; // unsigned distance
      tab[j] *= _inout; // make it signed

      // mark its neighbors if distance > 1 (triangle inequality)
      if ( dist>1 )
      {
        const int x = j%incy;          //col
        const int y = (j%incz) / incy; //row
        const int z = j/incz;          //slice
        // If not on boundary of volume and have not visited neighbor,
        // then push
        if (x>0 && !inout[j-incx]) filo.push(j-incx);
        if (y>0 && !inout[j-incy]) filo.push(j-incy);
        if (z>0 && !inout[j-incz]) filo.push(j-incz);
        if (x<dims[0]-1 && !inout[j+incx]) filo.push(j+incx);
        if (y<dims[1]-1 && !inout[j+incy]) filo.push(j+incy);
        if (z<dims[2]-1 && !inout[j+incz]) filo.push(j+incz);
      }
    }
  } // next i

  for (vtkIdType j=0; j<npoints; ++j)
  {
    const int x = j%incy;
    const int y = (j%incz) / incy;
    const int z = j/incz;

    MRIsetVoxVal(vol, x,y,z,0, tab[j] );
  }

  // free the objects that we New'd!  they want to be free!!!!!
  // according to vtkObject.h, the Delete() method should be used to delete.
  obb->Delete();
  implicit->Delete();
  sp->Delete();
  mesh->Delete();
  faces->Delete();
  points->Delete();
}
#endif

MRI*
CreateHemiMask(MRI* dpial,
               MRI* dwhite,
               const unsigned char lblWhite,
               const unsigned char lblRibbon,
               const unsigned char lblBackground)
{
  // allocate return volume
  MRI* mri = MRIalloc(dpial->width,
                      dpial->height,
                      dpial->depth,
                      MRI_UCHAR);

  for (unsigned int z(0), depth(dpial->depth);
       z<depth; ++z)
    for (unsigned int y(0), height(dpial->height);
         y<height; ++y)
      for (unsigned int x(0), width(dpial->width);
           x<width; ++x)
      {
        if (x == (unsigned int)Gx && 
            y == (unsigned int)Gy && 
            z == (unsigned int)Gz)
          DiagBreak() ;
#ifndef NO_VTK
        if ( MRIFvox(dwhite,x,y,z) > 0 )
          MRIsetVoxVal(mri, x,y,z,0, lblWhite);
        else if ( MRIFvox(dpial,x,y,z) > 0 )
          MRIsetVoxVal(mri, x,y,z,0, lblRibbon);
        else
          MRIsetVoxVal(mri,x,y,z,0, lblBackground);
#else
        if ( MRIgetVoxVal(dwhite,x,y,z,0) < 0 )
          MRIsetVoxVal(mri, x,y,z,0, lblWhite);
        else if ( MRIgetVoxVal(dpial,x,y,z,0) < 0 )
          MRIsetVoxVal(mri, x,y,z,0, lblRibbon);
        else
          MRIsetVoxVal(mri,x,y,z,0, lblBackground);
#endif
      } // next x,y,z

  return mri;
}

/*
  implicit assumption the masks do not overlap
  if this is the case, a message will be printed
*/
MRI* CombineMasks(MRI* maskOne,
                  MRI* maskTwo,
                  const unsigned char lblBackground)
{
  MRI* mri = MRIalloc(maskOne->width,
                      maskOne->height,
                      maskOne->depth,
                      MRI_UCHAR);
  unsigned char voxOne,voxTwo;

  unsigned int overlap = 0;
  for (unsigned int z(0), depth(maskOne->depth);
       z<depth; ++z)
    for (unsigned int y(0), height(maskOne->height);
         y<height; ++y)
      for (unsigned int x(0), width(maskOne->width);
           x<width; ++x)
      {
        voxOne = MRIvox(maskOne,x,y,z);
        voxTwo = MRIvox(maskTwo,x,y,z);
        if ( voxOne!=lblBackground && voxTwo!=lblBackground )
        {
          // overlap
          ++overlap;
          MRIsetVoxVal(mri,x,y,z,0,voxOne);
        }
        else if ( voxOne == lblBackground )
          MRIsetVoxVal(mri,x,y,z,0,voxTwo);
        else if ( voxTwo == lblBackground )
          MRIsetVoxVal(mri,x,y,z,0,voxOne);
      } // next x,y,z

  std::cout << " hemi masks overlap voxels = " << overlap << std::endl;
  return mri;
}


MRI* FilterLabel(MRI* initialMask,
                 const unsigned char lbl)
{
  MRI* mri = MRIalloc( initialMask->width,
                       initialMask->height,
                       initialMask->depth,
                       MRI_UCHAR);

  for (unsigned int z(0), depth(initialMask->depth);
       z<depth; ++z)
    for (unsigned int y(0), height(initialMask->height);
         y<height; ++y)
      for (unsigned int x(0), width(initialMask->width);
           x<width; ++x)
      {
        if ( MRIvox(initialMask,x,y,z) == lbl )
          MRIsetVoxVal(mri, x,y,z,0, (unsigned char)(1) );
        else
          MRIsetVoxVal(mri, x,y,z,0, (unsigned char)(0) );
      } // next x,y,z

  return mri;
}


