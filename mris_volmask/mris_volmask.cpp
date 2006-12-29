/**
 * @file  mris_volmask.cpp
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:11 $
 *    $Revision: 1.6 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */



// STL
#include <iostream>
#include <queue>
#include <string>

// VTK
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkStructuredPoints.h>
#include <vtkImplicitModeller.h>
#include <vtkOBBTree.h>

#include "cmd_line_interface.h"

// FS
extern "C" {
#include "fsenv.h"
#include "mrisurf.h"
#include "mri.h"
};
char *Progname;

// static function declarations

// forward declaration
struct IoParams;

class IoError : public std::exception {
public:
  IoError(const std::string& what)
      : m_what(what) {}
  virtual ~IoError() throw() {}
  virtual const char* what() {
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
          const std::string& b) {
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
               MRIS*& surfRightPial) throw(IoError) ;

/*

converts vertices to the index space

*/
void ConvertSurfaceVertexCoordinates(MRI_SURFACE* mris,
                                     MRI* vol);

void ComputeSurfaceDistanceFunction(MRIS* mris, //input surface
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

The user can either specify an input subject (thus implicitly setting the names of the in-out files)
or use an advanced mode and manually specify the files

*/
struct IoParams {
  typedef std::string StringType;


  StringType  subject;

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

  IoParams();
  void parse(int ac, char* av[]);
};

int
main(int ac, char* av[]) {
  // parse command-line
  IoParams params;
  try {
    params.parse(ac,av);
  } catch (std::exception& excp) {
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

  try {
    outputPath = LoadInputFiles(params,
                                mriTemplate,
                                surfLeftWhite,
                                surfLeftPial,
                                surfRightWhite,
                                surfRightPial);
  } catch (std::exception& e) {
    std::cerr << " Exception caught while processing input files \n"
    << e.what() << std::endl;
    exit(1);
  }

  /*
    Process LEFT hemisphere
  */


  //---------------------
  // proces white surface
  ConvertSurfaceVertexCoordinates(surfLeftWhite, mriTemplate);
  // allocate distance
  MRI* dLeftWhite = MRIalloc( mriTemplate->width,
                              mriTemplate->height,
                              mriTemplate->depth,
                              MRI_FLOAT );
  MRIcopyHeader(mriTemplate, dLeftWhite);
  ComputeSurfaceDistanceFunction(surfLeftWhite,
                                 params.capValue,
                                 dLeftWhite);
  // if the option is there, output distance
  if ( params.bSaveDistance )
    MRIwrite
    ( dLeftWhite,
      const_cast<char*>( (outputPath / "lh.dpial." + params.outRoot + ".mgz")
                         .c_str() )
    );

  //-----------------------
  // process pial surface
  ConvertSurfaceVertexCoordinates(surfLeftPial,  mriTemplate);
  MRI* dLeftPial = MRIalloc( mriTemplate->width,
                             mriTemplate->height,
                             mriTemplate->depth,
                             MRI_FLOAT);
  MRIcopyHeader(mriTemplate,dLeftPial);
  ComputeSurfaceDistanceFunction(surfLeftPial,
                                 params.capValue,
                                 dLeftPial);
  if ( params.bSaveDistance )
    MRIwrite
    ( dLeftPial,
      const_cast<char*>( (outputPath / "lh.dpial." + params.outRoot + ".mgz")
                         .c_str() )
    );

  // combine them and create a mask for the left hemi
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
  ConvertSurfaceVertexCoordinates(surfRightWhite,mriTemplate);
  MRI* dRightWhite = MRIalloc( mriTemplate->width,
                               mriTemplate->height,
                               mriTemplate->depth,
                               MRI_FLOAT);
  MRIcopyHeader(mriTemplate, dRightWhite);
  ComputeSurfaceDistanceFunction( surfRightWhite,
                                  params.capValue,
                                  dRightWhite);
  if ( params.bSaveDistance )
    MRIwrite
    ( dRightWhite,
      const_cast<char*>( (outputPath / "rh.dwhite." + params.outRoot + ".mgz")
                         .c_str() )
    );

  //--------------------
  // process pial
  ConvertSurfaceVertexCoordinates(surfRightPial, mriTemplate);
  MRI* dRightPial = MRIalloc( mriTemplate->width,
                              mriTemplate->height,
                              mriTemplate->depth,
                              MRI_FLOAT);
  MRIcopyHeader(mriTemplate, dRightPial);
  ComputeSurfaceDistanceFunction(surfRightPial,
                                 params.capValue,
                                 dRightPial);
  if ( params.bSaveDistance )
    MRIwrite
    ( dRightPial,
      const_cast<char*>( (outputPath / "rh.dpial." + params.outRoot + ".mgz")
                         .c_str() )
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

  finally combine the two created masks

  */
  MRI* finalMask = CombineMasks(maskLeftHemi, maskRightHemi,
                                params.labelBackground);
  MRIcopyHeader( mriTemplate, finalMask);
  // write final mask
  MRIwrite( finalMask,
            const_cast<char*>( (outputPath / (params.outRoot + ".mgz"))
                               .c_str() )
          );

  /*
     if present, also write the ribbon masks by themselves
  */
  if ( params.bSaveRibbon ) {
    std::cout << " writing ribbon files\n";
    // filter the mask of the left hemi
    MRI* ribbon = FilterLabel(maskLeftHemi,
                              params.labelLeftRibbon);
    MRIcopyHeader( mriTemplate, ribbon);
    MRIwrite( ribbon,
              const_cast<char*>( (outputPath / "lh." + params.outRoot + ".mgz")
                                 .c_str()  )
            );
    MRIfree(&ribbon);

    ribbon = FilterLabel(maskRightHemi,
                         params.labelRightRibbon);
    MRIcopyHeader( mriTemplate, ribbon);
    MRIwrite( ribbon,
              const_cast<char*>( (outputPath / "rh." + params.outRoot + ".mgz")
                                 .c_str() )
            );
    MRIfree(&ribbon);
  }

  return 0;
}




IoParams::IoParams() {
  labelBackground = 0;
  labelLeftWhite = 10;
  labelRightWhite = 20;
  labelLeftRibbon = 100;
  labelRightRibbon = 200;

  capValue = 3;
  bSaveDistance = false;
  bSaveRibbon = false;

  outRoot = "ribbon";
  surfWhiteRoot = "white";
  surfPialRoot = "pial";
}

void
IoParams::parse(int ac, char* av[]) {
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


  std::string strUse =  " surface root name (i.e. <subject>/surf/$hemi.<NAME>";
  CCmdLineInterface interface(av[0]);
  bool showHelp(false);


  interface.AddOptionBool( "help", &showHelp, " display help message");
  interface.AddOptionString( (ssurf+sw).c_str(), &surfWhiteRoot,
                             (strUse + " - default value is white").c_str()
                           );
  interface.AddOptionString( (ssurf+"pial").c_str(), &surfPialRoot,
                             (strUse + " - default value is pial").c_str()
                           );
  interface.AddOptionString( "out_root", &outRoot,
                             " default value is ribbon - output will then be mri/ribbon.mgz and  mri/lh.ribbon.mgz and mri/rh.ribbon.mgz (last 2 if -save_ribbon is used)"
                           );
  interface.AddOptionInt( (slbl+"background").c_str(), &iBackground,
                          " override default value for background label value (0)"
                        );
  interface.AddOptionInt( (slbl+sl+sw).c_str(), &iLeftWhite,
                          " override default value for left white matter label - 100"
                        );
  interface.AddOptionInt( (slbl+sl+srib).c_str(), &iLeftRibbon,
                          " override default value for left ribbon label - 10"
                        );
  interface.AddOptionInt( (slbl+sr+sw).c_str(), &iRightWhite,
                          " override default value for right white matter label - 200"
                        );
  interface.AddOptionInt( (slbl+sr+srib).c_str(), &iRightRibbon,
                          " override default value for right ribbon label - 20"
                        );
  interface.AddOptionFloat( "cap_distance", &capValue,
                            " maximum distance up to which the signed distance function computation is accurate"
                          );
  interface.AddOptionBool( "save_distance", &bSaveDistance,
                           " option to save the signed distance function as ?h.dwhite.mgz ?h.dpial.mgz in the mri directory"
                         );
  interface.AddOptionBool( "save_ribbon", &bSaveRibbon,
                           " option to save just the ribbon for the hemispheres - in the format ?h.outRoot.mgz"
                         );
  interface.AddIoItem(&subject, " subject");

  // if ac == 0, then print complete help
  if ( ac == 1 ) {
    std::cout << " Summary Description : \n\n"
    << " Computes a volume mask, at the same resolution as the <subject>/mri/brain.mgz.\n"
    << " The volume mask contains 4 values: LH_WM (default 10), LH_GM (default 100), RH_WM (default 20), RH_GM(default 200)\n"
    << " The algorithm uses the 4 surfaces situated in <subject>/surf/[lh|rh].[white|pial].surf and labels voxels based on the signed-distance function from the surface.\n";
    interface.PrintHelp();
    exit(0);
  }


  interface.Parse(ac,av);
  if ( showHelp ) {
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
               MRIS*& surfRightPial) throw(IoError) {
  // determine the mode of the application and infer the input file names
  // use the BOOST filesystem library to resolve paths (system-independent)

  // declare path objects
  std::string pathSurfLeftWhite,
  pathSurfLeftPial,
  pathSurfRightWhite,
  pathSurfRightPial,
  pathMriInput,
  pathOutput;

  FSENV* fsenv= FSENVgetenv();

  if ( !params.subject.empty() ) // application is in subject-mode
  {
    std::string subjDir = std::string( fsenv->SUBJECTS_DIR ) / params.subject;
    std::string pathSurf(subjDir / "surf");
    pathSurfLeftWhite = pathSurf / "lh.white";
    pathSurfLeftPial = pathSurf / "lh.pial";
    pathSurfRightWhite = pathSurf / "rh.white";
    pathSurfRightPial = pathSurf / "rh.pial";
    pathMriInput = subjDir / "mri" / "orig.mgz";

    pathOutput = subjDir / "mri";
  }

  FSENVfree(&fsenv);

  // load actual files now
  surfLeftWhite = MRISread
                  ( const_cast<char*>( pathSurfLeftWhite.c_str() )
                  );
  if ( !surfLeftWhite ) throw IoError( std::string("failed to read left white surface ")
                                         + pathSurfLeftWhite );

  surfLeftPial  = MRISread
                  ( const_cast<char*>( pathSurfLeftPial.c_str() )
                  );
  if ( !surfLeftPial ) throw IoError( std::string("failed to read left pial surface ")
                                        + pathSurfLeftPial );

  surfRightWhite = MRISread
                   ( const_cast<char*>( pathSurfRightWhite.c_str() )
                   );
  if ( !surfRightWhite ) throw IoError( std::string("failed to read right white surface ")
                                          + pathSurfRightWhite );

  surfRightPial = MRISread
                  ( const_cast<char*>( pathSurfRightPial.c_str() )
                  );
  if ( !surfRightPial ) throw IoError( std::string("failed to read right pial surface ")
                                         + pathSurfRightPial );

  mriTemplate = MRIread
                ( const_cast<char*>( pathMriInput.c_str() )
                );
  if ( !mriTemplate ) throw IoError( std::string("failed to read template mri ")
                                       + pathMriInput );

  return pathOutput;
}

void
ComputeSurfaceDistanceFunction(MRIS* mris,
                               float thickness,
                               MRI* vol) {
  vtkPoints* points;
  vtkCellArray* faces;
  vtkPolyData* mesh;

  points = vtkPoints::New();
  points->SetNumberOfPoints( mris->nvertices );
  VERTEX* pvtx = &( mris->vertices[0] );
  for (unsigned int ui(0), nvertices(mris->nvertices);
       ui < nvertices; ++ui, ++pvtx ) {
    points->SetPoint(ui, pvtx->x, pvtx->y, pvtx->z);
  } // next ui

  faces = vtkCellArray::New();
  faces->Allocate( mris->nfaces );
  FACE* pface = &( mris->faces[0] );
  for ( unsigned int ui(0), nfaces(mris->nfaces);
        ui < nfaces; ++ui, ++pface ) {
    faces->InsertNextCell(VERTICES_PER_FACE,
                          pface->v );
  } // next ui

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
  implicit->SetInput(mesh);
  double bounds[6] =
    { 0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1 };
  implicit->SetModelBounds(bounds);
  implicit->SetMaximumDistance( thickness );
  implicit->SetSampleDimensions(dims);
  implicit->AdjustBoundsOff();
  implicit->CappingOff();
  implicit->SetCapValue(thickness);
  implicit->SetProcessModeToPerVoxel();
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

  int incx,incy,incz;
  sp->GetIncrements(incx,incy,incz);

  for (vtkIdType i=0; i<npoints; ++i) {
    if (inout[i]) continue;

    double pt[3];
    sp->GetPoint(i,pt);
    const int _inout = obb->InsideOrOutside(pt);
    if ( _inout == 0 ) std::cerr << "Warning: voxel could not be classified\n";

    filo.push(i);
    while (!filo.empty()) {
      const vtkIdType j = filo.front();
      filo.pop();
      if (inout[j]) continue;

      inout[j] = _inout;
      const float dist = tab[j];
      tab[j] *= _inout;

      // mark its neighbors if distance > 1 (triangle inequality)
      if ( dist>1 ) {
        const int x = j%incy;
        const int y = (j%incz) / incy;
        const int z = j/incz;
        if (x>0 && !inout[j-incx]) filo.push(j-incx);
        if (y>0 && !inout[j-incy]) filo.push(j-incy);
        if (z>0 && !inout[j-incz]) filo.push(j-incz);
        if (x<dims[0]-1 && !inout[j+incx]) filo.push(j+incx);
        if (y<dims[1]-1 && !inout[j+incy]) filo.push(j+incy);
        if (z<dims[2]-1 && !inout[j+incz]) filo.push(j+incz);
      }
    }
  } // next i

  for (vtkIdType j=0; j<npoints; ++j) {
    const int x = j%incy;
    const int y = (j%incz) / incy;
    const int z = j/incz;

    MRIsetVoxVal(vol, x,y,z,0, tab[j] );
  }

}

MRI*
CreateHemiMask(MRI* dpial,
               MRI* dwhite,
               const unsigned char lblWhite,
               const unsigned char lblRibbon,
               const unsigned char lblBackground) {
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
           x<width; ++x) {
        if ( MRIFvox(dwhite,x,y,z) > 0 )
          MRIsetVoxVal(mri, x,y,z,0, lblWhite);
        else if ( MRIFvox(dpial,x,y,z) > 0 )
          MRIsetVoxVal(mri, x,y,z,0, lblRibbon);
        else
          MRIsetVoxVal(mri,x,y,z,0, lblBackground);
      } // next x,y,z

  return mri;
}

/*

implicit assumption the masks do not overlap
if this is the case, a message will be printed

*/
MRI* CombineMasks(MRI* maskOne,
                  MRI* maskTwo,
                  const unsigned char lblBackground) {
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
           x<width; ++x) {
        voxOne = MRIvox(maskOne,x,y,z);
        voxTwo = MRIvox(maskTwo,x,y,z);
        if ( voxOne!=lblBackground && voxTwo!=lblBackground ) {
          // overlap
          ++overlap;
          MRIsetVoxVal(mri,x,y,z,0,voxOne);
        } else if ( voxOne == lblBackground )
          MRIsetVoxVal(mri,x,y,z,0,voxTwo);
        else if ( voxTwo == lblBackground )
          MRIsetVoxVal(mri,x,y,z,0,voxOne);
      } // next x,y,z

  std::cout << " hemi masks overlap voxels = " << overlap << std::endl;
  return mri;
}

void
ConvertSurfaceVertexCoordinates(MRI_SURFACE* mris,
                                MRI* vol) {
  double cx, cy, cz;
  Real vx, vy, vz;

  VERTEX* pvtx = &( mris->vertices[0] );
  unsigned int nvertices = (unsigned int)mris->nvertices;

  for ( unsigned int ui=0;
        ui < nvertices;
        ++ui, ++pvtx ) {
    cx = pvtx->x;
    cy = pvtx->y;
    cz = pvtx->z;

    MRIsurfaceRASToVoxel( vol,
                          cx, cy, cz,
                          &vx, &vy, &vz);

    pvtx->x = vx;
    pvtx->y = vy;
    pvtx->z = vz;
  } // next ui, pvtx
}

MRI* FilterLabel(MRI* initialMask,
                 const unsigned char lbl) {
  MRI* mri = MRIalloc( initialMask->width,
                       initialMask->height,
                       initialMask->depth,
                       MRI_UCHAR);

  for (unsigned int z(0), depth(initialMask->depth);
       z<depth; ++z)
    for (unsigned int y(0), height(initialMask->height);
         y<height; ++y)
      for (unsigned int x(0), width(initialMask->width);
           x<width; ++x) {
        if ( MRIvox(initialMask,x,y,z) == lbl )
          MRIsetVoxVal(mri, x,y,z,0, (unsigned char)(1) );
        else
          MRIsetVoxVal(mri, x,y,z,0, (unsigned char)(0) );
      } // next x,y,z

  return mri;
}
