#include <stdlib.h>
#include <fstream>
#include "string_fixed.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include "SurfaceCollection.h"
#include "DataManager.h"
extern "C" {
#include "mrisurf.h"
}
#include "Scuba-impl.h"

#define Assert(x,e)   if(!(x)) { throw e; }

using namespace std;

char* Progname = "test_SurfaceCollection";

int main ( int argc, char** argv ) {

  cerr << "Beginning test..." << endl;

  try {

    string fnMRIS = "test_data/lh.white";
    SurfaceCollection surf;
    surf.SetSurfaceFileName( fnMRIS );
    MRIS* mris = surf.GetMRIS();
    
    DataManager dataMgr = DataManager::GetManager();
    MRISLoader mrisLoader = dataMgr.GetMRISLoader();
    Assert( 1 == mrisLoader.CountLoaded(), 
	    logic_error( "CountLoaded didn't return 1" ) );
    Assert( 1 == mrisLoader.CountReferences(mris),
	    logic_error( "CountReferences didn't return 1" ) );

    char* fnMRISC = strdup( fnMRIS.c_str() );
    MRIS* mrisComp = MRISread( fnMRISC );
    
    Assert( (mrisComp->nvertices == mris->nvertices),
	    logic_error( "Inequal nvertices" ));
    for( int nVertex = 0; nVertex < mris->nvertices; nVertex++ ) {
      VERTEX* v = &(mris->vertices[nVertex]);
      VERTEX* vComp = &(mrisComp->vertices[nVertex]);
      Assert( (v->x == vComp->x && v->y == vComp->y && v->z == vComp->z),
	      logic_error( "vertex comparison failed" ) );
    }
    
    MRISfree( &mrisComp );

  }
  catch( logic_error& e ) {
    cerr << "failed with exception: " << e.what() << endl;
    exit( 1 );
  }
  catch(...) {
    cerr << "failed." << endl;
    exit( 1 );
  }

  cerr << "Success" << endl;
  
  exit( 0 );
}
