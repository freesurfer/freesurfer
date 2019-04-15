#include "MRISDistanceField.h"
#include <string>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
 
#include "fsenv.h"

const char *Progname;
// the setting of SUBJECTS_DIR where bert resides is essential to this test.
class TestMRISDistanceField : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE ( TestMRISDistanceField );

    CPPUNIT_TEST( TestPointDist ); //tests a set of points for inclusion

  CPPUNIT_TEST_SUITE_END();
  public:
    TestMRISDistanceField()
    {
      std::string fullpath1( "/Users/krish/subjects/bert/surf/lh.white" );
      std::string fullpath2 ( "/Users/krish/subjects/bert/mri/aseg.mgz" );

      mris = MRISread(fullpath1.c_str() );	
      tmp = MRIread(fullpath2.c_str() );
      mridist = MRIcloneDifferentType(tmp, MRI_FLOAT);
      
      Math::ConvertSurfaceRASToVoxel(mris, tmp);
      distrunner = new MRISDistanceField(mris, mridist);
      distrunner->Generate();
    }

    ~TestMRISDistanceField()
    {
      delete distrunner;
      MRISfree(&mris);
      MRIfree(&tmp);
      MRIfree(&mridist);
    }


    void TestPointDist()
    {

    }

  private:
    MRIS *mris;
    MRI *tmp, *mridist;
    MRISDistanceField *distrunner;
};


int main ( int argc, char** argv ) {

  // this is needed by the freesurfer utils library
  Progname = argv[0];
  
  const int SUCCESS = 0;
  const int FAIL = 1;

  CppUnit::TextUi::TestRunner runner;
  runner.addTest( TestMRISDistanceField::suite() );

  if ( runner.run() ) {
    exit ( SUCCESS );
  }

  exit( FAIL );
}
