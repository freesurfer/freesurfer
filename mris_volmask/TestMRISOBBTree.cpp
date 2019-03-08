#include "MRISOBBTree.h"
#include <string>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/ui/text/TestRunner.h>
 
#include "fsenv.h"

const char *Progname;
// the setting of SUBJECTS_DIR where bert resides is essential to this test.
class TestMRISOBBTree : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE ( TestMRISOBBTree );

    CPPUNIT_TEST( TestComputeOBB ); //tests whether a single OBB Node is computed correctly
    CPPUNIT_TEST( TestBuildTree ); //tests whether tree gets built recursively
    CPPUNIT_TEST( TestPointInclusion ); //tests a set of points for inclusion

  CPPUNIT_TEST_SUITE_END();
  public:
    TestMRISOBBTree()
    {
      std::string fullpath1( FSENVgetSUBJECTS_DIR() );
      std::string fullpath2 ( fullpath1 );
      fullpath1.append("/bert/surf/lh.white");
      fullpath2.append("/bert/mri/aseg.mgz");

      mris = MRISread(fullpath1.c_str() );	
      tmp = MRIread(fullpath2.c_str() );
      OBBTree = new MRISOBBTree(mris);
      Math::ConvertSurfaceRASToVoxel(mris, tmp);
    }

    ~TestMRISOBBTree()
    {
      delete OBBTree;
    }

    void TestComputeOBB()
    {
      std::cerr << "TestComputeOBB (tests whether the root node created ( the biggest OBB) has the given orientation and size)." <<  std::endl;
      OBBNode *pOBB = new OBBNode;
      std::vector<int> faceidlist;
      faceidlist.reserve(mris->nfaces);
      for (int i=0; i< mris->nfaces; i++)
        faceidlist.push_back(i);
      OBBTree->ComputeOBB(faceidlist, pOBB);
      double delta = 0.001;
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pOBB->corner[0]  , 135.644  , delta);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pOBB->corner[1]  , 175.577  , delta);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pOBB->corner[2]  , 27.0459  , delta);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pOBB->axes[0][0] , -4.46    , delta);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pOBB->axes[1][1] , -106.193 , delta);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(pOBB->axes[2][2] , 1.0086   , delta);
      delete pOBB;
    }

    void TestBuildTree()
    {
      std::cerr << "TestBuildTree (tests the # of OBBNodes created for a 12 level tree). " << std::endl;
      OBBTree->SetMaxLevel(12);
      OBBTree->ConstructTree();
      CPPUNIT_ASSERT_EQUAL(OBBTree->GetNodesCount(), 8191);
    }

    void TestPointInclusion()
    {
      std::cerr << "TestPointInclusion ( tests a bunch of points whether they are inside the mesh)" << std::endl;
      OBBTree->ConstructTree();
      CPPUNIT_ASSERT_EQUAL( OBBTree->PointInclusionTest(152, 106, 128), 1);  // defn inside
      CPPUNIT_ASSERT_EQUAL( OBBTree->PointInclusionTest(179, 117, 128), 1);  // borderline inside

      CPPUNIT_ASSERT_EQUAL( OBBTree->PointInclusionTest(177, 117, 128), -1); // borderline outside
      CPPUNIT_ASSERT_EQUAL( OBBTree->PointInclusionTest(0, 0, 0), -1);       // definitely outside
      CPPUNIT_ASSERT_EQUAL( OBBTree->PointInclusionTest(300, 300, 300), -1); // embarassingly outside
    }

  private:
    MRISOBBTree* OBBTree;
    MRIS* mris;
    MRI* tmp;

};


int main ( int argc, char** argv ) {

  // this is needed by the freesurfer utils library
  Progname = argv[0];
  
  const int SUCCESS = 0;
  const int FAIL = 1;

  CppUnit::TextUi::TestRunner runner;
  runner.addTest( TestMRISOBBTree::suite() );

  if ( runner.run() ) {
    exit ( SUCCESS );
  }

  exit( FAIL );
}
