#define BOOST_TEST_MODULE My Test
#include <boost/test/included/unit_test.hpp>

#include <cstdio>

#include "mri.h"
#include "gcamorph.h"
#include "gca.h"

#include "gcalinearnode.hpp"
#include "gcalinearprior.hpp"

#include "chronometer.hpp"

void LoadData(MRI **mri, GCAM **gcam, GCA **gca, const int i) {
  const int MaxLength = 256;
  char fNameMRI[MaxLength], fNameGCAM[MaxLength], fNameGCA[MaxLength];

  snprintf(fNameMRI,MaxLength-1, "GCAMcompLabels%03i.mgz",i);
  snprintf(fNameGCAM,MaxLength-1, "GCAMcompLabels%03i.m3z",i);
  snprintf(fNameGCA,MaxLength-1, "GCAMcompLabels%03i.gca",i);

  *mri = MRIread(fNameMRI);
  *gcam = GCAMread(fNameGCAM);
  *gca = GCAread(fNameGCA);
}

BOOST_AUTO_TEST_CASE(LinearNode)
{
  SciGPU::Utilities::Chronometer tBasicCall, tTestCall;
  SciGPU::Utilities::Chronometer tExhumeNode;
  SciGPU::Utilities::Chronometer tExhumePrior;

  MRI *mri, *mriTest;
  GCAM *gcam, *gcamTest;
  GCA *gca, *gcaTest;

  const int i = 0;

  int nChanged, nChangedTest = 0;

  LoadData(&mri, &gcam, &gca, i);
  LoadData(&mriTest, &gcamTest, &gcaTest, i);

  BOOST_REQUIRE( mri != NULL );
  BOOST_REQUIRE( gcam != NULL );
  BOOST_REQUIRE( gca != NULL );

  BOOST_REQUIRE( mriTest != NULL );
  BOOST_REQUIRE( gcamTest != NULL );
  BOOST_REQUIRE( gcaTest != NULL );

  gcam->gca = gca;

  // Reference call
  tBasicCall.Start();
  nChanged = GCAMcomputeLabels(mri, gcam);
  tBasicCall.Stop();

  // Test call
  Freesurfer::GCAlinearNode gcaLN;
  Freesurfer::GCAlinearPrior gcaLP;

  BOOST_TEST_MESSAGE("Exhuming Nodes");
  tExhumeNode.Start();
  gcaLN.Exhume(gcaTest);
  tExhumeNode.Stop();

  BOOST_TEST_MESSAGE("Exhuming Priors");
  tExhumePrior.Start();
  gcaLP.Exhume(gcaTest);
  tExhumePrior.Stop();
  
  BOOST_TEST_MESSAGE("Running call");
  tTestCall.Start();
  // nChangedTest = GCAMcomputeLabelsLinearCPU(mriTest, gcamTest, gcaLN, gcaLP);
  tTestCall.Stop();
  
  BOOST_TEST_MESSAGE( "tBasicCall  : " << tBasicCall );
  BOOST_TEST_MESSAGE( "tExhumeNode : " << tExhumeNode );
  BOOST_TEST_MESSAGE( "tExhumePrior: " << tExhumePrior );
  BOOST_TEST_MESSAGE( "tTestCall   : " << tTestCall );

  BOOST_CHECK( nChanged == nChangedTest );

  // Clean up
  MRIfree( &mri );
  GCAMfree( &gcam );
  GCAfree( &gca );

  MRIfree( &mriTest );
  GCAMfree( &gcamTest );
  GCAfree( &gcaTest );
}
