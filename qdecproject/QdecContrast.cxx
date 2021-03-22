/**
 * @brief Stores a GLM contrast vector
 *
 * Stores a GLM contrast vector associated with a particular design (being
 * based on user selected factors, formulating a hypothesis to test).
 * Stores the name associated with it (to go in the .fsgd file), the
 * human-readable question (for use by GUI), and the vector itself (as an
 * int array).
 */
/*
 * Original Author: Nick Schmansky
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

#include <string.h>
#include <errno.h>
#include <assert.h>
#include <sys/stat.h>

#include "QdecContrast.h"
#include <stdio.h> // printf


// Constructors/Destructors
//

QdecContrast::QdecContrast ( vector< double > iaVector,
                             string isName,
                             string isQuestion )
{
  assert( iaVector.size() );
  this->maVector = iaVector;
  this->msName = isName;
  this->msQuestion = isQuestion;
}

QdecContrast::~QdecContrast ( )
{ }

//
// Methods
//


/**
 * @return string
 */
string QdecContrast::GetName ( )
{
  return this->msName;
}


/**
 * @return string
 */
string QdecContrast::GetQuestion ( )
{
  return this->msQuestion;
}


/**
 * @return string
 */
string QdecContrast::GetContrastStr ( )
{
  string contrast = "";
  for (unsigned int i=0; i < this->maVector.size();)
  {
    char tmpstr[1000];
    sprintf(tmpstr,"% 2.3f",this->maVector[i]);
    contrast += strdup(tmpstr);
    if (++i < this->maVector.size()) contrast += "  ";
  }
  contrast += ";\n";

  return contrast;
}


/**
 * Writes the contrast vector to a .mat file, which is readable by matlab, and
 * mri_glmfit.
 * @return int
 * @param string ifnWorkingDir
 */
int QdecContrast::WriteDotMatFile ( string ifnWorkingDir )
{
  string dirName = ifnWorkingDir + "/contrasts/";
  int err = mkdir( dirName.c_str(), 0777);
  if( err != 0 && errno != EEXIST )
  {
    fprintf( stderr,
             "ERROR: QdecContrast::WriteDotMatFile: "
             "could not create directory %s\n",
             dirName.c_str());
    return(-1);
  }

  this->mfnDotMatFileName = dirName; 
  this->mfnDotMatFileName += this->GetName();
  this->mfnDotMatFileName += ".mat";

  FILE* fp = fopen( this->mfnDotMatFileName.c_str(), "w");
  if( ! fp )
  {
    fprintf( stderr,
             "ERROR: QdecContrast::WriteDotMatFile: "
             "could not create file %s\n",
             this->mfnDotMatFileName.c_str());
    return(-2);
  }

  for(unsigned int i=0; i < this->maVector.size(); i++)
  {
    fprintf( fp, "%+4.5f ", this->maVector[i] );
  }
  fprintf( fp, "\n" );
  fclose( fp );

  return 0;
}


/**
 * @return string
 */
string QdecContrast::GetDotMatFileName ( )
{
  return this->mfnDotMatFileName;
}

