/**
 * @brief Container for a text-input file of data to QDEC.
 *
 * Implements loading/saving the white-space delimited data file containing
 * the list of subjects with their discrete and continuous factors.
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

#include <fstream>
#include <iostream>
#include <math.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "QdecDataTable.h"

// Constructors/Destructors
//

QdecDataTable::QdecDataTable() { mfnFileName = "not loaded"; }

QdecDataTable::~QdecDataTable() {
  while (this->mFactors.size() != 0) {
    delete this->mFactors.back();
    this->mFactors.pop_back();
  }
  while (this->mSubjects.size() != 0) {
    delete this->mSubjects.back();
    this->mSubjects.pop_back();
  }
}

//
// Methods
//

/**
 * Load white-space delimited file containing subject ids and their
 * discrete and continuous factors.  Upon exit, the mFactors member variable
 * contains a list (vector) of the factors read from the first line stored in
 * QdecFactor objects, and the mSubjects member variable contains a list
 * (vector) of all the subject data in QdecSubject objects (within which
 * each data item is stored in a QdecFactor object).
 * @return int
 * @param  isFileName
 * @param  osNewSubjDir
 * @param  isFsIdColName
 */
int QdecDataTable::Load(const char *isFileName, char *osNewSubjDir,
                        const char *isFsIdColName) {
  size_t tmpstrMaxSize = 200000; // maximum size of one line in the file
  char * tmpstr        = (char *)malloc(tmpstrMaxSize);
  assert(tmpstr);

  if (nullptr == isFileName) {
    std::cerr << "ERROR: QdecDataTable::Load: NULL filename!" << std::endl;
    return (-1);
  }

  // the file could contain a line beginning with SUBJECTS_DIR, to specify
  // where this data is located.  otherwise, assume enviro SUBJECTS_DIR
  if (NULL != osNewSubjDir)
    osNewSubjDir[0] = 0; // assume no new subj dir

  // delete any prior loaded data
  while (this->mFactors.size() != 0) {
    delete this->mFactors.back();
    this->mFactors.pop_back();
  }
  while (this->mSubjects.size() != 0) {
    delete this->mSubjects.back();
    this->mSubjects.pop_back();
  }

  this->mfnFileName = isFileName;

  std::cout << "\nLoading data table " << isFileName << "...\n";

  std::ifstream ifsDatFile;
  ifsDatFile.open(isFileName);
  if (!ifsDatFile.is_open()) {
    std::cerr << "ERROR: could not open " << isFileName << std::endl;
    return (-1);
  }

  // the file might be in Mac text format, which uses CR as the end-of-line,
  // so try getting a line, and if we reach end-of-file without getting a line,
  // then try using the Mac delimiter 'CR' (instead of the default 'LF')
#define UNIX_EOL 0x0A
#define MAC_EOL  0x0D
  char eol_delim = UNIX_EOL;
  ifsDatFile.getline(tmpstr, tmpstrMaxSize, eol_delim);
  if (ifsDatFile.eof()) {
    std::cout
        << "QdecDataTable::Load: will attempt to open as Mac text file...\n";
    eol_delim = MAC_EOL;
  }
  ifsDatFile.clear();
  ifsDatFile.seekg(0); // rewind

#undef WHITESPC
#define WHITESPC " ,\"\t\n\r"

  // Attempt to load the first non-commented line of the file,
  // which may fail if it cant detect end-of-line
  tmpstr[0] = '#';
  while (tmpstr[0] == '#') // ignore lines beginning with #
  {
    ifsDatFile.getline(tmpstr, tmpstrMaxSize, eol_delim);
    if (ifsDatFile.fail()) {
      std::cerr << "ERROR: QdecDataTable::Load failed to load first line of "
                << isFileName << "!\n";
      ifsDatFile.close();
      return (-1);
    }

    // while we're looking for that first line of factor names, look for
    // a line with SUBJECTS_DIR as the first string, and set it to whatever
    // follows it
    if (strncmp(tmpstr, "SUBJECTS_DIR", 12) == 0) {
      char *token = strtok(&tmpstr[13], WHITESPC);
      strcpy(osNewSubjDir, token);
      std::cout << "Setting SUBJECTS_DIR to " << osNewSubjDir << std::endl;
      tmpstr[0] = '#'; // continue trying find the first line with factor info
    }
  }

  /*
   * Count the number of columns in the first line from the file,
   * and, at the same time, find the column which contains the subject
   * names.  that column must have a special name (fsid, ID, Id...see below)
   * and must be column zero (the first column).
   * note: optionally the expected name of the fsid column can be passed-in
   * as a parameter (isFsIdColName)
   */
  int   ncols   = 0;
  int   fsidcol = -1;
  char *token   = strtok(tmpstr, WHITESPC); // get first token in this line
  while (token != NULL) {
    //cout << token << std::endl;
    if (!strcmp(token, "fsid"))
      fsidcol = ncols;
    else if (!strcmp(token, "ID"))
      fsidcol = ncols;
    else if (!strcmp(token, "Id"))
      fsidcol = ncols;
    else if (!strcmp(token, "subject_id"))
      fsidcol = ncols;
    else if (!strcmp(token, "subjid"))
      fsidcol = ncols;
    else if (!strcmp(token, "subject"))
      fsidcol = ncols;
    else if (!strcmp(token, "Subject"))
      fsidcol = ncols;
    else if (!strcmp(token, "Measure:volume"))
      fsidcol = ncols;
    else if (isFsIdColName && !strcmp(token, isFsIdColName))
      fsidcol = ncols;
    ncols++;
    token = strtok(NULL, WHITESPC); // get next token in this line
  }
  if (fsidcol == -1) {
    std::cerr << "ERROR: QdecDataTable::Load could not find column named "
                 "'fsid' or 'ID' in "
              << isFileName << "!\n";
    ifsDatFile.close();
    return (-1);
  }
  if (fsidcol != 0) {
    std::cerr << "ERROR: QdecDataTable::Load did not find a column named "
                 "'fsid', 'ID', or 'Subject' in the first column of "
              << isFileName << "!\n";
    ifsDatFile.close();
    return (-1);
  }
  int nFactors = ncols - 1;

  /*
   * Count the number of input rows (subjects)
   */
  int nInputs    = 0;
  int getLineRet = 0;
  while ((getLineRet =
              ifsDatFile.getline(tmpstr, tmpstrMaxSize, eol_delim).good())) {
    if ((strlen(tmpstr) > 2) && (tmpstr[0] != '#'))
      nInputs++;
  }

  // --------------------------------------------------
  std::cout << "Number of columns:  " << ncols << std::endl;
  std::cout << "fsid column:        " << fsidcol + 1 << std::endl;
  std::cout << "Number of factors:  " << nFactors << std::endl;
  std::cout << "Number of subjects: " << nInputs << std::endl;

  if (nInputs <= 0) {
    std::cerr << "ERROR: QdecDataTable::Load: number of subjects = " << nInputs
              << std::endl;
    ifsDatFile.close();
    return (-1);
  }
  if ((getLineRet == 0) && (strlen(tmpstr) > 2)) {
    std::cout << "ERROR: QdecDataTable::Load: problem parsing file "
              << isFileName << std::endl;
    std::cout << "This line did not appear to end with a newline: " << std::endl
              << tmpstr << std::endl;
    std::cout << "In your editor, place the cursor at the end of that"
              << std::endl
              << "line and press Enter once, then save the file." << std::endl;
    ifsDatFile.close();
    return (-1);
  }

  /*
   * Look for a file named ignore.factors which can list any factors the
   * user wishes to ignore (stuff like alternate subject ids which are 
   * neither continuous nor discrete data)
   */
  std::vector<std::string> sIgnoreFactors;
  std::string              fnDataTable = isFileName;
  std::string              fnPath;
  std::string::size_type   nPreLastSlash = fnDataTable.rfind('/');
  if (std::string::npos != nPreLastSlash)
    fnPath = fnDataTable.substr(0, nPreLastSlash);
  else
    fnPath = ".";
  std::stringstream fnIgnore;
  fnIgnore << fnPath << "/ignore.factors";
  std::ifstream ifsIgnoreFile(fnIgnore.str().c_str(), std::ios::in);
  if (ifsIgnoreFile.good()) {
    std::cout << "Reading factors to ignore from config file "
              << fnIgnore.str().c_str() << std::endl;
    char tmpstr2[2048];
    while (ifsIgnoreFile.getline(tmpstr2, 2000, UNIX_EOL).good()) {
      if (strlen(tmpstr2) >= 1) {
        std::cout << "\t" << tmpstr2 << std::endl;
        ;
        sIgnoreFactors.push_back(strdup(tmpstr2));
      }
    }
    ifsIgnoreFile.close();
    std::cout << "done." << std::endl;
  }

  /*
   * Look for a file named ordinal.factors which can list the continuous
   * factors that are known to be integers (not floating point).
   * useful only when saving the project file so that an integer doesnt get
   * written as a float.
   */
  std::vector<std::string> sOrdinalFactors;
  nPreLastSlash = fnDataTable.rfind('/');
  if (std::string::npos != nPreLastSlash)
    fnPath = fnDataTable.substr(0, nPreLastSlash);
  else
    fnPath = ".";
  std::stringstream fnOrdinal;
  fnOrdinal << fnPath << "/ordinal.factors";
  std::ifstream ifsOrdinalFile(fnOrdinal.str().c_str(), std::ios::in);
  if (ifsOrdinalFile.good()) {
    std::cout << "Reading factors to consider as ordinal from config file "
              << fnOrdinal.str().c_str() << std::endl;
    char tmpstr2[2048];
    while (ifsOrdinalFile.getline(tmpstr2, 2000, UNIX_EOL).good()) {
      if (strlen(tmpstr2) >= 1) {
        std::cout << "\t" << tmpstr2 << std::endl;
        sOrdinalFactors.push_back(strdup(tmpstr2));
      }
    }
    ifsOrdinalFile.close();
    std::cout << "done." << std::endl;
  }

  /*
   * read-in the factor names from the first non-commented line of input
   */
  ifsDatFile.clear();
  ifsDatFile.seekg(0); // rewind
  tmpstr[0] = '#';
  while (tmpstr[0] == '#') // ignore lines beginning with # and SUBJECTS_DIR
  {
    ifsDatFile.getline(tmpstr, tmpstrMaxSize, eol_delim);
    if (ifsDatFile.fail() || (NULL == tmpstr)) {
      std::cerr << "ERROR2: QdecDataTable::Load failed to load first line of "
                << isFileName << std::endl;
      ifsDatFile.close();
      return (-1);
    }
    if (strncmp(tmpstr, "SUBJECTS_DIR", 12) == 0)
      tmpstr[0] = '#'; // continue
  }
  token = strtok(tmpstr, WHITESPC);
  if (NULL == token) {
    std::cerr << "ERROR: QdecDataTable::Load failed to tokenize string: "
              << tmpstr << std::endl;
    ifsDatFile.close();
    return (-1);
  } //else std::cout << "token: %s\n",token);

  int nthfactor = 0;
  while ((nthfactor < nFactors) && (token)) {
    if (fsidcol == nthfactor) {
      // skip-past the fsid column
      token = strtok(NULL, WHITESPC);
      if (NULL == token) {
        std::cerr << "ERROR2: QdecDataTable::Load failed to tokenize string: "
                  << tmpstr << std::endl;
        ifsDatFile.close();
        return (-1);
      } //else std::cout << "token: %s\n",token);
    }

    char factor[2048];
    strncpy(factor, token, sizeof(factor) - 1);
    //cout << "factor: " << factor << std::endl;

    // determine if this factor should be ignored (by comparing against
    // what we may have read from the ignore.factors file parsed earlier)
    bool bIgnore = false;
    for (unsigned int n = 0; n < sIgnoreFactors.size(); n++) {
      //cout << this->mIgnoreFactors[n].c_str() << "  " << factor << std::endl;
      if (strcmp(sIgnoreFactors[n].c_str(), factor) == 0) {
        bIgnore = true;
      }
    }

    // if there exists a file called 'factor'.levels, where 'factor' is the
    // token read from the line, then it is a discrete factor, in which
    // case we'll read its valid levels, otherwise, assume its continuous or
    // ignore it

    // Extract the path of the data table from the data table file name.
    std::string            fnDataTable = isFileName;
    std::string            fnPath;
    std::string::size_type nPreLastSlash = fnDataTable.rfind('/');
    if (std::string::npos != nPreLastSlash)
      fnPath = fnDataTable.substr(0, nPreLastSlash);
    else
      fnPath = ".";

    // Build the levels file name.
    std::stringstream fnLevels;
    fnLevels << fnPath << "/" << factor << ".levels";
    // std::cout << fnLevels.str().c_str() << std::endl;

    // Try to open the levels file (implicitly means this is discrete)
    std::ifstream ifsLevelFile(fnLevels.str().c_str(), std::ios::in);
    if (ifsLevelFile.good()) {
      QdecFactor *qf = new QdecFactor(strdup(factor), // name
                                      QdecFactor::qdecDiscreteFactorType);
      std::cout << "Reading discrete factor levels from config file "
                << fnLevels.str().c_str() << std::endl;
      int  levelCount = 0;
      char tmpstr2[1000];
      while (ifsLevelFile.getline(tmpstr2, 1000, UNIX_EOL).good()) {
        if (strlen(tmpstr2) >= 1) {
          std::cout << "\t" << tmpstr2 << std::endl;
          qf->AddLevelName(tmpstr2);
          levelCount++;
        }
      }
      qf->SetHaveDotLevelsFile();
      ifsLevelFile.close();
      std::cout << "done." << std::endl;

      // error check: cannot have just one level (or no levels)
      if (levelCount < 2) {
        std::cerr << "ERROR: " << fnLevels.str().c_str() << " is invalid."
                  << " Must have at least two levels." << std::endl;
        ifsDatFile.close();
        return (-1);
      }
      // ok, so add it to our factors storage
      this->mFactors.push_back(qf);
    } else if (bIgnore) {
      // factor to ignore
      QdecFactor *qf = new QdecFactor(strdup(factor), // name
                                      QdecFactor::qdecIgnoreType);
      this->mFactors.push_back(qf);
    } else // must be continuous
    {
      QdecFactor *qf = new QdecFactor(strdup(factor), // name
                                      QdecFactor::qdecContinuousFactorType);

      // determine if this factor is ordinal (by comparing against
      // what we may have read from the ordinal.factors file parsed earlier)
      for (unsigned int n = 0; n < sOrdinalFactors.size(); n++) {
        if (strcmp(sOrdinalFactors[n].c_str(), factor) == 0) {
          qf->SetOrdinal();
          break;
        }
      }

      this->mFactors.push_back(qf);
    }

    nthfactor++;

    token = strtok(NULL, WHITESPC); // get next string in this line
    if (!token)
      break;
  } // end while (nthfactor < nFactors)

  // sanity checks
  if (nthfactor == 0) {
    std::cerr
        << "ERROR: QdecDataTable::Load failed to read any factors from the"
           " first line of "
        << isFileName << std::endl;
    ifsDatFile.close();
    return (-1);
  }
  if (nthfactor != nFactors) {
    std::cerr
        << "ERROR: QdecDataTable::Load failed to read all factors from the"
           " first line of "
        << isFileName << std::endl;
    ifsDatFile.close();
    return (-1);
  }

  /*
   * read-in each row of subject data from the data table
   */
  int numExpContFactors = 0; // these are used for sanity-checking
  int numExpDiscFactors = 0;
  for (int nthInput = 0; nthInput < nInputs; nthInput++) {
    std::string               subj_id;
    std::vector<QdecFactor *> theFactors;
    tmpstr[0] = '#';
    while (tmpstr[0] == '#') // ignore lines beginning with #
    {
      ifsDatFile.getline(tmpstr, tmpstrMaxSize, eol_delim);
      if (ifsDatFile.fail() || (NULL == tmpstr)) {
        std::cerr << "ERROR: QdecDataTable::Load failed to load line "
                  << nthInput + 2 << " of " << isFileName << std::endl;
        ifsDatFile.close();
        return (-1);
      }
    }
    token = strtok(tmpstr, WHITESPC); // a token is each column item
    if (NULL == token) {
      std::cerr << "ERROR3: QdecDataTable::Load failed to tokenize string: "
                << tmpstr << std::endl
                << "on line " << nthInput + 2 << " of " << isFileName
                << std::endl;
      ifsDatFile.close();
      return (-1);
    } // else std::cout << "token: " << token << std::endl;

    int numContFactors = 0; // these are used for sanity-checking
    int numDiscFactors = 0;
    int errs           = 0;
    nthfactor          = 0;
    for (int nthcol = 0; nthcol < ncols; nthcol++) {
      if (nthcol == fsidcol) // get subject id
      {
        subj_id = strdup(token);
      } else // get factor data
      {
        // start by assuming its continuous by trying to convert to a double
        double dtmp    = 0.0;
        int    retCode = sscanf(token, "%lf", &dtmp);
        if (retCode == 1 && !this->mFactors[nthfactor]->IsDiscrete() &&
            !this->mFactors[nthfactor]->HaveDotLevelsFile() &&
            !this->mFactors[nthfactor]
                 ->Ignore()) { // yes!  its a continuous factor
          //cout << "\t" << nthInput << ": " << dtmp << std::endl;
          QdecFactor *qf = new QdecFactor(
              this->mFactors[nthfactor]->GetFactorName().c_str(),
              QdecFactor::qdecContinuousFactorType, dtmp /* value */);
          if (this->mFactors[nthfactor]->IsOrdinal())
            qf->SetOrdinal();
          theFactors.push_back(qf); // save this factor data
          // for later sanity-checking:
          if (nthInput == 0)
            numExpContFactors++;
          numContFactors++;
        } else if (this->mFactors[nthfactor]->Ignore()) {
          // this column is being ignored (ie, its not a factor, but
          // we still need to save the data in case we save the file later)
          QdecFactor *qf =
              new QdecFactor(this->mFactors[nthfactor]->GetFactorName().c_str(),
                             QdecFactor::qdecIgnoreType,
                             (const char *)strdup(token) /* value */);
          theFactors.push_back(qf);
        } else // it must be a discrete factor
        {
          // if discrete, then check that its valid (a known level name, as
          // read from a factor.levels file that user optionally created)
          if (this->mFactors[nthfactor]->IsDiscrete() &&
              this->mFactors[nthfactor]->HaveDotLevelsFile() &&
              !this->mFactors[nthfactor]->ValidLevelName(token)) {
            std::cerr << std::endl
                      << "ERROR: Subject " << subj_id
                      << " has an invalid level '" << strdup(token)
                      << "' in the "
                      << this->mFactors[nthfactor]->GetFactorName().c_str()
                      << " column" << std::endl;
            std::cerr << "INFO: If '"
                      << this->mFactors[nthfactor]->GetFactorName().c_str()
                      << "' is a discrete factor, then create a file"
                      << std::endl
                      << "named '"
                      << this->mFactors[nthfactor]->GetFactorName().c_str()
                      << ".levels' containing the valid factor names"
                      << std::endl
                      << "one per line." << std::endl;
            ifsDatFile.close();
            return (-1);
          } else // we dont know about this discrete factor, so update mFactors
          {
            this->mFactors[nthfactor]->SetDiscrete();
            this->mFactors[nthfactor]->AddLevelName(token);
          }
          // and save-away this subjects discrete data
          QdecFactor *qf =
              new QdecFactor(this->mFactors[nthfactor]->GetFactorName().c_str(),
                             QdecFactor::qdecDiscreteFactorType,
                             (const char *)strdup(token), /* value */
                             this->mFactors[nthfactor]->GetLevelNames());
          theFactors.push_back(qf); // save this factor data
          // for later sanity-checking
          if (nthInput == 0)
            numExpDiscFactors++;
          numDiscFactors++;
        }
        nthfactor++;
      }

      token = strtok(NULL, WHITESPC); // get next factor in this line
      if (!token)
        break;

    } // end for (int nthcol = 0; nthcol < ncols; nthcol++)

    // sanity-check, if this row has more or less than expected data
    if (numExpContFactors != numContFactors) {
      std::cerr << "\nERROR: Subject " << subj_id
                << " seems to have a mismatched "
                   "number of continuous factors (expected "
                << numExpContFactors << ", found " << numContFactors << ")"
                << std::endl;
      errs++;
    }
    if (numExpDiscFactors != numDiscFactors) {
      std::cerr << "\nERROR: Subject " << subj_id
                << " seems to have a mismatched "
                   "number of discrete factors (expected "
                << numExpDiscFactors << ", found " << numDiscFactors << ")"
                << std::endl;
      errs++;
    }

    // check if this subj_id is already in the table (indicating an error)
    for (unsigned int m = 0; m < this->mSubjects.size(); m++) {
      if (subj_id == this->mSubjects[m]->GetId()) {
        std::cerr << "\nERROR: Subject " << subj_id << " already exists in "
                  << isFileName << std::endl;
        errs++;
      }
    }

    // add this subject data to our Subjects list (if no errors found)
    assert(theFactors.size());
    if (errs == 0) {
      QdecSubject *qsubj = new QdecSubject(subj_id.c_str(), theFactors);
      this->mSubjects.push_back(qsubj);
    }

    // continue to next subject...
  } // end for (int nthInput = 0; nthInput < nInputs; nthInput++)

  ifsDatFile.close();

  // check for discrete factors having only one level
  for (int n = 0; n < nFactors; n++) {
    if (this->mFactors[n]->IsDiscrete()) {
      if (1 == this->mFactors[n]->GetNumberOfLevels()) {
        std::cerr << "ERROR: " << this->mFactors[n]->GetFactorName()
                  << " is invalid. Must have at least two levels." << std::endl;
        return (-1);
      }
    }
  }

  // now remove any continuous factors that have both zero mean and std
  int purged = this->PurgeNullFactors();
  if (purged) {
    std::cout << std::endl
              << "Purged " << purged << " null (zero mean and std) factors."
              << std::endl;
  }

  std::cout << std::endl
            << "Data table " << isFileName << " loaded." << std::endl;

  free(tmpstr);

  return 0;
}

/**
 * @return int
 * @param  isFileName
 * 
 */
int QdecDataTable::Save(const char *isFileName) {
  int nFactors = this->GetContinuousFactorNames().size() +
                 this->GetDiscreteFactorNames().size();
  int nInputs = this->mSubjects.size();

  if ((nFactors == 0) || (nInputs == 0)) {
    std::cerr << "The data table is empty!  Nothing to save!" << std::endl;
    return 1;
  }

  FILE *fp = fopen(isFileName, "w");
  if (NULL == fp) {
    std::cerr << "ERROR: Unable to open file" << isFileName
              << " for writing!\n";
    return 1;
  }

  const char *fsid = "fsid";
  fprintf(fp, "%s", fsid); // first line/column: our file identifier
  for (int m = 0; m < nInputs; m++) {
    std::vector<QdecFactor *> subjectFactors = this->mSubjects[m]->GetFactors();

    // handle first line: the factor names
    if (m == 0) {
      for (unsigned int n = 0; n < subjectFactors.size(); n++) {
        fprintf(fp, " %s", subjectFactors[n]->GetFactorName().c_str());
      }
      fprintf(fp, "\n");
    }

    // now the data belonging to each subject, one subject per line
    fprintf(fp, "%s ", this->mSubjects[m]->GetId().c_str());
    for (unsigned int n = 0; n < subjectFactors.size(); n++) {
      if (subjectFactors[n]->IsDiscrete()) {
        fprintf(fp, "%s ", subjectFactors[n]->GetDiscreteValue().c_str());
      } else if (subjectFactors[n]->Ignore()) {
        fprintf(fp, "%s ", subjectFactors[n]->GetIgnoreValue().c_str());
      } else if (subjectFactors[n]->IsOrdinal()) {
        double value = subjectFactors[n]->GetContinuousValue();
        if (isnan(value)) {
          fprintf(fp, "NaN ");
        } else {
          fprintf(fp, "%d ", (int)value);
        }
      } else {
        fprintf(fp, "%.9lf ", subjectFactors[n]->GetContinuousValue());
      }
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  std::cout << "Saved file " << isFileName << std::endl;

  return 0;
}

/**
 * dumps factors and inputs to filepointer (stdout, or file)
 * @param  iFilePointer
 */
void QdecDataTable::Dump(FILE *fp) {
  int nFactors = this->GetContinuousFactorNames().size() +
                 this->GetDiscreteFactorNames().size();
  int nInputs = this->mSubjects.size();

  if (nFactors == 0)
    return;
  if (nInputs == 0)
    return;

  fprintf(fp, "Input table: %s\n", this->GetFileName().c_str());

  fprintf(fp, "Subj#, SubjID, Data...\n");
  for (int m = 0; m < nInputs; m++) {
    fprintf(fp, "%5d %s ", m + 1, this->mSubjects[m]->GetId().c_str());
    std::vector<QdecFactor *> subjectFactors = this->mSubjects[m]->GetFactors();
    for (unsigned int n = 0; n < subjectFactors.size(); n++) {
      if (subjectFactors[n]->IsDiscrete())
        fprintf(fp, "%s ", subjectFactors[n]->GetDiscreteValue().c_str());
      else if (subjectFactors[n]->Ignore())
        fprintf(fp, "%s ", subjectFactors[n]->GetIgnoreValue().c_str());
      else
        fprintf(fp, "%lf ", subjectFactors[n]->GetContinuousValue());
    }
    fprintf(fp, "\n");
  }

  for (int n = 0; n < nFactors; n++) {
    fprintf(fp, "%d  %s  %s %d\n", n + 1,
            this->mFactors[n]->GetFactorName().c_str(),
            this->mFactors[n]->GetFactorTypeName().c_str(),
            (int)this->mFactors[n]->GetLevelNames().size());
    if (this->mFactors[n]->IsDiscrete()) {
      std::vector<std::string> levelNames = this->mFactors[n]->GetLevelNames();
      for (unsigned int l = 0; l < levelNames.size(); l++) {
        fprintf(fp, "  %3d  %s\n", l + 1, levelNames[l].c_str());
      }
    }
  }

  // for testing GetMeanAndStdDev()...
  std::vector<std::string> contFactorNames = this->GetContinuousFactorNames();
  fprintf(fp,
          "                Continuous Factors:         Mean:       StdDev:\n"
          "                -------------------         -----       -------\n");
  for (unsigned int i = 0; i < this->GetContinuousFactorNames().size(); i++) {
    std::vector<double> vals =
        this->GetMeanAndStdDev(contFactorNames[i].c_str());
    char val1[1000];
    char val2[1000];
    sprintf(val1, "%5.3f", vals[0]);
    sprintf(val2, "%5.3f", vals[1]);
    fprintf(fp, "%35s %13s %13s\n", contFactorNames[i].c_str(), val1, val2);
  }

  fprintf(fp,
          "\n"
          "Number of subjects:   %d\n"
          "Number of factors:    %d (%d discrete, %d continuous)\n"
          "Number of classes:    %d\n"
          "Number of regressors: %d\n",
          nInputs, nFactors, (int)this->GetDiscreteFactorNames().size(),
          (int)this->GetContinuousFactorNames().size(),
          this->GetNumberOfClasses(), this->GetNumberOfRegressors());

  fprintf(fp, "============================================================\n");

  fflush(fp);
}

/**
 * @return string
 */
std::string QdecDataTable::GetFileName() { return mfnFileName; }

/**
 * @return std::vector< string >
 */
std::vector<std::string> QdecDataTable::GetSubjectIDs() {
  std::vector<std::string> ids;

  for (unsigned int i = 0; i < this->mSubjects.size(); i++) {
    ids.push_back(this->mSubjects[i]->GetId());
  }

  return ids;
}

/**
 * @return std::vector< QdecSubject* >
 */
std::vector<QdecSubject *> QdecDataTable::GetSubjects() {
  return this->mSubjects;
}

/**
 * @return QdecFactor*
 * @param isFactorName
 */
QdecFactor *QdecDataTable::GetFactor(const char *isFactorName) {
  unsigned int i;
  for (i = 0; i < this->mFactors.size(); i++) {
    if (0 == strcmp(isFactorName, this->mFactors[i]->GetFactorName().c_str())) {
      return this->mFactors[i];
    }
  }
  if (i == this->mFactors.size()) {
    std::cerr << "ERROR: QdecDataTable::GetFactor: '" << isFactorName
              << "' is not in datatable!" << std::endl;
  }

  return NULL;
}

/**
 * @return QdecFactor*
 * @param isFactorName
 */
QdecFactor *QdecDataTable::GetFactor(const char *isSubjectName,
                                     const char *isFactorName) {
  QdecFactor *qf = NULL;
  for (unsigned int i = 0; i < this->mSubjects.size(); i++) {
    if (0 == strcmp(isSubjectName, this->mSubjects[i]->GetId().c_str())) {
      qf = this->mSubjects[i]->GetFactor(isFactorName);
      break;
    }
  }
  if (NULL == qf) {
    std::cerr << "ERROR: QdecDataTable::GetFactor: '" << isFactorName
              << "' is not in datatable!" << std::endl;
  }

  return qf;
}

/**
 * @return std::vector< string >
 */
std::vector<std::string> QdecDataTable::GetDiscreteFactorNames() {
  std::vector<std::string> discreteFactorNames;
  for (unsigned int i = 0; i < this->mFactors.size(); i++) {
    if (this->mFactors[i]->IsDiscrete()) {
      discreteFactorNames.push_back(this->mFactors[i]->GetFactorName());
    }
  }
  return discreteFactorNames;
}

/**
 * @return std::vector< string >
 */
std::vector<std::string> QdecDataTable::GetContinuousFactorNames() {
  std::vector<std::string> continuousFactorNames;
  for (unsigned int i = 0; i < this->mFactors.size(); i++) {
    if (this->mFactors[i]->IsContinuous()) {
      continuousFactorNames.push_back(this->mFactors[i]->GetFactorName());
    }
  }
  return continuousFactorNames;
}

/**
 * GetNumberOfClasses( ) - returns the number of subjects in the table
 */
int QdecDataTable::GetNumberOfSubjects() { return this->mSubjects.size(); }

/**
 * GetNumberOfClasses( ) - returns the number of classes for the design.
 * The number of classes is just all the combinations of all
 * the levels for the discrete factors.
 */
int QdecDataTable::GetNumberOfClasses() {
  int nClasses = 1;

  for (unsigned int i = 0; i < this->mFactors.size(); i++) {
    if (this->mFactors[i]->IsDiscrete()) {
      nClasses *= this->mFactors[i]->GetLevelNames().size();
    }
  }
  return (nClasses);
}

/**
 * GetNumberOfRegressors() - returns the number of regressors for the
 * given design.
 */
int QdecDataTable::GetNumberOfRegressors() {
  int nReg = this->GetNumberOfClasses() *
             (this->GetContinuousFactorNames().size() + 1);
  return (nReg);
}

/**
 * GetMeanAndStdDev() - computes the average and stddev of continuous factor
 * @return std::vector< double > - first element is mean, second is the stddev
 * @param isFactorName
 */
std::vector<double> QdecDataTable::GetMeanAndStdDev(const char *isFactorName) {
  // can't find a mean on a discrete factor:
  assert(this->GetFactor(isFactorName)->IsContinuous());
  // or if there aren't any subjects
  assert(this->GetSubjects().size());

  double                     d        = 0.0;
  double                     Sum      = 0.0;
  double                     Sum2     = 0.0;
  long                       N        = 0;
  std::vector<QdecSubject *> subjects = this->GetSubjects();
  for (unsigned int i = 0; i < this->GetSubjects().size(); i++, N++) {
    d = subjects[i]->GetContinuousFactorValue(isFactorName);
    if (isnan(d)) {
      N--;
    } else {
      Sum += d;
      Sum2 += (d * d);
    }
  }

  double Avg    = Sum / N;
  double StdDev = sqrt(N * (Sum2 / N - Avg * Avg) / (N - 1));

  std::vector<double> tmp;
  tmp.push_back(Avg);
  tmp.push_back(StdDev);

  return (tmp);
}

/**
 * deletes all continuous factors that have a zero mean and zero stddev.
 * those are useless factors (probably originating from a stats table).
 * @return number of factors purged
 */
int QdecDataTable::PurgeNullFactors() {
  int purgeCount = 0;

  std::vector<QdecFactor *>::iterator iter = mFactors.begin();
  while (iter != mFactors.end()) {
    QdecFactor *factor = *iter;
    if (factor->IsContinuous()) {
      std::string         factorName = factor->GetFactorName();
      std::vector<double> vals = this->GetMeanAndStdDev(factorName.c_str());
      if ((vals[0] == 0.0) && (vals[1] == 0.0)) {
        // this factor has both zero mean and stddev, so get rid of it
        iter = mFactors.erase(iter);
        purgeCount++;

        // now must also delete this factor from each subject
        std::vector<QdecSubject *>::iterator iterSubj = mSubjects.begin();
        while (iterSubj != mSubjects.end()) {
          QdecSubject *subject = *iterSubj;
          subject->DeleteFactor(factorName.c_str());
          ++iterSubj;
        }
      } else
        ++iter;
    } else
      ++iter;
  }

  return purgeCount;
}

/**
 * merge a factor from a given data table into this data table
 *
 * @return int
 * @param  isFactorName
 * @param  iDataTable
 */
int QdecDataTable::MergeFactor(const char *   isFactorName,
                               QdecDataTable *iDataTable) {
  // first, some sanity checks on this alien data
  if ((NULL == isFactorName) || (0 == strlen(isFactorName))) {
    std::cerr << "ERROR: QdecDataTable::MergeFactor: invalid factor name!"
              << std::endl;
    return 1;
  }
  if (iDataTable->GetNumberOfSubjects() != this->GetNumberOfSubjects()) {
    std::cerr << "ERROR: QdecDataTable::MergeFactor: invalid input data table: "
                 "mismatch in number of subjects"
              << std::endl;
    return 1;
  }
  std::vector<std::string> theirIDs = iDataTable->GetSubjectIDs();
  std::vector<std::string> ourIDs   = this->GetSubjectIDs();
  for (unsigned int i = 0; i < theirIDs.size(); i++) {
    if (ourIDs[i].compare(theirIDs[i])) {
      std::cerr
          << "ERROR: QdecDataTable::MergeFactor: invalid input data table: "
             "mismatch in name of subjects: "
          << theirIDs[i].c_str() << " vs. " << ourIDs[i].c_str() << std::endl;
      return 1;
    }
  }

  // lets check if this factor is already in our table
  for (unsigned int i = 0; i < this->mFactors.size(); i++) {
    const char *theFactor = this->mFactors[i]->GetFactorName().c_str();
    if (0 == strcmp(isFactorName, theFactor))
      return 0; // already got it
  }

  // so far so good, so add this factor to our list of factors
  QdecFactor *newFactor = new QdecFactor(iDataTable->GetFactor(isFactorName));
  this->mFactors.push_back(newFactor);

  // and add the data for each subject
  for (unsigned int i = 0; i < this->mSubjects.size(); i++) {
    // get factor data for this subject from input table
    QdecFactor *theirFactor = new QdecFactor(iDataTable->GetFactor(
        this->mSubjects[i]->GetId().c_str(), isFactorName));
    // and stick in our table
    this->mSubjects[i]->AddFactor(theirFactor);
  }

  return 0;
}

/**
 * delete a factor from this data table
 *
 * @return int
 * @param  isFactorName
 */
int QdecDataTable::DeleteFactor(const char *isFactorName) {
  // first, some sanity checks on this alien data
  if ((NULL == isFactorName) || (0 == strlen(isFactorName))) {
    std::cerr << "ERROR: QdecDataTable::DeleteFactor: invalid factor name"
              << std::endl;
    return 1;
  }

  // search and destroy!  search and destroy!
  std::vector<QdecFactor *>::iterator iter = mFactors.begin();
  while (iter != mFactors.end()) {
    QdecFactor *factor     = *iter;
    std::string factorName = factor->GetFactorName();
    if (0 == strcmp(isFactorName, factorName.c_str())) {
      // found it, so get rid of it
      iter = mFactors.erase(iter);

      // now must also delete this factor from each subject
      std::vector<QdecSubject *>::iterator iterSubj = mSubjects.begin();
      while (iterSubj != mSubjects.end()) {
        QdecSubject *subject = *iterSubj;
        subject->DeleteFactor(factorName.c_str());
        ++iterSubj;
      }
    } else
      ++iter;
  }

  return 0;
}
