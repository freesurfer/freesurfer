/*
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


#ifndef ConfigFile_h
#define ConfigFile_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include "mri_ca_util.h"


/* [[]]
         The CConfigFile class needs to be enhanced so that a configFile can be written out to disk without seeking in the file for section names
          This means that all the entries in a section must be grouped and passed to a write command so that that section can be written out all at once and
          not revisited . This will simplify the code and make it possible to write out to cout just like a file.

         This has been started in CStatisticsVolumeHeader.h and other config file containing classes should be updated so the obsolete write() functions can
           be removed from CConfigFile (grep on obsolete in this file)

         Class CConfigFile should be a namespace of static functions which enable more operations on the iostream object
         -if only methods on the iostream class are used then a filter program could be written with a configFile read in from cin or written to  cout

         -dont want to replicate fstrean methods in this class, just want to augment it
         -Want to write a simple filter program by simply reading from cin and writing to cout
         -DONT WANT CConfigFile to be derived from iostream because then the cout and cin objects cannot be written to since they are not of the inherited type.

         1) rewrite classes which use CConfigFile:
         cin >> labelLUT;
           class labelLUT {
              friend ostream& operator>>(ostream&, labelLUT&);
           }

         2) verify that the operations in ConfigFile can work on an iostream object by replacing fstream with iostream and removing the open and close operations
         3) declare the config file functions static.
         4) invoke the functions in the classes which use CConfigFile with CConfigFile::get(istream, ...) and CConfigFile::write(ostream, ...)


         Currently:
         A config file Extends fstream with special extraction and insertion operators which pertain to the cfg file format
          Which is basically an exteneded .ini file format from windows which is easy to read and supports
          singletons, vectors and vectors of vectors (including matrices) of ints, doubles and strings

*/

#define debugAid 0
#define EOL 10

using namespace std ;

class CConfigFile
{
private:

  // Assuming the last character retrieved indicates the beginning of a string we are interested in
  // occurs after the current character
  // get the next sequence of characters up to the first 'ch' character or EOL or EOF
  // If EOL is found first then bEOLFoundFirst is set to true
  string getDelimitedString(fstream& ifs, char ch, bool& bEOLFoundFirst)
  {

    bEOLFoundFirst=false;
    string strDelimitedString;  // initially an empty string
    strDelimitedString.erase();


    char c;
    ifs.get(c);
    while (ifs.good() && (c!=ch) && (c!=EOL))
    {
      strDelimitedString.append(1,c);
      ifs.get(c);
    };


    if (c==EOL) bEOLFoundFirst=true;
    return(strDelimitedString);

  }

  string getDelimitedString(fstream& ifs, char ch, bool& bEOLFoundFirst, bool& bEQUALFound, bool& bSqBracketFound, bool& bEOFFound)
  {
    bEOLFoundFirst=false;
    bEQUALFound=false;
    bSqBracketFound=false;
    bEOFFound=false;
    string strDelimitedString;  // initially an empty string
    strDelimitedString.erase();


    char c;
    ifs.get(c);
    if (debugAid) cout << c;
    while (ifs.good() && (c!=ch) && (c!=EOL) && (c!='=') && (c!='[') && (c!='\n') )
    {
      strDelimitedString.append(1,c);
      ifs.get(c);
      if (debugAid) cout << c;
    };

    if (!ifs.good())  bEOFFound=true;

    if (c==EOL) bEOLFoundFirst=true;
    if (c=='=') bEQUALFound=true;
    if (c=='[') bSqBracketFound=true;

    return(strDelimitedString);

  }


  bool skipOver(fstream& ifs, char ch)
  {
    char c;

    ifs.get(c);
    while (ifs.good() && (c==ch))
    {
      ifs.get(c);
    };


    if (ifs.good())
    {
      ifs.putback(c);
      return true;
    }
    else
    {
      return false;
    }
  }

  bool skipToNoThrow(fstream& ifs, char ch)
  {
    char c;

    ifs.get(c);
    if (debugAid) cout << c;
    while (ifs.good() && (c!=ch))
    {
      ifs.get(c);
      if (debugAid) cout << c;
    };


    if (ifs.good())
    {
      return true;
    }
    else
    {
      return false;
    }
  }


  // skip to a '[' and then
  // get the sequence of characters up to the first ']'
  // If either cannot be found then the section name is left unchanged..
  bool getNextSectionName(fstream& ifs, string& strSectionName)
  {
    bool bFoundSection=false;
    bool bEOLFound;
    if  (skipToNoThrow(ifs,'['))
    {
      strSectionName=getDelimitedString(ifs, ']',bEOLFound);
      if (ifs.good())
        bFoundSection=true;
    }

    return(bFoundSection);
  }

  // skip to the end of the current line next non blank and get the sequence of characters up to the first '='
  string getNextValueName(fstream& ifs)
  {
    bool bEOLFound=true;
    string strNextValueName;

    while ((ifs.good()) && (bEOLFound))   // keep looking until we find a value which is terminated by an '=' char and not by EOL
      // because that could simply by part of a matrix definition and not a valid value in a config file
    {
      skipOver(ifs, ' ');
      strNextValueName=getDelimitedString(ifs,'=',bEOLFound);
    };

    return(strNextValueName);
  }

  bool getNextVectorElement(fstream& ifs, string& strNextElement, bool& bMoreCols, bool& bMoreRows)
  {
    if (!ifs.good()) return(false);

    bool bEQUALFound;
    bool bSqBracketFound;
    bool bEOLFound;
    bool bEOFFound;

    bMoreCols=true;
    bMoreRows=true;

    skipOver(ifs, ' ');
    strNextElement=getDelimitedString(ifs,',',bEOLFound, bEQUALFound, bSqBracketFound, bEOFFound);

    bMoreRows= (! ( (bEQUALFound) || (bSqBracketFound)  || (bEOFFound)  ));
    bMoreCols=((bMoreRows) && (!bEOLFound));

    if ( (bEQUALFound) || (bSqBracketFound) )
      return(false);
    else
    {
      if ((strNextElement.size()==0) && (bMoreCols==false)) // if this is the only element on a row and its empty then discount it
        return(false);
      else
        return(true);
    }

  }


  bool skipToSection(fstream& ifs, string& strSectionName)
  {
    bool bFoundSection=false;
    string strNextSectionName;
    while ((ifs.good()) && (getNextSectionName(ifs, strNextSectionName)))
    {
      if (strNextSectionName==strSectionName)
      {
        bFoundSection=true;
        break;
      }
    }

    if (!fs.good())
    {
      fs.clear();  // we read off the end of the file so clear the fs error bits
      fs.seekg(0,ios::end);  // move the file pointer to the end of the file

    }

    return(bFoundSection);
  }

  // Skip to the equal sign after the value strValueName but dont go beyond the current section
  bool skipToValue(fstream& ifs, string& strValueName)
  {
    bool bFoundValue=false;
    if (strValueName=="")
    {
      bFoundValue=true;
    }
    else
    {
      while (ifs.good())
      {
        string strNextValueName=getNextValueName(ifs);
        if (strNextValueName==strValueName)
        {
          bFoundValue=true;
          break;
        }
      }

      if (!fs.good())
      {
        fs.clear();  // we read off the end of the file so clear the fs error bits
        fs.seekg(0,ios::end);  // move the file pointer to the end of the file

      }



    }
    return(bFoundValue);
  }

  bool skipTo(string& strSectionName, string& strValueName)
  {
    bool bFoundSectionAndValue=false;
    if (skipToSection(fs,  strSectionName))
    {
      if (skipToNoThrow(fs,'\n'))
      {
        if (skipToValue(fs, strValueName))
          bFoundSectionAndValue=true;
      }
    }

    if (!fs.good())
    {
      fs.clear();  // we read off the end of the file so clear the fs error bits
      fs.seekg(0,ios::end);  // move the file pointer to the end of the file

    }

    return( bFoundSectionAndValue);
  }




public:
  // default ctor, copy, assign OK
  ~CConfigFile()
  {
    fs.close();

  }



  // Opens and reads in the file strFilepath, throws if there is an error
  // DEFAULT file mode is  CREATE (create the file if nec.), OVERWRITE_AT_BEGINNING (for all new data written to the file)
  // over-ride these default by passing in ios::nocreate  or ios::trunc   or ios::append as the nProt argument
  void init(string strTheFilepath, int nProt = 0)
  {
    strFilepath=strTheFilepath;

    fs.open(strFilepath.c_str(), ios::in | ios::out | nProt);


    if (!fs)
    {
      char* strErrorMsg=new char[500];
      sprintf(strErrorMsg,"Error: could not open config file %s", strFilepath.c_str());
      throw(strErrorMsg);
    }

  }

  // methods to retrieve values from a config file
  bool get(string& strValue, string strSectionName, string valueName)   // get string
  {
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    if (skipTo(strSectionName, valueName))
    {
      bool bEOLFound;
      strValue=getDelimitedString(fs, ' ',bEOLFound);
    }
  }


  bool get(TypeVectorString& vectStrValues, string strSectionName, string valueName)   // get 1D vector of strings
  {
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    if (skipTo(strSectionName, valueName))
    {
      // read in the comma separated strings into  vectIntValues
      //
      string strNextElement;
      bool bMoreCols=true;
      bool bMoreRowsDummy;
      while ((bMoreCols==true) && (getNextVectorElement(fs, strNextElement, bMoreCols, bMoreRowsDummy)) )
      {
        vectStrValues.push_back(strNextElement);
      }
    }
  }

  bool get(TypeMatrixString& matrixStrValues, string strSectionName, string valueName)   // get 2D matrix of strings
  {
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    if (skipTo(strSectionName, valueName))
    {

      // read in the comma separated strings into  matrixStrValues
      // until there is a line which contains a [, = or EOF
      // permit the number of columns in ea row to be different
      //
      string strNextElement;
      bool bMoreRows=true;
      while (bMoreRows==true)
      {
        TypeVectorString vectString;
        bool bMoreCols=true;
        while ( (bMoreCols==true) && (getNextVectorElement(fs, strNextElement, bMoreCols, bMoreRows)) )
        {
          vectString.push_back(strNextElement);
        }
        if (vectString.size()!=0) matrixStrValues.push_back(vectString);
      }
    }

  }



  bool get(int& nValue, string strSectionName, string valueName)   // get single int
  {
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    if (skipTo(strSectionName, valueName))
    {
      bool bEOLFound;
      string strValue=getDelimitedString(fs, ' ',bEOLFound);
      nValue=atoi(strValue.c_str());
    }
  }


  bool get(TypeVectorInt& vectIntValues, string strSectionName, string valueName)   // get vector of ints
  {
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    if (skipTo(strSectionName, valueName))
    {

      // read in the comma separated strings into  vectIntValues
      //
      string strNextElement;
      bool bMoreCols=true;
      bool bMoreRowsDummy;
      while ((bMoreCols==true) && (getNextVectorElement(fs, strNextElement, bMoreCols, bMoreRowsDummy)) )
      {
        if (strNextElement.size()>0)
        {
          int i=atoi(strNextElement.c_str());
          vectIntValues.push_back(i);
        }
      }
    }
  }



  bool get(TypeMatrixInt& matrixIntValues, string strSectionName, string valueName)   // get 2D matrix of ints
  {
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    if (skipTo(strSectionName, valueName))
    {

      // read in the comma separated strings into  matrixStrValues
      // until there is a line which contains a [, = or EOF
      // permit the number of columns in ea row to be different
      //
      string strNextElement;
      bool bMoreRows=true;
      while (bMoreRows==true)
      {
        TypeVectorInt vectInt;
        bool bMoreCols=true;
        while ( (bMoreCols==true) && (getNextVectorElement(fs, strNextElement, bMoreCols, bMoreRows)) )
        {
          if (strNextElement.size()>0)
          {
            int i=atoi(strNextElement.c_str());
            vectInt.push_back(i);
          }
        }
        if (vectInt.size()!=0)  matrixIntValues.push_back(vectInt);
      }
    }

  }

  bool get(double& dValue, string strSectionName, string valueName)   // get single int
  {
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    if (skipTo(strSectionName, valueName))
    {
      bool bEOLFound;
      string strValue=getDelimitedString(fs, ' ',bEOLFound);
      dValue=atof(strValue.c_str());
    }
  }

  bool get(float& fValue, string strSectionName, string valueName)
  {
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    if (skipTo(strSectionName, valueName))
    {
      bool bEOLFound;
      string strValue=getDelimitedString(fs, ' ',bEOLFound);
      fValue=atof(strValue.c_str());
    }
  }

  bool get(TypeVectorDouble& vectDouble, string strSectionName, string valueName)
  {
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    if (skipTo(strSectionName, valueName))
    {

      // read in the comma separated strings into  vectIntValues
      //
      string strNextElement;
      bool bMoreCols=true;
      bool bMoreRowsDummy;
      while ((bMoreCols==true) && (getNextVectorElement(fs, strNextElement, bMoreCols, bMoreRowsDummy)) )
      {
        if (strNextElement.size()>0)
        {
          double d=atof(strNextElement.c_str());
          vectDouble.push_back(d);
        }
      }
    }
  }


  bool get(TypeMatrixDouble& matrixDoubleValues, string strSectionName, string valueName)
  {
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    if (skipTo(strSectionName, valueName))
    {

      // read in the comma separated strings into  matrixStrValues
      // until there is a line which contains a [, = or EOF
      // permit the number of columns in ea row to be different
      //
      string strNextElement;
      bool bMoreRows=true;
      while (bMoreRows==true)
      {
        TypeVectorDouble vectDouble;
        bool bMoreCols=true;
        while ( (bMoreCols==true) && (getNextVectorElement(fs, strNextElement, bMoreCols, bMoreRows)) )
        {
          if (strNextElement.size()>0)
          {
            double d=atof(strNextElement.c_str());
            vectDouble.push_back(d);
          }
        }
        if (vectDouble.size()!=0)  matrixDoubleValues.push_back(vectDouble);
      }
    }
  }

  bool get(TypeVectorFloat& vectFloat, string strSectionName, string valueName)
  {
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    if (skipTo(strSectionName, valueName))
    {

      // read in the comma separated strings into  vectIntValues
      //
      string strNextElement;
      bool bMoreCols=true;
      bool bMoreRowsDummy;
      while ((bMoreCols==true) && (getNextVectorElement(fs, strNextElement, bMoreCols, bMoreRowsDummy)) )
      {
        if (strNextElement.size()>0)
        {
          float f=atof(strNextElement.c_str());
          vectFloat.push_back(f);
        }
      }
    }
  }


  bool get(TypeMatrixFloat& matrixFloatValues, string strSectionName, string valueName)
  {
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    if (skipTo(strSectionName, valueName))
    {

      // read in the comma separated strings into  matrixStrValues
      // until there is a line which contains a [, = or EOF
      // permit the number of columns in ea row to be different
      //
      string strNextElement;
      bool bMoreRows=true;
      while (bMoreRows==true)
      {
        TypeVectorFloat vectFloat;
        bool bMoreCols=true;
        while ( (bMoreCols==true) && (getNextVectorElement(fs, strNextElement, bMoreCols, bMoreRows)) )
        {
          if (strNextElement.size()>0)
          {
            float f=atof(strNextElement.c_str());
            vectFloat.push_back(f);
          }
        }
        if (vectFloat.size()!=0)  matrixFloatValues.push_back(vectFloat);
      }
    }
  }

  void writeSection(string strSectionName,bool bEchoToStdOut=false )
  {
    if (bEchoToStdOut)
    {
      cout << "  [" << strSectionName << "]\n";
    }
    else
    {
      fs << "  [" << strSectionName << "]\n";
    }
  }

  // Values can only be written out to a file once

  void writeValueName(string strValueName,bool bEchoToStdOut=false )
  {
    if (strValueName!="")
    {
      if (bEchoToStdOut)
      {
        cout << "    " << strValueName << "=";
      }
      else
      {
        fs << "    " << strValueName << "=";
      }
    }
  }
  // methods to be used by specific config file implementations to write out sections and values in a file

// Strings -----------
  void write(string strValue, string strValueName,bool bEchoToStdOut=false)
  {
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut)
    {
      cout << strValue << "\n";
    }
    else
    {
      fs << strValue << "\n";
    }
  }

  void write(TypeVectorString vectStrValues, string strValueName,bool bEchoToStdOut=false)
  {
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(vectStrValues);
    else printType(vectStrValues, fs);
  }

  void write(TypeMatrixString matrixStrValues, string strValueName,bool bEchoToStdOut=false)
  {
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(matrixStrValues);
    else printType(matrixStrValues, fs);
  }

  // Ints -----------
  void write(int nValue, string strValueName,bool bEchoToStdOut=false)
  {
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut)
    {
      cout << nValue << "\n";
    }
    else
    {
      fs <<  nValue << "\n";
    }
  }

  void write(TypeVectorInt vectIntValues, string strValueName,bool bEchoToStdOut=false)
  {
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(vectIntValues);
    else printType(vectIntValues, fs);
  }

  void write(TypeMatrixInt matrixIntValues, string strValueName,bool bEchoToStdOut=false)
  {
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(matrixIntValues);
    else printType(matrixIntValues, fs);
  }

// floats -----------
  void write(float fValue, string strValueName,bool bEchoToStdOut=false)
  {
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut)
    {
      cout << fValue << "\n";
    }
    else
    {
      fs << fValue << "\n";
    }
  }

  void write(TypeVectorFloat vectFloatValues, string strValueName,bool bEchoToStdOut=false)
  {
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(vectFloatValues);
    else printType(vectFloatValues, fs);
  }

  void write(TypeMatrixFloat matrixFloatValues, string strValueName,bool bEchoToStdOut=false)
  {
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(matrixFloatValues);
    else printType(matrixFloatValues, fs);
  }
// Doubles -----------
  void write(double dValue, string strValueName,bool bEchoToStdOut=false)
  {
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut)
    {
      cout << dValue << "\n";
    }
    else
    {
      fs << dValue << "\n";
    }
  }

  void write(TypeVectorDouble vectDoubleValues, string strValueName ,bool bEchoToStdOut=false)
  {
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(vectDoubleValues);
    else printType(vectDoubleValues, fs);
  }

  void write(TypeMatrixDouble matrixDoubleValues, string strValueName ,bool bEchoToStdOut=false)
  {
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(matrixDoubleValues);
    else printType(matrixDoubleValues, fs);

  }

#if 0
  // ============== obsolete ===================

// Strings -----------
  void write(string strValue, string strSectionName, string strValueName,bool bEchoToStdOut=false)
  {
    // if strSectionName is in the file move to it, ow write it to the file
    if ((bEchoToStdOut) || (!skipTo(strSectionName, strValueName)) )
    {
      writeSection(strSectionName,bEchoToStdOut);
    }

    // write the value
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut)
    {
      cout << strValue << "\n";
    }
    else
    {
      fs << strValue << "\n";
    }


  }

  void write(TypeVectorString vectStrValues, string strSectionName, string strValueName,bool bEchoToStdOut=false)
  {
    // if strSectionName is in the file move to it, ow write it to the file
    if ((bEchoToStdOut) || (!skipTo(strSectionName, strValueName)) )
    {
      writeSection(strSectionName,bEchoToStdOut);
    }

    // write the value
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(vectStrValues);
    else printType(vectStrValues, fs);
  }

  void write(TypeMatrixString matrixStrValues, string strSectionName, string strValueName,bool bEchoToStdOut=false)
  {
    // if strSectionName is in the file move to it, ow write it to the file
    if ((bEchoToStdOut) || (!skipTo(strSectionName, strValueName)) )
    {
      writeSection(strSectionName,bEchoToStdOut);
    }

    // write the value
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(matrixStrValues);
    else printType(matrixStrValues, fs);
  }

  // Ints -----------
  void write(int nValue, string strSectionName, string strValueName,bool bEchoToStdOut=false)
  {
    if (bEchoToStdOut)
    {
      writeSection(strSectionName,bEchoToStdOut);
      writeValueName(strValueName, bEchoToStdOut);
      cout << nValue << "\n";
    }
    else
    {
      fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
      // if strSectionName is in the file move to it, ow write it to the file
      if (!skipToSection(fs,strSectionName))
      {
        writeSection(strSectionName,bEchoToStdOut);
      }

      // write the value
      if (!skipToValue(fs,strValueName))
      {
        fs.seekg(0,ios::beg);
        skipToSection(fs,strSectionName);
        writeValueName(strValueName, bEchoToStdOut);
      }

      fs << nValue << "\n";

    }




  }

  void write(TypeVectorInt vectIntValues, string strSectionName, string strValueName,bool bEchoToStdOut=false)
  {
    // if strSectionName is in the file move to it, ow write it to the file
    if ((bEchoToStdOut) || (!skipTo(strSectionName, strValueName)) )
    {

      writeSection(strSectionName,bEchoToStdOut);
    }

    // write the value
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(vectIntValues);
    else printType(vectIntValues, fs);
  }

  void write(TypeMatrixInt matrixIntValues, string strSectionName, string strValueName,bool bEchoToStdOut=false)
  {
    // if strSectionName is in the file move to it, ow write it to the file
    if ((bEchoToStdOut) || (!skipTo(strSectionName, strValueName)) )
    {
      writeSection(strSectionName,bEchoToStdOut);
    }

    // write the value
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(matrixIntValues);
    else printType(matrixIntValues, fs);
  }

// floats -----------
  void write(float fValue, string strSectionName, string strValueName,bool bEchoToStdOut=false)
  {
    // if strSectionName is in the file move to it, ow write it to the file
    if ((bEchoToStdOut) || (!skipTo(strSectionName, strValueName)) )
    {
      writeSection(strSectionName,bEchoToStdOut);
    }

    // write the value
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut)
    {
      cout << fValue << "\n";
    }
    else
    {
      fs << fValue << "\n";
    }
  }

// Doubles -----------
  void write(double dValue, string strSectionName, string strValueName,bool bEchoToStdOut=false)
  {
    // if strSectionName is in the file move to it, ow write it to the file
    if ((bEchoToStdOut) || (!skipTo(strSectionName, strValueName)) )
    {
      writeSection(strSectionName,bEchoToStdOut);
    }

    // write the value
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut)
    {
      cout << dValue << "\n";
    }
    else
    {
      fs << dValue << "\n";
    }
  }

  void write(TypeVectorDouble vectDoubleValues, string strSectionName, string strValueName ,bool bEchoToStdOut=false)
  {
    // if strSectionName is in the file move to it, ow write it to the file
    if ((bEchoToStdOut) || (!skipTo(strSectionName, strValueName)) )
    {
      writeSection(strSectionName,bEchoToStdOut);
    }

    // write the value
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(vectDoubleValues);
    else printType(vectDoubleValues, fs);
  }

  void write(TypeMatrixDouble matrixDoubleValues, string strSectionName, string strValueName ,bool bEchoToStdOut=false)
  {
    // if strSectionName is in the file move to it, ow write it to the file
    if ((bEchoToStdOut) || (!skipTo(strSectionName, strValueName)) )
    {
      writeSection(strSectionName,bEchoToStdOut);
    }

    // write the value
    writeValueName(strValueName, bEchoToStdOut);
    if (bEchoToStdOut) printType(matrixDoubleValues);
    else printType(matrixDoubleValues, fs);

  }

  // =============end obsolete =================
#endif

  void print()
  {
    char strMsg[500];
    sprintf(strMsg,"ConfigFilename: %s",strFilepath.c_str());
    cout << strMsg;
  }

  // retrieve a list of [section names] containing a given prefix
  void find(TypeVectorString& vectString, string strPrefix)
  {
    vectString.clear();
    fs.clear();
    fs.seekg(0,ios::beg);  // move the file pointer to the beginning of the file
    string strNextSectionName;
    while ((fs.good()) && (getNextSectionName(fs,strNextSectionName)))
    {
      if (!strncmp(strNextSectionName.c_str(),strPrefix.c_str(), strlen(strPrefix.c_str()))  )
      {
        vectString.push_back(strNextSectionName);
      }
    }

  }


  // test the methods in CConfigFile
  // returns true if the test passed and false otherwise
  bool test()
  {
    bool bPassed=true;
    // [[ these tests should be self verifying tests !! ]]
    // int tests
    int singleInt;
    int singleIntActual=10;
    get(singleInt, "test", "singleInt");
    cout << "singleInt = " << singleInt << "\n";
    if (singleInt!=singleIntActual) bPassed=false;

    TypeVectorInt vectorInt;
    TypeVectorInt vectorIntActual;
    vectorIntActual.push_back(-2);
    vectorIntActual.push_back(3);
    vectorIntActual.push_back(4);
    vectorIntActual.push_back(5);
    vectorIntActual.push_back(-345);
    get(vectorInt,"test","vectorInt");
    cout << "vectorInt = ";
    printType(vectorInt);
    cout << "\n";
    if (vectorInt!=vectorIntActual) bPassed=false;

    TypeMatrixInt matrixInt,matrixInt2;
    TypeMatrixInt matrixIntActual;
    TypeVectorInt vectIntRow;
    vectIntRow.erase(vectIntRow.begin(),vectIntRow.end());
    vectIntRow.push_back(1);
    vectIntRow.push_back(2);
    vectIntRow.push_back(-3);
    matrixIntActual.push_back(vectIntRow);
    vectIntRow.erase(vectIntRow.begin(),vectIntRow.end());
    vectIntRow.push_back(-3);
    vectIntRow.push_back(0);
    vectIntRow.push_back(5);
    matrixIntActual.push_back(vectIntRow);
    get( matrixInt, "test","matrixInt");
    cout << "matrixInt = " ;
    printType(matrixInt);
    if (matrixInt!=matrixIntActual) bPassed=false;
    cout << "\n";

    TypeMatrixInt matrixInt2Actual;
    vectIntRow.erase(vectIntRow.begin(),vectIntRow.end());
    vectIntRow.push_back(2);
    vectIntRow.push_back(3);
    vectIntRow.push_back(4);
    vectIntRow.push_back(5);
    matrixInt2Actual.push_back(vectIntRow);
    vectIntRow.erase(vectIntRow.begin(),vectIntRow.end());
    vectIntRow.push_back(4);
    vectIntRow.push_back(5);
    vectIntRow.push_back(6);
    matrixInt2Actual.push_back(vectIntRow);
    vectIntRow.erase(vectIntRow.begin(),vectIntRow.end());
    vectIntRow.push_back(7);
    vectIntRow.push_back(8);
    vectIntRow.push_back(9);
    vectIntRow.push_back(0);
    matrixInt2Actual.push_back(vectIntRow);
    get( matrixInt2, "test","matrixInt2");
    cout << "matrixInt2 = " ;
    printType(matrixInt2);
    cout << "\n";
    if (matrixInt2!=matrixInt2Actual) bPassed=false;

    // string tests
    TypeVectorString vectString;
    TypeVectorString vectStringActual;
    vectStringActual.erase(vectStringActual.begin(),vectStringActual.end());
    vectStringActual.push_back("yata");
    vectStringActual.push_back("gobba");
    get(vectString,"test", "vectString");
    cout << "vectString = ";
    printType(vectString);
    cout << "\n";
    if (vectString!=vectStringActual) bPassed=false;

    double singleDouble;
    double singleDoubleActual=345.63;
    get(singleDouble, "test", "singleDouble");
    cout << "singleDouble = " << singleDouble << "\n";
    if (singleDouble!=singleDoubleActual) bPassed=false;


    TypeMatrixDouble matrixDouble;
    TypeMatrixDouble matrixDoubleActual;
    TypeVectorDouble vectDoubleRow;
    vectDoubleRow.erase(vectDoubleRow.begin(),vectDoubleRow.end());
    vectDoubleRow.push_back(45.0);
    vectDoubleRow.push_back(56);
    vectDoubleRow.push_back(23.232);
    matrixDoubleActual.push_back(vectDoubleRow);
    vectDoubleRow.erase(vectDoubleRow.begin(),vectDoubleRow.end());
    vectDoubleRow.push_back(12.0);
    vectDoubleRow.push_back(21);
    vectDoubleRow.push_back(-23.0);
    matrixDoubleActual.push_back(vectDoubleRow);
    get(matrixDouble,"test","matrixDouble");
    cout << "matrixDouble = ";
    printType(matrixDouble);
    cout << "\n";
    if (matrixDouble!=matrixDoubleActual) bPassed=false;


    /* Here's a section that can be placed into a config file
    [test]
    singleInt=10
    singleDouble=345.63
    vectorInt=-2, 3, 4, 5, -345
    matrixInt= 1, 2, -3
    -3, 0,  5
    matrixInt2=
    2, 3, 4, 5
    4, 5, 6
    7, 8, 9, 0
    vectString=yata, gobba
    matrixDouble= 45.0,  56,  23.232
     12.0,  21,  -23.0
     */

    if (bPassed)
    {
      cout << "Passed all tests.\n";
    }
    else
    {
      cout << "One or more tests were failed.\n";
    }
    return bPassed;
  }



protected:
  string strFilepath;
  fstream fs;  // an i/o stream  which writes to files




};

#endif



