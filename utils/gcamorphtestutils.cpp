/**
 * @brief Utilities to help with testing GCAmorph routines
 *
 * Reference:
 * "How to Stay Sane while Programming: Collected Wisdom from Broadmoor"
 */
/*
 * Original Author: Richard Edgar
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

#include <iostream>
#include <vector>
using namespace std;

#include "gcamorphtestutils.hpp"
#include "gcamorphtestutils.h"

// ======================================================================

GCAMorphUtils::GCAMorphUtils(void) : varTypeMap(), scalarTypeMap()
{
  std::cout << __FUNCTION__ << ": Will not write out saved, saved2 and saved original positions" << std::endl;
  std::cout << __FUNCTION__ << ": Will not write out spacing scalar value" << std::endl;

  // Sanity check
  if (this->nVars != 24) {
    cerr << __FUNCTION__ << ": Invalid nVars!" << endl;
    exit(EXIT_FAILURE);
  }
  if (this->nDims != 3) {
    cerr << __FUNCTION__ << ": Invalid nDims!" << endl;
    exit(EXIT_FAILURE);
  }
  if (this->nScalars != 2) {
    cerr << __FUNCTION__ << ": Invalid nScalars!" << endl;
  }

  // Create the variable type map
  this->varTypeMap["rx"] = NC_DOUBLE;
  this->varTypeMap["ry"] = NC_DOUBLE;
  this->varTypeMap["rz"] = NC_DOUBLE;

  this->varTypeMap["origx"] = NC_DOUBLE;
  this->varTypeMap["origy"] = NC_DOUBLE;
  this->varTypeMap["origz"] = NC_DOUBLE;

  this->varTypeMap["dx"] = NC_FLOAT;
  this->varTypeMap["dy"] = NC_FLOAT;
  this->varTypeMap["dz"] = NC_FLOAT;

  this->varTypeMap["odx"] = NC_FLOAT;
  this->varTypeMap["ody"] = NC_FLOAT;
  this->varTypeMap["odz"] = NC_FLOAT;

  this->varTypeMap["origArea"] = NC_FLOAT;
  this->varTypeMap["origArea1"] = NC_FLOAT;
  this->varTypeMap["origArea2"] = NC_FLOAT;
  this->varTypeMap["area"] = NC_FLOAT;
  this->varTypeMap["area1"] = NC_FLOAT;
  this->varTypeMap["area2"] = NC_FLOAT;

  this->varTypeMap["invalid"] = NC_CHAR;
  this->varTypeMap["label"] = NC_INT;
  this->varTypeMap["status"] = NC_INT;
  this->varTypeMap["labelDist"] = NC_FLOAT;

  this->varTypeMap["mean"] = NC_FLOAT;
  this->varTypeMap["variance"] = NC_FLOAT;

  // And another sanity check
  if (this->varTypeMap.size() != this->nVars) {
    cerr << __FUNCTION__ << ": Incorrect entries in varTypeMap" << endl;
    exit(EXIT_FAILURE);
  }

  // Create the scalar type map
  this->scalarTypeMap["exp_k"] = NC_DOUBLE;
  this->scalarTypeMap["neg"] = NC_INT;

  // Double check the size
  if (this->scalarTypeMap.size() != this->nScalars) {
    cerr << __FUNCTION__ << ": Incorrect entries in scalarTypeMap" << endl;
    exit(EXIT_FAILURE);
  }
}

void GCAMorphUtils::Write(const GCAM *src, string fName) const
{
  // Construct the filename, using fName passed by value
  fName += ".nc";

  std::cout << __FUNCTION__ << ": Writing file " << fName << " ... ";
  std::cout.flush();

  // Sanity check the GCAM
  if (src->ninputs != 1) {
    std::cerr << __FUNCTION__ << ": Must have ninputs=1" << std::endl;
  }

  // Reference for the file
  int ncid;

  // Open the file
  NC_SAFE_CALL(nc_create(fName.c_str(), NC_CLOBBER, &ncid));

  // Set up the dimensions
  int dimIDs[this->nDims];
  NC_SAFE_CALL(nc_def_dim(ncid, "x", src->width, &dimIDs[this->iX]));
  NC_SAFE_CALL(nc_def_dim(ncid, "y", src->height, &dimIDs[this->iY]));
  NC_SAFE_CALL(nc_def_dim(ncid, "z", src->depth, &dimIDs[this->iZ]));

  // Set up the variable IDs
  map< string, int > varIDmap;
  map< string, nc_type >::const_iterator myIt;

  for (myIt = this->varTypeMap.begin(); myIt != this->varTypeMap.end(); myIt++) {
    NC_SAFE_CALL(nc_def_var(ncid,
                            myIt->first.c_str(),  // Name of the variable
                            myIt->second,         // Type of the variable
                            this->nDims,
                            dimIDs,
                            &varIDmap[myIt->first]));
  }

  // Sanity check
  if (varTypeMap.size() != varIDmap.size()) {
    cerr << __FUNCTION__ << ": Failed to create varIDmap correctly" << endl;
    exit(EXIT_FAILURE);
  }

  // Set up the scalars

  map< string, int > scalarIDmap;
  map< string, nc_type >::const_iterator scalarVarIt;

  // Do the doubles
  for (scalarVarIt = this->scalarTypeMap.begin(); scalarVarIt != this->scalarTypeMap.end(); scalarVarIt++) {
    NC_SAFE_CALL(nc_def_var(ncid,
                            scalarVarIt->first.c_str(),  // Name of the variable
                            scalarVarIt->second,         // Type of variable
                            0,                           // Zero dimensions implies scalar
                            NULL,                        // No pointer to dimensions for scalar
                            &scalarIDmap[scalarVarIt->first]));
  }

  // Make the end of the 'definition' region
  NC_SAFE_CALL(nc_enddef(ncid));

  // Pack into contiguous arrays
  const size_t nElems = src->width * src->height * src->depth;

  vector< double > x(nElems), y(nElems), z(nElems);
  vector< double > origx(nElems), origy(nElems), origz(nElems);
  vector< float > dx(nElems), dy(nElems), dz(nElems);
  vector< float > odx(nElems), ody(nElems), odz(nElems);
  vector< float > area(nElems), area1(nElems), area2(nElems);
  vector< float > origArea(nElems), origArea1(nElems), origArea2(nElems);
  vector< char > invalid(nElems);
  vector< int > label(nElems), status(nElems);
  vector< float > labelDist(nElems);
  vector< float > mean(nElems), variance(nElems);

  // Ugly loop to do the writing element by element
  for (int i = 0; i < src->width; i++) {
    for (int j = 0; j < src->height; j++) {
      for (int k = 0; k < src->depth; k++) {
        const GCA_MORPH_NODE &gcamn = src->nodes[i][j][k];

        const size_t i1d = k + (src->depth * (j + (src->height * i)));

        x.at(i1d) = gcamn.x;
        y.at(i1d) = gcamn.y;
        z.at(i1d) = gcamn.z;

        origx.at(i1d) = gcamn.origx;
        origy.at(i1d) = gcamn.origy;
        origz.at(i1d) = gcamn.origz;

        dx.at(i1d) = gcamn.dx;
        dy.at(i1d) = gcamn.dy;
        dz.at(i1d) = gcamn.dz;

        odx.at(i1d) = gcamn.odx;
        ody.at(i1d) = gcamn.ody;
        odz.at(i1d) = gcamn.odz;

        area.at(i1d) = gcamn.area;
        area1.at(i1d) = gcamn.area1;
        area2.at(i1d) = gcamn.area2;

        origArea.at(i1d) = gcamn.orig_area;
        origArea1.at(i1d) = gcamn.orig_area1;
        origArea2.at(i1d) = gcamn.orig_area2;

        invalid.at(i1d) = gcamn.invalid;
        label.at(i1d) = gcamn.label;
        status.at(i1d) = gcamn.status;
        labelDist.at(i1d) = gcamn.label_dist;

        if (gcamn.gc != NULL) {
          mean.at(i1d) = gcamn.gc->means[0];
          variance.at(i1d) = gcamn.gc->covars[0];
        }
        else {
          mean.at(i1d) = -1;
          variance.at(i1d) = -1;
        }
      }
    }
  }

  // We use 'find' to get an exception if the name doesn't exist

  NC_SAFE_CALL(nc_put_var_double(ncid, varIDmap.find("rx")->second, &x[0]));
  NC_SAFE_CALL(nc_put_var_double(ncid, varIDmap.find("ry")->second, &y[0]));
  NC_SAFE_CALL(nc_put_var_double(ncid, varIDmap.find("rz")->second, &z[0]));

  NC_SAFE_CALL(nc_put_var_double(ncid, varIDmap.find("origx")->second, &origx[0]));
  NC_SAFE_CALL(nc_put_var_double(ncid, varIDmap.find("origy")->second, &origy[0]));
  NC_SAFE_CALL(nc_put_var_double(ncid, varIDmap.find("origz")->second, &origz[0]));

  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("dx")->second, &dx[0]));
  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("dy")->second, &dy[0]));
  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("dz")->second, &dz[0]));

  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("odx")->second, &odx[0]));
  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("ody")->second, &ody[0]));
  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("odz")->second, &odz[0]));

  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("origArea")->second, &origArea[0]));
  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("origArea1")->second, &origArea1[0]));
  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("origArea2")->second, &origArea2[0]));

  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("area")->second, &area[0]));
  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("area1")->second, &area1[0]));
  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("area2")->second, &area2[0]));

  NC_SAFE_CALL(nc_put_var_text(ncid, varIDmap.find("invalid")->second, &invalid[0]));

  NC_SAFE_CALL(nc_put_var_int(ncid, varIDmap.find("label")->second, &label[0]));
  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("labelDist")->second, &labelDist[0]));
  NC_SAFE_CALL(nc_put_var_int(ncid, varIDmap.find("status")->second, &status[0]));

  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("mean")->second, &mean[0]));
  NC_SAFE_CALL(nc_put_var_float(ncid, varIDmap.find("variance")->second, &variance[0]));

  // Sort out the scalars

  NC_SAFE_CALL(nc_put_var_double(ncid, scalarIDmap.find("exp_k")->second, &(src->exp_k)));
  NC_SAFE_CALL(nc_put_var_int(ncid, scalarIDmap.find("neg")->second, &(src->neg)));

  // Close the file
  NC_SAFE_CALL(nc_close(ncid));

  cout << "complete" << endl;
}

// ---------------

void GCAMorphUtils::Read(GCAM **dst, string fName) const
{
  // Make sure input pointer is NULL
  if (*dst != NULL) {
    cerr << __FUNCTION__ << ": dst pointer not NULL!" << endl;
    exit(EXIT_FAILURE);
  }

  // Construct the filename, using fName passed by value
  fName += ".nc";

  std::cout << __FUNCTION__ << ": Reading file " << fName << " ... ";
  std::cout.flush();

  // Reference for the file
  int ncid;

  // Open the file
  NC_SAFE_CALL(nc_open(fName.c_str(), NC_NOWRITE, &ncid));

  // Query number of variables and dimensions
  int nDimFile, nVarFile;
  NC_SAFE_CALL(nc_inq_ndims(ncid, &nDimFile));
  NC_SAFE_CALL(nc_inq_nvars(ncid, &nVarFile));

  if (nDimFile != static_cast< int >(this->nDims)) {
    cerr << "Invalid number of dimensions " << nDimFile << endl;
    exit(EXIT_FAILURE);
  }

  if (nVarFile != static_cast< int >(this->totalVars)) {
    std::cerr << "Invalid number of variables " << nVarFile << endl;
    exit(EXIT_FAILURE);
  }

  int dimIDs[3];
  size_t dimLen[3];

  NC_SAFE_CALL(nc_inq_dimid(ncid, "x", &dimIDs[this->iX]));
  NC_SAFE_CALL(nc_inq_dimid(ncid, "y", &dimIDs[this->iY]));
  NC_SAFE_CALL(nc_inq_dimid(ncid, "z", &dimIDs[this->iZ]));

  NC_SAFE_CALL(nc_inq_dimlen(ncid, dimIDs[this->iX], &dimLen[this->iX]));
  NC_SAFE_CALL(nc_inq_dimlen(ncid, dimIDs[this->iY], &dimLen[this->iY]));
  NC_SAFE_CALL(nc_inq_dimlen(ncid, dimIDs[this->iZ], &dimLen[this->iZ]));

  // Allocate the target
  *dst = GCAMalloc(dimLen[this->iX], dimLen[this->iY], dimLen[this->iZ]);
  if (*dst == NULL) {
    std::cerr << __FUNCTION__ << ": GCAMalloc failed!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Set this to have only one input
  (*dst)->ninputs = 1;

  // Get hold of the variable IDs
  map< string, int > varIDmap;
  map< string, nc_type >::const_iterator myIt;

  for (myIt = this->varTypeMap.begin(); myIt != this->varTypeMap.end(); myIt++) {
    // Get the variable ID
    NC_SAFE_CALL(nc_inq_varid(ncid,
                              myIt->first.c_str(),  // Name of the variable
                              &varIDmap[myIt->first]));

    // Do a simple sanity check
    nc_type varType;
    NC_SAFE_CALL(nc_inq_vartype(ncid, varIDmap.find(myIt->first)->second, &varType));
    if (varType != myIt->second) {
      cerr << __FUNCTION__ << ": Mismatched type for variable " << myIt->first << endl;
      exit(EXIT_FAILURE);
    }
  }

  // Get hold of the variable IDs for the scalars
  map< string, int > scalarIDmap;
  map< string, nc_type >::const_iterator scalarVarIt;

  for (scalarVarIt = this->scalarTypeMap.begin(); scalarVarIt != this->scalarTypeMap.end(); scalarVarIt++) {
    // Get the variable ID
    NC_SAFE_CALL(nc_inq_varid(ncid,
                              scalarVarIt->first.c_str(),  // Name of the variable
                              &scalarIDmap[scalarVarIt->first]));

    // Do a simple sanity check
    nc_type varType;
    NC_SAFE_CALL(nc_inq_vartype(ncid, scalarIDmap.find(scalarVarIt->first)->second, &varType));
    if (varType != scalarVarIt->second) {
      cerr << __FUNCTION__ << ": Mismatched type for scalar " << scalarVarIt->first << endl;
      exit(EXIT_FAILURE);
    }
  }

  // Read into contiguous arrays
  const size_t nElems = (*dst)->width * (*dst)->height * (*dst)->depth;

  vector< double > x(nElems), y(nElems), z(nElems);
  vector< double > origx(nElems), origy(nElems), origz(nElems);
  vector< float > dx(nElems), dy(nElems), dz(nElems);
  vector< float > odx(nElems), ody(nElems), odz(nElems);
  vector< float > area(nElems), area1(nElems), area2(nElems);
  vector< float > origArea(nElems), origArea1(nElems), origArea2(nElems);
  vector< char > invalid(nElems);
  vector< int > label(nElems), status(nElems);
  vector< float > labelDist(nElems);
  vector< float > mean(nElems), variance(nElems);

  // We use 'find' to get an exception if the name doesn't exist

  NC_SAFE_CALL(nc_get_var_double(ncid, varIDmap.find("rx")->second, &x[0]));
  NC_SAFE_CALL(nc_get_var_double(ncid, varIDmap.find("ry")->second, &y[0]));
  NC_SAFE_CALL(nc_get_var_double(ncid, varIDmap.find("rz")->second, &z[0]));

  NC_SAFE_CALL(nc_get_var_double(ncid, varIDmap.find("origx")->second, &origx[0]));
  NC_SAFE_CALL(nc_get_var_double(ncid, varIDmap.find("origy")->second, &origy[0]));
  NC_SAFE_CALL(nc_get_var_double(ncid, varIDmap.find("origz")->second, &origz[0]));

  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("dx")->second, &dx[0]));
  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("dy")->second, &dy[0]));
  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("dz")->second, &dz[0]));

  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("odx")->second, &odx[0]));
  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("ody")->second, &ody[0]));
  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("odz")->second, &odz[0]));

  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("origArea")->second, &origArea[0]));
  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("origArea1")->second, &origArea1[0]));
  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("origArea2")->second, &origArea2[0]));

  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("area")->second, &area[0]));
  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("area1")->second, &area1[0]));
  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("area2")->second, &area2[0]));

  NC_SAFE_CALL(nc_get_var_text(ncid, varIDmap.find("invalid")->second, &invalid[0]));

  NC_SAFE_CALL(nc_get_var_int(ncid, varIDmap.find("label")->second, &label[0]));
  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("labelDist")->second, &labelDist[0]));
  NC_SAFE_CALL(nc_get_var_int(ncid, varIDmap.find("status")->second, &status[0]));

  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("mean")->second, &mean[0]));
  NC_SAFE_CALL(nc_get_var_float(ncid, varIDmap.find("variance")->second, &variance[0]));

  // Unpack into the GCA_MORPH structure
  for (int i = 0; i < (*dst)->width; i++) {
    for (int j = 0; j < (*dst)->height; j++) {
      for (int k = 0; k < (*dst)->depth; k++) {
        // Save some typing
        GCA_MORPH_NODE &gcamn = (*dst)->nodes[i][j][k];
        // Compute 1D index
        const size_t i1d = k + ((*dst)->depth * (j + ((*dst)->height * i)));

        gcamn.x = x.at(i1d);
        gcamn.y = y.at(i1d);
        gcamn.z = z.at(i1d);

        gcamn.origx = origx.at(i1d);
        gcamn.origy = origy.at(i1d);
        gcamn.origz = origz.at(i1d);

        gcamn.dx = dx.at(i1d);
        gcamn.dy = dy.at(i1d);
        gcamn.dz = dz.at(i1d);

        gcamn.odx = odx.at(i1d);
        gcamn.ody = ody.at(i1d);
        gcamn.odz = odz.at(i1d);

        gcamn.orig_area = origArea.at(i1d);
        gcamn.orig_area1 = origArea1.at(i1d);
        gcamn.orig_area2 = origArea2.at(i1d);

        gcamn.area = area.at(i1d);
        gcamn.area1 = area1.at(i1d);
        gcamn.area2 = area2.at(i1d);

        gcamn.invalid = invalid.at(i1d);
        gcamn.label = label.at(i1d);
        gcamn.status = status.at(i1d);
        gcamn.label_dist = labelDist.at(i1d);

        // Mean and variance require special handling
        if (variance.at(i1d) >= 0) {
          gcamn.gc = alloc_gcs(1, 0, 1);
          gcamn.gc->means[0] = mean.at(i1d);
          gcamn.gc->covars[0] = variance.at(i1d);
        }
        else {
          gcamn.gc = NULL;
        }
      }
    }
  }

  // Fetch the scalars
  NC_SAFE_CALL(nc_get_var_double(ncid, scalarIDmap.find("exp_k")->second, &((*dst)->exp_k)));
  NC_SAFE_CALL(nc_get_var_int(ncid, scalarIDmap.find("neg")->second, &((*dst)->neg)));

  NC_SAFE_CALL(nc_close(ncid));

  cout << "complete" << endl;
}

// #########################################################################

static GCAMorphUtils myGCAMutils;

// ======================================================================

void WriteGCAMoneInput(const GCAM *src, const char *fName) { myGCAMutils.Write(src, fName); }
