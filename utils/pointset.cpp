#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>

#include "json.h"
#include "pointset.h"


/// Writes pointset to a JSON file
/// TODO: support the .dat output format as well and
/// use json.h instead of hard-coding this stuff

bool fsPointSet::save(std::string filename)
{
  std::ofstream jsonfile;
  jsonfile.open(filename);

  jsonfile << "{\n";
  jsonfile << "    \"data_type\": \"fs_pointset\",\n";
  jsonfile << "    \"points\": [\n";

  for (const Point& point : points) {
    jsonfile << "        {\n";
    jsonfile << "            \"coordinates\": {\n";
    jsonfile << "                \"x\": " << point.x << ",\n";
    jsonfile << "                \"y\": " << point.y << ",\n";
    jsonfile << "                \"z\": " << point.z << "\n";
    jsonfile << "            },\n";
    jsonfile << "            \"legacy_stat\": 1\n";
#if 0
    // Currently, this breaks loading into freeview
    if(point.labelname.length()>0)
      jsonfile << "            \"labelname\": \"" << point.labelname << "\"\n";
    jsonfile << "            \"index\": " << point.index << "\n";
    jsonfile << "            \"value\": " << point.value << "\n";
    jsonfile << "            \"count\": " << point.count << "\n";
#endif    
    // the Qt json loader requires EXACT comma syntax
    if (&point == &points.back()) {
      jsonfile << "        }\n";
    } else {
      jsonfile << "        },\n";
    }
  }

  jsonfile << "    ],\n";
  jsonfile << "    \"vox2ras\": \"" << vox2ras << "\"\n";
  jsonfile << "}\n";

  jsonfile.close();

  return true;
}

// Save in control point format, v6 compatible
bool fsPointSet::save_as_ctrlpoint(std::string filename)
{
  std::ofstream ctrpfile;
  ctrpfile.open(filename);

  for (const Point& point : points)
    ctrpfile << point.x << " " << point.y << " " << point.z << "\n";
  ctrpfile << "info\n";
  ctrpfile << "numpoints " << points.size() << "\n";
  ctrpfile << "useRealRAS 1\n";

  ctrpfile.close();

  return true;
}


fsPointSet loadfsPointSet(std::string filename)
{
  // read the JSON file from stream
  std::ifstream is(filename);
  nlohmann::json j;
  is >> j;

  // init pointset 
  fsPointSet ps = fsPointSet();

  // get space
  ps.vox2ras = j["vox2ras"];

  // get point data
  for (nlohmann::json point : j["points"]) {
    nlohmann::json coord = point["coordinates"];
    ps.add(coord["x"], coord["y"], coord["z"], point["legacy_stat"]);
  }

  return ps;
}

bool fsPointSet::writeCentroidTable(std::string outfile,std::string ltafile)
{  
  std::ofstream tablefile(outfile);
  if (!tablefile.is_open()) {
    std::cerr << "could not open writable file at " << outfile << std::endl;
    return(false);
  }
  std::cout << "writing results to " << outfile << std::endl;

  int max_id_chars = 2, max_label_chars = 10, haslabels=0;  // used for table formatting
  for (const Point& point : points) {
    int label_chars = point.labelname.length();
    if(label_chars > 0) haslabels=1;
    if(label_chars > max_label_chars) max_label_chars = label_chars;
    std::ostringstream ss;
    ss << point.index;
    int id_chars = ss.str().length();
    if(id_chars > max_id_chars) max_id_chars = id_chars;
  }

  int precision = 4;
  int cwidth = precision + 8;

  // table header
  tablefile << "#SegCentroid" << std::endl;
  if(ltafile.length()) tablefile << "# Reg " << ltafile << std::endl;
  tablefile << std::setw(max_id_chars) << "# ID";
  if(haslabels) {
    tablefile << std::left << "    " << std::setw(max_label_chars)
              << "StructName" << std::right;
  }
  tablefile << " " << std::setw(cwidth) << "R" + std::string(precision, ' ')
            << " " << std::setw(cwidth) << "A" + std::string(precision, ' ')
            << " " << std::setw(cwidth) << "S" + std::string(precision, ' ')
            << std::endl;

  int table_width = max_id_chars + (cwidth * 3) + 3;
  if(haslabels) table_width += max_label_chars + 4;
  tablefile << "# ";
  tablefile << std::string(table_width, '-') << std::endl;

  // table body
  tablefile << std::fixed << std::setprecision(precision);
  for (const Point& point : points) {
    tablefile << std::setw(max_id_chars) << point.index;
    if (haslabels) {
      tablefile << std::left << "    " << std::setw(max_label_chars) 
                << point.labelname << std::right;
    }
    tablefile << " " << std::setw(cwidth) << point.x
              << " " << std::setw(cwidth) << point.y
              << " " << std::setw(cwidth) << point.z
              << std::endl;
  }

  tablefile.close();

  return true;
}
