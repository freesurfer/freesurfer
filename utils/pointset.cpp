#include <fstream>

#include "json.h"
#include "pointset.h"


/// Writes pointset to a JSON file
/// TODO: support the .dat output format as well and
/// use json.h instead of hard-coding this stuff

bool PointSet::save(std::string filename)
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
bool PointSet::save_as_ctrlpoint(std::string filename)
{
  std::ofstream ctrpfile;
  ctrpfile.open(filename);

  for (const Point& point : points)
    ctrpfile << point.x << " " << point.y << " " << point.z << "\n";
  ctrpfile << "info\n";
  ctrpfile << "numpoints" << points.size() << "\n";
  ctrpfile << "useRealRAS 1\n";

  ctrpfile.close();

  return true;
}


PointSet loadPointSet(std::string filename)
{
  // read the JSON file from stream
  std::ifstream is(filename);
  nlohmann::json j;
  is >> j;

  // init pointset 
  PointSet ps = PointSet();

  // get space
  ps.vox2ras = j["vox2ras"];

  // get point data
  for (nlohmann::json point : j["points"]) {
    nlohmann::json coord = point["coordinates"];
    ps.add(coord["x"], coord["y"], coord["z"], point["legacy_stat"]);
  }

  return ps;
}
