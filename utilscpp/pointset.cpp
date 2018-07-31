#include <fstream>

#include "pointset.hpp"


/// Writes pointset to a JSON file
/// todo\ support the .dat output format as well and
/// use json.hpp instead of hard-coding this stuff

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
