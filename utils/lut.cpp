#include <fstream>
#include <sstream>

#include "lut.h"


LookupTable::LookupTable(std::string filename)
{
  importFromFile(filename);
}


void LookupTable::add(int label, std::string labelname)
{
  LabelInfo labelinfo = {true, labelname};
  labelmap[label] = labelinfo;
}


bool LookupTable::importFromFile(std::string filename)
{
  std::ifstream infile(filename);
  std::string line;

  while (std::getline(infile, line)) {
    std::istringstream ss(line);
    LabelInfo labelinfo;
    int label;
    if (!(ss >> label)) continue;
    ss >> labelinfo.name;
    labelmap[label] = labelinfo;
  }

  return true;
}


std::vector<int> LookupTable::labels()
{
  std::vector<int> labellist;
  for(auto const& l: labelmap) labellist.push_back(l.first);
  return labellist;
}


bool LookupTable::empty()
{
  return labelmap.empty();
}


bool LookupTable::hasNameInfo()
{
  bool has_names = false;
  for (auto const &label : labelmap) if (!label.second.name.empty()) has_names = true;
  return has_names;
}
