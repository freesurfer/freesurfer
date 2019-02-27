#ifndef LUT_H
#define LUT_H

#include <string>
#include <map>
#include <vector>


class LookupTable {

  struct LabelInfo {
    bool valid;
    std::string name;
  };

public:
  LookupTable() {};
  LookupTable(std::string filename);

  void add(int label, std::string labelname = "");
  bool importFromFile(std::string filename);

  std::vector<int> labels();

  bool hasNameInfo();
  bool hasColorInfo();
  bool empty();

  // extensions of std::map
  LabelInfo& operator[](int index) { return labelmap[index]; }
  typedef std::map<int, LabelInfo>::iterator iterator;
  typedef std::map<int, LabelInfo>::const_iterator const_iterator;
  iterator begin() { return labelmap.begin(); }
  iterator end() { return labelmap.end(); }

private:
  std::map<int, LabelInfo> labelmap;
};

#endif
