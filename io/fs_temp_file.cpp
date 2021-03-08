#include <iostream>
#include <string>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>

#include "utils.h"
#include "argparse.h"

#include "fs_temp_file.help.xml.h"


int main(int argc, char **argv) {

  ArgumentParser parser;
  parser.addHelp(fs_temp_file_help_xml, fs_temp_file_help_xml_len);
  parser.addArgument("-s", "--suffix", 1);
  parser.addArgument("--scratch");
  parser.parse(argc, argv);

  if (parser.exists("scratch")) setenv("FS_TMPDIR", "/scratch", 0);

  std::string suffix = parser.retrieve<std::string>("suffix");
  std::string tempfile = getTempFile(suffix);

  umask(077);

  std::fstream fs;
  fs.open(tempfile, std::ios::out);
  fs.close();

  std::cout << tempfile << "\n";
}
