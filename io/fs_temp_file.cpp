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
  parser.addArgument("-b", "--base", 1);
  parser.addArgument("--scratch");
  parser.parse(argc, argv);

  if (parser.exists("scratch")) setenv("FS_TMPDIR", "/scratch", 0);
  if (parser.exists("base")) {
    std::string basedir = parser.retrieve<std::string>("base");
    setenv("FS_TMPDIR", basedir.c_str(), 1);
  }

  std::string suffix = parser.retrieve<std::string>("suffix");
  std::string tempfile = makeTempFile(suffix);
  std::cout << tempfile << "\n";
}
