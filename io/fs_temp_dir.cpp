#include <iostream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "utils.h"
#include "argparse.h"

#include "fs_temp_dir.help.xml.h"


int main(int argc, char **argv) {

  ArgumentParser parser;
  parser.addHelp(fs_temp_dir_help_xml, fs_temp_dir_help_xml_len);
  parser.addArgument("-b", "--base", 1);
  parser.addArgument("--scratch");
  parser.parse(argc, argv);

  if (parser.exists("scratch")) setenv("FS_TMPDIR", "/scratch", 0);
  if (parser.exists("base")) {
    std::string basedir = parser.retrieve<std::string>("base");
    setenv("FS_TMPDIR", basedir.c_str(), 1);
  }

  std::string tempdir = makeTempDir();
  std::cout << tempdir << "\n";
}
