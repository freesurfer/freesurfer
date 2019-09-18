#include "argparse.h"
#include "log.h"
#include "annotation.h"

#include "mris_annot_diff.help.xml.h"
 

int main(int argc, const char *argv[])
{
  ArgumentParser parser;
  parser.addHelp(mris_annot_diff_help_xml, mris_annot_diff_help_xml_len);
  parser.addArgument("annot1");
  parser.addArgument("annot2");
  parser.parse(argc, argv);

  // load input
  std::vector<int> annot1 = readAnnotationIntoVector(parser.retrieve<std::string>("annot1"));
  std::vector<int> annot2 = readAnnotationIntoVector(parser.retrieve<std::string>("annot2"));

  // compare lengths
  int length = annot1.size();
  if (length != annot2.size()) logFatal(1) << "annotation sizes (" << length << " and " << annot2.size() << ") do not match";

  // diff labels
  int ndiffs = 0;
  for (int i = 0; i < length; i++) {
    if (annot1[i] != annot2[i]) ndiffs++;
  }

  if (ndiffs > 0) logFatal(1) << "found " << ndiffs << " differences between annotations";
  exit(0);
}
