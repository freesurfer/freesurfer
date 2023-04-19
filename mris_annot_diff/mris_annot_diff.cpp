#include "argparse.h"
#include "log.h"
#include "colortab.h"
#include "mrisurf.h"
#include "fio.h"

#include "mris_annot_diff.help.xml.h"
 
// the program only works with .annot
// --diff-ctab to compare ctab
int main(int argc, char *argv[])
{
  ArgumentParser parser;
  parser.addHelp(mris_annot_diff_help_xml, mris_annot_diff_help_xml_len);
  parser.addArgument("annot1");
  parser.addArgument("annot2");
  parser.addArgument("--verbose",   0, Bool, false);
  parser.addArgument("--diff-ctab", 0, Bool, false);
  parser.parse(argc, argv);

  const std::string fannot1 = parser.retrieve<std::string>("annot1");
  const std::string fannot2 = parser.retrieve<std::string>("annot2");

  // get nvertices from annot files
  FILE *fp1 = fopen(fannot1.c_str(), "r");
  int nvertices1 = freadInt(fp1);
  fclose(fp1);

  FILE *fp2 = fopen(fannot2.c_str(), "r");
  int nvertices2 = freadInt(fp2);
  fclose(fp2); 

  if (nvertices1 != nvertices2)
    fs::fatal() << "annotation sizes (" << nvertices1 << " and " << nvertices2 << ") do not match";

  bool verbose = parser.retrieve<bool>("verbose");

  // load annot files
  MRIS *surf1 = MRISoverAlloc(nvertices1, 0, nvertices1, 0);
  MRIS *surf2 = MRISoverAlloc(nvertices2, 0, nvertices2, 0);

  MRISreadAnnotation(surf1, fannot1.c_str());
  MRISreadAnnotation(surf2, fannot2.c_str());

  // diff labels
  int ndiffs = 0;
  for (int vtxno = 0; vtxno < nvertices1; vtxno++)
  {
    int ano1 = surf1->vertices[vtxno].annotation;
    int ano2 = surf2->vertices[vtxno].annotation;
    if (ano1 != ano2)
    {
      ndiffs++;

      if (verbose)
      {
        int r1, g1, b1;
        MRISAnnotToRGB(ano1, r1, g1, b1);

        int r2, g2, b2;
        MRISAnnotToRGB(ano2, r2, g2, b2);

        printf("[INFO] annotation differs at vtxno %d (%d and %d), (%d,%d,%d and %d,%d,%d)\n",
               vtxno, ano1, ano2, r1, g1, b1, r2, g2, b2);
      }
    }
  }

  int error = 0;
  if (ndiffs > 0)
  {
    fs::error() << "found " << ndiffs << " differences between annotations";
    error = 1;
  }
  else
    printf("No Annotation difference found!\n");

  bool diff_ctab = parser.retrieve<bool>("diff-ctab");
  if (!diff_ctab)
    exit(error);

  printf("\n");

  int nctabdiffs = 0;
  COLOR_TABLE *ctab1 = surf1->ct;
  COLOR_TABLE *ctab2 = surf2->ct;
  int todiff = ctab1->nentries;
  if (ctab1->nentries != ctab2->nentries)
  {
    nctabdiffs++;
    todiff = (ctab1->nentries < ctab2->nentries) ? ctab1->nentries : ctab2->nentries;
    error = 1;
    if (verbose)
      printf("[INFO] colortab number of entries (%d and %d) do not match\n", ctab1->nentries, ctab2->nentries);
  }

  for (int n = 0; n < todiff; n++)
  {
    COLOR_TABLE_ENTRY *cte1 = ctab1->entries[n];
    COLOR_TABLE_ENTRY *cte2 = ctab2->entries[n];
    if ((cte1 == NULL && cte2 != NULL) || (cte1 != NULL && cte2 == NULL))
    {
      nctabdiffs++;
      if (verbose)
        printf("[INFO] colortab entry %d, don't match\n", n);
    }
    else if (cte1 != NULL && cte2 != NULL)
    {
      if (strcmp(cte1->name, cte2->name) != 0)
      {
        nctabdiffs++;
        if (verbose)
          printf("[INFO] colortab entry %d, label name (%s and %s) don't match\n", n, cte1->name, cte2->name);
      }

      if (cte1->ri != cte2->ri ||
          cte1->gi != cte2->gi ||
          cte1->bi != cte2->bi ||
          cte1->ai != cte2->ai)
      {
        nctabdiffs++;
        if (verbose) 
          printf("[INFO] colortab entry %d, RGBi (%d,%d,%d,%d and %d,%d,%d,%d) don't match\n", n, cte1->ri, cte1->gi, cte1->bi, cte1->ai, cte2->ri, cte2->gi, cte2->bi, cte2->ai);
      }
      
      if (cte1->rf != cte2->rf ||
          cte1->gf != cte2->gf ||
          cte1->bf != cte2->bf ||
          cte1->af != cte2->af)
      {
        nctabdiffs++;
        if (verbose)
          printf("[INFO] colortab entry %d, RGBf (%.6f,%.6f,%.6f,%.6f and %.6f,%.6f,%.6f,%.6f) don't match\n", n, cte1->rf, cte1->gf, cte1->bf, cte1->af, cte2->rf, cte2->gf, cte2->bf, cte2->af);
      }
    }
  }

  if (nctabdiffs > 0)
    fs::fatal() << "found " << nctabdiffs << " differences between color tables";
  else
    printf("No color tables difference found!\n");

  exit(0);

#if 0
  // load input
  std::vector<int> annot1 = readAnnotationIntoVector(parser.retrieve<std::string>("annot1"));
  std::vector<int> annot2 = readAnnotationIntoVector(parser.retrieve<std::string>("annot2"));

  // compare lengths
  int length = annot1.size();
  if (length != annot2.size()) fs::fatal() << "annotation sizes (" << length << " and " << annot2.size() << ") do not match";

  // diff labels
  int ndiffs = 0;
  for (int i = 0; i < length; i++) {
    if (annot1[i] != annot2[i]) ndiffs++;
  }

  if (ndiffs > 0) fs::fatal() << "found " << ndiffs << " differences between annotations";
  exit(0);
#endif
}
