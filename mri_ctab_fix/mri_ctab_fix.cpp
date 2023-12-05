#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>

#include "colortab.h"
#include "log.h"

void parse_commandline(int argc, char **argv);
void print_usage(const char *prog);
std::map<std::string, int> buildLabelAnnotationMap(const COLOR_TABLE *ctab);
void checkColortab(const char *inctabfile, const char *outctabfile);
COLOR_TABLE *splitColortab(std::map<std::string, int>* label_map, COLOR_TABLE *ctab_merged);
int fix1(std::map<std::string, int>* label_map_1, std::map<std::string, int>* label_map_2, COLOR_TABLE *inctab_1, COLOR_TABLE *inctab_2);
COLOR_TABLE *mergeColortab(std::map<std::string, int>* label_map_1, std::map<std::string, int>* label_map_2, COLOR_TABLE *inctab_1, COLOR_TABLE *inctab_2);

const char *exception = NULL;
const char *mode = NULL;
int check = 0, fix = 0, merge = 0;
const char *inctab_1file = NULL, *inctab_2file = NULL, *mergedctabfile = NULL, *ctab_mergedfile = NULL;
const char *outfixedctab_1file = NULL, *outfixedctab_2file = NULL;

/* 
 * usage: mri_ctab_fix
 *             -c [-e <exception>] <inctab> [outfixedctab]
 *             -1 [-e <exception>] <inctab_1> <inctab_2> [outfixedctab_1 outfixedctab_2]
 *             -2 <inctab_1> <inctab_2> <ctab_merged> [outfixedctab_1 outfixedctab_2]
 *             -m <inctab_1> <inctab_2> [mergedctab]
 *
 *     mri_ctab_fix -c [-e <exception>] <inctab> [outfixedctab]
 *       check duplicated color scheme assignment in <inctab>
 *       output fixed LUT in outfixedctab if it is given
 *       use optional <exception> to specify label to skip
 *
 *     mri_ctab_fix -1 [-e <exception>] <inctab_1> <inctab_2> [outfixedctab_1 outfixedctab_2]
 *       inctab_1 is unchanged, inctab_2 is fixed so that
 *       same labels in inctab_1 and inctab_2 share the same annotations
 *       different labels across inctab_1 and inctab_2 have different annotations
 *       output fixed LUTs in outfixedctab_1 and outfixedctab_2 if they are given
 *       use optional <exception> to specify label to skip
 *
 *     mri_ctab_fix -m <inctab_1> <inctab_2> [mergedctab]
 *       make sure there are no duplicated annotations across two ctab files before running the merge
 *       merge inctab_1 and inctab_2, both inctab_1 and inctab_2 are unchanged;
 *       all labels in <inctab_1> are copied into mergedct, label ids are kept;
 *       any labels in <inctab_2> not in <inctab_1> are sequentially re-numbered after inctab_1
 *       merged ctab is saved as mergedctab, or printed to stdout
 *
 *     mri_ctab_fix -2 <inctab_1> <inctab_2> <ctab_merged> <outfixedctab_1> <outfixedctab_2>
 *       split <ctab_merged> into 2 LUTs with original labels found in <inctab_1>/<inctab_2>
 *       output split LUTs in outfixedctab_1 and outfixedctab_2
 *       after the split, same annotations across outfixedctab_1 and outfixedctab_2 share same label names and ids
 *
 */
int main(int argc, char *argv[])
{
  parse_commandline(argc, argv);

  if (strcmp(mode, "-c") == 0)
  {
    check = 1;
    printf("[INFO] check (-c %s)\n", inctab_1file);
    checkColortab(inctab_1file, inctab_2file);
    exit(0);
  }
  
  // build label name => annotation map for inctab_1
  COLOR_TABLE *inctab_1 = CTABreadASCII(inctab_1file);
  if (inctab_1 == NULL)
  {
    print_usage(argv[0]);
    exit(1);
  }
  std::map<std::string, int> label_map_1 = buildLabelAnnotationMap(inctab_1);

  // build label name => annotation map for inctab_2
  COLOR_TABLE *inctab_2 = CTABreadASCII(inctab_2file);  
  if (inctab_2 == NULL)
  {
    print_usage(argv[0]);
    exit(1);
  }
  std::map<std::string, int> label_map_2 = buildLabelAnnotationMap(inctab_2);

  if (merge)
  {
    printf("[INFO] merge (-m %s %s %s)\n", inctab_1file, inctab_2file, (mergedctabfile != NULL) ? mergedctabfile : "");
    COLOR_TABLE *mergedct = mergeColortab(&label_map_1, &label_map_2, inctab_1, inctab_2);
    if (mergedct == NULL)
    {
      print_usage(argv[0]);
      exit(1);
    }

    // output mergedct, print to stdout by default
    FILE *mergedfp = stdout;
    if (mergedctabfile != NULL)
    {
      mergedfp = fopen(mergedctabfile, "w");
      printf("output new ctab %s ...\n", mergedctabfile);
    }
    
    CTABprintASCII(mergedct, mergedfp);
    if (mergedctabfile != NULL)
      fclose(mergedfp);
    
    exit(0);
  }

  if (fix == 2)
  {
    printf("[INFO] fix2 (-2 %s %s %s)\n", inctab_1file, inctab_2file, ctab_mergedfile);
  
    COLOR_TABLE *ctab_merged = CTABreadASCII(ctab_mergedfile);
    if (ctab_merged == NULL)
    {
      print_usage(argv[0]);
      exit(1);
    }
    
    COLOR_TABLE *ct1 = splitColortab(&label_map_1, ctab_merged);
    // output fixed inctab_1      
    printf("output new ctab %s ...\n", outfixedctab_1file);
    FILE *newfp_ctab1 = fopen(outfixedctab_1file, "w");
    CTABprintASCII(ct1, newfp_ctab1);
    fclose(newfp_ctab1);
    
    COLOR_TABLE *ct2 = splitColortab(&label_map_2, ctab_merged);
    // output fixed inctab_2
    printf("output new ctab %s ...\n", outfixedctab_2file);
    FILE *newfp_ctab2 = fopen(outfixedctab_2file, "w");
    CTABprintASCII(ct2, newfp_ctab2);
    fclose(newfp_ctab2);
  }
  else if (fix == 1)
  {
    printf("[INFO] fix1 (-1 %s %s)\n", inctab_1file, inctab_2file);
    int fixcount = fix1(&label_map_1, &label_map_2, inctab_1, inctab_2);

    if (fixcount > 0)
    {
      if (outfixedctab_1file != NULL)
      {
       // inctab_1 is unchanged, output it anyway
        printf("\noutput new ctab %s ...\n", outfixedctab_1file);
        FILE *newfp_ctab1 = fopen(outfixedctab_1file, "w");
        CTABprintASCII(inctab_1, newfp_ctab1);
        fclose(newfp_ctab1);
      }

      if (outfixedctab_2file != NULL)
      {
        // inctab_2 is changed in place
        printf("\noutput new ctab %s ...\n", outfixedctab_2file);
        FILE *newfp_ctab2 = fopen(outfixedctab_2file, "w");
        CTABprintASCII(inctab_2, newfp_ctab2);
        fclose(newfp_ctab2);
      }
    }

    printf("\n\n Fixed Annotation Assignments: %d\n", fixcount);
  }
  
  exit(0);
}

// create a new COLOR_TABLE from ctab_merged with only labels found in label_map, keep the same label ids
// after the split: same annotations across two LUTs share same label names and ids
COLOR_TABLE *splitColortab(std::map<std::string, int>* label_map, COLOR_TABLE *ctab_merged)
{  
  /* Allocate our table. */
  COLOR_TABLE *splitct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (splitct == NULL) ErrorReturn(NULL, (ERROR_NO_MEMORY, "splitColortab: could not allocate table"));

  /* Allocate the array of NULL CTE ptrs. */
  splitct->nentries = ctab_merged->nentries;
  splitct->entries = (COLOR_TABLE_ENTRY **)calloc(splitct->nentries, sizeof(COLOR_TABLE_ENTRY *));
  if (splitct->entries == NULL)
    ErrorReturn(NULL, (ERROR_NO_MEMORY, "mergeColortab: could not allocate %d entries", splitct->nentries));

  /* Copy in the file name. */
  strncpy(splitct->fname, ctab_merged->fname, sizeof(splitct->fname));

  /* We'll write this version if we write to binary. */
  splitct->version = CTAB_VERSION_TO_WRITE;

  // copy the tissue type ctab if it is there
  splitct->ctabTissueType = ctab_merged->ctabTissueType;

  // build label name => annotation map for ctab_merged
  std::map<std::string, int> label_map_merged = buildLabelAnnotationMap(ctab_merged);

  // loop through label_map_merged
  // if it is found in label_map, output in splitct, keep the same label id and annotation
  std::map<std::string, int>::const_iterator label_map_merged_It;
  for (label_map_merged_It = label_map->begin(); label_map_merged_It != label_map->end(); label_map_merged_It++)
  {
    std::string label_merged = label_map_merged_It->first;
    //int annot_merged = label_map_merged_It->second;

    // try to find it in label_map
    std::map<std::string, int>::const_iterator label_map_It = label_map->find(label_merged);
    if (label_map_It != label_map->end())
    {
      // found it
      int ctab_merged_index = -1;
      CTABfindName(ctab_merged, label_merged.c_str(), &ctab_merged_index);
      if (ctab_merged_index < 0)
        printf("[ERROR] label %s not found in ctab_merged %s\n", label_merged.c_str(), ctab_merged->fname);

      if (splitct->entries[ctab_merged_index] != NULL)
        printf("[WARN] duplicated label id: %d %s found in ctab_merged %s\n", ctab_merged_index, label_merged.c_str(), ctab_merged->fname);
      else
      {
        splitct->entries[ctab_merged_index] = (CTE *)malloc(sizeof(CTE));
        // check malloc return ...

        // fill the entry
        strncpy(splitct->entries[ctab_merged_index]->name, ctab_merged->entries[ctab_merged_index]->name, sizeof(splitct->entries[ctab_merged_index]->name));
        splitct->entries[ctab_merged_index]->ri = ctab_merged->entries[ctab_merged_index]->ri;
        splitct->entries[ctab_merged_index]->gi = ctab_merged->entries[ctab_merged_index]->gi;
        splitct->entries[ctab_merged_index]->bi = ctab_merged->entries[ctab_merged_index]->bi;
        splitct->entries[ctab_merged_index]->ai = ctab_merged->entries[ctab_merged_index]->ai;

        /* Now calculate the float versions. */
        splitct->entries[ctab_merged_index]->rf = (float)splitct->entries[ctab_merged_index]->ri / 255.0;
        splitct->entries[ctab_merged_index]->gf = (float)splitct->entries[ctab_merged_index]->gi / 255.0;
        splitct->entries[ctab_merged_index]->bf = (float)splitct->entries[ctab_merged_index]->bi / 255.0;
        splitct->entries[ctab_merged_index]->af = (float)splitct->entries[ctab_merged_index]->ai / 255.0;
        splitct->entries[ctab_merged_index]->TissueType = ctab_merged->entries[ctab_merged_index]->TissueType;
        splitct->entries[ctab_merged_index]->count = 0;	
      }
    }
  }

  return splitct;
}


// inctab_1 is unchanged, inctab_2 is fixed so that
//   same labels in inctab_1 and inctab_2 share the same annotation
//   different labels across inctab_1 and inctab_2 have different annotations
//   (annotations across inctab_1 and inctab_2 are unique)
int fix1(std::map<std::string, int>* label_map_1, std::map<std::string, int>* label_map_2, COLOR_TABLE *inctab_1, COLOR_TABLE *inctab_2)
{
  int fixcount = 0;
  
  // loop through label_map_1, search them in label_map_2
  // same label, same annotation
  std::map<std::string, int>::const_iterator label_map_1_It;
  for (label_map_1_It = label_map_1->begin(); label_map_1_It != label_map_1->end(); label_map_1_It++)
  {
    std::string label_1 = label_map_1_It->first;
    int annot_1 = label_map_1_It->second;

    std::map<std::string, int>::iterator label_map_2_It = label_map_2->find(label_1);    
    if (label_map_2_It != label_map_2->end())
    {
      // found the same label in both inctab_1 and inctab_2
      // make sure they have the same annotation, fix inctab_2
      if (annot_1 != label_map_2_It->second)
      {
	printf("[INFO] replace inctab_2 label %s (%s) annotation %d => %d\n", label_1.c_str(), inctab_2->fname, label_map_2_It->second, annot_1);
	label_map_2_It->second = annot_1;

	// fix inctab_2
	int inctab_2_index = -1;
	CTABfindName(inctab_2, label_1.c_str(), &inctab_2_index);
	if (inctab_2_index < 0)
	{
	  printf("[ERROR] label %s not found in inctab_2 %s\n", label_1.c_str(), inctab_2->fname);
	  continue;
	}

	int ri, gi, bi;
	AnnotToRGB(annot_1, ri, gi, bi);
	
	// replace the entry with newly calculated annotation
	inctab_2->entries[inctab_2_index]->ri = ri;
	inctab_2->entries[inctab_2_index]->gi = gi;
	inctab_2->entries[inctab_2_index]->bi = bi;

	// the float versions
	inctab_2->entries[inctab_2_index]->rf = (float)ri / 255.0;
	inctab_2->entries[inctab_2_index]->gf = (float)gi / 255.0;
	inctab_2->entries[inctab_2_index]->bf = (float)bi / 255.0;

	fixcount++;
      }
    }
  }

  printf("\n");
  
  // check annotations, same annot must have same label
  for (label_map_1_It = label_map_1->begin(); label_map_1_It != label_map_1->end(); label_map_1_It++)
  {
    std::string label_1 = label_map_1_It->first;
    int annot_1 = label_map_1_It->second;

    std::map<std::string, int>::iterator label_map_2_It;
    for (label_map_2_It = label_map_2->begin(); label_map_2_It != label_map_2->end(); label_map_2_It++)
    {
      std::string label_2 = label_map_2_It->first;
      int annot_2 = label_map_2_It->second;
      if (annot_1 == annot_2 && label_1 != label_2)
      {
	printf("\n");
	
	// same annotation assigned to different label names
        printf("[WARN] the same annot %8d for label %s (%s) and label %s (%s)\n", annot_1, label_1.c_str(), inctab_1->fname, label_2.c_str(), inctab_2->fname);
	if (exception != NULL)
	{
	  std::size_t found = label_1.find(exception);
          if (found != std::string::npos)
            continue;
	}

	 // fix inctab_2
	int inctab_2_index = -1;
	CTABfindName(inctab_2, label_2.c_str(), &inctab_2_index);
	if (inctab_2_index < 0)
	{
	  printf("[ERROR] label %s not found in inctab_2 %s\n", label_2.c_str(), inctab_2->fname);
	  continue;
	}

	// generate a random annotation
	int annot_new = CTABgenUniqueAnnotation(inctab_2, inctab_1);
	label_map_2_It->second = annot_new;

	printf("[INFO] replace inctab_2 label %s (%s) annotation %d => %d\n", label_2.c_str(), inctab_2->fname, annot_2, annot_new);

	int ri, gi, bi;
	AnnotToRGB(annot_new, ri, gi, bi);
	
	// replace the entry with newly calculated annotation
	inctab_2->entries[inctab_2_index]->ri = ri;
	inctab_2->entries[inctab_2_index]->gi = gi;
	inctab_2->entries[inctab_2_index]->bi = bi;

	// the float versions
	inctab_2->entries[inctab_2_index]->rf = (float)ri / 255.0;
	inctab_2->entries[inctab_2_index]->gf = (float)gi / 255.0;
	inctab_2->entries[inctab_2_index]->bf = (float)bi / 255.0;

	fixcount++;	
      }
    }
  }

  return fixcount;
}

// make sure there are no duplicated annotations across two ctab files before running the merge
// all labels in inctab_1 are copied into mergedct, label ids are kept
// any labels in inctab_2 not in inctab_1 are sequentially re-numbered after inctab_1
COLOR_TABLE *mergeColortab(std::map<std::string, int>* label_map_1, std::map<std::string, int>* label_map_2, COLOR_TABLE *inctab_1, COLOR_TABLE *inctab_2)
{
  /* Allocate our table. */
  COLOR_TABLE *mergedct = (COLOR_TABLE *)calloc(1, sizeof(COLOR_TABLE));
  if (mergedct == NULL) ErrorReturn(NULL, (ERROR_NO_MEMORY, "mergeColortab: could not allocate table"));

  /* Allocate the array of NULL CTE ptrs. */
  mergedct->nentries = inctab_1->nentries + inctab_2->nentries;
  mergedct->entries = (COLOR_TABLE_ENTRY **)calloc(mergedct->nentries, sizeof(COLOR_TABLE_ENTRY *));
  if (mergedct->entries == NULL)
    ErrorReturn(NULL, (ERROR_NO_MEMORY, "mergeColortab: could not allocate %d entries", mergedct->nentries));

  /* Copy in the file name. */
  strncpy(mergedct->fname, inctab_1->fname, sizeof(mergedct->fname));
  
  /* We'll write this version if we write to binary. */
  mergedct->version = CTAB_VERSION_TO_WRITE;

  // copy the tissue type ctab if it is there
  mergedct->ctabTissueType = inctab_1->ctabTissueType;

  int maxlabel = 0;
  // loop through label_map_1, add all labels, keep the label ids
  std::map<std::string, int>::const_iterator label_map_1_It;
  for (label_map_1_It = label_map_1->begin(); label_map_1_It != label_map_1->end(); label_map_1_It++)
  {
    std::string label_1 = label_map_1_It->first;

    // keep the label ids in inctab_1
    int inctab_1_index = -1;
    CTABfindName(inctab_1, label_1.c_str(), &inctab_1_index);
    if (inctab_1_index < 0)
      printf("[ERROR] label %s not found in inctab_1 %s\n", label_1.c_str(), inctab_1->fname);

    if (mergedct->entries[inctab_1_index] != NULL)
      printf("[WARN] duplicated label id: %d %s found in inctab_1 %s\n", inctab_1_index, label_1.c_str(), inctab_1->fname);
    else
    {
      maxlabel = (inctab_1_index > maxlabel) ? inctab_1_index : maxlabel;
      mergedct->entries[inctab_1_index] = (CTE *)malloc(sizeof(CTE));
      // check malloc return ...

      // fill the entry
      strncpy(mergedct->entries[inctab_1_index]->name, inctab_1->entries[inctab_1_index]->name, sizeof(mergedct->entries[inctab_1_index]->name));
      mergedct->entries[inctab_1_index]->ri = inctab_1->entries[inctab_1_index]->ri;
      mergedct->entries[inctab_1_index]->gi = inctab_1->entries[inctab_1_index]->gi;
      mergedct->entries[inctab_1_index]->bi = inctab_1->entries[inctab_1_index]->bi;
      mergedct->entries[inctab_1_index]->ai = inctab_1->entries[inctab_1_index]->ai;

      /* Now calculate the float versions. */
      mergedct->entries[inctab_1_index]->rf = (float)mergedct->entries[inctab_1_index]->ri / 255.0;
      mergedct->entries[inctab_1_index]->gf = (float)mergedct->entries[inctab_1_index]->gi / 255.0;
      mergedct->entries[inctab_1_index]->bf = (float)mergedct->entries[inctab_1_index]->bi / 255.0;
      mergedct->entries[inctab_1_index]->af = (float)mergedct->entries[inctab_1_index]->ai / 255.0;
      mergedct->entries[inctab_1_index]->TissueType = inctab_1->entries[inctab_1_index]->TissueType;
      mergedct->entries[inctab_1_index]->count = 0;	
    }
  }

  // loop through label_map_2 labels, add any labels not in label_map_1
  std::map<std::string, int>::const_iterator label_map_2_It;
  for (label_map_2_It = label_map_2->begin(); label_map_2_It != label_map_2->end(); label_map_2_It++)
  {
    std::string label_2 = label_map_2_It->first;
    
    std::map<std::string, int>::iterator map_1_It = label_map_1->find(label_2);
    if (map_1_It == label_map_1->end())
    {
      // labels in ctab2, not in ctab1
      int inctab_2_index = -1;
      CTABfindName(inctab_2, label_2.c_str(), &inctab_2_index);
      if (inctab_2_index < 0)
        printf("[ERROR] label %s not found in inctab_2 %s\n", label_2.c_str(), inctab_2->fname);

      // re-number the label id, add at the end
      maxlabel++;
      int newlabel = maxlabel;

      mergedct->entries[newlabel] = (CTE *)malloc(sizeof(CTE));
      // check malloc return ...

      // fill the entry
      strncpy(mergedct->entries[newlabel]->name, inctab_2->entries[inctab_2_index]->name, sizeof(mergedct->entries[newlabel]->name));
      mergedct->entries[newlabel]->ri = inctab_2->entries[inctab_2_index]->ri;
      mergedct->entries[newlabel]->gi = inctab_2->entries[inctab_2_index]->gi;
      mergedct->entries[newlabel]->bi = inctab_2->entries[inctab_2_index]->bi;
      mergedct->entries[newlabel]->ai = inctab_2->entries[inctab_2_index]->ai;

      /* Now calculate the float versions. */
      mergedct->entries[newlabel]->rf = (float)mergedct->entries[newlabel]->ri / 255.0;
      mergedct->entries[newlabel]->gf = (float)mergedct->entries[newlabel]->gi / 255.0;
      mergedct->entries[newlabel]->bf = (float)mergedct->entries[newlabel]->bi / 255.0;
      mergedct->entries[newlabel]->af = (float)mergedct->entries[newlabel]->ai / 255.0;
      mergedct->entries[newlabel]->TissueType = inctab_2->entries[inctab_2_index]->TissueType;
      mergedct->entries[newlabel]->count = 0;	
    }
  }

  return mergedct;
}


// check duplicated color scheme assignment in inctabfile
// output fixed LUT in outctabfile if it is given
void checkColortab(const char *inctabfile, const char *outctabfile)
{
  COLOR_TABLE *ctab = CTABreadASCII(inctabfile);
  if (ctab == NULL)
  {
    print_usage("mri_ctab_fix");
    exit(1);
  }

  int dupCount = CTABprintAnnotationAssignment(ctab, true, NULL, exception);
  if (dupCount > 0)
  {
    FILE *newfp = stdout;
    if (outctabfile != NULL)
    {
      printf("\n\noutput fixed colortable in %s\n", outctabfile);
      newfp = fopen(outctabfile, "w");
    }
    CTABprintASCII(ctab, newfp);
    if (outctabfile != NULL)
      fclose(newfp);

    printf("\n\n Duplicated Annotation Assignments: %d\n", dupCount);
  }
  else
    printf("\n\n No Duplicated Annotation Assignments Found!\n");
}


// build lable-annotation for given COLOR_TABLE
std::map<std::string, int> buildLabelAnnotationMap(const COLOR_TABLE *ctab)
{
  std::map<std::string, int> label_map;
  for (int labelid = 0; labelid < ctab->nentries; labelid++)
  {
    if (ctab->entries[labelid] != NULL)
    {
      int r = ctab->entries[labelid]->ri;
      int g = ctab->entries[labelid]->gi;
      int b = ctab->entries[labelid]->bi;

      int annotation;
      RGBToAnnot(r, g, b, annotation);
    
      label_map[ctab->entries[labelid]->name] = annotation;
    }
  }

  return label_map;
}


// usage
void print_usage(const char *prog)
{
  printf("\n");  
  printf("%s can be run in 4 modes:\n", prog);
  printf("\n");
  printf("\t-c [-e <exception>] <inctab> [outfixedctab]\n");
  printf("\t   check duplicated color scheme assignment in <inctab>\n");
  printf("\t   output fixed LUT in outfixedctab if it is given\n");
  printf("\t   use optional <exception> to specify label to skip\n");

  printf("\n");  
  printf("\t-1 [-e <exception>] <inctab_1> <inctab_2> [outfixedctab_1 outfixedctab_2]\n");
  printf("\t   <inctab_1> is unchanged, <inctab_2> is fixed so that\n");
  printf("\t   same labels in inctab_1 and inctab_2 share the same annotations;\n");
  printf("\t   different labels across inctab_1 and inctab_2 have different annotations.\n");
  printf("\t   output fixed LUTs in outfixedctab_1 and outfixedctab_2 if they are given\n");
  printf("\t   use optional <exception> to specify label to skip\n");  

  printf("\n");
  printf("\t-2 <inctab_1> <inctab_2> <ctab_merged> <outfixedctab_1> <outfixedctab_2>\n");
  printf("\t   split <ctab_merged> into 2 LUTs with original labels found in <inctab_1>/<inctab_2>\n");
  printf("\t   output split LUTs in outfixedctab_1 and outfixedctab_2\n");
  printf("\t   after the split, same annotations across outfixedctab_1 and outfixedctab_2 share same label names and ids\n");

  printf("\n");  
  printf("\t-m <inctab_1> <inctab_2> [mergedctab]\n");
  printf("\t   make sure there are no duplicated annotations across two ctab files before running the merge\n");
  printf("\t   merge <inctab_1> and <inctab_2>, both <inctab_1> and <inctab_2> are unchanged\n");
  printf("\t   all labels in <inctab_1> are copied into mergedct, label ids are kept;\n");
  printf("\t   any labels in <inctab_2> not in <inctab_1> are re-numbered and put at the end\n");
  printf("\t   merged ctab is saved as mergedctab, or printed to stdout\n");
}


void parse_commandline(int argc, char **argv)
{
  if (argc < 2)
  {
    print_usage(argv[0]);
    exit(1);
  }

  int iargv = 1;
  mode = argv[iargv];
  if (strcmp(mode, "-c") != 0 && strcmp(mode, "-1") != 0 && strcmp(mode, "-2") != 0 && strcmp(mode, "-m") != 0)
  {
    print_usage(argv[0]);
    exit(1);
  }
  
  iargv++;
  if ((argc - iargv) >= 1 && strcmp(argv[iargv], "-e") == 0)
  {
    iargv++;
    if ((argc - iargv) < 1)
    {
      print_usage(argv[0]);
      exit(1);
    }
    exception = argv[iargv];
    iargv++;  // move to next argument
  }

  // need to have at least one argument left
  if ((argc - iargv) < 1)
  {
    print_usage(argv[0]);
    exit(1);
  }
  
  inctab_1file = argv[iargv++];
  if ((argc - iargv) >= 1)
    inctab_2file = argv[iargv++];

  if (strcmp(mode, "-c") == 0)
    return;

  if (strcmp(mode, "-1") == 0)
  {
    fix = 1;
    if (inctab_2file == NULL)
    {
      print_usage(argv[0]);
      exit(1);
    }

    if ((argc - iargv) >= 2)
    {
      outfixedctab_1file = argv[iargv++];
      outfixedctab_2file = argv[iargv++];
    }
  }
  else if (strcmp(mode, "-2") == 0)
  {
    fix = 2;
    if (inctab_2file == NULL || (argc - iargv) < 3)
    {
      print_usage(argv[0]);
      exit(1);
    }

    ctab_mergedfile = argv[iargv++];
    outfixedctab_1file = argv[iargv++];
    outfixedctab_2file = argv[iargv++];
  }
  else if (strcmp(mode, "-m") == 0)
  {
    merge = 1;
    if (inctab_2file == NULL)
    {
      print_usage(argv[0]);
      exit(1);
    }

    if ((argc - iargv) >= 1)
      mergedctabfile = argv[iargv++];
  }
}
