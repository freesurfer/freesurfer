/**
 * @brief utility to read an .xml file and output it in a readable format
 *
 * routines for printing of help text formated in xml
 */
/*
 * Original Author: Greg Terrono
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

//-------------------------------------------------------------------

#include <ctype.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xinclude.h>
#include <libxml/xmlIO.h>
#include <string.h>

static void printName(xmlNodePtr cur);
static void printContents(xmlDocPtr doc, xmlNodePtr cur);
static int tagNameIs(const char *c, xmlNodePtr cur);
static int wrdLen(char *c);

int outputHelpDoc(const xmlDocPtr doc)
// output the help text from the xml Doc
{
  xmlNodePtr cur = xmlDocGetRootElement(doc);

  // Catches the error of passing a blank document
  if (cur == NULL) {
    fprintf(stderr, "Empty document.\n");
    xmlFreeDoc(doc);
    return -1;
  }

  // Catches the error of passing a document with a different root tag
  if (!tagNameIs("help", cur)) {
    fprintf(stderr, "Document id of the wrong type!  The root node is not help.\n");
    xmlFreeDoc(doc);
    return -1;
  }

  cur = cur->xmlChildrenNode;

  printf("\t\t\t\tHelp");

  int exampleNum = 1;
  // Outputs the file
  while (cur != NULL) {
    // Outputs the arugments of the file
    if (tagNameIs("arguments", cur)) {
      xmlNodePtr argumentType;
      argumentType = cur->xmlChildrenNode;
      while (argumentType != NULL) {
        if (!tagNameIs("text", argumentType)) {
          printName(argumentType);
          printf(" ARGUMENTS");
          xmlNodePtr argumentElement;
          argumentElement = argumentType->xmlChildrenNode;
          int first = 1;
          while (argumentElement != NULL) {
            if (!tagNameIs("text", argumentElement)) {
              if (tagNameIs("intro", argumentElement) || tagNameIs("argument", argumentElement)) {
                if (first == 1) {
                  first = 0;
                }
                else {
                  printf("\n");
                }
              }
              printContents(doc, argumentElement);
            }
            argumentElement = argumentElement->next;
          }
        }
        argumentType = argumentType->next;
      }
    }
    // Outputs all the other tags
    else if (!tagNameIs("text", cur)) {
      printName(cur);

      if (tagNameIs("example", cur)) {
        printf(" %d", exampleNum++);
      }
      if (tagNameIs("outputs", cur)) {
        xmlNodePtr outputElement;
        outputElement = cur->xmlChildrenNode;
        int first = 1;
        while (outputElement != NULL) {
          if (!tagNameIs("text", outputElement)) {
            if (tagNameIs("intro", outputElement) || tagNameIs("output", outputElement)) {
              if (first == 1) {
                first = 0;
              }
              else {
                printf("\n");
              }
            }
            printContents(doc, outputElement);
          }
          outputElement = outputElement->next;
        }
      }
      else {
        printContents(doc, cur);
      }
    }
    cur = cur->next;
  }

  printf("\n\n\n");

  return 0;  // success
}

static int FileExists(const char *fname)
{
  FILE *fp = fopen(fname, "r");
  if (fp) {
    fclose(fp);
  }
  return (fp != NULL);
}

int outputHelp(const char *name)
// load and parse xml doc from file
{
  xmlDocPtr doc;
  char *fname = (char *)name;

  // Checks for the .xml file name
  if (name == NULL) {
    fprintf(stderr, "No filename passed.\n");
    return -1;
  }

  if (!FileExists(fname)) {
    // assume just the program name was specified, and try to find
    // its associated installed .help.xml file in the freesurfer install dir
    char *fshome = getenv("FREESURFER_HOME");
    if (NULL == fshome) {
      return -1;
    }
    char *fullname = (char *)malloc(strlen(name) + strlen(fshome) + strlen("/docs/xml/.help.xml") + 1);
    strcpy(fullname, fshome);
    strcat(fullname, "/docs/xml/");
    strcat(fullname, name);
    strcat(fullname, ".help.xml");
    fname = (char *)fullname;
  }

  // Opens the .xml file given as an argument
  doc = xmlParseFile(fname);
  // Catches an error in parsing the file
  if (doc == NULL) {
    fprintf(stderr, "Document not parsed successfully. \n");
    return -1;
  }

  return outputHelpDoc(doc);
}

int outputHelpXml(const unsigned char *text, unsigned int size)
// load and parse xml doc from memory
{
  xmlDocPtr doc;

  // Checks if really passed:
  if (text == NULL) {
    fprintf(stderr, "No xml string passed.\n");
    return -1;
  }

  // Parses the .xml file given as an argument
  doc = xmlParseMemory((char *)text, size);
  // Catches an error in parsing the file
  if (doc == NULL) {
    fprintf(stderr, "Document not parsed successfully. \n");
    return -1;
  }

  return outputHelpDoc(doc);
}

// Changes a string to upper case
static void toUpperCase(char *c)
{
  int i;
  for (i = 0; *(c + i) != '\0'; i++) {
    *(c + i) = toupper(*(c + i));
  }
}

// Replaces the first underscore or dash in a string with a space
static void replaceUnderscore(char *c)
{
  char *pos = strchr(c, '_');
  if (pos != NULL) {
    *pos = ' ';
  }
  pos = strchr(c, '-');
  if (pos != NULL) {
    *pos = ' ';
  }
}

// Prints the name of a tag in the correct format
static void printName(xmlNodePtr cur)
{
  char *n = (char *)malloc(strlen((char *)cur->name) + 1);
  strcpy(n, (char *)cur->name);
  toUpperCase(n);
  replaceUnderscore(n);
  printf("\n\n%s", n);
}

// Prints the text of a tag in the correct format
// at most FSPRINT_MAX_CHARS characters per line with the correct number of tabs
#define FSPRINT_MAX_CHARS 78
static void printContents(xmlDocPtr doc, xmlNodePtr cur)
{
  xmlChar *contents;
  contents = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
  unsigned int i = 0, j;
  while (i < strlen((char *)contents)) {
    int tabNum = 1;
    printf("\n\t");
    if (tagNameIs("explanation", cur)) {
      printf("\t");
      tabNum = 2;
    }
    if (*(contents + i) == ' ') {
      i++;
    }
    for (j = i; j < FSPRINT_MAX_CHARS - tabNum * 8 + i && j < strlen((char *)contents) &&
                (wrdLen((char *)(contents + j)) > FSPRINT_MAX_CHARS - tabNum * 8 ||
                 (unsigned)wrdLen((char *)(contents + j)) <= FSPRINT_MAX_CHARS - tabNum * 8 + i - j);
         j++) {
      if (*(contents + j) == '\n') {
        j++;
        break;
      }
      else if (wrdLen((char *)(contents + j)) > FSPRINT_MAX_CHARS - tabNum * 8) {
        for (; j < FSPRINT_MAX_CHARS - tabNum * 8 + i; j++) {
          printf("%c", *(contents + j));
        }
      }
      else {
        printf("%c", *(contents + j));
      }
    }
    i = j;
  }
  xmlFree(contents);
}

// Checks if the name of the tag is the same as the string parameter
static int tagNameIs(const char *c, xmlNodePtr cur) { return !(xmlStrcasecmp(cur->name, (const xmlChar *)c)); }

// Counts the number of chars until the next space, dash, underscore, slash,
// or end of the string
// Used to keep a word from running on two lines
static int wrdLen(char *c)
{
  unsigned int i;
  for (i = 0; i < strlen(c); i++) {
    if (*(c + i) == ' ' || *(c + i) == '_' || *(c + i) == '/') {
      return i;
    }
  }
  return i;
}

#ifdef BUILD_MAIN
//-------------------------------------------------------------------
int main(int argc, char *argv[])
{
  if (argv[1] == NULL) {
    // assuming input is being piped (ie 'cat somefile.help.xml | fsPrintHelp')
    return outputHelp("/dev/stdin");
  }
  return outputHelp(argv[1]);
  exit(0);
}
#endif
