/**
 * @brief Reads an .xml file and outputs it in html format
 *
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//-------------------------------------------------------------------
void toUpperCase(char *sPtr);
void printName(xmlNodePtr cur, FILE *f);
void replaceUnderscore(char *c);
void printContents(xmlDocPtr doc, xmlNodePtr cur, FILE *f);
void printTableContents(xmlDocPtr doc, xmlNodePtr cur, FILE *f);
int tagNameIs(char *c, xmlNodePtr cur);
int wrdLen(char *c);

//-------------------------------------------------------------------
int main(int argc, char *argv[])
{
  xmlDocPtr doc;
  xmlNodePtr cur;

  // Checks for the .xml file name
  if (argc == 1) {
    fprintf(stderr, "No file name passed.\n");
    return -1;
  }
  // Opens the .xml file given as an argument
  doc = xmlParseFile(argv[1]);
  // Catches an error in parsing the file
  if (doc == NULL) {
    fprintf(stderr, "Document not parsed successfully. \n");
    return -1;
  }

  cur = xmlDocGetRootElement(doc);

  // Catches the error of passing a blank document
  if (cur == NULL) {
    fprintf(stderr, "Empty document.\n");
    xmlFreeDoc(doc);
    return -1;
  }

  // Catches the error of passing a document with a different root tag
  if (!tagNameIs("help", cur)) {
    fprintf(stderr,
            "Document id of the wrong type, "
            "the root node is not help.\n");
    xmlFreeDoc(doc);
    return -1;
  }

  if (argc == 2) {
    printf("No output specified.\nOutput will be sent to output.html.\n");
  }
  FILE *outFile;
  outFile = fopen(argc > 2 ? argv[2] : "output.html", "w");

  cur = cur->xmlChildrenNode;

  fprintf(outFile, "<html>\n<body>\n<basefont face=\"sans-serif\"></basefont>\n");

  int exampleNum = 1;
  // Outputs the file
  while (cur != NULL) {
    // Outputs the arugments of the file
    if (tagNameIs("arguments", cur) || tagNameIs("outputs", cur)) {
      xmlNodePtr argumentType;
      argumentType = (tagNameIs("arguments", cur)) ? cur->xmlChildrenNode : cur;
      while (argumentType != NULL) {
        if (!tagNameIs("text", argumentType)) {
          printName(argumentType, outFile);
          if (tagNameIs("arguments", cur)) {
            fprintf(outFile, " ARGUMENTS");
          }
          fprintf(outFile, "</h1>");
          xmlNodePtr argumentElement;
          argumentElement = argumentType->xmlChildrenNode;
          int first = 1, justIntro = 0, justArg = 0;
          // fprintf (outFile,"\n<table border=\"3\">\n");
          while (1) {
            if (!tagNameIs("text", argumentElement)) {
              if ((first || justIntro) && !tagNameIs("intro", argumentElement))
                fprintf(outFile,
                        "\n<table border=\"3\">\n<tr><th>%s</th><th>"
                        "Explanation</th></tr>\n",
                        (tagNameIs("arguments", cur)) ? "Argument" : "Output");

              if (justArg) {
                if (tagNameIs("argument", argumentElement) || tagNameIs("output", argumentElement)) {
                  fprintf(outFile, "<td>&nbsp;</td>\n");
                }
                justArg = 0;
              }
              else if (tagNameIs("explanation", argumentElement)) {
                fprintf(outFile, "</tr>\n<tr>\n<td>&nbsp;</td>\n");
              }
              if (tagNameIs("argument", argumentElement) || tagNameIs("output", argumentElement)) {
                justArg = 1;
              }
              if (tagNameIs("intro", argumentElement)) {
                char *c, *intro;
                c = (char *)xmlNodeListGetString(doc, argumentElement->xmlChildrenNode, 1);
                if (first) {
                  intro = malloc(strlen("\n<h3>") + strlen(c) + 100);
                  strcpy(intro, "\n<h3>");
                  first = 0;
                }
                else {
                  intro = malloc(strlen("</tr>\n</table>\n<h3>") + strlen(c) + 100);
                  strcpy(intro, "</tr>\n</table>\n<h3>");
                }
                intro = strcat(intro, c);
                free(c);
                while (tagNameIs("intro", argumentElement->next) || tagNameIs("text", argumentElement->next)) {
                  if (!tagNameIs("text", argumentElement->next)) {
                    char *content;
                    content = (char *)xmlNodeListGetString(doc, argumentElement->next->xmlChildrenNode, 1);
                    char *secondHalf = malloc(strlen("<br />") + strlen(content) + 100);
                    strcpy(secondHalf, "<br />");
                    secondHalf = strcat(secondHalf, content);
                    char *temp = malloc(strlen(intro) + strlen(secondHalf) + 100);
                    strcpy(temp, intro);
                    temp = strcat(temp, secondHalf);
                    free(intro);
                    intro = temp;
                  }
                  argumentElement = argumentElement->next;
                }
                char *temp = malloc(strlen(intro) + strlen("</h3>\n") + 100);
                strcpy(temp, intro);
                temp = strcat(temp, "</h3>");
                // HACK this or the free below causes crash                free(intro);
                intro = temp;
                fprintf(outFile, "%s\n", intro);
                // HACK this causes crash, debug it:                free(intro);
                justIntro = 1;
              }
              else {
                if (first) {
                  first = 0;
                }
                else if ((tagNameIs("argument", argumentElement) || tagNameIs("output", argumentElement)) &&
                         !justIntro) {
                  fprintf(outFile, "</tr>\n");
                }

                justIntro = 0;
                printTableContents(doc, argumentElement, outFile);
              }
            }
            argumentElement = argumentElement->next;
            if (argumentElement == NULL) {
              if (!justIntro) {
                fprintf(outFile, "</tr>\n</table>\n");
              }
              break;
            }
          }
        }
        if (!tagNameIs("arguments", cur)) {
          break;
        }
        argumentType = argumentType->next;
      }
    }
    // Outputs all the other tags
    else if (!tagNameIs("text", cur)) {
      printName(cur, outFile);
      if (tagNameIs("EXAMPLE", cur)) {
        fprintf(outFile, " %d", exampleNum++);
      }
      fprintf(outFile, "</h1>");
      printContents(doc, cur, outFile);
    }
    cur = cur->next;
  }

  fprintf(outFile, "</html>\n</body>\n");

  return EXIT_SUCCESS;
}

// Changes a string to unpper case
void toUpperCase(char *c)
{
  int i;
  for (i = 0; *(c + i) != '\0'; i++) {
    *(c + i) = toupper(*(c + i));
  }
}

// Prints the name of a tag in the correct format
void printName(xmlNodePtr cur, FILE *f)
{
  xmlChar *n;
  n = (xmlChar *)cur->name;
  toUpperCase((char *)n);
  replaceUnderscore((char *)n);
  fprintf(f, "\n\n<h1>%s", n);
}

// Replaces the first underscore in a string with a space
void replaceUnderscore(char *c)
{
  char *pos = strchr(c, '_');
  if (pos != NULL) {
    *pos = ' ';
  }
}

// Prints the text of a tag in the corect format
void printContents(xmlDocPtr doc, xmlNodePtr cur, FILE *f)
{
  xmlChar *contents;
  contents = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
  unsigned int i;
  fprintf(f, "\n<p>");
  for (i = 0; i < strlen((char *)contents); i++) {
    if (*(contents + i) == '\n') {
      fprintf(f, "<br />");
    }
    else if (*(contents + i) == '<') {
      fprintf(f, "&lt;");
    }
    else if (*(contents + i) == '>') {
      fprintf(f, "&gt;");
    }
    else if (*(contents + i) == '&') {
      fprintf(f, "&amp;");
    }
    else if (*(contents + i) == '\"')
      fprintf(f, "\\\"");
    else {
      fprintf(f, "%c", *(contents + i));
    }
  }
  fprintf(f, "</p>");
  xmlFree(contents);
}

// Prints the text of a tag in the corect format
void printTableContents(xmlDocPtr doc, xmlNodePtr cur, FILE *f)
{
  xmlChar *contents;
  contents = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
  unsigned int i;
  if (tagNameIs("argument", cur) || tagNameIs("output", cur)) {
    fprintf(f, "<tr>\n");
  }
  if (tagNameIs("intro", cur)) {
    fprintf(f, "<h3><br />");
  }
  else {
    fprintf(f, "<td>");
  }
  for (i = 0; i < strlen((char *)contents); i++) {
    if (*(contents + i) == '\n') {
      fprintf(f, "<br />");
    }
    else if (*(contents + i) == '<') {
      fprintf(f, "&lt;");
    }
    else if (*(contents + i) == '>') {
      fprintf(f, "&gt;");
    }
    else if (*(contents + i) == '&') {
      fprintf(f, "&amp;");
    }
    else if (*(contents + i) == '\"')
      fprintf(f, "\\\"");
    else {
      fprintf(f, "%c", *(contents + i));
    }
  }
  if (!tagNameIs("intro", cur)) {
    fprintf(f, "</td>");
  }

  xmlFree(contents);
}

// Checks if the name of the tag is the same as the string parameter
int tagNameIs(char *c, xmlNodePtr cur)
{
  if (cur == NULL) {
    return 0;
  }
  return !(xmlStrcmp(cur->name, (const xmlChar *)c));
}

// Counts the number of chars until the next space, dash, underscore, slash,
// or end of the string
// Used to keep a word from running on two lines
int wrdLen(char *c)
{
  int i;
  for (i = 0; i < strlen(c); i++)
    if (*(c + i) == ' ' || *(c + i) == '_' || *(c + i) == '/') {
      return i;
    }
  return i;
}
