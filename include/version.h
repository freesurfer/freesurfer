#pragma once

#include <string>


std::string getVersion();
std::string getBuildStamp();
std::string getAllInfo(int argc, char **argv, const std::string& progname);

int handleVersionOption(int argc, char** argv, const char* progname);

char *argv2cmdline(int argc, char *argv[]);
char *VERuser(void);
char *VERfileTimeStamp(char *fname);
char *VERcurTimeStamp(void);
