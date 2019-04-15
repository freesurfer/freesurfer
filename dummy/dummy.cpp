#include <iostream>

#include "log.h"
#include "timer.h"

int main(int argc, const char **argv) {
  logDebug << currentDateTime();
}
