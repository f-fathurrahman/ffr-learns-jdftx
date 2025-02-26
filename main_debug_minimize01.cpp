#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ElecVars.h>
#include <electronic/ElecMinimizer.h>
#include <core/Util.h>
#include <commands/parser.h>

#include <iostream>
#include <cmath>

#include <string>

#include "my_jdftx.h"

using namespace std;

int main( int argc, char** argv ) {

  Everything e;
  InitParams ip("Performs Joint DFT calculations.", &e);
  initSystemCmdline(argc, argv, ip);
  parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
  e.setup();

  write_info(e);
  write_iGarr(e);

  initialize_Haux(e);

  e.eInfo.write(e.eVars.Hsub, "eVars_Hsub.bindat");
  e.eInfo.write(e.eVars.F, "eVars_F.bindat");
  e.eInfo.write(e.eVars.C, "eVars_C.bindat");

  MyElecMinimizer elecMin(e);
  //my_ElecMinimizer_minimize(e, elecMin, e.elecMinParams);
  e.elecMinParams.linePrefix = "SIMPLE_MIN: ";
  simple_minimize(e, elecMin, e.elecMinParams);

  e.ener.print();

  return 0;
}

