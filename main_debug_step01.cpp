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



  MyElecGradient g, gPrev, Kg;
  
  // get initial energy and gradient
  double E = my_ElecMinimizer_compute(e, elecMin, &g, &Kg);
  logPrintf("Initial energy E = %18.10f\n", E);
  e.ener.print();

  // step direction
  MyElecGradient d = clone(Kg);
  elecMin.constrain(d);

  double ET = 0.0;
  double alphaPrev = 0.0;

  // Try the test step:
  double alphaT = 1.0;
  my_ElecMinimizer_step(e, elecMin, d, alphaT-alphaPrev);
  alphaPrev = alphaT;
  ET = my_ElecMinimizer_compute(e, elecMin, 0, 0);
  logPrintf("ET = %le\n", ET);

  alphaT = 2.0; // total w.r.t to the first point
  // we will increment w.r.t previous point
  my_ElecMinimizer_step(e, elecMin, d, alphaT-alphaPrev);
  alphaPrev = alphaT;
  ET = my_ElecMinimizer_compute(e, elecMin, 0, 0);
  logPrintf("ET = %le\n", ET);

  alphaT = -2.0; // this will return to the first point
  my_ElecMinimizer_step(e, elecMin, d, alphaT);
  alphaPrev = alphaT;
  ET = my_ElecMinimizer_compute(e, elecMin, 0, 0);
  logPrintf("ET = %le\n", ET);


  return 0;
}

