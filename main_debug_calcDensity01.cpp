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

void initialize_Haux(Everything& e) {
  // This is required to properly initilize ElecVars for metallic system.
  if( !e.eVars.HauxInitialized &&
      e.eInfo.fillingsUpdate==ElecInfo::FillingsHsub ) {

    logPrintf("Haux is not yet initialized\n");
    logPrintf("Setting up things for FillingsHsub\n");

    // Using FillingsConst
    e.eInfo.fillingsUpdate = ElecInfo::FillingsConst;
    // XXX need to export this to file
    //
    // call the original member function of ElecVars to calculate energy and Hsub
    e.eVars.elecEnergyAndGrad( e.ener, 0, 0, true );
    // grad and Kgrad are not computed.
    // Hsub is computed
    e.ener.print();

    // Set again to FillingsHsub
    e.eInfo.fillingsUpdate = ElecInfo::FillingsHsub;
    
    // Update B:
    // Use Hsub for initial value of Haux
    logPrintf("Initial values of Haux_eigs\n");
    for(int q = e.eInfo.qStart; q < e.eInfo.qStop; q++) {
      // Set to Hsub_eig
      e.eVars.Haux_eigs[q] = e.eVars.Hsub_eigs[q];
    }
    e.eVars.HauxInitialized = true;
  }


}

int main( int argc, char** argv ) {

  Everything e;
  InitParams ip("Performs Joint DFT calculations.", &e);
  initSystemCmdline(argc, argv, ip);
  parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
  e.setup();

  write_info(e);

  initialize_Haux(e);

  // Using FillingsConst
  //e.eInfo.fillingsUpdate = ElecInfo::FillingsConst;

  e.eInfo.write(e.eVars.F, "eVars_F.bindat");
  e.eInfo.write(e.eVars.C, "eVars_C.bindat");
  e.eInfo.write(e.eVars.Hsub, "eVars_Hsub.bindat");

  e.eVars.n = my_calcDensity(e);
  write_eVars_n(e);

  return 0;
}
