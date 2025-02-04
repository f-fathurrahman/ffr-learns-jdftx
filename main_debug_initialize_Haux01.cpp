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

  logPrintf("---------------------\n");
  logPrintf("ENTER initialize_Haux\n");
  logPrintf("---------------------\n");

  // This is required to properly initilize ElecVars for metallic system.
  if( !e.eVars.HauxInitialized &&
      e.eInfo.fillingsUpdate==ElecInfo::FillingsHsub ) {

    logPrintf("Haux is not yet initialized\n");
    logPrintf("Setting up things for FillingsHsub\n");

    // Using FillingsConst
    e.eInfo.fillingsUpdate = ElecInfo::FillingsConst;
    // grad and Kgrad are not computed.
    double Etot = my_elecEnergyAndGrad( e, 0, 0, true );
    logPrintf("Etot = %18.10f\n", Etot);
    // Hsub is computed
    logPrintf("Energy components:");
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

  logPrintf("---------------------\n");
  logPrintf("EXIT initialize_Haux\n");
  logPrintf("---------------------\n");

}

int main( int argc, char** argv ) {

  Everything e;
  InitParams ip("Performs Joint DFT calculations.", &e);
  initSystemCmdline(argc, argv, ip);
  parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
  e.setup();

  write_info(e);

  initialize_Haux(e);

  e.eInfo.write(e.eVars.C, "eVars_C.dat");
  logPrintf("eVars.C is written to evars_C.dat\n");

  e.eInfo.write(e.eVars.Hsub, "eVars_Hsub.dat");
  logPrintf("eVars.Hsub is written to evars_Hsub.dat\n");

  return 0;
}

