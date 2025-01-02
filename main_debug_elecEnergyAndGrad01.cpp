#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ElecVars.h>
#include <electronic/ElecMinimizer.h>
#include <core/Util.h>
#include <commands/parser.h>

#include <iostream>
#include <cmath>

#include "my_jdftx.h"

using namespace std;

int main( int argc, char** argv )
{
  Everything e;
  InitParams ip("Performs Joint DFT calculations.", &e);
  initSystemCmdline(argc, argv, ip);
  parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
  e.setup();

  MyElecGradient g, Kg;
  
  g.init(e);
  Kg.init(e);

  // This is required to properly initilize ElecVars for metallic system.
  if(!e.eVars.HauxInitialized && e.eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
  {
    cout << "e.eVars.HauxInitialized = " << e.eVars.HauxInitialized << endl;
    cout << "e.eInfo.fillingsUpdate = " << e.eInfo.fillingsUpdate << endl;
    logPrintf("Haux is not yet initialized\n");
    logPrintf("Setting up things for FillingsHsub\n");

    // Using FillingsConst
    e.eInfo.fillingsUpdate = ElecInfo::FillingsConst;
    
    // call the original member function of ElecVars to calculate energy and Hsub
    e.eVars.elecEnergyAndGrad( e.ener, 0, 0, true );

    e.eInfo.fillingsUpdate = ElecInfo::FillingsHsub;
    
    // Update B:
    // Use Hsub for initial value of Haux
    for(int q = e.eInfo.qStart; q < e.eInfo.qStop; q++) {
      e.eVars.Haux_eigs[q] = e.eVars.Hsub_eigs[q];
      cout << e.eVars.Haux_eigs[q][0] << endl;
    }
    e.eVars.HauxInitialized = true;
  }

  double Etot = my_elecEnergyAndGrad( e, &g, &Kg, false ); // ener is included in e

  logPrintf("Etot = %18.10f\n", Etot);

  printf("\n");
  printf("%s is finished\n", argv[0]);
  printf("\n");

  return 0;
}