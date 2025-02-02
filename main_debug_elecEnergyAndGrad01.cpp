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

// formatted output
void export_Hsub(Everything& e) {
  for(int q = e.eInfo.qStart; q < e.eInfo.qStop; q++) {
    stringstream ss;
    ss << "Hsub_" << q << ".dat";
    FILE *fptr;
    fptr = fopen( ss.str().c_str(), "w");
    e.eVars.Hsub[q].print(fptr);
    fclose(fptr);
    cout << "Hsub written to " << ss.str() << endl;
  }
}


int main( int argc, char** argv )
{
  Everything e;
  InitParams ip("Performs Joint DFT calculations.", &e);
  initSystemCmdline(argc, argv, ip);
  parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
  e.setup();

  write_info(e);

  MyElecGradient g, Kg;
  // these will allocate memory?
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
    logPrintf("Initial values of Haux_eigs\n");
    for(int q = e.eInfo.qStart; q < e.eInfo.qStop; q++) {
      e.eVars.Haux_eigs[q] = e.eVars.Hsub_eigs[q];
      logPrintf("\nikspin index = %d\n", q);
      for(int ist=0; ist < e.eInfo.nStates; ist++) {
        logPrintf("%4d %18.10f\n",ist,  e.eVars.Haux_eigs[q][ist]);
      }
    }
    e.eVars.HauxInitialized = true;
  }

  e.eInfo.write(e.eVars.C, "eVars_C.dat");
  logPrintf("eVars.C is written to evars_C.dat\n");

  e.eInfo.write(e.eVars.Hsub, "eVars_Hsub.dat");
  logPrintf("eVars.Hsub is written to evars_Hsub.dat\n");

  // Not working
  /*
  std::vector<matrix> n;
  n.push_back(e.eVars.n[0]->toMatrix());
  e.eInfo.write(n, "eVars_n.dat");
  */

   matrix n = e.eVars.n[0]->toMatrix();
   cout << "n.nRows() = " << n.nRows() <<  " n.nCols() = " << n.nCols() << endl;
   cout << "Some n = " << n(0,0).x << " " << n(0,0).y << endl;
   cout << "Some n = " << n(1,0).x << " " << n(1,0).y << endl;

  //cout << "eVars.n = " << e.eVars.n[0](1) << endl;

  /*
  e.eInfo.write(e.eVars.n[0].to_matrix(), "eVars_n.dat");
  logPrintf("eVars.n is written to evars_n.dat\n");
  */

  //export_Hsub(e);

  double Etot = my_elecEnergyAndGrad( e, &g, &Kg, false ); // ener is included in e
  logPrintf("Etot = %18.10f\n", Etot);

  printf("\n");
  printf("%s is finished\n", argv[0]);
  printf("\n");

  return 0;
}