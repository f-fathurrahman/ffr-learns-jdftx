#include "my_jdftx.h"

// This is originally done in elecFluidMinimize before calling elecEnergyAndGrad
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
    e.eInfo.write(e.eVars.F, "eVars_F_step00.dat");
    e.eInfo.fillingsUpdate = ElecInfo::FillingsConst;
    e.eInfo.write(e.eVars.F, "eVars_F_step01.dat");
    //write_eVars_F(e);
    // grad and Kgrad are not computed.
    double Etot = my_elecEnergyAndGrad( e, 0, 0, true );
    logPrintf("Etot = %18.10f (can be ignored)\n", Etot);
    // Hsub is computed
    //logPrintf("Energy components:");
    //e.ener.print();
    // XXX: Etot is not used here, only Hsub

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
