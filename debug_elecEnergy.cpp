#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ElecVars.h>
#include <core/Util.h>
#include <commands/parser.h>

#include <iostream>
#include <cmath>

using namespace std;

void debug_elecEnergyAndGrad( Everything& e )
{

    // This is required to properly initilizing ElecVars for metallic system.
    if(!e.eVars.HauxInitialized && e.eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
    {   
        e.eInfo.fillingsUpdate = ElecInfo::FillingsConst;
        e.eVars.elecEnergyAndGrad(e.ener, 0, 0, true);
        e.eInfo.fillingsUpdate = ElecInfo::FillingsHsub;
        //Update B:
        for(int q = e.eInfo.qStart; q < e.eInfo.qStop; q++) e.eVars.Haux_eigs[q] = e.eVars.Hsub_eigs[q];
        e.eVars.HauxInitialized = true;
    }
    e.eVars.elecEnergyAndGrad( e.ener );

    logPrintf("KE     = %18.10f\n", e.ener.E["KE"]);
    logPrintf("Eloc   = %18.10f\n", e.ener.E["Eloc"]);
    logPrintf("EH     = %18.10f\n", e.ener.E["EH"]);
    logPrintf("Exc    = %18.10f\n", e.ener.E["Exc"]);
    logPrintf("Enl    = %18.10f\n", e.ener.E["Enl"]);    
    logPrintf("Epulay = %18.10f\n", e.ener.E["Epulay"]);
    logPrintf("Ewald  = %18.10f\n", e.ener.E["Eewald"]);
    //logPrintf("Total  = %18.10f\n", e.ener);
    cout << e.ener.E << endl;    
    
    return;
}

int main( int argc, char** argv )
{
    Everything e;
    InitParams ip("Performs Joint DFT calculations.", &e);
    initSystemCmdline(argc, argv, ip);

    parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
    e.setup();

    debug_elecEnergyAndGrad(e);

    return 0;
}