#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ElecVars.h>
#include <core/Util.h>
#include <commands/parser.h>

#include <iostream>
#include <cmath>

using namespace std;

// Dimensions of eVars and ColumnBundle
void debug_ColumnBundle_v01( Everything* e )
{
    unsigned long Nkspin = e->eVars.C.size();
    int Nstates = e->eVars.C[0].nCols();

    cout << "Nkspin  = " << Nkspin << endl;
    cout << "Nstates = " << Nstates << endl; // or nBands in JDFTX

    logPrintf("%lu\n", e->eVars.C.size()); // unsigned long
    logPrintf("%d\n", e->eVars.C[0].nCols()); // int
    logPrintf("%lu\n", e->eVars.C[0].colLength()); // unsigned long

    for(unsigned long i = 0; i < Nkspin; i++) {
        logPrintf("%4lu %4d %4lu\n", i+1, e->eVars.C[i].nCols(), e->eVars.C[i].colLength());
    }

    return;

}

void debug_ColumnBundle_v02( Everything* e )
{
    //logPrintf("eVars.C.size()     = %d\n", e->eVars.C.size());
    cout << "eVars.C.size()     = " << e->eVars.C.size() << endl; // Nkspin
    logPrintf("eVars.C[0].nCols() = %d\n", e->eVars.C[0].nCols()); // Nstates
    cout << "eVars.C[0].colLength = " << e->eVars.C[0].colLength() << endl; // Ngw

    // data
    //complex* dt = e->eVars.C[0].data();
    //for(unsigned int i=0; i < e->eVars.C[0].colLength(); i++) {
    //    logPrintf("%5d %18.10f %18.10f im\n", i+1, dt[i].x, dt[i].y);
    //}

    ColumnBundle Oc;
    Oc = O(e->eVars.C[0]);
    //dt = Oc.data();
    //for(unsigned int i=0; i < Oc.colLength(); i++) {
    //    logPrintf("%5d %18.10f %18.10f im\n", i+1, dt[i].x, dt[i].y);
    //}
    return;
}


int main( int argc, char** argv )
{
    Everything e;
    InitParams ip("Performs Joint DFT calculations.", &e);
    initSystemCmdline(argc, argv, ip);

    parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
    e.setup();
    debug_ColumnBundle_v01(&e);

    return 0;
}