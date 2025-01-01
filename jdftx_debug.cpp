#include "my_jdftx.h"

#include <iostream>
#include <cmath>

using namespace std;

void debug_dot( Everything* e )
{
    ColumnBundle c0 = e->eVars.C[0].getSub(0,1);
    ColumnBundle c1 = e->eVars.C[0].getSub(1,2);

    double detR = c0.basis->gInfo->detR;
    double d = dot(c0,c0);

    logPrintf("d = %18.10f\n", d*detR/2);

    return;
}

int main( int argc, char** argv )
{
    Everything e;
    InitParams ip("Performs Joint DFT calculations.", &e);
    initSystemCmdline(argc, argv, ip);

    parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
    e.setup();

    //debug_dot(&e);
    debug_elecEnergyAndGrad(e);

    return 0;
}