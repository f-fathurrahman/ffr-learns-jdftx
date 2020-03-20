#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ElecVars.h>
#include <electronic/ElecMinimizer.h>
#include <core/Util.h>
#include <commands/parser.h>

#include <iostream>
#include <cmath>

#include "MyElecGradient.cpp"
#include "my_elecEnergyAndGrad.cpp"

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

    double Etot = my_elecEnergyAndGrad( e, &g, &Kg, false ); // ener is included in e

    logPrintf("Etot = %18.10f\n", Etot);

    printf("\n");
    printf("%s is finished\n", argv[0]);
    printf("\n");

    return 0;
}