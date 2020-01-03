#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ElecVars.h>
#include <electronic/ElecMinimizer.h>
#include <core/Util.h>
#include <commands/parser.h>

#include <iostream>
#include <cmath>

using namespace std;

#include "my_elecEnergyAndGrad.cpp"

int main( int argc, char** argv )
{
    Everything e;
    InitParams ip("Performs Joint DFT calculations.", &e);
    initSystemCmdline(argc, argv, ip);
    parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
    e.setup();

    my_elecEnergyAndGrad( e, e.ener, 0, 0, false );

    printf("\n");
    printf("%s is finished\n", argv[0]);
    printf("\n");

    return 0;
}