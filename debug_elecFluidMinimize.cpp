#include <electronic/Everything.h>
#include <electronic/ElecMinimizer.h>
#include <electronic/ColumnBundle.h>
#include <commands/parser.h>

#include "my_elecFluidMinimize.cpp"

int main( int argc, char** argv )
{
    Everything e;
    InitParams ip("Performs Joint DFT calculations.", &e);
    initSystemCmdline( argc, argv, ip );
    parse( readInputFile(ip.inputFilename), e, ip.printDefaults );
    e.setup();

    // Call this for the original implementation
    // elecFluidMinimize( e );

    my_elecFluidMinimize( e );

    return 0;
}
