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

    ElecInfo& eInfo = e.eInfo;
    ElecVars& eVars = e.eVars;

    for(int q = e.eInfo.qStart; q < e.eInfo.qStop; q++)
    {
        const QuantumNumber& qnum = eInfo.qnums[q];
        diagMatrix& Fq = eVars.F[q];
        logPrintf("q = %d  weight = %f  tr(F) = %f\n", q, qnum.weight, trace(Fq));
        std::cout << Fq.data()[0] << std::endl;
        std::cout << Fq.data()[1] << std::endl;
        std::cout << Fq.data()[2] << std::endl;
        std::cout << Fq.data()[3] << std::endl;
    }

    return 0;

}