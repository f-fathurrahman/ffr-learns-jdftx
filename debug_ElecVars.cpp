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

    ElecVars& eVars = e.eVars;
    size_t Nkspin = eVars.C.size();
    std::cout << "Nkspin = " << eVars.C.size() << std::endl;

    //
    // Try to modify ColumnBundle data: not working
    //
    /*for(int q = e.eInfo.qStart; q < e.eInfo.qStop; q++) {
        complex* data = eVars.C[q].data();
        int Nstates = eVars.C[q].nCols();
        int Ngwk = eVars.C[q].colLength();
        //std::cout << Ngwk << std::endl;
        logPrintf("Ngwk = %d, Nstates = %d\n", Ngwk, Nstates);
        for(size_t i=0; i < Nstates*Ngwk; i++) {
            data[i].x = 1.0;
            data[i].y = 0.0;
        }
        eVars.orthonormalize(q);
        //cout << eVars.data()
    }*/

    MyElecGradient g, Kg;
    
    g.init(e);
    Kg.init(e);

    double Etot = my_elecEnergyAndGrad( e, &g, &Kg, false );
    logPrintf("Etot = %18.10f\n", Etot);

    logPrintf("\nTest dot:\n");
    for(size_t i = 0; i < Nkspin; i++)
    {
        std::cout << dot(eVars.C[i], O(eVars.C[i])) << std::endl;
    }

    logPrintf("\nTest dotc:\n");
    for(size_t i = 0; i < Nkspin; i++)
    {
        std::cout << dotc(eVars.C[i], O(eVars.C[i])).real() << std::endl;
    }

    logPrintf("\nTest dot eVars.C and eVars.C:\n");
    double ss = 0.0;
    for(size_t i = 0; i < Nkspin; i++)
    {
        ss = ss + dot(eVars.C[i], eVars.C[i]);
        std::cout << dot(eVars.C[i], eVars.C[i]) << std::endl;
    }
    std::cout << "ss = " << ss << std::endl;


    logPrintf("\nTest dot g and g:\n");
    ss = 0.0;
    for(size_t i = 0; i < Nkspin; i++)
    {
        ss = ss + dot(g.C[i], O(g.C[i]));
        std::cout << dot(g.C[i], O(g.C[i])) << std::endl;
    }
    std::cout << "ss = " << ss << std::endl;
    std::cout << "All: g,g  " << dot(g, g) << std::endl;
    std::cout << "All: g,Kg " << dot(g, Kg) << std::endl;

    matrix A = eVars.C[0]^O(eVars.C[0]);
    std::cout << A.nRows() << std::endl;
    std::cout << A.nCols() << std::endl;
    complex* dA = A.data();
    std::cout << dA[A.index(0,0)].x << std::endl;
    std::cout << dA[A.index(1,0)].x << std::endl;
    std::cout << dA[A.index(1,1)].x << std::endl;
    std::cout << dA[A.index(2,2)].x << std::endl;
    std::cout << dA[A.index(2,1)].x << std::endl;
    std::cout << dA[A.index(3,3)].x << std::endl;
    //ColumnBundle& C0 = g.C[0];

    //logPrintf("dot g Kg = %18.10f\n", dot(g,Kg));
    
    //int nCols = g.C[0].nCols();
    //logPrintf("nCols = %d\n", nCols);
    //std::cout << "Nkspin = " << g.C.size() << std::endl;
    //logPrintf("Nkspin = %d\n", g.C.size());
    //logPrintf("dot g0 Kg0 = %18.10f\n", dot(g.C[0], Kg.C[0]));
    
    //complex res = dotc(g.C[0], Kg.C[0]);
    //logPrintf("dotc g0 Kg0 = %18.10f %18.10f\n", res.x, res.y);
    //logPrintf("dot g1 Kg1 = %18.10f\n", dot(g.C[1], Kg.C[1]));

    //logPrintf("Pass here in main:27\n");

    return 0;

}