#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <core/Util.h>
#include <commands/parser.h>

#include <iostream>
#include <cmath>

using namespace std;

// important: use pass by reference
void debug_Everything(Everything* e)
{
    logPrintf("cntrl.ecut = %f\n", e->cntrl.Ecut);

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

    double detR = Oc.basis->gInfo->detR;
    cout << "detR   = " << endl;
    cout << "1/detR = " << 1.0/detR << endl;
    cout << "1/sqrt(detR) = " << 1.0/sqrt(detR) << endl;

    // Test dot
    ColumnBundle c1, c2;
    c1 = e->eVars.C[0].getSub(1,2);
    c2 = e->eVars.C[0].getSub(2,3);
    cout << "c1.colLength() = " << c1.colLength() << endl;
    cout << "c1 nCols()     = " << c1.nCols() << endl;
    cout << "dot(c1,c1)      = " << dot(c1, c1) << endl;
    cout << "dot(c1,c1)*detR = " << dot(c1, c1)*detR << endl;    
    cout << "dot(c1,c2)      = " << dot(c1, c2) << endl;    
    
    cout << "dot evars.C[1] = " << dot(e->eVars.C[1], e->eVars.C[1])/e->eVars.C[1].nCols() << endl;
    cout << "dot Oc         = " << dot(Oc, Oc) << endl;
}

int main( int argc, char** argv )
{
    Everything e;
    InitParams ip("Performs Joint DFT calculations.", &e);
    initSystemCmdline(argc, argv, ip);

    parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
    e.setup();
    debug_Everything(&e);

    return 0;
}