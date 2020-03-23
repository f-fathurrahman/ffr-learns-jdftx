#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ElecVars.h>
#include <core/Util.h>
#include <commands/parser.h>

#include <iostream>
#include <cmath>

using namespace std;


void debug_operator_O( Everything* e )
{
    size_t Ngw = e->eVars.C[0].colLength();
    ColumnBundle c0 = e->eVars.C[0].getSub(0,1);
    ColumnBundle c1 = O(c0);

    double detR = c0.basis->gInfo->detR;

    for(size_t i = 0; i < Ngw; i++) {
        logPrintf("%18.10f %18.10f %18.10f\n", c0.data()[i].x, c1.data()[i].x, c1.data()[i].x/detR);
    }

    return;
}


void debug_operators( Everything* e )
{
    ColumnBundle Oc;
    Oc = O(e->eVars.C[0]);

    double detR = Oc.basis->gInfo->detR;
    cout << "detR   = " << detR << endl;
    cout << "1/detR = " << 1.0/detR << endl;
    cout << "1/sqrt(detR) = " << 1.0/sqrt(detR) << endl;

    //e->eVars.orthonormalize(0);

    // Test dot between different column of C[0]
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

    c1 = e->eVars.C[0];
    cout << "c1.nCols() = " << c1.nCols() << endl;
    cout << "dot C0, C0 = " << dot(c1, c1)/c1.nCols() << endl;
    cout << "dot C0, C0 = " << dot(c1, c1) << endl;

    complex* data = c1.data();
    double ss = 0.0;
    for(size_t i = 0; i < c1.colLength(); i++) {
        size_t idx = c1.index(0,i);
        ss = ss + data[idx].x * data[idx].x + data[idx].y * data[idx].y;
        //logPrintf("dat = [%le,%le]\n", data[idx].x, data[idx].y);
    }
    cout << "ss = " << 2*ss << endl;

    return;
}

int main( int argc, char** argv )
{
    Everything e;
    InitParams ip("Performs Joint DFT calculations.", &e);
    initSystemCmdline(argc, argv, ip);

    parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
    e.setup();

    //debug_operator_O(&e);
    
    debug_operators(&e);

    return 0;
}