#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ElecVars.h>
#include <electronic/ElecMinimizer.h>
#include <core/Util.h>
#include <commands/parser.h>

#include <iostream>
#include <cmath>

#include <string>

#include "my_jdftx.h"

using namespace std;

int main( int argc, char** argv ) {

  Everything e;
  InitParams ip("Performs Joint DFT calculations.", &e);
  initSystemCmdline(argc, argv, ip);
  parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
  e.setup();

  write_info(e);

  cout << "Basis size = " << e.basis.size() << endl;
  int Nkspin = e.basis.size();
  for(int ikspin=0; ikspin < Nkspin; ikspin++) {
    cout << "Basis: " << e.basis[ikspin].nbasis << endl;
    std::stringstream ss;
    ss << "iGarr_" << ikspin+1 << ".bindat";
    e.basis[ikspin].iGarr.write(ss.str().c_str());
  }
  //
  auto data0 = e.basis[0].iGarr.data();
  for(size_t i=0; i<e.basis[0].nbasis; i++) {
    cout << (i+1) << " " << data0->x() << " " << data0->y() << " "  << data0->z() << endl;
    data0++; // increment by using the size of whatever this pointer points to
  }
  //
  //cout << e.basis[0].head.size() << endl;



  return 0;
}
