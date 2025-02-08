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

  // Write iGarr (Miller indices?)
  cout << "Basis size = " << e.basis.size() << endl;
  int Nkspin = e.basis.size();
  for(int ikspin=0; ikspin < Nkspin; ikspin++) {
    cout << "Basis: " << e.basis[ikspin].nbasis << endl;
    std::stringstream ss;
    ss << "iGarr_" << ikspin+1 << ".bindat";
    e.basis[ikspin].iGarr.write(ss.str().c_str());
  }

  e.eInfo.write(e.eVars.F, "eVars_F.bindat");
  e.eInfo.write(e.eVars.C, "eVars_C.bindat");

  e.eVars.n = my_calcDensity(e);
  write_eVars_n(e);

  return 0;
}
