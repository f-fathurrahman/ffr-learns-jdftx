#include <electronic/Everything.h>
#include <electronic/ColumnBundle.h>
#include <electronic/ElecVars.h>
#include <core/Util.h>
#include <commands/parser.h>

#include <iostream>
#include <cmath>

using namespace std;

// Dimensions of eVars and ColumnBundle
void debug_ColumnBundle_v01( Everything* e )
{
  unsigned long Nkspin = e->eVars.C.size();
  int Nstates = e->eVars.C[0].nCols();

  cout << "Nkspin  = " << Nkspin << endl;
  cout << "Nstates = " << Nstates << endl; // or nBands in JDFTX

  // Test formatted print
  logPrintf("%lu\n", e->eVars.C.size()); // unsigned long
  logPrintf("%d\n", e->eVars.C[0].nCols()); // int
  logPrintf("%lu\n", e->eVars.C[0].colLength()); // unsigned long

  // Print number of plane waves for each kpoint
  logPrintf("Number of plane waves for each kpoint:\n");
  for(unsigned long i = 0; i < Nkspin; i++) {
    logPrintf("%4lu %4d %4lu\n", i+1, e->eVars.C[i].nCols(), e->eVars.C[i].colLength());
  }

  return;
}

void test_print_data_01(Everything* e)
{
  // data
  int ikspin = 0;
  complex* dt = e->eVars.C[ikspin].data();
  int Nstates = e->eVars.C[ikspin].nCols();
  int Nbasis = e->eVars.C[ikspin].colLength();
  //
  logPrintf("ikspin = %d\n", ikspin);
  logPrintf("Nstates = %d\n", Nstates);
  logPrintf("Nbasis = %d\n", Nbasis);
  //
  for(unsigned int i=0; i < e->eVars.C[ikspin].colLength(); i++) {
    logPrintf("%5d %18.10f %18.10f im\n", i+1, dt[i].x, dt[i].y);
  }
  // access "last" data
  size_t idx_data = e->eVars.C[ikspin].index(Nstates, 0);
  cout << "N data total = " << Nstates*Nbasis << endl;
  cout << "idx_data = " << idx_data << endl;
  cout << dt[idx_data].x << endl;
  cout << dt[idx_data].y << endl;
}

void debug_ColumnBundle_v02( Everything* e )
{
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
  return;
}


void write_info(Everything &e)
{
  MPIUtil::File fp;
  mpiWorld->fopenWrite(fp, "INFO.dat");

  fprintf(fp, "detR (Cell Volume) = %18.10f\n", e.gInfo.detR);

  auto qStart = e.eInfo.qStart;
  auto qStop = e.eInfo.qStop;

  fprintf(fp, "qStart = %d\n", qStart);
  fprintf(fp, "qStop  = %d\n", qStop);

  auto Nkspin = e.eVars.C.size();
  auto Nstates = e.eVars.C[0].nCols();

  fprintf(fp, "Nkspin = %ld\n", Nkspin);
  fprintf(fp, "Nstates = %d\n", Nstates);

  for(unsigned long i = 0; i < Nkspin; i++) {
    fprintf(fp, "%ld \n", e.eVars.C[i].colLength());
  }
  mpiWorld->fclose(fp);
}

int main( int argc, char** argv )
{
  Everything e;
  InitParams ip("Performs Joint DFT calculations.", &e);
  initSystemCmdline(argc, argv, ip);

  parse(readInputFile(ip.inputFilename), e, ip.printDefaults);
  e.setup();
  
  //debug_ColumnBundle_v01(&e);
  //test_print_data_01(&e);

  write_info(e);

  e.eInfo.write(e.eVars.C, "eVars_C.dat");
  logPrintf("eVars.C is written to evars_C.dat\n");

  return 0;
}