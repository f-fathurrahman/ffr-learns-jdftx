#include "my_jdftx.h"

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