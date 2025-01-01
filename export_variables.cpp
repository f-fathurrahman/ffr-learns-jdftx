#include "my_jdftx.h"

void write_info(Everything &e)
{
  auto qStart = e.eInfo.qStart;
  auto qStop = e.eInfo.qStop;

  MPIUtil::File fp;
  mpiWorld->fopenWrite(fp, "INFO.dat");
  fprintf(fp, "qStart = %d\n", qStart);
  fprintf(fp, "qStop  = %d\n", qStop);
  mpiWorld->fclose(fp);
}

void export_eVars_C(Everything &e)
{
  auto qStart = e.eInfo.qStart;
  auto qStop = e.eInfo.qStop;
  MPIUtil::File fp;
  mpiWorld->fopenWrite(fp, "eVars_C.dat");
	for(int q=qStart; q<qStop; q++) {
		mpiWorld->fwriteData(e.eVars.C[q], fp);
  }
	mpiWorld->fclose(fp);

  return;
}