#include "my_jdftx.h"

ScalarFieldArray my_calcDensity(Everything& e)
{

  ScalarFieldArray density( e.eVars.n.size() );
  
  // Runs over all states and accumulates density to the corresponding spin channel
  // of the total density
  e.iInfo.augmentDensityInit(); // XXX what's this?
  for(int q=e.eInfo.qStart; q < e.eInfo.qStop; q++)
  {
    // For standard case
    density += e.eInfo.qnums[q].weight * diagouterI(
      e.eVars.F[q], e.eVars.C[q], density.size(), &e.gInfo
    );
    // additional contributions: pseudopotential contribution
    e.iInfo.augmentDensitySpherical(
      e.eInfo.qnums[q], e.eVars.F[q], e.eVars.VdagC[q]
    );

  }
  //
  e.iInfo.augmentDensityGrid(density);
  //
  for(ScalarField& ns: density)
  {
    nullToZero(ns, e.gInfo);
    ns->allReduceData(mpiWorld, MPIUtil::ReduceSum);
  }
  e.symm.symmetrize(density);
  
  return density;
}

