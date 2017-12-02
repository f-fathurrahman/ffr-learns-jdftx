#ifndef JDFTX_ELECTRONIC_SCFPARAMS_H
#define JDFTX_ELECTRONIC_SCFPARAMS_H

#include <core/Util.h>
#include <core/PulayParams.h>

//! @addtogroup ElecSystem
//! @{

//! Parameters controlling SCF iteration
struct SCFparams : public PulayParams
{
  int nEigSteps; //!< number of steps of the eigenvalue solver per iteration (use elecMinParams.nIterations if 0)
  double eigDiffThreshold; //!< convergence threshold on the RMS change of eigenvalues

  string historyFilename; //!< Read SCF history in order to resume a previous run
  
  enum MixedVariable
  {  MV_Density, //!< Mix electron density (n) and kinetic energy density (tau)
    MV_Potential //!< Mix the local electronic potential (Vscloc) and the kinetic energy potential (Vtau)
  }
  mixedVariable; //!< Whether we are mixing the density or the potential
  
  double qKerker; //!< Wavevector controlling Kerker preconditioning
  double qKappa; //!< wavevector controlling long-range damping (if negative, auto-set to zero or fluid Debye wave-vector as appropriate)
  
  bool verbose; //!< Whether the inner eigensolver will print progress
  double mixFractionMag;  //!< Mixing fraction for magnetization density / potential
  
  SCFparams()
  {  nEigSteps = 2; //for Davidson; the default for CG is 40 (and set by the command)
    eigDiffThreshold = 1e-8;
    mixedVariable = MV_Density;
    qKerker = 0.8;
    qKappa = -1.;
    verbose = false;
    mixFractionMag = 1.5;
  }
};

//! @}
#endif // JDFTX_ELECTRONIC_SCFPARAMS_H
