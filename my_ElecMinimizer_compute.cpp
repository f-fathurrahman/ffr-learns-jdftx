#include "my_jdftx.h"

double my_ElecMinimizer_compute(
  Everything& e,
  MyElecMinimizer& elecMin,
  MyElecGradient* grad,
  MyElecGradient* Kgrad )
{
  //
  if(grad) {
    grad->init(e);
  }
  if(Kgrad) {
    Kgrad->init(e);
  }
    
  double ener = my_elecEnergyAndGrad(e, grad, Kgrad, false);

  // some info if rotations are applied
  if(elecMin.rotExists) {
    logPrintf("my_ElecMinimizer_compute: rotations will be applied\n");
  } else {
    logPrintf("my_ElecMinimizer_compute: rotations will NOT be applied\n");
  }

  // onl
  if(grad) {
    for(int q = e.eInfo.qStart; q < e.eInfo.qStop; q++)
    {
      //Rotate wavefunction gradients if necessary:
      if(elecMin.rotExists) {
        grad->C[q] = grad->C[q] * elecMin.rotPrevCinv[q];
        Kgrad->C[q] = Kgrad->C[q] * elecMin.rotPrevCinv[q];
      }

      //Subspace gradient handling depends on mode:
      if( e.eInfo.fillingsUpdate == ElecInfo::FillingsHsub ) {
        //
        //Haux fillings: rotate gradient computed by ElecVars if necessary
        if(elecMin.rotExists) {
          grad->Haux[q] = elecMin.rotPrev[q] * grad->Haux[q] * dagger(elecMin.rotPrev[q]);
          Kgrad->Haux[q] = elecMin.rotPrev[q] * Kgrad->Haux[q] * dagger(elecMin.rotPrev[q]);
        }
      }
      else if(!e.eInfo.scalarFillings) {
        //Non-scalar fillings:
        grad->Haux[q] = dagger_symmetrize(complex(0,1) * ( e.eVars.F[q]*e.eVars.Hsub[q] - e.eVars.Hsub[q]*e.eVars.F[q]));
        // dagger_symmetrize is defined in core/matrixOperators.cpp
        Kgrad->Haux[q] = e.cntrl.subspaceRotationFactor * grad->Haux[q];
      }
      //else: constant scalar fillings (no subspace gradient)
    }
        
    //Cache gradient overlaps, if needed, for subspace rotation handling:
    if(elecMin.sra) {
      elecMin.sra->cacheGradientOverlaps(*grad, *Kgrad);
    }
    elecMin.KgradHaux = Kgrad->Haux;
  }
  return ener;
}