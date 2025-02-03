#include "my_jdftx.h"

double my_applyHamiltonian(
  Everything& e,
  int q,
  const diagMatrix& Fq,
  ColumnBundle& HCq,
  Energies& ener,
  bool need_Hsub
){

  //make sure wavefunction is available for this states
  assert( e.eVars.C[q] );

  const QuantumNumber& qnum = e.eInfo.qnums[q];
  
  std::vector<matrix> HVdagCq(e.iInfo.species.size());
  
  //Propagate grad_n (Vscloc) to HCq (which is grad_Cq upto weights and fillings) if required
  if(need_Hsub) {
    //
    HCq += Idag_DiagV_I( e.eVars.C[q], e.eVars.Vscloc); //Accumulate Idag Diag(Vscloc) I C
    //
    e.iInfo.augmentDensitySphericalGrad(qnum, e.eVars.VdagC[q], HVdagCq); //Contribution via pseudopotential density augmentation
    //
    // Contribution via orbital KE:
    if(e.exCorr.needsKEdensity() && e.eVars.Vtau[qnum.index()]) {
      for(int iDir=0; iDir<3; iDir++) {
        HCq -= (0.5 * e.gInfo.dV) * D(
          Idag_DiagV_I(
            D(
              e.eVars.C[q], iDir
            ),
            e.eVars.Vtau),
            iDir
          );
      }
    }
    // Contribution via atomic density matrix projections (DFT+U)
    if(e.eInfo.hasU) {
      e.iInfo.rhoAtom_grad(e.eVars.C[q], e.eVars.U_rhoAtom, HCq);
    }
  }

  //Kinetic energy:
  double KEq;
  { // why?
    ColumnBundle LCq = L( e.eVars.C[q] );
    if(HCq) {
      HCq += (-0.5) * LCq;
    }
    KEq = qnum.weight * (-0.5) * traceinner(Fq, e.eVars.C[q], LCq).real();
    ener.E["KE"] += KEq;
  }
  
  //Nonlocal pseudopotentials:
  ener.E["Enl"] += qnum.weight * e.iInfo.EnlAndGrad(qnum, Fq, e.eVars.VdagC[q], HVdagCq);
  if(HCq) e.iInfo.projectGrad(HVdagCq, e.eVars.C[q], HCq);
  
  //Compute subspace hamiltonian if needed:
  if(need_Hsub) {
    e.eVars.Hsub[q] = e.eVars.C[q] ^ HCq;
    e.eVars.Hsub[q].diagonalize( e.eVars.Hsub_evecs[q], e.eVars.Hsub_eigs[q] );
  }
  return KEq;
}