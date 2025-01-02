#include "my_jdftx.h"

double my_elecEnergyAndGrad( Everything& e, 
  MyElecGradient* grad, MyElecGradient* Kgrad, bool calc_Hsub )
{

  logPrintf("------------------------------------\n");
  logPrintf("**** ENTER my_elecEnergyAndGrad ****\n");
  logPrintf("------------------------------------\n");

  // Shortcuts
  const ElecInfo& eInfo = e.eInfo;
  ElecVars& eVars = e.eVars;
  Energies& ener = e.ener;

  //std::cout << "eInfo.mu = " << eInfo.mu << std::endl;

  // Cleanup old gradients:
  if( grad )
  {
    grad->C.assign(eInfo.nStates, ColumnBundle());
  }
  //
  if( Kgrad )
  {
    Kgrad->C.assign(eInfo.nStates, ColumnBundle());
  }

  // Determine whether Hsub and hence HC needs to be calculated:
  bool need_Hsub = calc_Hsub || grad;
  double mu = 0., Bz = 0.;

  //logPrintf("need_Hsub = %d\n", need_Hsub);

  if( eInfo.fillingsUpdate == ElecInfo::FillingsHsub )
  {    
    // Update nElectrons from mu, or mu from nElectrons as appropriate:
    if( std::isnan(eInfo.mu) )
    {
      //logPrintf("Finding mu: ");
      mu = eInfo.findMu( eVars.Haux_eigs, eInfo.nElectrons, Bz );
      //logPrintf(" mu = %18.10f\n", mu);
      //std::cout << "Bz = " << Bz << std::endl;
    }
    else
    {
      mu = eInfo.mu;
      ((ElecInfo&)eInfo).nElectrons = eInfo.nElectronsCalc( mu, eVars.Haux_eigs, Bz );
    }
    
    // Compute fillings from aux hamiltonian eigenvalues:
    for(int q=eInfo.qStart; q<eInfo.qStop; q++)
    {
      //logPrintf("q = %d\n", q);
      //logPrintf("Bz = %18.10f\n", Bz);
      //logPrintf("muEff = %18.10f\n", eInfo.muEff(mu, Bz, q));      
      eVars.F[q] = eInfo.smear( eInfo.muEff(mu, Bz, q), eVars.Haux_eigs[q] );
    }
    
    //logPrintf("\nBefore update fillings and energies:\n");
    //ener.print();

    // Update TS and muN:
    eInfo.updateFillingsEnergies( eVars.Haux_eigs, ener );

    //logPrintf("\nAfter update fillings and energies:\n");
    //ener.print();

    // Report for SCF (ElecMinimizer handles for minimize):
    if( e.cntrl.scf && eVars.n[0] ) eInfo.smearReport();
  }


  //
  // Update the density and density-dependent pieces if required:
  //
  eVars.n = eVars.calcDensity();
  
  if( e.exCorr.needsKEdensity() ) eVars.tau = eVars.KEdensity();
  
  // Atomic density matrix contributions for DFT+U (can be skipped for normal case)
  if( eInfo.hasU ) e.iInfo.rhoAtom_calc( eVars.F, eVars.C, eVars.rhoAtom);
  
  // Calculate density functional and its gradient
  eVars.EdensityAndVscloc( ener );

  // Update Vscloc projected onto spherical functions for ultrasoft psps
  if( need_Hsub ) e.iInfo.augmentDensityGridGrad( eVars.Vscloc );

  //logPrintf("\nAfter update density and density-dependent pieces:\n");
  //ener.print();


  //
  // update wavefunction dependent parts
  //

  //gradient w.r.t C (upto weights and fillings)
  std::vector<ColumnBundle> HC(eInfo.nStates);
  
  //Exact exchange if required:
  ener.E["EXX"] = 0.;
  //if( e.exCorr.exxFactor() )
  //{
  //  double aXX = e.exCorr.exxFactor();
  //  double omega = e.exCorr.exxRange();
  //  assert(e.exx);
  //  ener.E["EXX"] = (e.exx)(aXX, omega, eVars.F, eVars.C, need_Hsub ? &HC : 0); // not compiled ?
  //}
  
  //Do the single-particle contributions one state at a time to save memory (and for better cache warmth):
  ener.E["KE"] = 0.;
  ener.E["Enl"] = 0.;
  for(int q=eInfo.qStart; q < e.eInfo.qStop; q++)
  {
    double KEq = eVars.applyHamiltonian(q, eVars.F[q], HC[q], ener, need_Hsub);
    
    if(grad) //Calculate wavefunction gradients:
    {
      const QuantumNumber& qnum = eInfo.qnums[q];
      
      HC[q] -= O( eVars.C[q] ) * eVars.Hsub[q]; //Include orthonormality contribution
      
      grad->C[q] = HC[q] * ( eVars.F[q]*qnum.weight );
      
      if( Kgrad )
      {
        double Nq = qnum.weight*trace( eVars.F[q] );
        
        double KErollover = 2. * (Nq>1e-3 ? KEq/Nq : 1.);
        logPrintf("Nq = %f, KEq = %f, KErollover = %f\n", Nq, KEq, KErollover);
        
        precond_inv_kinetic( HC[q], KErollover ); //apply preconditioner
        
        std::swap( Kgrad->C[q], HC[q]); //this frees HC[q]
      }
    }
  }
  mpiWorld->allReduce(ener.E["KE"], MPIUtil::ReduceSum);
  mpiWorld->allReduce(ener.E["Enl"], MPIUtil::ReduceSum);

  //logPrintf("\nAfter calc KE and Enl:\n");
  //ener.print();

  double dmuContrib = 0.0;
  double dBzContrib = 0.0;

  // whether magnetization needs to be constrained
  bool Mconstrain = (eInfo.spinType==SpinZ) and std::isnan(eInfo.Bz);

  //logPrintf("Mconstrain  = %d\n", Mconstrain);
  
  // contribution due to N/M constraint via the mu/Bz gradient 
  if( grad and eInfo.fillingsUpdate==ElecInfo::FillingsHsub and (std::isnan(eInfo.mu) or Mconstrain) )
  {
    logPrintf("Pass here 170 in my_elecEnergyAndGrad\n");

    // numerator and denominator of dmuContrib resolved by spin channels (if any)
    double dmuNum[2] = {0.,0.}, dmuDen[2] = {0.,0.}; 
    
    logPrintf("Calculating dmuNum and dmuDen\n");
    double wsum = 0.0;
    for(int q=eInfo.qStart; q<eInfo.qStop; q++)
    {
      diagMatrix fprime = eInfo.smearPrime( eInfo.muEff(mu,Bz,q), eVars.Haux_eigs[q] );
      double w = eInfo.qnums[q].weight;
      int sIndex = eInfo.qnums[q].index();
      wsum = wsum + w;
      logPrintf("q = %d sIndex = %d w = %f\n", q, sIndex, w);
      dmuNum[sIndex] += w * trace(fprime * ( diag(eVars.Hsub[q]) - eVars.Haux_eigs[q]) );
      dmuDen[sIndex] += w * trace(fprime);
    }
    logPrintf("wsum = %f\n", wsum);
    
    mpiWorld->allReduce(dmuNum, 2, MPIUtil::ReduceSum);
    mpiWorld->allReduce(dmuDen, 2, MPIUtil::ReduceSum);
    
    if(std::isnan(eInfo.mu) and Mconstrain)
    {
      logPrintf("Pass here 189 in my_elecEnergyAndGrad\n");
      //Fixed N and M (effectively independent constraints on Nup and Ndn)
      double dmuContribUp = dmuNum[0]/dmuDen[0];
      double dmuContribDn = dmuNum[1]/dmuDen[1];
      dmuContrib = 0.5*(dmuContribUp + dmuContribDn);
      dBzContrib = 0.5*(dmuContribUp - dmuContribDn);
    }
    else if(Mconstrain)
    {
      logPrintf("Pass here 198 in my_elecEnergyAndGrad\n");
      //Fixed M only
      dmuContrib = 0.;
      dBzContrib = (dmuNum[0]-dmuNum[1])/(dmuDen[0]-dmuDen[1]);
    }
    else
    {
      logPrintf("Pass here 205 in my_elecEnergyAndGrad\n");
      //Fixed N only
      dmuContrib = (dmuNum[0]+dmuNum[1])/(dmuDen[0]+dmuDen[1]);
      dBzContrib = 0.;
    }
  }
  else {
    logPrintf("Pass here 209 in my_elecEnergyAndGrad.\n");
  }


  //Auxiliary hamiltonian gradient:
  if( grad && eInfo.fillingsUpdate==ElecInfo::FillingsHsub )
  {
    for(int q=eInfo.qStart; q < eInfo.qStop; q++)
    {

      logPrintf("q = %d\n", q);

      const QuantumNumber& qnum = eInfo.qnums[q];
      
      // gradient w.r.t fillings except for constraint contributions
      matrix gradF0 = eVars.Hsub[q] - eVars.Haux_eigs[q];
      
      // gradient w.r.t fillings
      matrix gradF = gradF0 - eye(eInfo.nBands)*eInfo.muEff(dmuContrib,dBzContrib,q);
      
      grad->Haux[q] = qnum.weight * dagger_symmetrize(eInfo.smearGrad(eInfo.muEff(mu,Bz,q), eVars.Haux_eigs[q], gradF));
      
      if( Kgrad ) //Drop the fermiPrime factors in preconditioned gradient:
        Kgrad->Haux[q] = (-e.cntrl.subspaceRotationFactor) * gradF0;
    }
  }

  logPrintf("eInfo.nBands = %d\n", eInfo.nBands);

  logPrintf("-----------------------------------\n");
  logPrintf("**** EXIT my_elecEnergyAndGrad ****\n");
  logPrintf("-----------------------------------\n");

  return relevantFreeEnergy(e);
}