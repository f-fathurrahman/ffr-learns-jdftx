void my_elecFluidMinimize( Everything& e )
{
  // Get variables
  //Control &cntrl = e.cntrl;
  ElecVars &eVars = e.eVars;
  ElecInfo& eInfo = e.eInfo;
  //IonInfo& iInfo = e.iInfo;
  Energies &ener = e.ener;

  // Haux will be initialized here
  // This section will not be executed for for FillingsConst (insulator case).

  if( !eVars.HauxInitialized &&
     eInfo.fillingsUpdate == ElecInfo::FillingsHsub )
  {

    //Constant nElectrons mode
    if(std::isnan(eInfo.mu))
    {   
      logPrintf("\nSetting the auxilliary hamiltonian equal to the subspace hamiltonian.\n");
      
      // calculate Hsub at current fillings:
      eInfo.fillingsUpdate = ElecInfo::FillingsConst;

      logPrintf("\nBefore calculating energy:\n");
      ener.print(); logPrintf("\n");

      eVars.elecEnergyAndGrad( e.ener, 0, 0, true );
      
      logPrintf("\nAfter calculating energy:\n");
      ener.print(); logPrintf("\n");

      eInfo.fillingsUpdate = ElecInfo::FillingsHsub;
      
      // Update B:
      for( int q=e.eInfo.qStart; q<e.eInfo.qStop; q++ )
      {
        eVars.Haux_eigs[q] = eVars.Hsub_eigs[q];
      }
      eVars.HauxInitialized = true;

      //print_Hsub_eigs( e );
    }
    // constant mu mode
    else
    {
      die("Constant chemical potential auxilliary hamiltonian fillings requires the band\n"
        "eigenvalues to either be read in using one of the two commands initial-state\n"
        "or elec-initial-eigenvals, or be automatically initialized during LCAO.");
    }

  }


  // This is only entered if we use fluid
  if( eVars.isRandom && eVars.fluidParams.fluidType!=FluidNone )
  {
    die("!!!! SHOULD NOT ENTER HERE !!!!\n");
  }

  //Prevent change in mu from abruptly changing electron count:
  if( eInfo.fillingsUpdate == ElecInfo::FillingsHsub &&
    !std::isnan(eInfo.mu) )
  { 
    double Bz, mu = eInfo.findMu( eVars.Haux_eigs, eInfo.nElectrons, Bz );
    
    logPrintf("Shifting auxilliary hamiltonian by %lf to set nElectrons=%lf\n", eInfo.mu-mu, eInfo.nElectrons);
    
    for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++) {
      eVars.Haux_eigs[q] += eye(eInfo.nBands)*(eInfo.mu-mu);
    }
  }
  std::cout << "fillingsUpdate = " << eInfo.fillingsUpdate << std::endl;
  std::cout << "mu = " << eInfo.mu << std::endl;

  double zz;

  zz = eVars.elecEnergyAndGrad(ener);
  logPrintf("zz = %18.10f\n", zz);
  ener.print();

  zz = eVars.elecEnergyAndGrad(ener);
  logPrintf("zz = %18.10f\n", zz);
  ener.print();

  logPrintf("\n-------- Electronic minimization -----------\n"); logFlush();
  
  //elecMinimize(e); driver: use scf, direct min, etc
  
  //ElecMinimizer emin(e);
  MyElecMinimizer emin(e);
  //emin.minimize(e.elecMinParams);
  //if (!e.ionDynamicsParams.tMax) {
  //  e.eVars.setEigenvectors(); //Don't spend time with this if running MD
  //  logPrintf("Pass here 93\n");
  //} 

  return;
}