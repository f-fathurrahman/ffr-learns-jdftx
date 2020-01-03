void my_elecEnergyAndGrad( Everything& e, 
    Energies& ener, ElecGradient* grad, ElecGradient* Kgrad, bool calc_Hsub )
{

    // This is required to properly initilizing ElecVars for metallic system.
    if(!e.eVars.HauxInitialized && e.eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
    {   
        logPrintf("\nSetting up things for FillingsHsub\n");

        // Using FillingsConst
        e.eInfo.fillingsUpdate = ElecInfo::FillingsConst;
        // call the original member function of ElecVars to calculate energy and Hsub
        e.eVars.elecEnergyAndGrad(e.ener, 0, 0, true);

        e.eInfo.fillingsUpdate = ElecInfo::FillingsHsub;
        
        // Update B:
        // Use Hsub for initial value of Haux
        for(int q = e.eInfo.qStart; q < e.eInfo.qStop; q++) e.eVars.Haux_eigs[q] = e.eVars.Hsub_eigs[q];
        e.eVars.HauxInitialized = true;
    }

    // Shortcuts
    const ElecInfo& eInfo = e.eInfo;
    ElecVars& eVars = e.eVars;

    cout << "eInfo.mu = " << eInfo.mu << endl;

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

    logPrintf("need_Hsub = %d\n", need_Hsub);

    if( eInfo.fillingsUpdate == ElecInfo::FillingsHsub )
    {        
        // Update nElectrons from mu, or mu from nElectrons as appropriate:
        if( std::isnan(eInfo.mu) )
        {
            logPrintf("Finding mu: ");
            mu = eInfo.findMu( eVars.Haux_eigs, eInfo.nElectrons, Bz );
            logPrintf(" mu = %18.10f\n", mu);
            cout << "Bz = " << Bz << endl;
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
        
        logPrintf("\nBefore update fillings and energies:\n");
        ener.print();

        // Update TS and muN:
        eInfo.updateFillingsEnergies( eVars.Haux_eigs, ener );

        logPrintf("\nAfter update fillings and energies:\n");
        ener.print();

        // Report for SCF (ElecMinimizer handles for minimize):
        if( e.cntrl.scf && eVars.n[0] ) eInfo.smearReport();
    }



    return;
}