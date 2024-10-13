// include?

struct MySubspaceRotationAdjust
{
  Everything& e;
  bool adjust; //whether adjustment is active (otherwise only check for Knorm indefiniteness)
  std::vector<matrix> KgPrevHaux; //preconditioned auxiliary gradient at previous step
  double gDotKgPrevHaux; //overlap of current auxiliary gradient with KgPrevHaux
  double cumulatedScale; //net change in subspace rotation factor since last direction reset
  double KnormTot, KnormAux; //latest preconditioned gradient norms: total and auxiliary alone
  
  MySubspaceRotationAdjust(Everything& e)
  : e(e), adjust(e.cntrl.subspaceRotationAdjust), gDotKgPrevHaux(0.), cumulatedScale(1.)
  {
  }
  
  void cacheGradientOverlaps(const MyElecGradient& grad, const MyElecGradient& Kgrad)
  {
    //Calculate overlaps of current gradient:
    KnormTot = dot(grad, Kgrad, &KnormAux);
    mpiWorld->bcast(KnormTot);
    mpiWorld->bcast(KnormAux);
    
    //Compute auxiliary overlap with previous preconditioned gradient:
    if(KgPrevHaux.size())
    {
      gDotKgPrevHaux = 0.;
      for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
        gDotKgPrevHaux += dotc(grad.Haux[q], KgPrevHaux[q]).real();
      mpiWorld->allReduce(gDotKgPrevHaux, MPIUtil::ReduceSum);
      mpiWorld->bcast(gDotKgPrevHaux);
    }
  }
  
  //Handle indefinite Knorm, adjust subspace rotation factor if necessary and report changes.
  //Returns true if CG needs to be reset
  bool report(const std::vector<matrix>& KgradHaux)
  {
    //Handle indefinite preconditioner issues:
    if(KnormTot <= 0.)
    {
      logPrintf("%s\tPreconditioner indefiniteness detected (grad_K will become NAN): ",
        e.elecMinParams.linePrefix);
      my_convergeEmptyStates(e);
      return true; //resets CG
    }
    //Adjust subspace rotation factor if necessary:
    if(adjust)
    {
      if(gDotKgPrevHaux)
      {
        double& kappa = e.cntrl.subspaceRotationFactor;
        double Emin_kappa = gDotKgPrevHaux / KnormTot; //dimensionless version of dEmin/d(log kappa)
        double kappaMax = 0.5*kappa*fabs(KnormTot)/fabs(KnormAux); //aux gradient magnitude should not exceed half that of total
        //Heuristic for adjusting subspace rotation factor (called kappa here):
        double scaleFactor = exp(Emin_kappa/hypot(1,Emin_kappa)); //saturated to range (1/e,e)
        scaleFactor = std::min(scaleFactor, kappaMax/kappa); //cap scale factor at maximum
        kappa *= scaleFactor;
        cumulatedScale *= scaleFactor;
        logPrintf("\tSubspaceRotationAdjust: set factor to %.3lg\n", kappa);
        if(fabs(log(cumulatedScale)) > 2.) //cumulated scale adjustment > e^2
        {
          logPrintf("\tSubspaceRotationAdjust: resetting CG because factor has changed by %lg\n", cumulatedScale);
          cumulatedScale = 1.;
          return true; //resets CG
        }
      }
      KgPrevHaux = KgradHaux; //Cached preconditioned gradient for next step
    }
    return false;
  }
};