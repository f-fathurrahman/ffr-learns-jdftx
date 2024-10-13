#include "my_jdftx.h"

void my_bandMinimize(Everything& e)
{
  bool fixed_H = true;
  std::swap(fixed_H, e.cntrl.fixed_H); //remember fixed_H flag and temporarily set it to true
  logPrintf("Minimization will be done independently for each quantum number.\n");
  e.ener.Eband = 0.;
  for(int q=e.eInfo.qStart; q<e.eInfo.qStop; q++)
  {
    logPrintf("\n---- Minimization of quantum number: ");
    e.eInfo.kpointPrint(globalLog, q, true);
    logPrintf(" ----\n");
    switch(e.cntrl.elecEigenAlgo)
    {
      case ElecEigenCG: { BandMinimizer(e, q).minimize(e.elecMinParams); break; }
      case ElecEigenDavidson: { BandDavidson(e, q).minimize(); break; }
    }
    e.ener.Eband += e.eInfo.qnums[q].weight * trace(e.eVars.Hsub_eigs[q]);
  }
  mpiWorld->allReduce(e.ener.Eband, MPIUtil::ReduceSum);
  if(e.cntrl.shouldPrintEigsFillings)
  {
    //Print the eigenvalues if requested
    print_Hsub_eigs(e);
    logPrintf("\n"); logFlush();
  }
  std::swap(fixed_H, e.cntrl.fixed_H); //restore fixed_H flag
  e.eVars.setEigenvectors();
}

void my_convergeEmptyStates(Everything& e)
{
  logPrintf("Converging empty states (this may take a while): "); logFlush();
  
  std::vector<diagMatrix> eigsPrev = e.eVars.Hsub_eigs;
  logSuspend();
  e.elecMinParams.fpLog = nullLog;
  
  my_bandMinimize(e); //this will also set state to eigenvectors
  logResume();
  e.elecMinParams.fpLog = globalLog;
  
  e.ener.Eband = 0.; //only affects printing (if non-zero Energies::print assumes band structure calc)
  logPrintf("|deigs|: %.3e\n", SCF::eigDiffRMS(e.eVars.Hsub_eigs, eigsPrev, e)); logFlush();
}

//
// Implementation
//

//
// Constructor
//
MyElecMinimizer::MyElecMinimizer(Everything& e)
: e(e),
  eVars(e.eVars),
  eInfo(e.eInfo),
  rotPrev(eInfo.nStates), rotPrevC(eInfo.nStates), rotPrevCinv(eInfo.nStates)
{

  for(int q=eInfo.qStart; q<eInfo.qStop; q++)
  {
    rotPrev[q] = eye(eInfo.nBands);
    rotPrevC[q] = eye(eInfo.nBands);
    rotPrevCinv[q] = eye(eInfo.nBands);
  }
  rotExists = false; //rotation is identity
  
  //Initialize subspace rotation adjuster if required:
  if( e.cntrl.subspaceRotationAdjust &&
      ( eInfo.fillingsUpdate==ElecInfo::FillingsHsub || !eInfo.scalarFillings) )
  {
    //logPrintf("make_shared MySubspaceRotationAdjust is entered\n");
    sra = std::make_shared<MySubspaceRotationAdjust>(e);
  }
}


void MyElecMinimizer::step(const MyElecGradient& dir, double alpha)
{
  assert(dir.eInfo == &eInfo);
  for(int q=eInfo.qStart; q<eInfo.qStop; q++)
  {
    axpy(alpha, rotExists ? dir.C[q]*rotPrevC[q] : dir.C[q], eVars.C[q]);
    if(eInfo.fillingsUpdate==ElecInfo::FillingsConst && eInfo.scalarFillings)
    {
      //Constant scalar fillings: no rotations required
      eVars.orthonormalize(q);
    }
    else
    {
      //Haux or non-scalar fillings: rotations required
      assert(dir.Haux[q]);
      matrix rot;
      if(eInfo.fillingsUpdate == ElecInfo::FillingsHsub)
      {
        //Haux fillings:
        matrix Haux = eVars.Haux_eigs[q];
        axpy(alpha, rotExists ? dagger(rotPrev[q])*dir.Haux[q]*rotPrev[q] : dir.Haux[q], Haux);
        //rotation chosen to diagonalize auxiliary matrix
        Haux.diagonalize(rot, eVars.Haux_eigs[q]);
      }
      else
      {
        //Non-scalar fillings:
        assert(!eInfo.scalarFillings);
        //auxiliary matrix directly generates rotations
        rot = cis(alpha * dir.Haux[q]);
        // cis is defined in core/matrixLinalg.cpp
      }
      matrix rotC = rot;
      eVars.orthonormalize(q, &rotC);
      rotPrev[q] = rotPrev[q] * rot;
      rotPrevC[q] = rotPrevC[q] * rotC;
      rotPrevCinv[q] = inv(rotC) * rotPrevCinv[q];
      rotExists = true; //rotation is no longer identity
    }
  }
}

double MyElecMinimizer::compute(MyElecGradient* grad, MyElecGradient* Kgrad)
{
  if(grad) grad->init(e);  
  if(Kgrad) Kgrad->init(e);
    
  //double ener = e.eVars.elecEnergyAndGrad(e.ener, grad, Kgrad);
  double ener = my_elecEnergyAndGrad(e, grad, Kgrad, false);
    
  if(grad)
  {
    for(int q=eInfo.qStart; q<eInfo.qStop; q++)
    {
      //Rotate wavefunction gradients if necessary:
      if(rotExists)
      {
        grad->C[q] = grad->C[q] * rotPrevCinv[q];
        Kgrad->C[q] = Kgrad->C[q] * rotPrevCinv[q];
      }

      //Subspace gradient handling depends on mode:
      if(eInfo.fillingsUpdate == ElecInfo::FillingsHsub)
      {
        //Haux fillings: rotate gradient computed by ElecVars if necessary
        if(rotExists)
        {
          grad->Haux[q] = rotPrev[q] * grad->Haux[q] * dagger(rotPrev[q]);
          Kgrad->Haux[q] = rotPrev[q] * Kgrad->Haux[q] * dagger(rotPrev[q]);
        }
      }
      else if(!eInfo.scalarFillings)
      {
        //Non-scalar fillings:
        grad->Haux[q] = dagger_symmetrize(complex(0,1) * (eVars.F[q]*eVars.Hsub[q] - eVars.Hsub[q]*eVars.F[q]));
        // dagger_symmetrize is defined in core/matrixOperators.cpp
        Kgrad->Haux[q] = e.cntrl.subspaceRotationFactor * grad->Haux[q];
      }
      //else: constant scalar fillings (no subspace gradient)
    }
        
    //Cache gradient overlaps, if needed, for subspace rotation handling:
    if(sra) sra->cacheGradientOverlaps(*grad, *Kgrad);
    KgradHaux = Kgrad->Haux;
  }
  return ener;
}

bool MyElecMinimizer::report(int iter)
{
  if(e.cntrl.shouldPrintEcomponents)
  {
    //Print the iteration header
    time_t timenow = time(0);
    logPrintf("------------------------------------------------------\n");
    logPrintf("Iteration %d   %s\n",iter, ctime(&timenow));
    //Print the energies
    e.ener.print(); logPrintf("\n"); logFlush();
  }
  
  if(e.cntrl.shouldPrintEigsFillings)
  {
    //Print the eigenvalues if requested
    print_Hsub_eigs(e);
    logPrintf("\n"); logFlush();
  }
    
  //Fillings update report:
  if(eInfo.fillingsUpdate==ElecInfo::FillingsHsub)
  {
    eInfo.smearReport();
  }
    
  //Dump:
  e.dump(DumpFreq_Electronic, iter);
    
  //Re-unitarize rotations:
  if(rotExists) {
    for(int q=eInfo.qStart; q<eInfo.qStop; q++)
    {
      rotPrevC[q] = rotPrev[q];
      rotPrevCinv[q] = dagger(rotPrev[q]);
    }
  }
    
  //Subspace rotation preconditioner handling:
  if(sra) return sra->report(KgradHaux);
  return false;
}



void MyElecMinimizer::constrain(MyElecGradient& dir)
{
  assert(dir.eInfo == &eInfo);
  //Project component of search direction along current wavefunctions:
  for(int q=eInfo.qStart; q<eInfo.qStop; q++) {
    dir.C[q] -= eVars.C[q] * (eVars.C[q]^O(dir.C[q]));
  }
}

//
// All processes minimize together; make sure scalars are in sync to round-off error
//
double MyElecMinimizer::sync(double x) const
{
  mpiWorld->bcast(x);
  return x;
}

// debug the minimize function
double MyElecMinimizer::my_minimize(const MinimizeParams& p)
{
  logPrintf("Entering MyElecMinimizer::my_minimize\n");

  MyElecGradient g, gPrev, Kg;

  // get initial energy and gradient
  double E = sync( compute(&g, &Kg) );
  logPrintf("E = %18.10f\n", E);

  // list of past energies
  EdiffCheck ediffCheck(p.nEnergyDiff, p.energyDiffThreshold);

  // step direction (will be reset in first iteration)
  MyElecGradient d = clone(Kg);
    
  //restrict search direction to allowed subspace
  constrain(d);

  // whether current direction is along the gradient
  bool forceGradDirection = true;
    
  // initially use the specified scheme, may switch to SD on trouble
  MinimizeParams::DirectionUpdateScheme currentDirUpdateScheme = p.dirUpdateScheme;
  bool gPrevUsed;
  switch(currentDirUpdateScheme)
  {
    case MinimizeParams::FletcherReeves:
    case MinimizeParams::SteepestDescent:
      gPrevUsed = false;
      break;
    default:
      gPrevUsed = true;
  }

  double alphaT = p.alphaTstart; //test step size
  double alpha = alphaT; //actual step size
  double beta = 0.0; //CG prev search direction mix factor
  double gKNorm = 0.0, gKNormPrev = 0.0; //current and previous norms of the preconditioned gradient

  // Select the linmin method:
  //Linmin linmin = getLinmin(p);

  //Iterate until convergence, max iteration count or kill signal
  int iter = 0;
  for( iter=0; !killFlag; iter++ )
  //for( iter=0; iter <= 1; iter++ )
  {

      if( report(iter) ) //optional reporting/processing
    {
      E = sync( compute(&g, &Kg) ); //update energy and gradient if state was modified
      logPrintf("%s\tState modified externally: resetting search direction.\n", p.linePrefix);
      forceGradDirection = true; //reset search direction
    }
        
    gKNorm = sync( dot(g,Kg) );
    
    logPrintf("\n");
    logPrintf("%sIter: %3d  %s: ", p.linePrefix, iter, p.energyLabel);
    logPrintf(p.energyFormat, E);
    logPrintf("  |grad|_K: %10.3le  alpha: %10.3le", sqrt(gKNorm/p.nDim), alpha);

    //Print prev step stats and set CG direction parameter if necessary
    beta = 0.0;
    
    if( !forceGradDirection )
    {
      double dotgd = sync(dot(g,d));
      double dotgPrevKg = gPrevUsed ? sync(dot(gPrev, Kg)) : 0.;

      logPrintf("  linmin: %10.3le", dotgd/sqrt(sync(dot(g,g))*sync(dot(d,d))));
      
      if( gPrevUsed )
      {
        logPrintf("  cgtest: %10.3le", dotgPrevKg/sqrt(gKNorm*gKNormPrev));
      }
      
      logPrintf("  t[s]: %9.2lf", clock_sec());

      //Update beta:
      switch(currentDirUpdateScheme)
      {
        case MinimizeParams::FletcherReeves:  beta = gKNorm/gKNormPrev; break;
        case MinimizeParams::PolakRibiere:    beta = (gKNorm-dotgPrevKg)/gKNormPrev; break;
        case MinimizeParams::HestenesStiefel: beta = (gKNorm-dotgPrevKg)/(dotgd-sync(dot(d,gPrev))); break;
        case MinimizeParams::SteepestDescent: beta = 0.0; break;
        //Should never encounter since LBFGS handled separately; just to eliminate compiler warnings
        case MinimizeParams::LBFGS: break;
      }

      if(beta < 0.0)
      {
        logPrintf("\n%sEncountered beta < 0, resetting CG.", p.linePrefix);
        beta = 0.0;
      }
    }
        
    forceGradDirection = false;
        
    logPrintf("\n");
        
    if( sqrt(gKNorm/p.nDim) < p.knormThreshold )
    {
      logPrintf("%sConverged (|grad|_K < %le).\n", p.linePrefix, p.knormThreshold);
      return E;
    }

    if( ediffCheck.checkConvergence(E) )
    {
      logPrintf("%sConverged (|Delta %s| < %le for %d iters).\n",
        p.linePrefix, p.energyLabel, p.energyDiffThreshold, p.nEnergyDiff);
      return E;
    }

    if(!std::isfinite(gKNorm))
    {
      logPrintf("%s|grad|_K = %le. Stopping ...\n", p.linePrefix, gKNorm);
      return E;
    }

    if(!std::isfinite(E))
    {
      logPrintf("%sE = %le. Stopping ...\n", p.linePrefix, E);
      return E;
    }
    
    if( iter >= p.nIterations) {
      break;
    }

    // We are not stopping the iterations ...

    if( gPrevUsed ) {
      gPrev = g;
    }
    
    gKNormPrev = gKNorm;

    // Update search direction
    d *= beta;
    axpy(-1.0, Kg, d);  // d = beta*d - Kg
        
    //restrict search direction to allowed subspace
    constrain(d);
    
    // Line minimization
    alphaT = std::min(alphaT, safeStepSize(d));
    if( my_linminQuad(*this, p, d, alphaT, alpha, E, g, Kg) )
    {
      // linmin succeeded:
      if(p.updateTestStepSize)
      {
        alphaT = alpha;

        if(alphaT < p.alphaTmin) {
          // bad step size: make sure next test step size is not too bad
          alphaT = p.alphaTstart; 
        }
      }
    }
    else
    {
      // linmin failed:
      logPrintf("%s\tUndoing step.\n", p.linePrefix);
      step(d, -alpha);
      
      E = sync(compute(&g, &Kg));
      
      if(beta)
      {
        // Failed, but not along the gradient direction:
        logPrintf("%s\tStep failed: resetting search direction.\n", p.linePrefix);
        forceGradDirection = true; //reset search direction
      }
      else
      {
        // Failed along the gradient direction
        logPrintf("%s\tStep failed along negative gradient direction.\n", p.linePrefix);
        logPrintf("%sProbably at roundoff error limit. (Stopping)\n", p.linePrefix);
        return E;
      }
    }

  } // end of iter

  logPrintf("Leaving MyElecMinimizer::my_minimize\n");

  return 0.0;
}
