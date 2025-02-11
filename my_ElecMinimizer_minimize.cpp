#include "my_jdftx.h"

double my_ElecMinimizer_minimize(
  Everything& e,
  MyElecMinimizer& elecMin,
  const MinimizeParams& p)
{
  logPrintf("-------------------------------\n");
  logPrintf("ENTER my_ElecMinimizer_minimize\n");
  logPrintf("-------------------------------\n");

  MyElecGradient g, gPrev, Kg;

  // get initial energy and gradient
  double E = elecMin.compute(&g, &Kg);
  logPrintf("Initial energy E = %18.10f\n", E);

  // list of past energies
  EdiffCheck ediffCheck(p.nEnergyDiff, p.energyDiffThreshold);

  // step direction (will be reset in first iteration)
  MyElecGradient d = clone(Kg);

  //restrict search direction to allowed subspace
  elecMin.constrain(d);

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

    if( elecMin.report(iter) ) //optional reporting/processing
    {
      E = elecMin.compute(&g, &Kg); //update energy and gradient if state was modified
      logPrintf("%s\tState modified externally: resetting search direction.\n", p.linePrefix);
      forceGradDirection = true; //reset search direction
    }
        
    gKNorm = dot(g,Kg);
    
    logPrintf("\n");
    logPrintf("%sIter: %3d  %s: ", p.linePrefix, iter, p.energyLabel);
    logPrintf(p.energyFormat, E);
    logPrintf("  |grad|_K: %10.3le  alpha: %10.3le", sqrt(gKNorm/p.nDim), alpha);

    //Print prev step stats and set CG direction parameter if necessary
    beta = 0.0;
    
    if( !forceGradDirection )
    {
      double dotgd = dot(g,d);
      double dotgPrevKg = gPrevUsed ? dot(gPrev, Kg) : 0.0;

      logPrintf("  linmin: %10.3le", dotgd/sqrt( dot(g,g) * dot(d,d)));
      
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
        case MinimizeParams::HestenesStiefel: beta = (gKNorm-dotgPrevKg)/(dotgd-dot(d,gPrev)); break;
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
    elecMin.constrain(d);
    
    // Line minimization
    alphaT = std::min(alphaT, elecMin.safeStepSize(d));
    if( my_linminQuad(elecMin, p, d, alphaT, alpha, E, g, Kg) )
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
      elecMin.step(d, -alpha);
      
      E = elecMin.compute(&g, &Kg);
      
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


  logPrintf("------------------------------\n");
  logPrintf("EXIT my_ElecMinimizer_minimize\n");
  logPrintf("------------------------------\n");

  return 0.0;
}


