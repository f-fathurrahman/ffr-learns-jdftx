#include "my_jdftx.h"

bool my_linminQuad(
  Everything& e,
  MyElecMinimizer& elecMin,
  const MinimizeParams& p,
  const MyElecGradient& d,
  double alphaT,
  double& alpha,
  double& E,
  MyElecGradient& g,
  MyElecGradient& Kg )
{
        
  double alphaPrev = 0.0; //the progress made so far along d
  double Eorig = E;
  double gdotd = dot(g,d); //directional derivative at starting point
  
  if( gdotd >= 0.0 ) {
    logPrintf("%s\tBad step direction: g.d > 0.\n", p.linePrefix);
    alpha = alphaPrev;
    logPrintf("%s\tUsing previous alpha = %f\n", p.linePrefix, alpha);
    return false;
  }

  //Test step and step size prediction:
  double ET = 0.0; //test step energy
  for(int s=0; s < p.nAlphaAdjustMax; s++) {
    //
    logPrintf("my_linMinQuad: trial step = %d alphaT=%le alphaPrev=%le\n", s, alphaT, alphaPrev);
    if(alphaT < p.alphaTmin) {
      logPrintf("%s\talphaT below threshold %le. Quitting step.\n", p.linePrefix, p.alphaTmin);
      alpha = alphaPrev;
      return false;
    }
    
    // Try the test step:
    my_ElecMinimizer_step(e, elecMin, d, alphaT-alphaPrev);
    alphaPrev = alphaT;
    ET = my_ElecMinimizer_compute(e, elecMin, 0, 0);
    
    // Check if step crossed domain of validity of parameter space:
    if( !std::isfinite(ET) ) {
      alphaT *= p.alphaTreduceFactor;
      logPrintf("%s\tTest step failed with %s = %le, reducing alphaT to %le.\n",
                p.linePrefix, p.energyLabel, ET, alphaT);
      continue;
    }
    //Predict step size:
    alpha = 0.5*pow(alphaT,2)*gdotd/(alphaT*gdotd + E - ET);

    //Check reasonableness of predicted step size:
    if(alpha < 0) {
      // Curvature has the wrong sign
      // That implies ET < E, so accept step for now, and try descending further next time
      alphaT *= p.alphaTincreaseFactor;
      
      logPrintf("\tWrong curvature in test step, increasing alphaT to %le.\n", alphaT);
      //
      // ffr: why we do early return here? alphaT is modified for the next iteration? 
      E = my_ElecMinimizer_compute(e, elecMin, &g, &Kg);
      return true;
    }
    
    if(alpha/alphaT > p.alphaTincreaseFactor) {
      alphaT *= p.alphaTincreaseFactor;
      logPrintf("%s\tPredicted alpha/alphaT>%lf, increasing alphaT to %le.\n",
                p.linePrefix, p.alphaTincreaseFactor, alphaT);
      continue; // continue iteration ?
    }
    
    if(alphaT/alpha < p.alphaTreduceFactor) {
      alphaT *= p.alphaTreduceFactor;
      logPrintf("%s\tPredicted alpha/alphaT < %lf, reducing alphaT to %le.\n",
                p.linePrefix, p.alphaTreduceFactor, alphaT);
      continue; // continue iteration
    }
    
    //Successful test step:
    logPrintf("my_linminQuad: successful trial step: using alpha = %le\n", alpha);
    break;
  }



  if(!std::isfinite(E)) {
    logPrintf("%s\tTest step failed %d times. Quitting step.\n", p.linePrefix, p.nAlphaAdjustMax);
    alpha = alphaPrev;
    return false;
  }

  // Actual step:
  for(int s=0; s < p.nAlphaAdjustMax; s++) {

    logPrintf("my_linminQuad: Actual step, iter=%d using alpha = %le\n", s, alpha);

    //Try the step:
    my_ElecMinimizer_step(e, elecMin, d, alpha-alphaPrev);
    alphaPrev = alpha;

    E = my_ElecMinimizer_compute(e, elecMin, &g, &Kg);
    
    if( !std::isfinite(E) ) {
      alpha *= p.alphaTreduceFactor;
      logPrintf("\tStep failed with %s = %le, reducing alpha to %le.\n", p.energyLabel, E, alpha);
      continue;
    }
    
    // Trial energy is higher, reduce alpha
    if(E > Eorig) {
      alpha *= p.alphaTreduceFactor;
      logPrintf("\tStep increased %s by %le, reducing alpha to %le.\n", p.energyLabel, E-Eorig, alpha);
      continue;
    }

    // Step successful:
    break;
  }
  
  if(!std::isfinite(E) || E>Eorig) {
    logPrintf("\tStep failed to reduce %s after %d attempts. Quitting step.\n", p.energyLabel, p.nAlphaAdjustMax);
    fflush(p.fpLog);
    return false;
  }
  
  return true;
}


