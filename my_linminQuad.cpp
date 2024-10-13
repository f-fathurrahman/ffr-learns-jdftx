#include "my_jdftx.h"

bool my_linminQuad(
  Minimizable<MyElecGradient>& obj,
  const MinimizeParams& p,
  const MyElecGradient& d,
  double alphaT, double& alpha, double& E, MyElecGradient& g, MyElecGradient& Kg )
{
        
  double alphaPrev = 0.0; //the progress made so far along d
  double Eorig = E;
  double gdotd = obj.sync(dot(g,d)); //directional derivative at starting point
  
  if( gdotd >= 0.0 )
  {
    logPrintf("%s\tBad step direction: g.d > 0.\n", p.linePrefix);
    alpha = alphaPrev;
    return false;
  }

  //Test step and step size prediction:
  double ET = 0.0; //test step energy
  for(int s=0; s<p.nAlphaAdjustMax; s++)
  {
    if(alphaT < p.alphaTmin)
    {
      logPrintf("%s\talphaT below threshold %le. Quitting step.\n", p.linePrefix, p.alphaTmin);
      alpha = alphaPrev;
      return false;
    }
    
    // Try the test step:
    obj.step(d, alphaT-alphaPrev); alphaPrev = alphaT;
    ET = obj.sync( obj.compute(0,0) );
    
    // Check if step crossed domain of validity of parameter space:
    if(!std::isfinite(ET))
    {
      alphaT *= p.alphaTreduceFactor;
      logPrintf("%s\tTest step failed with %s = %le, reducing alphaT to %le.\n",
                p.linePrefix, p.energyLabel, ET, alphaT);
      continue;
    }
    //Predict step size:
    alpha = 0.5*pow(alphaT,2)*gdotd/(alphaT*gdotd + E - ET);

    //Check reasonableness of predicted step size:
    if(alpha < 0)
    {
      // Curvature has the wrong sign
      // That implies ET < E, so accept step for now, and try descending further next time
      alphaT *= p.alphaTincreaseFactor;
      
      logPrintf("%s\tWrong curvature in test step, increasing alphaT to %le.\n", p.linePrefix, alphaT);
      
      E = obj.sync(obj.compute(&g, &Kg));
      return true;
    }
    
    if(alpha/alphaT > p.alphaTincreaseFactor)
    {
      alphaT *= p.alphaTincreaseFactor;
      logPrintf("%s\tPredicted alpha/alphaT>%lf, increasing alphaT to %le.\n",
                p.linePrefix, p.alphaTincreaseFactor, alphaT);
      continue;
    }
    
    if(alphaT/alpha < p.alphaTreduceFactor)
    {
      alphaT *= p.alphaTreduceFactor;
      logPrintf("%s\tPredicted alpha/alphaT < %lf, reducing alphaT to %le.\n",
                p.linePrefix, p.alphaTreduceFactor, alphaT);
      continue;
    }
    
    //Successful test step:
    break;
  }



  if(!std::isfinite(E))
  {
    logPrintf("%s\tTest step failed %d times. Quitting step.\n", p.linePrefix, p.nAlphaAdjustMax);
    alpha = alphaPrev;
    return false;
  }

  // Actual step:
  for(int s=0; s < p.nAlphaAdjustMax; s++)
  {
    //Try the step:
    obj.step(d, alpha-alphaPrev);
    alphaPrev = alpha;

    E = obj.sync(obj.compute(&g, &Kg));
    
    if(!std::isfinite(E))
    {
      alpha *= p.alphaTreduceFactor;
      logPrintf("%s\tStep failed with %s = %le, reducing alpha to %le.\n",
                p.linePrefix, p.energyLabel, E, alpha);
      continue;
    }
    
    if(E > Eorig)
    {
      alpha *= p.alphaTreduceFactor;
      logPrintf("%s\tStep increased %s by %le, reducing alpha to %le.\n",
                p.linePrefix, p.energyLabel, E-Eorig, alpha);
      continue;
    }

    // Step successful:
    break;
  }
  
  if(!std::isfinite(E) || E>Eorig)
  {
    logPrintf("%s\tStep failed to reduce %s after %d attempts. Quitting step.\n",
              p.linePrefix, p.energyLabel, p.nAlphaAdjustMax); fflush(p.fpLog);
    return false;
  }
  
  return true;
}


