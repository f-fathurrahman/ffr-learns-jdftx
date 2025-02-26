#include "my_jdftx.h"

void my_ElecMinimizer_step(
  Everything& e,
  MyElecMinimizer& elecMin,
  const MyElecGradient& dir,
  double alpha
)
{
  assert( dir.eInfo == &e.eInfo );

  bool is_rot_exist_prev = elecMin.rotExists;

  logPrintf("**** ENTER my_ElecMinimizer_step with alpha=%le\n", alpha);
  if(elecMin.rotExists) {
    logPrintf("my_ElecMinimizer_step: rotations will be applied\n");
  } else {
    logPrintf("my_ElecMinimizer_step: rotations will NOT be applied\n");
  }

  for(int q=elecMin.eInfo.qStart; q < elecMin.eInfo.qStop; q++) {
    //
    // Update step for wavefunctions
    //
    if(elecMin.rotExists) {
      axpy( alpha, dir.C[q] * elecMin.rotPrevC[q], e.eVars.C[q] );
    }
    else {
      axpy( alpha, dir.C[q], e.eVars.C[q] );
    }
    //
    // Now for Haux (if needed)
    //
    if( e.eInfo.fillingsUpdate==ElecInfo::FillingsConst && e.eInfo.scalarFillings) {
      // Constant scalar fillings: no rotations required
      e.eVars.orthonormalize(q);
    }
    else {
      // Haux or non-scalar fillings: rotations required
      assert( dir.Haux[q] );
      matrix rot; //?
      if( e.eInfo.fillingsUpdate == ElecInfo::FillingsHsub) {
        //Haux fillings:
        matrix Haux = e.eVars.Haux_eigs[q];
        if(elecMin.rotExists) {
          axpy(alpha, dagger( elecMin.rotPrev[q] ) * dir.Haux[q] * elecMin.rotPrev[q], Haux);
        }
        else {
          axpy(alpha, dir.Haux[q], Haux);
        }
        //rotation chosen to diagonalize auxiliary matrix
        Haux.diagonalize( rot, e.eVars.Haux_eigs[q] );
        // results are in rot and Haux_eigs?
      }
      else {
        //Non-scalar fillings:
        assert( !e.eInfo.scalarFillings );
        //auxiliary matrix directly generates rotations
        rot = cis(alpha * dir.Haux[q]);
        // cis is defined in core/matrixLinalg.cpp
      }
      matrix rotC = rot; // copy eigenvectors from Haux
      e.eVars.orthonormalize(q, &rotC); // orthonormalize wavefunc with extra rotations
      elecMin.rotPrev[q] = elecMin.rotPrev[q] * rot; // accumulate
      elecMin.rotPrevC[q] = elecMin.rotPrevC[q] * rotC; // rotPrev for wavefunc
      elecMin.rotPrevCinv[q] = inv(rotC) * elecMin.rotPrevCinv[q];
      elecMin.rotExists = true; //rotation is no longer identity
    }
  }

  if(!is_rot_exist_prev) {
    if(elecMin.rotExists) {
      logPrintf("!!!! my_ElecMinimizer_step: rotExists is set to true and it is set to false before\n");
    }
  }

  logPrintf(" **** EXIT my_ElecMinimizer_step\n");
}