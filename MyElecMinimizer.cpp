#include <core/Minimize.h>
#include <electronic/BandMinimizer.h>
#include <electronic/BandDavidson.h>
#include <electronic/Everything.h>
#include <electronic/ElecInfo.h>
#include <electronic/ColumnBundle.h>
#include <electronic/SCF.h>

#include "MyElecGradient.cpp"

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

#include "MySubspaceRotationAdjust.cpp"


//! Variational total energy minimizer for electrons
class MyElecMinimizer : public Minimizable<MyElecGradient>
{
public:
  MyElecMinimizer(Everything& e);
  
  //Virtual functions from Minimizable:
  void step(const MyElecGradient& dir, double alpha);
  double compute(MyElecGradient* grad, MyElecGradient* Kgrad);
  bool report(int iter);
  void constrain(MyElecGradient&);

  //!< All processes minimize together; make sure scalars are in sync to round-off error
  double sync(double x) const;
  
private:
  Everything& e;
  class ElecVars& eVars;
  const ElecInfo& eInfo;
  std::vector<matrix> KgradHaux; //!< latest preconditioned auxiliary gradient
  std::vector<matrix> rotPrev; //!< cumulated unitary rotations of subspace
  std::vector<matrix> rotPrevC; //!< cumulated transormation of wavefunctions (including non-unitary orthonormalization components)
  std::vector<matrix> rotPrevCinv; //!< inverse of rotPrevC (which is not just dagger, since these are not exactly unitary)
  
  bool rotExists; //!< whether rotPrev is non-trivial (not identity)
  std::shared_ptr<struct MySubspaceRotationAdjust> sra; //!< Subspace rotation adjustment helper
};

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
    logPrintf("make_shared MySubspaceRotationAdjust is entered\n");
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

#include "my_elecEnergyAndGrad.cpp"

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