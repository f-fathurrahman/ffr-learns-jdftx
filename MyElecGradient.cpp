#include "my_jdftx.h"

void MyElecGradient::init(const Everything& e)
{
  eInfo = &e.eInfo;
  C.resize(eInfo->nStates);
  Haux.resize(eInfo->nStates);
}

MyElecGradient& MyElecGradient::operator*=(double alpha)
{
  for(int q=eInfo->qStart; q<eInfo->qStop; q++)
  {
    if(C[q]) C[q] *= alpha;
    if(Haux[q]) Haux[q] *= alpha;
  }
  return *this;
}

void axpy(double alpha, const MyElecGradient& x, MyElecGradient& y)
{
  assert(x.eInfo == y.eInfo);
  for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
  {
    if(x.C[q]) {
      if(y.C[q]) {
        axpy(alpha, x.C[q], y.C[q]);
      }
      else {
        y.C[q] = alpha*x.C[q];
      }
    }
    
    if(x.Haux[q]) {
      if(y.Haux[q]) {
        axpy(alpha, x.Haux[q], y.Haux[q]);
      }
      else {
        y.Haux[q] = alpha*x.Haux[q];
      }
    }
  }
}

double dot(const MyElecGradient& x, const MyElecGradient& y, double* auxContrib)
{
  assert(x.eInfo == y.eInfo);
  
  std::vector<double> result(2, 0.); //calculate wavefunction and auxiliary contributions separately
  
  for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
  {
    if(x.C[q] && y.C[q]) {
      result[0] += dotc(x.C[q], y.C[q]).real()*2.0;
    }
    
    if(x.Haux[q] && y.Haux[q]) {
      result[1] += dotc(x.Haux[q], y.Haux[q]).real();
    }
  }
  mpiWorld->allReduceData(result, MPIUtil::ReduceSum);
  
  if(auxContrib) *auxContrib = result[1]; //store auxiliary contribution, if requested
  
  return result[0] + result[1]; //return total
}


// No aux contrib
double dot(const MyElecGradient& x, const MyElecGradient& y)
{
  assert(x.eInfo == y.eInfo);
  
  std::vector<double> result(2, 0.); //calculate wavefunction and auxiliary contributions separately
  
  //logPrintf("Initial: res = %18.10f, %18.10f\n", result[0], result[1]);

  for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
  {
    if(x.C[q] && y.C[q]) {
      result[0] += dotc(x.C[q], y.C[q]).real()*2.0;
    }
    
    if(x.Haux[q] && y.Haux[q]) {
      result[1] += dotc(x.Haux[q], y.Haux[q]).real();
    }
    //logPrintf("q = %d, res = %18.10f\n", q, dotc(x.C[q], y.C[q]).real()*2.0);
    //logPrintf("q = %d, res = %18.10f, %18.10f\n", q, result[0], result[1]);
  }
  mpiWorld->allReduceData(result, MPIUtil::ReduceSum);
  
//  if(auxContrib) *auxContrib = result[1]; //store auxiliary contribution, if requested
  
  return result[0] + result[1]; //return total
}

MyElecGradient clone(const MyElecGradient& x)
{
  return x;
}

void randomize(MyElecGradient& x)
{
  logPrintf("Calling MyElecGradient.randomize\n");
  randomize(x.C, *x.eInfo);
  for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++) {
    if(x.Haux[q])
    {
      randomize(x.Haux[q]);
      x.Haux[q] = dagger_symmetrize(x.Haux[q]); //make hermitian
    }
  }
}
