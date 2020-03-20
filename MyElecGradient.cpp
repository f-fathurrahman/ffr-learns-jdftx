struct MyElecGradient
{
  std::vector<ColumnBundle> C; //!< wavefunctions
  std::vector<matrix> Haux; //!< auxiliary Hamiltonian
  const ElecInfo* eInfo;
  
  void init(const Everything& e); //!< initialize C and Haux with the correct sizes for everything
  
  MyElecGradient& operator*=(double alpha); //!< scalar multiply
};

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
  
//  if(auxContrib) *auxContrib = result[1]; //store auxiliary contribution, if requested
  
  return result[0] + result[1]; //return total
}

MyElecGradient clone(const MyElecGradient& x)
{
  return x;
}

void randomize(MyElecGradient& x)
{
  randomize(x.C, *x.eInfo);
  for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++) {
    if(x.Haux[q])
    {
      randomize(x.Haux[q]);
      x.Haux[q] = dagger_symmetrize(x.Haux[q]); //make hermitian
    }
  }
}
