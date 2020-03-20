struct MyElecGradient
{
  std::vector<ColumnBundle> C; //!< wavefunctions
  std::vector<matrix> Haux; //!< auxiliary Hamiltonian
  const ElecInfo* eInfo;
  
  void init(const Everything& e); //!< initialize C and Haux with the correct sizes for everything
  
  MyElecGradient& operator*=(double alpha); //!< scalar multiply
};

double dot(const MyElecGradient& x, const MyElecGradient& y, double* auxContrib)
{
  assert(x.eInfo == y.eInfo);
  
  std::vector<double> result(2, 0.); //calculate wavefunction and auxiliary contributions separately
  
  for(int q=x.eInfo->qStart; q<x.eInfo->qStop; q++)
  {
    if(x.C[q] && y.C[q]) result[0] += dotc(x.C[q], y.C[q]).real()*2.0;
    if(x.Haux[q] && y.Haux[q]) result[1] += dotc(x.Haux[q], y.Haux[q]).real();
  }
  mpiWorld->allReduceData(result, MPIUtil::ReduceSum);
  
  if(auxContrib) *auxContrib=result[1]; //store auxiliary contribution, if requested
  
  return result[0]+result[1]; //return total
}