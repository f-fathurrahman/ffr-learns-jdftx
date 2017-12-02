#ifndef JDFTX_ELECTRONIC_SCF_H
#define JDFTX_ELECTRONIC_SCF_H

#include <core/Pulay.h>
#include <core/ScalarFieldArray.h>

//! @addtogroup ElecSystem
//! @{
//! @file SCF.h Class SCF and related definitions

//! @brief Variable that is mixed during SCF
//! Component names are density-like, but when mixing potential, they refer to corresponding gradient.
struct SCFvariable
{  ScalarFieldArray n; //!< electron density (or potential)
  ScalarFieldArray tau; //!< KE density (or potential) [mGGA only]
  std::vector<matrix> rhoAtom; //!< atomic density matrices (or corresponding potential) [DFT+U only]
};

//! @brief Self-Consistent Field method for converging electronic state
class SCF : public Pulay<SCFvariable>
{
public:
  SCF(Everything& e);
  
  //! Minimizes residual to achieve self-consistency
  void minimize();
  
  static double eigDiffRMS(const std::vector<diagMatrix>&, const std::vector<diagMatrix>&, const Everything& e); //!< weigted RMS difference between two sets of eigenvalues
  
protected:
  //---- Interface to Pulay ----
  double sync(double x) const;
  double cycle(double dEprev, std::vector<double>& extraValues);
  void report(int iter);
  void axpy(double alpha, const SCFvariable& X, SCFvariable& Y) const;
  double dot(const SCFvariable& X, const SCFvariable& Y) const;
  size_t variableSize() const;
  void readVariable(SCFvariable&, FILE*) const;
  void writeVariable(const SCFvariable&, FILE*) const;
  SCFvariable getVariable() const;
  void setVariable(const SCFvariable&);
  SCFvariable precondition(const SCFvariable&) const;
  SCFvariable applyMetric(const SCFvariable&) const;

private:
  Everything& e;
  bool mixTau; //!< whether KE needs to be mixed
  RealKernel kerkerMix, diisMetric; //!< convolution kernels for kerker preconditioning and the DIIS overlap metric
  
  double eigDiffRMS(const std::vector<diagMatrix>&, const std::vector<diagMatrix>&) const; //!< weighted RMS difference between two sets of eigenvalues
};

//! @}
#endif
