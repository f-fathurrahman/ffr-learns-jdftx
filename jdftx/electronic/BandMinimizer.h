#ifndef JDFTX_ELECTRONIC_BANDMINIMIZER_H
#define JDFTX_ELECTRONIC_BANDMINIMIZER_H

#include <core/Minimize.h>

class ColumnBundle;
class Everything;

//! @addtogroup ElecSystem
//! @{

//! Conjugate gradients eigen-solver
class BandMinimizer : public Minimizable<ColumnBundle>
{
public:
  BandMinimizer(Everything& e, int q); //!< Construct band-structure minimizer for quantum number q

  //Interface for Minimizable:
  double compute(ColumnBundle* grad, ColumnBundle* Kgrad);
  void step(const ColumnBundle& dir, double alpha);
  void constrain(ColumnBundle&);

private:
  Everything& e;
  class ElecVars& eVars;
  const class ElecInfo& eInfo;
  int q; //!< Current quantum number
};

//! @}
#endif
