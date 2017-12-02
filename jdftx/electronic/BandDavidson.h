#ifndef JDFTX_ELECTRONIC_BANDDAVIDSON_H
#define JDFTX_ELECTRONIC_BANDDAVIDSON_H

#include <core/Minimize.h>

class Everything;

//! @addtogroup ElecSystem
//! @{

//! Davidson eigensolver
class BandDavidson
{
public:
  BandDavidson(Everything& e, int q); //!< Construct Davidson eigenvalue solver for quantum number q
  void minimize(); //!< Converge eigenproblem with tolerance set by e.elecMinParams
  
private:
  Everything& e;
  class ElecVars& eVars;
  const class ElecInfo& eInfo;
  int q;  //!< Current quantum number
};

//! @}
#endif // JDFTX_ELECTRONIC_BANDDAVIDSON_H
