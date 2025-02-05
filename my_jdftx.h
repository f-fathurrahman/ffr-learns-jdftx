#ifndef MY_JDFTX_H
#define MY_JDFTX_H

#include <vector>

#include <core/Minimize.h>

#include <electronic/ElecInfo.h>
#include <electronic/ColumnBundle.h>
#include <electronic/Everything.h>
#include <electronic/BandMinimizer.h>
#include <electronic/BandDavidson.h>
#include <electronic/ElecInfo.h>
#include <electronic/ColumnBundle.h>
#include <electronic/SCF.h>

#include <core/ScalarFieldArray.h>
#include <core/Util.h>

#include <iostream>

//
// Definition of MyElecGradient
//
struct MyElecGradient
{
  std::vector<ColumnBundle> C; //!< wavefunctions
  std::vector<matrix> Haux; //!< auxiliary Hamiltonian
  const ElecInfo* eInfo;
  
  void init(const Everything& e); //!< initialize C and Haux with the correct sizes for everything
  
  MyElecGradient& operator*=(double alpha); //!< scalar multiply
};

//
// Some functions
//

void my_convergeEmptyStates(Everything& e);

void axpy(double alpha, const MyElecGradient& x, MyElecGradient& y);

double dot(const MyElecGradient& x, const MyElecGradient& y, double* auxContrib);

double dot(const MyElecGradient& x, const MyElecGradient& y);

void my_bandMinimize(Everything& e);

double my_elecEnergyAndGrad( Everything& e, 
  MyElecGradient* grad, MyElecGradient* Kgrad, bool calc_Hsub );

bool my_linminQuad(
  Minimizable<MyElecGradient>& obj,
  const MinimizeParams& p,
  const MyElecGradient& d,
  double alphaT, double& alpha, double& E, MyElecGradient& g, MyElecGradient& Kg );

MyElecGradient clone(const MyElecGradient& x);
void randomize(MyElecGradient& x);




//
// This is a struct
//
#include "MySubspaceRotationAdjust.cpp"
// Can we make this a class rather than a struct?
// Are there any reasons related to shared_pointer ?




//
//! Variational total energy minimizer for electrons
//
class MyElecMinimizer : public Minimizable<MyElecGradient>
{
public:
  MyElecMinimizer(Everything& e);
  
  //Virtual functions from Minimizable:
  void step(const MyElecGradient& dir, double alpha);
  double compute(MyElecGradient* grad, MyElecGradient* Kgrad);
  bool report(int iter);
  void constrain(MyElecGradient&);

  double my_minimize(const MinimizeParams& p);

  //!< All processes minimize together; make sure scalars are in sync to round-off error
  double sync(double x) const;

  // XXX They are originally private
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
// Functions
//
void my_elecFluidMinimize( Everything& e );

void write_info(Everything &e);

ScalarFieldArray my_calcDensity(Everything& e);

double my_applyHamiltonian(
  Everything& e,
  int q,
  const diagMatrix& Fq,
  ColumnBundle& HCq,
  Energies& ener,
  bool need_Hsub
);

void my_ElecMinimizer_step(
  Everything& e,
  MyElecMinimizer& elecMin,
  MyElecGradient& dir,
  double alpha
);

void write_eVars_F(Everything& e, const char* prefix="eVars_F_");
void write_eVars_Hsub(Everything& e, const char* prefix="eVars_Hsub_");
void write_eVars_n(Everything& e, const char* prefix="eVars_n_");

#endif

