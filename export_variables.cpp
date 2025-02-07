#include "my_jdftx.h"

// formatted output
void export_eVars_Hsub(Everything& e, const char* prefix) {
  for(int q = e.eInfo.qStart; q < e.eInfo.qStop; q++) {
    std::stringstream ss;
    ss << prefix << q+1 << ".dat";
    FILE *fptr;
    fptr = fopen( ss.str().c_str(), "w");
    e.eVars.Hsub[q].print(fptr);
    fclose(fptr);
    std::cout << "Hsub is written to " << ss.str() << std::endl;
  }
}

void write_eVars_F(Everything& e, const char* prefix) {
  for(int q = e.eInfo.qStart; q < e.eInfo.qStop; q++) {
    std::stringstream ss;
    ss << prefix << q+1 << ".dat";
    FILE *fptr;
    fptr = fopen( ss.str().c_str(), "w");
    e.eVars.F[q].print(fptr);
    fclose(fptr);
    std::cout << "eVars.F is written to " << ss.str() << std::endl;
  }
}

void write_eVars_n(Everything& e, const char* prefix) {
  int Nspin = e.eVars.n.size();
  for(int ispin = 0; ispin < Nspin; ispin++) {
    matrix n = e.eVars.n[ispin]->toMatrix();
    std::stringstream ss;
    ss << prefix << ispin+1 << ".dat";
    FILE *fptr;
    fptr = fopen( ss.str().c_str(), "w");
    n.print(fptr, "%18.10lg\t"); // only real parts ?
    //n.print(fptr, "%18.10lg%+18.10lgi\t"); // both real and imag parts
    fclose(fptr);
    std::cout << "eVars.n is written to " << ss.str() << std::endl;
  }
}
