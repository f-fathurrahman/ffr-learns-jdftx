/*-------------------------------------------------------------------
Copyright 2011 Ravishankar Sundararaman

This file is part of JDFTx.

JDFTx is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

JDFTx is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with JDFTx.  If not, see <http://www.gnu.org/licenses/>.
-------------------------------------------------------------------*/

#include <commands/command.h>
#include <electronic/Everything.h>


EnumStringMap<ExCorrType> exCorrTypeMap(
  ExCorrLDA_PZ, "lda", //Default LDA is PZ
  ExCorrLDA_PZ, "lda-PZ", 
  ExCorrLDA_PW, "lda-PW", 
  ExCorrLDA_PW_prec, "lda-PW-prec", 
  ExCorrLDA_VWN, "lda-VWN", 
  ExCorrLDA_Teter, "lda-Teter",
  ExCorrGGA_PBE, "gga", //Default GGA is PBE
  ExCorrGGA_PBE, "gga-PBE",
  ExCorrGGA_PBEsol, "gga-PBEsol",
  ExCorrGGA_PW91, "gga-PW91",
  ExCorrMGGA_TPSS, "mgga-TPSS",
  ExCorrMGGA_revTPSS, "mgga-revTPSS",
  ExCorrORB_GLLBsc, "orb-GLLBsc",
  ExCorrPOT_LB94, "pot-LB94",
  ExCorrHYB_PBE0, "hyb-PBE0",
  ExCorrHYB_HSE06, "hyb-HSE06",
  ExCorrHYB_HSE12, "hyb-HSE12",
  ExCorrHYB_HSE12s, "hyb-HSE12s",
  ExCorrHF, "Hartree-Fock"
);
EnumStringMap<ExCorrType> exCorrDescriptionMap(
  ExCorrLDA_PZ, "Perdew-Zunger LDA", 
  ExCorrLDA_PW, "Perdew-Wang LDA", 
  ExCorrLDA_PW_prec, "Perdew-Wang LDA with extended precision (used by PBE)",
  ExCorrLDA_VWN, "Vosko-Wilk-Nusair LDA", 
  ExCorrLDA_Teter, "Teter93 LSDA",
  ExCorrGGA_PBE, "Perdew-Burke-Ernzerhof GGA",
  ExCorrGGA_PBEsol, "Perdew-Burke-Ernzerhof GGA reparametrized for solids",
  ExCorrGGA_PW91, "Perdew-Wang GGA",
  ExCorrMGGA_TPSS, "Tao-Perdew-Staroverov-Scuseria meta GGA",
  ExCorrMGGA_revTPSS, "revised Tao-Perdew-Staroverov-Scuseria meta GGA",
  ExCorrORB_GLLBsc, "Orbital-dependent GLLB-sc potential (no total energy)",
  ExCorrPOT_LB94, "van Leeuwen-Baerends model potential (no total energy)",
  ExCorrHYB_PBE0, "Hybrid PBE with 1/4 exact exchange",
  ExCorrHYB_HSE06, "HSE06 'wPBEh' hybrid with 1/4 screened exact exchange",
  ExCorrHYB_HSE12, "Reparametrized screened exchange functional for accuracy (w=0.185 A^-1 and a=0.313)",
  ExCorrHYB_HSE12s, "Reparametrized screened exchange functional for k-point convergence (w=0.408 A^-1 and a=0.425)",
  ExCorrHF, "Full exact exchange with no correlation"
);

EnumStringMap<KineticType> kineticTypeMap(
  KineticTF, "lda-TF",
  KineticVW, "gga-vW", 
  KineticPW91, "gga-PW91k"
);
EnumStringMap<KineticType> kineticDescriptionMap(
  KineticTF, "Thomas-Fermi LDA kinetic energy", 
  KineticVW, "von Weisacker correction to LDA kinetic energy", 
  KineticPW91, "Perdew-Wang GGA kinetic energy parameterized by Lembarki and Chermette"
);


struct CommandElecExCorr : public Command
{
  CommandElecExCorr(const char* cmdName = "elec-ex-corr", string path="jdftx/Electronic/Functional") : Command(cmdName, path)
  {
    format = "<functional>";
    comments = "Specify the exchange-correlation functional, where <functional> is one of:"
      + addDescriptions(exCorrTypeMap.optionList(), linkDescription(exCorrTypeMap, exCorrDescriptionMap))
      + ".\n\nNote that lda is an alias for lda-pz, and gga for gga-pbe.\n\n";
    hasDefault = true;
    emptyParamError = "   eXchange/Correlation functional(s) must be specified.";
    
    #ifdef LIBXC_ENABLED
    format += "\n\t| <funcX> <funcC>\n\t| <funcXC>";
    comments +=
      "The second and third lines use eXchange/Correlation functionals from libXC \\cite LibXC.\n"
      "The exact entries below will depend on the version of LibXC linked against.\n"
      "Here, <funcX> is one of:"
      + addDescriptions(xcMap_X.optionList(), getLibXCdescription_X)
      + ",\n\n<funcC> is one of:"
      + addDescriptions(xcMap_C.optionList(), getLibXCdescription_C)
      + ",\n\nand <funcXC> is one of:"
      + addDescriptions(xcMap_XC.optionList(), getLibXCdescription_XC)
      + ".";
    #else
    comments += "Additional functionals can be enabled by compiling with LibXC support \\cite LibXC.";
    #endif
  }

  void process(ParamList& pl, Everything& e)
  {  process(pl, e.exCorr);
  }
  
  void printStatus(Everything& e, int iRep)
  {  printStatus(e.exCorr);
  }
  
protected:
  void process(ParamList& pl, ExCorr& exCorr)
  {  string key;
    pl.get(key, string(), "functional");
    if(key.length()) //Otherwise default functional set by ExCorr constructor
    {  if(exCorrTypeMap.getEnum(key.c_str(), exCorr.exCorrType)) //Found internal ExCorr functional
      {  //Set the functional name:
        exCorr.xcName = string(exCorrTypeMap.getString(exCorr.exCorrType));
      }
      #ifdef LIBXC_ENABLED
      else if(xcMap_X.getEnum(key.c_str(), exCorr.xcExchange)) {} //Found LibXC Exchange functional
      else if(xcMap_XC.getEnum(key.c_str(), exCorr.xcExcorr)) {} //Found LibXC ExCorr functional
      #endif
      else throw key + " is not a recognized exchange or exchange-correlation functional";
    }
    #ifdef LIBXC_ENABLED
    if((exCorr.xcExchange || exCorr.xcExcorr) && !exCorr.xcCorr)
    {  //LibXC will be used:
      exCorr.exCorrType = ExCorrLibXC;
      if(exCorr.xcExchange)
        pl.get(exCorr.xcCorr, 0, xcMap_C, "funcC", true); //required correlation functional
      //Set the short name:
      if(exCorr.xcExcorr)
        exCorr.xcName = string(xcMap_XC.getString(exCorr.xcExcorr));
      else
        exCorr.xcName = string(xcMap_X.getString(exCorr.xcExchange))
          + ':' + string(xcMap_C.getString(exCorr.xcCorr));
    }
    #endif
  }

  void printStatus(const ExCorr& exCorr)
  {  switch(exCorr.exCorrType)
    {
      #ifdef LIBXC_ENABLED
      case ExCorrLibXC:
      {  if(exCorr.xcExcorr)
          logPrintf("%s", xcMap_XC.getString(exCorr.xcExcorr));
        else
          logPrintf("%s %s", xcMap_X.getString(exCorr.xcExchange), xcMap_C.getString(exCorr.xcCorr));
        break;
      }
      #endif
      default:
        logPrintf("%s", exCorrTypeMap.getString(exCorr.exCorrType));
    }
  }
}
commandElecExCorr;


struct CommandElecExCorrCompare : public CommandElecExCorr
{
  CommandElecExCorrCompare() : CommandElecExCorr("elec-ex-corr-compare")
  {
    format = "<functional>";
    comments =
      "Compute total energies for other functionals at the final state for comparison.\n"
      "The available options for each parameter are identical to elec-ex-corr.\n"
      "\n"
      "This command may be specified multiple times. It invokes 'dump End ExcCompare'\n"
      "automatically, but the compute frequency can be controlled using dump explicitly.";
    hasDefault = false;
    allowMultiple = true;
    emptyParamError = "   eXchange/Correlation functional(s) must be specified.";
    
    #ifdef LIBXC_ENABLED
    format += "\n\t| <funcX> <funcC>\n\t| <funcXC>";
    #endif
    forbid("fix-electron-density");
    forbid("fix-electron-potential");
  }
  
  void process(ParamList& pl, Everything& e)
  {  e.exCorrDiff.push_back(std::shared_ptr<ExCorr>(new ExCorr));
    CommandElecExCorr::process(pl, *e.exCorrDiff.back());
    e.dump.insert(std::make_pair(DumpFreq_End, DumpExcCompare));
  }
  
  void printStatus(Everything& e, int iRep)
  {  CommandElecExCorr::printStatus(*e.exCorrDiff[iRep]);
  }
}
commandElecExCorrCompare;


struct CommandFluidExCorr : public CommandElecExCorr
{
  CommandFluidExCorr() : CommandElecExCorr("fluid-ex-corr", "jdftx/Fluid/Parameters")
  {
    format = "<kinetic> [<exchange-correlation>]";
    comments =
      "Kinetic energy functional for fluid convolution coupling where <kinetic> is one of:"
      + addDescriptions(kineticTypeMap.optionList(), linkDescription(kineticTypeMap, kineticDescriptionMap)) +
      #ifdef LIBXC_ENABLED
      addDescriptions(xcMap_K.optionList(), getLibXCdescription_K) +
      #endif
      ".\n\nThe available options for <exchange-correlation> are identical to elec-ex-corr\n"
      "and defaults to lda-pz.";
    hasDefault = true;
    emptyParamError = "   A kinetic energy functional must be specified.";
    require("elec-ex-corr");
  }
  
  void process(ParamList& pl, Everything& e)
  {  ExCorr& fluidExCorr = e.eVars.fluidParams.exCorr;
    //Get kinetic energy functional:
    string key; pl.get(key, string(), "kinetic");
    if(!key.length()) fluidExCorr.kineticType = KineticTF; //default: Thomas-Fermi
    else
    {  if(kineticTypeMap.getEnum(key.c_str(), fluidExCorr.kineticType)) {} //Found internal kinetic functional
      #ifdef LIBXC_ENABLED
      else if(xcMap_K.getEnum(key.c_str(), fluidExCorr.xcKinetic)) { fluidExCorr.kineticType = KineticLibXC; } //Found LibXC kinetic functional
      #endif
      else throw key + " is not a recognized kinetic energy functional";
    }
    //Set default exchange-correlation to be LDA
    fluidExCorr.exCorrType =  ExCorrLDA_PZ;
    fluidExCorr.xcName = exCorrTypeMap.getString(ExCorrLDA_PZ);
      
    CommandElecExCorr::process(pl, fluidExCorr);
  }
  
  void printStatus(Everything& e, int iRep)
  {  const ExCorr& fluidExCorr = e.eVars.fluidParams.exCorr;
    switch(fluidExCorr.exCorrType)
    {
      #ifdef LIBXC_ENABLED
      case KineticLibXC: logPrintf("%s ", xcMap_K.getString(fluidExCorr.xcKinetic)); break;
      #endif
      default: logPrintf("%s ", kineticTypeMap.getString(fluidExCorr.kineticType));
    }
    CommandElecExCorr::printStatus(fluidExCorr);
  }
}
commandFluidExCorr;


struct CommandVanDerWaals : public Command
{
  CommandVanDerWaals() : Command("van-der-waals", "jdftx/Electronic/Functional")
  {
    format = "[<scaleOverride>=0]";
    comments =
      "DFT+D2 pair-potential corrections for the long range Van der Waals\n"
      "interaction [S. Grimme, J. Comput. Chem. 27: 1787–1799 (2006)].\n"
      "\n"
      "Default scale factors are available for the gga-PBE, hyb-gga-xc-b3lyp\n"
      "and mgga-TPSS exchange-correlation functionals (see elec-ex-corr).\n"
      "Manually specify <scaleOverride> to use with other functionals";
  }

  void process(ParamList& pl, Everything& e)
  {  e.iInfo.vdWenable = true;
    pl.get(e.iInfo.vdWscale, 0., "scaleOverride");
  }

  void printStatus(Everything& e, int iRep)
  {  if(e.iInfo.vdWscale) logPrintf("%lg", e.iInfo.vdWscale);
  }
}
commandVanDerWaals;


struct CommandExchangeParameters : public Command
{
  CommandExchangeParameters() : Command("exchange-parameters", "jdftx/Electronic/Functional")
  {
    format = "<exxScale> [<exxOmega>=0]";
    comments =
      "Override exact-exchange parameters in a hybrid functional.\n"
      "Here <exxScale> is the scale fraction of exact exchange,\n"
      "and <exxOmega> is the screening parameter.\n"
      "This is only supported for internal hybrid functionals PBE0\n"
      "and HSExx, and not for the LibXC hybrid functionals.";
  }

  void process(ParamList& pl, Everything& e)
  {  pl.get(e.exCorr.exxScaleOverride, 0., "exxScale", true);
    pl.get(e.exCorr.exxOmegaOverride, 0., "exxOmega");
    if(e.exCorr.exxScaleOverride <= 0.) throw string("<exxScale> must be >= 0");
  }

  void printStatus(Everything& e, int iRep)
  {  logPrintf(" %lg", e.exCorr.exxScaleOverride);
    if(e.exCorr.exxOmegaOverride) logPrintf(" %lg", e.exCorr.exxOmegaOverride);
  }
}
commandExchangeParameters;
