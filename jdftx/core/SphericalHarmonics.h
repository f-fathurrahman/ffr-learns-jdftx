/*-------------------------------------------------------------------
Copyright 2012 Ravishankar Sundararaman

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

#ifndef JDFTX_CORE_SPHERICALHARMONICS_H
#define JDFTX_CORE_SPHERICALHARMONICS_H

#include <core/vector3.h>

//! @addtogroup Utilities
//! @{
//! @file SphericalHarmonics.h Real spherical Harmonics and spherical bessel functions

//! Spherical bessel function
inline double bessel_jl(int l, double x)
{	if(fabs(x) > 1.+0.1*l)
	{	double s, c; sincos(x, &s, &c);
		double xInv=1./x, xInvSq = xInv * xInv;
		switch(l)
		{	case 0: return xInv * s;
			case 1: return xInv * (xInv*s - c);
			case 2: return xInv * ((3*xInvSq-1)*s - 3*xInv*c);
			case 3: return xInv * ((15*xInvSq-6)*xInv*s + (1-15*xInvSq)*c);
			case 4: return xInv * ((1+xInvSq*(-45+xInvSq*105))*s + xInv*(10-105*xInvSq)*c);
			case 5: return xInv * (xInv*(15+xInvSq*(-420+xInvSq*945))*s + (-1+xInvSq*(105-945*xInvSq))*c);
			case 6: return xInv * ((-1+xInvSq*(210+xInvSq*(-4725 + 10395*xInvSq)))*s + xInv*(-21+xInvSq*(1260 - 10395*xInvSq))*c);
			default: return 0.; //unsupported l
		}
	}
	else //Series expansions about 0 to prevent roundoff errors (accurate to 1 in 1e15 at the crossover for each l)
	{	double term = 1.;
		for(int i=3; i<=2*l+1; i+=2)
			term *= x/i;
		double ret = term;
		double xSq = x*x;
		for(int i=2; i<=14; i+=2)
		{	term *= -xSq/(i*(i+2*l+1));
			ret += term;
		}
		return ret;
	}
}


//! Auto-generated code for computing real spherical harmonics and their Clebsch-Gordon coefficients
//! The external interface to these functions follows the end of this namespace block
namespace YlmInternal
{
	template<int lm> __hostanddev__ double Ylm(double x, double y, double z); //flat-indexed by lm := l*(l+1) + m
	#define Power pow //For code auto-generated by Mathematica
	#define DECLARE_Ylm(lm,code) \
		template<> __hostanddev__ double Ylm<lm>(double x, double y, double z) { return code; }
	DECLARE_Ylm(0, 0.28209479177387814)
	DECLARE_Ylm(1, 0.4886025119029199*y)
	DECLARE_Ylm(2, 0.4886025119029199*z)
	DECLARE_Ylm(3, 0.4886025119029199*x)
	DECLARE_Ylm(4, 1.0925484305920792*x*y)
	DECLARE_Ylm(5, 1.0925484305920792*y*z)
	DECLARE_Ylm(6, -0.31539156525252005*(Power(x,2) + Power(y,2) - 2.*Power(z,2)))
	DECLARE_Ylm(7, 1.0925484305920792*x*z)
	DECLARE_Ylm(8, 0.5462742152960396*(x - y)*(x + y))
	DECLARE_Ylm(9, -0.5900435899266435*y*(-3.*Power(x,2) + Power(y,2)))
	DECLARE_Ylm(10, 2.890611442640554*x*y*z)
	DECLARE_Ylm(11, -0.4570457994644658*y*(Power(x,2) + Power(y,2) - 4.*Power(z,2)))
	DECLARE_Ylm(12, 0.3731763325901154*z*(-3.*(Power(x,2) + Power(y,2)) + 2.*Power(z,2)))
	DECLARE_Ylm(13, -0.4570457994644658*x*(Power(x,2) + Power(y,2) - 4.*Power(z,2)))
	DECLARE_Ylm(14, 1.445305721320277*(x - y)*(x + y)*z)
	DECLARE_Ylm(15, 0.5900435899266435*x*(Power(x,2) - 3.*Power(y,2)))
	DECLARE_Ylm(16, 2.5033429417967046*x*(x - y)*y*(x + y))
	DECLARE_Ylm(17, -1.7701307697799304*y*(-3.*Power(x,2) + Power(y,2))*z)
	DECLARE_Ylm(18, -0.9461746957575601*x*y*(Power(x,2) + Power(y,2) - 6.*Power(z,2)))
	DECLARE_Ylm(19, -0.6690465435572892*y*z*(3.*(Power(x,2) + Power(y,2)) - 4.*Power(z,2)))
	DECLARE_Ylm(20, 0.03526184897173477*(9.*Power(Power(x,2) + Power(y,2),2) - 72.*(Power(x,2) + Power(y,2))*Power(z,2) + 24.*Power(z,4)))
	DECLARE_Ylm(21, -0.6690465435572892*x*z*(3.*(Power(x,2) + Power(y,2)) - 4.*Power(z,2)))
	DECLARE_Ylm(22, -0.47308734787878004*(x - y)*(x + y)*(Power(x,2) + Power(y,2) - 6.*Power(z,2)))
	DECLARE_Ylm(23, 1.7701307697799304*x*(Power(x,2) - 3.*Power(y,2))*z)
	DECLARE_Ylm(24, 0.6258357354491761*(Power(x,4) - 6.*Power(x,2)*Power(y,2) + Power(y,4)))
	DECLARE_Ylm(25, 0.6563820568401701*y*(5.*Power(x,4) - 10.*Power(x,2)*Power(y,2) + Power(y,4)))
	DECLARE_Ylm(26, 8.302649259524166*x*(x - y)*y*(x + y)*z)
	DECLARE_Ylm(27, 0.4892382994352504*y*(-3.*Power(x,2) + Power(y,2))*(Power(x,2) + Power(y,2) - 8.*Power(z,2)))
	DECLARE_Ylm(28, -4.793536784973324*x*y*z*(Power(x,2) + Power(y,2) - 2.*Power(z,2)))
	DECLARE_Ylm(29, 0.45294665119569694*y*(Power(Power(x,2) + Power(y,2),2) - 12.*(Power(x,2) + Power(y,2))*Power(z,2) + 8.*Power(z,4)))
	DECLARE_Ylm(30, 0.1169503224534236*z*(15.*Power(Power(x,2) + Power(y,2),2) - 40.*(Power(x,2) + Power(y,2))*Power(z,2) + 8.*Power(z,4)))
	DECLARE_Ylm(31, 0.45294665119569694*x*(Power(Power(x,2) + Power(y,2),2) - 12.*(Power(x,2) + Power(y,2))*Power(z,2) + 8.*Power(z,4)))
	DECLARE_Ylm(32, -2.396768392486662*(x - y)*(x + y)*z*(Power(x,2) + Power(y,2) - 2.*Power(z,2)))
	DECLARE_Ylm(33, -0.4892382994352504*x*(Power(x,2) - 3.*Power(y,2))*(Power(x,2) + Power(y,2) - 8.*Power(z,2)))
	DECLARE_Ylm(34, 2.0756623148810416*(Power(x,4) - 6.*Power(x,2)*Power(y,2) + Power(y,4))*z)
	DECLARE_Ylm(35, 0.6563820568401701*x*(Power(x,4) - 10.*Power(x,2)*Power(y,2) + 5.*Power(y,4)))
	DECLARE_Ylm(36, 1.3663682103838286*x*y*(3.*Power(x,4) - 10.*Power(x,2)*Power(y,2) + 3.*Power(y,4)))
	DECLARE_Ylm(37, 2.366619162231752*y*(5.*Power(x,4) - 10.*Power(x,2)*Power(y,2) + Power(y,4))*z)
	DECLARE_Ylm(38, -2.0182596029148967*x*(x - y)*y*(x + y)*(Power(x,2) + Power(y,2) - 10.*Power(z,2)))
	DECLARE_Ylm(39, 0.9212052595149236*y*(-3.*Power(x,2) + Power(y,2))*z*(3.*(Power(x,2) + Power(y,2)) - 8.*Power(z,2)))
	DECLARE_Ylm(40, 0.9212052595149236*x*y*(Power(Power(x,2) + Power(y,2),2) - 16.*(Power(x,2) + Power(y,2))*Power(z,2) + 16.*Power(z,4)))
	DECLARE_Ylm(41, 0.5826213625187314*y*z*(5.*Power(Power(x,2) + Power(y,2),2) - 20.*(Power(x,2) + Power(y,2))*Power(z,2) + 8.*Power(z,4)))
	DECLARE_Ylm(42, 0.06356920226762842*(-5.*Power(Power(x,2) + Power(y,2),3) + 90.*Power(Power(x,2) + Power(y,2),2)*Power(z,2) - 120.*(Power(x,2) + Power(y,2))*Power(z,4) + 16.*Power(z,6)))
	DECLARE_Ylm(43, 0.5826213625187314*x*z*(5.*Power(Power(x,2) + Power(y,2),2) - 20.*(Power(x,2) + Power(y,2))*Power(z,2) + 8.*Power(z,4)))
	DECLARE_Ylm(44, 0.4606026297574618*(x - y)*(x + y)*(Power(Power(x,2) + Power(y,2),2) - 16.*(Power(x,2) + Power(y,2))*Power(z,2) + 16.*Power(z,4)))
	DECLARE_Ylm(45, -0.9212052595149236*x*(Power(x,2) - 3.*Power(y,2))*z*(3.*(Power(x,2) + Power(y,2)) - 8.*Power(z,2)))
	DECLARE_Ylm(46, -0.5045649007287242*(Power(x,4) - 6.*Power(x,2)*Power(y,2) + Power(y,4))*(Power(x,2) + Power(y,2) - 10.*Power(z,2)))
	DECLARE_Ylm(47, 2.366619162231752*x*(Power(x,4) - 10.*Power(x,2)*Power(y,2) + 5.*Power(y,4))*z)
	DECLARE_Ylm(48, 0.6831841051919143*(Power(x,6) - 15.*Power(x,4)*Power(y,2) + 15.*Power(x,2)*Power(y,4) - Power(y,6)))
	#undef DECLARE_Ylm
	#undef Power
}

//! Index by combined lm := l*(l+1)+m index (useful when static-looping over all l,m)
template<int lm> __hostanddev__ double Ylm(const vector3<>& qhat)
{	return YlmInternal::Ylm<lm>(qhat[0],qhat[1],qhat[2]);
}

//! Index by l and m separately
template<int l, int m> __hostanddev__ double Ylm(const vector3<>& qhat)
{	return Ylm<l*(l+1)+m>(qhat);
}

//! Switch a function templated over l,m for all supported l,m with parenthesis enclosed argument list argList
#define SwitchTemplate_lm(l,m,fTemplate,argList) \
	switch(l*(l+1)+m) \
	{	case 0: fTemplate<0,0> argList; break; \
		case 1: fTemplate<1,-1> argList; break; \
		case 2: fTemplate<1,0> argList; break; \
		case 3: fTemplate<1,1> argList; break; \
		case 4: fTemplate<2,-2> argList; break; \
		case 5: fTemplate<2,-1> argList; break; \
		case 6: fTemplate<2,0> argList; break; \
		case 7: fTemplate<2,1> argList; break; \
		case 8: fTemplate<2,2> argList; break; \
		case 9: fTemplate<3,-3> argList; break; \
		case 10: fTemplate<3,-2> argList; break; \
		case 11: fTemplate<3,-1> argList; break; \
		case 12: fTemplate<3,0> argList; break; \
		case 13: fTemplate<3,1> argList; break; \
		case 14: fTemplate<3,2> argList; break; \
		case 15: fTemplate<3,3> argList; break; \
		case 16: fTemplate<4,-4> argList; break; \
		case 17: fTemplate<4,-3> argList; break; \
		case 18: fTemplate<4,-2> argList; break; \
		case 19: fTemplate<4,-1> argList; break; \
		case 20: fTemplate<4,0> argList; break; \
		case 21: fTemplate<4,1> argList; break; \
		case 22: fTemplate<4,2> argList; break; \
		case 23: fTemplate<4,3> argList; break; \
		case 24: fTemplate<4,4> argList; break; \
		case 25: fTemplate<5,-5> argList; break; \
		case 26: fTemplate<5,-4> argList; break; \
		case 27: fTemplate<5,-3> argList; break; \
		case 28: fTemplate<5,-2> argList; break; \
		case 29: fTemplate<5,-1> argList; break; \
		case 30: fTemplate<5,0> argList; break; \
		case 31: fTemplate<5,1> argList; break; \
		case 32: fTemplate<5,2> argList; break; \
		case 33: fTemplate<5,3> argList; break; \
		case 34: fTemplate<5,4> argList; break; \
		case 35: fTemplate<5,5> argList; break; \
		case 36: fTemplate<6,-6> argList; break; \
		case 37: fTemplate<6,-5> argList; break; \
		case 38: fTemplate<6,-4> argList; break; \
		case 39: fTemplate<6,-3> argList; break; \
		case 40: fTemplate<6,-2> argList; break; \
		case 41: fTemplate<6,-1> argList; break; \
		case 42: fTemplate<6,0> argList; break; \
		case 43: fTemplate<6,1> argList; break; \
		case 44: fTemplate<6,2> argList; break; \
		case 45: fTemplate<6,3> argList; break; \
		case 46: fTemplate<6,4> argList; break; \
		case 47: fTemplate<6,5> argList; break; \
		case 48: fTemplate<6,6> argList; break; \
	}

//! Use above macro to provide a non-templated version of the function
template<int l, int m> void set_Ylm(const vector3<> qHat, double& result) { result = Ylm<l,m>(qHat); }
inline double Ylm(int l, int m, const vector3<>& qHat) { double result=0.;  SwitchTemplate_lm(l,m, set_Ylm, (qHat, result)); return result; }

//! Term in real spherical harmonic expansion of a product of two real spherical harmonics
struct YlmProdTerm
{	int l, m; //!< angular quantum numbers of current term
	double coeff; //!< coefficient of Ylm with current l and m
	YlmProdTerm(int l, int m, double coeff) : l(l), m(m), coeff(coeff) {}
};

//! Real spherical harmonic expansion of product of two real spherical harmonics
//! (effectively returns list of non-zero Clebsch-Gordon coefficients in the modified real Ylm basis)
inline std::vector<YlmProdTerm> expandYlmProd(int lm1, int lm2)
{	if(lm2 > lm1) std::swap(lm1, lm2);
	std::vector<YlmProdTerm> result;
	#define ADD(l,m,coeff) result.push_back(YlmProdTerm(l,m,coeff)) //shorthand for use in Mathematica generated list below:
	switch(lm2 + (lm1*(lm1+1))/2)
	{	case 0: ADD(0,0,0.28209479177387814); break;
		case 1: ADD(1,-1,0.28209479177387814); break;
		case 2: ADD(0,0,0.28209479177387814); ADD(2,0,-0.126156626101008); ADD(2,2,-0.2185096861184158); break;
		case 3: ADD(1,0,0.28209479177387814); break;
		case 4: ADD(2,-1,0.2185096861184158); break;
		case 5: ADD(0,0,0.28209479177387814); ADD(2,0,0.252313252202016); break;
		case 6: ADD(1,1,0.28209479177387814); break;
		case 7: ADD(2,-2,0.2185096861184158); break;
		case 8: ADD(2,1,0.2185096861184158); break;
		case 9: ADD(0,0,0.28209479177387814); ADD(2,0,-0.126156626101008); ADD(2,2,0.2185096861184158); break;
		case 10: ADD(2,-2,0.28209479177387814); break;
		case 11: ADD(1,1,0.2185096861184158); ADD(3,1,-0.058399170081901854); ADD(3,3,-0.2261790131595403); break;
		case 12: ADD(3,-2,0.18467439092237178); break;
		case 13: ADD(1,-1,0.2185096861184158); ADD(3,-3,0.2261790131595403); ADD(3,-1,-0.058399170081901854); break;
		case 14: ADD(0,0,0.28209479177387814); ADD(2,0,-0.18022375157286857); ADD(4,0,0.04029925596769687); ADD(4,4,-0.23841361350444806); break;
		case 15: ADD(2,-1,0.28209479177387814); break;
		case 16: ADD(1,0,0.2185096861184158); ADD(3,0,-0.14304816810266882); ADD(3,2,-0.18467439092237178); break;
		case 17: ADD(1,-1,0.2185096861184158); ADD(3,-1,0.23359668032760741); break;
		case 18: ADD(3,-2,0.18467439092237178); break;
		case 19: ADD(2,1,0.15607834722743988); ADD(4,1,-0.06371871843402754); ADD(4,3,-0.16858388283618386); break;
		case 20: ADD(0,0,0.28209479177387814); ADD(2,0,0.09011187578643429); ADD(2,2,-0.15607834722743988); ADD(4,0,-0.1611970238707875); ADD(4,2,-0.18022375157286857); break;
		case 21: ADD(2,0,0.28209479177387814); break;
		case 22: ADD(1,-1,-0.126156626101008); ADD(3,-1,0.20230065940342062); break;
		case 23: ADD(1,0,0.252313252202016); ADD(3,0,0.2477666950834761); break;
		case 24: ADD(1,1,-0.126156626101008); ADD(3,1,0.20230065940342062); break;
		case 25: ADD(2,-2,-0.18022375157286857); ADD(4,-2,0.15607834722743988); break;
		case 26: ADD(2,-1,0.09011187578643429); ADD(4,-1,0.2207281154418226); break;
		case 27: ADD(0,0,0.28209479177387814); ADD(2,0,0.18022375157286857); ADD(4,0,0.24179553580618124); break;
		case 28: ADD(2,1,0.28209479177387814); break;
		case 29: ADD(3,-2,0.18467439092237178); break;
		case 30: ADD(1,1,0.2185096861184158); ADD(3,1,0.23359668032760741); break;
		case 31: ADD(1,0,0.2185096861184158); ADD(3,0,-0.14304816810266882); ADD(3,2,0.18467439092237178); break;
		case 32: ADD(2,-1,0.15607834722743988); ADD(4,-3,0.16858388283618386); ADD(4,-1,-0.06371871843402754); break;
		case 33: ADD(2,-2,0.15607834722743988); ADD(4,-2,0.18022375157286857); break;
		case 34: ADD(2,1,0.09011187578643429); ADD(4,1,0.2207281154418226); break;
		case 35: ADD(0,0,0.28209479177387814); ADD(2,0,0.09011187578643429); ADD(2,2,0.15607834722743988); ADD(4,0,-0.1611970238707875); ADD(4,2,0.18022375157286857); break;
		case 36: ADD(2,2,0.28209479177387814); break;
		case 37: ADD(1,-1,-0.2185096861184158); ADD(3,-3,0.2261790131595403); ADD(3,-1,0.058399170081901854); break;
		case 38: ADD(3,2,0.18467439092237178); break;
		case 39: ADD(1,1,0.2185096861184158); ADD(3,1,-0.058399170081901854); ADD(3,3,0.2261790131595403); break;
		case 40: ADD(4,-4,0.23841361350444806); break;
		case 41: ADD(2,-1,-0.15607834722743988); ADD(4,-3,0.16858388283618386); ADD(4,-1,0.06371871843402754); break;
		case 42: ADD(2,2,-0.18022375157286857); ADD(4,2,0.15607834722743988); break;
		case 43: ADD(2,1,0.15607834722743988); ADD(4,1,-0.06371871843402754); ADD(4,3,0.16858388283618386); break;
		case 44: ADD(0,0,0.28209479177387814); ADD(2,0,-0.18022375157286857); ADD(4,0,0.04029925596769687); ADD(4,4,0.23841361350444806); break;
		case 45: ADD(3,-3,0.28209479177387814); break;
		case 46: ADD(2,2,0.2261790131595403); ADD(4,2,-0.04352817137756816); ADD(4,4,-0.23032943298089034); break;
		case 47: ADD(4,-3,0.16286750396763996); break;
		case 48: ADD(2,-2,0.2261790131595403); ADD(4,-4,0.23032943298089034); ADD(4,-2,-0.04352817137756816); break;
		case 49: ADD(1,1,0.2261790131595403); ADD(3,1,-0.09403159725795937); ADD(5,1,0.01694331772935932); ADD(5,5,-0.2455320005465369); break;
		case 50: ADD(3,2,0.1486770096793976); ADD(5,2,-0.04482780509623635); ADD(5,4,-0.1552880720369528); break;
		case 51: ADD(3,-3,-0.21026104350168); ADD(5,-3,0.12679217987703037); break;
		case 52: ADD(3,-2,0.1486770096793976); ADD(5,-4,0.1552880720369528); ADD(5,-2,-0.04482780509623635); break;
		case 53: ADD(1,-1,0.2261790131595403); ADD(3,-1,-0.09403159725795937); ADD(5,-5,0.2455320005465369); ADD(5,-1,0.01694331772935932); break;
		case 54: ADD(0,0,0.28209479177387814); ADD(2,0,-0.21026104350168); ADD(4,0,0.07693494321105766); ADD(6,0,-0.011854396693264043); ADD(6,6,-0.25480059867297505); break;
		case 55: ADD(3,-2,0.28209479177387814); break;
		case 56: ADD(2,1,0.18467439092237178); ADD(4,1,-0.07539300438651343); ADD(4,3,-0.19947114020071635); break;
		case 57: ADD(2,-2,0.18467439092237178); ADD(4,-2,0.21324361862292307); break;
		case 58: ADD(2,-1,0.18467439092237178); ADD(4,-3,0.19947114020071635); ADD(4,-1,-0.07539300438651343); break;
		case 59: ADD(1,0,0.18467439092237178); ADD(3,0,-0.18806319451591874); ADD(5,0,0.05357947514468781); ADD(5,4,-0.19018826981554557); break;
		case 60: ADD(1,1,0.18467439092237178); ADD(3,1,0.11516471649044517); ADD(3,3,-0.1486770096793976); ADD(5,1,-0.08300496597356405); ADD(5,3,-0.1793112203849454); break;
		case 61: ADD(5,-2,0.19018826981554557); break;
		case 62: ADD(1,-1,0.18467439092237178); ADD(3,-3,0.1486770096793976); ADD(3,-1,0.11516471649044517); ADD(5,-3,0.1793112203849454); ADD(5,-1,-0.08300496597356405); break;
		case 63: ADD(5,-4,0.19018826981554557); break;
		case 64: ADD(2,1,0.1486770096793976); ADD(4,1,-0.09932258459927992); ADD(6,1,0.022177545476549994); ADD(6,5,-0.1801712311720527); break;
		case 65: ADD(0,0,0.28209479177387814); ADD(4,0,-0.1795148674924679); ADD(4,4,-0.15171775404828514); ADD(6,0,0.07112638015958425); ADD(6,4,-0.18818271355849853); break;
		case 66: ADD(3,-1,0.28209479177387814); break;
		case 67: ADD(2,0,0.20230065940342062); ADD(2,2,0.058399170081901854); ADD(4,0,-0.15078600877302686); ADD(4,2,-0.16858388283618386); break;
		case 68: ADD(2,-1,0.23359668032760741); ADD(4,-1,0.23841361350444806); break;
		case 69: ADD(2,-2,-0.058399170081901854); ADD(4,-2,0.16858388283618386); break;
		case 70: ADD(1,1,-0.058399170081901854); ADD(3,1,0.1456731240789439); ADD(3,3,0.09403159725795937); ADD(5,1,-0.0656211873953095); ADD(5,3,-0.14175796661021045); break;
		case 71: ADD(1,0,0.23359668032760741); ADD(3,0,0.05947080387175903); ADD(3,2,-0.11516471649044517); ADD(5,0,-0.1694331772935932); ADD(5,2,-0.17361734258475534); break;
		case 72: ADD(1,-1,0.20230065940342062); ADD(3,-1,0.126156626101008); ADD(5,-1,0.22731846124334898); break;
		case 73: ADD(3,-2,0.11516471649044517); ADD(5,-2,0.17361734258475534); break;
		case 74: ADD(1,-1,0.058399170081901854); ADD(3,-3,-0.09403159725795937); ADD(3,-1,-0.1456731240789439); ADD(5,-3,0.14175796661021045); ADD(5,-1,0.0656211873953095); break;
		case 75: ADD(2,2,-0.09403159725795937); ADD(4,2,0.13325523051897814); ADD(4,4,0.11752006695060024); ADD(6,2,-0.04435509095309999); ADD(6,4,-0.1214714192760309); break;
		case 76: ADD(2,1,0.11516471649044517); ADD(4,1,0.10257992428141023); ADD(4,3,-0.06785024228911189); ADD(6,1,-0.08589326429043577); ADD(6,3,-0.16297101049475005); break;
		case 77: ADD(0,0,0.28209479177387814); ADD(2,0,0.126156626101008); ADD(2,2,-0.1456731240789439); ADD(4,0,0.025644981070352558); ADD(4,2,-0.11468784191000729); ADD(6,0,-0.17781595039896067); ADD(6,2,-0.17178652858087154); break;
		case 78: ADD(3,0,0.28209479177387814); break;
		case 79: ADD(2,-1,-0.14304816810266882); ADD(4,-1,0.19466390027300617); break;
		case 80: ADD(2,0,0.2477666950834761); ADD(4,0,0.24623252122982908); break;
		case 81: ADD(2,1,-0.14304816810266882); ADD(4,1,0.19466390027300617); break;
		case 82: ADD(3,-2,-0.18806319451591874); ADD(5,-2,0.14175796661021045); break;
		case 83: ADD(1,-1,-0.14304816810266882); ADD(3,-1,0.05947080387175903); ADD(5,-1,0.21431790057875125); break;
		case 84: ADD(1,0,0.2477666950834761); ADD(3,0,0.168208834801344); ADD(5,0,0.23961469724456466); break;
		case 85: ADD(1,1,-0.14304816810266882); ADD(3,1,0.05947080387175903); ADD(5,1,0.21431790057875125); break;
		case 86: ADD(3,2,-0.18806319451591874); ADD(5,2,0.14175796661021045); break;
		case 87: ADD(4,-3,-0.20355072686733566); ADD(6,-3,0.10864734032983336); break;
		case 88: ADD(2,-2,-0.18806319451591874); ADD(4,-2,-0.04441841017299272); ADD(6,-2,0.17742036381239995); break;
		case 89: ADD(2,-1,0.05947080387175903); ADD(4,-1,0.09932258459927992); ADD(6,-1,0.22177545476549995); break;
		case 90: ADD(0,0,0.28209479177387814); ADD(2,0,0.168208834801344); ADD(4,0,0.15386988642211533); ADD(6,0,0.23708793386528085); break;
		case 91: ADD(3,1,0.28209479177387814); break;
		case 92: ADD(2,-2,-0.058399170081901854); ADD(4,-2,0.16858388283618386); break;
		case 93: ADD(2,1,0.23359668032760741); ADD(4,1,0.23841361350444806); break;
		case 94: ADD(2,0,0.20230065940342062); ADD(2,2,-0.058399170081901854); ADD(4,0,-0.15078600877302686); ADD(4,2,0.16858388283618386); break;
		case 95: ADD(1,-1,-0.058399170081901854); ADD(3,-3,-0.09403159725795937); ADD(3,-1,0.1456731240789439); ADD(5,-3,0.14175796661021045); ADD(5,-1,-0.0656211873953095); break;
		case 96: ADD(3,-2,0.11516471649044517); ADD(5,-2,0.17361734258475534); break;
		case 97: ADD(1,1,0.20230065940342062); ADD(3,1,0.126156626101008); ADD(5,1,0.22731846124334898); break;
		case 98: ADD(1,0,0.23359668032760741); ADD(3,0,0.05947080387175903); ADD(3,2,0.11516471649044517); ADD(5,0,-0.1694331772935932); ADD(5,2,0.17361734258475534); break;
		case 99: ADD(1,1,-0.058399170081901854); ADD(3,1,0.1456731240789439); ADD(3,3,-0.09403159725795937); ADD(5,1,-0.0656211873953095); ADD(5,3,0.14175796661021045); break;
		case 100: ADD(2,-2,-0.09403159725795937); ADD(4,-4,-0.11752006695060024); ADD(4,-2,0.13325523051897814); ADD(6,-4,0.1214714192760309); ADD(6,-2,-0.04435509095309999); break;
		case 101: ADD(2,-1,0.11516471649044517); ADD(4,-3,0.06785024228911189); ADD(4,-1,0.10257992428141023); ADD(6,-3,0.16297101049475005); ADD(6,-1,-0.08589326429043577); break;
		case 102: ADD(2,-2,0.1456731240789439); ADD(4,-2,0.11468784191000729); ADD(6,-2,0.17178652858087154); break;
		case 103: ADD(2,1,0.05947080387175903); ADD(4,1,0.09932258459927992); ADD(6,1,0.22177545476549995); break;
		case 104: ADD(0,0,0.28209479177387814); ADD(2,0,0.126156626101008); ADD(2,2,0.1456731240789439); ADD(4,0,0.025644981070352558); ADD(4,2,0.11468784191000729); ADD(6,0,-0.17781595039896067); ADD(6,2,0.17178652858087154); break;
		case 105: ADD(3,2,0.28209479177387814); break;
		case 106: ADD(2,-1,-0.18467439092237178); ADD(4,-3,0.19947114020071635); ADD(4,-1,0.07539300438651343); break;
		case 107: ADD(2,2,0.18467439092237178); ADD(4,2,0.21324361862292307); break;
		case 108: ADD(2,1,0.18467439092237178); ADD(4,1,-0.07539300438651343); ADD(4,3,0.19947114020071635); break;
		case 109: ADD(5,-4,0.19018826981554557); break;
		case 110: ADD(1,-1,-0.18467439092237178); ADD(3,-3,0.1486770096793976); ADD(3,-1,-0.11516471649044517); ADD(5,-3,0.1793112203849454); ADD(5,-1,0.08300496597356405); break;
		case 111: ADD(5,2,0.19018826981554557); break;
		case 112: ADD(1,1,0.18467439092237178); ADD(3,1,0.11516471649044517); ADD(3,3,0.1486770096793976); ADD(5,1,-0.08300496597356405); ADD(5,3,0.1793112203849454); break;
		case 113: ADD(1,0,0.18467439092237178); ADD(3,0,-0.18806319451591874); ADD(5,0,0.05357947514468781); ADD(5,4,0.19018826981554557); break;
		case 114: ADD(2,-1,0.1486770096793976); ADD(4,-1,-0.09932258459927992); ADD(6,-5,0.1801712311720527); ADD(6,-1,0.022177545476549994); break;
		case 115: ADD(4,-4,0.15171775404828514); ADD(6,-4,0.18818271355849853); break;
		case 116: ADD(2,-1,-0.11516471649044517); ADD(4,-3,0.06785024228911189); ADD(4,-1,-0.10257992428141023); ADD(6,-3,0.16297101049475005); ADD(6,-1,0.08589326429043577); break;
		case 117: ADD(2,2,-0.18806319451591874); ADD(4,2,-0.04441841017299272); ADD(6,2,0.17742036381239995); break;
		case 118: ADD(2,1,0.11516471649044517); ADD(4,1,0.10257992428141023); ADD(4,3,0.06785024228911189); ADD(6,1,-0.08589326429043577); ADD(6,3,0.16297101049475005); break;
		case 119: ADD(0,0,0.28209479177387814); ADD(4,0,-0.1795148674924679); ADD(4,4,0.15171775404828514); ADD(6,0,0.07112638015958425); ADD(6,4,0.18818271355849853); break;
		case 120: ADD(3,3,0.28209479177387814); break;
		case 121: ADD(2,-2,-0.2261790131595403); ADD(4,-4,0.23032943298089034); ADD(4,-2,0.04352817137756816); break;
		case 122: ADD(4,3,0.16286750396763996); break;
		case 123: ADD(2,2,0.2261790131595403); ADD(4,2,-0.04352817137756816); ADD(4,4,0.23032943298089034); break;
		case 124: ADD(1,-1,-0.2261790131595403); ADD(3,-1,0.09403159725795937); ADD(5,-5,0.2455320005465369); ADD(5,-1,-0.01694331772935932); break;
		case 125: ADD(3,-2,-0.1486770096793976); ADD(5,-4,0.1552880720369528); ADD(5,-2,0.04482780509623635); break;
		case 126: ADD(3,3,-0.21026104350168); ADD(5,3,0.12679217987703037); break;
		case 127: ADD(3,2,0.1486770096793976); ADD(5,2,-0.04482780509623635); ADD(5,4,0.1552880720369528); break;
		case 128: ADD(1,1,0.2261790131595403); ADD(3,1,-0.09403159725795937); ADD(5,1,0.01694331772935932); ADD(5,5,0.2455320005465369); break;
		case 129: ADD(6,-6,0.25480059867297505); break;
		case 130: ADD(2,-1,-0.1486770096793976); ADD(4,-1,0.09932258459927992); ADD(6,-5,0.1801712311720527); ADD(6,-1,-0.022177545476549994); break;
		case 131: ADD(2,-2,0.09403159725795937); ADD(4,-4,-0.11752006695060024); ADD(4,-2,-0.13325523051897814); ADD(6,-4,0.1214714192760309); ADD(6,-2,0.04435509095309999); break;
		case 132: ADD(4,3,-0.20355072686733566); ADD(6,3,0.10864734032983336); break;
		case 133: ADD(2,2,-0.09403159725795937); ADD(4,2,0.13325523051897814); ADD(4,4,-0.11752006695060024); ADD(6,2,-0.04435509095309999); ADD(6,4,0.1214714192760309); break;
		case 134: ADD(2,1,0.1486770096793976); ADD(4,1,-0.09932258459927992); ADD(6,1,0.022177545476549994); ADD(6,5,0.1801712311720527); break;
		case 135: ADD(0,0,0.28209479177387814); ADD(2,0,-0.21026104350168); ADD(4,0,0.07693494321105766); ADD(6,0,-0.011854396693264043); ADD(6,6,0.25480059867297505); break;
	}
	#undef ADD
	return result;
}
//! Wrapper function expandYlmProd with individual indices
inline std::vector<YlmProdTerm> expandYlmProd(int l1, int m1, int l2, int m2)
{	int lm1 = l1*(l1+1) + m1;
	int lm2 = l2*(l2+1) + m2;
	return expandYlmProd(lm1, lm2);
}

//! @}
#endif // JDFTX_CORE_SPHERICALHARMONICS_H

