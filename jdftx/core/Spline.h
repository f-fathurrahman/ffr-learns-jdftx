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

#ifndef JDFTX_CORE_SPLINE_H
#define JDFTX_CORE_SPLINE_H

//! @addtogroup Utilities
//! @{

//! @file Spline.h Spline interpolation routines

#include <core/scalar.h>
#include <cstdio>

//! C4-continuous interpolation using quintic splines
namespace QuinticSpline
{
	//! @brief Generate quintic spline coefficients for a set of uniformly-spaced samples.
	//! @param samples An array of uniformly-spaced samples
	//! @param oddExtension The boundary condition at the first sample is odd/even for oddExtension=true/false respectively.
	//! Natural boundary conditions (third and fourth derivatives zero) are imposed on the last sample.
	std::vector<double> getCoeff(const std::vector<double>& samples, bool oddExtension=false);
	
	//! @brief Compute value of quintic spline.
	//! Warning: x is not range-checked
	//! @param coeff pointer to coefficient array generated by getCoeff
	//! @param x location to evaluate spline in the continuous range [0, nCoeff-1)
	__hostanddev__ double value(const double* coeff, double x);
	
	//! @brief Compute derivative (w.r.t x) of quintic spline.
	//! Warning: x is not range-checked
	//! @param coeff pointer to coefficient array generated by getCoeff
	//! @param x location to evaluate spline in the continuous range [0, nCoeff-1)
	__hostanddev__ double deriv(const double* coeff, double x);
	
	//! @brief Gradient propagation corresponding to value()
	//! @param E_value Input derivative with respect to value
	//! @param E_coeff
	__hostanddev__ void valueGrad(double E_value, double* E_coeff, double x);
}

//! @}

//--------------  Implementation -----------------
//!@cond

#include <cmath>

namespace QuinticSpline
{
	__hostanddev__ void getBernsteinCoeffs(const double* coeff, double x, double& tR, double& tL, double (&b)[6])
	{	int j = (int)x;
		tR = x - j; //right weight for interval
		tL = 1.-tR; //left weight for interval
		//Load blip coefficients:
		double c[6];
		for(int i=0; i<6; i++) c[i] = coeff[j+i];
		//Convert to Bernstein polynomial coefficients:
		b[0] = (1./66) * (c[0] + 26*c[1] + 66*c[2] + 26*c[3] + c[4]);
		b[1] = (1./33) * (8*c[1] + 33*c[2] + 18*c[3] + c[4]);
		b[2] = (2./33) * (2*c[1] + 15*c[2] + 12*c[3] + c[4]);
		b[3] = (2./33) * (c[1] + 12*c[2] + 15*c[3] + 2*c[4]);
		b[4] = (1./33) * (c[1] + 18*c[2] + 33*c[3] + 8*c[4]);
		b[5] = (1./66) * (c[1] + 26*c[2] + 66*c[3] + 26*c[4] + c[5]);
	}
	
	//Compute value of quintic spline
	__hostanddev__ double value(const double* coeff, double x)
	{	double tR, tL, b[6]; getBernsteinCoeffs(coeff, x, tR, tL, b);
		//Evaluate Bernstein polynomial by de Casteljau's reduction
		double c[5], d[4];
		for(int i=0; i<5; i++) c[i] = tL*b[i] + tR*b[i+1]; //5->4
		for(int i=0; i<4; i++) d[i] = tL*c[i] + tR*c[i+1]; //4->3
		for(int i=0; i<3; i++) c[i] = tL*d[i] + tR*d[i+1]; //3->2
		for(int i=0; i<2; i++) d[i] = tL*c[i] + tR*c[i+1]; //2->1
		return tL*d[0] + tR*d[1]; //1->0
	}
	
	//Compute derivative of quintic spline
	__hostanddev__ double deriv(const double* coeff, double x)
	{	double tR, tL, b[6]; getBernsteinCoeffs(coeff, x, tR, tL, b);
		//Derivative by by de Casteljau's reduction
		double c[5], d[4];
		for(int i=0; i<5; i++) c[i] = b[i+1] - b[i]; //5->4
		for(int i=0; i<4; i++) d[i] = tL*c[i] + tR*c[i+1]; //4->3
		for(int i=0; i<3; i++) c[i] = tL*d[i] + tR*d[i+1]; //3->2
		for(int i=0; i<2; i++) d[i] = tL*c[i] + tR*c[i+1]; //2->1
		return 5.*(tL*d[0] + tR*d[1]); //1->0
	}
	
	//Gradient propagation corresponding to value
	//On the GPU, final results are accumulated using shared memory. (Uses 6 doubles per thread of the thread-block)
	#ifdef __CUDA_ARCH__
	extern __shared__ double shared_E_coeff[];
	#endif
	__hostanddev__ void valueGrad(double E_value, double* E_coeff, double x)
	{	int j = (int)x;
		double tR = x - j; //right weight for interval
		double tL = 1.-tR; //left weight for interval
		double b[6], c[6];
		//Backtrace de Casteljau's reduction:
		b[0]=tL; b[1]=tR; //0->1
		c[0]=0.; for(int i=0; i<2; i++) { c[i] += tL*b[i]; c[i+1] = tR*b[i]; } //1->2
		b[0]=0.; for(int i=0; i<3; i++) { b[i] += tL*c[i]; b[i+1] = tR*c[i]; } //2->3
		c[0]=0.; for(int i=0; i<4; i++) { c[i] += tL*b[i]; c[i+1] = tR*b[i]; } //3->4
		b[0]=0.; for(int i=0; i<5; i++) { b[i] += tL*c[i]; b[i+1] = tR*c[i]; } //4->5
		//Propagate from Bernstein coefficients 'b' to blip coefficients 'c':
		c[0] = (1./66) * (b[0]);
		c[1] = (1./66) * (26*b[0] + 16*b[1] + 8*b[2] + 4*b[3] + 2*b[4] + b[5]);
		c[2] = (1./33) * (33*b[0] + 33*b[1] + 30*b[2] + 24*b[3] + 18*b[4] + 13*b[5]);
		c[3] = (1./33) * (13*b[0] + 18*b[1] + 24*b[2] + 30*b[3] + 33*b[4] + 33*b[5]);
		c[4] = (1./66) * (b[0] + 2*b[1] + 4*b[2] + 8*b[3] + 16*b[4] + 26*b[5]);
		c[5] = (1./66) * (b[5]);
		//Accumulate E_coeff:
		#ifndef __CUDA_ARCH__
			for(int i=0; i<6; i++) E_coeff[j+i] += E_value * c[i];
		#else
			int iThread = threadIdx.x;
			for(int i=0; i<6; i++) shared_E_coeff[iThread*6+i] = E_value * c[i];
			//Accumulate results to first thread:
			int extent = blockDim.x/2;
			int stride = (blockDim.x+1)/2;
			while(extent)
			{	__syncthreads();
				if(iThread<extent)
					for(int i=0; i<6; i++) shared_E_coeff[iThread*6+i] += shared_E_coeff[(stride+iThread)*6+i];
				extent = stride/2;
				stride = (stride+1)/2;
			}
			__syncthreads();
			//Save to the global memory:
			if(iThread==0)
				for(int i=0; i<6; i++) E_coeff[j+i] += shared_E_coeff[i];
		#endif
	}
}
//! @endcond

#endif // JDFTX_CORE_SPLINE_H
