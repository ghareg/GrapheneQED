#include "opermask.h"
#include <iostream>
#include <gsl/gsl_sf_laguerre.h>

int constrDif(Occup const* const* qstates, Count lftNum, Count rgtNum, Index* lft, Index* rgt)
{
	if (lftNum == rgtNum) {
		return 0;
	}
	const Occup* lftOccup = qstates[lftNum];
	const Occup* rgtOccup = qstates[rgtNum];
	int lftCount = 0;
	int rgtCount = 0;
	Index lftP = 0;
	Index rgtP = 0;
	while(lftP < Ne && rgtP < Ne) {
		if (lftOccup[lftP] == rgtOccup[rgtP]) {
			++lftP;
			++rgtP;
		}
		else if (lftOccup[lftP] > rgtOccup[rgtP]) {
			rgtCount++;
			if (rgtCount == 3) {
				return -1;
			}
			rgt[rgtCount - 1] = rgtP;
			++rgtP;
		}
		else {
			lftCount++;
			if (lftCount == 3) {
				return -1;
			}
			lft[lftCount - 1] = lftP;
			++lftP;
		}
	}
	if (lftP < Ne) {
		while(lftP < Ne) {
			lftCount++;
			if (lftCount == 3) {
				return -1;
			}
			lft[lftCount - 1] = lftP;
			++lftP;
		}
	}
	else {
		while(rgtP < Ne) {
			rgtCount++;
			if (rgtCount == 3) {
				return -1;
			}
			rgt[rgtCount - 1] = rgtP;
			++rgtP;
		}
	}
	
	return lftCount;
}

complex<double> MatElement(Count lftNum, Count rgtNum, const Occup qnummap[][2], Occup const* const* qstates, double tem, double det, const complex<double>* mee, double* sasNumState)
{
	Index lft[2];
	Index rgt[2];
	int diffNum = constrDif(qstates, lftNum, rgtNum, lft, rgt);
	if (diffNum == -1) {
	   return complex<double>(0.0, 0.0);
	}

	Occup const* lftOccup = qstates[lftNum];
	Occup const* rgtOccup = qstates[rgtNum];
	double OneMatElem = 0.0;
	if (diffNum == 0) {
		for (int i = 0; i < Ne; ++i) {
			OneMatElem += OneBodyElement(qnummap[lftOccup[i]][0], qnummap[lftOccup[i]][1], qnummap[rgtOccup[i]][0], qnummap[rgtOccup[i]][1], tem, det);
		}
		sasNumState[2] = SzState(lftOccup, qnummap);
	}
	else if (diffNum == 1) {
		if (qnummap[lftOccup[lft[0]]][1] == qnummap[rgtOccup[rgt[0]]][1]) {
			OneMatElem += sign(lft[0] + rgt[0]) * OneBodyElement(qnummap[lftOccup[lft[0]]][0], qnummap[lftOccup[lft[0]]][1], qnummap[rgtOccup[rgt[0]]][0], qnummap[rgtOccup[rgt[0]]][1], tem, det);
			sasNumState[0] = sign(lft[0] + rgt[0]);
			sasNumState[1] = qnummap[lftOccup[lft[0]]][0] == ninit ? sasNumState[0] : -sasNumState[0];
		}
	}
	
	complex<double> TwoMatElem(0.0, 0.0);
	if (diffNum == 0) {
		for (int i = 0; i < Ne; ++i) {
			for (int j = i + 1; j < Ne; ++j) {
				int n1 = qnummap[lftOccup[j]][0];
				int j1 = qnummap[lftOccup[j]][1];
				int n2 = qnummap[lftOccup[i]][0];
				int j2 = qnummap[lftOccup[i]][1];
				TwoMatElem += TwoBodyElement(n1, j1, n2, j2, n2, j2, n1, j1, mee) - TwoBodyElement(n1, j1, n2, j2, n1, j1, n2, j2, mee);
			}
		}
	}
	else if (diffNum == 1) {
		for (int i  = 0; i < lft[0]; ++i) {
			if (i < rgt[0]) {
				int n1 = qnummap[lftOccup[lft[0]]][0];
				int j1 = qnummap[lftOccup[lft[0]]][1];
				int n2 = qnummap[lftOccup[i]][0];
				int j2 = qnummap[lftOccup[i]][1];
				int n3 = qnummap[rgtOccup[i]][0];
				int j3 = qnummap[rgtOccup[i]][1];
				int n4 = qnummap[rgtOccup[rgt[0]]][0];
				int j4 = qnummap[rgtOccup[rgt[0]]][1];
				TwoMatElem += sign(lft[0] + rgt[0]) * 
					(TwoBodyElement(n1, j1, n2, j2, n3, j3, n4, j4, mee) - TwoBodyElement(n1, j1, n2, j2, n4, j4, n3, j3, mee));
			}
			else {
				int n1 = qnummap[lftOccup[lft[0]]][0];
				int j1 = qnummap[lftOccup[lft[0]]][1];
				int n2 = qnummap[lftOccup[i]][0];
				int j2 = qnummap[lftOccup[i]][1];
				int n3 = qnummap[rgtOccup[rgt[0]]][0];
				int j3 = qnummap[rgtOccup[rgt[0]]][1];
				int n4 = qnummap[rgtOccup[i + 1]][0];
				int j4 = qnummap[rgtOccup[i + 1]][1];
				TwoMatElem += sign(lft[0] + rgt[0] + 1) * 
					(TwoBodyElement(n1, j1, n2, j2, n3, j3, n4, j4, mee) - TwoBodyElement(n1, j1, n2, j2, n4, j4, n3, j3, mee));
			}
		}
		for (int i = lft[0] + 1; i < Ne; ++i) {
			if (i <= rgt[0]) {
				int n1 = qnummap[lftOccup[i]][0];
				int j1 = qnummap[lftOccup[i]][1];
				int n2 = qnummap[lftOccup[lft[0]]][0];
				int j2 = qnummap[lftOccup[lft[0]]][1];
				int n3 = qnummap[rgtOccup[i - 1]][0];
				int j3 = qnummap[rgtOccup[i - 1]][1];
				int n4 = qnummap[rgtOccup[rgt[0]]][0];
				int j4 = qnummap[rgtOccup[rgt[0]]][1];
				TwoMatElem += sign(lft[0] + rgt[0] + 1) * 
					(TwoBodyElement(n1, j1, n2, j2, n3, j3, n4, j4, mee) - TwoBodyElement(n1, j1, n2, j2, n4, j4, n3, j3, mee));
			}
			else {
				int n1 = qnummap[lftOccup[i]][0];
				int j1 = qnummap[lftOccup[i]][1];
				int n2 = qnummap[lftOccup[lft[0]]][0];
				int j2 = qnummap[lftOccup[lft[0]]][1];
				int n3 = qnummap[rgtOccup[rgt[0]]][0];
				int j3 = qnummap[rgtOccup[rgt[0]]][1];
				int n4 = qnummap[rgtOccup[i]][0];
				int j4 = qnummap[rgtOccup[i]][1];
			TwoMatElem += sign(lft[0] + rgt[0]) * 
					(TwoBodyElement(n1, j1, n2, j2, n3, j3, n4, j4, mee) - TwoBodyElement(n1, j1, n2, j2, n4, j4, n3, j3, mee));
			}
		}
	}
	else {
		int n1 = qnummap[lftOccup[lft[1]]][0];
		int j1 = qnummap[lftOccup[lft[1]]][1];
		int n2 = qnummap[lftOccup[lft[0]]][0];
		int j2 = qnummap[lftOccup[lft[0]]][1];
		int n3 = qnummap[rgtOccup[rgt[0]]][0];
		int j3 = qnummap[rgtOccup[rgt[0]]][1];
		int n4 = qnummap[rgtOccup[rgt[1]]][0];
		int j4 = qnummap[rgtOccup[rgt[1]]][1];
		TwoMatElem += sign(lft[0] + lft[1] + rgt[0] + rgt[1]) * 
			(TwoBodyElement(n1, j1, n2, j2, n3, j3, n4, j4, mee) - TwoBodyElement(n1, j1, n2, j2, n4, j4, n3, j3, mee));
	}

	return OneMatElem + TwoMatElem;
}

int SzState(const Occup* occup, const Occup qnummap[][2])
{
	int upCount = 0;
	for (int i = 0; i < Ne; ++i) {
		if (qnummap[occup[i]][0] == ninit) {
			upCount++;
		}
	}
	return 2 * upCount - Ne;
}	

double OneBodyElement(int n, int j, int np, int jp, double tem, double det)
{
	if (n == np) {
		if (n == ninit) {
			return -det;
		}
		else {
			return det;
		}
	}
	else if (n - np == 1 || np - n == 1) {
		return tem;
	}

	return 0.0;
}

double factorial(int n)
{
	double sum = 1;
	for(int i = 1; i <=n; ++i) {
		sum *= i;
	}
	return sum;
}

complex<double> Bs(int n1, int n2, int kx, int ky)
{
	double ellbLy=sqrt(lambda/(2*Pi*Ns));
	complex<double> result(0.0,0.0);
	if(n1 <= n2) {
  		result = sqrt(factorial(n1)/factorial(n2))*pow(sqrt(2) * Pi * ellbLy, (n2-n1)) * pow(kx/lambda + (double)ky * I, n2-n1)*gsl_sf_laguerre_n (n1, n2-n1,2 * Pi * Pi * ellbLy * ellbLy *(kx*kx/(lambda*lambda) + ky * ky));
	}
  	else {
		result = sqrt(factorial(n2)/factorial(n1))*pow(sqrt(2) * Pi * ellbLy, (n1-n2)) * pow(-kx/lambda+ (double)ky * I, n1-n2)*gsl_sf_laguerre_n (n2, n1-n2,2 * Pi * Pi * ellbLy * ellbLy *(kx*kx/(lambda*lambda) + ky * ky));
	}
  
	return result;
}

complex<double> Bf(int n1, int n2, int kx, int ky)
{
	
	if (n1 == 0 || n2 == 0) {
	   return Cnp(n1) * Cnp(n2) * Bs(n1, n2, kx, ky);
	}
	else {
		return Cnm(n1) * Cnm(n2) * Bs(n1 - 1, n2 -1, kx, ky) + Cnp(n1) * Cnp(n2) * Bs(n1, n2, kx, ky);
	}
}

complex<double> Mee(int n1, int j1, int n2, int j2, int n3, int j3, int n4, int j4)
{
  complex<double> result=0;
  double ellbLy=sqrt(lambda/(2*Pi*Ns));
  for(int kx=-lim; kx<=lim; kx++)
  {
    for(int ky=-lim; ky<=lim; ky++)
    {
      if((j1-j4-kx)%Ns==0 && (kx!=0 || ky!=0))
		  result += (ellbLy / lambda) * (1/sqrt(kx*kx/(lambda*lambda)+ky*ky)) * exp(-2*Pi*Pi*ellbLy*ellbLy*(kx*kx/(lambda*lambda)+ky*ky))*Bf(n1, n4, kx, ky)*Bf(n2, n3, -kx,-ky)*exp((2*Pi*ky*(j1-j3)/Ns)*I);
    }
  }
  
  return result;
}

complex<double> MeeMat(int n1, int j1, int n2, int j2, int n3, int j3, int n4, int j4, const complex<double>* mee)
{
  if(mee)
  {
  const int nsmax = 2 * Ns;;
  int lft1 = (n1 - ninit) * Ns + j1;
  int lft2 = (n2 - ninit) * Ns + j2;
  int rgt1 = (n3 - ninit) * Ns + j3;
  int rgt2 = (n4 - ninit) * Ns + j4;
  
  return mee[lft1*nsmax*nsmax*nsmax+lft2*nsmax*nsmax+rgt1*nsmax+rgt2];
  }

  return complex<double>(0.0,0.0);
}

complex<double> TwoBodyElement(int n1, int j1, int n2, int j2, int n3, int j3, int n4, int j4, const complex<double>* mee)
{
	if ((j1 + j2 - j3 - j4)%Ns==0) {
		return MeeMat(n1, j1, n2, j2, n3, j3, n4, j4, mee);
	}
	
	return complex<double>(0.0, 0.0);
}

double calculateElNum(Count i,  const Occup qnummap[][2], Occup const* const* qstates)
{
	const Occup* state = qstates[i];
	int elNum = 0;
	for (int n = 0; n < Ne; ++n) {
		if (qnummap[state[n]][0] == ninit) {
			++elNum;
		}
	}
	return elNum;
}

