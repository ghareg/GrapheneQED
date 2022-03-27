#ifndef OPERMASK_H_
#define OPERMASK_H_

#include <complex> 
#include "basis.h"

typedef int Index;
using std::complex;

extern const double Pi;
extern const int Ne;
extern const int Ns;
extern const int ninit;
const double lambda=0.98;
const complex<double> I(0,1);
const int lim=80;

inline double sign(int i) {return i % 2 ? -1.0: 1.0;}
inline double Cnp (int n) {return n == 0 ? 1.0 : sqrt(0.5);}
inline double Cnm (int n) {return n == 0 ? 0 : sqrt(0.5);}

int constrDif(Occup const* const* qstates, Count lftNum, Count rgtNum, Index* lft, Index* rgt);
complex<double> MatElement(Count lftNum, Count rgtNum, const Occup qnummap[][2], Occup const* const* qstates, double tem, double det, const complex<double>* mee, double* sasNumState);
int SzState(const Occup* occup, const Occup qnummap[][2]);

double factorial(int n);
complex<double> Bs(int n1, int n2, int kx, int ky);
complex<double> Bf(int n1, int n2, int kx, int ky);
complex<double> Mee(int n1, int j1, int n2, int j2, int n3, int j3, int n4, int j4);
complex<double> MeeMat(int n1, int j1, int n2, int j2, int n3, int j3, int n4, int j4, const complex<double>* mee);
double OneBodyElement(int n, int j, int np, int jp, double tem, double det);
complex<double> TwoBodyElement(int n1, int j1, int n2, int j2, int n3, int j3, int n4, int j4, const complex<double>* mee);

double calculateElNum(Count i,  const Occup qnummap[][2], Occup const* const* qstates);

template <typename arrayT>
void resize(arrayT*& mat, Count currSize, Count newSize)
{
	arrayT* oldmat = mat;
	mat = new arrayT[newSize];
	for (Count i = 0; i < currSize; ++i) {
		mat[i] = oldmat[i];
	}
	delete[] oldmat;
}

#endif
